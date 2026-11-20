"""Mutation testing: plant one realistic bug at a time and check that the test suite catches it.

Owner: Yan Xin. Standard library only.

    python tools/mutation_test.py [--quick] [--keep]

For each mutant: copy include/ and tests/ into mutants_tmp/<id>/, apply one text replacement (it must
match exactly once, otherwise SETUP-ERROR), build run_tests with the project flags and run it.
A failing run (or a crash or timeout) means KILLED; a passing run means SURVIVED.
The unmutated copy runs first: if it fails, the script stops with BASELINE-FAILED, because then
every mutant would look killed.
Mutants whose file does not exist yet are reported as NOT-YET and do not fail the run.
Writes results/mutation.csv and exits 1 if any mutant survived or could not be set up.
"""

import argparse
import concurrent.futures
import csv
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
WORK = ROOT / "mutants_tmp"

# Keep in sync with: Makefile and build_all.bat.
CXX = os.environ.get("CXX", "g++")
FLAGS = ["-std=c++17", "-O2", "-Wall", "-Wextra", "-pthread"]
LDFLAGS = ["-static"] if os.name == "nt" else []
EXE = "run_tests.exe" if os.name == "nt" else "run_tests"
RUN_TIMEOUT_S = 300

# The "old" text must appear exactly once in the file. The tests that kill each mutant are listed in results/mutation.csv.
MUTANTS = [
    dict(id="M1", file="include/dp.h", desc="parent window misses the right neighbour",
         old="    if (up_right < best) best = up_right, d = +1;\n",
         new=""),
    dict(id="M2", file="include/dp.h", desc="tie rule uses <=, so the rightmost parent wins",
         old="    if (up < best) best = up, d = 0;\n    if (up_right < best) best = up_right, d = +1;",
         new="    if (up <= best) best = up, d = 0;\n    if (up_right <= best) best = up_right, d = +1;"),
    dict(id="M3", file="include/dp.h", desc="sentinel columns are 0 instead of +inf",
         old="M = Grid<int64_t>(h, w + 2, INF);",
         new="M = Grid<int64_t>(h, w + 2, 0);"),
    dict(id="M4", file="include/dp.h", desc="backtrack starts from the last column instead of the argmin",
         old="return detail::backtrack(P, detail::argmin_first(&M.at(e.h - 1, 1), e.w));",
         new="return detail::backtrack(P, e.w - 1);"),
    dict(id="M5", file="include/dp.h", desc="rolling DP updates one row in place",
         old="std::vector<int64_t>& row = buf[i & 1];",
         new="std::vector<int64_t>& row = buf[(i - 1) & 1];"),
    dict(id="M6", file="include/energy.h", desc="local update band is [lo, hi] instead of [lo-1, hi]",
         old="int c0 = std::max(lo - 1, 0)",
         new="int c0 = std::max(lo, 0)"),
    dict(id="M7", file="include/energy.h", desc="separable smoothing uses weights [1,1,1]",
         old="s[c] = left + 2 * mid + right;",
         new="s[c] = left + mid + right;"),
    dict(id="M8", file="include/greedy.h", desc="greedy starts at column 0 instead of the top-row argmin",
         old="s[0] = int(std::min_element(top, top + e.w) - top);",
         new="s[0] = 0;"),
    dict(id="M9", file="include/graphcut.h", desc="no padding column",
         old="sg.W = pad ? e.w + 1 : e.w;",
         new="sg.W = e.w;"),
    dict(id="M10", file="include/graphcut.h", desc="missing INF backward edge in each row chain",
         old="            g.add_edge(sg.node(i, j + 1), sg.node(i, j), INF);\n",
         new=""),
    dict(id="M11", file="include/graphcut.h", desc="residual update forgets cap[e^1] += x",
         old="        g.cap[e ^ 1] += x;\n",
         new=""),
    dict(id="M12", file="include/graphcut.h", desc="seam = leftmost-in-T instead of rightmost-in-S",
         old="        for (int j = 0; j < sg.W; ++j)\n            if (S[sg.node(i, j)]) s[i] = j;",
         new="        for (int j = sg.W - 1; j >= 0; --j)\n            if (!S[sg.node(i, j)]) s[i] = j;"),
    dict(id="M13", file="include/carve.h", desc="horizontal pass forgets to transpose idx back",
         old="        img = transpose(img);\n        idx = transpose(idx);\n    }",
         new="        img = transpose(img);\n    }"),
    dict(id="M14", file="include/core.h", desc="remove_seam keeps the seam pixel and drops its right neighbour",
         old="        std::copy(src, src + s[r], dst);\n        std::copy(src + s[r] + 1, src + g.w, dst + s[r]);",
         new="        std::copy(src, src + s[r] + 1, dst);\n        std::copy(src + s[r] + 2, src + g.w, dst + s[r] + 1);"),
    dict(id="M15", file="include/carve.h", desc="Gradient energy is transposed instead of recomputed",
         old="E = o.energy == EnergyKind::Gradient ? compute_energy(G, o.energy, o.threads) : transpose(E);",
         new="E = transpose(E);"),
]


def prepare(name, mutant=None):
    """Copies the sources into WORK/name and applies the mutant. Returns an error string or None."""
    d = WORK / name
    shutil.rmtree(d, ignore_errors=True)
    shutil.copytree(ROOT / "include", d / "include")
    shutil.copytree(ROOT / "tests", d / "tests")
    if mutant:
        path = d / mutant["file"]
        text = path.read_text(encoding="utf-8")
        count = text.count(mutant["old"])
        if count != 1:
            return f"'old' text matches {count} times in {mutant['file']}"
        path.write_text(text.replace(mutant["old"], mutant["new"]), encoding="utf-8")
    return None


def build_and_run(name):
    """Returns (status, detail): status is PASS, FAIL, BUILD-ERROR or TIMEOUT."""
    d = WORK / name
    exe = d / EXE
    cmd = [CXX, *FLAGS, "-I", str(d / "include"), "-I", str(ROOT / "third_party"),
           str(d / "tests" / "run_tests.cpp"), "-o", str(exe), *LDFLAGS]
    b = subprocess.run(cmd, capture_output=True, text=True)
    if b.returncode != 0:
        return "BUILD-ERROR", b.stderr.strip().splitlines()[0] if b.stderr.strip() else "compiler failed"
    try:
        # cwd = the mutant's folder, so parallel runs do not share tests_tmp/.
        r = subprocess.run([str(exe), "--csv", str(d / "tests.csv")], cwd=d, capture_output=True,
                           text=True, timeout=RUN_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        return "TIMEOUT", f"no result after {RUN_TIMEOUT_S} s"
    if r.returncode == 0:
        return "PASS", ""
    failed = re.findall(r"^FAIL\s+(\S+)", r.stdout, flags=re.M)
    return "FAIL", " ".join(failed) if failed else f"exit code {r.returncode} (crash)"


def run_mutant(m):
    if not (ROOT / m["file"]).exists():
        return m, "NOT-YET", "file not present"
    err = prepare(m["id"], m)
    if err:
        return m, "SETUP-ERROR", err
    status, detail = build_and_run(m["id"])
    if status == "PASS":
        return m, "SURVIVED", ""
    if status == "BUILD-ERROR":
        return m, "SETUP-ERROR", "mutant does not compile: " + detail
    return m, "KILLED", detail


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--quick", action="store_true", help="only the first 5 mutants (CI)")
    ap.add_argument("--keep", action="store_true", help="keep mutants_tmp/ for inspection")
    args = ap.parse_args()
    mutants = MUTANTS[:5] if args.quick else MUTANTS

    WORK.mkdir(exist_ok=True)
    prepare("baseline")
    status, detail = build_and_run("baseline")
    if status != "PASS":
        print(f"BASELINE-FAILED: the unmutated code does not pass ({status}: {detail})")
        return 1
    print("baseline: all tests pass")

    with concurrent.futures.ThreadPoolExecutor(max_workers=os.cpu_count() or 2) as pool:
        results = list(pool.map(run_mutant, mutants))

    out = ROOT / "results" / "mutation.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["id", "file", "description", "status", "killed_by"])
        for m, st, detail in results:
            w.writerow([m["id"], m["file"], m["desc"], st, detail if st == "KILLED" else ""])

    print(f"{'id':<5} {'file':<14} {'status':<12} detail")
    for m, st, detail in results:
        print(f"{m['id']:<5} {m['file']:<14} {st:<12} {detail or m['desc']}")
    counts = {s: sum(1 for _, st, _ in results if st == s) for s in ("KILLED", "SURVIVED", "SETUP-ERROR", "NOT-YET")}
    active = len(results) - counts["NOT-YET"]
    print(f"{counts['KILLED']}/{active} KILLED" + (f", {counts['NOT-YET']} not yet present" if counts["NOT-YET"] else ""))

    if not args.keep:
        shutil.rmtree(WORK, ignore_errors=True)
    return 1 if counts["SURVIVED"] or counts["SETUP-ERROR"] else 0


if __name__ == "__main__":
    sys.exit(main())
