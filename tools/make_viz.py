"""Embeds the generated results into viz/index.html, so the page works offline as a single file.

Owner: Shahana. Standard library only.

    python tools/make_viz.py          # rewrite the data block in viz/index.html
    python tools/make_viz.py --check  # exit 1 if the data block is out of date

Reads (anything missing is skipped with a warning, so this works before the bench has run):
  results/tests.csv
  results/bench_scaling.csv, bench_flow.csv, bench_improve.csv, bench_sobel.csv, bench_carve.csv,
  results/seam_counts.csv            (lines starting with '#' are kept as the machine note)
  results/demo/<image>/seams.json     {"w","h","minW","seams"} for the resize slider (packed seams)
  results/demo/<image>/original.*     the picture the seams were cut from, embedded unchanged
  results/demo/<image>/summary.csv    method,input_w,input_h,out_w,out_h,removed_energy,ms
  results/demo/<image>/<method>.jpg   one per summary row (.png for a .png input), plus input and energy
and writes them as one JSON object between /*DATA-BEGIN*/ and /*DATA-END*/ in viz/index.html.
Images are ordered from easiest to hardest by their difficulty: the average energy of the pixels
the DP carve removed (removed energy / removed pixels in the dp row of summary.csv).
"""

import argparse
import base64
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
HTML = ROOT / "viz" / "index.html"
BEGIN, END = "/*DATA-BEGIN*/", "/*DATA-END*/"
SIZE_WARN_BYTES = 8_000_000  # each photo adds about 1 MB (its JPEGs and packed seams, base64)
BENCH = ["bench_scaling", "bench_flow", "bench_improve", "bench_sobel", "bench_carve", "seam_counts"]
# Every bench CSV starts its note with the same machine line; anything after it is a caveat for that file.
NOTE_LABEL = {"bench_flow": "max flow", "bench_carve": "carve table"}
EXTRAS = ["input", "energy"]
MIME = {".jpg": "image/jpeg", ".jpeg": "image/jpeg", ".png": "image/png"}


def warn(msg):
    print(f"warning: {msg}", file=sys.stderr)


def number(v):
    """CSV cells become int or float when they look like one; everything else stays a string."""
    for cast in (int, float):
        try:
            return cast(v)
        except ValueError:
            pass
    return v


def read_csv(path):
    """Returns (rows, notes): rows as dicts with numeric cells converted, notes = '#' comment lines."""
    notes, lines = [], []
    with open(path, encoding="utf-8", newline="") as f:
        for line in f:
            (notes if line.startswith("#") else lines).append(line)
    rows = [{k: number(v) for k, v in row.items()} for row in csv.DictReader(lines)]
    return rows, [n.lstrip("#").strip() for n in notes]


def data_uri(path):
    return f"data:{MIME[path.suffix.lower()]};base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def picture(d, stem):
    """results/demo/<image>/<stem>.jpg or .png, whichever the demo wrote; None if neither."""
    return next((d / f"{stem}{ext}" for ext in MIME if (d / f"{stem}{ext}").exists()), None)


def merge_notes(notes):
    """[(csv name, note)] -> the shared machine line once, then each file's extra caveat with its label."""
    if not notes:
        return []
    machine = min((n for _, n in notes), key=len)
    out = [machine]
    for name, n in notes:
        extra = n[len(machine):].lstrip(", ") if n.startswith(machine) else n
        line = f"{NOTE_LABEL.get(name, name)}: {extra}" if extra else ""
        if line and line not in out:
            out.append(line)
    return out


def collect():
    data = {"tests": [], "bench": {}, "notes": [], "images": []}
    tests = RESULTS / "tests.csv"
    if tests.exists():
        data["tests"] = [{k: r[k] for k in ("id", "name", "status", "ms")} for r in read_csv(tests)[0]]
    else:
        warn("results/tests.csv not found (run the tests)")
    for name in BENCH:
        path = RESULTS / f"{name}.csv"
        if path.exists():
            rows, notes = read_csv(path)
            data["bench"][name] = rows
            data["notes"] += [(name, n) for n in notes]
        else:
            warn(f"results/{name}.csv not found (run the bench)")
    demo = RESULTS / "demo"
    for d in sorted(p for p in demo.glob("*") if p.is_dir()) if demo.exists() else []:
        seams, summary, original = d / "seams.json", d / "summary.csv", picture(d, "original")
        if not seams.exists() or not summary.exists() or original is None:
            warn(f"{d.relative_to(ROOT)}: seams.json, summary.csv or original.* missing, skipped")
            continue
        slider = {**json.loads(seams.read_text(encoding="utf-8")), "src": data_uri(original)}
        image = {"name": d.name, "slider": slider, "gallery": [], "img": {}, "difficulty": None}
        for row in read_csv(summary)[0]:
            if row["method"] == "dp" and row["out_w"] < row["input_w"]:
                image["difficulty"] = round(row["removed_energy"] / ((row["input_w"] - row["out_w"]) * row["input_h"]), 1)
            pic = picture(d, row["method"])
            if pic:
                image["gallery"].append({**row, "src": data_uri(pic)})
            else:
                warn(f"{d.relative_to(ROOT)}/{row['method']}.jpg or .png missing, skipped")
        for name in EXTRAS:
            if picture(d, name):
                image["img"][name] = data_uri(picture(d, name))
        data["images"].append(image)
    data["notes"] = merge_notes(data["notes"])
    if not data["images"]:
        warn("no results/demo/<image>/ folders found (run seamcarve demo)")
    # Easiest first; an image without a dp row goes last, then by name so the order is stable.
    data["images"].sort(key=lambda im: (im["difficulty"] is None, im["difficulty"] or 0, im["name"]))
    return data


def render(html, data):
    if html.count(BEGIN) != 1 or html.count(END) != 1:
        sys.exit(f"{HTML}: expected exactly one {BEGIN} ... {END} block")
    head, rest = html.split(BEGIN)
    _, tail = rest.split(END)
    block = "const DATA = " + json.dumps(data, separators=(",", ":"), sort_keys=True) + ";"
    return head + BEGIN + block + END + tail


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="exit 1 if viz/index.html is out of date")
    args = ap.parse_args()
    if not (ROOT / "results").is_dir():
        # Rewriting the page now would replace its embedded data with an empty set.
        print("results/ not found: run 'make results' first (the page keeps the data it already has)")
        return 1
    html = HTML.read_text(encoding="utf-8")
    new = render(html, collect())
    if args.check:
        if new != html:
            print("viz/index.html is out of date: run python tools/make_viz.py")
            return 1
        print("viz/index.html is up to date")
        return 0
    HTML.write_text(new, encoding="utf-8", newline="\n")
    size = len(new.encode("utf-8"))
    print(f"wrote {HTML.relative_to(ROOT)} ({size / 1e6:.2f} MB)")
    if size > SIZE_WARN_BYTES:
        warn(f"viz/index.html is over {SIZE_WARN_BYTES / 1e6:.0f} MB; shrink the demo inputs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
