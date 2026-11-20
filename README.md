# Seam Carving - CSD3131 Assignment 2

Content-aware image resizing by seam carving, in C++17 with no libraries beyond image reading and writing.
Seams are found by dynamic programming (the main solution), by a greedy rule, and as a minimum cut with three max-flow algorithms, and all of them are compared on the same inputs.

## Marking this? Start here

| Step | Command or link | What you should see |
|---|---|---|
| Build and run every test (Windows, MinGW) | `build_all.bat` (runs `mingw32-make test`) | `ALL 46 TESTS PASSED`, then `DONE: make test succeeded` |
| Same on Linux or macOS | `make` (or `make test`) | same |
| Interactive visualiser (resize slider, method gallery, charts, tests) | open `viz/index.html` in a browser (works offline) | the page |
| Generated results quoted in the report (not stored in the repo) | `make results` (several minutes), plus `make mutation` for the mutation table | `results/`: CSVs and images, listed below |
| Demo video | TODO: link | |

## Quick usage

```bash
seamcarve carve images/starry_night.jpg -W 50% -o out.png --seams seams.png
seamcarve carve images/synthetic.png -W 60% -H 80% -m dinic -o dinic.png
seamcarve baseline images/synthetic.png -W 50% --kind scale -o scaled.png
seamcarve dumpdp --demo
seamcarve bench --quick
```

`seamcarve --help` lists every command and option.
Sizes are pixels (`300`) or percentages (`50%`); methods are `dp` (default), `rolling`, `memo`, `mt`, `greedy`, `brute`, `ff`, `ek` and `dinic`.
Exit codes: 0 success, 1 runtime error, 2 bad usage.

## Requirement traceability

| Brief item | Where it is done | Proof |
|---|---|---|
| 1a Energy map (gradient / Sobel) | `include/energy.h` | EN1-EN7 |
| 1b DP minimum seam, vertical and horizontal | `include/dp.h` + transpose in `include/carve.h` | DP1-DP9, CV6 |
| 1c Remove seam, update energy, repeat to target size | `include/carve.h` + `remove_and_update` in `include/energy.h` | CV1-CV5, EN5 |
| 2a(i) Number of seams is exponential | `count_seams` in `include/dp.h` | DP8, `results/seam_counts.csv` |
| 2a(ii) Recurrence over subproblems | SRTBOT comment in `include/dp.h`, `seamcarve dumpdp --demo` | DP6 (hand-filled table) |
| 2b Divide and conquer for the convolution | `sobel_direct`, `sobel_separable`, `sobel_dc` and threads in `include/energy.h` | EN4, EN7, `results/bench_sobel.csv` |
| 3 Greedy, compared with DP | `include/greedy.h` | GR1-GR4, `results/bench_carve.csv`, gallery |
| 4 Graph cut (min cut / max flow), compared | `include/graphcut.h` (Ford-Fulkerson, Edmonds-Karp, Dinic) | GC1-GC9, `results/bench_flow.csv` |
| 5b Time and space complexity | work counters + `seamcarve bench` | `results/bench_scaling.csv` (fitted slopes in the visualiser) |
| 5c Comparison with greedy and graph cut | bench `carve` suite + baselines | `results/bench_carve.csv`, `results/demo/` |
| 5d Improvements | rolling DP rows, local energy update, threads, Dinic | `results/bench_improve.csv` |
| 7a Bonus: visualise all seams, step by step | overlay and frames in `include/carve.h`, `viz/index.html` | CV5 |
| 7b Bonus: interactive resizing | removal-order map, packed seams + slider (original and resized side by side, next seam drawn) | CV4, CV8 |

## Architecture

```
            core.h      types, image I/O, seam utilities, transpose, timer, CSV, test framework
           /      \
     energy.h     dp.h          energy map and updates | brute force, memo, bottom-up, rolling, threads, counting
         |   \   /   \
         |  greedy.h  graphcut.h    greedy + baselines | seam graph, FF / EK / Dinic, certificate
          \     |      /
            carve.h             carving driver, both directions, overlay, frames, removal-order map, packed seams
               |
            bench.h             benchmark suites
               |
            main.cpp            the seamcarve CLI
```

Every module is one header holding its implementation and its own tests (`#ifdef SC_TESTS`).
There are two executables: `seamcarve` (from `src/main.cpp`) and `run_tests` (from `tests/run_tests.cpp`, which picks up whichever modules exist). `make` compiles each into the system temp folder, runs it and deletes it, so neither file ever sits in the project.

| Files | Committed by | Presents |
|---|---|---|
| project setup, `include/core.h`, `tests/run_tests.cpp`, `Makefile`, `build_all.bat`, `viz/index.html` | Henrey | slides 10-11, demo and conclusion |
| `include/energy.h`, `tools/make_viz.py`, `README.md` | Shahana | slides 1-2, introduction and seam carving |
| `include/dp.h`, `src/main.cpp` | Vera | slides 3-4, dynamic programming |
| `include/greedy.h`, `include/bench.h` | Asha | slides 5-6, divide and conquer and greedy |
| `include/graphcut.h` | Zafir | slide 7, graph cut |
| `include/carve.h`, `tools/mutation_test.py` | Yan Xin | slides 8-9, improvements and testing |

## Tests

```bash
build_all.bat                  # or make test: build + all tests (build_all.bat is a wrapper around make)
make test ARGS="--filter DP"   # one module (CORE, EN, DP, GR, GC, CV)
mingw32-make mutation          # or make mutation: mutation testing, a few minutes
```

There are 46 tests: CORE 7, EN 7, DP 9, GR 6, GC 9 and CV 8.
Each one checks its result by an independent route: brute force, a second algorithm, a hand-computed table, or a max-flow min-cut certificate.
Random tests use fixed seeds, so every failure is reproducible.
`tools/mutation_test.py` plants 15 realistic bugs one at a time and checks that the suite catches each one; `results/mutation.csv` lists which tests killed which mutant.

## Benchmarks and results

```bash
mingw32-make results           # or make results: tests, demo images, full benchmarks, visualiser data
seamcarve bench --quick        # about a minute
```

| File | Content |
|---|---|
| `results/tests.csv` | every test, status and time |
| `results/mutation.csv` | every mutant and the tests that killed it |
| `results/bench_scaling.csv` | per-seam time against n for every method |
| `results/bench_carve.csv` | full 50% carve per image and method: time and removed energy |
| `results/bench_flow.csv` | Ford-Fulkerson vs Edmonds-Karp vs Dinic |
| `results/bench_improve.csv` | memory, energy update modes, threads, max-flow solvers |
| `results/bench_sobel.csv` | direct vs separable vs divide-and-conquer vs threaded energy |
| `results/seam_counts.csv` | number of seams against 2^(h-1) |
| `results/demo/<image>/` | every method on the same input (320 px), energy map, slider data (full resolution) |

Timing protocol: one warm-up run, then the median of five timed runs, built with `-O2`; each CSV names the machine in its first line.

## Credits

- Image reading and writing: [stb_image and stb_image_write](https://github.com/nothings/stb) by Sean Barrett (public domain / MIT).
- Images: see `images/SOURCES.md`.
- Algorithms follow the CSD3131 lectures (dynamic programming, max flow) and Avidan and Shamir, *Seam Carving for Content-Aware Image Resizing* (2007).
- The MIT seam carving lecture used for the greedy example and the baselines: [mitmath/computational-thinking](https://github.com/mitmath/computational-thinking).
- The report contains the AI-usage declaration.
