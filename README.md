# bmssp-vs-dijkstra

Implementations, experiments and data for the paper

> **A Comprehensive Theoretical and Experimental Comparison of Dijkstra's Algorithm and
> Duan et al.'s SSSP Algorithms**

The repository contains independent C++17 implementations of

* Dijkstra's algorithm (binary heap with lazy deletion, and indexed heap with DecreaseKey),
* the O(m log^(2/3) n) algorithm of Duan, Mao, Mao, Shu and Yin (STOC 2025) — `bmssp25`,
* the O(m sqrt(log n) + sqrt(mn log n loglog n)) algorithm of Duan, Mao, Shu and Yin (ICALP 2026) — `bmssp26`,

together with the drivers that produce every measurement, table and figure of the paper.

## Layout

| Path | Content |
|------|---------|
| `research-paper.pdf` | the paper (IEEE conference format) |
| `explanation.pdf` | a plain-language guide to the algorithms, the code, the experiments and the results |
| `figure_data.txt` | every value plotted in the paper's figures, tab-separated, for redrawing them |
| `paper/` | LaTeX source of the paper; `figures/`, `tables/` and `numbers.tex` are generated |
| `explanation/` | LaTeX source of the guide |
| `scripts/src/` | C++ sources (`dijkstra.hpp`, `bmssp25.hpp`, `bmssp26.hpp`, `graphs.hpp`, `common.hpp`, `main.cpp`) |
| `scripts/Makefile` | builds `bin/sssp` (timing), `bin/sssp_count` (operation counts), `bin/sssp_debug` (sanitizers) |
| `scripts/fetch_data.sh` | downloads the five SSTD 2005 road networks |
| `scripts/run_experiments.py` | runs all experiments and writes `results/*.csv` |
| `scripts/make_figures.py` | turns `results/*.csv` into figures, tables and the numbers quoted in the text |
| `scripts/build_pdfs.sh` | runs `make_figures.py` and typesets both PDFs |
| `results/` | the raw measurements behind the paper, `summary.json`, `environment.json`, validation witnesses |

## Reproducing the paper

Requirements: a C++17 compiler, `make`, `git`, Python 3 (with `matplotlib` and `numpy` for the figures)
and, for the PDFs, [Tectonic](https://tectonic-typesetting.github.io) or any LaTeX distribution with `IEEEtran`.

```sh
sh scripts/fetch_data.sh snap                 # road networks -> data/sstd2005, SNAP graphs -> data/snap
make -C scripts                               # bin/sssp and bin/sssp_count
make -C scripts check                         # sanitizer build + randomized validation
python3 scripts/run_experiments.py            # all measurements (about one hour) -> results/
python3 -m venv .venv && .venv/bin/pip install matplotlib numpy
sh scripts/build_pdfs.sh                      # figures, tables, numbers.tex, then both PDFs
```

`build_pdfs.sh` runs `scripts/make_figures.py` (writes `paper/figures`, `paper/tables`, `paper/numbers.tex`,
`figure_data.txt` and `results/summary.json`) and typesets `paper/main.tex` and `explanation/explanation.tex`.

`run_experiments.py` accepts experiment names (`info validate real variants reduce params scaling
density weights memory hubreduce snap`), `--force` to overwrite existing results and `--counts-only` to redo
only the operation counts. Timings depend on the machine;
operation counts and validation results are deterministic for the fixed seeds.

## Using the solver directly

```sh
scripts/bin/sssp info  "sstd:data/sstd2005/Datasets - SSTD 2005/OL Road Network's Edges (Edge ID, Start Node ID, End Node ID, L2 Distance).txt"
scripts/bin/sssp run   rand:65536:3:1 --algo bmssp26 --sources 3 --reps 5
scripts/bin/sssp_count run grid:256:1 --algo bmssp25 --k 2 --t 6 --sources 1
scripts/bin/sssp_count validate --trials 20000 --maxn 100 --seed 1
```

Graph specifications: `sstd:<edge file>`, `snap:<edge list[.gz]>:<unit|rand>`, `rand:<n>:<deg>:<seed>[:<wmax>]`, `grid:<side>:<seed>`,
`ba:<n>:<d>:<seed>`, `hub:<n>:<deg>:<seed>`. Options: `--k`, `--t`, `--delta` (parameters),
`--reduce <delta>` (run on the degree-reduced graph), `--nofix` and `--literal-base` (pseudocode as
printed), `--source0`, `--min-reach <fraction>` (only sources that reach that share of the vertices). Every run is checked against Dijkstra's algorithm; the `ok` column must be 1.

## Data

The road networks are the SSTD 2005 collection (Li, Cheng, Hadjieleftheriou, Kollios and Teng,
"On Trip Planning Queries in Spatial Databases", SSTD 2005), downloaded from
<https://github.com/mehdi-ansari/Real-dataset-shortest-path-sample>. They are not redistributed here;
`scripts/fetch_data.sh` fetches the exact commit used.

Six further graphs (web-Stanford, web-Google, roadNet-PA/TX/CA, soc-Pokec) come from the SNAP collection
(J. Leskovec and A. Krevl, <http://snap.stanford.edu/data>); `scripts/fetch_data.sh snap` downloads them
(about 200 MB). Without them `run_experiments.py` simply skips the `snap` experiment.
