#!/usr/bin/env python3
"""Turn the raw measurements in results/ into the figures, tables and numbers of the paper.

    .venv/bin/python scripts/make_figures.py

Reads   results/*.csv
Writes  paper/figures/*.pdf      (the four figures of the paper)
        explanation/figures/*.pdf (three further charts used only by the guide)
        figure_data.txt          (every plotted value as plain text, for redrawing the figures)
        paper/tables/*.tex       (table bodies)
        paper/numbers.tex        (every number quoted in the running text, as LaTeX macros)
        results/summary.json     (the same aggregates in machine-readable form)

Aggregation rule used everywhere: for each (graph, algorithm, source) take the median
over the repetitions; a ratio to Dijkstra is formed per source (same source, same graph)
and summarised over the sources by its geometric mean.
"""
import csv
import json
import math
import os
import statistics as st
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
FIG = os.path.join(ROOT, "paper", "figures")
TAB = os.path.join(ROOT, "paper", "tables")
GUIDEFIG = os.path.join(ROOT, "explanation", "figures")
os.makedirs(FIG, exist_ok=True)
os.makedirs(GUIDEFIG, exist_ok=True)
os.makedirs(TAB, exist_ok=True)

# Categorical colours (validated for colour-vision deficiency); every series also has
# its own marker / hatch so that the figures survive greyscale printing.
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d2"
STYLE = {
    "dijkstra":    dict(label="Dijkstra (binary heap)", color="#2a78d6", marker="o", ls="-", hatch=""),
    "dijkstra-dk": dict(label="Dijkstra (DecreaseKey)", color="#8a8983", marker="s", ls="--", hatch="...."),
    "bmssp25":     dict(label="BMSSP-25", color="#eb6834", marker="^", ls="-", hatch="////"),
    "bmssp26":     dict(label="BMSSP-26", color="#1baf7a", marker="D", ls="-", hatch="\\\\\\\\"),
}
ROADS = ["OL", "TG", "CAL", "SF", "NA"]

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "axes.linewidth": 0.6,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "text.color": INK, "grid.color": GRID, "grid.linewidth": 0.5, "lines.linewidth": 1.3,
    "lines.markersize": 4, "pdf.fonttype": 42, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "hatch.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False,
})
COL1, COL2 = 3.45, 7.1  # IEEE column / page width in inches


def read(name):
    path = os.path.join(RES, name + ".csv")
    if not os.path.exists(path):
        return []
    with open(path) as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k, v in r.items():
            if k in ("graph", "algo", "variant", "family", "name", "snapgraph", "weights"):
                continue
            try:
                r[k] = int(v)
            except ValueError:
                try:
                    r[k] = float(v)
                except ValueError:
                    pass
    return rows


def gmean(xs):
    return math.exp(sum(math.log(x) for x in xs) / len(xs))


def per_source_time(rows, key):
    """key(row) -> hashable configuration; returns {config: {source: median ns}}."""
    acc = defaultdict(lambda: defaultdict(list))
    for r in rows:
        assert r["ok"] == 1, "a run failed the correctness check: %r" % r
        acc[key(r)][r["source"]].append(r["ns"])
    return {c: {s: st.median(v) for s, v in d.items()} for c, d in acc.items()}


def ratio(tab, cfg, base):
    """Per-source ratios cfg/base and their summary."""
    # both configurations must have been run from exactly the same start vertices
    assert set(tab[cfg]) == set(tab[base]), "source sets differ: %r vs %r" % (cfg, base)
    rs = [tab[cfg][s] / tab[base][s] for s in tab[cfg]]
    return dict(gmean=gmean(rs), lo=min(rs), hi=max(rs), n=len(rs))


def per_source_count(rows, key, fields):
    acc = defaultdict(lambda: defaultdict(list))
    for r in rows:
        assert r["ok"] == 1
        for f in fields:
            acc[key(r)][f].append(r[f] / r["m"])
    return {c: {f: st.mean(v) for f, v in d.items()} for c, d in acc.items()}


def fmt(x, nd=2):
    return ("%." + str(nd) + "f") % x


def save(fig, name, guide=False):
    """Figures of the paper go to paper/figures; extra charts used only by the guide to explanation/figures."""
    fig.savefig(os.path.join(GUIDEFIG if guide else FIG, name + ".pdf"))
    plt.close(fig)
    print("  figure", name, "(guide only)" if guide else "")


FIGDATA = []  # (title, column names, rows): every plotted value, written to figure_data.txt


def figdata(title, header, rows):
    FIGDATA.append((title, header, rows))


# Plotting instructions written above each block of figure_data.txt (matched by the start of the block title).
HOWTO = {
    "FIGURE 1": [
        "Chart type: stacked bar chart (as in the paper, Fig. 1).",
        "X axis: road network, in the order OL, TG, CAL, SF, NA. Two bars per network, side by side:",
        "        left bar = bmssp25, right bar = bmssp26.",
        "Y axis: 'Edge scans per arc', linear scale, start at 0 (largest bar is about 5.1).",
        "Each bar is stacked from three columns, bottom to top: FindPivots, base_case, relaxation_loops",
        "        (their sum is the column 'total'). Use one colour per part and a legend with the three parts.",
        "Reference line: draw a horizontal line at y = 1 and label it 'Dijkstra' (Dijkstra scans every arc once).",
    ],
    "FIGURE 2": [
        "Chart type: heat map, one for bmssp25 and one for bmssp26, side by side (as in the paper, Fig. 2).",
        "X axis: t, the columns t=2 ... t=10.   Y axis: k, the rows k = 1 ... 6 (k = 1 at the bottom).",
        "Cell value: running time relative to Dijkstra (1.0 would mean 'as fast as Dijkstra').",
        "        Colour the cells by value, light = small = faster, dark = large = slower; use the same colour",
        "        scale for both heat maps and print the value (one decimal) inside each cell.",
        "Mark two cells per heat map: the prescribed parameters (given in the block title) and the smallest value.",
        "Alternative without a heat map: line chart with t on the x axis, the ratio on the y axis, one line per k.",
    ],
    "FIGURE 3": [
        "Chart type: line chart with markers, four panels in a row, one per value of 'family'",
        "        (rand = random graphs, grid = grids, ba = preferential attachment, hub = hub source).",
        "X axis: log2_n, i.e. the graph has n = 2^(log2_n) vertices; label it 'log2 n'. Linear scale.",
        "Y axis: 'Time relative to Dijkstra', linear scale, start at 0; use the same range (0 to about 9) in all panels.",
        "Lines: one per column dijkstra-dk, bmssp25, bmssp26 (three lines per panel).",
        "Reference line: horizontal line at y = 1 = Dijkstra with a binary heap (the baseline).",
        "Reading: a point at y = 2 means 'twice as slow as Dijkstra'; every bmssp25 and bmssp26 point is above 1.",
    ],
    "FIGURE 4": [
        "Chart type: line chart with markers, four panels in a row, one per value of 'family' (as for Figure 3).",
        "X axis: log2_n (n = 2^(log2_n) vertices), label 'log2 n'. Linear scale.",
        "Y axis: 'Comparisons per arc', linear scale, start at 0. Use one common range (0 to about 18) for",
        "        rand, grid and ba, and a separate range (0 to about 50) for hub.",
        "Lines: the four columns cmp_dijkstra, cmp_dijkstra-dk, cmp_bmssp25, cmp_bmssp26 (four lines per panel).",
        "The columns add_* (additions per arc) are NOT plotted in the paper; they are quoted in the text.",
        "        To plot 'comparisons + additions per arc' instead, add cmp_X + add_X for each algorithm X.",
    ],
    "EXTRA (table V": [
        "Not a figure in the paper. To plot: grouped bar chart, x = network, one bar per algorithm, y = ratio",
        "        (linear, from 0); the two last columns can be drawn as error bars; reference line at y = 1.",
    ],
    "EXTRA (table VI": [
        "Not a figure in the paper. To plot: two grouped bar charts, x = network, one bar per algorithm;",
        "        y = comparisons_per_arc in the first chart and additions_per_arc in the second (linear, from 0).",
    ],
    "EXTRA (table IX": [
        "Not a figure in the paper. To plot: grouped bar chart, x = graph, y = time_ratio (linear, from 0),",
        "        one bar per algorithm; make one chart for weights = unit and one for weights = random;",
        "        reference line at y = 1 = Dijkstra. 'comparisons_vs_dijkstra' below 1 means fewer comparisons.",
    ],
    "EXTRA (section VI-I of the paper, not a figure): random graphs with n = 2^18, time ratio by average": [
        "Not a figure in the paper. To plot: line chart, x = average_out_degree (logarithmic scale, base 2),",
        "        y = time ratio to Dijkstra (linear, from 0), one line per algorithm, reference line at y = 1.",
    ],
    "EXTRA (section VI-I of the paper, not a figure): random graphs with n = 2^18, time ratio by weight": [
        "Not a figure in the paper. To plot: line chart, x = log2_of_largest_weight (0 means all weights are 1),",
        "        y = time ratio to Dijkstra (linear, from 0), one line per algorithm, reference line at y = 1.",
    ],
}


def write_figdata():
    path = os.path.join(ROOT, "figure_data.txt")
    with open(path, "w") as f:
        f.write("# Values plotted in the figures of research-paper.pdf (tab-separated, one block per figure).\n"
                "# Generated by scripts/make_figures.py from results/*.csv; do not edit by hand.\n"
                "# 'ratio' = running time divided by the time of Dijkstra with a binary heap on the same graph and\n"
                "# start vertex, geometric mean over start vertices. 'per arc' = count divided by the number of arcs m.\n"
                "# Algorithms: dijkstra = binary heap with lazy deletion, dijkstra-dk = heap with DecreaseKey,\n"
                "# bmssp25 / bmssp26 = the 2025 / 2026 algorithm of Duan et al. (BMSSP-25 / BMSSP-26 in the paper).\n"
                "# Layout of a block: a title line starting with '##', plotting instructions in lines starting\n"
                "# with '#', then one header line with the column names, then the data rows.\n"
                "# In a spreadsheet: paste the header line and the rows; the columns separate at the tabs.\n")
        for title, header, rows in sorted(FIGDATA, key=lambda b: (not b[0].startswith("FIGURE"), b[0][:8] if b[0].startswith("FIGURE") else "")):
            f.write("\n## %s\n" % title)
            how = [v for k, v in HOWTO.items() if title.startswith(k)]
            assert how, "no plotting instructions for block: " + title
            f.write("# HOW TO PLOT\n")
            for line in how[0]:
                f.write("#   " + line + "\n")
            f.write("\t".join(header) + "\n")
            for r in rows:
                f.write("\t".join(("%.4g" % x) if isinstance(x, float) else str(x) for x in r) + "\n")
    print("  wrote figure_data.txt (%d blocks)" % len(FIGDATA))


def write(name, text):
    with open(os.path.join(TAB, name + ".tex"), "w") as f:
        f.write(text)
    print("  table", name)


ROADPTS = {}  # per road network: frontier width and ratios, reused by the frontier figure
NUM = {}      # macro name -> value (string)
SUMMARY = {}  # json summary


def macro(name, value):
    assert name.isalpha(), name
    NUM[name] = str(value)


# ---------------------------------------------------------------- road networks
def road_tables_and_figures():
    info = {r["graph"]: r for r in read("graph_info")}
    T = per_source_time(read("real_time"), lambda r: (r["graph"], r["algo"]))
    cnt_rows = read("real_count")
    fields = ["cmp", "add", "scan", "relaxed", "scan_pivot", "scan_base", "scan_main", "ds_insert", "ds_batch",
              "ds_pulled", "heap_push", "sel_elems", "pivots", "frontier", "stale_erase"]
    Cn = per_source_count(cnt_rows, lambda r: (r["graph"], r["algo"]), fields)
    par = {}
    for r in cnt_rows:
        par[(r["graph"], r["algo"])] = (r["k"], r["t"], r["levels"])
    peak = defaultdict(list)
    for r in cnt_rows:
        if r["algo"] == "dijkstra-dk":
            peak[r["graph"]].append(r["peak"] / r["n"])

    # Table: data sets
    lines = []
    for g in ROADS:
        i = info[g]
        k25, t25, l25 = par[(g, "bmssp25")]
        k26, t26, l26 = par[(g, "bmssp26")]
        lines.append("%s & %s & %s & %d & %.2f & %.1f & %d/%d/%d & %d/%d/%d \\\\" % (
            g, "{:,}".format(i["n"]).replace(",", "\\,"), "{:,}".format(i["m"]).replace(",", "\\,"),
            i["max_out"], i["m"] / i["n"], 100 * st.mean(peak[g]), k25, t25, l25, k26, t26, l26))
    write("datasets", "\n".join(lines) + "\n")

    # Table: running times
    lines = []
    summ = {}
    for g in ROADS:
        base = (g, "dijkstra")
        med = {a: st.median(T[(g, a)].values()) / 1e6 for a in STYLE}
        rr = {a: ratio(T, (g, a), base) for a in STYLE}
        summ[g] = {"median_ms": med, "ratio": rr}
        lines.append("%s & %s & %s & %s & %s & %s & %s & %s \\\\" % (
            g, fmt(med["dijkstra"], 3), fmt(med["dijkstra-dk"], 3), fmt(rr["dijkstra-dk"]["gmean"]),
            fmt(med["bmssp25"], 3), fmt(rr["bmssp25"]["gmean"]), fmt(med["bmssp26"], 3), fmt(rr["bmssp26"]["gmean"])))
    write("real_time", "\n".join(lines) + "\n")
    SUMMARY["real_time"] = summ
    r25 = [summ[g]["ratio"]["bmssp25"]["gmean"] for g in ROADS]
    r26 = [summ[g]["ratio"]["bmssp26"]["gmean"] for g in ROADS]
    macro("RealRatioMinTwoFive", fmt(min(r25), 1)); macro("RealRatioMaxTwoFive", fmt(max(r25), 1))
    macro("RealRatioMinTwoSix", fmt(min(r26), 1)); macro("RealRatioMaxTwoSix", fmt(max(r26), 1))
    macro("RealRatioMin", fmt(min(r25 + r26), 1)); macro("RealRatioMax", fmt(max(r25 + r26), 1))
    macro("RealSources", len(next(iter(T.values()))))
    # share of degree-two vertices (chains) and the slowdown of BMSSP-25 go together on the road networks
    by_deg = sorted(ROADS, key=lambda g: -info[g]["deg2_pct"])
    by_ratio = sorted(ROADS, key=lambda g: -summ[g]["ratio"]["bmssp25"]["gmean"])
    assert by_deg[:2] == by_ratio[:2] or set(by_deg[:2]) == set(by_ratio[:2])
    top = sorted(by_deg[:2], key=ROADS.index)
    macro("ChainNets", " and ".join(top))
    macro("ChainPctA", fmt(info[top[0]]["deg2_pct"], 1)); macro("ChainPctB", fmt(info[top[1]]["deg2_pct"], 1))
    rest = [g for g in ROADS if g not in top]
    macro("ChainRestMin", fmt(min(info[g]["deg2_pct"] for g in rest), 0)); macro("ChainRestMax", fmt(max(info[g]["deg2_pct"] for g in rest), 0))
    macro("ChainRatioMin", fmt(min(summ[g]["ratio"]["bmssp25"]["gmean"] for g in top), 1))
    macro("ChainRatioMax", fmt(max(summ[g]["ratio"]["bmssp25"]["gmean"] for g in top), 1))
    macro("RestRatioMin", fmt(min(summ[g]["ratio"]["bmssp25"]["gmean"] for g in rest), 1))
    macro("RestRatioMax", fmt(max(summ[g]["ratio"]["bmssp25"]["gmean"] for g in rest), 1))
    faster26 = [g for g in ROADS if summ[g]["ratio"]["bmssp26"]["gmean"] < summ[g]["ratio"]["bmssp25"]["gmean"]]
    macro("RealTwoSixFasterCount", len(faster26))
    macro("RealTwoSixFasterWord", ["none", "one", "two", "three", "four", "five"][len(faster26)])
    macro("RealTwoSixSlowerWord", ["none", "one", "two", "three", "four", "five"][len(ROADS) - len(faster26)])
    macro("RealTwoSixFasterList", ", ".join(faster26))
    SUMMARY["real_26_faster_on"] = faster26

    # Table: operation counts per arc
    lines = []
    for g in ROADS:
        assert abs(Cn[(g, "dijkstra")]["add"] - 1) < 1e-9 and abs(Cn[(g, "dijkstra")]["scan"] - 1) < 1e-9
        row = [g, fmt(Cn[(g, "dijkstra")]["cmp"])]
        for a in ("bmssp25", "bmssp26"):
            for f in ("cmp", "add", "scan"):
                row.append(fmt(Cn[(g, a)][f]))
        lines.append(" & ".join(row) + " \\\\")
    write("real_count", "\n".join(lines) + "\n")
    SUMMARY["real_count"] = {"%s/%s" % k: v for k, v in Cn.items()}
    for tag, a in (("TwoFive", "bmssp25"), ("TwoSix", "bmssp26")):
        cr = [Cn[(g, a)]["cmp"] / Cn[(g, "dijkstra")]["cmp"] for g in ROADS]
        ar = [Cn[(g, a)]["add"] / Cn[(g, "dijkstra")]["add"] for g in ROADS]
        sr = [Cn[(g, a)]["scan"] for g in ROADS]
        macro("CmpRatioMin" + tag, fmt(min(cr))); macro("CmpRatioMax" + tag, fmt(max(cr)))
        macro("AddRatioMin" + tag, fmt(min(ar), 1)); macro("AddRatioMax" + tag, fmt(max(ar), 1))
        macro("ScanMin" + tag, fmt(min(sr), 1)); macro("ScanMax" + tag, fmt(max(sr), 1))
    for tag, a in (("TwoFive", "bmssp25"), ("TwoSix", "bmssp26")):
        pv = [Cn[(g, a)]["scan_pivot"] for g in ROADS]
        macro("PivotScanMin" + tag, fmt(min(pv), 2)); macro("PivotScanMax" + tag, fmt(max(pv), 2))
        mn = [Cn[(g, a)]["scan_main"] for g in ROADS]
        macro("MainScanMin" + tag, fmt(min(mn), 1)); macro("MainScanMax" + tag, fmt(max(mn), 1))
    # cost per comparison-or-addition, relative to Dijkstra (time ratio / operation ratio)
    for tag, a in (("TwoFive", "bmssp25"), ("TwoSix", "bmssp26")):
        v = []
        for g in ROADS:
            ops = (Cn[(g, a)]["cmp"] + Cn[(g, a)]["add"]) / (Cn[(g, "dijkstra")]["cmp"] + Cn[(g, "dijkstra")]["add"])
            v.append(summ[g]["ratio"][a]["gmean"] / ops)
        macro("PerOpMin" + tag, fmt(min(v), 1)); macro("PerOpMax" + tag, fmt(max(v), 1))
        o = [(Cn[(g, a)]["cmp"] + Cn[(g, a)]["add"]) / (Cn[(g, "dijkstra")]["cmp"] + Cn[(g, "dijkstra")]["add"])
             for g in ROADS]
        macro("OpsRatioMin" + tag, fmt(min(o), 1)); macro("OpsRatioMax" + tag, fmt(max(o), 1))
        SUMMARY["per_op_cost_" + a] = dict(zip(ROADS, v))
    words = ["none", "one", "two", "three", "four", "five"]
    macro("RoadsMoreAddsTwoSix", words[sum(1 for g in ROADS if Cn[(g, "bmssp26")]["add"] > Cn[(g, "bmssp25")]["add"])])
    macro("RoadsMoreCmpTwoSix", words[sum(1 for g in ROADS if Cn[(g, "bmssp26")]["cmp"] > Cn[(g, "bmssp25")]["cmp"])])
    SUMMARY["frontier_peak_over_n"] = {g: st.mean(v) for g, v in peak.items()}
    macro("RoadPeakMin", fmt(100 * min(st.mean(v) for v in peak.values()), 1))
    macro("RoadPeakMax", fmt(100 * max(st.mean(v) for v in peak.values()), 1))
    ROADPTS.update({g: dict(peak=st.mean(peak[g]),
                            cmp={a: Cn[(g, a)]["cmp"] / Cn[(g, "dijkstra")]["cmp"] for a in ("bmssp25", "bmssp26")},
                            ops={a: (Cn[(g, a)]["cmp"] + Cn[(g, a)]["add"]) / (Cn[(g, "dijkstra")]["cmp"] + 1.0)
                                 for a in ("bmssp25", "bmssp26")},
                            time={a: summ[g]["ratio"][a]["gmean"] for a in ("bmssp25", "bmssp26")}) for g in ROADS})

    # Figure: running time relative to Dijkstra
    fig, ax = plt.subplots(figsize=(COL1, 2.05))
    algs = ["dijkstra-dk", "bmssp25", "bmssp26"]
    w = 0.26
    for j, a in enumerate(algs):
        xs = np.arange(len(ROADS)) + (j - 1) * w
        ys = [summ[g]["ratio"][a]["gmean"] for g in ROADS]
        lo = [summ[g]["ratio"][a]["gmean"] - summ[g]["ratio"][a]["lo"] for g in ROADS]
        hi = [summ[g]["ratio"][a]["hi"] - summ[g]["ratio"][a]["gmean"] for g in ROADS]
        ax.bar(xs, ys, w * 0.92, color=STYLE[a]["color"], hatch=STYLE[a]["hatch"], edgecolor="white", linewidth=0.4,
               label=STYLE[a]["label"], zorder=3)
        ax.errorbar(xs, ys, yerr=[lo, hi], fmt="none", ecolor=INK, elinewidth=0.6, capsize=1.5, zorder=4)
        for x, y, h in zip(xs, ys, hi):
            ax.text(x, y + h + 0.12, fmt(y, 1), ha="center", va="bottom", fontsize=6, color=INK)
    ax.axhline(1, color=STYLE["dijkstra"]["color"], lw=1.2, zorder=2, label="Dijkstra (binary heap) = 1")
    ax.set_xticks(range(len(ROADS)))
    ax.set_xticklabels(["%s\n$n$=%s" % (g, "{:,}".format(info[g]["n"])) for g in ROADS])
    ax.set_ylabel("Running time relative to Dijkstra")
    ax.set_ylim(0, max(summ[g]["ratio"][a]["hi"] for g in ROADS for a in algs) * 1.22)
    ax.yaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    h, l = ax.get_legend_handles_labels()
    order = [1, 2, 3, 0]
    ax.legend([h[i] for i in order], [l[i] for i in order], loc="upper left", ncol=2, frameon=False,
              handlelength=1.4, columnspacing=0.9, bbox_to_anchor=(-0.02, 1.27))
    save(fig, "real_ratio", guide=True)
    figdata("EXTRA (table V of the paper, not a figure): road networks, time ratio to Dijkstra",
            ["network", "algorithm", "ratio", "smallest_per_start_vertex", "largest_per_start_vertex"],
            [(g, a, summ[g]["ratio"][a]["gmean"], summ[g]["ratio"][a]["lo"], summ[g]["ratio"][a]["hi"])
             for g in ROADS for a in algs])
    for tag, a in (("TwoFive", "bmssp25"), ("TwoSix", "bmssp26")):
        macro("RealSrcMin" + tag, fmt(min(summ[g]["ratio"][a]["lo"] for g in ROADS), 1))
        macro("RealSrcMax" + tag, fmt(max(summ[g]["ratio"][a]["hi"] for g in ROADS), 1))

    # Figure: operation counts and where the edge scans happen
    fig, axs = plt.subplots(1, 3, figsize=(COL2, 1.95), gridspec_kw=dict(wspace=0.3))
    for ax, f, title in ((axs[0], "cmp", "Comparisons per arc"), (axs[1], "add", "Additions per arc")):
        for j, a in enumerate(["dijkstra", "bmssp25", "bmssp26"]):
            xs = np.arange(len(ROADS)) + (j - 1) * w
            ys = [Cn[(g, a)][f] for g in ROADS]
            ax.bar(xs, ys, w * 0.92, color=STYLE[a]["color"], hatch=STYLE[a]["hatch"], edgecolor="white",
                   linewidth=0.4, label=STYLE[a]["label"] if a != "dijkstra" else "Dijkstra", zorder=3)
            for x, y in zip(xs, ys):
                ax.text(x, y + 0.01 * ax.get_ylim()[1], fmt(y, 1), ha="center", va="bottom", fontsize=5.5, rotation=90)
        ax.set_xticks(range(len(ROADS)))
        ax.set_xticklabels(ROADS)
        ax.set_ylabel(title)
        ax.yaxis.grid(True, zorder=0)
        ax.set_axisbelow(True)
        ax.set_ylim(0, ax.get_ylim()[1] * 1.18)
    axs[0].legend(loc="upper left", ncol=3, frameon=False, handlelength=1.4, columnspacing=0.8,
                  bbox_to_anchor=(-0.05, 1.2))
    ax = axs[2]
    parts = [("scan_pivot", "FindPivots", "#b9d4f2"), ("scan_base", "base case", "#5b9ae0"),
             ("scan_main", "relaxation loops", "#1c4f91")]
    for j, a in enumerate(["bmssp25", "bmssp26"]):
        xs = np.arange(len(ROADS)) + (j - 0.5) * 0.38
        bottom = np.zeros(len(ROADS))
        for f, lab, col in parts:
            ys = np.array([Cn[(g, a)][f] for g in ROADS])
            ax.bar(xs, ys, 0.34, bottom=bottom, color=col, edgecolor="white", linewidth=0.6,
                   label=lab if j == 0 else None, zorder=3)
            bottom += ys
        for x, y in zip(xs, bottom):
            ax.text(x, y + 0.08, "B-25" if a == "bmssp25" else "B-26", ha="center", va="bottom", fontsize=5, color=MUTED)
    ax.set_xticks(range(len(ROADS)))
    ax.set_xticklabels(ROADS)
    ax.set_ylabel("Edge scans per arc")
    ax.yaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.12)
    ax.legend(loc="upper left", ncol=3, frameon=False, handlelength=1.0, columnspacing=0.7,
              bbox_to_anchor=(-0.12, 1.2))
    save(fig, "real_counts", guide=True)
    figdata("EXTRA (table VI of the paper, not a figure): road networks, comparisons and additions per arc",
            ["network", "algorithm", "comparisons_per_arc", "additions_per_arc"],
            [(g, a, Cn[(g, a)]["cmp"], Cn[(g, a)]["add"]) for g in ROADS for a in ("dijkstra", "bmssp25", "bmssp26")])

    # Figure of the paper: where the edge scans happen (single panel)
    fig, ax = plt.subplots(figsize=(COL1, 1.9))
    for j, a in enumerate(["bmssp25", "bmssp26"]):
        xs = np.arange(len(ROADS)) + (j - 0.5) * 0.38
        bottom = np.zeros(len(ROADS))
        for f, lab, col in parts:
            ys = np.array([Cn[(g, a)][f] for g in ROADS])
            ax.bar(xs, ys, 0.34, bottom=bottom, color=col, edgecolor="white", linewidth=0.6,
                   label=lab if j == 0 else None, zorder=3)
            bottom += ys
        for x, y in zip(xs, bottom):
            ax.text(x, y + 0.08, "B-25" if a == "bmssp25" else "B-26", ha="center", va="bottom", fontsize=5.5, color=MUTED)
    ax.axhline(1, color=STYLE["dijkstra"]["color"], lw=1.2, zorder=4, label="Dijkstra")
    ax.set_xticks(range(len(ROADS)))
    ax.set_xticklabels(ROADS)
    ax.set_ylabel("Edge scans per arc")
    ax.yaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.1)
    ax.legend(loc="upper center", ncol=4, frameon=False, handlelength=1.0, columnspacing=0.8,
              bbox_to_anchor=(0.5, 1.2))
    save(fig, "edge_scans")
    figdata("FIGURE 1: road networks, edge scans per arc by part of the algorithm (Dijkstra = 1.0 everywhere)",
            ["network", "algorithm", "FindPivots", "base_case", "relaxation_loops", "total"],
            [(g, a, Cn[(g, a)]["scan_pivot"], Cn[(g, a)]["scan_base"], Cn[(g, a)]["scan_main"], Cn[(g, a)]["scan"])
             for g in ROADS for a in ("bmssp25", "bmssp26")])


# ------------------------------------------------------------ variants / defects
def variants_tables():
    rows = read("variants_count")
    acc = defaultdict(lambda: defaultdict(list))
    for r in rows:
        assert r["ok"] == 1
        key = (r["graph"], r["algo"], r["variant"])
        for f in ("stale_erase", "stale_pull", "dup_u", "bound_viol", "full_nonempty", "scan", "cmp"):
            acc[key][f].append(r[f])
        acc[key]["scan_m"].append(r["scan"] / r["m"])
        acc[key]["cmp_m"].append(r["cmp"] / r["m"])
        acc[key]["n"].append(r["n"])
    T = per_source_time(read("variants_time"), lambda r: (r["graph"], r["algo"], r["variant"]))
    lines = []
    out = {}
    order = [("bmssp25", "as-printed", "BMSSP-25 as printed"), ("bmssp25", "repaired", "BMSSP-25 with deletion"),
             ("bmssp26", "as-printed", "BMSSP-26 as printed"), ("bmssp26", "no-deletion", "\\quad base case repaired only"),
             ("bmssp26", "printed-base-case", "\\quad deletion only"), ("bmssp26", "repaired", "BMSSP-26 with both repairs")]
    for g in ("SF", "NA"):
        for a, v, label in order:
            d = acc[(g, a, v)]
            ref = (g, a, "repaired")
            tr = ratio(T, (g, a, v), ref)["gmean"]
            out["%s/%s/%s" % (g, a, v)] = dict(stale=st.mean(d["stale_erase"]), dup=st.mean(d["dup_u"]),
                                               bound=st.mean(d["bound_viol"]), scan=st.mean(d["scan_m"]),
                                               cmp=st.mean(d["cmp_m"]), time_vs_repaired=tr)
            lines.append("%s & %s & %s & %s & %s & %s & %s & %s \\\\" % (
                g if (a, v) == order[0][:2] else "", label, "{:,.0f}".format(st.mean(d["stale_erase"])).replace(",", "\\,"),
                "{:,.0f}".format(st.mean(d["dup_u"])).replace(",", "\\,"),
                "{:,.0f}".format(st.mean(d["bound_viol"])).replace(",", "\\,"),
                fmt(st.mean(d["scan_m"])), fmt(st.mean(d["cmp_m"])), fmt(tr)))
        if g == "SF":
            lines.append("\\hline")
    write("variants", "\n".join(lines) + "\n")
    SUMMARY["variants"] = out
    allv = {}
    for g in ROADS:
        for a, v, _ in order:
            d = acc[(g, a, v)]
            allv["%s/%s/%s" % (g, a, v)] = dict(stale=st.mean(d["stale_erase"]), dup=st.mean(d["dup_u"]),
                                                bound=st.mean(d["bound_viol"]),
                                                dup_frac_n=st.mean(d["dup_u"]) / d["n"][0],
                                                time_vs_repaired=ratio(T, (g, a, v), (g, a, "repaired"))["gmean"])
    SUMMARY["variants_all"] = allv
    d25 = [allv["%s/bmssp25/as-printed" % g]["dup_frac_n"] * 100 for g in ROADS]
    d26 = [allv["%s/bmssp26/as-printed" % g]["dup_frac_n"] * 100 for g in ROADS]
    macro("DupPctMinTwoFive", fmt(min(d25), 1)); macro("DupPctMaxTwoFive", fmt(max(d25), 1))
    macro("DupPctMinTwoSix", fmt(min(d26), 1)); macro("DupPctMaxTwoSix", fmt(max(d26), 1))
    t25 = [allv["%s/bmssp25/as-printed" % g]["time_vs_repaired"] for g in ROADS]
    t26 = [allv["%s/bmssp26/as-printed" % g]["time_vs_repaired"] for g in ROADS]
    macro("PrintedTimeMinTwoFive", fmt(min(t25))); macro("PrintedTimeMaxTwoFive", fmt(max(t25)))
    macro("PrintedTimeMinTwoSix", fmt(min(t26))); macro("PrintedTimeMaxTwoSix", fmt(max(t26)))

    # validation summary
    for name in ("small_n16", "medium_n400"):
        path = os.path.join(RES, "validation", name + ".csv")
        if not os.path.exists(path):
            continue
        with open(path) as f:
            v = {r["variant"]: {k: int(x) for k, x in r.items() if k != "variant"} for r in csv.DictReader(f)}
        SUMMARY["validation_" + name] = v
        tag = "Small" if name == "small_n16" else "Medium"
        macro("Val%sTrials" % tag, "{:,}".format(v["dijkstra"]["runs"]).replace(",", "\\,"))
        macro("Val%sRuns" % tag, "{:,}".format(sum(x["runs"] for x in v.values())).replace(",", "\\,"))
        macro("Val%sWrong" % tag, sum(x["wrong"] for x in v.values()))
        macro("Val%sDupTwoFive" % tag, "{:,}".format(v["bmssp25/literal/paper-params"]["runs_with_duplicates"]).replace(",", "\\,"))
        macro("Val%sDupTwoSix" % tag, "{:,}".format(v["bmssp26/literal/paper-params"]["runs_with_duplicates"]).replace(",", "\\,"))
        macro("Val%sBoundTwoSix" % tag, "{:,}".format(v["bmssp26/literal/paper-params"]["runs_with_bound_violation"]).replace(",", "\\,"))
        macro("Val%sStaleTwoSix" % tag, "{:,}".format(v["bmssp26/fix+basefix/paper-params"]["runs_with_stale_keys"]).replace(",", "\\,"))
        macro("Val%sStaleTwoFive" % tag, "{:,}".format(v["bmssp25/fix/paper-params"]["runs_with_stale_keys"]).replace(",", "\\,"))
        if name == "small_n16":
            macro("ValMinNDupTwoFive", v["bmssp25/literal/paper-params"]["min_n_dup"])
            macro("ValMinNDupTwoSix", v["bmssp26/literal/paper-params"]["min_n_dup"])
            macro("ValMinNBoundTwoSix", v["bmssp26/literal/paper-params"]["min_n_bound"])


# -------------------------------------------------------------- degree reduction
def reduce_table():
    rows = read("reduce_time")
    if not rows:
        return
    T = per_source_time(rows, lambda r: (r["graph"], r["algo"], r["variant"]))
    size = {(r["graph"], r["variant"]): (r["n"], r["m"]) for r in rows}
    lines = []
    out = {}
    for g in ROADS:
        base = (g, "dijkstra", "original")
        n3, m3 = size[(g, "delta3")]
        n0, m0 = size[(g, "original")]
        vals = [g, fmt(n3 / n0), fmt(m3 / m0)]
        for a, v in (("dijkstra", "delta3"), ("bmssp25", "original"), ("bmssp25", "delta3"), ("bmssp26", "original"),
                     ("bmssp26", "delta3"), ("bmssp26", "delta4")):
            r = ratio(T, (g, a, v), base)["gmean"]
            out["%s/%s/%s" % (g, a, v)] = r
            vals.append(fmt(r, 1))
        lines.append(" & ".join(vals) + " \\\\")
    write("reduce", "\n".join(lines) + "\n")
    SUMMARY["reduce"] = out
    s25 = [out["%s/bmssp25/delta3" % g] / out["%s/bmssp25/original" % g] for g in ROADS]
    s26 = [out["%s/bmssp26/delta3" % g] / out["%s/bmssp26/original" % g] for g in ROADS]
    macro("ReduceSlowMinTwoFive", fmt(min(s25), 1)); macro("ReduceSlowMaxTwoFive", fmt(max(s25), 1))
    macro("ReduceSlowMinTwoSix", fmt(min(s26), 1)); macro("ReduceSlowMaxTwoSix", fmt(max(s26), 1))
    macro("ReduceNGrowMin", fmt(min(size[(g, "delta3")][0] / size[(g, "original")][0] for g in ROADS), 1))
    macro("ReduceNGrowMax", fmt(max(size[(g, "delta3")][0] / size[(g, "original")][0] for g in ROADS), 1))


# --------------------------------------------------------------- parameter sweep
def params_figure():
    rows = read("params_time")
    if not rows:
        return
    T = per_source_time(rows, lambda r: (r["graph"], r["algo"], r["k"], r["t"]))
    info = {(r["graph"], r["algo"], r["k"], r["t"]): r["levels"] for r in rows}
    defaults = {}
    for r in read("real_count"):
        defaults[(r["graph"], r["algo"])] = (r["k"], r["t"])
    best = {}
    for g in ("SF", "NA"):
        for a in ("bmssp25", "bmssp26"):
            cells = {}
            for (gg, aa, k, t) in T:
                if gg == g and aa == a:
                    cells[(k, t)] = ratio(T, (g, a, k, t), (g, "dijkstra", 0, 0))["gmean"]
            kb, tb = min(cells, key=cells.get)
            dk, dt = defaults[(g, a)]
            best[(g, a)] = dict(best_k=kb, best_t=tb, best_ratio=cells[(kb, tb)], default_k=dk, default_t=dt,
                                default_ratio=cells[(dk, dt)], worst_ratio=max(cells.values()),
                                levels_best=info[(g, a, kb, tb)], cells={"%d,%d" % c: v for c, v in cells.items()})
    SUMMARY["params"] = {"%s/%s" % k: v for k, v in best.items()}
    # the text says "on both networks": make sure that stays true if the data are regenerated
    assert best[("SF", "bmssp26")]["best_t"] == best[("NA", "bmssp26")]["best_t"]
    for g in ("SF", "NA"):
        for tag, a in (("TwoFive", "bmssp25"), ("TwoSix", "bmssp26")):
            b = best[(g, a)]
            macro("Par%s%sBestK" % (g, tag), b["best_k"]); macro("Par%s%sBestT" % (g, tag), b["best_t"])
            macro("Par%s%sBest" % (g, tag), fmt(b["best_ratio"])); macro("Par%s%sDefault" % (g, tag), fmt(b["default_ratio"]))
            macro("Par%s%sWorst" % (g, tag), fmt(b["worst_ratio"], 1))
            macro("Par%s%sGain" % (g, tag), "%.0f" % (100 * (1 - b["best_ratio"] / b["default_ratio"])))
    fig, axs = plt.subplots(1, 2, figsize=(COL2, 2.25), gridspec_kw=dict(wspace=0.12))
    g = "SF"
    ks = list(range(1, 7))
    ts = list(range(2, 11))
    vmax = max(max(float(v) for v in best[(g, a)]["cells"].values()) for a in ("bmssp25", "bmssp26"))
    vmin = min(min(float(v) for v in best[(g, a)]["cells"].values()) for a in ("bmssp25", "bmssp26"))
    for ax, a in zip(axs, ("bmssp25", "bmssp26")):
        M = np.array([[best[(g, a)]["cells"]["%d,%d" % (k, t)] for t in ts] for k in ks])
        im = ax.imshow(M, cmap="Blues", vmin=vmin * 0.9, vmax=vmax, aspect="auto", origin="lower")
        for i, k in enumerate(ks):
            for j, t in enumerate(ts):
                v = M[i, j]
                ax.text(j, i, fmt(v, 1), ha="center", va="center", fontsize=6,
                        color="white" if (v - vmin * 0.9) / (vmax - vmin * 0.9) > 0.55 else INK)
        b = best[(g, a)]
        for (k, t, style) in ((b["default_k"], b["default_t"], "-"), (b["best_k"], b["best_t"], ":")):
            ax.add_patch(plt.Rectangle((ts.index(t) - 0.5, ks.index(k) - 0.5), 1, 1, fill=False, ec="#eb6834", lw=1.4,
                                       ls=style, zorder=5))
        ax.set_xticks(range(len(ts))); ax.set_xticklabels(ts)
        ax.set_yticks(range(len(ks))); ax.set_yticklabels(ks)
        ax.set_xlabel("$t$")
        ax.set_title("%s on SF" % STYLE[a]["label"], fontsize=8)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0)
    axs[0].set_ylabel("$k$")
    cb = fig.colorbar(im, ax=axs, fraction=0.025, pad=0.015)
    cb.set_label("Time relative to Dijkstra", fontsize=7)
    cb.outline.set_visible(False)
    save(fig, "params_heatmap")
    for gg in ("SF", "NA"):
        figdata("FIGURE 2%s: parameter sweep on %s, time ratio to Dijkstra for every k and t "
                "(prescribed: bmssp25 k=%d t=%d, bmssp26 k=%d t=%d)" % (
                    "" if gg == "SF" else " (same sweep on NA, quoted in the text only)", gg,
                    best[(gg, "bmssp25")]["default_k"], best[(gg, "bmssp25")]["default_t"],
                    best[(gg, "bmssp26")]["default_k"], best[(gg, "bmssp26")]["default_t"]),
                ["algorithm", "k"] + ["t=%d" % t for t in ts],
                [(a, k) + tuple(best[(gg, a)]["cells"]["%d,%d" % (k, t)] for t in ts)
                 for a in ("bmssp25", "bmssp26") for k in ks])


# --------------------------------------------------------------------- scaling
def scaling_figures():
    rows = read("scaling_time")
    if not rows:
        return
    T = per_source_time(rows, lambda r: (r["family"], r["log2n"], r["algo"]))
    fams = [("rand", "Random, out-degree 3"), ("grid", "Grid, random weights"), ("ba", "Preferential attachment"),
            ("hub", "Hub source")]
    fig, axs = plt.subplots(1, 4, figsize=(COL2, 1.85), sharey=True, gridspec_kw=dict(wspace=0.08))
    out = {}
    for ax, (fam, title) in zip(axs, fams):
        es = sorted({e for (f, e, a) in T if f == fam})
        for a in ("dijkstra-dk", "bmssp25", "bmssp26"):
            ys = []
            for e in es:
                r = ratio(T, (fam, e, a), (fam, e, "dijkstra"))
                ys.append(r["gmean"])
                out["%s/%d/%s" % (fam, e, a)] = r["gmean"]
            ax.plot(es, ys, color=STYLE[a]["color"], marker=STYLE[a]["marker"], ls=STYLE[a]["ls"],
                    label=STYLE[a]["label"], markeredgecolor="white", markeredgewidth=0.4, zorder=3)
        ax.axhline(1, color=STYLE["dijkstra"]["color"], lw=1.2, zorder=2,
                   label="Dijkstra (binary heap) = 1" if fam == "rand" else None)
        ax.set_title(title, fontsize=7.5)
        ax.set_xlabel("$\\log_2 n$")
        ax.set_xticks([e for e in es if e % 2 == 0])
        ax.yaxis.grid(True, zorder=0)
        ax.set_axisbelow(True)
    axs[0].set_ylabel("Time relative to Dijkstra")
    axs[0].set_ylim(0, None)
    h, l = axs[0].get_legend_handles_labels()
    order = [3, 0, 1, 2]
    fig.legend([h[i] for i in order], [l[i] for i in order], loc="upper center", ncol=4, frameon=False,
               handlelength=2.0, columnspacing=1.6, bbox_to_anchor=(0.5, 1.09))
    save(fig, "scaling_ratio")
    SUMMARY["scaling_ratio"] = out
    figdata("FIGURE 3: synthetic families, time ratio to Dijkstra as n grows",
            ["family", "log2_n", "dijkstra-dk", "bmssp25", "bmssp26"],
            [(fam, e) + tuple(out["%s/%d/%s" % (fam, e, a)] for a in ("dijkstra-dk", "bmssp25", "bmssp26"))
             for fam, _ in fams for e in sorted({e for (f, e, a) in T if f == fam})])
    pts = sorted({(f, e) for (f, e, a) in T})
    macro("SynthPoints", len(pts))
    allr = [out["%s/%d/%s" % (f, e, a)] for f, e in pts for a in ("bmssp25", "bmssp26")]
    allr += [d["time"][a] for d in ROADPTS.values() for a in ("bmssp25", "bmssp26")]
    macro("AllRatioMin", fmt(min(allr), 1)); macro("AllRatioMax", fmt(max(allr), 1))
    macro("SynthTwoSixFasterCount", sum(1 for f, e in pts if out["%s/%d/bmssp26" % (f, e)] < out["%s/%d/bmssp25" % (f, e)]))
    macro("TimeRandTwentyTwoFive", fmt(out["rand/20/bmssp25"], 1)); macro("TimeRandTwentyTwoSix", fmt(out["rand/20/bmssp26"], 1))
    # absolute times on random graphs, for the comparison with published numbers
    for e, word in ((16, "Sixteen"), (20, "Twenty"), (22, "TwentyTwo")):
        for tag, a in (("Dij", "dijkstra"), ("TwoFive", "bmssp25")):
            macro("Ms%s%s" % (word, tag), "%.0f" % (st.median(T[("rand", e, a)].values()) / 1e6))
        macro("RatioRand%sTwoFive" % word, fmt(out["rand/%d/bmssp25" % e], 1))
    late = [out["grid/%d/%s" % (e, a)] for f, e in pts if f == "grid" and e >= 14 for a in ("bmssp25", "bmssp26")]
    macro("GridLateMin", fmt(min(late), 1)); macro("GridLateMax", fmt(max(late), 1))
    dk = [out["%s/%d/dijkstra-dk" % (f, e)] for f, e in pts if e >= 21]
    macro("DkLargeMax", fmt(max(dk), 1))
    for fam, _ in fams:
        es = sorted({e for (f, e, a) in T if f == fam})
        for tag, a in (("TwoFive", "bmssp25"), ("TwoSix", "bmssp26")):
            vals = [out["%s/%d/%s" % (fam, e, a)] for e in es]
            name = {"rand": "Rand", "grid": "Grid", "ba": "Ba", "hub": "Hub"}[fam]
            macro("Scale%s%sFirst" % (name, tag), fmt(vals[0], 1)); macro("Scale%s%sLast" % (name, tag), fmt(vals[-1], 1))
            macro("Scale%s%sMin" % (name, tag), fmt(min(vals), 2)); macro("Scale%s%sMax" % (name, tag), fmt(max(vals), 1))
            macro("Scale%s%sArgmin" % (name, tag), es[vals.index(min(vals))])
        macro("Scale%sMaxExp" % {"rand": "Rand", "grid": "Grid", "ba": "Ba", "hub": "Hub"}[fam], es[-1])

    # machine-independent growth: comparisons and comparisons + additions per arc
    rows = read("scaling_count")
    if not rows:
        return
    Cn = per_source_count(rows, lambda r: (r["family"], r["log2n"], r["algo"]), ["cmp", "add", "scan"])
    pk = defaultdict(list)
    for r in rows:
        if r["algo"] == "dijkstra-dk":
            pk[(r["family"], r["log2n"])].append(r["peak"] / r["n"])
    pk = {k: st.mean(v) for k, v in pk.items()}
    fig, axs = plt.subplots(1, 4, figsize=(COL2, 1.9), gridspec_kw=dict(wspace=0.2))
    cout = {}
    algs = ("dijkstra", "dijkstra-dk", "bmssp25", "bmssp26")
    for col, (fam, title) in enumerate(fams):
        es = sorted({e for (f, e, a) in Cn if f == fam})
        for a in algs:
            for e in es:
                c = Cn[(fam, e, a)]
                cout["%s/%d/%s" % (fam, e, a)] = dict(cmp=c["cmp"], add=c["add"], scan=c["scan"])
            axs[col].plot(es, [Cn[(fam, e, a)]["cmp"] for e in es], color=STYLE[a]["color"],
                          marker=STYLE[a]["marker"], ls=STYLE[a]["ls"], label=STYLE[a]["label"],
                          markeredgecolor="white", markeredgewidth=0.4, zorder=3)
        axs[col].set_title(title, fontsize=7.5)
        axs[col].set_xlabel("$\\log_2 n$")
        axs[col].set_xticks([e for e in es if e % 2 == 0])
        axs[col].yaxis.grid(True, zorder=0)
        axs[col].set_axisbelow(True)
    axs[0].set_ylabel("Comparisons per arc")
    top = max(axs[c].get_ylim()[1] for c in range(3))  # common scale for the first three families
    for c in range(3):
        axs[c].set_ylim(0, top)
    axs[3].set_ylim(0, None)  # the hub family has its own scale
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=4, frameon=False, handlelength=2.0, columnspacing=1.6,
               bbox_to_anchor=(0.5, 1.09))
    save(fig, "scaling_cmp")
    figdata("FIGURE 4: synthetic families, comparisons per arc as n grows "
            "(additions per arc are given as well; they are quoted in the text)",
            ["family", "log2_n"] + ["cmp_" + a for a in algs] + ["add_" + a for a in algs],
            [(fam, e) + tuple(Cn[(fam, e, a)]["cmp"] for a in algs) + tuple(Cn[(fam, e, a)]["add"] for a in algs)
             for fam, _ in fams for e in sorted({e for (f, e, a) in Cn if f == fam})])
    SUMMARY["scaling_count"] = cout
    SUMMARY["scaling_frontier_peak_over_n"] = {"%s/%d" % k: v for k, v in pk.items()}
    names = {"rand": "Rand", "grid": "Grid", "ba": "Ba", "hub": "Hub"}
    tags = {"dijkstra": "Dij", "dijkstra-dk": "DijDK", "bmssp25": "TwoFive", "bmssp26": "TwoSix"}
    slopes = {}
    cross = {}
    for fam, _ in fams:
        es = sorted({e for (f, e, a) in Cn if f == fam})
        hi = [e for e in es if e >= 16]
        macro("Peak%s" % names[fam], fmt(100 * pk[(fam, es[-1])], 1 if pk[(fam, es[-1])] < 0.1 else 0))
        for a in algs:
            ys = [Cn[(fam, e, a)]["cmp"] for e in hi]
            slope, icpt = np.polyfit(hi, ys, 1)
            slopes["%s/%s" % (fam, a)] = dict(slope=float(slope), intercept=float(icpt))
            macro("Slope%s%s" % (names[fam], tags[a]), fmt(slope, 2))
            macro("CmpLast%s%s" % (names[fam], tags[a]), fmt(ys[-1], 1))
            macro("AddLast%s%s" % (names[fam], tags[a]), fmt(Cn[(fam, es[-1], a)]["add"], 1))
            macro("ScanLast%s%s" % (names[fam], tags[a]), fmt(Cn[(fam, es[-1], a)]["scan"], 1))
        # smallest size from which the new algorithm needs fewer operations for all larger sizes
        for a in ("bmssp25", "bmssp26"):
            for base in ("dijkstra", "dijkstra-dk"):
                for what, f in (("Cmp", lambda c: c["cmp"]), ("Ops", lambda c: c["cmp"] + c["add"])):
                    first = None
                    for e in reversed(es):
                        if f(Cn[(fam, e, a)]) < f(Cn[(fam, e, base)]):
                            first = e
                        else:
                            break
                    cross["%s/%s/%s/%s" % (fam, a, base, what)] = first
                    macro("Cross%s%s%s%s" % (what, names[fam], tags[a], tags[base]), first if first else "none")
    SUMMARY["cmp_slopes_per_doubling"] = slopes
    for tag, a in (("Dij", "dijkstra"), ("TwoFive", "bmssp25"), ("TwoSix", "bmssp26")):
        c = Cn[("rand", 20, a)]
        macro("OpsRandTwenty" + tag, fmt(c["cmp"] + c["add"], 1))
        macro("CmpRandTwenty" + tag, fmt(c["cmp"], 1))
    for fam in ("rand", "ba"):
        for tag, a in (("TwoFive", "bmssp25"), ("TwoSix", "bmssp26")):
            ys = [Cn[(fam, e, a)]["cmp"] for e in sorted({e for (f, e, x) in Cn if f == fam})]
            macro("CmpRange%s%sMin" % (names[fam], tag), fmt(min(ys), 0)); macro("CmpRange%s%sMax" % (names[fam], tag), fmt(max(ys), 0))
    hubes = sorted({e for (f, e, x) in Cn if f == "hub"})
    for tag, a in (("TwoFive", "bmssp25"), ("TwoSix", "bmssp26")):
        r = [Cn[("hub", e, a)]["cmp"] / Cn[("hub", e, "dijkstra")]["cmp"] for e in hubes]
        macro("HubCmpRatioMin" + tag, fmt(min(r), 1)); macro("HubCmpRatioMax" + tag, fmt(max(r), 1))
    ges = sorted({e for (f, e, x) in Cn if f == "grid"})
    macro("PeakGridSmall", fmt(100 * pk[("grid", ges[0])], 1))
    SUMMARY["fewer_operations_from_log2n"] = cross
    # the text says "in both families"
    assert cross["rand/bmssp26/dijkstra/Cmp"] == cross["ba/bmssp26/dijkstra/Cmp"]


# ----------------------------------------------------------- density and weights
def sweep_figure():
    dt, wt = read("density_time"), read("weights_time")
    if not dt or not wt:
        return
    Td = per_source_time(dt, lambda r: (r["degree"], r["algo"]))
    Tw = per_source_time(wt, lambda r: (r["log2wmax"], r["algo"]))
    fig, axs = plt.subplots(1, 2, figsize=(COL1, 1.75), sharey=True, gridspec_kw=dict(wspace=0.08))
    out = {"density": {}, "weights": {}}
    for ax, T, xs, name, xlabel in ((axs[0], Td, sorted({k[0] for k in Td}), "density", "Average out-degree $m/n$"),
                                    (axs[1], Tw, sorted({k[0] for k in Tw}), "weights", "$\\log_2$ of largest weight")):
        for a in ("bmssp25", "bmssp26"):
            ys = [ratio(T, (x, a), (x, "dijkstra"))["gmean"] for x in xs]
            out[name][a] = dict(zip(map(str, xs), ys))
            ax.plot(xs, ys, color=STYLE[a]["color"], marker=STYLE[a]["marker"], ls=STYLE[a]["ls"],
                    label=STYLE[a]["label"], markeredgecolor="white", markeredgewidth=0.4, zorder=3)
        ax.axhline(1, color=STYLE["dijkstra"]["color"], lw=1.2, zorder=2, label="Dijkstra = 1")
        ax.set_xlabel(xlabel)
        ax.yaxis.grid(True, zorder=0)
        ax.set_axisbelow(True)
    axs[0].set_xscale("log", base=2)
    axs[0].set_xticks([2, 4, 8, 16, 32]); axs[0].set_xticklabels([2, 4, 8, 16, 32])
    axs[0].set_ylabel("Time relative to Dijkstra")
    axs[0].set_ylim(0, None)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=3, frameon=False, handlelength=1.8, columnspacing=1.2,
               bbox_to_anchor=(0.5, 1.1))
    save(fig, "sweeps", guide=True)
    SUMMARY["sweeps"] = out
    figdata("EXTRA (section VI-I of the paper, not a figure): random graphs with n = 2^18, time ratio by average out-degree",
            ["average_out_degree", "bmssp25", "bmssp26"],
            [(int(x), out["density"]["bmssp25"][x], out["density"]["bmssp26"][x])
             for x in sorted(out["density"]["bmssp25"], key=int)])
    figdata("EXTRA (section VI-I of the paper, not a figure): random graphs with n = 2^18, time ratio by weight range",
            ["log2_of_largest_weight", "bmssp25", "bmssp26"],
            [(int(x), out["weights"]["bmssp25"][x], out["weights"]["bmssp26"][x])
             for x in sorted(out["weights"]["bmssp25"], key=int)])
    for tag, a in (("TwoFive", "bmssp25"), ("TwoSix", "bmssp26")):
        d = out["density"][a]; w = out["weights"][a]
        ks = sorted(d, key=int)
        macro("DensLow" + tag, fmt(d[ks[0]], 1)); macro("DensHigh" + tag, fmt(d[ks[-1]], 1))
        macro("DensMin" + tag, fmt(min(d.values()), 1)); macro("DensMax" + tag, fmt(max(d.values()), 1))
        kw = sorted(w, key=int)
        macro("WeightUnit" + tag, fmt(w[kw[0]], 1)); macro("WeightWide" + tag, fmt(w[kw[-1]], 1))
        macro("WeightMin" + tag, fmt(min(w.values()), 1)); macro("WeightMax" + tag, fmt(max(w.values()), 1))


# ------------------------------------------------- hub graphs, degree reduction
def hub_table():
    trows, crows = read("hubreduce_time"), read("hubreduce_count")
    if not trows or not crows:
        return
    e = 18
    T = defaultdict(list)
    for r in trows:
        assert r["ok"] == 1
        if r["log2n"] == e:
            T[(r["variant"], r["algo"])].append(r["ns"])
    T = {k: st.median(v) for k, v in T.items()}
    M = (1 << e) * 4 - 1  # arcs of the graph as given
    C, N = {}, {}
    for r in crows:
        if r["log2n"] == e:
            C[(r["variant"], r["algo"])] = r["cmp"] / M
            N[r["variant"]] = r["n"]
    base = T[("original", "dijkstra")]
    lines = []
    out = {}
    for v, label in (("original", "as given"), ("delta3", "$\\delta=3$")):
        vals = [label, "{:,}".format(N[v]).replace(",", "\\,")]
        for a in ("dijkstra", "bmssp25", "bmssp26"):
            vals.append(fmt(C[(v, a)], 1))
            vals.append(fmt(T[(v, a)] / base, 1))
        lines.append(" & ".join(vals) + " \\\\")
        for a in ("dijkstra", "bmssp25", "bmssp26"):
            out["%s/%s" % (v, a)] = dict(cmp=C[(v, a)], time_vs_dijkstra_as_given=T[(v, a)] / base,
                                         time_vs_dijkstra_same_graph=T[(v, a)] / T[(v, "dijkstra")])
    write("hub", "\n".join(lines) + "\n")
    SUMMARY["hub_reduce_n2e18"] = out
    macro("HubCmpDij", fmt(C[("original", "dijkstra")], 0)); macro("HubCmpTwoFive", fmt(C[("original", "bmssp25")], 0))
    macro("HubCmpTwoSix", fmt(C[("original", "bmssp26")], 0))
    macro("HubRedCmpDij", fmt(C[("delta3", "dijkstra")], 0)); macro("HubRedCmpTwoFive", fmt(C[("delta3", "bmssp25")], 0))
    macro("HubRedCmpTwoSix", fmt(C[("delta3", "bmssp26")], 0))
    macro("HubRedSameTwoFive", fmt(out["delta3/bmssp25"]["time_vs_dijkstra_same_graph"], 1))
    macro("HubRedSameTwoSix", fmt(out["delta3/bmssp26"]["time_vs_dijkstra_same_graph"], 1))
    macro("HubRedDijSlow", fmt(out["delta3/dijkstra"]["time_vs_dijkstra_as_given"], 1))


# -------------------------------------------------------- SNAP graphs (real)
def snap_table():
    trows, crows = read("snap_time"), read("snap_count")
    if not trows or not crows:
        return
    T = per_source_time(trows, lambda r: (r["snapgraph"], r["weights"], r["algo"]))
    Cn = per_source_count(crows, lambda r: (r["snapgraph"], r["weights"], r["algo"]), ["cmp", "add", "scan"])
    size, peak, reach = {}, defaultdict(list), defaultdict(list)
    for r in crows:
        size[r["snapgraph"]] = (r["n"], r["m"])
        if r["algo"] == "dijkstra-dk":
            peak[(r["snapgraph"], r["weights"])].append(r["peak"] / r["n"])
            reach[r["snapgraph"]].append(r["reached"] / r["n"])
    names = []
    for r in crows:
        if r["snapgraph"] not in names:
            names.append(r["snapgraph"])
    lines, out = [], {}
    for g in names:
        n, m = size[g]
        vals = [g, "%.0f" % (n / 1e3), "%.2f" % (m / 1e6), "%.1f" % (100 * st.mean(peak[(g, "rand")]))]
        for w in ("unit", "rand"):
            for a in ("bmssp25", "bmssp26"):
                r = ratio(T, (g, w, a), (g, w, "dijkstra"))["gmean"]
                c = Cn[(g, w, a)]["cmp"] / Cn[(g, w, "dijkstra")]["cmp"]
                cdk = Cn[(g, w, a)]["cmp"] / Cn[(g, w, "dijkstra-dk")]["cmp"]
                out["%s/%s/%s" % (g, w, a)] = dict(time=r, cmp_vs_lazy=c, cmp_vs_dk=cdk,
                                                   cmp=Cn[(g, w, a)]["cmp"], add=Cn[(g, w, a)]["add"])
        for a in ("bmssp25", "bmssp26"):
            vals += [fmt(out["%s/unit/%s" % (g, a)]["time"], 1), fmt(out["%s/rand/%s" % (g, a)]["time"], 1),
                     fmt(out["%s/rand/%s" % (g, a)]["cmp_vs_lazy"], 2)]
        lines.append(" & ".join(vals) + " \\\\")
        out[g] = dict(n=n, m=m, frontier_unit=st.mean(peak[(g, "unit")]), frontier_rand=st.mean(peak[(g, "rand")]),
                      reach=st.mean(reach[g]),
                      dijkstra_cmp_unit=Cn[(g, "unit", "dijkstra")]["cmp"], dijkstra_cmp_rand=Cn[(g, "rand", "dijkstra")]["cmp"])
    write("snap", "\n".join(lines) + "\n")
    SUMMARY["snap"] = out
    rows = []
    for g in names:
        for w in ("unit", "rand"):
            for a in ("bmssp25", "bmssp26"):
                r = ratio(T, (g, w, a), (g, w, "dijkstra"))
                o = out["%s/%s/%s" % (g, w, a)]
                rows.append((g, size[g][0], size[g][1], "unit" if w == "unit" else "random", a, r["gmean"], r["lo"],
                             r["hi"], o["cmp_vs_lazy"], o["cmp_vs_dk"], 100 * st.mean(peak[(g, w)])))
    figdata("EXTRA (table IX of the paper, not a figure): SNAP graphs (web: Stanford, Google; social: Pokec; "
            "roads: roadPA, roadTX, roadCA), ten start vertices each",
            ["graph", "vertices", "arcs", "weights", "algorithm", "time_ratio", "smallest_per_start_vertex",
             "largest_per_start_vertex", "comparisons_vs_dijkstra", "comparisons_vs_dijkstra-dk",
             "dijkstra_peak_frontier_percent_of_n"], rows)
    inst = [(g, w) for g in names for w in ("unit", "rand")]
    words = ["none", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve"]
    macro("SnapInstances", words[len(inst)])
    macro("SnapGraphs", words[len(names)])
    for tag, a in (("TwoFive", "bmssp25"), ("TwoSix", "bmssp26")):
        t = [out["%s/%s/%s" % (g, w, a)]["time"] for g, w in inst]
        macro("SnapRatioMin" + tag, fmt(min(t), 1)); macro("SnapRatioMax" + tag, fmt(max(t), 1))
        fewer = [(g, w) for g, w in inst if out["%s/%s/%s" % (g, w, a)]["cmp_vs_lazy"] < 1]
        macro("SnapFewerCmp" + tag, words[len(fewer)])
        fewer_dk = [(g, w) for g, w in inst if out["%s/%s/%s" % (g, w, a)]["cmp_vs_dk"] < 1]
        macro("SnapFewerCmpDK" + tag, words[len(fewer_dk)])
        SUMMARY["snap_fewer_cmp_" + a] = ["%s/%s" % x for x in fewer]
        SUMMARY["snap_fewer_cmp_vs_dk_" + a] = ["%s/%s" % x for x in fewer_dk]
    for tag, a in (("TwoFive", "bmssp25"), ("TwoSix", "bmssp26")):
        for w, wn in (("unit", "Unit"), ("rand", "Rand")):
            t = [out["%s/%s/%s" % (g, w, a)]["time"] for g in names]
            macro("Snap%sTimeMin%s" % (wn, tag), fmt(min(t), 1)); macro("Snap%sTimeMax%s" % (wn, tag), fmt(max(t), 1))
            c = [out["%s/%s/%s" % (g, w, a)]["cmp_vs_lazy"] for g in names]
            macro("Snap%sCmpMin%s" % (wn, tag), fmt(min(c), 2)); macro("Snap%sCmpMax%s" % (wn, tag), fmt(max(c), 2))
    fl = [g for g, w in inst if w == "rand" and out["%s/rand/bmssp25" % g]["cmp_vs_lazy"] < 1]
    macro("SnapFewerListTwoFive", " and ".join(fl) if fl else "none")
    # statements made in the text; fail loudly if regenerated data no longer support them
    assert all(out["%s/unit/%s" % (g, a)]["cmp_vs_lazy"] > 1 for g in names for a in ("bmssp25", "bmssp26"))
    assert all(w == "rand" for g, w in inst if out["%s/%s/bmssp25" % (g, w)]["cmp_vs_lazy"] < 1)
    best = min(inst, key=lambda x: out["%s/%s/bmssp25" % x]["time"])
    macro("SnapBestInstance", "%s with %s weights" % (best[0], "random" if best[1] == "rand" else "unit"))
    macro("SnapTwoSixFaster", words[sum(1 for g, w in inst if out["%s/%s/bmssp26" % (g, w)]["time"] < out["%s/%s/bmssp25" % (g, w)]["time"])])
    roads = [g for g in names if g.startswith("road")]
    others = [g for g in names if not g.startswith("road")]
    macro("SnapRoadPeakMax", fmt(100 * max(out[g]["frontier_rand"] for g in roads), 1))
    macro("SnapOtherPeakMin", fmt(100 * min(out[g]["frontier_rand"] for g in others), 0))
    macro("SnapOtherPeakMax", fmt(100 * max(out[g]["frontier_rand"] for g in others), 0))
    for grp, gs in (("Road", roads), ("Other", others)):
        for tag, a in (("TwoFive", "bmssp25"), ("TwoSix", "bmssp26")):
            c = [out["%s/%s/%s" % (g, w, a)]["cmp_vs_lazy"] for g in gs for w in ("unit", "rand")]
            macro("SnapCmp%sMin%s" % (grp, tag), fmt(min(c), 2)); macro("SnapCmp%sMax%s" % (grp, tag), fmt(max(c), 2))
            t = [out["%s/%s/%s" % (g, w, a)]["time"] for g in gs for w in ("unit", "rand")]
            macro("SnapTime%sMin%s" % (grp, tag), fmt(min(t), 1)); macro("SnapTime%sMax%s" % (grp, tag), fmt(max(t), 1))


# ---------------------------------------------------------------------- memory
def memory_table():
    rows = read("memory")
    if not rows:
        return
    M = {(r["graph"], r["algo"]): r["peak_rss_bytes"] / 2**20 for r in rows}
    lines = []
    for g in ("SF", "NA", "rand-20", "rand-22"):
        label = {"rand-20": "Random, $n=2^{20}$", "rand-22": "Random, $n=2^{22}$"}.get(g, g)
        lines.append("%s & %s & %s & %s & %s \\\\" % (label, fmt(M[(g, "dijkstra")], 0), fmt(M[(g, "dijkstra-dk")], 0),
                                                    fmt(M[(g, "bmssp25")], 0), fmt(M[(g, "bmssp26")], 0)))
    write("memory", "\n".join(lines) + "\n")
    SUMMARY["memory_mb"] = {"%s/%s" % k: v for k, v in M.items()}
    macro("MemDij", fmt(M[("rand-22", "dijkstra")], 0))
    macro("MemTwoFive", fmt(M[("rand-22", "bmssp25")], 0))
    macro("MemTwoSix", fmt(M[("rand-22", "bmssp26")], 0))
    macro("MemRatioTwoFive", fmt(M[("rand-22", "bmssp25")] / M[("rand-22", "dijkstra")], 1))
    macro("MemRatioTwoSix", fmt(M[("rand-22", "bmssp26")] / M[("rand-22", "dijkstra")], 1))


# ------------------------------------------------- parameters at feasible sizes
def theory_table():
    """Parameter values prescribed by the two papers, evaluated for concrete n (delta = 3)."""
    lines = []
    out = {}
    for e in (10, 20, 30, 40, 64):
        lg = float(e)
        k25 = max(1, math.floor(lg ** (1 / 3) + 1e-9))
        t25 = max(1, math.floor(lg ** (2 / 3) + 1e-9))
        l25 = math.ceil(lg / t25 - 1e-9)
        t26 = max(2, math.ceil(math.sqrt(lg * math.log2(lg) / 3.0) - 1e-9))
        k26 = max(1, math.ceil(t26 / math.log2(t26) - 1e-9)) if t26 > 2 else t26
        l26 = math.ceil(lg / t26 - 1e-9)
        g25 = lg ** (1 / 3)
        g26 = math.sqrt(lg / math.log2(lg))
        out[e] = dict(k25=k25, t25=t25, l25=l25, t26=t26, k26=k26, l26=l26, gain25=g25, gain26=g26)
        lines.append("$2^{%d}$ & %d & %d & %d & %s & %d & %d & %d & %s & %s \\\\" % (
            e, k25, t25, l25, fmt(g25), k26, t26, l26, "{:,}".format(t26 ** 3).replace(",", "\\,"), fmt(g26)))
    write("theory_params", "\n".join(lines) + "\n")
    SUMMARY["theory_params"] = out
    # log n beyond which sqrt(log n loglog n) < log^{2/3} n, i.e. (loglog n)^3 < log n
    x = 16.0
    while math.log2(x) ** 3 >= x:
        x += 1.0
    macro("FaceValueCrossLog", int(x))
    macro("GainTwoFiveAtTwenty", fmt(20 ** (1 / 3), 1)); macro("GainTwoSixAtTwenty", fmt(math.sqrt(20 / math.log2(20)), 1))
    macro("GainTwoFiveAtSixtyFour", fmt(64 ** (1 / 3), 1)); macro("GainTwoSixAtSixtyFour", fmt(math.sqrt(64 / 6), 1))


def validation_table():
    sets = []
    for name, label in (("small_n16", "$n\\le 16$"), ("medium_n400", "$n\\le 400$")):
        path = os.path.join(RES, "validation", name + ".csv")
        if not os.path.exists(path):
            return
        with open(path) as f:
            sets.append((label, {r["variant"]: {k: int(x) for k, x in r.items() if k != "variant"}
                                 for r in csv.DictReader(f)}))
    rows = [("BMSSP-25, as printed", "bmssp25/literal/paper-params"),
            ("BMSSP-25, with deletion", "bmssp25/fix/paper-params"),
            ("BMSSP-26, as printed", "bmssp26/literal/paper-params"),
            ("BMSSP-26, both repairs", "bmssp26/fix+basefix/paper-params")]
    lines = []
    num = lambda x: "{:,}".format(x).replace(",", "\\,")
    for label, v in sets:
        for i, (name, key) in enumerate(rows):
            d = v[key]
            lines.append("%s & %s & %s & %s & %s & %s & %s \\\\" % (
                label if i == 0 else "", name, num(d["runs"]), num(d["wrong"]), num(d["runs_with_stale_keys"]),
                num(d["runs_with_duplicates"]), num(d["runs_with_bound_violation"])))
        lines.append("\\hline")
    write("validation", "\n".join(lines[:-1]) + "\n")


if __name__ == "__main__":
    theory_table()
    validation_table()
    road_tables_and_figures()
    variants_tables()
    reduce_table()
    params_figure()
    scaling_figures()
    sweep_figure()
    hub_table()
    snap_table()
    memory_table()
    write_figdata()
    with open(os.path.join(RES, "environment.json")) as f:
        SUMMARY["environment"] = json.load(f)
    with open(os.path.join(RES, "summary.json"), "w") as f:
        json.dump(SUMMARY, f, indent=1, sort_keys=True)
    with open(os.path.join(ROOT, "paper", "numbers.tex"), "w") as f:
        f.write("% Generated by scripts/make_figures.py from results/*.csv. Do not edit by hand.\n")
        for k in sorted(NUM):
            f.write("\\newcommand{\\%s}{%s}\n" % (k, NUM[k]))
    print("  wrote paper/numbers.tex (%d macros) and results/summary.json" % len(NUM))
