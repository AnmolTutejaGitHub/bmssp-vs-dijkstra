#!/usr/bin/env python3
"""Run every experiment of the study and store the raw measurements as CSV.

    python3 scripts/run_experiments.py              # everything
    python3 scripts/run_experiments.py real params  # selected experiments
    python3 scripts/run_experiments.py --force real # overwrite existing results
    python3 scripts/run_experiments.py --force --counts-only   # redo only the operation counts

Experiments (each writes results/<name>_time.csv and/or results/<name>_count.csv):
    info      properties of the road networks
    validate  randomized differential testing of all implementation variants
    real      four algorithms on the five road networks (main comparison)
    variants  pseudocode-as-printed versus repaired variants on the road networks
    reduce    effect of the degree-reduction transformation
    params    sweep of the parameters k and t on the two largest road networks
    scaling   synthetic families, n = 2^10 ... 2^22
    density   random digraphs with n = 2^18 and growing average degree
    weights   random digraphs with n = 2^18 and growing weight range
    memory    peak resident set size on the largest instances
    hubreduce hub graphs with and without degree reduction
    snap      six graphs of the SNAP collection (needs `sh scripts/fetch_data.sh snap`)

Timing uses bin/sssp, operation counts use bin/sssp_count (see the Makefile).
Only the Python standard library is needed. Run on an otherwise idle machine.
"""
import json
import os
import platform
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "scripts", "bin")
RES = os.path.join(ROOT, "results")
DATA = os.path.join(ROOT, "data", "sstd2005", "Datasets - SSTD 2005")
SEED = 20261012

# short name -> (full name, edge file)
ROADS = {
    "OL": ("Oldenburg", "OL Road Network's Edges (Edge ID, Start Node ID, End Node ID, L2 Distance).txt"),
    "TG": ("San Joaquin County", "TG Road Network's Edges (Edge ID, Start Node ID, End Node ID, L2 Distance).txt"),
    "CAL": ("California", "California Road Network's Edges (Edge ID, Start Node ID, End Node ID, L2 Distance).txt"),
    "SF": ("San Francisco", "San Francisco Road Network's Edges (Edge ID, Start Node ID, End Node ID, L2 Distance).txt"),
    "NA": ("North America", "North America Road Network's Edges (Edge ID, Start Node ID, End Node ID, L2 Distance).txt"),
}
ALGOS = ["dijkstra", "dijkstra-dk", "bmssp25", "bmssp26"]


def road_spec(short):
    return "sstd:" + os.path.join(DATA, ROADS[short][1])


def call(binary, args):
    out = subprocess.run([os.path.join(BIN, binary)] + [str(a) for a in args], check=True,
                         capture_output=True, text=True).stdout
    return [l for l in out.splitlines() if l.strip()]


class Sink:
    """Collects the CSV rows printed by several invocations into one file."""

    def __init__(self, name, force):
        self.path = os.path.join(RES, name + ".csv")
        self.skip = os.path.exists(self.path) and not force
        self.header = None
        self.rows = []

    def add(self, lines, extra=None):
        hdr, rows = lines[0], lines[1:]
        if extra:
            hdr += "," + ",".join(extra.keys())
            rows = [r + "," + ",".join(str(v) for v in extra.values()) for r in rows]
        assert self.header in (None, hdr), "inconsistent csv header"
        self.header = hdr
        self.rows += rows

    def close(self):
        if self.skip or self.header is None:
            return
        with open(self.path, "w") as f:
            f.write(self.header + "\n" + "\n".join(self.rows) + "\n")
        print("  wrote %s (%d rows)" % (os.path.relpath(self.path, ROOT), len(self.rows)), flush=True)


COUNTS_ONLY = False  # set by --counts-only: redo the operation counts, keep the timings


def both(name, force, jobs, time_args, count_args, timing=True, counting=True):
    """jobs: list of (tag, spec, algo, extra args, extra columns)."""
    for binary, base, on in (("sssp", time_args, timing), ("sssp_count", count_args, counting)):
        if not on or (COUNTS_ONLY and binary == "sssp"):
            continue
        sink = Sink(name + ("_time" if binary == "sssp" else "_count"), force)
        if sink.skip:
            print("  %s exists, skipping (use --force to redo)" % os.path.relpath(sink.path, ROOT))
            continue
        for i, (tag, spec, algo, extra, cols) in enumerate(jobs):
            t0 = time.time()
            merged = list(extra) + base
            sink.add(call(binary, ["run", spec, "--algo", algo, "--tag", tag, "--seed", SEED] + merged), cols)
            print("    [%s %d/%d] %s %s %s  %.1fs" % (binary, i + 1, len(jobs), tag, algo, " ".join(map(str, extra)),
                                                     time.time() - t0), flush=True)
        sink.close()


def exp_info(force):
    path = os.path.join(RES, "graph_info.csv")
    if os.path.exists(path) and not force:
        print("  graph_info.csv exists, skipping")
        return
    rows = []
    for short in ROADS:
        kv = dict(p.split("=") for p in call("sssp", ["info", road_spec(short)])[0].split(","))
        kv = {"graph": short, "name": ROADS[short][0], **kv}
        rows.append(kv)
    with open(path, "w") as f:
        f.write(",".join(rows[0].keys()) + "\n")
        for r in rows:
            f.write(",".join(r.values()) + "\n")
    print("  wrote results/graph_info.csv")


def exp_validate(force):
    os.makedirs(os.path.join(RES, "validation"), exist_ok=True)
    for name, args in (("small_n16", ["--trials", 200000, "--maxn", 16, "--seed", SEED, "--dump",
                                      os.path.join(RES, "validation")]),
                       ("medium_n400", ["--trials", 30000, "--maxn", 400, "--seed", 777])):
        path = os.path.join(RES, "validation", name + ".csv")
        if os.path.exists(path) and not force:
            print("  validation/%s.csv exists, skipping" % name)
            continue
        lines = call("sssp_count", ["validate"] + args)
        with open(path, "w") as f:
            f.write("\n".join(lines) + "\n")
        print("  wrote results/validation/%s.csv" % name)


def exp_real(force):
    jobs = [(g, road_spec(g), a, [], None) for g in ROADS for a in ALGOS]
    both("real", force, jobs, ["--sources", 20, "--reps", 9], ["--sources", 20])


def exp_variants(force):
    variants = [
        ("bmssp25", "repaired", []),
        ("bmssp25", "as-printed", ["--nofix"]),
        ("bmssp26", "repaired", []),
        ("bmssp26", "no-deletion", ["--nofix"]),
        ("bmssp26", "printed-base-case", ["--literal-base"]),
        ("bmssp26", "as-printed", ["--nofix", "--literal-base"]),
    ]
    jobs = [(g, road_spec(g), a, flags + ["--variant", v], None) for g in ROADS for a, v, flags in variants]
    both("variants", force, jobs, ["--sources", 20, "--reps", 9], ["--sources", 20])


def exp_reduce(force):
    jobs = []
    for g in ROADS:
        for a in ("dijkstra", "bmssp25", "bmssp26"):
            jobs.append((g, road_spec(g), a, ["--variant", "original"], None))
            jobs.append((g, road_spec(g), a, ["--reduce", 3, "--variant", "delta3"], None))
            if a != "bmssp25":
                jobs.append((g, road_spec(g), a, ["--reduce", 4, "--delta", 4, "--variant", "delta4"], None))
    both("reduce", force, jobs, ["--sources", 10, "--reps", 7], ["--sources", 10])


def exp_params(force):
    jobs = []
    for g in ("SF", "NA"):
        for k in range(1, 7):
            for t in range(2, 11):
                jobs.append((g, road_spec(g), "bmssp25", ["--k", k, "--t", t], None))
                jobs.append((g, road_spec(g), "bmssp26", ["--k", k, "--t", t], None))
        jobs.append((g, road_spec(g), "dijkstra", [], None))
    both("params", force, jobs, ["--sources", 5, "--reps", 5], ["--sources", 5])


def synth_specs():
    specs = []
    for e in range(10, 23):
        specs.append(("rand", e, "rand:%d:3:%d" % (1 << e, e), []))
    for e in range(10, 23, 2):
        specs.append(("grid", e, "grid:%d:%d" % (1 << (e // 2), e), []))
    for e in range(10, 22):
        specs.append(("ba", e, "ba:%d:2:%d" % (1 << e, e), []))
    for e in range(10, 22):
        specs.append(("hub", e, "hub:%d:3:%d" % (1 << e, e), ["--source0"]))
    return specs


def exp_scaling(force):
    jobs = []
    for fam, e, spec, extra in synth_specs():
        for a in ALGOS:
            src = extra if extra else ["--sources", 3]
            jobs.append(("%s-%d" % (fam, e), spec, a, src, {"family": fam, "log2n": e}))
    both("scaling", force, jobs, ["--reps", 5], [])


def exp_density(force):
    jobs = []
    for deg in (2, 3, 4, 6, 8, 12, 16, 24, 32):
        for a in ALGOS:
            jobs.append(("rand-18-d%d" % deg, "rand:%d:%d:%d" % (1 << 18, deg, deg), a, ["--sources", 3],
                         {"degree": deg}))
    both("density", force, jobs, ["--reps", 5], [])


def exp_weights(force):
    jobs = []
    for wexp in (0, 1, 4, 10, 20, 30):
        for a in ALGOS:
            jobs.append(("rand-18-w%d" % wexp, "rand:%d:3:18:%d" % (1 << 18, 1 << wexp), a, ["--sources", 3],
                         {"log2wmax": wexp}))
    both("weights", force, jobs, ["--reps", 5], [])


def exp_hubreduce(force):
    """Hub graphs as given (source of out-degree n-1) and after degree reduction with delta = 3."""
    jobs = []
    for e in (14, 16, 18):
        spec = "hub:%d:3:%d" % (1 << e, e)
        for a in ALGOS:
            jobs.append(("hub-%d" % e, spec, a, ["--source0", "--variant", "original"], {"log2n": e}))
            jobs.append(("hub-%d" % e, spec, a, ["--source0", "--reduce", 3, "--variant", "delta3"], {"log2n": e}))
    both("hubreduce", force, jobs, ["--reps", 5], [])


SNAP = [  # short name, file in data/snap (downloaded by scripts/fetch_data.sh)
    ("Stanford", "web-Stanford.txt.gz"), ("Google", "web-Google.txt.gz"), ("roadPA", "roadNet-PA.txt.gz"),
    ("roadTX", "roadNet-TX.txt.gz"), ("roadCA", "roadNet-CA.txt.gz"), ("Pokec", "soc-pokec-relationships.txt.gz"),
]


def exp_snap(force):
    """Web, road and social graphs of the SNAP collection, with unit and with random weights."""
    jobs = []
    for name, fname in SNAP:
        path = os.path.join(ROOT, "data", "snap", fname)
        if not os.path.exists(path):
            sys.exit("missing %s: run `sh scripts/fetch_data.sh snap` first" % path)
        for mode in ("unit", "rand"):
            for a in ALGOS:
                jobs.append(("%s-%s" % (name, mode), "snap:%s:%s" % (path, mode), a,
                             ["--sources", 10, "--min-reach", 0.25], {"snapgraph": name, "weights": mode}))
    both("snap", force, jobs, ["--reps", 3], [])


def exp_memory(force):
    """Peak resident set size of one complete run (graph + reference Dijkstra + algorithm)."""
    path = os.path.join(RES, "memory.csv")
    if os.path.exists(path) and not force:
        print("  memory.csv exists, skipping")
        return
    graphs = [("SF", road_spec("SF")), ("NA", road_spec("NA")), ("rand-20", "rand:%d:3:20" % (1 << 20)),
              ("rand-22", "rand:%d:3:22" % (1 << 22))]
    mac = platform.system() == "Darwin"
    rows = ["graph,algo,peak_rss_bytes"]
    for tag, spec in graphs:
        for a in ALGOS:
            cmd = ["/usr/bin/time", "-l" if mac else "-v", os.path.join(BIN, "sssp"), "run", spec, "--algo", a,
                   "--sources", "1", "--reps", "1", "--seed", str(SEED)]
            err = subprocess.run(cmd, check=True, capture_output=True, text=True).stderr
            rss = None
            for line in err.splitlines():
                if "maximum resident set size" in line:            # macOS: bytes
                    rss = int(line.split()[0])
                elif "Maximum resident set size" in line:          # GNU time: kilobytes
                    rss = int(line.split()[-1]) * 1024
            rows.append("%s,%s,%d" % (tag, a, rss))
            print("    %s %s %.1f MB" % (tag, a, rss / 2**20), flush=True)
    with open(path, "w") as f:
        f.write("\n".join(rows) + "\n")
    print("  wrote results/memory.csv")


def environment():
    def sh(cmd):
        try:
            return subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip()
        except Exception:
            return ""
    env = {
        "date": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu": sh("sysctl -n machdep.cpu.brand_string") or platform.processor(),
        "cores": os.cpu_count(),
        "memory_bytes": sh("sysctl -n hw.memsize"),
        "compiler": sh("c++ --version | head -1"),
        "flags": "-std=c++17 -O3 -DNDEBUG",
        "python": platform.python_version(),
        "dataset_commit": sh("git -C '%s' rev-parse HEAD" % os.path.join(ROOT, "data", "sstd2005")),
        "seed": SEED,
    }
    with open(os.path.join(RES, "environment.json"), "w") as f:
        json.dump(env, f, indent=2)


EXPERIMENTS = {
    "info": exp_info, "validate": exp_validate, "real": exp_real, "variants": exp_variants,
    "reduce": exp_reduce, "params": exp_params, "scaling": exp_scaling, "density": exp_density,
    "weights": exp_weights, "memory": exp_memory, "hubreduce": exp_hubreduce, "snap": exp_snap,
}

if __name__ == "__main__":
    args = sys.argv[1:]
    force = "--force" in args
    COUNTS_ONLY = "--counts-only" in args
    names = [a for a in args if not a.startswith("--")]
    if not names:  # everything; the SNAP graphs only if they have been downloaded
        have_snap = all(os.path.exists(os.path.join(ROOT, "data", "snap", f)) for _, f in SNAP)
        names = [n for n in EXPERIMENTS if n != "snap" or have_snap]
    os.makedirs(RES, exist_ok=True)
    for b in ("sssp", "sssp_count"):
        if not os.path.exists(os.path.join(BIN, b)):
            sys.exit("missing scripts/bin/%s: run `make -C scripts` first" % b)
    if not os.path.isdir(DATA):
        sys.exit("missing data set: run `sh scripts/fetch_data.sh` first")
    environment()
    for n in names:
        if n not in EXPERIMENTS:
            sys.exit("unknown experiment %s (choose from %s)" % (n, ", ".join(EXPERIMENTS)))
        print("== %s" % n, flush=True)
        t0 = time.time()
        EXPERIMENTS[n](force)
        print("   done in %.1f s" % (time.time() - t0), flush=True)
