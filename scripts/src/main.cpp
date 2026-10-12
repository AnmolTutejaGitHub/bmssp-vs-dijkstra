// Driver: timing runs, operation counts and randomized validation.
//
//   sssp info <graph>
//   sssp run <graph> --algo <name> [--k K] [--t T] [--delta D] [--reduce D]
//            [--nofix] [--literal-base] [--sources N] [--min-reach F] [--source0] [--seed S] [--reps R]
//   sssp validate [--trials N] [--seed S] [--maxn N]
//
// Algorithms: dijkstra (lazy binary heap), dijkstra-dk (indexed heap with
// DecreaseKey), bmssp25 (Duan et al. 2025), bmssp26 (Duan et al. 2026).
// Every run is checked against Dijkstra before its result line is printed.
#include <cstring>
#include <iostream>
#include <map>
#include <memory>
#include "bmssp25.hpp"
#include "bmssp26.hpp"
#include "dijkstra.hpp"
#include "graphs.hpp"

namespace sp {
Counters C;
}
using namespace sp;

struct Args {
  std::map<std::string, std::string> kv;
  std::vector<std::string> pos;
  Args(int argc, char** argv, int from) {
    for (int i = from; i < argc; i++) {
      std::string a = argv[i];
      if (a.rfind("--", 0) == 0) {
        std::string key = a.substr(2);
        if (i + 1 < argc && std::strncmp(argv[i + 1], "--", 2) != 0)
          kv[key] = argv[++i];
        else
          kv[key] = "1";
      } else {
        pos.push_back(a);
      }
    }
  }
  bool has(const std::string& k) const { return kv.count(k) > 0; }
  long long num(const std::string& k, long long def) const {
    auto it = kv.find(k);
    return it == kv.end() ? def : std::stoll(it->second);
  }
  double real(const std::string& k, double def) const {
    auto it = kv.find(k);
    return it == kv.end() ? def : std::stod(it->second);
  }
  std::string str(const std::string& k, const std::string& def) const {
    auto it = kv.find(k);
    return it == kv.end() ? def : it->second;
  }
};

// Uniform interface over the four algorithms.
struct Solver {
  virtual ~Solver() {}
  virtual void run(vid s) = 0;
  virtual Len dist(vid v) const = 0;
  virtual std::string params() const { return "k=0,t=0,levels=0"; }
};
struct SolverDijkstra : Solver {
  DijkstraLazy a;
  explicit SolverDijkstra(const Graph& g) : a(g) {}
  void run(vid s) override { a.run(s); }
  Len dist(vid v) const override { return a.dist[v]; }
};
struct SolverDijkstraDK : Solver {
  DijkstraIndexed a;
  explicit SolverDijkstraDK(const Graph& g) : a(g) {}
  void run(vid s) override { a.run(s); }
  Len dist(vid v) const override { return a.dist[v]; }
};
struct Solver25 : Solver {
  Bmssp25 a;
  Solver25(const Graph& g, int k, int t, bool fix) : a(g, k, t) { a.fix = fix; }
  void run(vid s) override { a.run(s); }
  Len dist(vid v) const override { return a.dist(v); }
  std::string params() const override {
    return "k=" + std::to_string(a.k) + ",t=" + std::to_string(a.t) + ",levels=" + std::to_string(a.top);
  }
};
struct Solver26 : Solver {
  Bmssp26 a;
  Solver26(const Graph& g, int t, int k, double delta, bool fix, bool baseFix) : a(g, t, k, delta) {
    a.fix = fix;
    a.baseFix = baseFix;
  }
  void run(vid s) override { a.run(s); }
  Len dist(vid v) const override { return a.dist(v); }
  std::string params() const override {
    return "k=" + std::to_string(a.k) + ",t=" + std::to_string(a.t) + ",levels=" + std::to_string(a.top);
  }
};

static std::unique_ptr<Solver> makeSolver(const std::string& algo, const Graph& g, const Args& A) {
  int k = (int)A.num("k", 0), t = (int)A.num("t", 0);
  bool fix = !A.has("nofix"), baseFix = !A.has("literal-base");
  if (algo == "dijkstra") return std::unique_ptr<Solver>(new SolverDijkstra(g));
  if (algo == "dijkstra-dk") return std::unique_ptr<Solver>(new SolverDijkstraDK(g));
  if (algo == "bmssp25") return std::unique_ptr<Solver>(new Solver25(g, k, t, fix));
  if (algo == "bmssp26") return std::unique_ptr<Solver>(new Solver26(g, t, k, A.real("delta", 3.0), fix, baseFix));
  throw std::runtime_error("unknown algorithm " + algo);
}

static void graphStats(const Graph& g, uint32_t& maxOut, uint32_t& maxIn) {
  std::vector<uint32_t> in(g.n, 0);
  maxOut = maxIn = 0;
  for (vid u = 0; u < g.n; u++) {
    maxOut = std::max(maxOut, g.outdeg(u));
    for (uint32_t e = g.head[u]; e < g.head[u + 1]; e++) in[g.to[e]]++;
  }
  for (vid u = 0; u < g.n; u++) maxIn = std::max(maxIn, in[u]);
}

static int cmdInfo(const Args& A) {
  Graph g = makeGraph(A.pos.at(0));
  uint32_t mo, mi;
  graphStats(g, mo, mi);
  Len wmin = INF_LEN, wmax = 0;
  for (Len w : g.w) {
    wmin = std::min(wmin, w);
    wmax = std::max(wmax, w);
  }
  // eccentricity of vertex 0 in hops (BFS) as a rough depth indicator
  std::vector<uint32_t> lev(g.n, NONE);
  std::vector<vid> q{0};
  lev[0] = 0;
  uint32_t depth = 0;
  size_t reached = 1;
  for (size_t i = 0; i < q.size(); i++) {
    vid u = q[i];
    for (uint32_t e = g.head[u]; e < g.head[u + 1]; e++)
      if (lev[g.to[e]] == NONE) {
        lev[g.to[e]] = lev[u] + 1;
        depth = std::max(depth, lev[u] + 1);
        q.push_back(g.to[e]);
        reached++;
      }
  }
  size_t deg2 = 0;  // vertices with exactly two outgoing arcs (interior points of chains in a symmetric graph)
  for (vid u = 0; u < g.n; u++) deg2 += g.outdeg(u) == 2;
  int k25, t25, t26, k26;
  Bmssp25::defaults(g.n, k25, t25);
  Bmssp26::defaults(g.n, 3.0, t26, k26);
  std::printf("n=%u,m=%llu,max_out=%u,max_in=%u,wmin=%llu,wmax=%llu,bfs_depth0=%u,reach0=%zu,"
              "k25=%d,t25=%d,t26=%d,k26=%d,deg2_pct=%.1f\n",
              g.n, (unsigned long long)g.m, mo, mi, (unsigned long long)wmin, (unsigned long long)wmax, depth,
              reached, k25, t25, t26, k26, 100.0 * deg2 / g.n);
  return 0;
}

static int cmdRun(const Args& A) {
  std::string spec = A.pos.at(0);
  std::string algo = A.str("algo", "dijkstra");
  std::string tag = A.str("tag", spec);
  Graph orig = makeGraph(spec);
  int reduce = (int)A.num("reduce", 0);
  Graph red;
  std::vector<vid> rep;
  uint64_t prep_ns = 0;
  (void)prep_ns;
  if (reduce) {
    uint64_t t0 = now_ns();
    red = reduceDegree(orig, (unsigned)reduce, rep);
    prep_ns = now_ns() - t0;
  }
  const Graph& g = reduce ? red : orig;
  int nsrc = (int)A.num("sources", 5), reps = (int)A.num("reps", 5);
  uint64_t seed = (uint64_t)A.num("seed", 1);
  DijkstraLazy ref(orig);
  std::vector<vid> sources;
  if (A.has("source0")) {
    sources.push_back(0);
  } else {
    // Uniformly random sources. With --min-reach F a candidate is kept only if
    // it reaches at least the fraction F of all vertices (for directed graphs
    // in which a random vertex may reach almost nothing).
    double minReach = A.real("min-reach", 0.0);
    Rng r(seed * 0x9e3779b97f4a7c15ULL + 12345);
    for (int tries = 0; (int)sources.size() < nsrc && tries < 2000; tries++) {
      vid c = (vid)r.below(orig.n);
      if (minReach > 0) {
        ref.run(c);
        size_t reached = 0;
        for (vid v = 0; v < orig.n; v++) reached += ref.dist[v] != INF_LEN;
        if ((double)reached < minReach * orig.n) continue;
      }
      sources.push_back(c);
    }
  }
  std::unique_ptr<Solver> S = makeSolver(algo, g, A);
#ifdef COUNT_OPS
  reps = 1;
  std::printf("graph,algo,variant,n,m,k,t,levels,source,ok,reached,cmp,add,scan,relaxed,heap_push,heap_pop,"
              "ds_insert,ds_batch,ds_pulled,ds_pull,ds_split,sel_elems,scan_pivot,scan_base,scan_main,calls,"
              "base_calls,partial,pivots,frontier,stale_erase,stale_pull,dup_u,u_total,bound_viol,full_nonempty,"
              "stale_parent,peak\n");
#else
  std::printf("graph,algo,variant,n,m,k,t,levels,source,rep,ns,ok,prep_ns\n");
#endif
  std::string variant = A.str("variant", "default");
  // parameter triple "k=..,t=..,levels=.." -> three csv fields
  std::string ps = S->params();
  int pk = 0, pt = 0, pl = 0;
  std::sscanf(ps.c_str(), "k=%d,t=%d,levels=%d", &pk, &pt, &pl);
  int bad = 0;
  for (vid s : sources) {
    ref.run(s);
    size_t reached = 0;
    (void)reached;
    for (vid v = 0; v < orig.n; v++) reached += ref.dist[v] != INF_LEN;
    vid src = reduce ? rep[s] : s;
    for (int r = 0; r < reps; r++) {
      C.clear();
      uint64_t t0 = now_ns();
      S->run(src);
      uint64_t ns = now_ns() - t0;
      Counters snap = C;
      bool ok = true;
      for (vid v = 0; v < orig.n && ok; v++) ok = S->dist(reduce ? rep[v] : v) == ref.dist[v];
      bad += !ok;
#ifdef COUNT_OPS
      (void)ns;
      std::printf("%s,%s,%s,%u,%llu,%d,%d,%d,%u,%d,%zu", tag.c_str(), algo.c_str(), variant.c_str(), g.n,
                  (unsigned long long)g.m, pk, pt, pl, s, (int)ok, reached);
      for (uint64_t x : {snap.cmp, snap.add, snap.scan, snap.relaxed, snap.heap_push, snap.heap_pop, snap.ds_insert,
                         snap.ds_batch, snap.ds_pulled, snap.ds_pull, snap.ds_split, snap.sel_elems, snap.scan_pivot,
                         snap.scan_base, snap.scan_main, snap.calls, snap.base_calls, snap.partial, snap.pivots,
                         snap.frontier, snap.stale_erase, snap.stale_pull, snap.dup_u, snap.u_total, snap.bound_viol,
                         snap.full_nonempty, snap.stale_parent, snap.peak})
        std::printf(",%llu", (unsigned long long)x);
      std::printf("\n");
#else
      (void)snap;
      std::printf("%s,%s,%s,%u,%llu,%d,%d,%d,%u,%d,%llu,%d,%llu\n", tag.c_str(), algo.c_str(), variant.c_str(), g.n,
                  (unsigned long long)g.m, pk, pt, pl, s, r, (unsigned long long)ns, (int)ok,
                  (unsigned long long)prep_ns);
#endif
    }
  }
  std::fflush(stdout);
  return bad ? 2 : 0;
}

// ------------------------------------------------------------- validation
static std::vector<Len> bellmanFord(const Graph& g, vid s) {
  std::vector<Len> d(g.n, INF_LEN);
  d[s] = 0;
  for (vid round = 0; round < g.n; round++) {
    bool ch = false;
    for (vid u = 0; u < g.n; u++) {
      if (d[u] == INF_LEN) continue;
      for (uint32_t e = g.head[u]; e < g.head[u + 1]; e++)
        if (d[u] + g.w[e] < d[g.to[e]]) {
          d[g.to[e]] = d[u] + g.w[e];
          ch = true;
        }
    }
    if (!ch) break;
  }
  return d;
}

static Graph randomSmall(Rng& r, vid maxn) {
  vid n = 1 + (vid)r.below(maxn);
  int shape = (int)r.below(6);
  uint64_t wr;
  switch (r.below(5)) {
    case 0: wr = 1; break;      // all weights 0: every path ties
    case 1: wr = 2; break;      // weights 0/1
    case 2: wr = 4; break;
    case 3: wr = 1000; break;
    default: wr = 1ULL << 30;
  }
  std::vector<Edge> es;
  auto add = [&](vid u, vid v) { es.push_back({u, v, r.below(wr)}); };
  if (shape == 0) {  // sparse random
    uint64_t m = r.below(3 * (uint64_t)n + 1);
    for (uint64_t i = 0; i < m; i++) add((vid)r.below(n), (vid)r.below(n));
  } else if (shape == 1) {  // dense random
    uint64_t m = r.below((uint64_t)n * n + 1);
    for (uint64_t i = 0; i < m; i++) add((vid)r.below(n), (vid)r.below(n));
  } else if (shape == 2) {  // path with chords
    for (vid i = 0; i + 1 < n; i++) add(i, i + 1);
    for (vid i = 0; i < n / 3; i++) add((vid)r.below(n), (vid)r.below(n));
  } else if (shape == 3) {  // star out of vertex 0 plus random arcs
    for (vid i = 1; i < n; i++) add(0, i);
    for (vid i = 0; i < n; i++) add((vid)r.below(n), (vid)r.below(n));
  } else if (shape == 4) {  // undirected tree plus a few arcs
    for (vid i = 1; i < n; i++) {
      vid p = (vid)r.below(i);
      Len w = r.below(wr);
      es.push_back({p, i, w});
      es.push_back({i, p, w});
    }
    for (vid i = 0; i < n / 4; i++) add((vid)r.below(n), (vid)r.below(n));
  } else {  // grid-like
    vid side = 1;
    while ((side + 1) * (side + 1) <= n) side++;
    n = side * side;
    for (vid y = 0; y < side; y++)
      for (vid x = 0; x < side; x++) {
        vid a = y * side + x;
        if (x + 1 < side) {
          add(a, a + 1);
          add(a + 1, a);
        }
        if (y + 1 < side) {
          add(a, a + side);
          add(a + side, a);
        }
      }
  }
  return Graph::fromEdges(n, es);
}

struct Tally {
  uint64_t runs = 0, wrong = 0, stale_erase = 0, stale_pull = 0, dup_u = 0, runs_with_stale = 0, runs_with_dup = 0;
  uint64_t bound_viol = 0, runs_with_bound = 0, full_nonempty = 0, minN_stale = ~0ULL, minN_dup = ~0ULL,
           minN_bound = ~0ULL, stale_parent = 0;
};

// Smallest instance (by vertices, then arcs) on which a phenomenon was observed.
struct Witness {
  bool have = false;
  Graph g;
  vid s = 0;
  void offer(const Graph& h, vid src) {
    if (have && (h.n > g.n || (h.n == g.n && h.m >= g.m))) return;
    have = true;
    g = h;
    s = src;
  }
  void write(const std::string& path, const std::string& what) const {
    if (!have) return;
    FILE* f = std::fopen(path.c_str(), "w");
    if (!f) return;
    std::fprintf(f, "# %s\n# n m source, then one arc per line: tail head weight\n%u %llu %u\n", what.c_str(), g.n,
                 (unsigned long long)g.m, s);
    for (vid u = 0; u < g.n; u++)
      for (uint32_t e = g.head[u]; e < g.head[u + 1]; e++)
        std::fprintf(f, "%u %u %llu\n", u, g.to[e], (unsigned long long)g.w[e]);
    std::fclose(f);
  }
};

static int cmdValidate(const Args& A) {
  uint64_t trials = (uint64_t)A.num("trials", 2000), seed = (uint64_t)A.num("seed", 1);
  vid maxn = (vid)A.num("maxn", 60);
  Rng r(seed);
  std::map<std::string, Tally> tally;
  Witness wDup25, wDup26, wStale26, wBound26;
  auto check = [&](const std::string& name, Solver& S, vid src, const std::vector<Len>& want,
                   const std::vector<vid>* rep) {
    C.clear();
    S.run(src);
    Tally& T = tally[name];
    T.runs++;
    bool ok = true;
    for (vid v = 0; v < (vid)want.size() && ok; v++) ok = S.dist(rep ? (*rep)[v] : v) == want[v];
    T.wrong += !ok;
    T.stale_erase += C.stale_erase;
    T.stale_pull += C.stale_pull;
    T.dup_u += C.dup_u;
    T.runs_with_stale += C.stale_erase > 0;
    T.runs_with_dup += C.dup_u > 0;
    T.bound_viol += C.bound_viol;
    T.full_nonempty += C.full_nonempty;
    T.stale_parent += C.stale_parent;
    T.runs_with_bound += C.bound_viol > 0;
    uint64_t nn = want.size();
    if (C.stale_erase > 0) T.minN_stale = std::min(T.minN_stale, nn);
    if (C.dup_u > 0) T.minN_dup = std::min(T.minN_dup, nn);
    if (C.bound_viol > 0) T.minN_bound = std::min(T.minN_bound, nn);
    return ok;
  };
  for (uint64_t it = 0; it < trials; it++) {
    Graph g = randomSmall(r, maxn);
    vid s = (vid)r.below(g.n);
    std::vector<Len> want = bellmanFord(g, s);
    {
      SolverDijkstra d(g);
      check("dijkstra", d, s, want, nullptr);
      SolverDijkstraDK dk(g);
      check("dijkstra-dk", dk, s, want, nullptr);
    }
    // paper parameters and random parameters (small k, t force partial executions)
    int rk = 1 + (int)r.below(4), rt = 1 + (int)r.below(4);
    for (int mode = 0; mode < 2; mode++) {
      int k = mode ? rk : 0, t = mode ? rt : 0;
      std::string sfx = mode ? "/random-params" : "/paper-params";
      {
        Solver25 a(g, k, t, true);
        check("bmssp25/fix" + sfx, a, s, want, nullptr);
        Solver25 b(g, k, t, false);
        check("bmssp25/literal" + sfx, b, s, want, nullptr);
        if (!mode && C.dup_u) wDup25.offer(g, s);
      }
      {
        int t6 = mode ? 2 + (int)r.below(4) : 0, k6 = mode ? 1 + (int)r.below(4) : 0;
        Solver26 a(g, t6, k6, 3.0, true, true);
        check("bmssp26/fix+basefix" + sfx, a, s, want, nullptr);
        Solver26 b(g, t6, k6, 3.0, false, true);
        check("bmssp26/basefix-only" + sfx, b, s, want, nullptr);
        Solver26 c(g, t6, k6, 3.0, true, false);
        check("bmssp26/fix-only" + sfx, c, s, want, nullptr);
        Solver26 d(g, t6, k6, 3.0, false, false);
        check("bmssp26/literal" + sfx, d, s, want, nullptr);
        if (!mode && C.dup_u) wDup26.offer(g, s);
        if (!mode && C.stale_parent) wStale26.offer(g, s);
        if (!mode && C.bound_viol) wBound26.offer(g, s);
      }
    }
    // degree-reduced inputs
    {
      unsigned delta = 3 + (unsigned)r.below(3);
      std::vector<vid> rep;
      Graph h = reduceDegree(g, delta, rep);
      Solver25 a(h, 0, 0, true);
      check("bmssp25/fix/reduced", a, rep[s], want, &rep);
      Solver26 b(h, 0, 0, (double)delta, true, true);
      check("bmssp26/fix+basefix/reduced", b, rep[s], want, &rep);
      SolverDijkstra d(h);
      check("dijkstra/reduced", d, rep[s], want, &rep);
    }
  }
  if (A.has("dump")) {
    std::string dir = A.str("dump", ".");
    wDup25.write(dir + "/witness_2025_duplicate_U.txt",
                 "2025 pseudocode without deletion: a vertex is reported complete twice in one call");
    wDup26.write(dir + "/witness_2026_duplicate_U.txt",
                 "2026 pseudocode without deletion: a vertex is reported complete twice in one call");
    wStale26.write(dir + "/witness_2026_stale_key.txt",
                   "2026 pseudocode: a vertex of U_i is still stored in the caller's D when the sub-call returns");
    wBound26.write(dir + "/witness_2026_base_case_bound.txt",
                   "2026 base case as printed: returned B' exceeds the smallest value left in D");
  }
  std::printf("variant,runs,wrong,runs_with_stale_keys,stale_keys,stale_in_parent,stale_pulls,runs_with_duplicates,"
              "duplicates,runs_with_bound_violation,bound_violations,full_with_nonempty_D,min_n_stale,min_n_dup,"
              "min_n_bound\n");
  uint64_t wrongFixed = 0;
  auto mn = [](uint64_t x) { return x == ~0ULL ? 0ULL : (unsigned long long)x; };
  for (auto& kv : tally) {
    const Tally& T = kv.second;
    std::printf("%s,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu\n", kv.first.c_str(),
                (unsigned long long)T.runs, (unsigned long long)T.wrong, (unsigned long long)T.runs_with_stale,
                (unsigned long long)T.stale_erase, (unsigned long long)T.stale_parent,
                (unsigned long long)T.stale_pull,
                (unsigned long long)T.runs_with_dup, (unsigned long long)T.dup_u,
                (unsigned long long)T.runs_with_bound, (unsigned long long)T.bound_viol,
                (unsigned long long)T.full_nonempty, mn(T.minN_stale), mn(T.minN_dup), mn(T.minN_bound));
    if (kv.first.find("literal") == std::string::npos && kv.first.find("-only") == std::string::npos)
      wrongFixed += T.wrong;
  }
  return wrongFixed ? 2 : 0;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: sssp info|run|validate ... (see the header of main.cpp)\n");
    return 1;
  }
  std::string cmd = argv[1];
  Args A(argc, argv, 2);
  try {
    if (cmd == "info") return cmdInfo(A);
    if (cmd == "run") return cmdRun(A);
    if (cmd == "validate") return cmdValidate(A);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  std::fprintf(stderr, "unknown command %s\n", cmd.c_str());
  return 1;
}
