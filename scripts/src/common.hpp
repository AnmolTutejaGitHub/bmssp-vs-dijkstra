// Shared types for the SSSP study: graph (CSR), tie-broken distance labels,
// operation counters and a worst-case linear selection routine.
#pragma once
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

namespace sp {

using Len = uint64_t;  // path length; edge weights are non-negative integers
using vid = uint32_t;  // vertex id
constexpr Len INF_LEN = ~0ULL;
constexpr vid NONE = 0xffffffffu;

// ---------------------------------------------------------------- counters
// Machine-independent cost measures. They are compiled in only with
// -DCOUNT_OPS so that the timing binary carries no instrumentation.
struct Counters {
  // comparison-addition model
  uint64_t cmp = 0;   // comparisons between path lengths / labels
  uint64_t add = 0;   // additions d[u] + w(u,v)
  // algorithmic work
  uint64_t scan = 0;        // edge scans (an edge looked at in any loop)
  uint64_t relaxed = 0;     // scans that changed or confirmed a label
  uint64_t heap_push = 0;   // priority-queue insertions (Dijkstra, local searches, base cases)
  uint64_t heap_pop = 0;
  uint64_t ds_insert = 0;   // Insert on the partial-sorting structure
  uint64_t ds_batch = 0;    // items passed to BatchPrepend (2025) / Merge (2026)
  uint64_t ds_pulled = 0;   // keys returned by Pull
  uint64_t ds_pull = 0;     // Pull calls
  uint64_t ds_split = 0;    // block splits
  uint64_t sel_elems = 0;   // elements handed to linear selection
  // where the edge scans happen
  uint64_t scan_pivot = 0;  // inside FindPivots
  uint64_t scan_base = 0;   // inside base cases
  uint64_t scan_main = 0;   // relaxation loops of BMSSP proper
  // recursion shape
  uint64_t calls = 0;       // BMSSP calls with l > 0
  uint64_t base_calls = 0;
  uint64_t partial = 0;     // calls returning B' < B
  uint64_t pivots = 0;      // keys inserted as pivots
  uint64_t frontier = 0;    // sum of |S| over calls with l > 0
  // diagnostics for the stale-key defect
  uint64_t stale_erase = 0; // completed vertices that were still stored in D
  uint64_t stale_parent = 0;  // 2026: vertices of U_i present in the caller's D before Merge
  uint64_t stale_pull = 0;  // already-completed vertices returned by Pull
  uint64_t dup_u = 0;       // vertices reported complete twice inside one call
  uint64_t u_total = 0;     // total size of all returned sets U at the top level
  uint64_t peak = 0;        // largest number of entries held by Dijkstra's heap
  // diagnostics for the returned boundary (2026 paper, Lemma 3.7(2) and Lemma 3.1)
  uint64_t bound_viol = 0;  // calls returning B' larger than the smallest value left in D
  uint64_t full_nonempty = 0;  // calls returning B' = B together with a non-empty D
  void clear() { *this = Counters(); }
};
extern Counters C;

#ifdef COUNT_OPS
#define CNT(x) (++::sp::C.x)
#define CNTN(x, n) (::sp::C.x += (uint64_t)(n))
#else
#define CNT(x) ((void)0)
#define CNTN(x, n) ((void)0)
#endif

// ------------------------------------------------------------------- graph
struct Edge {
  vid u, v;
  Len w;
};

struct Graph {
  vid n = 0;
  uint64_t m = 0;
  std::vector<uint32_t> head;  // size n+1
  std::vector<vid> to;
  std::vector<Len> w;

  static Graph fromEdges(vid n, const std::vector<Edge>& es) {
    Graph g;
    g.n = n;
    g.m = es.size();
    g.head.assign((size_t)n + 1, 0);
    for (const Edge& e : es) g.head[e.u + 1]++;
    for (vid i = 0; i < n; i++) g.head[i + 1] += g.head[i];
    g.to.resize(g.m);
    g.w.resize(g.m);
    std::vector<uint32_t> at(g.head.begin(), g.head.end() - 1);
    for (const Edge& e : es) {
      uint32_t p = at[e.u]++;
      g.to[p] = e.v;
      g.w[p] = e.w;
    }
    return g;
  }
  uint32_t outdeg(vid u) const { return head[u + 1] - head[u]; }
};

// ------------------------------------------------------------------ labels
// A distance label is the tuple <length, #edges, vertex, predecessor>, compared
// lexicographically. This is the tie-breaking rule of Duan et al. that makes
// all path lengths distinct (2025 paper, Assumption 2.1; 2026 paper, Section 2.3).
struct Label {
  Len len;
  uint32_t hops;
  vid v;
  vid pred;
};
constexpr Label INF_LABEL = {INF_LEN, 0xffffffffu, NONE, NONE};

inline bool lt(const Label& a, const Label& b) {
  CNT(cmp);
  if (a.len != b.len) return a.len < b.len;
  if (a.hops != b.hops) return a.hops < b.hops;
  if (a.v != b.v) return a.v < b.v;
  return a.pred < b.pred;
}
inline bool rawLt(const Label& a, const Label& b) {  // uncounted, for diagnostics only
  if (a.len != b.len) return a.len < b.len;
  if (a.hops != b.hops) return a.hops < b.hops;
  if (a.v != b.v) return a.v < b.v;
  return a.pred < b.pred;
}
inline bool same(const Label& a, const Label& b) {
  return a.len == b.len && a.hops == b.hops && a.v == b.v && a.pred == b.pred;
}
struct LabelLess {
  bool operator()(const Label& a, const Label& b) const { return lt(a, b); }
};

// Per-vertex part of a label (the vertex id is implicit).
struct VLabel {
  Len len;
  uint32_t hops;
  vid pred;
};

struct LabelArray {
  std::vector<VLabel> d;
  void reset(vid n, vid s) {
    d.assign(n, VLabel{INF_LEN, 0xffffffffu, NONE});
    d[s] = VLabel{0, 0, NONE};
  }
  Label get(vid v) const { return Label{d[v].len, d[v].hops, v, d[v].pred}; }
  // Candidate label of v through edge (u,v) of weight w: one addition.
  Label through(vid u, vid v, Len w) const {
    CNT(add);
    return Label{d[u].len + w, d[u].hops + 1, v, u};
  }
  // cand <= d[v] ?  (one comparison; both labels belong to the same vertex)
  bool notWorse(const Label& c) const {
    CNT(cmp);
    const VLabel& x = d[c.v];
    if (c.len != x.len) return c.len < x.len;
    if (c.hops != x.hops) return c.hops < x.hops;
    return c.pred <= x.pred;
  }
  void set(const Label& c) { d[c.v] = VLabel{c.len, c.hops, c.pred}; }
};

// --------------------------------------------------------------- selection
// select_kth rearranges a[0..n) so that a[k] is the element of rank k, every
// element before it is smaller and every element after it is larger.
// Quickselect with median-of-three pivots; once the work spent exceeds a
// linear budget the pivot rule switches to median-of-medians (Blum et al.),
// so the running time is O(n) in the worst case (introselect, Musser 1997).
// All keys handled here are distinct labels.
template <class T, class Less>
void insertion_sort(T* a, size_t n, Less less) {
  for (size_t i = 1; i < n; i++) {
    T x = a[i];
    size_t j = i;
    while (j > 0 && less(x, a[j - 1])) {
      a[j] = a[j - 1];
      j--;
    }
    a[j] = x;
  }
}

template <class T, class Less>
void select_kth(T* a, size_t n, size_t k, Less less);

template <class T, class Less>
size_t mom_pivot(T* a, size_t lo, size_t hi, Less less) {
  size_t g = 0;
  for (size_t i = lo; i < hi; i += 5) {
    size_t len = std::min<size_t>(5, hi - i);
    insertion_sort(a + i, len, less);
    std::swap(a[lo + g], a[i + len / 2]);
    g++;
  }
  select_kth(a + lo, g, g / 2, less);
  return lo + g / 2;
}

template <class T, class Less>
void select_kth(T* a, size_t n, size_t k, Less less) {
  CNTN(sel_elems, n);
  size_t lo = 0, hi = n;
  size_t budget = 4 * n + 64, spent = 0;
  bool guarded = false;
  while (hi - lo > 12) {
    size_t len = hi - lo, p;
    if (!guarded) {
      size_t x = lo, y = lo + len / 2, z = hi - 1;
      if (less(a[y], a[x])) std::swap(x, y);
      if (less(a[z], a[y])) {
        std::swap(y, z);
        if (less(a[y], a[x])) std::swap(x, y);
      }
      p = y;
    } else {
      p = mom_pivot(a, lo, hi, less);
    }
    std::swap(a[p], a[hi - 1]);
    size_t store = lo;
    for (size_t i = lo; i + 1 < hi; i++)
      if (less(a[i], a[hi - 1])) {
        std::swap(a[i], a[store]);
        store++;
      }
    std::swap(a[store], a[hi - 1]);
    if (k == store) return;
    if (k < store)
      hi = store;
    else
      lo = store + 1;
    spent += len;
    if (spent > budget) guarded = true;
  }
  insertion_sort(a + lo, hi - lo, less);
}

// ------------------------------------------------------------------- misc
inline uint64_t now_ns() {
  return (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

// 2^e capped so that products with small factors never overflow.
inline uint64_t pow2cap(uint64_t e) { return e >= 40 ? (1ULL << 40) : (1ULL << e); }

struct Rng {  // splitmix64
  uint64_t s;
  explicit Rng(uint64_t seed) : s(seed) {}
  uint64_t next() {
    uint64_t z = (s += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }
  uint64_t below(uint64_t n) { return next() % n; }
};

}  // namespace sp
