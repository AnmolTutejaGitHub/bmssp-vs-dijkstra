// Graph sources: the SSTD 2005 road networks, synthetic generators and the
// degree-reduction transformation used by both Duan et al. papers.
#pragma once
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "common.hpp"

namespace sp {

// Road network edge file of the SSTD 2005 collection:
//   "<edge id> <start node> <end node> <L2 distance>"   (one road segment per line)
// A segment is an undirected road, so it yields two opposite arcs. Distances
// carry six decimals and are scaled by 10^6 to exact integers.
inline Graph loadSSTD(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  std::vector<Edge> es;
  std::string line;
  vid n = 0;
  while (std::getline(in, line)) {
    std::istringstream ls(line);
    long long id, u, v;
    double d;
    if (!(ls >> id >> u >> v >> d)) continue;
    if (u < 0 || v < 0 || d < 0) throw std::runtime_error("bad line in " + path);
    Len w = (Len)std::llround(d * 1e6);
    es.push_back({(vid)u, (vid)v, w});
    es.push_back({(vid)v, (vid)u, w});
    n = std::max<vid>(n, (vid)std::max(u, v) + 1);
  }
  return Graph::fromEdges(n, es);
}

constexpr Len WMAX_DEFAULT = 1u << 20;  // synthetic weights are uniform in [1, 2^20]

// Edge list of the SNAP collection: lines "<from> <to>", comments start with '#'.
// The files carry no weights. mode "unit": every arc has weight 1; mode "rand":
// arc i gets a weight uniform in [1, 2^20] derived from a fixed seed and i.
// The arcs are used exactly as listed (the roadNet files already contain both
// directions of every road); vertices are renumbered consecutively.
// A ".gz" file is read through `gunzip -c`.
inline Graph loadSNAP(const std::string& path, const std::string& mode) {
  if (mode != "unit" && mode != "rand") throw std::runtime_error("snap weight mode must be unit or rand");
  bool gz = path.size() > 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
  if (path.find('\'') != std::string::npos) throw std::runtime_error("bad path");
  FILE* f = gz ? popen(("gunzip -c '" + path + "'").c_str(), "r") : std::fopen(path.c_str(), "r");
  if (!f) throw std::runtime_error("cannot open " + path);
  std::vector<Edge> es;
  vid n = 0;
  char line[256];
  Rng r(0x5eed5eedULL);
  while (std::fgets(line, sizeof line, f)) {
    if (line[0] == '#' || line[0] == '\n') continue;
    char* end;
    unsigned long u = std::strtoul(line, &end, 10);
    unsigned long v = std::strtoul(end, &end, 10);
    Len w = mode == "unit" ? 1 : 1 + r.below(WMAX_DEFAULT);
    es.push_back({(vid)u, (vid)v, w});
    n = std::max<vid>(n, (vid)std::max(u, v) + 1);
  }
  if (gz) pclose(f); else std::fclose(f);
  if (es.empty()) throw std::runtime_error("no arcs read from " + path);
  // The ids in the files have gaps; renumber the vertices that occur, in id order.
  std::vector<vid> id(n, NONE);
  for (const Edge& e : es) id[e.u] = id[e.v] = 0;
  vid used = 0;
  for (vid i = 0; i < n; i++)
    if (id[i] != NONE) id[i] = used++;
  for (Edge& e : es) {
    e.u = id[e.u];
    e.v = id[e.v];
  }
  return Graph::fromEdges(used, es);
}

// Random sparse digraph: a random Hamiltonian cycle (strong connectivity)
// plus (deg-1)*n uniformly random arcs. Weights are uniform in [1, WMAX].
inline Graph genRandom(vid n, unsigned deg, uint64_t seed, Len WMAX = WMAX_DEFAULT) {
  Rng r(seed);
  std::vector<vid> perm(n);
  for (vid i = 0; i < n; i++) perm[i] = i;
  for (vid i = n; i > 1; i--) std::swap(perm[i - 1], perm[r.below(i)]);
  std::vector<Edge> es;
  es.reserve((size_t)n * deg);
  for (vid i = 0; i < n; i++) es.push_back({perm[i], perm[(i + 1) % n], 1 + r.below(WMAX)});
  for (uint64_t i = 0; i < (uint64_t)n * (deg - 1); i++) {
    vid u = (vid)r.below(n), v = (vid)r.below(n);
    es.push_back({u, v, 1 + r.below(WMAX)});
  }
  return Graph::fromEdges(n, es);
}

// side x side grid, 4-neighbourhood, both directions, independent weights.
inline Graph genGrid(vid side, uint64_t seed) {
  const Len WMAX = WMAX_DEFAULT;
  Rng r(seed);
  vid n = side * side;
  std::vector<Edge> es;
  es.reserve((size_t)n * 4);
  for (vid y = 0; y < side; y++)
    for (vid x = 0; x < side; x++) {
      vid a = y * side + x;
      if (x + 1 < side) {
        es.push_back({a, a + 1, 1 + r.below(WMAX)});
        es.push_back({a + 1, a, 1 + r.below(WMAX)});
      }
      if (y + 1 < side) {
        es.push_back({a, a + side, 1 + r.below(WMAX)});
        es.push_back({a + side, a, 1 + r.below(WMAX)});
      }
    }
  return Graph::fromEdges(n, es);
}

// Preferential attachment (Barabasi-Albert): every new vertex attaches to d
// earlier vertices chosen proportionally to degree; arcs in both directions.
// Produces a heavy-tailed degree distribution.
inline Graph genBA(vid n, unsigned d, uint64_t seed) {
  const Len WMAX = WMAX_DEFAULT;
  Rng r(seed);
  std::vector<Edge> es;
  std::vector<vid> ends;  // every arc endpoint once: sampling it is degree-proportional
  es.reserve((size_t)n * d * 2);
  ends.reserve((size_t)n * d * 2);
  for (vid v = 1; v < n; v++) {
    unsigned links = std::min<unsigned>(d, v);
    for (unsigned j = 0; j < links; j++) {
      vid u = ends.empty() ? 0 : ends[r.below(ends.size())];
      if (j == 0 && v <= d) u = v - 1;  // keep the seed part connected
      Len w = 1 + r.below(WMAX);
      es.push_back({v, u, w});
      es.push_back({u, v, w});
      ends.push_back(u);
      ends.push_back(v);
    }
  }
  return Graph::fromEdges(n, es);
}

// Wide-frontier instance: a random sparse digraph in which vertex 0 also has an
// arc to every other vertex with a weight above every ordinary edge weight.
// A search from vertex 0 starts with all n vertices in its frontier.
inline Graph genHub(vid n, unsigned deg, uint64_t seed) {
  const Len WMAX = WMAX_DEFAULT;
  Rng r(seed);
  std::vector<vid> perm(n);
  for (vid i = 0; i < n; i++) perm[i] = i;
  for (vid i = n; i > 1; i--) std::swap(perm[i - 1], perm[r.below(i)]);
  std::vector<Edge> es;
  es.reserve((size_t)n * (deg + 1));
  for (vid i = 0; i < n; i++) es.push_back({perm[i], perm[(i + 1) % n], 1 + r.below(WMAX)});
  for (uint64_t i = 0; i < (uint64_t)n * (deg - 1); i++) {
    vid u = (vid)r.below(n), v = (vid)r.below(n);
    es.push_back({u, v, 1 + r.below(WMAX)});
  }
  for (vid v = 1; v < n; v++) es.push_back({0, v, WMAX + 1 + r.below(WMAX)});
  return Graph::fromEdges(n, es);
}

// Degree reduction (Frederickson-style, as described by Duan et al.): vertex v
// with Delta_v incident arcs becomes a zero-weight directed cycle of
// ceil(Delta_v / (delta-2)) vertices, each carrying at most delta-2 arc
// endpoints, so every vertex of the result has in-degree + out-degree <= delta.
// delta = 3 gives the constant-degree graph of the 2025 paper (in/out-degree <= 2).
// rep[v] is a vertex of the cycle of v; its distance equals the distance of v.
inline Graph reduceDegree(const Graph& g, unsigned delta, std::vector<vid>& rep) {
  if (delta < 3) throw std::runtime_error("delta must be at least 3");
  const unsigned slots = delta - 2;
  std::vector<uint32_t> deg(g.n, 0);
  for (vid u = 0; u < g.n; u++)
    for (uint32_t e = g.head[u]; e < g.head[u + 1]; e++) {
      deg[u]++;
      deg[g.to[e]]++;
    }
  rep.assign(g.n, 0);
  std::vector<uint32_t> csize(g.n);
  uint64_t total = 0;
  for (vid v = 0; v < g.n; v++) {
    csize[v] = std::max<uint32_t>(1, (deg[v] + slots - 1) / slots);
    rep[v] = (vid)total;
    total += csize[v];
  }
  if (total >= NONE) throw std::runtime_error("reduced graph too large");
  std::vector<Edge> es;
  es.reserve(g.m + total);
  for (vid v = 0; v < g.n; v++)
    if (csize[v] >= 2)
      for (uint32_t i = 0; i < csize[v]; i++) es.push_back({rep[v] + i, rep[v] + (i + 1) % csize[v], 0});
  std::vector<uint32_t> used(g.n, 0);
  for (vid u = 0; u < g.n; u++)
    for (uint32_t e = g.head[u]; e < g.head[u + 1]; e++) {
      vid v = g.to[e];
      vid a = rep[u] + used[u]++ / slots;
      vid b = rep[v] + used[v]++ / slots;
      es.push_back({a, b, g.w[e]});
    }
  return Graph::fromEdges((vid)total, es);
}

// "sstd:<edge file>", "snap:<edge list>:<unit|rand>", "rand:<n>:<deg>:<seed>[:<wmax>]", "grid:<side>:<seed>",
// "ba:<n>:<d>:<seed>", "hub:<n>:<deg>:<seed>"
inline Graph makeGraph(const std::string& spec) {
  std::vector<std::string> p;
  size_t c = spec.find(':');
  if (c == std::string::npos) throw std::runtime_error("bad graph spec " + spec);
  std::string kind = spec.substr(0, c);
  if (kind == "sstd") return loadSSTD(spec.substr(c + 1));
  if (kind == "snap") {  // snap:<file>:<unit|rand>
    std::string rest = spec.substr(c + 1);
    size_t c2 = rest.rfind(':');
    if (c2 == std::string::npos) throw std::runtime_error("use snap:<file>:<unit|rand>");
    return loadSNAP(rest.substr(0, c2), rest.substr(c2 + 1));
  }
  std::stringstream ss(spec.substr(c + 1));
  std::string tok;
  while (std::getline(ss, tok, ':')) p.push_back(tok);
  auto num = [&](size_t i) { return (uint64_t)std::stoull(p.at(i)); };
  if (kind == "rand")
    return genRandom((vid)num(0), (unsigned)num(1), num(2), p.size() > 3 ? (Len)num(3) : WMAX_DEFAULT);
  if (kind == "grid") return genGrid((vid)num(0), num(1));
  if (kind == "ba") return genBA((vid)num(0), (unsigned)num(1), num(2));
  if (kind == "hub") return genHub((vid)num(0), (unsigned)num(1), num(2));
  throw std::runtime_error("unknown graph kind " + kind);
}

}  // namespace sp
