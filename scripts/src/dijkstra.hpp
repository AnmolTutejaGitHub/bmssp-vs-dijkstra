// Dijkstra's algorithm, two binary-heap variants.
//   DijkstraLazy  : heap of (key, vertex) pairs with lazy deletion (the usual
//                   practical implementation; the heap may hold up to m entries)
//   DijkstraIndexed: indexed heap with DecreaseKey (at most n entries)
#pragma once
#include "common.hpp"

namespace sp {

struct DijkstraLazy {
  const Graph& g;
  std::vector<Len> dist;
  struct Item {
    Len key;
    vid v;
  };
  std::vector<Item> h;

  explicit DijkstraLazy(const Graph& g_) : g(g_) { h.reserve(g.n); }

  static bool less(const Item& a, const Item& b) {
    CNT(cmp);
    return a.key < b.key;
  }
  void push(Item x) {
    CNT(heap_push);
    size_t i = h.size();
    h.push_back(x);
#ifdef COUNT_OPS
    if (h.size() > C.peak) C.peak = h.size();
#endif
    while (i > 0) {
      size_t p = (i - 1) / 2;
      if (!less(x, h[p])) break;
      h[i] = h[p];
      i = p;
    }
    h[i] = x;
  }
  Item pop() {
    CNT(heap_pop);
    Item top = h[0], x = h.back();
    h.pop_back();
    size_t n = h.size(), i = 0;
    if (n) {
      for (;;) {
        size_t c = 2 * i + 1;
        if (c >= n) break;
        if (c + 1 < n && less(h[c + 1], h[c])) c++;
        if (!less(h[c], x)) break;
        h[i] = h[c];
        i = c;
      }
      h[i] = x;
    }
    return top;
  }

  void run(vid s) {
    dist.assign(g.n, INF_LEN);
    h.clear();
    dist[s] = 0;
    push({0, s});
    while (!h.empty()) {
      Item it = pop();
      vid u = it.v;
      CNT(cmp);
      if (it.key != dist[u]) continue;  // stale entry
      for (uint32_t e = g.head[u]; e < g.head[u + 1]; e++) {
        CNT(scan);
        CNT(add);
        Len nd = it.key + g.w[e];
        vid v = g.to[e];
        CNT(cmp);
        if (nd < dist[v]) {
          CNT(relaxed);
          dist[v] = nd;
          push({nd, v});
        }
      }
    }
  }
};

struct DijkstraIndexed {
  const Graph& g;
  std::vector<Len> dist;
  std::vector<vid> h;         // heap of vertices keyed by dist
  std::vector<uint32_t> pos;  // position in heap, NONE if absent

  explicit DijkstraIndexed(const Graph& g_) : g(g_) { h.reserve(g.n); }

  bool less(vid a, vid b) const {
    CNT(cmp);
    return dist[a] < dist[b];
  }
  void up(size_t i) {
    vid x = h[i];
    while (i > 0) {
      size_t p = (i - 1) / 2;
      if (!less(x, h[p])) break;
      h[i] = h[p];
      pos[h[i]] = (uint32_t)i;
      i = p;
    }
    h[i] = x;
    pos[x] = (uint32_t)i;
  }
  vid pop() {
    CNT(heap_pop);
    vid top = h[0], x = h.back();
    h.pop_back();
    pos[top] = NONE;
    size_t n = h.size(), i = 0;
    if (n) {
      for (;;) {
        size_t c = 2 * i + 1;
        if (c >= n) break;
        if (c + 1 < n && less(h[c + 1], h[c])) c++;
        if (!less(h[c], x)) break;
        h[i] = h[c];
        pos[h[i]] = (uint32_t)i;
        i = c;
      }
      h[i] = x;
      pos[x] = (uint32_t)i;
    }
    return top;
  }

  void run(vid s) {
    dist.assign(g.n, INF_LEN);
    pos.assign(g.n, NONE);
    h.clear();
    dist[s] = 0;
    h.push_back(s);
    pos[s] = 0;
    while (!h.empty()) {
      vid u = pop();
      Len du = dist[u];
      for (uint32_t e = g.head[u]; e < g.head[u + 1]; e++) {
        CNT(scan);
        CNT(add);
        Len nd = du + g.w[e];
        vid v = g.to[e];
        CNT(cmp);
        if (nd < dist[v]) {
          CNT(relaxed);
          dist[v] = nd;
          if (pos[v] == NONE) {
            CNT(heap_push);
            h.push_back(v);
#ifdef COUNT_OPS
            if (h.size() > C.peak) C.peak = h.size();
#endif
            up(h.size() - 1);
          } else {
            up(pos[v]);
          }
        }
      }
    }
  }
};

}  // namespace sp
