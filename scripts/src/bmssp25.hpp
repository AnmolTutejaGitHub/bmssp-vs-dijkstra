// The O(m log^{2/3} n) algorithm of Duan, Mao, Mao, Shu and Yin (STOC 2025),
// written from the pseudocode of the paper:
//   Algorithm 1  FindPivots   (k rounds of Bellman-Ford style relaxation)
//   Algorithm 2  BaseCase     (bounded Dijkstra that settles k+1 vertices)
//   Algorithm 3  BMSSP        (recursion over the partial-sorting structure)
//   Lemma 3.3    block structure with Insert / BatchPrepend / Pull
// (numbering of the full version, arXiv:2504.17033v2)
// One switch deviates from the paper on purpose: with `fix` set, a vertex is
// deleted from D as soon as a recursive call reports it complete (the repair
// of Castro, Clementino and de Freitas for the stale-key defect).
#pragma once
#include <cmath>
#include <map>
#include "common.hpp"

namespace sp {

// Binary heap of labels with lazy deletion, used by the bounded local searches.
struct LabelHeap {
  std::vector<Label> h;
  bool empty() const { return h.empty(); }
  void clear() { h.clear(); }
  void push(const Label& x) {
    CNT(heap_push);
    size_t i = h.size();
    h.push_back(x);
    while (i > 0) {
      size_t p = (i - 1) / 2;
      if (!lt(x, h[p])) break;
      h[i] = h[p];
      i = p;
    }
    h[i] = x;
  }
  Label pop() {
    CNT(heap_pop);
    Label top = h[0], x = h.back();
    h.pop_back();
    size_t n = h.size(), i = 0;
    if (n) {
      for (;;) {
        size_t c = 2 * i + 1;
        if (c >= n) break;
        if (c + 1 < n && lt(h[c + 1], h[c])) c++;
        if (!lt(h[c], x)) break;
        h[i] = h[c];
        i = c;
      }
      h[i] = x;
    }
    return top;
  }
};

// ------------------------------------------------------------ Lemma 3.3
// Two sequences of blocks. D1 receives single insertions; its blocks hold at
// most M pairs and are indexed by a red-black tree over their upper bounds.
// D0 receives batch prepends. A direct-address table gives O(1) lookup and
// deletion of a key (one table per recursion level, as in Castro et al.).
class BlockSeq25 {
  struct Block {
    std::vector<Label> items;
    Label ub;
    bool inD1;
  };
  static constexpr int ABSENT = -1, TMP = -2;
  size_t M = 1;
  Label B = INF_LABEL;
  size_t cnt = 0;
  std::vector<Block> pool;
  std::vector<int> freeIds;
  std::map<Label, int, LabelLess> d1;
  std::vector<int> d0;  // front of the sequence is d0.back()
  std::vector<int> blk;
  std::vector<uint32_t> pos;
  std::vector<Label> scratch, tmp;
  std::vector<std::pair<size_t, size_t>> ranges;

  int newBlock(bool inD1, const Label& ub) {
    int b;
    if (!freeIds.empty()) {
      b = freeIds.back();
      freeIds.pop_back();
    } else {
      b = (int)pool.size();
      pool.emplace_back();
    }
    pool[b].items.clear();
    pool[b].inD1 = inD1;
    pool[b].ub = ub;
    return b;
  }
  void release(int b) {
    pool[b].items.clear();
    freeIds.push_back(b);
  }
  void place(int b, const Label& x) {
    blk[x.v] = b;
    pos[x.v] = (uint32_t)pool[b].items.size();
    pool[b].items.push_back(x);
  }
  void reindex(int b) {
    std::vector<Label>& it = pool[b].items;
    for (uint32_t i = 0; i < it.size(); i++) {
      blk[it[i].v] = b;
      pos[it[i].v] = i;
    }
  }
  // Split an over-full block of D1 around its median.
  void split(std::map<Label, int, LabelLess>::iterator it) {
    CNT(ds_split);
    int b = it->second;
    size_t n = pool[b].items.size(), h = n / 2;
    select_kth(pool[b].items.data(), n, h - 1, LabelLess());
    Label ub = pool[b].items[h - 1];
    int nb = newBlock(true, ub);
    pool[nb].items.assign(pool[b].items.begin(), pool[b].items.begin() + h);
    pool[b].items.erase(pool[b].items.begin(), pool[b].items.begin() + h);
    reindex(nb);
    reindex(b);
    d1.emplace_hint(it, ub, nb);
  }

 public:
  void resize(vid n) {
    blk.assign(n, ABSENT);
    pos.assign(n, 0);
  }
  bool empty() const { return cnt == 0; }
  size_t size() const { return cnt; }
  bool contains(vid v) const { return blk[v] >= 0; }

  void init(size_t M_, const Label& B_) {
    assert(cnt == 0);
    M = M_;
    B = B_;
    d0.clear();
    d1.clear();
    d1.emplace(B, newBlock(true, B));
  }
  void clear() {
    for (int b : d0) {
      for (const Label& x : pool[b].items) blk[x.v] = ABSENT;
      release(b);
    }
    for (auto& kv : d1) {
      for (const Label& x : pool[kv.second].items) blk[x.v] = ABSENT;
      release(kv.second);
    }
    d0.clear();
    d1.clear();
    cnt = 0;
  }

  bool erase(vid v) {
    int b = blk[v];
    if (b < 0) return false;
    Block& bl = pool[b];
    uint32_t p = pos[v];
    Label last = bl.items.back();
    bl.items[p] = last;
    pos[last.v] = p;
    bl.items.pop_back();
    blk[v] = ABSENT;
    cnt--;
    if (bl.items.empty() && bl.inD1 && !same(bl.ub, B)) {
      d1.erase(bl.ub);
      release(b);
    }
    return true;
  }

  void insert(const Label& x) {
    CNT(ds_insert);
    int b = blk[x.v];
    if (b >= 0) {
      if (!lt(x, pool[b].items[pos[x.v]])) return;
      erase(x.v);
    }
    auto it = d1.lower_bound(x);  // block with the smallest upper bound >= x
    assert(it != d1.end());
    place(it->second, x);
    cnt++;
    if (pool[it->second].items.size() > M) split(it);
  }

  // L may contain several pairs with the same key; the smallest is kept.
  void batchPrepend(std::vector<Label>& L) {
    CNTN(ds_batch, L.size());
    tmp.clear();
    for (const Label& x : L) {
      int b = blk[x.v];
      if (b == TMP) {
        Label& y = tmp[pos[x.v]];
        if (lt(x, y)) y = x;
        continue;
      }
      if (b >= 0) {
        if (!lt(x, pool[b].items[pos[x.v]])) continue;
        erase(x.v);
      }
      blk[x.v] = TMP;
      pos[x.v] = (uint32_t)tmp.size();
      tmp.push_back(x);
    }
    if (tmp.empty()) return;
    cnt += tmp.size();
    size_t half = (M + 1) / 2;
    ranges.clear();
    if (tmp.size() <= M) {
      ranges.push_back({0, tmp.size()});
    } else {
      // repeated median splitting until every piece has at most ceil(M/2) pairs
      std::vector<std::pair<size_t, size_t>> st;
      st.push_back({0, tmp.size()});
      while (!st.empty()) {
        auto [lo, hi] = st.back();
        st.pop_back();
        if (hi - lo <= half) {
          ranges.push_back({lo, hi});
          continue;
        }
        size_t mid = lo + (hi - lo) / 2;
        select_kth(tmp.data() + lo, hi - lo, mid - lo, LabelLess());
        st.push_back({mid, hi});  // popped second: ranges come out ascending
        st.push_back({lo, mid});
      }
    }
    for (size_t r = ranges.size(); r-- > 0;) {  // largest first, smallest ends in front
      int b = newBlock(false, INF_LABEL);
      pool[b].items.assign(tmp.begin() + ranges[r].first, tmp.begin() + ranges[r].second);
      reindex(b);
      d0.push_back(b);
    }
  }

  // Returns at most M keys with the smallest values and a separator x.
  void pull(std::vector<vid>& out, Label& x) {
    CNT(ds_pull);
    out.clear();
    if (cnt <= M) {
      for (int b : d0)
        for (const Label& y : pool[b].items) out.push_back(y.v);
      for (auto& kv : d1)
        for (const Label& y : pool[kv.second].items) out.push_back(y.v);
      clear();
      d1.emplace(B, newBlock(true, B));
      x = B;
      CNTN(ds_pulled, out.size());
      return;
    }
    while (!d0.empty() && pool[d0.back()].items.empty()) {
      release(d0.back());
      d0.pop_back();
    }
    // Prefixes holding more than M pairs (or everything) from each sequence:
    // together they contain the M+1 smallest pairs of the structure.
    scratch.clear();
    size_t c = 0;
    for (size_t i = d0.size(); i-- > 0 && c <= M;) {
      const std::vector<Label>& it = pool[d0[i]].items;
      scratch.insert(scratch.end(), it.begin(), it.end());
      c += it.size();
    }
    c = 0;
    for (auto it = d1.begin(); it != d1.end() && c <= M; ++it) {
      const std::vector<Label>& its = pool[it->second].items;
      scratch.insert(scratch.end(), its.begin(), its.end());
      c += its.size();
    }
    assert(scratch.size() > M);
    select_kth(scratch.data(), scratch.size(), M, LabelLess());
    x = scratch[M];
    for (size_t i = 0; i < M; i++) {
      out.push_back(scratch[i].v);
      erase(scratch[i].v);
    }
    CNTN(ds_pulled, out.size());
  }
};

// ---------------------------------------------------------------- BMSSP
struct Bmssp25 {
  const Graph& g;
  LabelArray L;
  int k = 1, t = 1, top = 1;
  bool fix = true;  // delete completed vertices from D (Castro et al.)

  struct Res {
    Label Bp;
    size_t ubeg, uend;
  };

  std::vector<BlockSeq25> Ds;            // one structure per level
  std::vector<uint64_t> Mof, capOf;      // 2^{(l-1)t} and k 2^{lt}
  std::vector<std::vector<vid>> Sbuf, Wbuf;
  std::vector<std::vector<uint32_t>> ust;  // membership stamps of U per level
  std::vector<uint32_t> uEpoch;
  std::vector<vid> Ubuf, P, cur, nxt, path;
  std::vector<Label> K;
  std::vector<size_t> dupIdx;
  std::vector<uint32_t> wst, rst, bst;   // stamps: in W, in W_i, extracted in base case
  uint32_t wEpoch = 0, rEpoch = 0, bEpoch = 0;
  std::vector<vid> par, rootOf;
  std::vector<uint32_t> tsize;
  LabelHeap heap;

  // Default parameters of the paper: k = floor(log^{1/3} n), t = floor(log^{2/3} n).
  static void defaults(vid n, int& k, int& t) {
    double lg = std::log2((double)std::max<vid>(n, 2));
    k = std::max(1, (int)std::floor(std::cbrt(lg) + 1e-9));
    t = std::max(1, (int)std::floor(std::cbrt(lg) * std::cbrt(lg) + 1e-9));
  }

  Bmssp25(const Graph& g_, int k_ = 0, int t_ = 0) : g(g_) {
    defaults(g.n, k, t);
    if (k_ > 0) k = k_;
    if (t_ > 0) t = t_;
    double lg = std::log2((double)std::max<vid>(g.n, 2));
    top = std::max(1, (int)std::ceil(lg / t - 1e-9));
    Ds.resize(top + 1);
    Mof.assign(top + 1, 1);
    capOf.assign(top + 1, 1);
    Sbuf.resize(top + 1);
    Wbuf.resize(top + 1);
    ust.resize(top + 1);
    uEpoch.assign(top + 1, 0);
    for (int l = 1; l <= top; l++) {
      Ds[l].resize(g.n);
      Mof[l] = pow2cap((uint64_t)(l - 1) * t);
      capOf[l] = (uint64_t)k * pow2cap((uint64_t)l * t);
      ust[l].assign(g.n, 0);
    }
    wst.assign(g.n, 0);
    rst.assign(g.n, 0);
    bst.assign(g.n, 0);
    par.assign(g.n, NONE);
    rootOf.assign(g.n, NONE);
    tsize.assign(g.n, 0);
  }

  // Algorithm 1.
  void findPivots(const Label& B, const std::vector<vid>& S, std::vector<vid>& W) {
    W.clear();
    P.clear();
    ++wEpoch;
    for (vid x : S)
      if (wst[x] != wEpoch) {
        wst[x] = wEpoch;
        W.push_back(x);
      }
    cur.assign(S.begin(), S.end());
    const size_t limit = (size_t)k * S.size();
    for (int i = 1; i <= k; i++) {
      nxt.clear();
      ++rEpoch;
      for (vid u : cur) {
        for (uint32_t e = g.head[u]; e < g.head[u + 1]; e++) {
          CNT(scan);
          CNT(scan_pivot);
          vid v = g.to[e];
          Label c = L.through(u, v, g.w[e]);
          if (L.notWorse(c)) {
            CNT(relaxed);
            L.set(c);
            if (lt(c, B)) {
              if (rst[v] != rEpoch) {
                rst[v] = rEpoch;
                nxt.push_back(v);
              }
              if (wst[v] != wEpoch) {
                wst[v] = wEpoch;
                W.push_back(v);
              }
            }
          }
        }
      }
      if (W.size() > limit) {
        P.assign(S.begin(), S.end());
        return;
      }
      cur.swap(nxt);
    }
    // F = {(u,v) : u,v in W, d[v] = d[u] + w(u,v)}; with tie-breaking the tight
    // edge into v is the one from its recorded predecessor.
    for (vid v : W) {
      par[v] = NONE;
      rootOf[v] = NONE;
      tsize[v] = 0;
    }
    for (vid u : W) {
      for (uint32_t e = g.head[u]; e < g.head[u + 1]; e++) {
        CNT(scan);
        CNT(scan_pivot);
        vid v = g.to[e];
        if (wst[v] != wEpoch || L.d[v].pred != u) continue;
        CNT(add);
        CNT(cmp);
        if (L.d[u].len + g.w[e] == L.d[v].len && L.d[u].hops + 1 == L.d[v].hops) par[v] = u;
      }
    }
    for (vid v : W) {
      path.clear();
      vid x = v;
      while (rootOf[x] == NONE && par[x] != NONE) {
        path.push_back(x);
        x = par[x];
      }
      vid r = rootOf[x] == NONE ? x : rootOf[x];
      rootOf[x] = r;
      for (vid y : path) rootOf[y] = r;
      tsize[r]++;
    }
    for (vid x : S)
      if (par[x] == NONE && tsize[x] >= (uint32_t)k) P.push_back(x);
  }

  // Algorithm 2.
  Res baseCase(const Label& B, vid x) {
    CNT(base_calls);
    size_t ubeg = Ubuf.size();
    heap.clear();
    ++bEpoch;
    heap.push(L.get(x));
    size_t found = 0;
    while (!heap.empty() && found < (size_t)k + 1) {
      Label top_ = heap.pop();
      vid u = top_.v;
      if (bst[u] == bEpoch) continue;  // duplicate entry
      CNT(cmp);                        // stale-entry test, counted like Dijkstra's
      if (!same(top_, L.get(u))) continue;
      bst[u] = bEpoch;
      Ubuf.push_back(u);
      found++;
      for (uint32_t e = g.head[u]; e < g.head[u + 1]; e++) {
        CNT(scan);
        CNT(scan_base);
        vid v = g.to[e];
        Label c = L.through(u, v, g.w[e]);
        if (L.notWorse(c) && lt(c, B)) {
          CNT(relaxed);
          L.set(c);
          heap.push(c);
        }
      }
    }
    if (found <= (size_t)k) return Res{B, ubeg, Ubuf.size()};
    // k+1 vertices found: the last one extracted has the largest label
    Label Bp = L.get(Ubuf.back());
    Ubuf.pop_back();
    CNT(partial);
    return Res{Bp, ubeg, Ubuf.size()};
  }

  // Remove the entries of Ubuf listed in dupIdx (ascending positions >= from).
  void dropDuplicates(size_t from) {
    size_t w = from, d = 0;
    for (size_t i = from; i < Ubuf.size(); i++) {
      if (d < dupIdx.size() && dupIdx[d] == i) {
        d++;
        continue;
      }
      Ubuf[w++] = Ubuf[i];
    }
    Ubuf.resize(w);
    dupIdx.clear();
  }

  // Algorithm 3. The frontier S is Sbuf[l]; U is returned as a range of Ubuf.
  Res bmssp(int l, const Label& B) {
    if (l == 0) return baseCase(B, Sbuf[0][0]);
    CNT(calls);
    CNTN(frontier, Sbuf[l].size());
    findPivots(B, Sbuf[l], Wbuf[l]);
    BlockSeq25& D = Ds[l];
    D.init(Mof[l], B);
    Label Bp = B;
    for (vid x : P) {
      Label lx = L.get(x);
      CNT(pivots);
      D.insert(lx);
      if (lt(lx, Bp)) Bp = lx;
    }
    const size_t ubeg = Ubuf.size();
    const uint64_t cap = capOf[l];
    const uint32_t ue = ++uEpoch[l];
    std::vector<uint32_t>& us = ust[l];
    std::vector<vid>& Si = Sbuf[l - 1];
    while (Ubuf.size() - ubeg < cap && !D.empty()) {
      Label Bi;
      D.pull(Si, Bi);
#ifdef COUNT_OPS
      for (vid x : Si)
        if (us[x] == ue) CNT(stale_pull);
#endif
      Res r = bmssp(l - 1, Bi);
      dupIdx.clear();
      for (size_t i = r.ubeg; i < r.uend; i++) {
        vid u = Ubuf[i];
        if (us[u] == ue) {  // already in U: only possible without the deletion fix
          CNT(dup_u);
          dupIdx.push_back(i);
        }
        us[u] = ue;
        if (D.contains(u)) {
          CNT(stale_erase);
          if (fix) D.erase(u);
        }
      }
      K.clear();
      for (size_t i = r.ubeg; i < r.uend; i++) {
        vid u = Ubuf[i];
        for (uint32_t e = g.head[u]; e < g.head[u + 1]; e++) {
          CNT(scan);
          CNT(scan_main);
          vid v = g.to[e];
          Label c = L.through(u, v, g.w[e]);
          if (L.notWorse(c)) {
            CNT(relaxed);
            L.set(c);
            if (!lt(c, Bi)) {
              if (lt(c, B)) D.insert(c);
            } else if (!lt(c, r.Bp)) {
              K.push_back(c);
            }
          }
        }
      }
      for (vid x : Si) {
        Label lx = L.get(x);
        if (!lt(lx, r.Bp) && lt(lx, Bi)) K.push_back(lx);
      }
      D.batchPrepend(K);
      Bp = r.Bp;
      if (!dupIdx.empty()) dropDuplicates(r.ubeg);  // U <- U u U_i is a set union
    }
    if (lt(B, Bp)) Bp = B;
#ifdef COUNT_OPS
    if (lt(Bp, B)) CNT(partial);
#endif
    for (vid x : Wbuf[l])
      if (us[x] != ue && lt(L.get(x), Bp)) {
        us[x] = ue;
        Ubuf.push_back(x);
      }
    D.clear();
    return Res{Bp, ubeg, Ubuf.size()};
  }

  void run(vid s) {
    L.reset(g.n, s);
    Ubuf.clear();
    Sbuf[top].assign(1, s);
    Res r = bmssp(top, INF_LABEL);
    (void)r;
    CNTN(u_total, r.uend - r.ubeg);
  }
  Len dist(vid v) const { return L.d[v].len; }
};

}  // namespace sp
