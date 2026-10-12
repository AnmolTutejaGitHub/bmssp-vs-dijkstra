// The O(m sqrt(log n) + sqrt(mn log n loglog n)) algorithm of Duan, Mao, Shu
// and Yin (ICALP 2026), written from the pseudocode of the paper:
//   Algorithm 2  FindPivots  (local Dijkstra searches + spanning forest)
//   Algorithm 5  Partition   (forest -> edge-disjoint subtrees of Theta(k) vertices)
//   Algorithm 3  BMSSP
//   Algorithm 4  BaseCase    (Dijkstra on a search tree, |U| <= t^3 + 1)
//   Lemma 3.4    search tree of blocks with Insert / Merge / Pull
// Algorithm, line and lemma numbers are those of the full version,
// arXiv:2602.07868v2 (the proceedings version numbers its lemmas 3, 6, 9, 10
// where the full version has 3.1, 3.4, 3.7, 3.8; algorithms and lines agree).
// Two switches deviate from the paper on purpose (see the paper text):
//   fix      delete a vertex from D when a sub-call reports it complete
//   baseFix  return B' = d_B[D] measured after the last relaxation of the base case
#pragma once
#include <cmath>
#include <map>
#include <set>
#include "bmssp25.hpp"  // LabelHeap
#include "common.hpp"

namespace sp {

// Child structures are handed to the parent as a flat list of pairs plus the
// end index and upper bound of every block, in ascending block order.
struct BlockList {
  std::vector<Label> items;
  std::vector<std::pair<size_t, Label>> blocks;  // (end index, upper bound)
  void clear() {
    items.clear();
    blocks.clear();
  }
};

// ------------------------------------------------------------ Lemma 3.4
class BlockTree26 {
  struct Block {
    std::vector<Label> items;
    Label ub;
  };
  typedef std::map<Label, int, LabelLess> Tree;
  size_t M = 1;
  Label B = INF_LABEL;
  size_t cnt = 0;
  std::vector<Block> pool;
  std::vector<int> freeIds;
  Tree tr;  // upper bound -> block
  std::vector<int> blk;
  std::vector<uint32_t> pos;

  int newBlock(const Label& ub) {
    int b;
    if (!freeIds.empty()) {
      b = freeIds.back();
      freeIds.pop_back();
    } else {
      b = (int)pool.size();
      pool.emplace_back();
    }
    pool[b].items.clear();
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
  // Split blocks larger than M into two equal halves around the median.
  void normalize(Tree::iterator it) {
    int b = it->second;
    size_t n = pool[b].items.size();
    if (n <= M) return;
    CNT(ds_split);
    size_t h = n / 2;
    select_kth(pool[b].items.data(), n, h - 1, LabelLess());
    Label ub = pool[b].items[h - 1];
    int nb = newBlock(ub);
    pool[nb].items.assign(pool[b].items.begin(), pool[b].items.begin() + h);
    pool[b].items.erase(pool[b].items.begin(), pool[b].items.begin() + h);
    reindex(nb);
    reindex(b);
    Tree::iterator lowIt = tr.emplace_hint(it, ub, nb);
    normalize(lowIt);
    normalize(it);
  }
  // 0: reject x (the key is stored with a smaller value); 1: store x;
  // 2: store x, an older pair with the same key was deleted (blocks may have gone)
  int admit(const Label& x) {
    int b = blk[x.v];
    if (b < 0) return 1;
    if (!lt(x, pool[b].items[pos[x.v]])) return 0;
    erase(x.v);
    return 2;
  }

 public:
  void resize(vid n) {
    blk.assign(n, -1);
    pos.assign(n, 0);
  }
  bool empty() const { return cnt == 0; }
  size_t size() const { return cnt; }
  bool contains(vid v) const { return blk[v] >= 0; }

  void init(size_t M_, const Label& B_) {
    assert(cnt == 0);
    M = M_;
    B = B_;
    tr.clear();
    tr.emplace(B, newBlock(B));
  }
  void clear() {
    for (auto& kv : tr) {
      for (const Label& x : pool[kv.second].items) blk[x.v] = -1;
      release(kv.second);
    }
    tr.clear();
    cnt = 0;
  }
  // Move the whole content out (ascending block order) and empty the structure.
  void exportTo(BlockList& out) {
    out.clear();
    for (auto& kv : tr) {
      const std::vector<Label>& it = pool[kv.second].items;
      if (it.empty()) continue;
      out.items.insert(out.items.end(), it.begin(), it.end());
      out.blocks.push_back({out.items.size(), kv.first});
    }
    clear();
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
    blk[v] = -1;
    cnt--;
    if (bl.items.empty() && !same(bl.ub, B)) {
      tr.erase(bl.ub);
      release(b);
    }
    return true;
  }

  void insert(const Label& x) {
    CNT(ds_insert);
    if (!admit(x)) return;
    Tree::iterator it = tr.lower_bound(x);
    assert(it != tr.end());
    place(it->second, x);
    cnt++;
    normalize(it);
  }

  // Merge a structure whose values are all smaller than the values stored here.
  void merge(const BlockList& c) {
    const size_t total = c.items.size();
    if (total == 0) return;
    CNTN(ds_batch, total);
    if (total < M) {  // one block, appended to the first block
      for (const Label& x : c.items) {
        if (!admit(x)) continue;
        place(tr.begin()->second, x);
        cnt++;
      }
      normalize(tr.begin());
      return;
    }
    const size_t third = (M + 2) / 3;
    size_t start = 0;
    for (size_t j = 0; j < c.blocks.size(); j++) {
      size_t end = c.blocks[j].first;
      bool last = (j + 1 == c.blocks.size());
      if (!last && !(end - start >= third && total - end >= third)) continue;
      const Label& ub = c.blocks[j].second;
      // every admitted pair of the chunk goes to the block with upper bound ub
      int b = -1;
      Tree::iterator it = tr.end();
      for (size_t i = start; i < end; i++) {
        const Label& x = c.items[i];
        int a = admit(x);
        if (a == 0) continue;
        if (a == 2) b = -1;  // admit may have deleted the block: look it up again
        if (b < 0) {
          it = tr.lower_bound(ub);
          if (it == tr.end() || !same(it->first, ub)) it = tr.emplace_hint(it, ub, newBlock(ub));
          b = it->second;
        }
        place(b, x);
        cnt++;
      }
      if (b >= 0) normalize(tr.find(ub));
      start = end;
    }
  }

  // Returns the M keys with the smallest values and a separator x
  // (all keys and x = B if at most M pairs are stored).
  void pull(std::vector<vid>& out, Label& x) {
    CNT(ds_pull);
    out.clear();
    if (cnt <= M) {
      for (auto& kv : tr)
        for (const Label& y : pool[kv.second].items) out.push_back(y.v);
      clear();
      tr.emplace(B, newBlock(B));
      x = B;
      CNTN(ds_pulled, out.size());
      return;
    }
    size_t need = M;
    Tree::iterator it = tr.begin();
    for (;;) {
      int b = it->second;
      std::vector<Label>& its = pool[b].items;
      size_t s = its.size();
      if (s <= need) {  // the whole block is among the M smallest
        for (const Label& y : its) {
          out.push_back(y.v);
          blk[y.v] = -1;
        }
        cnt -= s;
        need -= s;
        release(b);
        it = tr.erase(it);
        continue;
      }
      select_kth(its.data(), s, need, LabelLess());
      x = its[need];
      for (size_t i = 0; i < need; i++) {
        out.push_back(its[i].v);
        blk[its[i].v] = -1;
      }
      its.erase(its.begin(), its.begin() + need);
      cnt -= need;
      reindex(b);
      break;
    }
    // normalisation: an under-sized remainder joins the next block
    {
      int b = it->second;
      Tree::iterator nx = std::next(it);
      if (pool[b].items.size() * 3 < M && nx != tr.end()) {
        int nb = nx->second;
        for (const Label& y : pool[b].items) place(nb, y);
        release(b);
        tr.erase(it);
        normalize(nx);
      }
    }
    CNTN(ds_pulled, out.size());
  }
};

// Search tree with M = 1 for the base case (Algorithm 4, line 1).
class BaseSet26 {
  std::set<Label, LabelLess> s;
  std::vector<uint8_t> in;
  std::vector<VLabel> key;
  Label B = INF_LABEL;

 public:
  void resize(vid n) {
    in.assign(n, 0);
    key.resize(n);
  }
  bool empty() const { return s.empty(); }
  size_t size() const { return s.size(); }
  void init(const Label& B_) {
    assert(s.empty());
    B = B_;
  }
  void insert(const Label& x) {
    CNT(ds_insert);
    if (in[x.v]) {
      Label old{key[x.v].len, key[x.v].hops, x.v, key[x.v].pred};
      if (!lt(x, old)) return;
      s.erase(old);
    }
    s.insert(x);
    in[x.v] = 1;
    key[x.v] = VLabel{x.len, x.hops, x.pred};
  }
  // smallest pair and the separator (next smallest value, or B)
  void pull(Label& u, Label& x) {
    CNT(ds_pull);
    CNT(ds_pulled);
    auto it = s.begin();
    u = *it;
    s.erase(it);
    in[u.v] = 0;
    x = s.empty() ? B : *s.begin();
  }
  Label minOrB() const { return s.empty() ? B : *s.begin(); }
  void exportTo(BlockList& out) {
    out.clear();
    for (const Label& x : s) {
      out.items.push_back(x);
      out.blocks.push_back({out.items.size(), x});
      in[x.v] = 0;
    }
    s.clear();
  }
};

// ---------------------------------------------------------------- BMSSP
struct Bmssp26 {
  const Graph& g;
  LabelArray L;
  int t = 2, k = 2, top = 1;
  bool fix = true;      // delete completed vertices from D
  bool baseFix = true;  // base case returns d_B[D] after its last relaxation

  struct Res {
    Label Bp;
    size_t ubeg, uend;
  };
  struct Groups {  // the sets P_j of one call
    std::vector<vid> mem;
    std::vector<uint32_t> off;
    std::vector<vid> pivot;
    std::vector<uint32_t> live;
    std::vector<uint8_t> inJ;
    size_t count() const { return pivot.size(); }
  };

  std::vector<BlockTree26> Ds;
  BaseSet26 base;
  std::vector<uint64_t> Mof, capOf;  // t 2^{(l-1)t} and t^3 2^{lt}
  uint64_t t3 = 8;
  std::vector<std::vector<vid>> Sbuf, Wbuf;
  std::vector<Groups> Gs;
  std::vector<std::vector<int32_t>> grpOf;  // index j of the set P_j containing v, per level
  std::vector<std::vector<uint32_t>> ust;
  std::vector<uint32_t> uEpoch;
  std::vector<vid> Ubuf;
  std::vector<uint32_t> J;
  std::vector<size_t> dupIdx;
  BlockList handoff;
  // FindPivots scratch
  std::vector<uint32_t> inF, inQ, inS, inW, inK, ext, sst;
  uint32_t fEpoch = 0, kEpoch = 0, sEpoch = 0;
  std::vector<vid> kpar, Kv, FV;
  std::vector<uint32_t> lid;
  std::vector<std::pair<vid, vid>> FE;
  std::vector<uint32_t> foff, fadj, gmemLocal, goff, ustack;
  std::vector<uint8_t> fvis;
  struct Frame {
    uint32_t v, ei, ustart;
  };
  std::vector<Frame> frames;
  LabelHeap heap;

  // Parameters of the paper: t = ceil(sqrt(log n loglog n / delta)),
  // k = ceil(t / log t). delta is the degree bound of the input graph.
  static void defaults(vid n, double delta, int& t, int& k) {
    double lg = std::log2((double)std::max<vid>(n, 4));
    double lglg = std::log2(lg);
    t = std::max(2, (int)std::ceil(std::sqrt(lg * lglg / delta) - 1e-9));
    k = kOf(t);
  }
  static int kOf(int t) {
    if (t <= 2) return std::max(1, t);
    return std::max(1, (int)std::ceil(t / std::log2((double)t) - 1e-9));
  }

  Bmssp26(const Graph& g_, int t_ = 0, int k_ = 0, double delta = 3.0) : g(g_) {
    defaults(g.n, delta, t, k);
    if (t_ > 0) {
      t = t_;
      k = kOf(t);
    }
    if (k_ > 0) k = k_;
    double lg = std::log2((double)std::max<vid>(g.n, 2));
    top = std::max(1, (int)std::ceil(lg / t - 1e-9));
    t3 = (uint64_t)t * t * t;
    Ds.resize(top + 1);
    Mof.assign(top + 1, 1);
    capOf.assign(top + 1, 1);
    Sbuf.resize(top + 1);
    Wbuf.resize(top + 1);
    Gs.resize(top + 1);
    grpOf.resize(top + 1);
    ust.resize(top + 1);
    uEpoch.assign(top + 1, 0);
    for (int l = 1; l <= top; l++) {
      Ds[l].resize(g.n);
      Mof[l] = (uint64_t)t * pow2cap((uint64_t)(l - 1) * t);
      capOf[l] = t3 * pow2cap((uint64_t)l * t);
      grpOf[l].assign(g.n, -1);
      ust[l].assign(g.n, 0);
    }
    base.resize(g.n);
    for (auto* a : {&inF, &inQ, &inS, &inW, &inK, &ext, &sst, &lid}) a->assign(g.n, 0);
    kpar.assign(g.n, NONE);
  }

  // Relax(u, v, B) of Algorithm 1 in the paper. On success c holds the new label.
  bool relax(vid u, uint32_t e, const Label& B, Label& c) {
    c = L.through(u, g.to[e], g.w[e]);
    if (!L.notWorse(c) || !lt(c, B)) return false;
    CNT(relaxed);
    L.set(c);
    return true;
  }

  // Algorithm 2 followed by Algorithm 5. Fills Gs[l] (the sets P_j) and Wbuf[l];
  // membership in Q is recorded in inQ.
  void findPivots(int l, const Label& B) {
    const std::vector<vid>& S = Sbuf[l];
    std::vector<vid>& W = Wbuf[l];
    Groups& G = Gs[l];
    std::vector<int32_t>& grp = grpOf[l];
    W.clear();
    FV.clear();
    FE.clear();
    ++fEpoch;
    for (vid x : S) inS[x] = fEpoch;
    for (vid x : S) {
      if (inF[x] == fEpoch) continue;  // already inside some tree
      ++kEpoch;
      Kv.clear();
      Kv.push_back(x);
      inK[x] = kEpoch;
      heap.clear();
      heap.push(L.get(x));
      bool merged = false;
      while (!heap.empty() && Kv.size() < (size_t)k && !merged) {
        Label top_ = heap.pop();
        vid u = top_.v;
        if (ext[u] == kEpoch) continue;  // duplicate entry
        CNT(cmp);                        // stale-entry test, counted like Dijkstra's
        if (!same(top_, L.get(u))) continue;
        ext[u] = kEpoch;
        for (uint32_t e = g.head[u]; e < g.head[u + 1]; e++) {
          CNT(scan);
          CNT(scan_pivot);
          Label c;
          if (!relax(u, e, B, c)) continue;
          vid v = c.v;
          if (inK[v] != kEpoch) {
            inK[v] = kEpoch;
            Kv.push_back(v);
          }
          kpar[v] = u;  // the edge (u,v) replaces any older edge into v
          if (inF[v] == fEpoch) {
            merged = true;  // the search ran into an existing tree
            break;
          }
          heap.push(c);
        }
      }
      if (merged || Kv.size() >= (size_t)k) {
        for (size_t i = 0; i < Kv.size(); i++) {
          vid y = Kv[i];
          if (i > 0) FE.push_back({kpar[y], y});
          if (inF[y] != fEpoch) {
            inF[y] = fEpoch;
            lid[y] = (uint32_t)FV.size();
            FV.push_back(y);
          }
        }
      } else {  // fewer than k vertices reachable: an arborescence W_j with root x
        inQ[x] = fEpoch;
        for (vid y : Kv)
          if (inW[y] != fEpoch) {
            inW[y] = fEpoch;
            W.push_back(y);
          }
      }
    }
    // --- partition the forest into subtrees with at least k vertices
    G.mem.clear();
    G.off.assign(1, 0);
    G.pivot.clear();
    G.live.clear();
    G.inJ.clear();
    if (FV.empty()) return;
    const uint32_t nf = (uint32_t)FV.size();
    foff.assign(nf + 1, 0);
    for (auto& e : FE) {
      foff[lid[e.first] + 1]++;
      foff[lid[e.second] + 1]++;
    }
    for (uint32_t i = 0; i < nf; i++) foff[i + 1] += foff[i];
    fadj.resize(foff[nf]);
    {
      std::vector<uint32_t>& at = ustack;  // borrowed as a cursor array
      at.assign(foff.begin(), foff.end() - 1);
      for (auto& e : FE) {
        uint32_t a = lid[e.first], b = lid[e.second];
        fadj[at[a]++] = b;
        fadj[at[b]++] = a;
      }
    }
    fvis.assign(nf, 0);
    ustack.clear();
    gmemLocal.clear();
    goff.assign(1, 0);
    const uint32_t s = (uint32_t)k;
    for (uint32_t r0 = 0; r0 < nf; r0++) {
      if (fvis[r0]) continue;
      size_t firstGroup = goff.size() - 1;
      fvis[r0] = 1;
      frames.clear();
      frames.push_back({r0, foff[r0], 0});
      ustack.push_back(r0);
      while (!frames.empty()) {
        Frame& f = frames.back();
        if (f.ei < foff[f.v + 1]) {
          uint32_t c = fadj[f.ei++];
          if (!fvis[c]) {
            fvis[c] = 1;
            uint32_t us_ = (uint32_t)ustack.size();
            frames.push_back({c, foff[c], us_});
            ustack.push_back(c);
          }
          continue;
        }
        Frame done = f;
        frames.pop_back();
        if (!frames.empty()) {
          Frame& p = frames.back();
          if (ustack.size() - p.ustart >= s) {  // report a group, keep the root
            gmemLocal.insert(gmemLocal.end(), ustack.begin() + p.ustart, ustack.end());
            goff.push_back((uint32_t)gmemLocal.size());
            ustack.resize(p.ustart);
            ustack.push_back(p.v);
          }
        } else {  // the remainder joins the last group of this tree
          gmemLocal.insert(gmemLocal.end(), ustack.begin() + done.ustart, ustack.end());
          if (goff.size() - 1 == firstGroup)
            goff.push_back((uint32_t)gmemLocal.size());
          else
            goff.back() = (uint32_t)gmemLocal.size();
          ustack.resize(done.ustart);
        }
      }
    }
    // --- P_j = vertices of S \ Q in group j that are in no earlier group
    for (size_t j = 0; j + 1 < goff.size(); j++) {
      size_t before = G.mem.size();
      int32_t id = (int32_t)G.pivot.size();
      for (uint32_t i = goff[j]; i < goff[j + 1]; i++) {
        vid v = FV[gmemLocal[i]];
        if (inS[v] == fEpoch && inQ[v] != fEpoch && grp[v] < 0) {
          grp[v] = id;
          G.mem.push_back(v);
        }
      }
      if (G.mem.size() == before) continue;
      G.off.push_back((uint32_t)G.mem.size());
      G.pivot.push_back(NONE);
      G.live.push_back((uint32_t)(G.mem.size() - before));
      G.inJ.push_back(0);
    }
  }

  // p_j <- argmin of d[] over the vertices still in P_j
  Label selectPivot(int l, uint32_t j) {
    Groups& G = Gs[l];
    std::vector<int32_t>& grp = grpOf[l];
    Label best = INF_LABEL;
    for (uint32_t i = G.off[j]; i < G.off[j + 1]; i++) {
      vid v = G.mem[i];
      if (grp[v] != (int32_t)j) continue;
      Label lv = L.get(v);
      if (lt(lv, best)) best = lv;
    }
    G.pivot[j] = best.v;
    return best;
  }

  // Algorithm 4. The remaining frontier stays in `base`.
  Res baseCase(const Label& B) {
    CNT(base_calls);
    base.init(B);
    for (vid x : Sbuf[0]) base.insert(L.get(x));
    const size_t ubeg = Ubuf.size();
    Label Bp = B;
    while (!base.empty() && Ubuf.size() - ubeg <= t3) {
      Label ul;
      base.pull(ul, Bp);
      vid u = ul.v;
      Ubuf.push_back(u);
      for (uint32_t e = g.head[u]; e < g.head[u + 1]; e++) {
        CNT(scan);
        CNT(scan_base);
        Label c;
        if (relax(u, e, B, c)) base.insert(c);
      }
    }
    if (baseFix) Bp = base.minOrB();
#ifdef COUNT_OPS
    if (lt(Bp, B)) CNT(partial);
#endif
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

  // Algorithm 3. Frontier in Sbuf[l]; U returned as a range of Ubuf; the
  // returned structure D is Ds[l] (or `base` for l = 0).
  Res bmssp(int l, const Label& B) {
    if (l == 0) return baseCase(B);
    const std::vector<vid>& S = Sbuf[l];
    CNT(calls);
    CNTN(frontier, S.size());
    BlockTree26& D = Ds[l];
    D.init(Mof[l], B);
    Groups& G = Gs[l];
    std::vector<int32_t>& grp = grpOf[l];
    findPivots(l, B);
    Label Bp = B;
    for (uint32_t j = 0; j < G.count(); j++) {
      Label lp = selectPivot(l, j);
      CNT(pivots);
      D.insert(lp);
      if (lt(lp, Bp)) Bp = lp;
    }
    const size_t ubeg = Ubuf.size();
    const uint64_t cap = capOf[l];
    const uint32_t ue = ++uEpoch[l];
    std::vector<uint32_t>& us = ust[l];
    std::vector<vid>& Si = Sbuf[l - 1];
    while (Ubuf.size() - ubeg <= cap && !D.empty()) {
      Label Bi;
      D.pull(Si, Bi);
      const size_t np = Si.size();
      ++sEpoch;
      for (size_t i = 0; i < np; i++) sst[Si[i]] = sEpoch;
      for (size_t i = 0; i < np; i++) {
        vid x = Si[i];
#ifdef COUNT_OPS
        if (us[x] == ue) CNT(stale_pull);
#endif
        int32_t j = grp[x];
        if (j < 0 || G.pivot[j] != x) continue;
        for (uint32_t q = G.off[j]; q < G.off[j + 1]; q++) {  // bring in {v in P_j : d[v] < B_i}
          vid v = G.mem[q];
          if (grp[v] != j || sst[v] == sEpoch) continue;
          if (lt(L.get(v), Bi)) {
            sst[v] = sEpoch;
            Si.push_back(v);
          }
        }
      }
      Res r = bmssp(l - 1, Bi);
      if (l == 1)
        base.exportTo(handoff);
      else
        Ds[l - 1].exportTo(handoff);
#ifdef COUNT_OPS
      for (size_t i = r.ubeg; i < r.uend; i++)
        if (D.contains(Ubuf[i])) CNT(stale_parent);
      if (!handoff.items.empty()) {
        bool below = false;
        for (const Label& y : handoff.items) below = below || rawLt(y, r.Bp);
        if (below) CNT(bound_viol);
        if (same(r.Bp, Bi)) CNT(full_nonempty);
      }
#endif
      D.merge(handoff);
      J.clear();
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
        int32_t j = grp[u];
        if (j >= 0) {
          grp[u] = -1;
          G.live[j]--;
          if (G.pivot[j] == u) {
            G.pivot[j] = NONE;
            if (G.live[j] > 0) {
              J.push_back((uint32_t)j);
              G.inJ[j] = 1;
            }
          }
        }
      }
      for (size_t i = r.ubeg; i < r.uend; i++) {
        vid u = Ubuf[i];
        for (uint32_t e = g.head[u]; e < g.head[u + 1]; e++) {
          CNT(scan);
          CNT(scan_main);
          Label c;
          if (!relax(u, e, B, c) || lt(c, Bi)) continue;
          D.insert(c);
          int32_t j = grp[c.v];
          if (j >= 0 && !G.inJ[j] && lt(c, L.get(G.pivot[j]))) G.pivot[j] = c.v;
        }
      }
      for (uint32_t j : J) {  // "for each non-empty P_j from j in J" (Section 3.3)
        G.inJ[j] = 0;
        if (G.live[j] == 0) continue;  // emptied by later vertices of U_i
        Label lp = selectPivot(l, j);
        CNT(pivots);
        D.insert(lp);
      }
      Bp = r.Bp;
      if (!dupIdx.empty()) dropDuplicates(r.ubeg);  // U <- U u U_i is a set union
    }
#ifdef COUNT_OPS
    if (lt(Bp, B)) CNT(partial);
#endif
    for (vid x : S) {
      Label lx = L.get(x);
      if (!lt(lx, Bp) && lt(lx, B)) D.insert(lx);
    }
    const size_t wbeg = Ubuf.size();
    for (vid x : Wbuf[l])
      if (us[x] != ue && lt(L.get(x), Bp)) {  // W' = {x in W \ U : d[x] < B'}
        us[x] = ue;
        Ubuf.push_back(x);
        if (D.contains(x)) {
          CNT(stale_erase);
          if (fix) D.erase(x);
        }
      }
    for (size_t i = wbeg; i < Ubuf.size(); i++) {
      vid u = Ubuf[i];
      for (uint32_t e = g.head[u]; e < g.head[u + 1]; e++) {
        CNT(scan);
        CNT(scan_main);
        Label c;
        if (relax(u, e, B, c) && !lt(c, Bp)) D.insert(c);
      }
    }
    for (vid v : G.mem) grp[v] = -1;
    return Res{Bp, ubeg, Ubuf.size()};
  }

  void run(vid s) {
    L.reset(g.n, s);
    Ubuf.clear();
    Sbuf[top].assign(1, s);
    Res r = bmssp(top, INF_LABEL);
    (void)r;
    CNTN(u_total, r.uend - r.ubeg);
    Ds[top].clear();
  }
  Len dist(vid v) const { return L.d[v].len; }
};

}  // namespace sp
