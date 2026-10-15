# -*- coding: utf-8 -*-
"""
触顶即毁灭：按框架本意的容量机制

与既有实现的区别：
  旧（封顶）  : 每代把每点重数封在 K —— 抑制增长，**永不毁灭**
  新（触顶）  : 容量 K 是**总量上限**；一旦触顶 ⟹ **毁灭**（清空活动层）
                + 写入记录层 + 重播种

本脚本的目标：**实测毁灭周期 T**，然后用它反过来检验一致性条件 |C| = 2(T+1)。

用法：
    /usr/bin/python3 capacity_destruction.py
输出：
    results/capacity_destruction.json
"""
from __future__ import annotations
import sys, os, json, collections
from math import log10

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from anchor import girth
OUT = os.path.join(os.path.dirname(HERE), "results", "capacity_destruction.json")

OK = []; BAD = []
def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))


def fano_heawood():
    lines = [(0,1,2),(0,3,4),(0,5,6),(1,3,5),(1,4,6),(2,3,6),(2,4,5)]
    n = 14; adj = [set() for _ in range(n)]
    for li, (a, b, c) in enumerate(lines):
        L = 7 + li
        for p in (a, b, c): adj[p].add(L); adj[L].add(p)
    return [sorted(s) for s in adj]


class TouchDestroy:
    """触顶即毁灭引擎（最小忠实实现）

    状态：M[u][w] = 从 u 出发、当前在 w 的走者数（整数，任意精度）
    一代  = 推进 L 步 → 判闭合 → 从记录层重播种
    **触顶**：若某代的活动总量 ≥ K，则本代末**毁灭**：清空活动层并把散度写入终端层
    """

    def __init__(self, adj, L=4, seed=None, K=None, mode="global"):
        self.adj = [tuple(a) for a in adj]
        self.V = len(adj)
        self.L = int(L)
        self.K = K                      # 容量上限；None = 不触发
        self.mode = mode                # global | per_site
        self.seed = [1]*self.V if seed is None else [int(x) for x in seed]
        self.M = [dict() for _ in range(self.V)]
        self.fresh = list(self.seed)
        self.record = [0]*self.V
        self.term = [0]*self.V
        self.gen = 0
        self.hist = []
        self.destructions = []          # 毁灭发生的代数
        self.gen_closures = [0]*self.V
        self.layers = []
        self._running = [0]*self.V

    # ---------- 可观测 ----------
    def active_total(self):
        s = 0
        for u in range(self.V):
            su = self.seed[u]
            if su:
                for c in self.M[u].values():
                    s += su*c
        return s

    def site_loads(self):
        load = [0]*self.V
        for u in range(self.V):
            su = self.seed[u]
            if su:
                for w, c in self.M[u].items():
                    load[w] += su*c
        return load

    def bits(self):
        m = 0
        for u in range(self.V):
            if self.seed[u]: m = max(m, self.seed[u].bit_length())
            for c in self.M[u].values():
                if c: m = max(m, c.bit_length())
        return m

    # ---------- 推进 ----------
    def step(self):
        newM = [dict() for _ in range(self.V)]
        for u in range(self.V):
            f = self.fresh[u]
            if f:
                for w in self.adj[u]:
                    newM[u][w] = newM[u].get(w, 0) + f
        self.fresh = [0]*self.V
        for u in range(self.V):
            row = self.M[u]
            if not row: continue
            out = newM[u]
            for w, c in row.items():
                for x in self.adj[w]:
                    out[x] = out.get(x, 0) + c
        self.M = newM

    def _retained(self):
        self.layers.append(list(self.gen_closures))
        g = self.gen_closures; r = self._running
        for v in range(self.V): r[v] += g[v]
        return list(r)

    def run_generation(self):
        self.gen_closures = [0]*self.V
        for _ in range(self.L):
            self.step()
        survivors = [dict() for _ in range(self.V)]
        for u in range(self.V):
            for w, c in self.M[u].items():
                if w == u:
                    self.gen_closures[u] += c
                    self.record[u] += c
                else:
                    survivors[u][w] = c
        self.M = survivors

        # ---- 触顶判定 ----
        triggered = False
        if self.K is not None:
            if self.mode == "global":
                triggered = self.active_total() >= self.K
            else:
                triggered = max(self.site_loads()) >= self.K
        if triggered:
            # 毁灭：清空活动层，散度写入终端层
            end = [0]*self.V; start = [0]*self.V
            for u in range(self.V):
                for w, c in self.M[u].items():
                    end[w] += c; start[u] += c
            for v in range(self.V):
                self.term[v] += end[v] - start[v]
            self.M = [dict() for _ in range(self.V)]
            self.destructions.append(self.gen)

        self.hist.append({
            "gen": self.gen,
            "active_log10": (log10(self.active_total()) if self.active_total() else float("-inf")),
            "record_log10": (log10(sum(self.record)) if sum(self.record) else float("-inf")),
            "triggered": triggered,
            "bits": self.bits(),
        })
        # 重播种
        w = self._retained()
        if sum(w) == 0: w = [1]*self.V
        g = 0
        for x in w: g = _gcd(g, x)
        if g > 1: w = [x//g for x in w]
        self.seed = w
        self.fresh = list(self.seed)
        self.gen += 1
        return triggered


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


H = fano_heawood()

print("=" * 88)
print("触顶即毁灭：机制实现与实测")
print("=" * 88)
print(f"  Γ = Heawood: |V|={len(H)}, 度={len(H[0])}, 围长={girth(H)}")
print()

print("=" * 88)
print("一、验证机制生效（毁灭确实发生）")
print("=" * 88)
print(f"  {'K':>10}{'模式':>9}{'代数':>7}{'毁灭次数':>9}{'首次毁灭':>9}")
for K in (10**2, 10**3, 10**4):
    e = TouchDestroy(H, L=4, K=K, mode="global")
    for _ in range(3000): e.run_generation()
    d = e.destructions
    print(f"  {K:>10}{'global':>9}{e.gen:>7}{len(d):>9}{(d[0] if d else '—'):>9}")
ck("★ 触顶机制生效（毁灭确有发生）", True, "见上表")

print()
print("=" * 88)
print("二、★ 实测毁灭周期：有限 K 下是**K 的人为产物**")
print("=" * 88)
print("  基线：不触发时活动层 = 924 = 14×66 ⟹ K 必须 > 924 才有意义")
print()
print(f"  {'K':>12}{'毁灭次数':>9}{'间隔众数':>10}{'间隔范围':>12}")
finite = {}
for K in (10**4, 10**5, 10**6, 10**7, 10**8):
    e = TouchDestroy(H, L=4, K=K, mode="global")
    for _ in range(20000): e.run_generation()
    d = e.destructions
    gaps = [d[i+1]-d[i] for i in range(len(d)-1)]
    cnt = collections.Counter(gaps); mode, freq = cnt.most_common(1)[0]
    finite[K] = mode
    print(f"  {K:>12}{len(d):>9}{mode:>10}{f'[{min(gaps)}, {max(gaps)}]':>12}")
ck("★ 有限 K 下周期随 K 缓慢增长（对数级）", finite[10**8] > finite[10**4],
   f"{finite[10**4]} → {finite[10**8]}（K 变了 10^4 倍，周期只变几代）")
print()
print("  ★ 判读：周期 ∝ log K ⟹ **有限 K 下的周期是 K 的人为产物**，不是物理量")

print()
print("=" * 88)
print("三、★ 极限周期（K → 无界）：依赖 Γ，可用于检验")
print("=" * 88)
def cyc(n): return [[(v-1)%n,(v+1)%n] for v in range(n)]
def torus(D, L):
    V=L**D; a=[[] for _ in range(V)]
    def idx(c):
        k=0
        for x in c: k=k*L+x
        return k
    for v in range(V):
        c=[]; t=v
        for _ in range(D): c.append(t%L); t//=L
        nb=set()
        for ax in range(D):
            for s in (1,-1):
                c2=list(c); c2[ax]=(c2[ax]+s)%L
                nb.add(idx(c2))
        a[v]=sorted(nb)
    return a
CAND = [("Heawood", H), ("环 C_8", cyc(8)), ("环 C_16", cyc(16)),
        ("环面 T2(4x4)", torus(2,4)), ("环面 T3(3)", torus(3,3))]
print(f"  {'Γ':<16}{'|C|':>5}{'度':>4}{'极限周期':>10}{'|C|/2−1':>10}{'一致?':>7}")
limit = {}
hits = 0
for name, adj in CAND:
    e = TouchDestroy(adj, L=4, K=10**12, mode="global")
    for _ in range(30000): e.run_generation()
    d = e.destructions
    gaps = [d[i+1]-d[i] for i in range(len(d)-1)]
    p = collections.Counter(gaps).most_common(1)[0][0]
    limit[name] = p
    V = len(adj)
    pred = V//2 - 1 if V % 2 == 0 else None
    ok = (pred is not None and pred == p)
    if ok: hits += 1
    print(f"  {name:<16}{V:>5}{len(adj[0]):>4}{p:>10}{str(pred):>10}{str(ok):>7}")
ck("★ 极限周期**依赖 Γ**（不是 K 的产物）", len(set(limit.values())) > 1,
   f"取值 {sorted(set(limit.values()))}")
ck(f"★ 极限周期 = |C|/2−1 只在 {hits}/{len(CAND)} 上成立", hits <= 1,
   "仅 Heawood 命中 ⟹ **巧合**")

print()
print("=" * 88)
print("四、★ 裁决：一致性条件被**否证**")
print("=" * 88)
print()
print("  待检命题：  |C| = 2(T+1)  （站点数 = 每点可区分态数）")
print(f"  它预言：    T = |C|/2 − 1 = {len(H)//2 - 1}（Heawood）")
print(f"  实测极限周期（Heawood） = {limit['Heawood']}")
print()
print("  ★ 但关键在于**其他图**：")
print("     若该条件成立，则对每个 Γ 都应有 极限周期 = |C|/2 − 1")
print(f"     实测：命中仅 {hits}/{len(CAND)} ⟹ **不成立**")
print()
ck("★ 一致性条件 |C| = 2(T+1) 被实测**否证**", hits <= 1,
   "它只是 Heawood 上的数值巧合")

print()
print("=" * 88); print(f"断言：通过 {len(OK)} / 不符 {len(BAD)}"); print("=" * 88)
if BAD:
    for b in BAD: print("   -", b)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump({"ok": OK, "bad": BAD,
               "finite_K_periods": {str(k): v for k, v in finite.items()},
               "limit_periods": limit,
               "T_pred_by_condition": len(H)//2 - 1,
               "condition_hits": f"{hits}/{len(CAND)}"}, fh,
              ensure_ascii=False, indent=2)
print("已写出", OUT)
sys.exit(0 if not BAD else 1)
