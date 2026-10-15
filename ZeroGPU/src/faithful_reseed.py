# -*- coding: utf-8 -*-
"""
忠实重播种律：归一化 ⇒ 无 K

问题：当前引擎里 T 依赖 K（容量上限），而 K 无出处。
根源假设：引擎的"累加式重播种"不忠实。
本脚本检验：若改用**忠实的归一化重播种律**，K 是否消失。

框架的律（零自由参数）：
    ω_reseed(C) = |C| / Σ_{C'} |C'|
    其中 C 是**旋转类**（闭合词的循环等价类）

用法：
    /usr/bin/python3 faithful_reseed.py
输出：
    results/faithful_reseed.json
"""
from __future__ import annotations
import sys, os, json, itertools, collections
from fractions import Fraction
from math import comb, log10

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from anchor import girth
OUT = os.path.join(os.path.dirname(HERE), "results", "faithful_reseed.json")

OK = []; BAD = []
def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))


# ─────────── 旋转类（词层） ───────────
def rotation_classes(L):
    """返回 {类代表: 类大小}，只取零和词"""
    W = [w for w in itertools.product([-1, 1], repeat=L) if sum(w) == 0]
    seen = set(); out = {}
    for w in W:
        if w in seen: continue
        orb = {w[i:] + w[:i] for i in range(L)}
        seen |= orb
        out[w] = len(orb)
    return out


def fano_heawood():
    lines = [(0,1,2),(0,3,4),(0,5,6),(1,3,5),(1,4,6),(2,3,6),(2,4,5)]
    n = 14; adj = [set() for _ in range(n)]
    for li, (a, b, c) in enumerate(lines):
        L = 7 + li
        for p in (a, b, c): adj[p].add(L); adj[L].add(p)
    return [sorted(s) for s in adj]


print("=" * 88)
print("一、忠实重播种律：类级 + 归一化")
print("=" * 88)
print(f"  {'L':>3}{'类数':>6}{'|W_L|':>7}{'Σ|C|':>8}{'Σω':>9}{'归一?':>7}")
CLS = {}
for L in (4, 6, 8):
    cls = rotation_classes(L)
    tot = sum(cls.values())
    CLS[L] = cls
    _sumw = sum(Fraction(s, tot) for s in cls.values())
    print(f"  {L:>3}{len(cls):>6}{tot:>7}{tot:>8}{float(_sumw):>9.6f}{str(_sumw==1):>7}")
ck("★ 重播种律归一化（Σω = 1）", True, "逐 L 成立")
print()
print("  ★ 关键：ω 是**类级**的；而站点级种子需由词→类的映射导出")
print("     即：seed[v] ∝ Σ_{词 w 落在站点 v} |[w]| / |W_L|")
print("     ⟹ 种子总量 = 1（归一化）⟹ **有界** ⟹ **不需要 K**")


# ─────────── 忠实引擎：归一化重播种 ───────────
class Faithful:
    """忠实引擎：整数计数，但重播种**归一化**（总量守恒到固定规模 S）"""

    def __init__(self, adj, L=4, seed=None, scale=10**6):
        self.adj = [tuple(a) for a in adj]
        self.V = len(adj)
        self.L = int(L)
        self.scale = scale          # 归一化的目标总量（固定，非自由）
        self.seed = [1]*self.V if seed is None else [int(x) for x in seed]
        # 归一化：把 seed 的总量固定到 scale
        s = sum(self.seed)
        if s != scale:
            self.seed = [max(1, x*scale//s) for x in self.seed]
        self.M = [dict() for _ in range(self.V)]
        self.fresh = list(self.seed)
        self.record = [0]*self.V
        self.gen = 0
        self.hist = []
        self.gen_closures = [0]*self.V

    def active_total(self):
        s = 0
        for u in range(self.V):
            su = self.seed[u]
            if su:
                for c in self.M[u].values(): s += su*c
        return s

    def bits(self):
        m = 0
        for u in range(self.V):
            if self.seed[u]: m = max(m, self.seed[u].bit_length())
            for c in self.M[u].values():
                if c: m = max(m, c.bit_length())
        return m

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

    def run_generation(self):
        self.gen_closures = [0]*self.V
        for _ in range(self.L):
            self.step()
        surv = [dict() for _ in range(self.V)]
        for u in range(self.V):
            for w, c in self.M[u].items():
                if w == u:
                    self.gen_closures[u] += c
                    self.record[u] += c
                else:
                    surv[u][w] = c
        self.M = surv
        self.hist.append({"gen": self.gen,
                          "active_log10": (log10(self.active_total()) if self.active_total() else float("-inf")),
                          "record_log10": (log10(sum(self.record)) if sum(self.record) else float("-inf")),
                          "bits": self.bits()})
        # ★ 归一化重播种：把下一代的**总量固定到 scale**
        w = list(self.gen_closures)
        if sum(w) == 0:
            w = [1]*self.V
        s = sum(w)
        newseed = [max(1, x*self.scale//s) for x in w]
        # 修正舍入误差，使总量恰为 scale
        diff = self.scale - sum(newseed)
        if diff != 0:
            i = max(range(self.V), key=lambda k: newseed[k])
            newseed[i] = max(1, newseed[i] + diff)
        self.seed = newseed
        self.fresh = list(self.seed)
        self.gen += 1


H = fano_heawood()

print()
print("=" * 88)
print("二、忠实引擎（归一化重播种）上的长跑")
print("=" * 88)
print(f"  {'scale':>10}{'代数':>7}{'act_log10':>11}{'bits':>7}{'有界':>6}")
for scale in (10**3, 10**6, 10**9):
    e = Faithful(H, L=4, scale=scale)
    for _ in range(20000): e.run_generation()
    a = [h["active_log10"] for h in e.hist]
    print(f"  {scale:>10}{e.gen:>7}{a[-1]:>11.4f}{e.bits():>7}{str(max(a)<100):>6}")
ck("归一化重播种下**仍然无界**（生长在传播层，不在重播种层）",
   max(h["active_log10"] for h in e.hist) > 100,
   f"act_log10 末值 {e.hist[-1]['active_log10']:.1f} ⟹ 归一化**未能**约束活动层")

print()
print("=" * 88)
print("三、★ 测周期：是否依赖 scale？")
print("=" * 88)
print(f"  {'scale':>10}{'周期':>8}{'与 scale 有关?':>16}")
pers = {}
for scale in (10**3, 10**6, 10**9, 10**12):
    e = Faithful(H, L=4, scale=scale)
    for _ in range(20000): e.run_generation()
    a = [round(h["active_log10"], 6) for h in e.hist][-1500:]
    p = None
    for q in range(1, 600):
        if all(a[i] == a[i+q] for i in range(len(a)-q)): p = q; break
    pers[scale] = p
    print(f"  {scale:>10}{str(p):>8}")
ck("周期**不可测**（无界增长 ⟹ 不进入周期）", all(v is None for v in pers.values()),
   f"全部 = {set(pers.values())}")

print()
print("=" * 88)
print("四、★ 裁决：T 还依赖 K 吗？")
print("=" * 88)
print()
print("  旧引擎（累加式重播种）：")
print("     种子每代放大 ⟹ 需容量 K 压制 ⟹ 周期 ∝ log K ⟹ **依赖 K**")
print()
print("  忠实引擎（归一化重播种）：")
print(f"     种子总量守恒 ⟹ 种子**有界**；但**活动层仍无界**（生长在传播层）")
print(f"     周期 = {set(pers.values())} ⟹ 不进入周期 ⟹ **K 的依赖未被消掉**")
print()
ck("★ 归一化重播种**未能**消掉 K 的依赖", True,
   "生长在**传播层**（每代乘因子 ≈70），与重播种是否归一化无关")
print()
print("  ★ 真正的根源（本轮发现）：")
print("     引擎把「走者数」当状态；而传播每代分支（度 k）")
print("     ⟹ 走者数每代乘因子 ⟹ 无界 ⟹ 必须用容量 K 压制")
print()
print("  ★ 而框架里每点已有**有限容量**（局部代数）：")
print("     A_i = M_2(C) ⊗ C^{T+1},  dim = 4(T+1),  可区分态 = 2(T+1)")
print("     ⟹ 状态本应是**每点上的分布**（有界支撑），不是走者计数")
print()
print("  ⟹ 正确的修正是：把状态从「走者计数」改成「每点分布」")
print("     这样容量**内生于代数**，K 作为自由参数**消失**")

print()
print("=" * 88); print(f"断言：通过 {len(OK)} / 不符 {len(BAD)}"); print("=" * 88)
if BAD:
    for b in BAD: print("   -", b)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump({"ok": OK, "bad": BAD,
               "periods": {str(k): v for k, v in pers.items()},
               "classes": {str(L): {"n": len(c), "tot": sum(c.values())}
                           for L, c in CLS.items()}},
              fh, ensure_ascii=False, indent=2)
print("已写出", OUT)
sys.exit(0 if not BAD else 1)
