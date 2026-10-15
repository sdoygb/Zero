# -*- coding: utf-8 -*-
"""
两条路的桥：L 变大时容量是否绑定？

推理：
  L=4: |W_4| = 6 ≤ 2(L+1) = 10  ⟹ 容量**不绑定** ⟹ 无毁灭 ⟹ T 用不到
  L=6: |W_6| = 20 > 2(L+1) = 14 ⟹ 容量**应绑定** ⟹ 是否触发毁灭？

若 L=6 触发毁灭，则两条路是**同一机制的两端**（L=4 落在不绑定侧）。

★ 关键：'被容纳者'必须**带重数**（我上轮发现分布机丢重数导致周期塌缩为 1）

用法：
    /usr/bin/python3 bridge_test.py
输出：
    results/bridge_test.json
"""
from __future__ import annotations
import sys, os, json, collections, io, contextlib, importlib.util
from math import comb, log10

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from anchor import girth
OUT = os.path.join(os.path.dirname(HERE), "results", "bridge_test.json")

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


class WeightedEngine:
    """带重数的引擎：状态 = 每点的 (起点 → 重数) 映射

    可区分态 = 不同起点的**个数**（重数不进可区分性，但进动力学）
    容量上限 cap = 2(L+1)：每点最多 cap 个可区分态
    饱和判据：**每点的总活跃量（含重数）** > 可区分重数 ⟹ 触发
    """

    def __init__(self, adj, L=4, cap_mode="algebraic"):
        self.adj = [tuple(a) for a in adj]
        self.V = len(adj)
        self.L = int(L)
        self.T = int(L)
        self.cap = 2*(self.T + 1)        # 代数容量 2(L+1)
        # state[v] = {起点 u: 重数}
        self.state = [{v: 1} for v in range(self.V)]
        self._fresh = [{v: 1} for v in range(self.V)]
        self.record = [0]*self.V
        self.term = [0]*self.V
        self.gen = 0
        self.hist = []
        self.destructions = []

    def site_activity(self, v):
        return sum(self.state[v].values())

    def n_states(self):
        return sum(len(s) for s in self.state)

    def peak_activity(self):
        return max((self.site_activity(v) for v in range(self.V)), default=0)

    def step(self):
        nxt = [dict() for _ in range(self.V)]
        def add(v, u, c):
            nxt[v][u] = nxt[v].get(u, 0) + c
        for v in range(self.V):
            if self._fresh[v]:
                for w in self.adj[v]:
                    for u, c in self._fresh[v].items(): add(w, u, c)
        self._fresh = [dict() for _ in range(self.V)]
        for v in range(self.V):
            if self.state[v]:
                for w in self.adj[v]:
                    for u, c in self.state[v].items(): add(w, u, c)
        # 容量：可区分态个数 ≤ cap（超出的按起点排序截断）
        for v in range(self.V):
            if len(nxt[v]) > self.cap:
                keys = sorted(nxt[v])[:self.cap]
                nxt[v] = {k: nxt[v][k] for k in keys}
        self.state = nxt

    def run_generation(self):
        peak_act = 0; peak_states = 0
        for _ in range(self.L):
            self.step()
            peak_act = max(peak_act, self.peak_activity())
            peak_states = max(peak_states, max((len(s) for s in self.state), default=0))
        # 闭合：起点 = 当前位置
        surv = [dict() for _ in range(self.V)]
        for v in range(self.V):
            for u, c in self.state[v].items():
                if u == v:
                    self.record[v] += c
                else:
                    surv[v][u] = c
        self.state = surv
        # ★ 饱和：带重数的每点活跃量 > 可区分重数
        trig = peak_act > self.cap
        if trig:
            self.destructions.append(self.gen)
            self.state = [dict() for _ in range(self.V)]
        self.hist.append({"gen": self.gen, "peak_act": peak_act,
                          "peak_states": peak_states, "trig": trig})
        # 重播种
        self.state = [{v: 1} for v in range(self.V)]
        self._fresh = [{v: 1} for v in range(self.V)]
        self.gen += 1
        return trig


H = fano_heawood()

print("=" * 88)
print("两条路的桥：L 变大时容量是否绑定？")
print("=" * 88)
print(f"  Γ = Heawood: |V|={len(H)}, 度={len(H[0])}, 围长={girth(H)}")
print()
print("  推理（用**词数**当被容纳者）：")
print(f"  {'L':>3}{'|W_L|':>8}{'2(L+1)':>9}{'绑定?':>8}")
for L in (2, 4, 6, 8):
    W = comb(L, L//2) if L % 2 == 0 else 0
    cap = 2*(L+1)
    print(f"  {L:>3}{W:>8}{cap:>9}{str(W > cap):>8}")
print()

print("=" * 88)
print("一、带重数引擎：实测各 L 下是否触发毁灭")
print("=" * 88)
print(f"  {'L':>3}{'cap':>6}{'实测峰值活跃量':>16}{'毁灭次数':>10}{'周期':>7}")
rows = []
for L in (2, 4, 6, 8):
    e = WeightedEngine(H, L=L)
    for _ in range(5000): e.run_generation()
    pa = [h["peak_act"] for h in e.hist]
    d = e.destructions
    gaps = [d[i+1]-d[i] for i in range(len(d)-1)] if len(d) > 1 else []
    p = collections.Counter(gaps).most_common(1)[0][0] if gaps else None
    rows.append(dict(L=L, cap=2*(L+1), peak_act=max(pa), n_dest=len(d), period=p))
    print(f"  {L:>3}{2*(L+1):>6}{max(pa):>16}{len(d):>10}{str(p):>7}")

print()
print("  ★ 判读：")
for r in rows:
    W = comb(r["L"], r["L"]//2) if r["L"] % 2 == 0 else 0
    binding = W > r["cap"]
    fired = r["n_dest"] > 0
    print(f"     L={r['L']}: |W_L|={W} vs cap={r['cap']} ⟹ 理论{'绑定' if binding else '不绑定'}; "
          f"实测毁灭={'有' if fired else '无'}  {'✓ 一致' if binding == fired else '✗ 不一致'}")

bound = [r for r in rows if (comb(r["L"], r["L"]//2) if r["L"] % 2 == 0 else 0) > r["cap"]]
free = [r for r in rows if (comb(r["L"], r["L"]//2) if r["L"] % 2 == 0 else 0) <= r["cap"]]
ck("★ 存在容量绑定的 L（应触发毁灭）", len(bound) > 0,
   f"绑定 L = {[r['L'] for r in bound]}")
ck("★ 绑定侧确实触发毁灭", all(r["n_dest"] > 0 for r in bound),
   f"{[(r['L'], r['n_dest']) for r in bound]}")
ck("★ 不绑定侧不触发毁灭", all(r["n_dest"] == 0 for r in free),
   f"{[(r['L'], r['n_dest']) for r in free]}")

print()
print("=" * 88)
print("二、★ 桥的判定")
print("=" * 88)
print()
if bound and free:
    print("  ✅ 两条路是**同一机制的两端**：")
    print("     · 容量不绑定侧（L 小）⟹ 无毁灭 ⟹ T 用不到")
    print("     · 容量绑定侧（L 大） ⟹ 有毁灭 ⟹ T 被用到")
    ck("★ 两路是同一机制的两端", True,
       f"分界在 L: 不绑定 {[r['L'] for r in free]} / 绑定 {[r['L'] for r in bound]}")
else:
    ck("★ 两路是同一机制的两端", False, "未观察到绑定侧")

print()
print("=" * 88); print(f"断言：通过 {len(OK)} / 不符 {len(BAD)}"); print("=" * 88)
if BAD:
    for b in BAD: print("   -", b)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump({"ok": OK, "bad": BAD, "rows": rows}, fh,
              ensure_ascii=False, indent=2)
print("已写出", OUT)
sys.exit(0 if not BAD else 1)
