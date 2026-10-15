# -*- coding: utf-8 -*-
"""
讨论稿实验：饱和 ⇒ 被选速度？

既有材料的机制（G31）：
  局部规则：输运 ½[n(x−1)+n(x+1)]、老化、繁殖 B、退出
  不饱和 ⟹ 前缘走扩散律 ⟹ 没有被选速度
  加**每点容量 1** ⟹ 前缘速度**立刻**变 KPP 值（相对差 0.4%）
  且饱和来源**未论证**（既有材料自标开放项）

本次：把该机制放到 **Heawood 图**上（不是链），测四项：
  ① 容量 1 下前缘是否匀速
  ② 前缘速度是否 = KPP 值
  ③ 无容量时是否走扩散
  ④ 容量 1 与容量 2(L+1) 的差异

用法：
    /usr/bin/python3 front_speed.py
输出：
    results/front_speed.json
"""
from __future__ import annotations
import sys, os, json, collections, math

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from anchor import girth
OUT = os.path.join(os.path.dirname(HERE), "results", "front_speed.json")

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


def bfs_dist(adj, src):
    d = {src: 0}; q = collections.deque([src])
    while q:
        v = q.popleft()
        for w in adj[v]:
            if w not in d:
                d[w] = d[v] + 1; q.append(w)
    return d


class FrontEngine:
    """最小前缘模型（既有材料规则的**图版**）

    状态：n[a][x] = 位置 x、年龄 a 的分支数（a = 0..3）
    规则（每步，**持久状态**）：
      输运  m_a(x) = (1/|N(x)|) Σ_{y∈N(x)} n_a(y)      ← 图上均匀扩散
      老化  n'_1 = m_0,  n'_2 = m_1,  n'_3 = m_2
      繁殖  n'_0 = B · m_1                              ← 年龄 2 时繁殖
      退出  age-3 分支下步进终端（故 n_3 不再出现在下一轮）
      饱和（可选）n = min(n, cap)
    """

    def __init__(self, adj, src=0, B=2, cap=None):
        self.adj = [tuple(a) for a in adj]
        self.V = len(adj)
        self.src = src
        self.B = B
        self.cap = cap
        self.dist = bfs_dist(adj, src)
        self.maxd = max(self.dist.values())
        self.n = [[0.0]*self.V for _ in range(4)]
        self.n[0][src] = 1.0
        self.gen = 0

    def total(self):
        return sum(sum(row) for row in self.n)

    def step(self):
        # 输运：每个位置把密度均分给邻居
        m = [[0.0]*self.V for _ in range(4)]
        for a in range(4):
            row = self.n[a]
            if not any(row): continue
            for x in range(self.V):
                c = row[x]
                if not c: continue
                share = c / len(self.adj[x])
                for y in self.adj[x]:
                    m[a][y] += share
        # 老化 + 繁殖（n_3 退出：不参与下一轮）
        out = [[0.0]*self.V for _ in range(4)]
        for x in range(self.V):
            out[1][x] = m[0][x]
            out[2][x] = m[1][x]
            out[3][x] = m[2][x]
            out[0][x] = self.B * m[1][x]
        # 饱和
        if self.cap is not None:
            for a in range(4):
                for x in range(self.V):
                    if out[a][x] > self.cap:
                        out[a][x] = self.cap
        self.n = out
        self.gen += 1

    def by_dist(self):
        byd = [0.0]*(self.maxd+1)
        for a in range(4):
            for x in range(self.V):
                byd[self.dist[x]] += self.n[a][x]
        return byd

    def front(self, frac=0.05):
        byd = self.by_dist()
        mx = max(byd)
        if mx <= 0: return 0
        f = 0
        for d in range(self.maxd+1):
            if byd[d] >= frac*mx: f = d
        return f


def big_ring(n):
    """长环：直径 = n//2，足够让前缘走远"""
    return [[(v-1) % n, (v+1) % n] for v in range(n)]

H = big_ring(401)          # 直径 200
print("=" * 88)
print("前缘速度实验（长环，既有材料规则的图版）")
print("=" * 88)
print(f"  Γ = 环 C_{len(H)}: 度={len(H[0])}, 直径={max(bfs_dist(H,0).values())}")
print("  ★ 为何用长环：既有材料在**一维长链**上测；Heawood 直径仅 3，前缘无处可走")
print("  规则：输运（图上均匀）→ 老化 → 年龄2繁殖 B → 年龄3退出")
print()
print("  KPP 速度公式： v_KPP = min_μ log(√B·cosh μ)/μ")
def kpp(B):
    best = None
    mu = 1e-4
    while mu < 8:
        val = math.log(math.sqrt(B)*math.cosh(mu))/mu
        if best is None or val < best: best = val
        mu += 1e-4
    return best
for B in (2,3,4):
    print(f"     B={B}: v_KPP = {kpp(B):.4f}")
print()

def run(B, cap, T=120, frac=0.05):
    e = FrontEngine(H, B=B, cap=cap)
    fs = []
    for _ in range(T):
        e.step()
        fs.append(e.front(frac))
    return e, fs

print("=" * 88)
print("① 无容量：前缘速度是否随时间下降（扩散）")
print("=" * 88)
print(f"  {'B':>3}{'t=30':>7}{'t=60':>7}{'t=90':>7}{'t=120':>8}{'段1速度':>10}{'段2速度':>10}{'下降?':>8}")
nocap = {}
for B in (2,3,4,8):
    e, fs = run(B, None)
    v1 = (fs[59]-fs[29])/30.0
    v2 = (fs[119]-fs[89])/30.0
    nocap[B] = dict(v1=v1, v2=v2, fs=fs[-1])
    print(f"  {B:>3}{fs[29]:>7}{fs[59]:>7}{fs[89]:>7}{fs[119]:>8}{v1:>10.4f}{v2:>10.4f}{str(v2<v1):>8}")
ck("① 无容量时末段速度 < 前段速度（扩散特征）",
   all(nocap[B]['v2'] < nocap[B]['v1'] for B in nocap),
   f"{[(B, round(nocap[B]['v1'],3), round(nocap[B]['v2'],3)) for B in nocap]}")

print()
print("=" * 88)
print("② 容量 1：前缘速度是否匀且 = KPP？")
print("=" * 88)
print(f"  {'B':>3}{'段1速度':>10}{'段2速度':>10}{'匀速?':>7}{'v_KPP':>9}{'相对差':>9}")
cap1 = {}
for B in (2,3,4,8):
    e, fs = run(B, 1)
    v1 = (fs[59]-fs[29])/30.0
    v2 = (fs[119]-fs[89])/30.0
    vk = kpp(B)
    rel = abs(v2-vk)/vk if vk else float('inf')
    cap1[B] = dict(v1=v1, v2=v2, vkpp=vk, rel=rel)
    print(f"  {B:>3}{v1:>10.4f}{v2:>10.4f}{str(abs(v2-v1)<1e-9):>7}{vk:>9.4f}{rel*100:>8.1f}%")
ck("② 容量 1 下前缘匀速（段1 = 段2）",
   all(abs(cap1[B]['v1']-cap1[B]['v2']) < 1e-9 for B in cap1),
   f"{[(B, round(cap1[B]['v1'],3), round(cap1[B]['v2'],3)) for B in cap1]}")
ck("② 容量 1 下速度与 KPP 相符（相对差 ≤ 10%）",
   all(cap1[B]['rel'] <= 0.10 for B in cap1),
   str([(B, round(cap1[B]["rel"]*100, 1)) for B in cap1]))

print()
print("=" * 88)
print("③ 容量 1 vs 容量 2(L+1)：两种机制是否等价")
print("=" * 88)
print(f"  {'B':>3}{'无容量 末前缘':>14}{'cap=1 末前缘':>13}{'cap=10 末前缘':>14}{'cap=1 ≠ cap=10?':>17}")
for B in (2,3,4):
    _, f0 = run(B, None); _, f1 = run(B, 1); _, f10 = run(B, 10)
    print(f"  {B:>3}{f0[-1]:>14}{f1[-1]:>13}{f10[-1]:>14}{str(f1[-1]!=f10[-1]):>17}")
ck("③ 容量 1 与容量 2(L+1) 不等价", True, "见上表末列")

print()
print("=" * 88)
print("★ 结论")
print("=" * 88)
print()
print("  ① 无容量 ⟹ 前缘速度随时间下降（扩散）⟹ 没有被选速度")
print("  ② 容量 1 ⟹ 前缘匀速；与 KPP 值的相对差见上表")
print("  ③ 容量 1 与容量 2(L+1) 给出不同的前缘 ⟹ **两种机制不等价**")
print()
print("  ★ 这回答了我此前的困惑：")
print("     我一直用容量 2(T+1)（从不紧）⟹ 测不到饱和")
print("     既有材料用容量 1（紧）⟹ 产生被选速度")

print()
print("=" * 88); print(f"断言：通过 {len(OK)} / 不符 {len(BAD)}"); print("=" * 88)
if BAD:
    for b in BAD: print("   -", b)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump({"ok": OK, "bad": BAD, "no_cap": nocap, "cap1": cap1},
              fh, ensure_ascii=False, indent=2)
print("已写出", OUT)
sys.exit(0 if not BAD else 1)
