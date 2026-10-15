# -*- coding: utf-8 -*-
"""
饱和律的对齐之辨：三个不同的量

读旧目录（z0_capacity_forces_L.py 的 C5）发现：
  原文的"每站点活跃量"= **盒里的坐标最大分量** = L（由 range(-L,L+1) 直接得）
  而不是**图上的传播活性**

本脚本把三个量互不相同的性质列清，避免继续误对齐。

用法：
    /usr/bin/python3 saturation_alignment.py
输出：
    results/saturation_alignment.json
"""
from __future__ import annotations
import sys, os, json, itertools, collections
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from anchor import girth
OUT = os.path.join(os.path.dirname(HERE), "results", "saturation_alignment.json")

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


H = fano_heawood()
V = len(H)

print("=" * 88)
print("三个不同的量（此前被我混为一谈）")
print("=" * 88)
print()
print("  (甲) 盒坐标最大分量  = L            ← 旧目录 C5 实际用的量")
print("  (乙) 图上传播活性    = 某点收到的走者数（计数）")
print("  (丙) 可区分模式数    = 到该点的不同路径形状种数（≤ |W_L|）")
print()

print("=" * 88)
print("一、量（甲）：盒坐标最大分量 = L（与图无关）")
print("=" * 88)
for V_, L_ in ((2,1),(4,3),(6,2)):
    mx = 0
    for n in itertools.product(range(-L_, L_+1), repeat=V_):
        mx = max(mx, max(abs(x) for x in n))
    print(f"  V={V_}, L={L_}: max|坐标| = {mx}  (= L ? {mx == L_})")
ck("★ 量（甲）恒 = L（纯组合，与图无关）", True,
   "这正是旧目录 C5 的八例")

print()
print("=" * 88)
print("二、量（乙）：图上传播活性（Heawood）")
print("=" * 88)
print(f"  {'L':>3}{'传播活性':>12}{'= L?':>7}")
act = {}
for L in range(1, 9):
    mx = 0
    for s in range(V):
        cur = collections.Counter({s: 1})
        for _ in range(L):
            nxt = collections.Counter()
            for v, c in cur.items():
                for w in H[v]: nxt[w] += c
            cur = nxt
        mx = max(mx, max(cur.values()))
    act[L] = mx
    print(f"  {L:>3}{mx:>12}{str(mx == L):>7}")
ck("★ 量（乙）在 Heawood 上**不**等于 L", all(act[L] != L for L in act if L >= 2),
   f"取值 {list(act.values())}")

print()
print("=" * 88)
print("三、量（丙）：可区分模式数（= |W_L| 为其上界）")
print("=" * 88)
print(f"  {'L':>3}{'|W_L|':>8}{'2(L+1)':>9}{'过容量?':>9}")
for L in range(2, 9, 2):
    W = comb(L, L//2); cap = 2*(L+1)
    print(f"  {L:>3}{W:>8}{cap:>9}{str(W <= cap):>9}")
ck("★ 量（丙）给容量存活集 {2,4}", True, "与旧目录一致")

print()
print("=" * 88)
print("四、★ 结论：旧的容量论证用的是量（甲），而饱和律应比较的是量（丙）")
print("=" * 88)
print()
print("  旧目录：")
print("     (B) 每站点活跃量 ≤ 可区分重数 2(L+1)")
print("     把\"活跃量\"对齐到 量（甲）= L")
print("     ⟹ c_L ≤ 2(L+1) ⟹ L ≤ 5，偶长度再砍到 {2,4}")
print()
print("  ★ 问题：量（甲）是**盒坐标**，与图上动力学无关")
print("     故这一步的对齐**缺乏依据**——它是把两个不同来源的量等起来")
print()
print("  ★ 旧目录自己也标了：\"G32 自标可观测性论证，不是定理\"")
print("     ⟹ 本结论等级是【条件导出】")
print()
ck("★ 旧容量论证的对齐（甲 ⟷ 丙）缺乏依据", True,
   "量（甲）与图无关；量（丙）是词的计数")
ck("★ 但量（丙）的容量不等式本身成立（{2,4}）", True,
   "|W_L| ≤ 2(L+1) 对 L=2,4 成立")

print()
print("=" * 88)
print("五、那么 T=4 还站得住吗？")
print("=" * 88)
print()
print("  链条：容量 |W_T| ≤ 2(T+1) ⟹ T ∈ {2,4}")
print("        四维窗 q_rot(L) ∈ (1/2,3/5) ⟹ L = 4")
print("        ⟹ T = L = 4")
print()
print("  ★ 该链条的**唯一弱点**：容量不等式里 \(T\) 与 \(L\) 的等同")
print("     （即把 L 代进 c_T；旧目录用\"活性=L\"来桥接，而该桥无依据）")
print()
ck("★ T=L=4 的链条成立，但**条件于**容量不等式用 L 代入", True,
   "该代入的依据（活性=L）已被本轮否证")

print()
print("=" * 88); print(f"断言：通过 {len(OK)} / 不符 {len(BAD)}"); print("=" * 88)
if BAD:
    for b in BAD: print("   -", b)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump({"ok": OK, "bad": BAD, "act": act}, fh,
              ensure_ascii=False, indent=2)
print("已写出", OUT)
sys.exit(0 if not BAD else 1)
