#!/usr/bin/env python3
"""
Z002 §6.7／§6.8 维数账的算术核验（不依赖 Lean，独立复算）。

口径（与 §6.7 的定义句一致，B 取 1-骨架的关联映射，边流量 ↦ 顶点上的边界）：
    ker B = 圈空间 Z_1，dim = |E| - |V| + 1（连通图，欧拉公式）
    im B  = 零和超平面 H_Q，dim = |V| - 1
"闭合＝运动"读法比的是这两个数；相等**仅当** |E| = 2|V| - 2。
"""
from itertools import combinations

def dims(nv, ne):
    """(圈空间维数, 零和维数)——要求图连通。"""
    return ne - nv + 1, nv - 1

# ---- ① K_m：m = 2..12 ----
print("① 完全图 K_m（|V|=m, |E|=C(m,2)）")
print(f"{'m':>3} | {'|E|':>4} | {'dim Z_1':>7} | {'dim H_Q':>7} | 相等?")
for m in range(2, 13):
    z, h = dims(m, m * (m - 1) // 2)
    print(f"{m:>3} | {m*(m-1)//2:>4} | {z:>7} | {h:>7} | {'✔ 相等' if z == h else '否'}")
    assert (z == h) == (m == 4), f"m={m}"

# ---- ② §6.8 的五个实例（|E| 逐个数出来）----
instances = [
    ("K_4",              4,  6),   # K4: C(4,2)=6
    ("C_4",              4,  4),   # 4-cycle
    ("4x4 方格",         16, 24),  # 4*3*2 = 24
    ("Heawood",          14, 21),  # 3-正则、围长 6
    ("K_4 口 C_4 (16槽)", 16, 40),  # K4 方向 6*4 + C4 方向 16
]
print("\n② §6.8 的实例")
print(f"{'图':>18} | {'|V|':>4} | {'|E|':>4} | {'dim Z_1':>7} | {'dim H_Q':>7} | §6.8 写的")
doc = {"K_4": "3=3", "C_4": "1≠3", "4x4 方格": "9≠15", "Heawood": "8≠13",
       "K_4 口 C_4 (16槽)": "25≠15"}
for name, nv, ne in instances:
    z, h = dims(nv, ne)
    got = f"{z}≠{h}" if z != h else f"{z}={h}"
    print(f"{name:>18} | {nv:>4} | {ne:>4} | {z:>7} | {h:>7} | {doc[name]:>10}"
          + ("  ✔ 与文档一致" if got == doc[name] else f"  ✘ 文档={doc[name]} 算={got}"))
    assert got == doc[name], name
    assert (z == h) == (name == "K_4"), name

# ---- ③ 判别式的两个等价算术形式 ----
print("\n③ 判别式")
for m in range(2, 13):
    a = (m * (m - 1) // 2 - (m - 1) == m - 1)          # 圈空间 = 零和
    b = (m * (m - 1) // 2 == 2 * (m - 1))              # |E| = 2|V| - 2
    c = ((m - 2) * (m - 3) // 2 == 1) if m >= 3 else False   # 配对唯一性 §6.8(一)
    if a or b or c:
        print(f"  m={m}: 圈空间=零和 {a} ; |E|=2|V|-2 {b} ; 每边不相邻边数=1 {c}")
    assert a == b == (m == 4) or (not a and not b), m
print("\n全部断言通过。")
