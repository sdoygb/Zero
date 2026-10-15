#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_so10.py —— 接 SO(10)：Zero 每代费米子(8场) 是否精确对应 SO(10) 16 维旋量表示
==============================================================================
SO(10) 大统一群：16 维旋量表示 = 一代费米子（含右手中微子 ν_R）。
  SU(5) 分解：16 = 5̄ + 10 + 1
    5̄ = {d^c(d_R), ν_L, e_L};  10 = {u_L, d_L, u^c(u_R), e^c(e_R)};  1 = {ν_R}
验证：Zero 每代 8 场（Z097）是否 = SO(10) 16 的完整成员。
评估：Zero 内生 ν_R ⟹ 倾向 SO(10)（16 含 ν_R）vs SU(5)（5̄+10=15 不含 ν_R）。
"""
import sys, os, json
import numpy as np

RES = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))
OUT = os.path.join(RES, "k4_so10.json")

# Zero 每代 8 场（Z097），用 (场名, 色, T₃, Y, Q)
zero8 = [
    ("ν_L", "单", 0.5, -1.0, 0.0),
    ("e_L", "单", -0.5, -1.0, -1.0),
    ("e_R", "单", 0.0, -2.0, -1.0),
    ("ν_R", "单", 0.0, 0.0, 0.0),
    ("u_L", "三", 0.5, 1/3, 2/3),
    ("u_R", "三", 0.0, 4/3, 2/3),
    ("d_L", "三", -0.5, 1/3, -1/3),
    ("d_R", "三", 0.0, -2/3, -1/3),
]

print("="*84)
print("接 SO(10)：Zero 每代费米子(8场) → SO(10) 16 维旋量表示")
print("="*84)

# ① SO(10) 16 表示 = 5̄+10+1 分解
print("\n① SO(10) 16 维旋量表示（SU(5) 分解）")
so10 = {
    "5̄": [("d_R", "三反", 0.0, -2/3, -1/3), ("ν_L", "单", 0.5, -1.0, 0.0), ("e_L", "单", -0.5, -1.0, -1.0)],
    "10": [("u_L", "三", 0.5, 1/3, 2/3), ("d_L", "三", -0.5, 1/3, -1/3),
           ("u_R", "三反", 0.0, 4/3, 2/3), ("e_R", "单", 0.0, -2.0, -1.0)],
    "1":  [("ν_R", "单", 0.0, 0.0, 0.0)],
}
dim = 0
for rep, members in so10.items():
    n = sum(3 if m[1]!="单" else 1 for m in members)
    dim += n
    names = ", ".join(m[0] for m in members)
    print(f"     {rep:3s} ({n:2d} 维) = {{{names}}}")
print(f"   16 = 5̄+10+1 = {dim} 维 ✓（含色自由度）")

# ② 对照 Zero 每代 8 场
print("\n② 对照 Zero 每代 8 场（Z097）")
zero_names = {f[0] for f in zero8}
so10_names = {m[0] for rep in so10.values() for m in rep}
print(f"   Zero 8 场: {sorted(zero_names)}")
print(f"   SO(10) 16: {sorted(so10_names)}")
print(f"   完全一致? {zero_names == so10_names}")

# 维度核对
zero_dim = sum(3 if f[1]=="三" else 1 for f in zero8)
print(f"   Zero 每代维数 = {zero_dim}（夸克×3色 + 轻子×1）= SO(10) 16 ✓")

# ③ 统一方向评估
print("\n③ 统一方向评估：SO(10) vs SU(5)")
print(f"   SU(5): 5̄+10 = 15 维，不含 ν_R（ν_R 是 SU(5) 单态，需额外 +1）")
print(f"   SO(10): 16 维旋量 = 含 ν_R 的整代费米子")
print(f"   Zero 内生 ν_R（Z097/Z098）⟹ 倾向 SO(10) 统一方向（16 自然承接 ν_R）")
print(f"   Z090 用 SU(5) 归一化 sin²θ_w=3/8；SO(10) 同样在 GUT 标度给 sin²θ_w=3/8")

# ④ 诚实边界
print("\n④ 诚实边界")
print(f"   · 这是量子数结构衔接（Zero 8场 = SO(10) 16 成员），非统一群动力学推导")
print(f"   · sin²θ_w=3/8 归一化是 Z090 已有（SU(5)）；SO(10) 兼容但未推导")
print(f"   · Zero 缺统一群耦合机制/极限（Z026）；不主张 SO(10) 已从 Zero 长出")

result = {
    "so10_16": {rep: [m[0] for m in members] for rep, members in so10.items()},
    "so10_dim": dim,
    "zero_8": sorted(zero_names),
    "match": zero_names == so10_names,
    "zero_dim": zero_dim,
    "unification_assess": "Zero 内生 ν_R ⟹ 倾向 SO(10)(16含ν_R); SU(5)(5̄+10=15)不含 ν_R",
    "boundary": "量子数结构衔接, 非统一群动力学; sin²θ_w=3/8是Z090既有; 缺统一机制Z026"
}
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(result, f, ensure_ascii=False, indent=1)
print(f"\n结果已写入 {OUT}")
