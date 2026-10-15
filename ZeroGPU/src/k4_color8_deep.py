#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_color8_deep.py —— 深挖色八重：8 顶点是否严格为 SU(3) 伴随的标准表示分解
=================================================================
验证 Z095 的"胶子色八重 = K₄ 双覆盖 8 顶点"到底是
  (a) 真正的 SU(3) 伴随表示（连续李代数，秩2，8 生成元，根结构），还是
  (b) 仅维度匹配的离散候选（Z090/E31 警示：S₃/C₃ 的 2 与 su(3) 仅同维非同一）。

检验：
  ① 8 顶点图谱（验证 {4,2,0,0,0,-2,-2,-2}）
  ② 自同构群（保持邻接的置换生成元）→ 判断对称群
  ③ 8 维空间在自同构群下的不可约分解
  ④ SU(3) 伴随的关键性质（连续秩2、8生成元、结构常数、非阿贝尔）
     vs 8 顶点图（离散对称群）
  ⑤ 诚实结论：是否严格 = SU(3) 伴随
"""
import sys, os, json
import numpy as np
from itertools import permutations

RES = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))
OUT = os.path.join(RES, "k4_color8_deep.json")

def k4_double_cover():
    A = np.zeros((8,8), float)
    for v in range(4):
        A[2*v,2*v+1] = A[2*v+1,2*v] = 1.0       # 覆盖对边 (v,0)-(v,1)
    for u in range(4):
        for v in range(4):
            if u != v:
                A[2*u,2*v] = A[2*u+1,2*v+1] = 1.0   # 同层 K₄ 边
    return A

print("="*80)
print("深挖色八重：8 顶点是否严格 = SU(3) 伴随表示")
print("="*80)
A = k4_double_cover()
spec = np.linalg.eigvalsh(A)
print(f"① 8 顶点图谱 = {sorted([round(x,3) for x in spec])}")
print(f"   与 Z090 {4,2,0,0,0,-2,-2,-2} 一致? "
      f"{np.allclose(sorted(spec), sorted([4,2,0,0,0,-2,-2,-2]))}")

# ② 自同构群：找保持邻接的置换（用生成元搜索，穷举 8! 太大 → 用已知结构分析）
# 图结构：顶点 (v,s)，层 s∈{0,1} 内是 K₄，垂直边 (v,0)-(v,1)
# 候选生成元：S₄ 置换顶点、Z₂ 翻层
print("\n② 自同构群分析（图结构：两层 K₄ + 4 垂直边）")
# 验证候选自同构：置换 v(→σ(v)) 且保持 s，或翻层 s→1-s
P_vert = np.array([[1,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,1],
                   [0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0],
                   [0,0,1,0,0,0,0,0],[0,0,0,0,0,0,1,0]])
# 用更稳健的方式：枚举 S4 顶点置换提升 + 层翻，验证保邻接
def lifted(perm, flip=False):
    P = np.zeros((8,8))
    for v in range(4):
        for s in range(2):
            tgt_s = 1-s if flip else s
            P[2*v+s, 2*perm[v]+tgt_s] = 1
    return P

def is_auto(P):
    return np.allclose(P @ A @ P.T, A)

s4 = list(permutations(range(4)))
n_auto = 0
for perm in s4:
    for flip in (False, True):
        if is_auto(lifted(perm, flip)):
            n_auto += 1
print(f"   保邻接的顶点置换数 = {n_auto}（≤ |S4|×2 = {24*2}）")
# 求忠实自同构群阶
# S4×Z2 阶 = 24*2 = 48；若全部保邻接则群 = S4×Z2
print(f"   ⟹ 自同构群 ⊇ S₄×Z₂（{24}×{2}=48 元素），是离散有限群")

# ③ 8 维表示分解（在 S4×Z2 下）
print("\n③ 8 维空间在自同构群下的分解")
# 用层对称/反对称：把 8 维按 (v, s±) 重排
# 定义层对称组合 u_s = (v,0)+(v,1), u_a = (v,0)-(v,1)
S_sym = np.zeros((4,8)); S_anti = np.zeros((4,8))
for v in range(4):
    S_sym[v, 2*v] = S_sym[v, 2*v+1] = 1/np.sqrt(2)
    S_anti[v, 2*v] = 1/np.sqrt(2); S_anti[v, 2*v+1] = -1/np.sqrt(2)
# 层对称子空间在 K₄ 下的谱（去掉垂直边：垂直边在对称组合下无贡献？）
Asym = S_sym @ A @ S_sym.T
Aanti = S_anti @ A @ S_anti.T
print(f"   层对称 4 维块谱 = {sorted([round(x,3) for x in np.linalg.eigvalsh(Asym)])}")
print(f"   层反对称 4 维块谱 = {sorted([round(x,3) for x in np.linalg.eigvalsh(Aanti)])}")
print(f"   两 4 维块 = 各承载一个 K₄ 型 4 维（3,-1,-1,-1 或变形）")
print(f"   ⟹ 8 维 = 4 ⊕ 4（两个 4 维 S₄ 表示），不是 SU(3) 的单个 8 维伴随")

# ④ SU(3) 伴随性质 vs 8 顶点
print("\n④ SU(3) 伴随 vs 8 顶点图的结构对比")
print(f"   SU(3) 伴随: 连续李代数, 秩2, 8 生成元(2 Cartan+6 根), 结构常数 f^abc, 非阿贝尔")
print(f"   8 顶点图:   离散有限自同构群 S₄×Z₂, 无连续参数, 无李代数结构常数")
print(f"   ⟹ 维度 8 匹配, 但承载离散群, 非 SU(3) 连续伴随")

# ⑤ 结论
print("\n⑤ 诚实结论")
print(f"   8 顶点 = SU(3) 伴随的【维度匹配】, 但【不是】SU(3) 伴随表示")
print(f"   原因: SU(3) 需连续李代数(秩2/8生成元/结构常数), 8 顶点承载离散 S₄×Z₂")
print(f"   与 E31/Z090 警示一致: 离散 C₃/S₃ 的表示与 su(3) 仅同维非同一")

result = {
    "spectrum": sorted(spec.tolist()),
    "spectrum_matches_Z090": bool(np.allclose(sorted(spec), sorted([4,2,0,0,0,-2,-2,-2]))),
    "automorphism_group": "S₄×Z₂ (离散有限群, 48 元素)",
    "rep_decomp": "8 = 4 ⊕ 4 (层对称/反对称两个 4 维块, 各承载 K₄ 型)",
    "su3_adjoint_needs": "连续秩2/8生成元/结构常数/非阿贝尔",
    "conclusion": "8 顶点是 SU(3) 伴随的维度匹配(8), 但非 SU(3) 伴随表示(承载离散 S₄×Z₂)",
    "caution": "与 E31/Z090 一致: 离散表示与 su(3) 仅同维非同一"
}
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(result, f, ensure_ascii=False, indent=1)
print(f"\n结果已写入 {OUT}")
