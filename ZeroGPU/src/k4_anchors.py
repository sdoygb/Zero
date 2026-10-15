# -*- coding: utf-8 -*-
"""
ZeroGPU · 量子特征锚定原型：从 K₄ 结构提取模式的 (自旋, 色) 模拟量

不预设"三代轻子"。用"锚"对照粒子种类：从 K₄/双覆盖结构提取每个模式的
量子特征（自旋、色），与 PDG 粒子对照。

锚定来源：
  自旋(费米/玻色)：模式在双覆盖结构(8顶点 spinor)下的变换性质
     - 单覆盖(K₄ 4顶)特征向量 C₃ 下整数表示 → 玻色子候选
     - 双覆盖(8顶) → spinor(半整数自旋) → 费米子候选
  色(轻子/夸克)：模式在 C₃(3轴旋转) 下的表示
     - C₃ 平凡表示(不变) → 色单态 → 轻子候选
     - C₃ 复表示(ω/ω²) → 色三重态 → 夸克候选
  代：3 频率槽（频率序）
  电荷：U(1) 相（Z086 输入，本脚本标记为输入位）

C₃ 旋转：四面体 3 重轴（固定顶点 0，循环置换顶点 1→2→3→1）。
用法：/usr/local/bin/python3 k4_anchors.py
"""
import numpy as np

def K4():
    return np.array([[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]], float)

def C3_matrix():
    """固定顶点0，循环 1→2→3→1 的置换矩阵"""
    P = np.zeros((4,4))
    P[0,0]=1; P[2,1]=1; P[3,2]=1; P[1,3]=1
    return P

def k4_double_cover():
    A = np.zeros((8,8), float)
    for v in range(4):
        A[2*v,2*v+1] = A[2*v+1,2*v] = 1.0
    for u in range(4):
        for v in range(4):
            if u!=v:
                A[2*u,2*v] = A[2*u+1,2*v+1] = 1.0
    return A

def C3_double(P4):
    """把 C3 置换提升到双覆盖（覆盖对同作用）"""
    P = np.zeros((8,8))
    for v in range(4):
        w = int(np.argmax(P4[v]))     # v -> w
        P[2*v,2*w]=1; P[2*v+1,2*w+1]=1
    return P

print("=" * 86)
print("量子特征锚定原型：K₄ 结构的模式变换性质 → 自旋 / 色")
print("=" * 86)

A = K4(); P = C3_matrix()
w, V = np.linalg.eigh(A)
print(f"\n[K₄ 单覆盖 4 顶点] 邻接谱: {sorted(np.round(w,4))}")
# 每个特征向量在 C3 下的变换性质
print("  特征模式在 C₃ 旋转下的分类：")
for i in range(4):
    v = V[:,i]
    # 作用 C3 三次，看是否回到自身（C3 表示特征值）
    pv = P @ v
    ratio = pv / v if np.max(np.abs(v))>1e-9 else None
    # 判断：C3 特征值 = 1(平凡) / ω,ω²(复=三重/旋量)
    p3 = np.linalg.matrix_power(P,3) @ v
    cyc3 = np.allclose(p3, v)
    p1 = np.allclose(P@v, v)
    print(f"    λ={w[i]:.3f} 模式{i}: C₃-不动(p1)={p1}  C₃³循环回自身={cyc3}")

# C3 作用在 -1 特征空间上的表示
print("\n  C₃ 作用在 -1 特征空间（3 维）的表示分解：")
neg_idx = [i for i in range(4) if abs(w[i]+1)<1e-6]
if neg_idx:
    E = V[:,neg_idx]  # 3x3 基
    # P 投影到 E 空间
    PE = E.T @ P @ E
    e_c3 = np.linalg.eigvals(PE)
    print(f"    C₃ 表示特征值 = {sorted(np.round(e_c3,3))}")
    print("    → 含复单位根(ω/ω²) 则=色三重/旋量；全±1 则=整数表示(玻色)")

print("\n" + "=" * 86)
print("[K₄ 双覆盖 8 顶点] spinor 结构：费米子(半整数自旋)候选")
print("=" * 86)
B = k4_double_cover(); P8 = C3_double(P)
wB, VB = np.linalg.eigh(B)
print(f"  双覆盖谱: {sorted(np.round(wB,4))}")
print(f"  → 双覆盖给出半整数自旋(spinor)承载；模式=费米子候选")
print(f"  → 8 顶点 = 色 SU(3) 八重态候选（费米子按色 8 变换）")

print("\n" + "=" * 86)
print("锚定对照（提案，非结论）：")
print("  - 双覆盖 spinor 模式 + 色三重(C₃复表示) → 夸克候选（带色荷）")
print("  - 双覆盖 spinor 模式 + 色单态(C₃平凡)  → 轻子候选（无色）")
print("  - 单覆盖整数表示 + 色单态                → 玻色子候选（自旋 0/1）")
print("  - 代 = 3 频率槽序号；电荷 = U(1) 相（Z086 输入位）")
print("=" * 86)
