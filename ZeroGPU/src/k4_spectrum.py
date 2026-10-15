# -*- coding: utf-8 -*-
"""
ZeroGPU · 质量比内生来源候选①：图结构多尺度谱

Z089 定义：粒子质量 = 稳态线性化算子的激发模式本征值（质量=频率）。
本脚本检验一个内生机理候选——质量不来自 seed 相邻比，而来自 Zero 承载图
（K₄ 及其复合/覆盖结构）自身的谱（特征值 = 模流速率 = 频率）。

若某个 Zero 图结构的谱给出 3 个悬殊的特征值（时间尺度 1/|λ|），则 3 代
质量比是内生的（不依赖 seed）。检验 K₄ 及复合结构的谱。

候选结构：
  ① K₄（Γ 原生，4 顶）
  ② K₄ 双覆盖（8 顶，girth=4）
  ③ K₄ 笛卡尔积 K₄□K₄（16 顶）
  ④ 超立方体 Q3(8)/Q4(16)
  ⑤ K₄ 幂（邻接幂 K₄^2：步长-2 连通）—— 谱幂放大

质量映射候选：m_i ∝ 1/|λ_i|（模流时间）或 m_i ∝ |λ_i|（模式能量）。
三代 = 3 个最大/最显著模式。检验谱比是否极端（观测轻子 1:207:3477）。

用法：/usr/local/bin/python3 k4_spectrum.py
"""
import numpy as np

def eigvals(A, dec=4):
    return sorted(round(float(x), dec) for x in np.linalg.eigvalsh(A))

def time_scale(eigs):
    """时间尺度 ∝ 1/|λ|（非零）；质量 ∝ 模流速率 = |λ|"""
    nz = [abs(e) for e in eigs if abs(e) > 1e-9]
    return sorted(set(round(1.0 / x, 4) for x in nz))

def K4():
    return np.array([[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]], float)

def k4_double_cover():
    """K₄ 双覆盖：8 顶，每原顶点复制为 ±，边跨覆盖对 + 覆盖内按 K₄ 结构"""
    A = np.zeros((8,8), float)
    for v in range(4):          # 覆盖对 (2v, 2v+1)
        A[2*v, 2*v+1] = A[2*v+1, 2*v] = 1.0
    for u in range(4):
        for v in range(4):
            if u != v:
                A[2*u, 2*v] = A[2*u+1, 2*v+1] = 1.0  # 同符号跨覆盖
    return A

def cartesian(A, B):
    n, m = A.shape[0], B.shape[0]
    C = np.zeros((n*m, n*m), float)
    for i in range(n):
        for j in range(m):
            for k in range(n):
                if A[i,k]:
                    C[i*m+j, k*m+j] = 1.0
            for l in range(m):
                if B[j,l]:
                    C[i*m+j, i*m+l] = 1.0
    return C

def hypercube(d):
    n = 2**d
    A = np.zeros((n,n), float)
    for i in range(n):
        for b in range(d):
            A[i, i ^ (1<<b)] = 1.0
    return A

print("=" * 84)
print("质量比内生来源候选①：图结构多尺度谱（特征值=模流速率=频率）")
print("=" * 84)
print("观测三代（轻子）质量比 = 1 : 207 : 3477（e:μ:τ），Koide Q=2/3")

graphs = {
    "K₄(4)":           K4(),
    "K₄双覆盖(8)":     k4_double_cover(),
    "K₄□K₄笛卡尔(16)": cartesian(K4(), K4()),
    "超立方Q3(8)":     hypercube(3),
    "超立方Q4(16)":    hypercube(4),
}
for name, A in graphs.items():
    ev = eigvals(A)
    nz = [e for e in ev if abs(e) > 1e-9]
    maxr = max(abs(e) for e in nz) / min(abs(e) for e in nz) if nz else 0
    print(f"\n  {name}: 谱={ev}")
    print(f"      非零|λ|={[abs(round(e,4)) for e in nz]}")
    print(f"      最大/最小|λ| 比={maxr:.3f}  时间尺度(1/|λ|)={time_scale(ev)}")
    print(f"      → 正则图谱为整数/小比，不提供悬殊多尺度")

print("\n" + "=" * 84)
print("结论：图自身谱（正则图）为整数或小整数比，不给出 1:207:3477 的极端")
print("质量比。质量比极端性需非线性层间动力学（回灌不动点），非图谱。")
print("=" * 84)
