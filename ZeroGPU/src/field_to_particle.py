#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
field_to_particle.py —— 场→粒子的结构性连接核验
=================================================
假说：场的 S₃ 对称性（F=1⊕2，壁处 S₃→S₂ 破缺）内禀地含 Z₃ 子群，
     而 S₃ 的 2 维不可约表示限制到 Z₃ 分裂为 ω⊕ω² —— 这正是 Koide
     三代质量形（Z₃ 正则表示）的"非迹部分"。
目的：精确群论核验该连接成立与否，并划清"结构性来源（内生）"与
     "数值参数（具名输入）"的边界。
核验：
  A. S₃ 的特征标表 + 类分解
  B. χ₂|_{{Z₃}} = ω ⊕ ω²（2 维不可约 → Z₃ 分裂为共轭一维特征）
  C. 场 F=1⊕2，方向部分 = ω⊕ω² = 三代相位（Koide 结构）来源
  D. 与 K₄ 色 C₃ 分解（Z090）对齐："三"在 Zero 两处独立内生
  E. 诚实边界：√2 振幅、δ=2/9 不内生（具名输入）
"""
import json
import numpy as np

OUT = "/Users/oygb/Downloads/Zero/ZeroGPU/results/field_to_particle.json"

# ---------------------------------------------------------------
# A. S₃ 特征标表（精确，符号解）
# ---------------------------------------------------------------
# S₃ 三个共轭类：{e} 阶1, {3-cycles} 阶2, {transpositions} 阶3
# 类的大小：1, 2, 3（S₃ 的类 1,2,3）
# 特征标（按类 [e, 3-cyc, trans]）：
chi_1   = np.array([1, 1, 1])    # 平凡
chi_sgn = np.array([1, 1, -1])   # 符号
chi_2   = np.array([2, -1, 0])   # 2 维不可约（旋量，dim2）
classes = [1, 2, 3]              # 各类大小
order   = 6

def inner(a, b):
    return sum(classes[i]*a[i]*np.conj(b[i]) for i in range(3)) / order

rep_table = {
    "chi_1":  chi_1, "chi_sgn": chi_sgn, "chi_2": chi_2
}
print("=== A. S₃ 特征标内积（正交性）===")
for k, v in rep_table.items():
    for k2, v2 in rep_table.items():
        print(f"  <{k},{k2}> = {inner(v, v2):+.4f}")

# 场的 3 维置换表示：χ_F(σ) = #Fix(σ)
chi_F = np.array([3, 0, 1])   # e→fix3, 3-cyc→fix0, trans→fix1
print(f"  χ_F = {chi_F}, <χ_F,1>={inner(chi_F,chi_1):.3f}, "
      f"<χ_F,sgn>={inner(chi_F,chi_sgn):.3f}, <χ_F,2>={inner(chi_F,chi_2):.3f}")
print("  ⇒  F = 1 ⊕ 2   (3 = 1 + 2)")

# ---------------------------------------------------------------
# B. χ₂ 限制到 Z₃ ⊂ S₃
# ---------------------------------------------------------------
# Z₃ = {e, (123), (132)} 3 阶循环子群，只含类 {e} 与 {3-cycles}
# χ₂ 在 Z₃ 上：χ₂(e)=2, χ₂(3-cyc)=-1（两个 3-cycle 同特征标）
chi2_Z3 = np.array([2, -1, -1])          # 在 Z₃ 的 [e, a, a²] 上
# Z₃ 的 3 个一维特征（在 [e,a,a²] 上）：
z1   = np.array([1, 1, 1])               # 平凡
z_om = np.array([1, np.exp(2j*np.pi/3), np.exp(4j*np.pi/3)])  # ω
z_om2= np.array([1, np.exp(4j*np.pi/3), np.exp(2j*np.pi/3)])  # ω²

def inner3(a, b):
    return (a[0]*np.conj(b[0]) + a[1]*np.conj(b[1]) + a[2]*np.conj(b[2])) / 3

print("\n=== B. χ₂|_{{Z₃}} 分裂（3 阶循环子群内）===")
c1  = inner3(chi2_Z3, z1)
com = inner3(chi2_Z3, z_om)
com2= inner3(chi2_Z3, z_om2)
print(f"  ⟨χ₂,1⟩={c1:+.3f}   ⟨χ₂,ω⟩={com:+.3f}   ⟨χ₂,ω²⟩={com2:+.3f}")
print(f"  ⇒  χ₂|_{{Z₃}} = ω ⊕ ω²   (共轭一维特征, 实系数)")

# ---------------------------------------------------------------
# C. 场的方向部分 = 三代相位（Koide Z₃ 结构）
# ---------------------------------------------------------------
print("\n=== C. 场→粒子 结构对应 ===")
print("  场 F = 1 ⊕ 2 ：1=强度(迹), 2=方向(无迹)")
print("  2 维方向部分限制到 Z₃ = ω ⊕ ω²")
print("  Koide 三代质量形用 Z₃ 正则表示 1 ⊕ ω ⊕ ω² (Z014)")
print("  ⇒ 场的方向信息(ω⊕ω²) = 三代质量形的非迹部分")

# ---------------------------------------------------------------
# D. 与 K₄ 色 C₃ 分解对齐（Z090）
# ---------------------------------------------------------------
print("\n=== D. '三'在 Zero 两处独立内生 ===")
print("  粒子侧: K₄ 单覆盖谱 {3,-1,-1,-1}, -1 空间在 C₃ 下 = {1,ω,ω²} (Z090)")
print("  场  侧: S₃ 方向表示 2, 限制到 Z₃ = {ω,ω²} (本文)")
print("  两者都含 3 阶循环对称 ⇒ 三代是 Zero 内生的循环结构")

# ---------------------------------------------------------------
# E. 诚实边界：√2 振幅、δ 角不内生
# ---------------------------------------------------------------
print("\n=== E. 诚实边界 ===")
print("  Koide 数值参数: √2 振幅、δ=2/9 不从 S₃/Z₃ 群论给出")
print("  ⇒ 结构性'三'(Z₃)内生; 数值质量比仍具名输入(U(1), Z086)")

result = {
    "s3_inner": {f"{k},{k2}": inner(v, v2).real for k, v in rep_table.items() for k2, v2 in rep_table.items()},
    "chi_F": chi_F.tolist(),
    "chi_F_decomp": "1 ⊕ 2",
    "chi2_restrict_Z3": {
        "⟨χ2,1⟩": c1.real, "⟨χ2,ω⟩": com.real, "⟨χ2,ω²⟩": com2.real,
        "decomp": "ω ⊕ ω²"
    },
    "field_to_particle": "2维方向 = Z₃ 三代非迹部分 (结构对应)",
    "three_is_endogenous": "粒子侧 C₃(K₄) ∥ 场侧 Z₃(S₃): 3 阶循环内生",
    "boundary": "√2振幅/δ=2/9 数值不内生, 仍具名输入(Z086)"
}
with open(OUT, "w") as f:
    json.dump(result, f, ensure_ascii=False, indent=2)
print(f"\n结果已写入 {OUT}")
