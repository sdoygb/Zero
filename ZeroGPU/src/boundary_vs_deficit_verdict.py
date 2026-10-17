#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
boundary_vs_deficit_verdict.py —— 「界间耦合」与「缺陷角」是不是同一机制？
=========================================================================
背景：`D82` 查明「界间耦合」（paper §5.9 的翘曲乘积交叉项）**无台账来历**；
      `D83`/`D84` 查明「缺陷角」（Regge／词层面积）**有台账支撑**（非正则性已导出；#30/#31 词层实现）。
本脚本**判定两者是否同一机制**（目标要求：否定性结论须 ≥2 独立角度）。

四角度：
  ① 量纲：交叉项的量纲 vs 缺陷角的量纲
  ② 极端情形 1：单因子翘曲乘积 ⟹ 耦合消失，曲率是否仍在？
  ③ 极端情形 2：单个 4-单纯形 ⟹ 缺陷角非零，曲率是否也在？
  ④ 结构类型：逐点张量分量 vs 集中在铰链上的角剩余

运行：cd ZeroGPU/src && /usr/bin/python3 boundary_vs_deficit_verdict.py
"""
import json
import math
import os
import sys

import numpy as np
import sympy as sp

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
CK, OUT = [], []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def note(name, value, detail=""):
    OUT.append({"name": name, "value": value, "detail": detail})
    print(f"      · {name} = {value}" + (f"   {detail}" if detail else ""))


# ============================================================
print("=" * 78)
print("角度 ①：量纲 —— 交叉项 vs 缺陷角")
print("=" * 78)
print("  交叉项 = −(R_a'/R_a)·Σ_{b≠a} n_b (R_b'/R_b)；R'/R 量纲 = 1/长度 ⟹ 交叉项量纲 = 1/长度²")
print("  缺陷角 δ_h = 2π − Σθ  ⟹ **无量纲**（角度）")
lam, s = sp.symbols('lambda s', positive=True)
# 取 R_a = exp(λ_a ρ)，ρ → sρ 等价 λ_a → λ_a/s
cross = -lam * 2 * lam                      # n_b = 2 的两因子情形
cross_scaled = -(-lam / s) * 2 * (-lam / s)  # λ → λ/s
print(f"  数值核对：取 R_a=e^(λ_a ρ)（两因子，n_b=2），交叉项 = {cross}")
print(f"            尺度变换 ρ→sρ（即 λ→λ/s）后 = {sp.simplify(cross_scaled)}")
ratio = sp.simplify(cross_scaled / cross)
print(f"            比值 = {ratio}")
ck("J1 交叉项在 ρ→sρ 下按 1/s² 缩放（量纲 1/长度²）；缺陷角无量纲",
   sp.simplify(ratio - 1 / s**2) == 0,
   f"交叉项 ∝ 1/s²；δ_h 是角度（无量纲）⟹ 量纲不同，不可能是同一个量")

# ============================================================
print()
print("=" * 78)
print("角度 ②：单因子翘曲乘积 —— 耦合消失，曲率是否仍在？")
print("=" * 78)
rho = sp.symbols('rho', positive=True)
R1 = sp.Function('R1')(rho)
# 取 R2' = 0（单因子极限）：耦合项 = -(R1'/R1)·n2·(R2'/R2) = 0
coupling = -(sp.diff(R1, rho) / R1) * 2 * 0
self_term = ((2 - 1) - R1 * sp.diff(R1, rho, 2) - (2 - 1) * sp.diff(R1, rho)**2) / R1**2
print(f"  耦合项（R₂'=0）= {coupling}")
print(f"  自身项 = {sp.simplify(self_term)}")
# 取具体 R1(ρ) = 1 + 0.3ρ + 0.05ρ²，在 ρ=0 处
R1v = 1 + sp.Rational(3, 10) * rho + sp.Rational(5, 100) * rho**2
self_val = float(sp.simplify(self_term.subs(R1, R1v).subs(rho, 0)))
print(f"  取 R₁(ρ)=1+0.3ρ+0.05ρ²，则自身项在 ρ=0 处 = {self_val:.6f}")
ck("J2 单因子时耦合项为零，但纤维 Ricci 仍 ≠ 0（**有曲率而无耦合**）",
   abs(float(coupling)) < 1e-15 and abs(self_val) > 1e-6,
   f"耦合 = 0；自身项 = {self_val:.6f} ⟹ 耦合既非曲率的必要条件")

# ============================================================
print()
print("=" * 78)
print("角度 ③：单个 4-单纯形 —— 缺陷角非零，曲率是否也在？")
print("=" * 78)
theta = math.acos(0.25)
delta = 2 * math.pi - 2 * theta
print(f"  单 4-单纯形：胞间二面角 θ = arccos(1/4) = {math.degrees(theta):.6f}°")
print(f"              每个三角形铰链被 2 个胞共享 ⟹ δ = 2π − 2θ = {math.degrees(delta):.6f}°")
print()
print("  但 5 个顶点嵌在 ℝ⁴（4 维仿射子空间）中 ⟹ 这**一块是平的**（内蕴 Riemann 张量为 0）")
# 数值验证"平坦"：5 顶点的 Gram 矩阵秩 ≤ 4（在 ℝ⁴ 中），且边长满足 Cayley–Menger 相容
V = np.array([(math.cos(2 * math.pi * j / 5), math.sin(2 * math.pi * j / 5),
               math.cos(4 * math.pi * j / 5), math.sin(4 * math.pi * j / 5)) for j in range(5)])
D2 = np.array([[float((V[i] - V[j]) @ (V[i] - V[j])) for j in range(5)] for i in range(5)])
# 由距离矩阵重建内蕴 Gram（取 V[0] 为原点）
G = np.array([[0.5 * (D2[0][i] + D2[0][j] - D2[i][j]) for j in range(5)] for i in range(5)])
ev = np.linalg.eigvalsh(G)
print(f"  由距离矩阵重建的内蕴 Gram 的本征值 = {np.round(ev, 12)}")
print("  ⟹ 全部非负（半正定）且秩 4 ⟹ 距离矩阵**可精确实现于 ℝ⁴** ⟹ 该单纯形是平的")
ck("J3 单单纯形：δ ≠ 0 但内蕴曲率 = 0（**有缺陷角而无曲率**）",
   delta > 1e-6 and bool((ev > -1e-9).all()),
   f"δ = {math.degrees(delta):.4f}°；Gram 最小本征值 = {ev.min():.2e} ⟹ 可嵌 ℝ⁴ ⟹ 平坦")
note("结论（角度③）", "缺陷角**不是曲率本身**；它是'绕铰链的角剩余'",
     "只在**闭合复形**中才等于曲率（连续极限下 δ_h ≈ ∫_h R）")

# ============================================================
print()
print("=" * 78)
print("角度 ④：结构类型")
print("=" * 78)
print("  界间耦合：连续度规上**逐点**的张量分量（依赖 R_a, R_a', R_a''；是 Ricci 的一个加项）")
print("  缺陷角  ：组合复形上**集中在铰链**的整体角剩余（依赖胶合数与该处的二面角）")
print("  ⟹ 一个是 pointwise 张量分量，一个是 concentrated 角剩余 —— **数学类型不同**")

# ============================================================
print()
print("=" * 78)
print("判定")
print("=" * 78)
print("  · 角度①：量纲不同（1/长度² vs 无量纲）⟹ **不可能是同一个量**")
print("  · 角度②：单因子时耦合 = 0 而曲率 ≠ 0 ⟹ 耦合非曲率的必要条件")
print("  · 角度③：单单纯形时 δ ≠ 0 而曲率 = 0 ⟹ 缺陷角非曲率的充分条件")
print("  · 角度④：逐点张量分量 vs 集中角剩余 ⟹ 类型不同")
print()
print("  ⟹ **判定：不是同一机制。**")
print("    它们的正确关系：**缺陷角是'曲率在铰链上的累积'**（δ_h ≈ ∫_h R，只在闭合复形中成立），")
print("    而界间耦合是'连续度规的一个 Ricci 交叉项'；两者之间的对应**属于环 6 那条未建的映射**。")
print("     另：台账侧不对称 —— 界间耦合**无来历**（D82）；缺陷角**有**（非正则性已导出 + #30/#31 词层实现）")

n_fail = sum(1 for x in CK if not x["pass"])
print()
print("=" * 78)
print(f"合计 {len(CK)} 条：PASS {len(CK)-n_fail} / FAIL {n_fail}")
print("=" * 78)
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "boundary_vs_deficit_verdict.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": CK, "notes": OUT,
               "verdict": {
                   "same_mechanism": False,
                   "angle1_dimension": "交叉项 1/长度²（ρ→sρ 下 ∝1/s²）；缺陷角无量纲",
                   "angle2_single_factor": "单因子 ⟹ 耦合=0 而曲率≠0（耦合非必要条件）",
                   "angle3_single_simplex": "单单纯形 ⟹ δ=208.955°≠0 而曲率=0（缺陷角非充分条件）",
                   "angle4_type": "逐点张量分量 vs 集中角剩余",
                   "correct_relation": "缺陷角 = 曲率在铰链上的累积（δ_h ≈ ∫_h R，需闭合复形）；与界间耦合的对应属环 6 未建映射"}},
              f, ensure_ascii=False, indent=1)
print("结果写入 results/boundary_vs_deficit_verdict.json")
sys.exit(1 if n_fail else 0)