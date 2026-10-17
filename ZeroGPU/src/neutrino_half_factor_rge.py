#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""neutrino_half_factor_rge.py —— 检验「1/2 是否来自 RG 跑动」（Z163 §3 的旁证）。

问题：m_D=m_τ/2 的 1/2 会不会是「m_τ 从低能跑到 seesaw 标度」的结果？
做法：SM 一圈 RGE（g1,g2,g3,y_t,y_τ）从 M_Z 积分到 M_R=√(M_I M_X)，
看 y_τ(M_R)/y_τ(M_Z) 是否 ≈1/2。

若比值 ≈0.8（而不是 0.5），则 1/2 不是跑动，只能是结构（两槽轨道归一化）。

输出：ZeroGPU/results/neutrino_half_factor_rge.json
用法：python3 ZeroGPU/src/neutrino_half_factor_rge.py
"""
from __future__ import annotations

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


M_Z = 91.1876
M_R = 1.624300e10
V = 246.22
M_TAU = 1.77686

# 初始值（M_Z，SM，近似）
g1 = math.sqrt(4 * math.pi * 0.01017)     # U(1)_Y，GUT 归一
g2 = math.sqrt(4 * math.pi * 0.0338)
g3 = math.sqrt(4 * math.pi * 0.1179)
y_t = math.sqrt(2) * 172.7 / V
y_b = math.sqrt(2) * 4.18 / V
y_tau = math.sqrt(2) * M_TAU / V


def dln_y(y, g1, g2, g3, y_t, y_b, y_tau, which):
    """一圈：d ln y / d ln μ。"""
    # 系数（SM，Yukawa 对角近似；上型/下型/轻子）
    if which == "t":
        return (1 / (16 * math.pi**2)) * (4.5 * y_t**2 + 1.5 * y_b**2 + y_tau**2
                                          - 8 * g3**2 - 2.25 * g2**2 - 17 / 12 * g1**2)
    if which == "tau":
        return (1 / (16 * math.pi**2)) * (1.5 * y_tau**2 + 3 * y_b**2
                                          - 2.25 * g2**2 - 15 / 4 * g1**2)
    raise ValueError


def dln_g(g, b, g3=None, g2=None, g1=None):
    """一圈规范耦合：d ln g / d ln μ = b g²/(16π²)。"""
    return b * g**2 / (16 * math.pi**2)


def run(mu0, mu1, n=4000):
    global g1, g2, g3, y_t, y_b, y_tau
    t0, t1 = math.log(mu0), math.log(mu1)
    h = (t1 - t0) / n
    state = [g1, g2, g3, y_t, y_b, y_tau]
    for _ in range(n):
        g1, g2, g3, y_t, y_b, y_tau = state
        # b 系数（SM）：g1: 41/10，g2: −19/6，g3: −7
        k1 = [dln_g(g1, 41 / 10) * g1]
        k2 = [dln_g(g2, -19 / 6) * g2]
        k3 = [dln_g(g3, -7) * g3]
        kt = [dln_y(y_t, g1, g2, g3, y_t, y_b, y_tau, "t") * y_t]
        kb = [dln_y(y_b, g1, g2, g3, y_t, y_b, y_tau, "b") * y_b if False else 0.0]
        # y_b 简化：忽略（量级小）
        ktau = [dln_y(y_tau, g1, g2, g3, y_t, y_b, y_tau, "tau") * y_tau]
        # 简单欧拉（h 小）
        state = [g1 + h * k1[0], g2 + h * k2[0], g3 + h * k3[0],
                 y_t + h * kt[0], y_b, y_tau + h * ktau[0]]
    g1, g2, g3, y_t, y_b, y_tau = state


y_tau0 = y_tau
run(M_Z, M_R)
ratio = y_tau / y_tau0
print(f"   y_τ(M_Z)={y_tau0:.5f}  y_τ(M_R)={y_tau:.5f}  比值 {ratio:.4f}")
print(f"   ⟹ m_D/m_τ = {ratio:.4f}（若 1/2 来自跑动，应 ≈0.5）")
check("跑动比值在 0.75–0.9（远大于 1/2）", 0.75 < ratio < 0.90, f"{ratio:.4f}")
check("比值与 1/2 的差距 >0.25", abs(ratio - 0.5) > 0.25, f"{abs(ratio-0.5):.3f}")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "y_tau_MZ": y_tau0, "y_tau_MR": y_tau,
       "ratio": ratio, "M_Z": M_Z, "M_R": M_R}
with open(os.path.join(OUT, "neutrino_half_factor_rge.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/neutrino_half_factor_rge.json")
sys.exit(0 if not FAIL else 1)
