#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""neutrino_amplitude_winding.py —— 振幅拟合的约定核对（更正版，推导稿 §15）。

背景：相位律是 sqrt(m) ∝ M（Z121），即 m ∝ M^2。于是
    R = m2/m3 = (M2/M3)^2,   X = R^2 = (M2/M3)^4.
早期版本误用 X = (M2/M3)^2（等于把质量律当成 m ∝ M），因此"卷绕振幅
A = 24 ln2/(2π)" 的 0.045% 好结果作废。

本脚本显式算两种约定：
  · 约定 C（正确，m ∝ M^2）：R=(M2/M3)^2, X=R^2；
  · 约定 W（错误，m ∝ M）  ：R=M2/M3,      X=R^2。
并在约定 C 下解出满足 R_obs 的振幅 A_correct，再检验卷绕振幅是否也满足。

输出：ZeroGPU/results/neutrino_amplitude_winding.json
用法：python3 ZeroGPU/src/neutrino_amplitude_winding.py
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


X_obs = 7.42e-5 / 2.510e-3
R_obs = math.sqrt(X_obs)


def masses(A: float):
    th0 = math.acos(-1 / A)
    M = [abs(1 + A * math.cos(th0 + 2 * math.pi * k / 3)) for k in range(3)]
    return th0, M


def r_phase(A: float) -> float:
    """约定 C（相位律 m ∝ M^2）：R = (M2/M3)^2。"""
    _, M = masses(A)
    return (M[1] / M[2]) ** 2


def r_linear(A: float) -> float:
    """约定 W（错误，m ∝ M）：R = M2/M3。"""
    _, M = masses(A)
    return M[1] / M[2]


# 约定 C 下解出 A：R(A) = R_obs（取 A>2 的分支）
lo, hi = 3.0, 10.0
for _ in range(300):
    mid = (lo + hi) / 2
    if r_phase(mid) < R_obs:
        lo = mid
    else:
        hi = mid
A_correct = (lo + hi) / 2
th0, M = masses(A_correct)
R_c = r_phase(A_correct)
X_c = R_c * R_c
print(f"  A_correct = {A_correct:.6f}   A^2 = {A_correct*A_correct:.6f}")
print(f"  零点 θ0 = {math.degrees(th0):.4f}°   M = {[round(x,6) for x in M]}")
print(f"  R = {R_c:.6f}   X = {X_c:.6f}   (观测 R={R_obs:.6f}, X={X_obs:.6f})")

A_wind = 24 * math.log(2) / (2 * math.pi)
thw, Mw = masses(A_wind)
R_w_c = r_phase(A_wind)
print(f"  卷绕振幅 A_wind = 24 ln2/2π = {A_wind:.6f}")
print(f"    约定 C 下 R = {R_w_c:.6f}（应为 {R_obs:.6f}，差 {100*(R_w_c/R_obs-1):+.2f}%）")
print(f"    约定 W 下 R = {r_linear(A_wind):.6f}（这才是它 0.02% 拟合的原因）")

R_oct = 3 - 2 * math.sqrt(2)
print(f"  八分圆 R = 3−2√2 = {R_oct:.6f}（偏差 {100*(R_oct/R_obs-1):+.3f}%）")

check("约定 C 下解出的 A 满足 R_obs（<1e-6）", abs(R_c / R_obs - 1) < 1e-6, f"{100*(R_c/R_obs-1):+.2e}%")
check("A_correct 与 3√2 差 >1%", abs(A_correct / (3 * math.sqrt(2)) - 1) > 0.01,
      f"{100*(A_correct/(3*math.sqrt(2))-1):+.2f}%")
check("卷绕振幅在约定 C 下失败（R 偏离 >50%）", abs(R_w_c / R_obs - 1) > 0.5,
      f"{100*(R_w_c/R_obs-1):+.1f}%")
check("卷绕振幅在约定 W 下成立（R 偏离 <0.2%）", abs(r_linear(A_wind) / R_obs - 1) < 0.002,
      f"{100*(r_linear(A_wind)/R_obs-1):+.3f}%")
check("八分圆在约定 C 下成立（<0.5%）", abs(R_oct / R_obs - 1) < 0.005, f"{100*(R_oct/R_obs-1):+.3f}%")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {
    "pass": PASS, "fail": FAIL,
    "A_correct_phase_law": A_correct, "A_wind_24ln2_2pi": A_wind,
    "R_obs": R_obs, "X_obs": X_obs,
    "R_phase_law_fit": R_c, "X_phase_law_fit": X_c,
    "R_wind_under_phase_law": R_w_c, "R_wind_under_linear": r_linear(A_wind),
    "R_octagon": R_oct,
}
with open(os.path.join(OUT, "neutrino_amplitude_winding.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/neutrino_amplitude_winding.json")
sys.exit(0 if not FAIL else 1)
