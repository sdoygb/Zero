#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_angle_from_closure.py —— 从 Zero 闭合路径(2π)逼出质量机制的"角"
======================================================================================
用户提示: 闭合的路径在 Zero 是 2π. 从这条路把"角"逼出来.

Zero 已知: 稳态不对称 2:1 (三项, 分母3); seed=[8,4,1.5,1] (首项8).
质量机制(参考 cosmos, 底座要 Zero): M(θ)=1+√2cosθ, 三代 θ_i=δ+2πi/3,
      振幅 √2 的角 = π/4, 零点 = 3π/4, 三代间隔 = 2π/3.

本脚本验证: 这些角能否从 闭合2π 的"分数份"逼出 (分数来自 Zero 的 3 和 8).
"""
from __future__ import annotations
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

tau = 2 * np.pi  # 闭合路径 = 2π (一圈)

# Zero 数: 3 (三项/稳态分母), 8 (seed 首项), 4 (seed 次项), 2 (相邻比)
P_("=" * 96)
P_("从 Zero 闭合路径(2π)逼出质量机制的角 · k4_angle_from_closure.py")
P_("=" * 96)
P_(f"闭合路径 τ=2π = {np.degrees(tau):.1f}°  (Zero 一圈)")
P_(f"Zero 数: 稳态不对称三项(3), seed 首项(8), seed 次项(4), 相邻比(2)")
P_("")

out = {}
def angle(name, frac, note=""):
    a = frac * tau
    P_(f"  {name:22s} = {frac}·2π = {np.degrees(a):10.4f}°   {note}")
    return a

# ============ 1. 质量机制需要的角 (目标) ============
P_("【质量机制需要的角 (参考)】")
P_(f"  振幅 √2 的角 π/4     = 45.000°    (2cos(π/4)=√2={np.sqrt(2):.8f})")
P_(f"  零点 3π/4            = 135.000°   (M=1+√2cos(3π/4)=0)")
P_(f"  三代间隔 2π/3        = 120.000°   (Z₃ 轨道)")
P_(f"  电子相角 δ(反解)      = 132.73°    (靠近零点, 差 δ₀=2.27°)")
P_("")

# ============ 2. 从 Zero 分数逼角 ============
P_("【Zero 闭合 2π 分份 → 角 (分数来自 3 与 8)】")
# seed 8 → 八分; 稳态 3 → 三分之一; 组合
fracs = {
    "√2 角 π/4   = (1/8)·2π  [seed 8]": 1/8,
    "零点 3π/4   = (3/8)·2π  [3·1/8: 三项×八分]": 3/8,
    "三代 2π/3   = (1/3)·2π  [三项]": 1/3,
    "π/2         = (2/8)·2π  [2 相邻比]": 2/8,
    "π           = (4/8)·2π  [4 seed 次项]": 4/8,
}
results = {}
for name, f in fracs.items():
    a = angle(name, f)
    results[name] = {"frac": f, "deg": round(float(np.degrees(a)), 4)}
P_("")

# ============ 3. 校验: 这些角能否重建质量机制 ============
P_("【校验: 用 Zero 角重建质量算子 M(θ)=1+√2cosθ】")
P_("  假设 (Zero 内生):")
P_("    √2 = 2cos(π/4), π/4 = (1/8)·2π  ← seed 首项 8 (闭合分 8 份)")
P_("    零点 = 3π/4 = (3/8)·2π           ← 三项(3) × 八分(1/8)")
P_("    三代间隔 = 2π/3 = (1/3)·2π       ← 三项(3)")
P_("  则 M(θ_i)=1+2cos(π/4)·cos(θ_i), θ_i=δ+2πi/3, δ 扫向零点附近 → 电子轻")
P_("")

# 校验 √2 恒等式
c8 = np.cos(2*np.pi/8)
P_(f"  ✓ 2cos(2π/8) = {2*c8:.10f}  (=√2={np.sqrt(2):.10f}?)")
out["sqrt2_check"] = float(2*c8)
# 零点
P_(f"  ✓ M(3π/4) = 1+√2cos(3π/4) = {1+np.sqrt(2)*np.cos(3*np.pi/4):.2e} (≈0)")
out["zero_check"] = float(1+np.sqrt(2)*np.cos(3*np.pi/4))
# 角差: 零点 − √2角 = π/2 = θ_物质
diff = 3*np.pi/4 - np.pi/4
P_(f"  ✓ 零点(3π/4) − √2角(π/4) = π/2 = θ_物质  ({np.degrees(diff):.1f}°  = 2·(1/8)·2π 或 2/8 圈)")
out["theta_matter"] = float(diff)

# ============ 4. 用 Zero 角重建三代质量, δ 反解 ============
from scipy.optimize import least_squares
obs = np.array([0.511, 105.66, 1776.86])
def Mth(d): return np.array([1+np.sqrt(2)*np.cos(d+2*np.pi*i/3) for i in range(3)])
def resid(d):
    v = Mth(d[0]); v = v/v.min(); o = np.sqrt(obs)/np.sqrt(obs).min()
    return (v-o)/o
sol = least_squares(resid, [2.3], bounds=([-np.pi],[np.pi]))
m = Mth(sol.x[0])**2
koide = float(m.sum()/np.sqrt(m).sum()**2)
P_("\n【δ 反解 + 质量重建 (形状用 Zero 角, δ 自由)】")
P_(f"  δ = {np.degrees(sol.x[0])%360:.4f}°  (反解)")
P_(f"  Koide = {koide:.5f}  (观测 0.6667)")
P_(f"  质量比 = {m/m.min()}")
P_(f"  δ 与零点 3π/4 差 δ₀ = {np.degrees(3*np.pi/4-sol.x[0])%360:.4f}° = {float(3*np.pi/4-sol.x[0])%tau:.6f} rad")
out["delta_recon"] = {"delta_deg": float(np.degrees(sol.x[0])%360),
                      "koide": round(koide,5),
                      "ratio": (m/m.min()).round(1).tolist(),
                      "d0_rad": float(3*np.pi/4-sol.x[0])%tau}

os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_angle_from_closure.json"), "w", encoding="utf-8") as f_:
    json.dump(out, f_, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "k4_angle_from_closure.txt"), "w", encoding="utf-8") as f_:
    f_.write("\n".join(LOG))
P_("\n已存 k4_angle_from_closure.json / .txt")
