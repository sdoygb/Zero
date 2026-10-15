#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_confinement_closure.py —— Zero 版禁闭: 色闭合 = Berry 2π 判据, 单色荷不闭合
======================================================================================
武器(最老理论 0.4/7.9): 闭合 ⟺ Berry 相位 ∮A = 2πn; 单色荷相位 2π/3 不闭合; 三色荷=2π 闭合.
Zero 底座: 闭合 2π 分 8 (seed 8), 三代间隔 = 1/3 圈 = 2π/3.
核心: "三" 同时给三代间隔与色三重态闭合 → 质子 = 三色闭合 = Berry 2π.
"""
from __future__ import annotations
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

P_("=" * 96)
P_("Zero 版禁闭 · k4_confinement_closure.py  (色闭合 = Berry 2π 判据)")
P_("=" * 96)

tau = 2*np.pi
third = tau/3            # 1/3 圈 = 三代间隔 / 色闭合间隔
eighth = tau/8           # 闭合 2π 分 8 (seed 8) 单格

P_("闭合 2π 分 8 (seed 8): 单格 = π/4 = %.4f° = 45°" % np.degrees(eighth))
P_("三代间隔 = 1/3 圈 = 2π/3 = %.4f° = 120°" % np.degrees(third))
P_("")

# ============ 1. Berry 相位闭合判据 ============
def closed(bp, tol=1e-9):
    return abs(bp/tau - round(bp/tau)) < tol

P_("【1. 闭合判据 ∮A = 2πn】")
rows = [
    ("单色荷 (夸克)", third, 1),
    ("二色荷", 2*third, 2),
    ("三色荷 (色单态/质子)", 3*third, 3),
    ("四色荷", 4*third, 4),
    ("δ⁸ 最大圆满回路 (对照)", tau, 1),
    ("单步 δ (0.4 对照)", tau/8, 1),
]
P_("  %-22s %10s %12s %8s" % ("结构", "Berry 相位", "rad/2π", "闭合?"))
for nm, bp, k in rows:
    ratio = bp/tau
    c = "✓ 闭合" if closed(bp) else "✗ 不闭合"
    P_("  %-22s %8.4f° %10.4f %10s" % (nm, np.degrees(bp), ratio, c))

P_("")
P_("  结论: 单色荷 2π/3 不闭合 → 禁闭; 三色荷 2π 闭合 → 色单态(质子).")
P_("        四色荷 = 4×2π/3 = 8π/3, 与 2π 差 2π/3, 也不闭合(4/3 圈).")

# ============ 2. 三色相位叠加 ============
P_("\n" + "=" * 96)
P_("【2. 三色荷相位叠加: α+β+γ=0 (色单态)】")
ph = np.array([0, third, 2*third])     # 0, 2π/3, 4π/3
S = np.exp(1j*ph).sum()
P_("  三色相位: 0, 2π/3, 4π/3 rad")
P_("  叠加 Σe^{iφ} = %.4f %+.4fi  |Σ| = %.2e  → %s"
  % (S.real, S.imag, abs(S), "色单态闭合 (和为0)" if abs(S)<1e-12 else "非单态"))

# 单色荷叠加
S1 = np.exp(1j*np.array([0]))
P_("  单色荷 |Σ|=%.4f ≠ 0 → 净色荷 ≠ 0 → 不闭合" % float(np.abs(S1[0])))

# ============ 3. "三"的双重角色 ============
P_("\n" + "=" * 96)
P_("【3. '三' 的双重角色: 三代间隔 = 色闭合间隔】")
P_("  三代 (Z121): 相角差 2π/3 (seed 8 的三项) → 三代质量槽")
P_("  色闭合 (7.9): 三色荷相位差 2π/3 → 色单态")
P_("  同一个 1/3 圈: 三代的'三' 与 色荷的'三' 同源 (Λ=3 / seed 三项)")
P_("  单色荷 = 2π/3 (不圆满) ⇔ 单代 = 部分结构: 都需整体闭合")

# ============ 4. 质子 = 三色闭合: 质量机制待接 ============
P_("\n" + "=" * 96)
P_("【4. 质子单粒子路径 (下一步)】")
P_("  质子 = K₄ 上三色闭合结构 (Berry 2π); 直接套质量=偏离零点² (Z121)")
P_("  m_p = 938.27 MeV 需由某稳态/角无参数给出 — 待攻")

out = {
  "single_quark_berry": third, "three_quark_berry": 3*third,
  "single_not_closed": not closed(third), "triplet_closed": closed(3*third),
  "color_singlet_sum": abs(S), "third_of_circle_deg": np.degrees(third),
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_confinement_closure.json"), "w", encoding="utf-8") as f_:
    json.dump(out, f_, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "k4_confinement_closure.txt"), "w", encoding="utf-8") as f_:
    f_.write("\n".join(LOG))
P_("\n已存 k4_confinement_closure.json / .txt")
