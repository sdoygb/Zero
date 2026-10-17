#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""phase_circle.py —— 粒子槽的相位圆图（24 格 = Z8 × Z3）。

画：
  · 单位圆与 24 条格线（每格 15°）；
  · Z8 八分圆顶点（0°,45°,…）；
  · Z3 三代间隔（120° 的虚线）；
  · 带电轻子三相位（90°,210°,330°）与零点（135°，=9 格）；
  · K4 三槽的模角（ln(3/2), ln2, ln(8/3)）——注意：槽是模能，不是相位，
    这里只把它们按数值当角度画在圆上作参照。

输出：ZeroGPU/results/phase_circle.png
用法：python3 ZeroGPU/src/phase_circle.py
"""
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as fm

for cand in ["PingFang SC", "Arial Unicode MS", "Heiti SC", "STHeiti", "Songti SC"]:
    try:
        fm.findfont(cand, fallback_to_default=False)
        plt.rcParams["font.sans-serif"] = [cand]
        break
    except Exception:
        pass
plt.rcParams["axes.unicode_minus"] = False

fig, ax = plt.subplots(figsize=(9, 9))
ax.set_aspect("equal")
th = np.linspace(0, 2 * np.pi, 400)
ax.plot(np.cos(th), np.sin(th), "k-", lw=1.2)
for k in range(24):
    a = k * np.pi / 12
    ax.plot([0, np.cos(a)], [0, np.sin(a)], color="0.85", lw=0.8, zorder=0)
for k in range(4):
    a = k * np.pi / 2
    ax.plot([0, np.cos(a)], [0, np.sin(a)], color="0.6", lw=1.2, zorder=0)

oct_ang = np.arange(0, 2 * np.pi, np.pi / 4)
ax.scatter(1.06 * np.cos(oct_ang), 1.06 * np.sin(oct_ang), s=70, marker="o",
           facecolors="none", edgecolors="tab:orange", lw=1.6, zorder=5, label="Z8 八分圆顶点")
ax.plot([], [], color="tab:green", ls=":", label="Z3 三代间隔 120°")
for k in range(3):
    a = math.radians(90) + k * 2 * np.pi / 3
    ax.plot([0, np.cos(a)], [0, np.sin(a)], color="tab:green", lw=0.9, ls=":", zorder=4)

for name, d in {"θe=90°": 90, "θμ=210°": 210, "θτ=330°": 330}.items():
    a = math.radians(d)
    ax.scatter([math.cos(a)], [math.sin(a)], s=160, color="tab:red", zorder=6)
    ax.annotate(name, (1.14 * math.cos(a), 1.14 * math.sin(a)), fontsize=11,
                color="tab:red", ha="center", va="center")
a = math.radians(135)
ax.scatter([math.cos(a)], [math.sin(a)], s=340, marker="*", color="purple", zorder=7)
ax.annotate("轻子零点 135°（9 格）", (1.26 * math.cos(a), 1.26 * math.sin(a)),
            fontsize=11, color="purple", ha="center", va="center")

for name, v in [("ln(3/2)", math.log(1.5)), ("ln2", math.log(2)), ("ln(8/3)", math.log(8 / 3))]:
    ax.scatter([math.cos(v)], [math.sin(v)], s=90, marker="^", color="tab:blue", zorder=6)
    ax.annotate(f"K4槽 {name}", (1.30 * math.cos(v), 1.30 * math.sin(v)),
                fontsize=10, color="tab:blue", ha="center")

for d in range(0, 360, 30):
    a = math.radians(d)
    ax.text(0.92 * math.cos(a), 0.92 * math.sin(a), str(d) + "°", fontsize=7,
            color="0.45", ha="center", va="center")

ax.set_title("粒子槽的相位圆：24 格（每格 15°）＝ Z8 八分圆 × Z3 三代", fontsize=14)
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.axis("off")
ax.legend(loc="lower left", fontsize=10, frameon=False)
plt.tight_layout()
out = "ZeroGPU/results/phase_circle.png"
plt.savefig(out, dpi=160, bbox_inches="tight")
print("saved", out)
