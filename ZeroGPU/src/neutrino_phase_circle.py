#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""neutrino_phase_circle.py —— 中微子的相位圆图：左右手位置、零点、偏心相位。

画在 24 格（=Z8×Z3）的相位圆上：
  · 轻子零点 135°（灰星）与三代相位 90°/210°/330°（红点）；
  · 中微子零点 144.74°（紫星，A_ν=√(3/2) 的零点）；
  · ν_L 在零点（蓝实心点）、ν_R 在对径 324.74°（蓝空心圈）——这是模型选择；
  · 偏心 2/9 rad = 12.73°（橙弧）；
  · δ0 = π/12 − 2/9 = 2.27°（绿弧）。

注意：ν_L/ν_R 的位置与"对径"是本文的建模选择，不是框架结论。

输出：ZeroGPU/results/neutrino_phase_circle.png
用法：python3 ZeroGPU/src/neutrino_phase_circle.py
"""
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as fm

for cand in ["PingFang SC", "Arial Unicode MS", "Heiti SC", "STHeiti"]:
    try:
        fm.findfont(cand, fallback_to_default=False)
        plt.rcParams["font.sans-serif"] = [cand]
        break
    except Exception:
        pass
plt.rcParams["axes.unicode_minus"] = False

fig, ax = plt.subplots(figsize=(9.5, 9.5))
ax.set_aspect("equal")
th = np.linspace(0, 2 * np.pi, 500)
ax.plot(np.cos(th), np.sin(th), "k-", lw=1.2)
for k in range(24):
    a = k * np.pi / 12
    ax.plot([0, np.cos(a)], [0, np.sin(a)], color="0.88", lw=0.7, zorder=0)
for k in range(4):
    a = k * np.pi / 2
    ax.plot([0, np.cos(a)], [0, np.sin(a)], color="0.65", lw=1.0, zorder=0)

a = math.radians(135)
ax.scatter([math.cos(a)], [math.sin(a)], s=300, marker="*", color="0.55", zorder=6)
ax.annotate("轻子零点 135°", (1.20 * math.cos(a), 1.20 * math.sin(a)), fontsize=10, color="0.4", ha="center")
for name, d in {"θe=90°": 90, "θμ=210°": 210, "θτ=330°": 330}.items():
    a = math.radians(d)
    ax.scatter([math.cos(a)], [math.sin(a)], s=60, color="tab:red", zorder=5)
    ax.annotate(name, (1.12 * math.cos(a), 1.12 * math.sin(a)), fontsize=9, color="tab:red", ha="center")

th0 = 180 - math.degrees(math.acos(1 / math.sqrt(1.5)))
a0 = math.radians(th0)
ax.scatter([math.cos(a0)], [math.sin(a0)], s=420, marker="*", color="purple", zorder=8)
ax.annotate("中微子零点 144.74°", (1.30 * math.cos(a0), 1.30 * math.sin(a0)), fontsize=11, color="purple", ha="center")
ax.scatter([math.cos(a0)], [math.sin(a0)], s=170, color="tab:blue", zorder=9)
ax.annotate("ν_L（在零点）", (1.14 * math.cos(a0 + 0.42), 1.14 * math.sin(a0 + 0.42)), fontsize=11, color="tab:blue", ha="center")
aR = a0 + math.pi
ax.scatter([math.cos(aR)], [math.sin(aR)], s=170, facecolors="none", edgecolors="tab:blue", lw=2.0, zorder=9)
ax.annotate(f"ν_R（对径 {math.degrees(aR) % 360:.2f}°）", (1.16 * math.cos(aR), 1.16 * math.sin(aR)), fontsize=11, color="tab:blue", ha="center")
ax.plot([math.cos(a0), math.cos(aR)], [math.sin(a0), math.sin(aR)], color="tab:blue", lw=1.0, ls="--", alpha=0.6, zorder=4)

oc = 2 / 9
ts = np.linspace(a0, a0 + oc, 60)
ax.plot(1.02 * np.cos(ts), 1.02 * np.sin(ts), color="tab:orange", lw=3, zorder=7)
ax.annotate(f"偏心 2/9 rad = {math.degrees(oc):.2f}°", (1.24 * math.cos(a0 + oc / 2), 1.24 * math.sin(a0 + oc / 2)),
            fontsize=11, color="tab:orange", ha="center")
d0 = math.pi / 12 - 2 / 9
ts = np.linspace(a0, a0 - d0, 40)
ax.plot(0.97 * np.cos(ts), 0.97 * np.sin(ts), color="tab:green", lw=2, zorder=7)
ax.annotate(f"δ₀ = π/12−2/9 = {math.degrees(d0):.2f}°", (0.78 * math.cos(a0 - 0.10), 0.78 * math.sin(a0 - 0.10)),
            fontsize=10, color="tab:green", ha="center")

for d in range(0, 360, 30):
    a = math.radians(d)
    ax.text(0.90 * math.cos(a), 0.90 * math.sin(a), str(d) + "°", fontsize=7, color="0.5", ha="center", va="center")

ax.set_title("中微子：左右手位置、零点、偏心相位（24 格=Z8×Z3）", fontsize=14)
ax.set_xlim(-1.55, 1.55)
ax.set_ylim(-1.55, 1.55)
ax.axis("off")
plt.tight_layout()
out = "ZeroGPU/results/neutrino_phase_circle.png"
plt.savefig(out, dpi=160, bbox_inches="tight")
print("saved", out)
