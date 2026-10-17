#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""neutrino_prediction_band.py —— 中微子预言的「中心值＋带」报告（Z163 §7.3）。

中心值：Z162 的形状（A²=3/2、m1=0、R）＋ Z163 的标度（M_R=√(M_I M_X)、m_D=m_τ/2）。
带：M_X 的 ±5%／±10%（Z155 一环值的阈值示意）与形状的 A²∈[3/2, 1.513633]（Z164 的残差带）。
逐项给出 m2、m3、Σmν、mβ、mββ、X 的中心值与带，并检查它们是否覆盖数据。

输出：ZeroGPU/results/neutrino_prediction_band.json
用法：python3 ZeroGPU/src/neutrino_prediction_band.py
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


V, M_Z, M_TAU = 246.22, 91.1876, 1.77686
M_I = M_Z * math.exp(4 * math.pi)
M_X0 = 1.009e13
M3_OBS = math.sqrt(2.510e-3)
D21, D31 = 7.42e-5, 2.510e-3
Ue2, Ue3 = 0.5443, 0.1492


def shape(A2=1.5):
    s = (math.sqrt(3) / 2) * math.sqrt(A2 - 1)
    f2, f3 = 1.5 - s, 1.5 + s
    R = (f2 / f3) ** 2
    X = R**2
    return R, X


def obs(MX, A2=1.5):
    R, X = shape(A2)
    m3 = (M_TAU / 2) ** 2 / math.sqrt(M_I * MX) * 1e9
    m2 = R * m3
    Sigma = m2 + m3
    m_beta = math.sqrt(Ue2**2 * m2**2 + Ue3**2 * m3**2) * 1e3
    m_bb = (Ue2**2 * m2 + Ue3**2 * m3) * 1e3
    return {"m3": m3, "m2": m2, "Sigma": Sigma, "m_beta": m_beta, "m_bb": m_bb, "X": X}


A2_DATA = 1.513633
base = obs(M_X0)
print("== 中心值（M_X=1.009e13，A²=3/2） ==")
print(f"   m2={base['m2']*1e3:.2f}  m3={base['m3']*1e3:.2f}  Σ={base['Sigma']*1e3:.2f}  "
      f"mβ={base['m_beta']:.2f}  mββ={base['m_bb']:.2f}  X={base['X']:.6f}（meV）")

bands = {}
for name, f in (("±5%", 0.05), ("±10%", 0.10)):
    lo = obs(M_X0 * (1 - f))
    hi = obs(M_X0 * (1 + f))
    def srt(a, b):
        return (min(a, b), max(a, b))
    bands[name] = {"m3": srt(lo["m3"], hi["m3"]), "Sigma": srt(lo["Sigma"], hi["Sigma"]),
                   "m_bb": srt(lo["m_bb"], hi["m_bb"])}
    print(f"   M_X{name}: m3 {lo['m3']*1e3:.2f}–{hi['m3']*1e3:.2f}；"
          f"Σ {lo['Sigma']*1e3:.2f}–{hi['Sigma']*1e3:.2f}；"
          f"mββ {lo['m_bb']:.2f}–{hi['m_bb']:.2f} meV")

shape_band = obs(M_X0, A2_DATA)
print(f"   A²=1.513633: X={shape_band['X']:.6f}（数据 X={D21/D31:.6f}）")

print("\n== 覆盖检查 ==")
b10 = bands["±10%"]
check("M_X ±10% 的带覆盖 m3 的实测值", b10["m3"][0] <= M3_OBS <= b10["m3"][1],
      f"{b10['m3'][0]*1e3:.2f}–{b10['m3'][1]*1e3:.2f} vs {M3_OBS*1e3:.2f}")
check("Σ 的中心值与带都在 DESI DR2 上界（64 meV）内", b10["Sigma"][1] < 0.064,
      f"≤{b10['Sigma'][1]*1e3:.2f}")
check("mββ 的带远低于下一代灵敏度（10 meV）", b10["m_bb"][1] < 10,
      f"≤{b10['m_bb'][1]:.2f}")
check("形状带（A² 到 1.5136）把 X 移到数据值附近（<0.5%）",
      abs(shape_band["X"] / (D21 / D31) - 1) < 0.005,
      f"{shape_band['X']:.6f} vs {D21/D31:.6f}")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "center": base, "bands": bands,
       "shape_band_X": shape_band["X"]}
with open(os.path.join(OUT, "neutrino_prediction_band.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/neutrino_prediction_band.json")
sys.exit(0 if not FAIL else 1)
