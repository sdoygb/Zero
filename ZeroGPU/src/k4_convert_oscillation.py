#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_convert_oscillation.py —— 测相互转化振荡：物质↔时空双向驱动量化
====================================================================================
用户洞察: 引力波 = 物质能量↔时空的相互转化动力学(不是横波传播)
量化模流中物质扇区(粒子源 m) ↔ 时空扇区(势 φ) 的互相驱动:
  同代相关 corr(m_t, φ_t): 物质直接决定时空
  滞后相关 corr(m_t, φ_{t−lag}): 时空→物质 (前代势引导当代物质)
           corr(φ_t, m_{t−lag}): 物质→时空 (前代物质弯曲当代时空)
  相互转化 = 两个方向滞后相关都显著; 单向只有物质→时空
对比: 单向(feedback=0) vs 双向(feedback=0.5)
"""
from __future__ import annotations
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from k4_inhomogeneity_scan import InhomogEngine, laplacian, solve_potential

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
K4 = [[1,2,3],[0,2,3],[0,1,3],[0,1,2]]
L4 = laplacian(K4); n = 4

def analyze(feedback, N=2500, p=1.0):
    eng = InhomogEngine(K4, seed=[8.0,4.0,1.5,1.0], T=4, retention="top2",
                        p=p, beta=0.2, feedback=feedback)
    M = []; P = []
    for _ in range(N):
        eng.run_generation()
        m = eng.mass_dist(eng.record)
        xs = m - 1.0/n
        P.append(solve_potential(L4, xs))
        M.append(m)
    M = np.array(M); P = np.array(P)     # (N, n)
    # 去趋势(消除慢漂移), 聚焦振荡
    Md = M - M.mean(axis=0); Pd = P - P.mean(axis=0)
    for v in range(n):
        Md[:,v] -= np.polyval(np.polyfit(np.arange(N), Md[:,v], 1), np.arange(N))
        Pd[:,v] -= np.polyval(np.polyfit(np.arange(N), Pd[:,v], 1), np.arange(N))
    def corr(a, b):
        a0 = a.ravel() - a.ravel().mean(); b0 = b.ravel() - b.ravel().mean()
        den = np.linalg.norm(a0)*np.linalg.norm(b0)
        return float(a0 @ b0 / den) if den > 0 else 0.0
    res = {}
    for lag in [0, 1, 2, 4, 8]:
        if lag == 0:
            res[0] = corr(Md, Pd)                      # 同代: 物质↔时空同步
        else:
            # corr(m_t, φ_{t−lag}): 时空→物质(前代势引导当代物质)
            m_to_phi = corr(Pd[lag:], Md[:-lag])
            # corr(φ_t, m_{t−lag}): 物质→时空(前代物质弯曲当代时空)
            phi_to_m = corr(Md[lag:], Pd[:-lag])
            res[lag] = {"时空→物质": round(m_to_phi, 4), "物质→时空": round(phi_to_m, 4)}
    return M, P, res

print("="*86)
print("测相互转化振荡：物质↔时空双向驱动量化")
print("="*86)

results = {}
for fb, label in [(0.0, "单向(feedback=0)"), (0.5, "双向(feedback=0.5)")]:
    M, P, res = analyze(feedback=fb)
    print(f"\n【{label}】")
    print(f"  同代相关 corr(m_t,φ_t) = {res[0]:.4f}  (物质↔时空同步度)")
    for lag in [1,2,4,8]:
        r = res[lag]
        print(f"  lag={lag}: 时空→物质={r['时空→物质']:.4f}  物质→时空={r['物质→时空']:.4f}")
    results[label] = res

print("\n" + "="*86)
print("相互转化判定")
print("="*86)
for label in ["单向(feedback=0)", "双向(feedback=0.5)"]:
    r = results[label]
    stm = max(abs(r[lag]["时空→物质"]) for lag in [1,2,4,8])
    mts = max(abs(r[lag]["物质→时空"]) for lag in [1,2,4,8])
    both = stm > 0.05 and mts > 0.05
    verdict = "双向相互转化 ✓ (物质↔时空互驱)" if both else \
              ("仅单向 (物质→时空)" if mts > 0.05 else "无显著相互驱动")
    print(f"  {label}: 时空→物质 max|r|={stm:.4f}  物质→时空 max|r|={mts:.4f} → {verdict}")

os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_convert_oscillation.json"), "w", encoding="utf-8") as f:
    json.dump({"config": {"graph": "K4", "N": 2500, "p": 1.0, "beta": 0.2},
               "results": {k: {str(l): v for l,v in results[k].items()} for k in results}},
              f, ensure_ascii=False, indent=1)
print(f"\n已存 {os.path.join(RES,'k4_convert_oscillation.json')}")
