#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_selfenergy.py —— 时空源升级：物质+时空自身能量共同拉伸时空（GR 自参考）
====================================================================================
用户洞察: 物质=能量, 物质和能量都能引起时空拉伸。
GR 核心: 时空曲率由总能量-动量决定(不只质量), 且引力场自身有能量(自参考非线性)。
升级: 时空源 x = (1−s)·物质能量 + s·时空自身能量(曲率能量), 自洽迭代。
  物质能量 = 粒子频率槽质量分布 m
  时空自身能量 = 曲率能量 κ² (κ=−Lφ, Z033 曲率代理), 归一化
  s = 自参考强度(时空自身能量参与弯曲时空的比例)
自洽: φ 由源 x 决定, 源含 φ 的能量 → 迭代求解(GR 自洽场精神)
测: s 扫描 → 时空拉伸(φ范围/曲率) 是否被自参考放大(非线性)
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

def solve_selfconsistent(m, s, L, n_iter=30, tol=1e-9):
    """自洽求解: x=(1−s)·m + s·(κ²归一化), φ=solve(L,x), 迭代到自洽。"""
    m = m / m.sum()
    phi = solve_potential(L, m - 1.0/n)
    for it in range(n_iter):
        kappa = -L @ phi                       # 曲率代理 κ=−Lφ
        sel = kappa**2
        sel = sel / sel.sum() if sel.sum() > 0 else np.ones(n)/n
        x = (1-s) * m + s * sel
        x -= x.mean()
        phi_new = solve_potential(L, x)
        if np.linalg.norm(phi_new - phi) < tol:
            phi = phi_new; break
        phi = phi_new
    return phi

def engine_phi(fb, p=1.0, N=1500, s=None):
    """跑耦合引擎, 收集势(可选自参考源)。"""
    eng = InhomogEngine(K4, seed=[8.0,4.0,1.5,1.0], T=4, retention="top2",
                        p=p, beta=0.2, feedback=fb)
    for _ in range(100): eng.run_generation()   # 预热
    tail = []
    for _ in range(N):
        eng.run_generation()
        m = eng.mass_dist(eng.record)
        if s is None:
            phi = solve_potential(L4, m - 1.0/n)
        else:
            phi = solve_selfconsistent(m, s, L4)
        tail.append(phi)
    return np.array(tail)

print("="*90)
print("时空源升级：物质+时空自身能量共同拉伸时空（GR 自参考）")
print("="*90)

out = {"config": {"graph": "K4", "selfenergy": "kappa^2 (kappa=-Lphi)", "coup": "x=(1-s)*m+s*kappa^2_norm"}}
# s 扫描: 自参考强度 → 时空拉伸(φ范围/曲率)
print("\n【1】自参考 s 扫描 → 时空拉伸（单向引擎 fb=0）")
out["s_scan_fb0"] = {}
for s in [0.0, 0.3, 0.6, 0.9]:
    phis = engine_phi(fb=0.0, s=s, N=500)
    phi_range = float(np.mean(np.max(phis,axis=1) - np.min(phis,axis=1)))
    phi_norm = float(np.mean(np.linalg.norm(phis,axis=1)))
    tail_std = float(np.std(phis[-200:], axis=0).max())
    out["s_scan_fb0"][s] = {"phi_range": round(phi_range,5), "phi_norm": round(phi_norm,5),
                            "tail_std": round(tail_std,6)}
    print(f"  s={s:.1f}: 时空拉伸(φ范围)={phi_range:.5f}  |φ|={phi_norm:.5f}  "
          f"末段波动={tail_std:.6f}")

# 自参考下的相互转化(fb): 时空动态
print("\n【2】自参考 + 双向耦合 → 时空动态")
out["coupling"] = {}
for s, fb in [(0.0, 0.8), (0.6, 0.8), (0.6, 0.4), (0.3, 0.8)]:
    phis = engine_phi(fb=fb, s=s, N=1500)
    tail_std = float(np.std(phis[-500:], axis=0).max())
    phi_range = float(np.mean(np.max(phis[-500:],axis=1) - np.min(phis[-500:],axis=1)))
    out["coupling"][f"s{s}_fb{fb}"] = {"phi_dyn": round(tail_std,6), "phi_range": round(phi_range,5)}
    print(f"  s={s:.1f} fb={fb:.1f}: 时空动态(末段波动)={tail_std:.6f}  拉伸(φ范围)={phi_range:.5f}")

os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_selfenergy.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("\n已存 k4_selfenergy.json")
