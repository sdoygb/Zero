#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_binary_gw.py —— 观测对应: 双星轨道衰减 + 引力波 chirp (可证伪预言)
====================================================================================
把 Zero 四极辐射定律(Z110: P∝(d³I/dt³)²) 应用到双星轨道(牛顿力学+辐射能量损失):
  双星: 质量 M=Σm, 约化 μ=m1m2/M, 圆周轨道半径 r, 角速度 Ω=√(GM/r³)
  四极矩:  I = μ·r² ;  d³I/dt³ = 3μ r² Ω³ (圆周轨道三阶导)
  辐射功率: P = (1/5)·(d³I/dt³)²        ← Zero 四极定律
  轨道能量损失:  dE/dt = −P,  E=−GMμ/(2r) → 轨道收缩 dr/dt<0
预言(可证伪):
  1) 轨道周期衰减:  dP_orb/dt < 0 (轨道收缩)
  2) 引力波 chirp:  f_gw(t)=Ω(t)/π 随时间上升 (四极主谐波, 频率爬升)
  3) 双星寿命:      从 r0 到并合的时间 τ
对比 GR(Hulse-Taylor PSR B1913+16): m≈1.4M☉, P_orb≈7.75h, dP/dt≈−2.4×10⁻¹² s/s
"""
from __future__ import annotations
import sys, os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))

def evolve_binary(m1, m2, r0, steps=200000, dt=1e-3):
    """牛顿圆周轨道 + 四极辐射能量损失. 返回轨道演化历史."""
    G = 1.0
    M = m1 + m2; mu = m1*m2/M
    r = r0
    E = -G*M*mu/(2*r)
    hist = []
    merged = False
    for t in range(steps):
        Om = np.sqrt(G*M/r**3)
        # Zero 四极辐射功率: P=(1/5)(d³I/dt³)², d³I/dt³=3μr²Ω³
        d3I = 3*mu*r**2*Om**3
        P = (1.0/5.0)*d3I**2
        # 轨道能量损失: dE/dt=−P → dr/dt
        drdt = -(2*r**2/(G*M*mu))*P
        r = r + drdt*dt
        if r <= 1e-4:
            merged = True
            break
        hist.append((t*dt, r, Om, P, 2*np.pi/Om, Om/np.pi))
        E = E - P*dt
    hist = np.array(hist)
    return hist, M, mu, merged

def main():
    out = {"config": {"model": "binary + Zero quadrupole radiation P=(1/5)(d³I)²",
                      "gw_freq": "f_gw=2·f_orb=Ω/π (quadrupole harmonic)",
                      "GR_ref": "Hulse-Taylor PSR B1913+16: dP/dt≈−2.4e-12 s/s, P_orb≈7.75h"}}
    print("="*90)
    print("观测对应: 双星轨道衰减 + 引力波 chirp (可证伪预言)")
    print("="*90)

    # 用 Hulse-Taylor 类参数(无量纲): m1=m2=1, M=2, mu=0.5
    m1=m2=1.0; M=m1+m2; mu=m1*m2/M
    for label, r0 in [("r0=5", 5.0), ("r0=8", 8.0), ("r0=10", 10.0)]:
        hist, Mv, muv, merged = evolve_binary(m1, m2, r0)
        t = hist[:,0]; r = hist[:,1]; Om = hist[:,2]; P = hist[:,3]
        Porb = hist[:,4]; fgw = hist[:,5]
        # 轨道周期衰减率 dP_orb/dt (末段拟合)
        dPdt = float(np.polyfit(t[-2000:], Porb[-2000:], 1)[0])
        # chirp: f_gw 从初始到末段的变化
        f0 = fgw[0]; f_end = fgw[-1]
        chirp_rate = float((f_end - f0)/(t[-1]-t[0]))
        lifetime = float(t[-1])
        dI3_0 = float(3*muv*r[0]**2*Om[0]**3)
        P0 = float(P[0])
        # 无量纲周期衰减率
        dPdt_norm = float(dPdt / Porb[0])
        out[label] = {"r0": float(r[0]), "M": Mv, "mu": muv, "merged": merged,
                      "P_orb_initial": round(float(Porb[0]),4),
                      "f_gw_initial": round(f0,6), "f_gw_final": round(f_end,6),
                      "chirp_rate_dfdt": round(chirp_rate,6),
                      "dP_orb_dt": round(dPdt,8), "dP_orb_dt_norm": round(dPdt_norm,8),
                      "d3I_initial": round(dI3_0,6), "P_initial": round(P0,8),
                      "lifetime_to_merge": round(lifetime,2)}
        print(f"\n[{label}] m1=m2=1, M={Mv}, mu={muv:.3f}, r0={r[0]:.1f}, 并合={merged}")
        print(f"  初始轨道周期 P_orb = {Porb[0]:.4f}")
        print(f"  引力波频率: 初始 f_gw={f0:.5f} → 末段 f_gw={f_end:.5f}  (chirp 上升={chirp_rate>0})")
        print(f"  chirp 率 df_gw/dt = {chirp_rate:.2e}  (>0 ⇒ 频率爬升)")
        print(f"  轨道周期衰减率 dP_orb/dt = {dPdt:.2e}  (<0 ⇒ 轨道收缩)")
        print(f"  归一化衰减率 dP/dt÷P_orb = {dPdt_norm:.2e}")
        print(f"  初始辐射功率 P0={P0:.2e}, d³I/dt³={dI3_0:.3e}")
        print(f"  双星寿命(到并合) τ = {lifetime:.1f}")

    # 验证 P∝(d³I)² 在轨道演化中: 采样历史点
    print("\n[验证] 轨道演化中 P∝(d³I)² 保持")
    hist, Mv, muv, merged = evolve_binary(m1, m2, 8.0)
    t=hist[:,0]; r=hist[:,1]; Om=hist[:,2]; P=hist[:,3]
    d3I = 3*muv*r**2*Om**3
    sel = np.arange(0, len(t), max(1,len(t)//20))
    ratio = P[sel]/d3I[sel]**2
    ratio_std = float(np.std(ratio)/np.mean(ratio))
    out["P_quad_hold"] = {"ratio_mean": round(float(np.mean(ratio)),6),
                          "ratio_std_rel": round(ratio_std,4)}
    print(f"  P/(d³I)² 在轨道演化中比值 {np.mean(ratio):.4f} (相对离散 {ratio_std:.4f})")

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_binary_gw.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已存 k4_binary_gw.json")

if __name__ == "__main__":
    main()
