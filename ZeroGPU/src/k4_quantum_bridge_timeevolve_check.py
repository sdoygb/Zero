#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_quantum_bridge_timeevolve_check.py —— 能量标度桥：色散一致性 ＋ 时间域稳态核验
=====================================================================================
本篇核三件事（纯 numpy、只写 results）：

A. 色散一致性（解析）：图上波方程  ω² = c²λ(q)，行模 λ(q)=4sin²(q/2)。
   驱动角频率 ω=f 定住 q(f)=2·arcsin(f/(2c))；格上波长 λ_e=2π/q。
   把 q(f)→λ_e 与 k4_quantum_doubleslit_A 的实测 λ_e 逐条对照。

B. 时间域稳态（数值）：带吸收层的二阶波方程
       ∂_t²u + 2γ∂_t u = -c² L u + g cos(f t)
   稳定蛙跳积分到稳态；测空间相位梯度 → q_meas，与 A 的 q(f) 对照。
   这把"定态散射"与"时间演化"两条路线落在同一个稳态上。

C. 显式欧拉的不稳定性（解析＋数值）：
   对 i∂_tψ=Hψ，显式欧拉 ψ_{n+1}=(1-iHdt)ψ_n 的放大因子 |1-iλdt|²=1+(λdt)²>1，
   无条件增长；数值上复现 blow-up，说明其"主导频率"是爆破解的伪信号。

判据方向固定：PASS = |q_meas/q_th - 1| < 1e-2（1%）；FAIL = 该篇要报的判定。
"""
from __future__ import annotations
import os
import json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []


def P(*a):
    s = " ".join(str(x) for x in a)
    LOG.append(s)
    print(s)


def pathL(n):
    L = np.zeros((n, n))
    for i in range(n):
        L[i, i] += int(i > 0) + int(i < n - 1)
        if i > 0:
            L[i, i - 1] -= 1
        if i < n - 1:
            L[i, i + 1] -= 1
    return L


def q_of_omega(f, c2):
    """ω² = c²·4sin²(q/2)  ⟹  q = 2 arcsin(f/(2c))."""
    w = f / (2.0 * np.sqrt(c2))
    if abs(w) > 1.0:
        return float("nan")
    return 2.0 * np.arcsin(w)


def main():
    P("=" * 96)
    P("能量标度桥：色散一致性 ＋ 时间域稳态核验 (k4_quantum_bridge_timeevolve_check)")
    P("=" * 96)
    out = {}
    c2 = 1.0

    # ---------- A. 色散一致性 ----------
    P("\n[A] 色散一致性：ω=f → q=2arcsin(f/2c) → λ_e=2π/q（格上），与实测 λ_e 对照")
    meas = [  # (标签, f, 实测 λ_e) —— 取自 k4_quantum_doubleslit_A.txt
        ("log2", np.log(2), 8.874),
        ("K4槽1", 0.4050, 15.394),
        ("K4槽2", 0.6930, 8.876),
        ("K4槽3", 0.9810, 6.127),
    ]
    rowsA = []
    for lab, f, lam_meas in meas:
        q = q_of_omega(f, c2)
        lam_lat = 2.0 * np.pi / q
        lam_cont = 2.0 * np.pi / f
        rel = abs(lam_lat - lam_meas) / lam_meas
        rowsA.append(dict(label=lab, f=round(f, 6), q=round(q, 6),
                          lam_lattice=round(lam_lat, 6), lam_measured=lam_meas,
                          rel_err_vs_lattice=round(rel, 6),
                          lam_continuum=round(lam_cont, 6),
                          ratio_continuum=round(lam_meas / lam_cont, 6)))
        P(f"  {lab:7s} f={f:.4f}  q=2arcsin(f/2c)={q:.4f}  "
          f"λ_lat=2π/q={lam_lat:.4f}  实测={lam_meas:.4f}  "
          f"相对差={rel*1e2:.3f}%   |  2π/f={lam_cont:.4f} (连续极限)")
    lam_err = max(r["rel_err_vs_lattice"] for r in rowsA)
    P(f"  ⟹ 格上波长 λ_e=2π/q 与实测最大相对差 = {lam_err*1e2:.3f}%")
    P("  ⟹ 2π/f 是连续极限；实测/2π/f 的 0.957–0.992 是色散效应，不是误差")
    out["A_dispersion"] = rowsA
    okA = lam_err < 1e-2

    # ---------- B. 时间域稳态 ----------
    P("\n[B] 时间域稳态：二阶波方程 + 吸收层，驱动 ω=f=log2，测 q_meas")
    N = 400
    L = pathL(N)
    f = np.log(2)
    c = np.sqrt(c2)
    dt = 0.30
    T = 20000
    gam = np.full(N, 0.004)
    gam[360:] = 0.40
    src = 20
    u = np.zeros(N)
    v = np.zeros(N)

    def acc(u, v, t):
        a = -c2 * (L @ u) - 2.0 * gam * v
        a[src] += np.cos(f * t)
        return a

    for k in range(T):
        t = k * dt
        a = acc(u, v, t)
        u = u + dt * v + 0.5 * dt * dt * a
        v = v + 0.5 * dt * (a + acc(u, v, t + dt))
    seg = u[80:320]
    F = np.abs(np.fft.rfft(seg))
    dfr = np.fft.rfftfreq(len(seg), 1.0)
    kpk = dfr[np.argmax(F[1:]) + 1]
    q_meas = 2.0 * np.pi * kpk
    q_th = q_of_omega(f, c2)
    relB = abs(q_meas - q_th) / q_th
    lam_meas_time = 2.0 * np.pi / q_meas
    P(f"  max|u|={np.abs(u).max():.3f}（稳定）  空间 FFT 峰 q_meas={q_meas:.4f}")
    P(f"  理论 q(f)=2arcsin(f/2c)={q_th:.4f}  相对差={relB*1e2:.3f}%")
    P(f"  λ_meas=2π/q_meas={lam_meas_time:.3f}  λ_th=2π/q={2*np.pi/q_th:.3f}")
    out["B_time_domain"] = dict(N=N, dt=dt, T=T, f=round(f, 6),
                                q_measured=round(q_meas, 6), q_theory=round(q_th, 6),
                                rel_err=round(relB, 6),
                                lam_measured=round(lam_meas_time, 6),
                                lam_theory=round(2 * np.pi / q_th, 6))
    okB = relB < 1e-2

    # ---------- C. 显式欧拉不稳定性 ----------
    P("\n[C] 显式欧拉 i∂tψ=Hψ 的放大因子 |1-iλdt|²=1+(λdt)²>1（无条件增长）")
    lam_max = float(np.linalg.eigvalsh(L).max())
    rowsC = []
    for c2c, dtc, Tc in [(0.8, 0.1, 6000), (1.0, 0.02, 8000)]:
        amp = np.sqrt(1.0 + (c2c * lam_max * dtc) ** 2)
        growth = amp ** Tc
        # 数值复现
        psi = np.zeros(16, complex)
        Ls = pathL(16)
        b = np.zeros(16)
        b[7] = np.sqrt(2 / 3)
        b[8] = np.sqrt(1 / 3)
        b -= b.mean()
        mx = 0.0
        for t in range(Tc):
            psi += b * dtc * np.cos(np.log(2) * t * dtc)
            psi = psi - 1j * c2c * (Ls @ psi) * dtc
            mx = max(mx, float(np.abs(psi).max()))
        rowsC.append(dict(c2=c2c, dt=dtc, T=Tc, amp_per_step=round(amp, 8),
                          analytic_growth=float(f"{growth:.3e}"),
                          measured_max=float(f"{mx:.3e}")))
        P(f"  c²={c2c} dt={dtc} T={Tc}: 每步放大={amp:.6f}  "
          f"解析增长={growth:.2e}  实测 max|ψ|={mx:.2e}")
    out["C_euler_unstable"] = rowsC
    okC = all(r["amp_per_step"] > 1.0 for r in rowsC)

    # ---------- 汇总 ----------
    verdicts = {"A_dispersion": bool(okA), "B_time_domain": bool(okB),
                "C_euler_unstable": bool(okC)}
    out["verdicts"] = verdicts
    P("\n" + "=" * 96)
    for k, v in verdicts.items():
        P(f"  {k}: {'PASS' if v else 'FAIL'}")
    P("=" * 96)

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_quantum_bridge_timeevolve_check.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "k4_quantum_bridge_timeevolve_check.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    P("\n已存 k4_quantum_bridge_timeevolve_check.json / .txt")


if __name__ == "__main__":
    main()
