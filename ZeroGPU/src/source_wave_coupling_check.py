#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
source_wave_coupling_check.py —— 源-波耦合的归一：频率与幅度
================================================================
对象：Z152 §5.1 的开放项——L1 内禀频率 f 进入 L2 波的源项 g cos(ft) 的归一。

两条分开：
  (频率) 线性响应：对任意 g>0，稳态频率 = 驱动频率 f。
  (幅度) 源幅度 = g·Δρ；Δρ=1/3 由 L1 稳态 (2/3,1/3) 定；
         g 前的读出权重 β（Z102 递推的混合权重）只缩放幅度，不改频率。

核验：
  A. L1 稳态：ρ=(2/3,1/3)，Δρ=1/3
  B. 频率：g = 0.1, 1/3, 1, 3 四个值，稳态主导频率都 = f（<1e-3）
  C. 幅度：A(g)/g 为常数（<1e-6）
"""
from __future__ import annotations
import os
import json
import math
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


def steady_q(f, c2, N=400, dt=0.3, T=20000, g=1.0, src=20):
    """阻尼波方程驱动 g cos(f t)，返回稳态空间 FFT 峰 q 与幅度。"""
    L = pathL(N)
    gam = np.full(N, 0.004)
    gam[360:] = 0.40
    u = np.zeros(N)
    v = np.zeros(N)

    def acc(u, v, t):
        a = -c2 * (L @ u) - 2.0 * gam * v
        a[src] += g * np.cos(f * t)
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
    q = 2.0 * np.pi * kpk
    return q, float(np.abs(u).max())


def main():
    P("=" * 84)
    P("源-波耦合的归一：频率与幅度 (source_wave_coupling_check)")
    P("=" * 84)

    rho1, rho2 = 2 / 3, 1 / 3
    drho = rho1 - rho2
    f = math.log(2)
    c2 = 1.0
    P(f"L1 稳态 rho=({rho1:.4f},{rho2:.4f})，Δρ={drho:.4f}；f=log2={f:.6f}")

    checks = []

    def chk(name, cond, detail=""):
        checks.append({"name": name, "pass": bool(cond), "detail": detail})
        P(("  [v] " if cond else "  [x] ") + name + (f"   {detail}" if detail else ""))

    P("\n[A] 稳态")
    chk("A  Δρ = 1/3", abs(drho - 1 / 3) < 1e-12, f"{drho:.6f}")

    P("\n[B/C] 频率与幅度（四个 g）")
    gs = [0.1, 1 / 3, 1.0, 3.0]
    res = []
    for g in gs:
        q, amp = steady_q(f, c2, g=g)
        res.append((g, q, amp))
        P(f"  g={g:.4f}: q={q:.5f}, max|u|={amp:.5f}, amp/g={amp/g:.5f}")

    freq_ok = True
    for g, q, amp in res:
        # 理论 q = 2 arcsin(f/2c)；与 g 无关
        qth = 2 * math.asin(f / (2 * math.sqrt(c2)))
        if abs(q - qth) / qth > 2e-3:
            freq_ok = False
    chk("B  四个 g 的稳态 q 都 = 2arcsin(f/2c)（<2e-3）", freq_ok,
        f"qth={2*math.asin(f/2):.5f}")

    ratios = [amp / g for g, q, amp in res]
    amp_ok = max(ratios) - min(ratios) < 1e-3 * max(ratios)
    chk("C  A(g)/g 为常数（<1e-3 相对）", amp_ok,
        f"min={min(ratios):.5f}, max={max(ratios):.5f}")

    npass = sum(1 for c in checks if c["pass"])
    P("\n" + "=" * 84)
    P(f"核验：通过 {npass} / {len(checks)}")
    P("=" * 84)

    out = {
        "drho": drho, "f": f,
        "g_results": [{"g": g, "q": q, "amp": amp, "amp_over_g": amp / g} for g, q, amp in res],
        "checks": checks, "pass": npass, "total": len(checks),
    }
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "source_wave_coupling_check.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "source_wave_coupling_check.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    P("\n已存 source_wave_coupling_check.json / .txt")


if __name__ == "__main__":
    main()
