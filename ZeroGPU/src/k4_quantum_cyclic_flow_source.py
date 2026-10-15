#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_quantum_cyclic_flow_source.py —— 循环流稳态作源: 模流边流 J(t) 以 f 振荡驱动 L2 波
======================================================================================
验证链(与 B' 的唯一区别: 源是时变流, 不是静止稳态):
  Z087 稳态 ω=(2/3,1/3)  →  模流边流 J(t)=(ω_A−ω_B)·e^{if t},  f=log2 (Z027 §5.1 内禀时间)
  → 源(边流)以 f 振荡注入 L2:  i∂tψ = c²Lψ + J(t)(e_A − e_B)
  → 测波主导频率 vs f.  若波频率=f, 则"循环流→源以f振荡→波受迫f"整链成立.

物理: 粒子=流(权重经边流动, 分布锁定), 不是态(静止). 模流非平凡=边流时变.
"""
from __future__ import annotations
import os, json
import numpy as np
from numpy.fft import rfft, rfftfreq

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

def path_laplacian(n):
    L = np.zeros((n, n))
    for i in range(n):
        L[i, i] += (1 if i > 0 else 0) + (1 if i < n - 1 else 0)
        if i > 0: L[i, i-1] -= 1
        if i < n-1: L[i, i+1] -= 1
    return L

def run(n, c2, f, J0, T, dt, srcA, srcB, probe):
    """一阶量子行走 + 循环流边流源. 返回源处与探测处主导频率(FFT)."""
    L = path_laplacian(n)
    psi = np.zeros(n, complex)
    pa, pp = np.zeros(T, complex), np.zeros(T, complex)
    for t in range(T):
        tt = t * dt
        J = J0 * np.cos(f * tt)              # 模流边流: 权重差振荡, 节奏 f
        psi[srcA] += J * dt                   # 源 A 注入(权重流入 A 端)
        psi[srcB] -= J * dt                   # 源 B 抽走(权重流出 B 端) —— 守恒边流
        psi = psi - 1j * c2 * (L @ psi) * dt  # 一阶量子行走
        pa[t] = psi[srcA]
        pp[t] = psi[probe]
    def dom(x):
        fr = rfftfreq(T, dt); F = np.abs(rfft(x.real))
        return float(fr[np.argmax(F[1:]) + 1])
    return dom(pa), dom(pp)

def main():
    P_("=" * 100)
    P_("循环流稳态作源: 模流边流 J(t)=Δω·e^{if t} 驱动 L2 波, 测波频率 vs f")
    P_("=" * 100)
    out = {}

    f = np.log(2)                             # 内禀频率 log2 (Z087)
    wA, wB = 2/3, 1/3
    J0_nat = wA - wB                          # 稳态权重差(内生源强度)
    P_(f"稳态 ω_A={wA:.3f}, ω_B={wB:.3f}; 模流边流振幅 J0=Δω={J0_nat:.4f}")
    P_(f"内禀频率 f=log2={f:.4f} (Z087/Z027 §5.1); 源 J(t)=Δω·cos(f·t) 注入边 A→B")
    P_("")

    n, c2, T, dt = 16, 0.8, 6000, 0.1
    srcA, srcB = 7, 8
    probe = 12
    out["setup"] = {"steady": [wA, wB], "f_intrinsic": round(f, 5), "J0_natural": J0_nat,
                    "engine": f"i∂tψ=c²Lψ+J(e_A−e_B), path{n}, c²={c2}"}

    # 源耦合强度扫描: J0_nat(1/3) 到强耦合
    P_(f"图 path{n}, c²={c2}, 源边 ({srcA},{srcB}), 探测 {probe}; T={T}, dt={dt}")
    results = []
    for label, J0 in [("内生 J0=Δω", J0_nat), ("×3", J0_nat*3), ("×10", J0_nat*10)]:
        fA, fP = run(n, c2, f, J0, T, dt, srcA, srcB, probe)
        P_(f"  源耦合 {label:12s}: 源处主导频率={fA:.4f} (f_A/f={fA/f:.3f}), 探测处={fP:.4f} (f_P/f={fP/f:.3f})")
        results.append({"coupling": label, "J0": round(J0, 4),
                        "freq_at_source": round(fA, 4), "ratio_src": round(fA/f, 3),
                        "freq_at_probe": round(fP, 4), "ratio_probe": round(fP/f, 3)})
    out["sweep"] = results

    # 更大图/更长确认
    P_("\n【更大图确认】path32, 源中间两顶点")
    n2, srcA2, srcB2, probe2 = 32, 15, 16, 26
    fA, fP = run(n2, c2, f, J0_nat, T, dt, srcA2, srcB2, probe2)
    P_(f"  path32 源处={fA:.4f} (f_A/f={fA/f:.3f}), 探测处={fP:.4f} (f_P/f={fP/f:.3f})")
    out["path32"] = {"freq_at_source": round(fA, 4), "ratio_src": round(fA/f, 3),
                     "freq_at_probe": round(fP, 4), "ratio_probe": round(fP/f, 3)}

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_quantum_cyclic_flow_source.json"), "w", encoding="utf-8") as f_:
        json.dump(out, f_, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "k4_quantum_cyclic_flow_source.txt"), "w", encoding="utf-8") as f_:
        f_.write("\n".join(LOG))
    P_("\n已存 k4_quantum_cyclic_flow_source.json / .txt")

if __name__ == "__main__":
    main()
