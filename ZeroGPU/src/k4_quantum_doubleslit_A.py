#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_quantum_doubleslit_A.py —— 桥选项 A: L1 内禀频率内生驱动 → λ_e 实测
======================================================================================
Z117 判定: Z115/Z116 的 λ=π / k=√λ_max(K₄)=2 都是拿图谱补"能量从哪来"的洞, 都补错了.
选项 A(用户选定): 让 L1 记录层稳态内禀频率 f 成为 L2 图上波的驱动源,
                波长 λ_e 由图上传播实测决定(不预设 λ=2π/f 的解析关系, 关系是发现不是假设).

桥结构(Zero 自洽):
  L1 稳态 ω_v=record_v/Σ → 内禀频率 f=log(p/q)=log2 (Z087, geom seed, path8-32 普适)
    → 假设 L1 内禀频率 f 进入 L2 波为角频率 ω=2c·sin(q/2) (图上行模色散)
    → q = 2·arcsin(f/2c)  (c=1, 能量单位选择)
    → 单频散射 (E−H+iΣ)ψ=b 驱动, λ_e 实测(源行相位梯度)

关键检验: λ_e(实测) 与 2π/f 的关系 —— 若 λ_e≈2π/f 则"物质波"直觉成立;
         若否, 关系由图色散决定(λ_e=2π/q_from_f), 是图上特有的桥.
"""
from __future__ import annotations
import os, json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

# ---------------- 复用 v2 引擎 ----------------
def build_H(nr, nc, xb=90, apertures=(90, 150), no_barrier=False, c2=1.0):
    N = nr * nc; rows = []; cols = []; vals = []; deg = np.zeros(N)
    def ae(a, b):
        rows.extend([a, b]); cols.extend([b, a]); vals.extend([-1.0, -1.0]); deg[a] += 1; deg[b] += 1
    for y in range(nr):
        for x in range(nc):
            v = y * nc + x
            if x + 1 < nc and (no_barrier or not (x == xb - 1 and y not in apertures)):
                ae(v, v + 1)
            if y + 1 < nr:
                ae(v, v + nc)
    rows += list(range(N)); cols += list(range(N)); vals += list(deg)
    return sp.csr_matrix((vals, (rows, cols)), shape=(N, N)) * c2

def absorber(nr, nc, m=24, sig0=1.5):
    yy, xx = np.mgrid[0:nr, 0:nc]
    de = np.minimum(np.minimum(yy, nr - 1 - yy), np.minimum(xx, nc - 1 - xx))
    return np.where(de < m, sig0 * ((m - de) / m) ** 2, 0.0)

def solve_scatter(nr, nc, q, xb=90, apertures=(90, 150), no_barrier=False,
                  src=(120, 20), c2=1.0, m=24, sig0=1.5):
    H = build_H(nr, nc, xb, apertures, no_barrier, c2)
    N = nr * nc
    E = 2 * c2 * (1 - np.cos(q))
    b = np.zeros(N, complex); b[src[0] * nc + src[1]] = 1.0
    A = (E * sp.eye(N, format="csr") - H + 1j * sp.diags(absorber(nr, nc, m, sig0).ravel())).tocsc()
    return spla.splu(A).solve(b).reshape(nr, nc)

def q_from_freq(f, c2=1.0):
    """L1 内禀频率 f → 驱动波数 q: 假设 f 为角频率 ω, 色散 ω=c·2sin(q/2) (行模 λ=4sin²(q/2))"""
    w = f / c2
    return float(np.arcsin(w / 2) * 2) if w <= 2 else float("nan")

def q_of_E(E, c2=1.0):
    lam = E / c2
    return float(np.arccos(1 - lam / 2)) if 0 <= lam <= 4 else float("nan")

def measure_lambda(nr, nc, q):
    """对照跑(无挡板)测源行相位梯度 → λ_e 实测"""
    pc = solve_scatter(nr, nc, q, no_barrier=True)
    y0 = nr // 2
    xw = np.arange(30, 91)
    qm = float(np.polyfit(xw, np.unwrap(np.angle(pc[y0, xw])), 1)[0])
    return abs(qm), 2 * np.pi / abs(qm)

def main():
    P_("=" * 96)
    P_("桥选项A: L1 内禀频率 f 内生驱动 → λ_e 实测 (k4_quantum_doubleslit_A)")
    P_("=" * 96)
    nr, nc = 240, 260; xb = 90; ap = (90, 150)
    out = {"bridge": "L1 f -> L2 wave, λ_e measured (not assumed 2π/f)"}

    # L1 内禀频率: Z087 geom seed → log2 (跨 path8-32 普适); K₄ 三槽 Z089
    f_log2 = np.log(2)
    k4_3slots = [0.405, 0.693, 0.981]      # Z089 trip seed: [log1.5, log2, log(1.5*~)...]
    candidates = {"log2(geom, L1实测)": f_log2,
                  "K4槽1(0.405)": 0.405, "K4槽2(0.693)": 0.693, "K4槽3(0.981)": 0.981}
    P_(f"L1 内禀频率候选: log2={f_log2:.4f}; K₄ 三槽 {k4_3slots}")
    P_("")

    rows = {}
    for name, f in candidates.items():
        q_in = q_from_freq(f)
        if np.isnan(q_in):
            P_(f"[{name}] f={f:.3f} 超出色散范围 (f/2c>1), 跳过"); continue
        qm, lam_meas = measure_lambda(nr, nc, q_in)
        lam_2pi_f = 2 * np.pi / f
        ratio = lam_meas / lam_2pi_f
        P_(f"[{name}] f={f:.4f}")
        P_(f"    ω=f 假设 → q_in=2·arcsin(f/2)={q_in:.4f}")
        P_(f"    实测: q_meas={qm:.4f} (误差 {abs(qm-q_in)/q_in*100:.2f}%), λ_e,实测={lam_meas:.3f} 格")
        P_(f"    2π/f={lam_2pi_f:.3f} → λ_e/2πf={ratio:.3f}  (≈1 ⟹ 物质波直觉; 否则图上色散桥)")
        rows[name] = {"f": round(f, 4), "q_in": round(q_in, 4), "q_measured": round(qm, 4),
                      "lambda_measured": round(lam_meas, 3), "lambda_2pi_over_f": round(lam_2pi_f, 3),
                      "ratio_lambda_vs_2pif": round(ratio, 3)}
        P_("")
    out["candidates"] = rows

    # 双缝验证: 用 log2 的 q 驱动, 测条纹间距
    P_("【双缝】用 log2 驱动 (q=" + f"{q_from_freq(f_log2):.4f}" + ") 测条纹间距")
    p12 = solve_scatter(nr, nc, q_from_freq(f_log2), xb, ap)
    p1 = solve_scatter(nr, nc, q_from_freq(f_log2), xb, (ap[0],))
    p2 = solve_scatter(nr, nc, q_from_freq(f_log2), xb, (ap[1],))
    y0 = nr // 2
    dsl = {}
    for L in [40, 80, 120]:
        x = xb + L
        ykw = np.arange(y0 - 60, y0 + 61)
        ph = np.unwrap(np.angle(p1[ykw, x] / p2[ykw, x]))
        A = np.vstack([ykw, np.ones_like(ykw)]).T
        cf, *_ = np.linalg.lstsq(A, ph, rcond=None)
        per = 2 * np.pi / abs(cf[0])
        dsl[L] = round(float(per), 3)
        P_(f"  L={L:3d}: 实测条纹间距={per:.3f} 格")
    out["doubleslit_log2"] = {"q_drive": round(q_from_freq(f_log2), 4), "fringe_spacing": dsl}
    P_(f"  λ_e,实测={rows['log2(geom, L1实测)']['lambda_measured']:.3f} → 条纹间距 ∝ λ_e·L/d 检验")

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_quantum_doubleslit_A.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "k4_quantum_doubleslit_A.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(LOG))
    P_("\n已存 k4_quantum_doubleslit_A.json / .txt")

if __name__ == "__main__":
    main()
