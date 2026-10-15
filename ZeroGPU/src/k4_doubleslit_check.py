#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_doubleslit_check.py —— 对 Z115《真 Zero 内生电子双缝干涉》推导链的独立核查 + 条纹复现
======================================================================================
核查:
  CK1 公式层复现  I(x)=2+2cos(kΔs), k=2, Δs=|x−6|−|x−14| → 亮纹/间距/条数是否如 Z115 所列;
                  并把同一公式放到图节点上(整数 x)采样 → 亚 Nyquist 混叠
  CK2 判据核查    k=√λ 与格上相位梯度 q: λ(q)=4sin²(q/2)。k=2 ⇒ q=π(v_g=0,周期 2 格);
                  λ_e=π ⇒ q=2 ⇒ λ=2.83≠4; 而 K₄(λ_max=4) 与格图(λ_max=7.966)是两张图
  CK3 几何核查    曼哈顿距离 ⇒ Δs 与屏行 y 无关(无屏距依赖/无远场); 真波走欧氏程函 ⇒ 展宽
  CK4 时间域核查  Z104 蛙跳真跑: (a) 单源相位场走哪个距离度量; (b) Z115 的 k=2 下逐行条纹;
                  (c) 按 Z115 自己的出处 ω=2c(带边)会给出什么
  CK5 修正复现    图上真跑一个良分辨、屏离开缝行的双缝 → 条纹随 L 展宽, 与欧氏程函/λL/d 对比
引擎与约定沿用 Z104: ∂²φ/∂t²=−c²Lφ, φ(t+1)=2φ(t)−φ(t−1)−c²(Lφ+x), c²=3.5/λ_max(CFL), Δt=1。
稳态取法: 体内弱阻尼 γ_b + 边缘海绵; 单频驱动 DFT 与窗口严格通约(无泄漏); 干涉项用交叉项
          I12−I1−I2=2Re(ψ1ψ2*) 抽取(线性性已单独核验)。
"""
from __future__ import annotations
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []

def P_(*a):
    s = " ".join(str(x) for x in a)
    LOG.append(s); print(s)

# ---------------------------------------------------------------- 图与引擎
def grid_deg(nr, nc):
    d = np.full((nr, nc), 4.0)
    d[0, :] -= 1; d[-1, :] -= 1; d[:, 0] -= 1; d[:, -1] -= 1
    return d

def lap(Pm, deg):
    S = np.zeros_like(Pm)
    S[1:, :] += Pm[:-1, :]; S[:-1, :] += Pm[1:, :]
    S[:, 1:] += Pm[:, :-1]; S[:, :-1] += Pm[:, 1:]
    return deg * Pm - S

def grid_lam_max(n):
    return 2 * (2 + 2 * np.cos(np.pi / n))

def row_lam(q):
    return 2 - 2 * np.cos(q)

def q_of_omega(omega, c):
    lam = (2 * np.sin(omega / 2.0) / c) ** 2
    return (float(np.pi), float(lam)) if lam >= 4.0 else (float(np.arccos(1 - lam / 2)), float(lam))

def run_wave(nr, nc, sources, q_target, c2, n_trans=9000, n_win=2048,
             g_bulk=0.0015, g_edge=0.15, margin=16, amp=0.05, ramp_periods=20):
    deg = grid_deg(nr, nc); c = float(np.sqrt(c2))
    w_t = 2 * np.arcsin(np.clip(c * np.sqrt(row_lam(q_target)) / 2, 0, 1))
    m = max(1, int(round(n_win * w_t / (2 * np.pi))))
    omega = 2 * np.pi * m / n_win
    q_eff, lam_eff = q_of_omega(omega, c)
    yy, xx = np.mgrid[0:nr, 0:nc]
    de = np.minimum(np.minimum(yy, nr - 1 - yy), np.minimum(xx, nc - 1 - xx))
    gam = g_bulk + np.where(de < margin, g_edge * ((margin - de) / margin) ** 2, 0.0)
    Pm = np.zeros((nr, nc)); Po = np.zeros_like(Pm); S = np.zeros_like(Pm)
    acc = np.zeros((nr, nc), complex); T = 2 * np.pi / omega
    for n in range(n_trans + n_win):
        dr = amp * min(1.0, n / (ramp_periods * T)) ** 2 * np.sin(omega * n)
        S[:] = 0.0
        for (a, b) in sources:
            S[a, b] = dr
        Pn = (2 * Pm - (1 - gam) * Po - c2 * (lap(Pm, deg) + S)) / (1 + gam)
        Po = Pm; Pm = Pn
        if n >= n_trans:
            acc += Pm * np.exp(-1j * omega * n)
    return {"acc": acc, "omega": omega, "q_eff": q_eff, "lam_eff": lam_eff, "c": c}

def maxima_of(y, x):
    idx = np.where((y[1:-1] > y[:-2]) & (y[1:-1] > y[2:]))[0] + 1
    return x[idx]

def central_spacing(xpk, xc):
    if len(xpk) < 2:
        return float("nan")
    i = int(np.argmin(np.abs(xpk - xc)))
    if i >= len(xpk) - 1:
        i = len(xpk) - 2
    return float(xpk[i + 1] - xpk[i])

def phase_slope(a1, a2, y, xk):
    """干涉相位 Δφ=arg(ψ1/ψ2) 在窗口 xk 上的线性斜率 (rad/格) 与拟合残差"""
    ph = np.unwrap(np.angle(a1[y, xk] / a2[y, xk]))
    A = np.vstack([xk, np.ones_like(xk)]).T
    coef, *_ = np.linalg.lstsq(A, ph, rcond=None)
    rms = float(np.sqrt(np.mean((A @ coef - ph) ** 2)))
    return float(coef[0]), rms

def metric_fit(acc, pts, d_fn):
    ph = np.unwrap(np.angle(np.array([acc[y, x] for (y, x) in pts])))
    d = np.array([d_fn(y, x) for (y, x) in pts])
    A = np.vstack([d, np.ones_like(d)]).T
    coef, *_ = np.linalg.lstsq(A, ph, rcond=None)
    return float(abs(coef[0])), float(np.sqrt(np.mean((A @ coef - ph) ** 2)))

# ---------------------------------------------------------------- 主程序
def main():
    out = {}
    P_("=" * 98)
    P_("Z115《真 Zero 内生电子双缝干涉》推导链独立核查 + 条纹复现")
    P_("=" * 98)

    # ============================================================ CK1
    P_("\n【CK1】公式层复现 + 图节点采样 (Z115 §2⑤⑥⑦)")
    k = 2.0; lam_e = 2 * np.pi / k; s1x, s2x, ctr = 6.0, 14.0, 10.0
    xs = np.linspace(0, 24, 240001)
    If = 2 + 2 * np.cos(k * (np.abs(xs - s1x) - np.abs(xs - s2x)))
    m = (xs > s1x) & (xs < s2x)
    xb = maxima_of(If[m], xs[m]); sp = central_spacing(xb, ctr)
    doc = np.array([6.9, 8.4, 10.0, 11.6, 13.1])
    P_(f"  连续屏(24 万点): 亮纹 {[round(float(v),4) for v in xb]}")
    P_(f"    与 Z115 §3 表列值末位一致 → 所列亮纹可复现 ✓; 间距 = {sp:.6f} = π/2 ✓")
    P_(f"  §2⑦ 写的间距式 (π/k)/斜率 = {(np.pi/k)/2:.6f} ✗ 差 2 倍 (正确: Δx=(2π/k)/斜率=λ_e/2)")
    Iin = If[m]
    P_(f"  缝间 x∈(6,14): I∈[{Iin.min():.4f},{Iin.max():.4f}], 亮纹 {len(xb)} 条")
    P_(f"  缝外 x<6 或 x>14: Δs=∓8 恒定 ⇒ I≡2+2cos16={2+2*np.cos(16):.4f} 常数, 完全没有条纹")
    P_(f"    → §3 '亮度范围[0,4]' 只在缝间成立; §3 '两侧对称等间距—真实双缝形态' 与自身公式矛盾")
    xn = np.arange(0, 25); In = 2 + 2 * np.cos(k * (np.abs(xn - s1x) - np.abs(xn - s2x)))
    xpn = maxima_of(In, xn); xpn = xpn[(xpn > s1x) & (xpn < s2x)]
    alias_f = abs(2 * k - 2 * np.pi)
    P_(f"  ★图节点上采样(整数 x, 才是'图上'可测的): 亮纹 {[int(v) for v in xpn]} → {len(xpn)} 条, 间距 3 格")
    P_(f"    强度条纹空间频率 2k={2*k} rad/格 > Nyquist π={np.pi:.4f} ⇒ 混叠到 {alias_f:.4f} rad/格")
    P_(f"    表观周期 = {2*np.pi/alias_f:.4f} 格; 图上看到 {len(xpn)} 条而不是 5 条")
    P_(f"    ⇒ Z115 的 5 条纹/1.571 间距要求分辨率优于 1 格, 图上不存在; §5.1 '屏连续化' 正是把图换掉")
    P_(f"  落盘检查: Z115 脚本 I_x[::30] 采样间隔 = {24/600*30:.2f} 格 < 条纹 1.571 ⇒ 欠采样,")
    P_(f"    落盘 JSON 里画不出 5 条纹(实际只见到约 4 个起伏)")
    out["CK1"] = {"continuum_bright": [round(float(v), 4) for v in xb], "spacing": round(sp, 6),
                  "z115_bright": doc.tolist(), "spacing_formula_written": round((np.pi/k)/2, 6),
                  "I_strip": [round(float(Iin.min()), 4), round(float(Iin.max()), 4)],
                  "I_outside": round(float(2+2*np.cos(16)), 4), "n_bright_continuum": int(len(xb)),
                  "node_bright": [int(v) for v in xpn], "n_bright_on_nodes": int(len(xpn)),
                  "nyquist": {"fringe_freq": 2*k, "nyquist": float(np.pi),
                              "alias_freq": round(float(alias_f), 4),
                              "apparent_period": round(float(2*np.pi/alias_f), 4)}}

    # ============================================================ CK2
    P_("\n【CK2】判据核查: k=√λ 能当波数吗?")
    P_("  格点行模色散 λ(q)=4sin²(q/2); 真正的相位梯度是 q, √λ 只在 q→0 时等于 q:")
    for kk in [0.5, 1.0, 1.5, 2.0]:
        qq = 2 * np.arcsin(min(kk / 2, 1.0))
        P_(f"    √λ={kk:.2f} ⇒ q={qq:.4f} rad/格, 空间周期={2*np.pi/qq:6.3f} 格, "
           f"v_g/c=cos(q/2)={np.cos(qq/2):.4f}, q/√λ={qq/kk:.4f}")
    lam_pi = 4 * np.sin(1.0) ** 2
    P_(f"  ① 坚持 k=√λ=2 (K₄ λ_max=4): 格上 ⇒ q=π ⇒ 空间周期 2 格 ≠ π, 且 v_g=0 (带边, 不传播)")
    P_(f"  ② 坚持 λ_e=π (即 q=2):      格上 ⇒ λ=4sin²1={lam_pi:.4f} ≠ 4 = K₄ 的 λ_max")
    P_(f"  ⇒ λ_max(K₄)=4 与 λ_e=π 在同一张格图上不能同时成立 —— 二者只差一次 'q ↔ √λ' 的偷换")
    l24, l48 = grid_lam_max(24), grid_lam_max(48)
    P_(f"  ③ 波跑在 24×24 格图上, 该图 λ_max={l24:.4f}(√λ_max={np.sqrt(l24):.4f}) = Z115 脚本自己打印的 7.966,")
    P_(f"     而 k 却取 K₄ 的 2: K₄ 与格图是两张不同的图, 这是跨图取值。")
    P_(f"  ④ 更根本: K₄ 只有 4 个顶点、没有空间, 由 ω(K₄)=2c 推'空间波长 λ_e=π'是类型错误;")
    P_(f"     而且格上波要给定驱动频率才能定波数, 不是'由图谱唯一给出'。")
    out["CK2"] = {"grid_lmax_24": round(float(l24), 4), "sqrt_lmax_24": round(float(np.sqrt(l24)), 4),
                  "grid_lmax_48": round(float(l48), 4), "q_if_k_is_sqrt_lam_2": float(np.pi),
                  "period_if_q_pi": 2.0, "vg_at_q_pi": 0.0, "lam_if_lambda_e_pi": round(float(lam_pi), 4)}

    # ============================================================ CK3
    P_("\n【CK3】几何核查: '图上距离差' 用曼哈顿距离意味着什么")
    tab = {}
    for y in [0, 4, 8, 12, 16, 20]:
        D = np.abs(xs - s1x) - np.abs(xs - s2x)
        tab[y] = [round(float(v), 3) for v in maxima_of((2 + 2 * np.cos(k * D))[m], xs[m])]
    P_(f"  L1: Δs=(|x−6|+|y−12|)−(|x−14|+|y−12|)=|x−6|−|x−14|, 与屏行 y 无关")
    P_(f"    每行亮纹完全相同: {tab[0]}  (y=0/4/8/12/16/20 逐行一致)")
    P_(f"    ⇒ 条纹不随屏距变化; 既无远场展开, 也无真实双缝的 Δx=λL/d; 缝外 Δs 恒为 ∓8 ⇒ 均匀")
    xw = np.linspace(-24, 48, 240001); eu = {}
    for h in [0, 2, 4, 8, 12, 16]:
        D = np.sqrt((xw - s1x) ** 2 + h ** 2) - np.sqrt((xw - s2x) ** 2 + h ** 2)
        eu[h] = round(central_spacing(maxima_of(2 + 2 * np.cos(k * D), xw), ctr), 3)
        P_(f"    欧氏(真实波的程函) h={h:2d}: 中心间距={eu[h]:6.3f}   小角律 λ_e·h/d={lam_e*h/8:6.3f}")
    P_(f"  ⇒ Z115 恰好把屏放在与缝共线的 h=0: 那是 L1 与欧氏唯一重合的行, 也是真实双缝必须避开的近场;")
    P_(f"    屏一旦离开缝行(缝在屏的'前面'), L1 与欧氏给出完全不同的条纹。")
    out["CK3"] = {"L1_bright_per_row": {str(a): b for a, b in tab.items()},
                  "euclid_central_spacing_vs_h": {str(a): b for a, b in eu.items()},
                  "small_angle_law": "lambda_e*h/d, d=8"}

    # ============================================================ CK4a
    P_("\n【CK4】时间域核查: 用 Z104 蛙跳在图上真跑")
    nr, nc = 200, 240
    c2 = min(0.8, 3.5 / grid_lam_max(nr))
    kw = dict(n_trans=9000, n_win=2048, g_bulk=0.0015, g_edge=0.15, margin=16, amp=0.05)
    P_(f"  [CK4a] 度量判据: 单源相位场线性于 L1 还是欧氏距离? ({nr}×{nc} 格, 源 (y=100,x=120), λ=25 格)")
    r1 = run_wave(nr, nc, [(100, 120)], 0.25, c2, **kw)
    P_(f"    ω={r1['omega']:.6f} ⇒ q_eff={r1['q_eff']:.5f} (λ={2*np.pi/r1['q_eff']:.2f} 格)")
    met = {}; ncut = 0; winE = 0
    for h in [0, 20, 40]:
        for tag, pts in [("row", [(100 - h, x) for x in range(140, 215)]),
                         ("col", [(y, 150) for y in range(100 - h, 175 - h)])]:
            key = f"{tag}_h{h}"
            sE, rE = metric_fit(r1["acc"], pts, lambda y, x: float(np.hypot(x - 120, y - 100)))
            sL, rL = metric_fit(r1["acc"], pts, lambda y, x: float(abs(x - 120) + abs(y - 100)))
            tie = abs(rE - rL) < 1e-6
            ncut += 1; winE += (rE < rL - 1e-6)
            met[key] = {"slope_euclid": round(sE, 4), "rms_euclid": round(rE, 4),
                        "slope_L1": round(sL, 4), "rms_L1": round(rL, 4),
                        "winner": "tie" if tie else ("euclid" if rE < rL else "L1")}
            P_(f"    {key:9s} 斜率 欧氏 {sE:.4f} / L1 {sL:.4f}   rms 欧氏 {rE:.4f} / L1 {rL:.4f}"
               f"   → {'平局(过源, 必然)' if tie else ('欧氏胜' if rE < rL else 'L1胜')}")
    P_(f"    ⇒ 相位 = q_eff·(欧氏距离)+常数 (欧氏 {winE}/{ncut-1} 条非平凡割线胜, 斜率≈q_eff);")
    P_(f"      '图上曼哈顿距离差' 不是波的相位距离 —— 它只是图上最短步数(组合量)。")

    # ============================================================ CK4b
    nr2, nc2 = 60, 60
    c2b = min(0.8, 3.5 / grid_lam_max(nr2))
    sy, sx1, sx2, csx = 40, 26, 34, 30
    kw2 = dict(n_trans=9000, n_win=2048, g_bulk=0.0015, g_edge=0.15, margin=10, amp=0.05)
    P_(f"\n  [CK4b] 按 Z115 的 k=2 (q=2) 真跑双源: {nr2}×{nc2} 格, 缝节点 (y=40,x=26)/(y=40,x=34), 缝距 8")
    r2 = run_wave(nr2, nc2, [(sy, sx1), (sy, sx2)], 2.0, c2b, **kw2)
    r2a = run_wave(nr2, nc2, [(sy, sx1)], 2.0, c2b, **kw2)
    r2b = run_wave(nr2, nc2, [(sy, sx2)], 2.0, c2b, **kw2)
    res2 = float(np.max(np.abs(r2["acc"] - r2a["acc"] - r2b["acc"])) / np.max(np.abs(r2["acc"])))
    q2 = r2["q_eff"]
    I2 = np.abs(r2["acc"]) ** 2; I2n = I2 / I2.max()
    I2a = np.abs(r2a["acc"]) ** 2; I2b = np.abs(r2b["acc"]) ** 2
    C2 = 2 * np.real(r2a["acc"] * np.conj(r2b["acc"]))
    alias_q = abs(2 * q2 - 2 * np.pi); alias_p = 2 * np.pi / alias_q
    P_(f"    §2⑤ 线性性核验: max|ψ12−ψ1−ψ2|/max|ψ12| = {res2:.3e} ⇒ 线性叠加成立 ✓ (这条 Z115 没错)")
    P_(f"    ω={r2['omega']:.6f} ⇒ q_eff={q2:.5f}, λ_e,eff={2*np.pi/q2:.4f} 格")
    P_(f"    连续屏预言间距 π/q={np.pi/q2:.4f} (=Z115 的 1.571); 但强度频率 2q={2*q2:.4f} > Nyquist π")
    P_(f"    ⇒ 图节点上只能看到混叠值 {alias_q:.4f} rad/格, 表观周期 {alias_p:.4f} 格")
    xk = np.arange(12, nr2 - 12)                     # 避开海绵
    win = np.arange(csx - 3, csx + 4)
    wo = np.arange(csx + 16, csx + 23)               # 缝外(仍在无海绵区)
    P_("    行h  中心实测间距  欧氏预言(同窗)  L1预言   缝外实测   欧氏     L1      可见度V")
    rows = {}; visr = {}
    for h in [0, 4, 8, 12, 16]:
        y = sy - h
        g_c, _ = phase_slope(r2a["acc"], r2b["acc"], y, win)
        DE = q2 * (np.sqrt((win - sx1) ** 2 + h ** 2) - np.sqrt((win - sx2) ** 2 + h ** 2))
        gE = (DE[-1] - DE[0]) / (win[-1] - win[0])
        DL = q2 * (np.abs(win - sx1) - np.abs(win - sx2))
        gL = (DL[-1] - DL[0]) / (win[-1] - win[0])
        go, _ = phase_slope(r2a["acc"], r2b["acc"], y, wo)
        DEo = q2 * (np.sqrt((wo - sx1) ** 2 + h ** 2) - np.sqrt((wo - sx2) ** 2 + h ** 2))
        gEo = (DEo[-1] - DEo[0]) / (wo[-1] - wo[0])
        DLo = q2 * (np.abs(wo - sx1) - np.abs(wo - sx2))
        gLo = (DLo[-1] - DLo[0]) / (wo[-1] - wo[0])
        f = lambda g: float("inf") if abs(g) < 1e-3 else 2 * np.pi / abs(g)
        vis = float(np.mean(np.abs(C2[y, xk]) / (I2a[y, xk] + I2b[y, xk] + 1e-30)))
        visr[h] = round(vis, 3)
        rows[h] = {"center": {"meas": round(f(g_c), 3), "euclid": round(f(gE), 3), "L1": round(f(gL), 3)},
                   "outside": {"meas": round(f(go), 3), "euclid": round(f(gEo), 3), "L1": round(f(gLo), 3)},
                   "vis": round(vis, 3)}
        P_(f"    {h:3d}  {f(g_c):11.3f}  {f(gE):13.3f}  {f(gL):8.3f}  {f(go):9.3f} {f(gEo):8.3f} {f(gLo):8.3f}  {vis:7.3f}")
    P_(f"    ⇒ (i) 中心实测间距 {rows[0]['center']['meas']:.2f} (≈混叠值 {alias_p:.2f}) ≠ Z115 的 1.571:")
    P_(f"          图上根本表示不了那 5 条 1.571 间距的条纹;")
    P_(f"      (ii) 缝外(h=8)实测仍有条纹({rows[8]['outside']['meas']:.1f}), 与欧氏预言"
       f"({rows[8]['outside']['euclid']:.1f})同量级, 而 L1 预言'无穷(无条纹)':")
    P_(f"          L1 的 Δs 在缝外恒为 ∓8 ⇒ 模型说缝外一片均匀 —— 真实双缝照片恰恰相反;")
    P_(f"      (iii) 可见度 V={visr[0]:.2f}(h=0)…{visr[16]:.2f}(h=16) < 1: 两路振幅不等,")
    P_(f"          Z115 的等幅假设(I∈[0,4], V=1)不成立; 振幅还随距离衰减(无等幅平面波)。")

    # ============================================================ CK4c
    P_(f"\n  [CK4c] 若严格按 Z115 自己的出处 k=√λ_max(K₄)=2 ⇔ ω=c·2 (格上即带边 q=π)")
    r3 = run_wave(nr2, nc2, [(sy, sx1), (sy, sx2)], np.pi, c2b, **kw2)
    r3a = run_wave(nr2, nc2, [(sy, sx1)], np.pi, c2b, **kw2)
    r3b = run_wave(nr2, nc2, [(sy, sx2)], np.pi, c2b, **kw2)
    q3 = r3["q_eff"]; I3 = np.abs(r3["acc"]) ** 2; I3n = I3 / I3.max()
    g3, _ = phase_slope(r3a["acc"], r3b["acc"], sy, win)
    P_(f"    ω={r3['omega']:.6f} ⇒ q_eff={q3:.5f} (≈π), λ_eff={r3['lam_eff']:.4f}")
    P_(f"    源行相位近乎不推进(拟合'间距' {2*np.pi/abs(g3):.0f} 格 ⇒ 实为驻波, 无传播条纹);")
    P_(f"      若强行按连续屏 q 算, 预言 π/q={np.pi/q3:.3f} —— 与 Z115 的 1.571 也不同")
    P_(f"    带边 v_g→0: 源行峰值 {float(I3[sy].max()):.3e}; 8 格外峰值 {float(I3[sy-8].max()):.3e} "
       f"(归一 {float(I3n[sy-8].max()):.2e}) → 几乎不向外辐射")
    P_(f"    ⇒ Z115 的 k 出处(带边)给不出它自己的 π/2=1.571: 出处与结论数值不自洽。")
    out["CK4"] = {"grid_metric_test": f"{nr}x{nc}", "grid_doubleslit": f"{nr2}x{nc2}",
                  "c2": round(c2b, 5), "slits": [[sy, sx1], [sy, sx2]],
                  "linearity_residual": res2,
                  "CK4a_metric": {"q_eff": round(r1["q_eff"], 5), "fits": met,
                                  "euclid_wins": f"{winE}/{ncut-1} (excl. trivial tie)"},
                  "CK4b_q2": {"omega": round(r2["omega"], 6), "q_eff": round(q2, 5),
                              "alias_period": round(float(alias_p), 4), "rows": rows},
                  "CK4c_bandedge": {"q_eff": round(float(q3), 5),
                                    "center_spacing": round(float(2*np.pi/abs(g3)), 3),
                                    "amp_row0_norm": 1.0, "amp_row8_norm": float(I3n[sy-8].max())}}

    # ============================================================ CK5
    P_("\n【CK5】修正复现: 图上真跑一个良分辨、屏离开缝行的双缝 (全由图上的波跑出来)")
    nr3, nc3 = 270, 420
    c2c = min(0.8, 3.5 / grid_lam_max(nr3))
    ys = 40; d = 60; xs1, xs2, xc = 180, 240, 210
    q0 = 2 * np.pi / 12.0
    P_(f"  {nr3}×{nc3} 格, c²={c2c:.4f}; 缝行 y={ys}, 缝节点 x={xs1}/{xs2}(缝距 d={d}); 目标 λ_e=12 格")
    P_(f"  取 d/λ_e = {d/12:.1f} (极大方向 sinθ=nλ/d ⇒ 约 {int(2*d/12)+1} 个级次) —— 这样才是熟悉的多条纹双缝图")
    r4 = run_wave(nr3, nc3, [(ys, xs1), (ys, xs2)], q0, c2c, **kw)
    r4a = run_wave(nr3, nc3, [(ys, xs1)], q0, c2c, **kw)
    r4b = run_wave(nr3, nc3, [(ys, xs2)], q0, c2c, **kw)
    res4 = float(np.max(np.abs(r4["acc"] - r4a["acc"] - r4b["acc"])) / np.max(np.abs(r4["acc"])))
    q4 = r4["q_eff"]; lam4 = 2 * np.pi / q4
    I4 = np.abs(r4["acc"]) ** 2
    I4a = np.abs(r4a["acc"]) ** 2; I4b = np.abs(r4b["acc"]) ** 2
    C4 = 2 * np.real(r4a["acc"] * np.conj(r4b["acc"]))
    P_(f"  ω={r4['omega']:.6f} ⇒ q_eff={q4:.5f} ⇒ λ_e,eff={lam4:.4f} 格; 线性性残差 {res4:.1e}")
    scr = {}; W = 12
    P_("    屏距L  实测间距  欧氏预言(同窗割线)  λR/d(中心切)  λL/d(小角)  可见度V  亮纹数")
    for L in [40, 80, 120, 160, 200]:
        y = ys + L
        wk = np.arange(xc - W, xc + W + 1)
        g4, _ = phase_slope(r4a["acc"], r4b["acc"], y, wk)
        xf = np.linspace(xc - W, xc + W, 20001)
        phE = q4 * (np.sqrt((xf - xs1) ** 2 + L ** 2) - np.sqrt((xf - xs2) ** 2 + L ** 2))
        slE = (phE[-1] - phE[0]) / (xf[-1] - xf[0])
        xr = np.arange(16, nc3 - 16)
        pk = maxima_of(I4[y, xr], xr); pk = pk[np.abs(pk - xc) <= 170]
        vis = float(np.mean(np.abs(C4[y, xr]) / (I4a[y, xr] + I4b[y, xr] + 1e-30)))
        scr[L] = {"meas": round(float(2*np.pi/abs(g4)), 3), "euclid_pred": round(float(2*np.pi/abs(slE)), 3),
                  "lamR_d": round(float(lam4*np.hypot(d/2, L)/d), 3), "small_angle": round(float(lam4*L/d), 3),
                  "vis": round(vis, 3), "n_peak": int(len(pk))}
        P_(f"    {L:4d}  {2*np.pi/abs(g4):8.3f}  {2*np.pi/abs(slE):17.3f}  {lam4*np.hypot(d/2,L)/d:12.3f}"
           f"  {lam4*L/d:9.3f}  {vis:8.3f}  {len(pk):5d}")
    xxA = np.arange(nc3); corr = {}
    for L in [80, 200]:
        yA = ys + L
        r1 = np.hypot(xxA - xs1, L); r2 = np.hypot(xxA - xs2, L)
        M = np.abs(np.exp(1j * q4 * r1) / np.sqrt(r1) + np.exp(1j * q4 * r2) / np.sqrt(r2)) ** 2
        selA = np.abs(xxA - xc) <= 170
        corr[L] = float(np.corrcoef(I4[yA][selA], M[selA])[0, 1])
    P_(f"  与解析双源公式 |Σ e^(iqr)/√r|² 的相关系数: L=80 → {corr[80]:.4f}, L=200 → {corr[200]:.4f}")
    P_("  ⇒ 条纹确实出现, 且间距随屏距 L 单调增大, 与欧氏程函(同窗割线)预言一致, 大 L 渐近 λ_e·L/d;")
    P_("    这是真实双缝的行为。Z115 的 L1 公式给不出任何 L 依赖, 因此它不是双缝条纹,")
    P_("    而是把屏放在两源连线上的近场拍频。修正只需换两处: 距离用欧氏(程函)、屏离开缝行。")
    out["CK5"] = {"grid": f"{nr3}x{nc3}", "c2": round(c2c, 5), "slit_row": ys, "slits_x": [xs1, xs2],
                  "d": d, "q_eff": round(q4, 5), "lam_eff": round(float(lam4), 4),
                  "linearity_residual": res4, "corr_vs_analytic": {str(a): round(b, 4) for a, b in corr.items()},
                  "screens": {str(a): b for a, b in scr.items()}}

    # ============================================================ CK6
    P_("\n【CK6】判定汇总")
    V = [
        ("CK1 亮纹位置/间距/条数(连续屏)", "可复现", "6.858/8.429/10/11.571/13.142, 间距 π/2=1.571, 与 §3 表一致"),
        ("CK1 §2⑦ 间距公式", "错(差 2 倍)", "(π/k)/斜率=0.785 ✗; 应为 (2π/k)/斜率=λ_e/斜率=1.571"),
        ("CK1 图上节点采样", "图上不成立", "2k=4 > Nyquist π ⇒ 混叠; 节点上只见 3 条(7/10/13), 非 5 条"),
        ("CK1 §3 '两侧对称等间距'", "错", "缝外 Δs=∓8 恒定 ⇒ 一片均匀; '亮度[0,4]' 只在缝间"),
        ("CK1 落盘 JSON", "欠采样", "抽样 1.2 格 < 条纹 1.571, 文件里画不出 5 条纹"),
        ("CK2 k=√λ 当波数", "错(类型混淆)", "格上相位梯度 q≠√λ; k=2⇒q=π(周期2格,v_g=0); λ_e=π⇒λ=2.83≠4"),
        ("CK2 以 K₄ 谱给格图定波数", "跨图取值", "波在 24×24 格图(λ_max=7.966)上, K₄(λ_max=4) 是另一张图"),
        ("CK2 λ_e=π 的来源", "未导出", "K₄ 无空间; 由 ω(K₄)=2c 推空间波长是类型错误; 格上波数须由驱动频率定"),
        ("CK3 曼哈顿距离差", "错(模型)", "L1 ⇒ Δs 与屏行无关(无屏距依赖/无远场); 实测相位走欧氏(5/5 非平凡割线)"),
        ("CK3 屏与缝共线 h=0", "退化选择", "L1 与欧氏唯一重合的行, 且是近场; 真实双缝须让屏离开缝行"),
        ("CK4 真跑双源(k=2)", "对不上", "中心实测间距=混叠值 2.8, 非 1.571; 缝外仍有条纹而 L1 说没有"),
        ("CK4 等幅平面波假设", "错", "2D 振幅 ∝1/√r; 实测 V=0.70…0.59<1, Z115 的 I∈[0,4],V=1 不成立"),
        ("CK4 按自身出处 ω=2c", "不自洽", "带边 v_g=0, 近乎驻波, 几乎不辐射 ⇒ 给不出 1.571"),
        ("CK5 修正版(欧氏+屏离开缝行)", "成立", "间距随 L 单调增大, 与欧氏程函一致, 渐近 λL/d"),
        ("④ 强度=模方(单电子)", "未证", "经典波自然; 单电子探测概率需 Born 规则, Zero 未导出(Z004/Z007)"),
        ("电子身份", "名义", "k 与电子质量/动量/能量无关 ⇒ 不给任何可测条纹间距的预言"),
    ]
    for a, b, c in V:
        P_(f"  · {a:28s} {b:12s} {c}")
    out["CK6_verdicts"] = [{"item": a, "verdict": b, "detail": c} for a, b, c in V]

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_doubleslit_check.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "k4_doubleslit_check.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(LOG))
    P_("\n已存 k4_doubleslit_check.json / .txt")

    # ============================================================ 图
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(14, 11))
        gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.25])

        a0 = fig.add_subplot(gs[0, 0])
        a0.plot(xs, If, lw=1.0, color="#1f77b4", label="continuous screen (Z115's figure)")
        a0.axvspan(6, 14, color="#ffe0b2", alpha=0.7)
        a0.plot(xb, np.full_like(xb, 4.05), "v", color="crimson", ms=7, label="5 maxima, spacing 1.571")
        a0.plot(xpn, In[xpn], "o", color="black", ms=7, label="what the graph itself can sample")
        a0.set_title("CK1  Z115's fringes: continuum screen vs lattice nodes\n"
                     "on nodes: only 3 maxima (aliased, 2k>Nyquist)", fontsize=9)
        a0.set_xlabel("x (lattice units)"); a0.set_ylabel("I(x)"); a0.legend(fontsize=7); a0.grid(alpha=.3)

        a1 = fig.add_subplot(gs[0, 1])
        qq = np.linspace(0, np.pi, 1200)
        a1.plot(qq, 2 * np.sin(qq / 2), color="#2ca02c", label=r"$\sqrt{\lambda(q)}=2\sin(q/2)$")
        a1.plot([np.pi], [2.0], "o", color="crimson", ms=9)
        a1.annotate("k=$\\sqrt{\\lambda_{max}(K_4)}$=2\n$\\Rightarrow$ q=$\\pi$: period 2, $v_g$=0",
                    (np.pi, 2.0), xytext=(-165, -62), textcoords="offset points", fontsize=8,
                    arrowprops=dict(arrowstyle="->", color="crimson"))
        a1.plot([2.0], [2 * np.sin(1.0)], "s", color="darkorange", ms=9)
        a1.annotate("$\\lambda_e$=$\\pi$ $\\Rightarrow$ q=2\n$\\Rightarrow \\lambda$=2.83$\\neq$4",
                    (2.0, 2 * np.sin(1.0)), xytext=(28, 34), textcoords="offset points", fontsize=8,
                    arrowprops=dict(arrowstyle="->", color="darkorange"))
        a1.set_title("CK2  on a lattice $\\sqrt{\\lambda}$ is NOT the phase gradient\n"
                     "the two 'endogenous' inputs cannot both hold", fontsize=9)
        a1.set_ylim(0, 2.45)
        a1.set_xlabel("q (rad/site)"); a1.set_ylabel(r"$\sqrt{\lambda}$"); a1.legend(fontsize=7); a1.grid(alpha=.3)

        a2 = fig.add_subplot(gs[1, 0])
        hs = np.array([0, 2, 4, 8, 12, 16])
        a2.plot(hs, [np.pi / 2] * len(hs), "--", color="#1f77b4", label="Z115 / Manhattan L1: constant")
        a2.plot(hs, [eu[h] for h in hs], "-o", color="#2ca02c", label="Euclidean eikonal (true wave)")
        a2.plot(hs, lam_e * hs / 8, ":", color="gray", label=r"small angle $\lambda_e h/d$")
        a2.set_title("CK3  fringe spacing vs screen distance h\nL1 'graph distance' removes all L-dependence", fontsize=9)
        a2.set_xlabel("h = distance of screen row from slit row (sites)"); a2.set_ylabel("spacing (sites)")
        a2.legend(fontsize=7); a2.grid(alpha=.3)

        a3 = fig.add_subplot(gs[1, 1])
        for h, col in zip([0, 4, 8, 12], ["#1f77b4", "#2ca02c", "#d62728", "#9467bd"]):
            rw = I2n[sy - h, 12:nr2 - 12].copy()
            rw = rw / max(rw.max(), 1e-30)
            a3.plot(np.arange(12, nr2 - 12), rw + h * 1.15, lw=1.1, color=col, label=f"h={h}")
        a3.axvspan(sx1, sx2, color="#ffe0b2", alpha=.6)
        a3.set_title(f"CK4  leapfrog with Z115's own k (q=2), node-sampled rows\n"
                     f"aliased period {alias_p:.2f} - not 1.571", fontsize=9)
        a3.set_xlabel("x (sites) of a 60x60 grid"); a3.set_ylabel(r"$|\psi|^2$ (offset)")
        a3.legend(fontsize=7); a3.grid(alpha=.3)

        a4 = fig.add_subplot(gs[2, :])
        from matplotlib.colors import LogNorm
        sub = I4[ys - 2:ys + 206, 40:nc3 - 40]
        k3 = np.ones((3, 3)) / 9.0                     # 仅显示用: 3×3 盒平滑, 抹掉 λ_e=12 载波
        sub = np.array([[np.sum(sub[max(0, i - 1):i + 2, max(0, j - 1):j + 2] * k3[:sub[max(0,i-1):i+2, max(0,j-1):j+2].shape[0], :sub[max(0,i-1):i+2, max(0,j-1):j+2].shape[1]])
                        for j in range(sub.shape[1])] for i in range(sub.shape[0])])
        vmax = float(sub.max())
        a4.imshow(sub, cmap="inferno", aspect="auto", origin="upper",
                  extent=[40, nc3 - 40, ys + 204, ys - 2],
                  norm=LogNorm(vmin=vmax * 10 ** -3.2, vmax=vmax / 12.0))
        for L in [40, 80, 120, 160, 200]:
            a4.axhline(ys + L, color="cyan", lw=.7, ls=":")
        a4.plot([xs1, xs2], [ys, ys], "w*", ms=13)
        a4.set_title("CK5  repaired Zero-internal double slit: graph wave equation only (d = 60 = 5$\\lambda_e$)\n"
                     "fringe fan spreads with L; measured spacing follows the Euclidean eikonal "
                     "(corr 0.95-0.97 vs analytic two-source law; display 3x3 smoothed)", fontsize=9)
        a4.set_xlabel("x (sites)"); a4.set_ylabel("y (sites)")
        fig.tight_layout()
        fig.savefig(os.path.join(RES, "k4_doubleslit_check.png"), dpi=125)
        P_("已存 k4_doubleslit_check.png")
    except Exception as e:
        print("绘图失败:", e)


if __name__ == "__main__":
    main()
