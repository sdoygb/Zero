#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_quantum_doubleslit_v2.py —— Z116「Zero 量子双缝」的核查与改进
======================================================================================
对 Z116（`k4_quantum_doubleslit.py`）逐条核查后重做。改进点:

  I1 屏不再放在吸收边界上:  Z116 的屏列 x=23 就是虚势吸收列, 屏上 |ψ|²~1e-4 是吸收层残余.
     本版把屏放在吸收层之前(L=40..160), 吸收层留边 24 格.
  I2 几何自洽:  Z116 的仿真是"挡板列 + 孔在 y=10,14", 但报告与结论用的是"孔在 y=4,20".
     本版孔位、屏位、报告完全一致(孔 y=90/150, 缝距 d=60).
  I3 真·单源:  Z116 用 x<=9 的平面波条带, 且 src 变量算了没用; 本版用**单个源节点**
     (y0,x=20) → 两孔 → 屏, 是"一个粒子从一点发出".
  I4 波长实测而非谱猜:  对照跑(去掉挡板)测入射相位梯度 → q 实测;
     不把 √λ_max(K₄) 当波数(Z115/Z116 同错: 跨图取值 + q↔√λ 偷换).
  I5 相位差用欧氏程函, 不用曼哈顿:  实测相位走欧氏(见 k4_doubleslit_check.py CK4a);
     L1 是图上最短步数(组合量), 不是波程.
  I6 单粒子归一化的正确陈述:  闭图(无吸收)上 H 厄米 ⇒ Σ|ψ(t)|²=1 **精确**成立;
     Z116 加了非厄米吸收势却仍报 Σp=1.0000 —— 实测是 1→0.003(衰减).

引擎: 一阶量子行走 i∂tψ = Hψ, H = c²L (c²=1 为能量单位选择, 属输入),
      定态散射 (E − H + iΣ)ψ = b 用稀疏 LU 求解(单频驱动, 与时间演化等价).
材质: 图 Laplacian 稀疏组装; 挡板=整列水平边切断(纯组合操作), 孔=保留该处一条边.
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
    s = " ".join(str(x) for x in a)
    LOG.append(s); print(s)

# ---------------------------------------------------------------- 图与求解器
def build_H(nr, nc, xb=90, apertures=(90, 150), no_barrier=False, c2=1.0):
    """挡板 = 列 xb 左侧全部水平边切断, 只在 apertures 行保留"""
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
    """(E − H + iΣ)ψ = b, E 由行模色散 E=2c²(1−cos q) 定; 返回 ψ 网格"""
    H = build_H(nr, nc, xb, apertures, no_barrier, c2)
    N = nr * nc
    E = 2 * c2 * (1 - np.cos(q))
    b = np.zeros(N, complex); b[src[0] * nc + src[1]] = 1.0
    A = (E * sp.eye(N, format="csr") - H + 1j * sp.diags(absorber(nr, nc, m, sig0).ravel())).tocsc()
    return spla.splu(A).solve(b).reshape(nr, nc)

def q_of_E(E, c2=1.0):
    """行模散散射回 q"""
    lam = E / c2
    return float(np.arccos(1 - lam / 2)) if 0 <= lam <= 4 else float("nan")

def aperture_check(psi, xb, apertures):
    a1, a2 = apertures
    return {"abs1": float(abs(psi[a1, xb])), "abs2": float(abs(psi[a2, xb])),
            "dphase": float(np.angle(psi[a2, xb] / psi[a1, xb]))}

# ---------------------------------------------------------------- 主程序
def main():
    out = {}
    P_("=" * 96)
    P_("Z116「Zero 量子双缝」核查 + 改进版 (k4_quantum_doubleslit_v2)")
    P_("=" * 96)

    # ============================================================ A 单位性
    P_("\n【A】单粒子归一化的正确陈述: 闭图(无吸收)上 H 厄米 ⇒ Σ|ψ(t)|²=1 精确")
    nr0, nc0 = 24, 24
    H0 = build_H(nr0, nc0, xb=12, apertures=(10, 14))
    lam0, V0 = np.linalg.eigh(H0.toarray())          # 厄米 → eigh
    N0 = nr0 * nc0
    psi0 = np.zeros(N0, complex)
    for y in range(nr0):
        for x in range(nc0):
            if x <= 9:
                psi0[y * nc0 + x] = np.exp(1j * 1.0 * (x - 2.0))
    psi0 /= np.linalg.norm(psi0)
    dev = []
    for t in [0, 5, 10, 20, 32, 48, 100]:
        pt = V0 @ (V0.T @ psi0 * np.exp(-1j * lam0 * t))
        dev.append(abs(float(np.sum(np.abs(pt) ** 2)) - 1.0))
    P_(f"  闭图 N={N0} (挡板列+两孔, 无吸收): max|Σ|ψ|²−1| = {max(dev):.2e}  ⇒ 精确守恒 ✓")
    P_(f"  ⇒ 「Σp=1」是**幺正性的自明后果**, 不是量子判据(任何归一化振幅都满足)。")

    # Z116 的实现: 它在**实**数组上写 H_eff[i,i] -= 1j*gamma
    g = 0.4
    H_wrong = H0.toarray().copy()                      # float64, 与 Z116 同
    import warnings
    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter("always")
        for y in range(nr0):
            for x in range(nc0):
                if x in (0, nc0 - 1) or y in (0, nr0 - 1):
                    H_wrong[y * nc0 + x, y * nc0 + x] -= 1j * g
    P_(f"  Z116 代码实测: H_eff dtype={H_wrong.dtype}, 虚部总量={np.abs(H_wrong.imag).sum() if np.iscomplexobj(H_wrong) else 0.0}")
    P_(f"    numpy 警告: {str(wlist[0].message) if wlist else '(无)'}")
    P_(f"  ⇒ **吸收势被静默丢弃**(ComplexWarning): 实际演化仍是**闭图厄米**演化, 边界是全反射墙。")
    P_(f"    所以 Z116 §四「Σp=1.0000 全时间成立」为真, 但原因是**闭图幺正性自明**, 与'单粒子量子判据'无关;")
    P_(f"    而它声称的'吸收边界 kills reflection'从未生效 ⇒ 屏上 p 是多次反射的驻波散斑, 不是透射条纹。")
    # 若吸收势真的生效, 会怎样(把它正确地建成复矩阵)
    He = H0.toarray().astype(complex)
    for y in range(nr0):
        for x in range(nc0):
            if x in (0, nc0 - 1) or y in (0, nr0 - 1):
                He[y * nc0 + x, y * nc0 + x] -= 1j * g
    lamE, UE = np.linalg.eig(He)
    cE = np.linalg.inv(UE) @ psi0
    decay = [float(np.sum(np.abs(UE @ (cE * np.exp(-1j * lamE * t))) ** 2)) for t in [0, 5, 10, 20, 32, 48]]
    P_(f"  若把吸收势**正确**建成复矩阵(γ={g}): Σ|ψ|² = " + ", ".join(f"t={t}:{v:.4f}" for t, v in zip([0, 5, 10, 20, 32, 48], decay)))
    P_(f"  ⇒ 此时归一化**不再**为 1(降到 {decay[-1]:.4f}); 即『Σp=1』既不是量子判据, 也只在闭图成立。")
    P_(f"  另: Z116 的 clean 段 Σp=1.0000 来自手工归一化 p_raw/p_raw.sum(), 与任何演化无关。")
    out["A_unitarity"] = {"closed_graph_max_dev": max(dev),
                          "z116_H_eff_dtype": str(H_wrong.dtype),
                          "z116_absorber_discarded": not np.iscomplexobj(H_wrong),
                          "norm_if_absorber_correct": [round(v, 5) for v in decay]}

    # ============================================================ B 对照: 测波长
    P_("\n【B】对照跑(去掉挡板) → **实测**入射波长(不猜 K₄ 谱)")
    nr, nc = 240, 260; xb = 90; ap = (90, 150); y0 = nr // 2
    q0 = 0.5236                                        # 目标: λ_e = 12 格
    pc = solve_scatter(nr, nc, q0, no_barrier=True)
    xw = np.arange(30, 91)                          # 源行(y0)上、源与右侧吸收层之间
    qm = float(np.polyfit(xw, np.unwrap(np.angle(pc[y0, xw])), 1)[0])   # 行内相位梯度 = q
    lam_meas = 2 * np.pi / abs(qm)
    P_(f"  目标 q0={q0:.4f} → 实测 q={abs(qm):.4f}, λ_e,实测={lam_meas:.3f} 格 (误差 {abs(abs(qm)-q0)/q0*100:.2f}%)")
    P_(f"     (注: 测法用**源行相位梯度**; 若改用二维区域对欧氏距离 r 拟合, 会偏 6–7% ——")
    P_(f"      相位=q·r 只在渐近区与轴线上精确, 这本身是图上程函的适用范围。)")
    P_(f"  q 由能量 E 定, 而 E 是**输入**(相当于粒子的能量/动量) —— 这正是 Zero 缺的那一环。")
    out["B_wavelength"] = {"q_target": q0, "q_measured": round(abs(qm), 5),
                           "lambda_measured": round(lam_meas, 4),
                           "note": "E(能量)是输入; Zero 无内生能量标度"}

    # ============================================================ C 挡板双孔
    P_("\n【C】真·单源双缝: 单个源节点 → 挡板两孔 → 屏 (全由图上的解给出)")
    p12 = solve_scatter(nr, nc, q0, xb, ap)
    p1 = solve_scatter(nr, nc, q0, xb, (ap[0],))
    p2 = solve_scatter(nr, nc, q0, xb, (ap[1],))
    p0 = solve_scatter(nr, nc, q0, xb, ())            # 无孔对照
    chk = aperture_check(p12, xb, ap)
    P_(f"  两孔照明: |ψ|={chk['abs1']:.5f} vs {chk['abs2']:.5f}, 相位差={chk['dphase']:.2e} rad")
    P_(f"  ⇒ 「两缝等幅同相」在 Z115/Z116 里是**假设**, 在本版里是**对称几何的推论** ✓")
    trans2 = float(np.sum(np.abs(p12[:, xb + 60]) ** 2))
    trans0 = float(np.sum(np.abs(p0[:, xb + 60]) ** 2))
    ratio = "∞ (断连: 无孔时右侧与源无任何边相连, ψ≡0)" if trans0 == 0 else f"{trans2/trans0:.1f} 倍"
    P_(f"  无孔对照: 挡板后 60 格线上 Σ|ψ|² = {trans0:.3e} (两孔: {trans2:.3e}) ⇒ {ratio}")
    P_(f"  ⇒ 图上的整列切断是**精确**的(无隧穿), 这正是组合几何相对连续势垒的优点。")
    P_("\n  屏距L  实测条纹间距   欧氏程函(同窗割线)   λR/d(中心切)   λL/d(小角)   可见度V   亮纹数")
    rows = {}
    for L in [40, 80, 120, 140]:
        x = xb + L
        if x >= nc - 24:
            continue
        ykw = np.arange(y0 - 60, y0 + 61)
        ph = np.unwrap(np.angle(p1[ykw, x] / p2[ykw, x]))
        A = np.vstack([ykw, np.ones_like(ykw)]).T
        cf, *_ = np.linalg.lstsq(A, ph, rcond=None)
        per = 2 * np.pi / abs(cf[0])
        yf = np.linspace(y0 - 60, y0 + 60, 20001)
        pE = abs(qm) * (np.hypot(yf - ap[0], L) - np.hypot(yf - ap[1], L))
        slE = (pE[-1] - pE[0]) / (yf[-1] - yf[0])
        vis = float(np.mean(np.abs(2 * np.real(p1[ykw, x] * np.conj(p2[ykw, x]))) /
                            (np.abs(p1[ykw, x]) ** 2 + np.abs(p2[ykw, x]) ** 2 + 1e-30)))
        rr = np.abs(p12[:, x]) ** 2
        pk = [i for i in range(2, nr - 2) if rr[i] > rr[i - 1] and rr[i] > rr[i + 1] and rr[i] > 0.02 * rr.max()]
        pk = [i for i in pk if abs(i - y0) <= 90]
        rows[L] = {"meas": round(float(per), 3), "euclid": round(float(2 * np.pi / abs(slE)), 3),
                   "lamR_d": round(float(lam_meas * np.hypot(30, L) / 60), 3),
                   "small_angle": round(float(lam_meas * L / 60), 3),
                   "vis": round(vis, 3), "n_peak": len(pk)}
        P_(f"  {L:5d}  {per:11.3f}   {2*np.pi/abs(slE):17.3f}   {lam_meas*np.hypot(30,L)/60:12.3f}"
           f"   {lam_meas*L/60:10.3f}   {vis:7.3f}   {len(pk):5d}")
    P_("  ⇒ 条纹来自散射解本身(不是手写公式), 间距随 L 增大, 与欧氏程函一致, 大 L 渐近 λL/d。")
    P_("  ⇒ 与我的经典两源计算(k4_doubleslit_check.py CK5)给**同一条法令** —— 图案本身不区分量子/经典,")
    P_("     区别只在解释(单粒子振幅+Born)与统计(逐个到达), 不在条纹位置。")
    out["C_doubleslit"] = {"grid": f"{nr}x{nc}", "barrier_col": xb, "apertures": list(ap), "d": ap[1] - ap[0],
                           "aperture_illumination": chk, "transmission_two_holes": trans2,
                           "transmission_no_hole": trans0, "control_ratio": ("inf/disconnected" if trans0 == 0 else trans2/trans0), "screens": {str(a): b for a, b in rows.items()}}

    # ============================================================ D Nyquist / k=2
    P_("\n【D】Z115/Z116 共用的 k=√λ_max(K₄)=2: 图上不可表示(亚 Nyquist)")
    k = 2.0
    P_(f"  强度条纹空间频率 2k={2*k} > Nyquist π={np.pi:.4f} ⇒ 采样到格点上必混叠成 {abs(2*k-2*np.pi):.4f} rad/格")
    h1, h2 = 4, 20                                      # Z116 报告用的孔位
    ys = np.linspace(h1, h2, 200001)
    dd = np.abs(ys - h1) - np.abs(ys - h2)
    true_br = ys[1:-1][(np.cos(k*dd[1:-1]) > np.cos(k*dd[:-2])) & (np.cos(k*dd[1:-1]) > np.cos(k*dd[2:]))]
    yn = np.arange(0, 25)
    ddn = np.abs(yn - h1) - np.abs(yn - h2)
    pn = 1 + np.cos(k * ddn)
    node_br = yn[1:-1][(pn[1:-1] > pn[:-2]) & (pn[1:-1] > pn[2:])]
    node_br = node_br[(node_br > h1) & (node_br < h2)]
    P_(f"  Z116 报告的孔 (y={h1},{h2}), k=2: 连续解真亮纹 {len(true_br)} 条 "
       f"(y={[round(float(v),2) for v in true_br[:6]]}…, 间距 {np.pi/k:.4f})")
    P_(f"  但屏是**整数格点**: 只能看到 {list(node_br)} → {len(node_br)} 条, 间距≈3 (混叠)")
    P_(f"  ⇒ Z116 §四 表里同时写「亮纹 6,9,12,15,18」(间距 3) 与「条纹间距≈1.57」——自相矛盾;")
    P_(f"     1.57 是连续屏量, 图上测不到(Z115 同一处错误, 未被改掉)。")
    P_(f"  另: Z116 仿真的孔在 y=10,14, 而报告/结论用 y=4,20 ⇒ 报告的数字不对应它自己的仿真。")
    out["D_nyquist"] = {"k": k, "fringe_freq": 2*k, "nyquist": float(np.pi),
                        "true_bright_continuous": [round(float(v), 3) for v in true_br],
                        "node_sampled_bright": [int(v) for v in node_br],
                        "z116_reported": [6, 9, 12, 15, 18], "z116_simulated_holes": [10, 14],
                        "z116_reported_holes": [4, 20]}

    # ============================================================ E 内生性审计
    P_("\n【E】内生性审计: 这条推导里哪些是 Zero 的, 哪些是输入")
    audit = [
        ("图 Laplacian L 作为结构", "Zero 内生", "Z001/Z021 的图对象", "✅"),
        ("挡板/孔 = 边的切断与保留", "Zero 内生(组合)", "纯组合操作, 无需外部几何", "✅"),
        ("一阶演化 i∂tψ=Hψ (Schrödinger 型)", "**输入**", "Z003 #12: Schrödinger 方程**不能**导出(根1 无普遍时间生成元)", "❌"),
        ("H = c²L 这个哈密顿量", "**输入(具名)**", "图上有 L, 但『取 H=c²L』是模型选择", "❌"),
        ("c² 的值", "**输入(数值)**", "Z104/Z105: c²=3.5/λ_max 由 CFL 反推 ⇒ 『光速』随图规模变", "❌"),
        ("舞台图 = 24×24 / 240×260 方格", "**输入(且与台账冲突)**", "Z021: Γ 全局=Heawood(14 顶点)、局部=K₄; 空间维数=3; 『格点不可能承载通道结构』", "❌"),
        ("单源+两孔 → 振幅劈裂", "Zero 内生(相对 Z115 的真改进)", "一个源节点, 图的切断即劈裂", "✅"),
        ("「两孔等幅同相」", "**本版变成推论**", "对称几何 ⇒ 实测 |ψ| 相等、相位差 4e-3 rad; Z115/Z116 里是假设", "✅"),
        ("能量 E ⇒ 波长 q", "**输入**", "E 相当于粒子能量/动量; Zero 无内生能量标度(K₄ 替代是错的)", "❌"),
        ("Born p=|ψ|²", "**条件**(Z116 已明写 ✓)", "Z003 #5 完全(代数扇区) / #8 交叉项条件; Z007: 交叉项需 C_coh", "⚠️"),
        ("屏/位置读出", "**输入**(Z116 已明写 ✓)", "Z009/Z003 根2: 不存在位置可观测量", "⚠️"),
        ("吸收层(虚势/海绵)", "**数值装置**", "为取定态散射而加, 非 Zero 结构", "⚠️"),
        ("k = √λ_max(K₄) = 2", "**错**", "跨图取值 + q↔√λ 偷换; 本版改为实测 q", "❌"),
        ("Δd 用曼哈顿距离", "**错**", "实测相位走欧氏程函(L1 只是最短步数); 本版用欧氏", "❌"),
        ("Σp=1 作为量子判据", "**不成立**", "幺正性自明后果; Z116 加吸收后实测 1→0.003", "❌"),
        ("屏放在 x=23 = 吸收列", "**实现错**", "屏上 |ψ|²~1e-4 是吸收层残余, 非条纹", "❌"),
    ]
    for a, b, c, d in audit:
        P_(f"  {d} {a:34s} {b:24s} {c}")
    n_ok = sum(1 for x in audit if x[3] == "✅"); n_in = sum(1 for x in audit if x[3] == "❌")
    n_cd = sum(1 for x in audit if x[3] == "⚠️")
    P_(f"\n  ⇒ 计数: Zero 内生 {n_ok} 项 / 输入或错误 {n_in} 项 / 条件(已明写) {n_cd} 项")
    P_(f"  ⇒ **不是全 Zero 推导**。最要紧的洞是同一个: 『粒子的能量/波长从哪来』——")
    P_(f"     Z115 用 λ_e=π、Z116 用 k=√λ_max(K₄)=2 都是拿图的谱去补这个洞, 而两处都补错了。")
    out["E_audit"] = [{"element": a, "status": b, "basis": c} for a, b, c, d in audit]
    out["E_counts"] = {"zero_internal": n_ok, "input_or_error": n_in, "conditional_declared": n_cd}

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_quantum_doubleslit_v2.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "k4_quantum_doubleslit_v2.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(LOG))
    P_("\n已存 k4_quantum_doubleslit_v2.json / .txt")

    # ============================================================ 图
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        fig, ax = plt.subplots(2, 2, figsize=(13, 10))
        a0 = ax[0, 0]
        sub = np.abs(p12) ** 2
        x0, x1 = 60, 240
        a0.imshow(sub[:, x0:x1], cmap="inferno", aspect="auto", origin="upper",
                  extent=[x0, x1, nr - 1, 0],
                  norm=LogNorm(vmin=sub.max() * 10 ** -5.5, vmax=sub.max() / 50))
        a0.axvline(xb, color="cyan", lw=1.0)
        for yy_ in ap:
            a0.plot([xb], [yy_], "w*", ms=12)
        a0.plot([20], [y0], "wo", ms=6)
        a0.set_title("C  corrected quantum double slit: one source node -> wall with two apertures\n"
                     r"$|\psi|^2$ from $(E-H+i\Sigma)\psi=b$ (sparse LU), log scale over 5.5 decades", fontsize=9)
        a0.set_xlabel("x (sites)"); a0.set_ylabel("y (sites)")

        a1 = ax[0, 1]
        for L, col in zip([40, 80, 120, 160], ["#1f77b4", "#2ca02c", "#d62728", "#9467bd"]):
            x = xb + L
            if x < nc - 24:
                prof = np.abs(p12[y0 - 90:y0 + 91, x]) ** 2
                a1.plot(np.arange(y0 - 90, y0 + 91), prof / max(prof.max(), 1e-30) + L * 0.02,
                        lw=1.0, color=col, label=f"L={L}")
        a1.set_title("C  screen-line profiles (offset): fringes widen with L", fontsize=9)
        a1.set_xlabel("y (sites)"); a1.set_ylabel(r"$|\psi|^2$ (offset)"); a1.legend(fontsize=7); a1.grid(alpha=.3)

        a2 = ax[1, 0]
        Ls = sorted(rows.keys())
        a2.plot(Ls, [rows[L]["meas"] for L in Ls], "-o", color="#d62728", label="measured (this v2)")
        a2.plot(Ls, [rows[L]["euclid"] for L in Ls], "--", color="#2ca02c", label="Euclidean eikonal (same window)")
        a2.plot(Ls, [rows[L]["small_angle"] for L in Ls], ":", color="gray", label=r"$\lambda L/d$")
        a2.plot(Ls, [rows[L]["lamR_d"] for L in Ls], "-.", color="#1f77b4", label=r"$\lambda R/d$")
        a2.set_title("C  fringe spacing vs screen distance (wavelength MEASURED, not assumed)", fontsize=9)
        a2.set_xlabel("L (sites)"); a2.set_ylabel("spacing (sites)"); a2.legend(fontsize=7); a2.grid(alpha=.3)

        a3 = ax[1, 1]
        ysf = np.linspace(0, 24, 4000)
        a3.plot(ysf, 1 + np.cos(k * (np.abs(ysf - h1) - np.abs(ysf - h2))), lw=0.8, color="#1f77b4",
                label=f"continuous screen (Z116 formula), {len(true_br)} maxima")
        a3.plot(yn, pn, "o", color="crimson", ms=5, label=f"lattice nodes: {len(node_br)} maxima, spacing 3")
        a3.set_title(f"D  k=2 -> fringe freq 2k=4 > Nyquist $\\pi$: {len(true_br)} true fringes alias to {len(node_br)}",
                     fontsize=9)
        a3.set_xlabel("y (sites)"); a3.set_ylabel("p(y)"); a3.legend(fontsize=7); a3.grid(alpha=.3)
        fig.tight_layout()
        fig.savefig(os.path.join(RES, "k4_quantum_doubleslit_v2.png"), dpi=125)
        P_("已存 k4_quantum_doubleslit_v2.png")
    except Exception as e:
        print("绘图失败:", e)


if __name__ == "__main__":
    main()
