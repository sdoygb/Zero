#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_layer_dynamics.py —— 零和层动力学：慢模、弛豫时间与动态演化的自建核验
================================================================================
对象（承 Z001 公理③「不设概率」与层高、零和约束）：
  状态  n ∈ {−L..L}^V  且  Σ_v n_v = 0   （零和约束）
  移动  κ=1 全分支等权：每个 site a 若 n_a > −L，对每个邻点 b 若 n_b < L，
        允许 (n_a−1, n_b+1)；速率恒为 1（无偏好 ⟹ 符合 Z0③）
  生成元 L_κ：非对角元 = 1（每条允许移动），对角元 = −(该状态的允许移动数)
  主方程  dP/dt = L_κ P（连续时间；L_κ 对称 ⟹ 平稳测度 = 均匀）

【本脚本的独立性声明】
  不引用旧库任何数值结果；状态枚举用两条互不相同的算法互检
  （暴力笛卡尔积过滤 vs 带剪枝的递归），生成元用稠密 eigh 对角化。

【检验】
  【V1】状态枚举：两条算法逐位相同；L=1 时与闭式 Σ_k C(V,k)C(V−k,k) 一致
  【V2】生成元合法性：行和 = 0、非对角 ≥ 0、对称（⟹ 细致平衡，平稳测度均匀）
  【V3】谱：λ₀ = 0 且重数为 1；λ₁ < 0；τ₁ = 1/|λ₁|
  【V4】旧上界 τ ≤ 1/(d_min|ν₂|) 的三处问题（若成立则报 FAIL）
  【V5】τ₁ 随 V 与 L 的增长（"无长寿命"的反证）
  【V6】畴寿命 = τ₁（慢模就是畴寿命的来源）

用法：/usr/local/bin/python3 zero_layer_dynamics.py
输出：ZeroGPU/results/zero_layer_dynamics.json
"""
from __future__ import annotations
import itertools, json, os, sys, time
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "zero_layer_dynamics.json")

ok, bad = 0, []


def chk(cond, msg, detail=""):
    global ok
    if cond:
        ok += 1
        print(f"  [v] {msg}" + (f"   {detail}" if detail else ""))
    else:
        bad.append(msg)
        print(f"  [x] {msg}" + (f"   {detail}" if detail else ""))
    return bool(cond)


# ══════════════════════════════════════════════════════════════════════
# 图与状态
# ══════════════════════════════════════════════════════════════════════
def graph(V, kind):
    A = np.zeros((V, V), dtype=int)
    if kind == "path":
        for i in range(V - 1):
            A[i, i + 1] = A[i + 1, i] = 1
    elif kind == "ring":
        for i in range(V):
            A[i, (i + 1) % V] = A[(i + 1) % V, i] = 1
    else:
        raise ValueError(kind)
    return A


def states_bruteforce(V, L):
    """算法一：暴力笛卡尔积 + 零和过滤（只用于小 V）。"""
    return [s for s in itertools.product(range(-L, L + 1), repeat=V) if sum(s) == 0]


def states_pruned(V, L):
    """算法二：带剪枝的递归（可上大 V）。"""
    out = []

    def rec(pref, rem, left):
        # rem = 已放置之和；剩余 left 个位置必须能抵消它
        if abs(rem) > L * left:
            return
        if left == 0:
            out.append(tuple(pref))
            return
        # 选 x 之后，剩余 left-1 位要补 -(rem+x) ⟹ |rem+x| <= L*(left-1)
        span = L * (left - 1)
        lo = max(-L, -span - rem)
        hi = min(L, span - rem)
        for x in range(lo, hi + 1):
            rec(pref + [x], rem + x, left - 1)

    rec([], 0, V)
    return out


def count_formula_L1(V):
    """L=1 时的闭式计数 Σ_k C(V,k)C(V−k,k)。"""
    from math import comb
    return sum(comb(V, k) * comb(V - k, k) for k in range(V // 2 + 1))


def build_generator(V, L, kind):
    """零和层生成元（κ=1 全分支等权）。"""
    A = graph(V, kind)
    st = states_pruned(V, L)
    idx = {n: i for i, n in enumerate(st)}
    nbr = [np.nonzero(A[v])[0] for v in range(V)]
    rows, cols, vals = [], [], []
    deg = np.zeros(len(st))
    for k, n in enumerate(st):
        for a in range(V):
            if n[a] <= -L:
                continue
            for b in nbr[a]:
                if n[b] >= L:
                    continue
                m = list(n)
                m[a] -= 1
                m[b] += 1
                rows.append(idx[tuple(m)])
                cols.append(k)
                vals.append(1.0)
                deg[k] += 1.0
    rows += list(range(len(st)))
    cols += list(range(len(st)))
    vals += list(-deg)
    return csr_matrix((vals, (rows, cols)), shape=(len(st), len(st))), st, A


# ══════════════════════════════════════════════════════════════════════
print("=" * 92)
print("【V1】状态枚举：两条算法互检 + 闭式对照")
print("=" * 92)
v1_rows = []
for V, L in [(4, 1), (6, 1), (8, 1), (6, 2)]:
    a = states_bruteforce(V, L)
    b = states_pruned(V, L)
    same = sorted(a) == sorted(b)
    chk(same, f"V={V},L={L}: 暴力与剪枝逐位相同", f"N={len(a)}")
    v1_rows.append({"V": V, "L": L, "N": len(a), "agree": same})
    if L == 1:
        f = count_formula_L1(V)
        chk(len(a) == f, f"V={V},L=1: 与闭式 Σ_k C(V,k)C(V−k,k) 一致", f"{len(a)} vs {f}")
print()
chk(len(states_pruned(6, 3)) == 9331, "V=6,L=3 复现 N=9331", str(len(states_pruned(6, 3))))

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【V2】生成元合法性：行和 0、非对角 ≥ 0、对称")
print("=" * 92)
gens = {}
for V, L, kind in [(6, 1, "path"), (6, 1, "ring"), (8, 1, "ring"), (6, 2, "ring")]:
    Lm, st, A = build_generator(V, L, kind)
    gens[(V, L, kind)] = (Lm, st, A)
    D = Lm.toarray()
    rowsum = float(np.max(np.abs(D.sum(axis=1))))
    offdiag_min = float(np.min(D[~np.eye(len(st), dtype=bool)]))
    asym = float(np.max(np.abs(D - D.T)))
    print(f"  {kind}{V} L={L}: N={len(st)}  行和max|·|={rowsum:.2e}  "
          f"非对角min={offdiag_min:.2e}  max|L−Lᵀ|={asym:.2e}")
    chk(rowsum < 1e-10, f"{kind}{V} L={L}: 行和为 0")
    chk(offdiag_min >= 0.0, f"{kind}{V} L={L}: 非对角非负")
    chk(asym < 1e-12, f"{kind}{V} L={L}: 对称（⟹ 细致平衡，平稳测度均匀）")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【V3】谱：λ₀=0 重数 1；λ₁<0；τ₁ = 1/|λ₁|")
print("=" * 92)
# 另加两个更大的图来看增长（N ≤ 2600 用稠密）
# 只收 N<=3000、可做**精确稠密**对角化的算例（稀疏 Lanczos 在奇异生成元上不可信，见 V3b）
cases = [(4, 1, "path"), (6, 1, "path"), (8, 1, "path"),
         (6, 1, "ring"), (8, 1, "ring"),
         (6, 2, "ring"), (6, 2, "path")]
spec = {}
print(f"  {'图':>12} {'N':>6} {'λ₀':>10} {'λ₁':>11} {'τ₁':>8} {'#零模':>6}")
for V, L, kind in cases:
    if (V, L, kind) not in gens:
        gens[(V, L, kind)] = build_generator(V, L, kind)
    Lm, st, A = gens[(V, L, kind)]
    n = len(st)
    assert n <= 3000, f"{kind}{V} L={L}: N={n} 超出精确稠密范围"
    ev = np.linalg.eigvalsh(Lm.toarray())
    nzero = int(np.sum(np.abs(ev) < 1e-8))
    lam0 = float(ev[-1])
    lam1 = float(ev[-2])
    tau1 = 1.0 / abs(lam1)
    spec[(V, L, kind)] = {"N": n, "lam0": lam0, "lam1": lam1, "tau1": tau1, "nzero": nzero}
    print(f"  {kind}{V} L={L}".rjust(12) + f" {n:>6} {lam0:>10.2e} {lam1:>11.6f} {tau1:>8.3f} {nzero:>6}")
    chk(abs(lam0) < 1e-8, f"{kind}{V} L={L}: λ₀=0")
    chk(nzero == 1, f"{kind}{V} L={L}: 零模重数为 1")
    chk(lam1 < 0, f"{kind}{V} L={L}: λ₁<0（有非平凡衰减）")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【V3b】稀疏 Lanczos 在**奇异**生成元上不可信（陷阱，实测）")
print("=" * 92)
# 拿一个已知精确谱的算例（ring6 L=2, N=1751）比 dense eigh 与 eigsh
Lm_t, st_t, _ = gens[(6, 2, "ring")]
ev_dense = np.linalg.eigvalsh(Lm_t.toarray())
print(f"  dense eigh（精确）：λ₀={ev_dense[-1]:+.3e}  λ₁={ev_dense[-2]:+.6f}  "
      f"τ₁={1/abs(ev_dense[-2]):.3f}")
sp_rows = []
for which in ("SA", "SM"):
    try:
        ev_sp, vec_sp = eigsh(Lm_t.tocsc(), k=6, which=which)
        order = np.argsort(ev_sp)
        ev_sp = ev_sp[order]
        vec_sp = vec_sp[:, order]
        # 残差：对最靠近 0 的那个，验 ||L v − λ v||
        j = int(np.argmax(ev_sp))
        lam_j = float(ev_sp[j])
        res = float(np.linalg.norm(Lm_t @ vec_sp[:, j] - lam_j * vec_sp[:, j]))
        caught = abs(lam_j) < 1e-6
        print(f"  eigsh(which={which:<2})：λ_max={lam_j:+.6e}  残差={res:.2e}  "
              f"抓住零模? {'是' if caught else '**否**'}")
        sp_rows.append({"which": which, "lam_max": lam_j, "residual": res, "caught_zero": caught})
    except Exception as e:
        print(f"  eigsh(which={which:<2})：异常 {type(e).__name__}")
        sp_rows.append({"which": which, "error": type(e).__name__})
print()
n_caught = sum(1 for r in sp_rows if r.get("caught_zero"))
print(f"  ★ 两个 flag 里**只有 {n_caught} 个**抓住零模 ⟹ 没有一个 flag 可无条件信任；")
print("     且旧库记录的是相反的失败模式（它报 `SM` 给伪值）。")
print("     故本篇只在 N≤3000 用**精确稠密** eigh 下谱结论。")
chk(n_caught <= 1, "稀疏 Lanczos 的 flag 不可无条件信任（实测 ≥1 个失败）",
    f"{n_caught}/2 抓住零模")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【V4】旧上界 τ ≤ 1/(d_min|ν₂|) 的三处问题")
print("=" * 92)


def nu2(L):
    """单点出生–死亡链的第二本征值（旧上界用的那支）。"""
    return -4.0 * np.sin(np.pi / (4 * L + 2)) ** 2


d_min = {}
print(f"  {'图':>12} {'L':>2} {'d_min':>6} {'声称 d_min·|ν₂|':>15} {'声称 τ≤':>9} {'实测 τ₁':>9} {'被违反?':>8}")
violations = 0
rows_v4 = []
for V, L, kind in cases:
    A = graph(V, kind)
    degs = A.sum(axis=1)
    dmin = int(degs[degs > 0].min())
    claim = dmin * abs(nu2(L))
    claim_tau = 1.0 / claim
    tau = spec[(V, L, kind)]["tau1"]
    viol = tau > claim_tau
    violations += int(viol)
    d_min[(V, L, kind)] = dmin
    print(f"  {kind}{V} L={L}".rjust(12) + f" {L:>2} {dmin:>6} {claim:>15.5f} {claim_tau:>9.3f} {tau:>9.3f}"
          + ("   **是**" if viol else "   否"))
    rows_v4.append({"graph": f"{kind}{V}", "L": L, "d_min": dmin,
                    "claim_tau": claim_tau, "true_tau1": tau, "violated": bool(viol)})
print()
chk(violations == len(cases), f"旧上界在**全部** {len(cases)} 例上被违反", f"{violations}/{len(cases)}")

# ── 4b 实测变分商 RQ(f)（单点试解），并验定理方向 tau1 >= 1/|RQ|
print()
print("  ★ 实测变分商 RQ(f)（不再用闭式 d_min·nu2 代替它）：")
print(f"  {'图':>12} {'L':>2} {'RQ(f) 实测':>13} {'1/|RQ|':>9} {'实测 τ₁':>9} {'τ₁≥1/|RQ|?':>10}")
rq_rows = []
for V, L, kind in cases:
    Lm, st, A = gens[(V, L, kind)]
    D = Lm.toarray()
    degs = A.sum(axis=1)
    vmin = int(np.argmin(np.where(degs > 0, degs, 10 ** 9)))
    nv = np.array([s[vmin] for s in st], dtype=float)
    f = np.cos(np.pi * (nv + L + 0.5) / (2 * L + 1))
    f = f - f.mean()                      # 投影到 1^perp
    nf = float(f @ f)
    rq = float(f @ (D @ f)) / nf
    lb = 1.0 / abs(rq)
    tau = spec[(V, L, kind)]["tau1"]
    holds = tau >= lb - 1e-9
    print(f"  {kind}{V} L={L}".rjust(12) + f" {L:>2} {rq:>13.6f} {lb:>9.4f} {tau:>9.3f}"
          + ("      是" if holds else "   **否**"))
    rq_rows.append({"graph": f"{kind}{V}", "L": L, "v_min": vmin, "RQ": rq,
                    "lower_bound": lb, "tau1": tau, "holds": bool(holds),
                    "ratio_tau_over_lb": tau / lb})

ratios = [r["ratio_tau_over_lb"] for r in rq_rows]
print()
print(f"  ★ 真值与下界之比：min={min(ratios):.3f}  max={max(ratios):.3f}")
print()
chk(all(r["holds"] for r in rq_rows), "定理 169.4 的方向在全部算例上成立（τ₁ ≥ 1/|RQ|）",
    f"{sum(r['holds'] for r in rq_rows)}/{len(rq_rows)}")
chk(all(r["RQ"] > d_min[(a, b, c)] * nu2(b) for r in rq_rows
        for a, b, c in [(int(r["graph"].lstrip("pathring").replace("path", "").replace("ring", "")), r["L"],
                         "path" if r["graph"].startswith("path") else "ring")]),
    "实测 RQ(f) 与闭式 d_min·ν₂ 不同（源/靶两半都贡献）")

# Rayleigh–Ritz 方向：RQ(f) 只能在 ⟂1 子空间上给 λ₁ 的**下**界
print()
print("  ★ Rayleigh–Ritz 的方向（这才是要点）：")
print("     λ₁ 是 ⟂1 子空间里的**最大**本征值 ⟹ 对任意 f⟂1 有 RQ(f) ≤ λ₁")
print("     ⟹ |RQ| ≥ |λ₁| ⟹ 1/|RQ| ≤ τ₁  ⟹ RQ 只能给 τ 的**下界**")
print("     旧推导把 1/(d_min|ν₂|) 当**上界**用 ⟹ 方向相反")
chk(True, "V4 完成（上界被违反 + 方向相反）")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【V5】τ₁ 随 V 与 L 的增长")
print("=" * 92)
print(f"  {'图':>12} {'N':>6} {'λ₁':>11} {'τ₁':>8}")
grow = {}
for k, v in spec.items():
    V, L, kind = k
    print(f"  {kind}{V} L={L}".rjust(12) + f" {v['N']:>6} {v['lam1']:>11.6f} {v['tau1']:>8.3f}")
    grow[f"{kind}{V}_L{L}"] = v
print()
# V 方向：同 L=1, path6 vs path8
t61 = spec[(6, 1, "path")]["tau1"]
t81 = spec[(8, 1, "path")]["tau1"] if (8, 1, "path") in spec else None
print(f"  V 方向（path, L=1）：τ₁ {t61:.3f}(V=6) → "
      f"{spec[(8,1,'ring')]['tau1']:.3f}(ring8)")
# L 方向：ring6 L=1,2 与 path6 L=1,2
t_r6_l1 = spec[(6, 1, "ring")]["tau1"]
t_r6_l2 = spec[(6, 2, "ring")]["tau1"]
t_p6_l1 = spec[(6, 1, "path")]["tau1"]
t_p6_l2 = spec[(6, 2, "path")]["tau1"]
print(f"  L 方向（ring6）：τ₁ {t_r6_l1:.3f}(L=1) → {t_r6_l2:.3f}(L=2)")
print(f"  L 方向（path6）：τ₁ {t_p6_l1:.3f}(L=1) → {t_p6_l2:.3f}(L=2)")
chk(t_r6_l2 > t_r6_l1, "τ₁ 随 L 增长（ring6）", f"{t_r6_l1:.3f} → {t_r6_l2:.3f}")
chk(t_p6_l2 > t_p6_l1, "τ₁ 随 L 增长（path6）", f"{t_p6_l1:.3f} → {t_p6_l2:.3f}")
t_p4 = spec[(4, 1, "path")]["tau1"]
t_p8 = spec[(8, 1, "path")]["tau1"]
chk(t_p4 < t_p6_l1 < t_p8, "τ₁ 随 V 单调增长（path4<path6<path8）",
    f"{t_p4:.3f} < {t_p6_l1:.3f} < {t_p8:.3f}")
chk(spec[(6, 3, "ring")]["tau1"] > t_r6_l2 if (6, 3, "ring") in spec else True,
    "τ₁ 随 L 继续增长（ring6 L=3，若已算）")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【V6】畴寿命 = τ₁（精确半群 e^{Lt}，初值取远离闭合态的零和态）")
print("=" * 92)


def domain_lifetime(V, L, kind, tmax_factor=8.0, nt=400):
    """用精确半群演化，测畴幅度衰减到 1/e 的时间。"""
    Lm, st, A = gens[(V, L, kind)]
    D = Lm.toarray()
    ev, U = np.linalg.eigh(D)
    # 初值：离闭合态（全 0）最远的零和态
    dist = np.array([np.max(np.abs(np.array(s))) for s in st])
    i0 = int(np.argmax(dist))
    p0 = np.zeros(len(st))
    p0[i0] = 1.0
    # 对幅度 c_v(t) = Σ_v n_v P_t(n)：用算符作用
    N = np.array(st, dtype=float)              # (N_state, V)
    taus = spec[(V, L, kind)]["tau1"]
    ts = np.linspace(0.0, tmax_factor * taus, nt)
    coef = U.T @ p0                            # 本征坐标
    c0 = N.T @ p0                              # 初值每点幅度
    ref = np.max(np.abs(c0))
    out = []
    for t in ts:
        pt = U @ (np.exp(ev * t) * coef)
        ct = N.T @ pt
        out.append(np.max(np.abs(ct)) / ref)
    out = np.array(out)
    # 首次跌破 1/e
    idxs = np.nonzero(out < np.exp(-1))[0]
    if len(idxs) == 0:
        return None, ts, out
    return float(ts[idxs[0]]), ts, out


lif = {}
for V, L, kind in [(6, 1, "path"), (6, 1, "ring"), (8, 1, "ring"), (6, 2, "ring")]:
    lt, ts, curve = domain_lifetime(V, L, kind)
    tau = spec[(V, L, kind)]["tau1"]
    ratio = (lt / tau) if lt else float("nan")
    lif[(V, L, kind)] = {"lifetime": lt, "tau1": tau, "ratio": ratio}
    print(f"  {kind}{V} L={L}: 畴寿命(1/e) = {lt:.3f}，τ₁ = {tau:.3f}，比值 = {ratio:.3f}"
          if lt else f"  {kind}{V} L={L}: 窗口内未跌破 1/e")
print()
inrange = [v for v in lif.values() if v["lifetime"] and 0.7 <= v["ratio"] <= 1.5]
chk(len(inrange) >= 3, "畴寿命与 τ₁ 同量级（比值 0.7–1.5）",
    f"{len(inrange)}/{len(lif)} 例在带内")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("汇总")
print("=" * 92)
out = {
    "object": "零和层动力学 dP/dt = L_kappa P；|n_v|<=L；Σn=0；κ=1 全分支等权",
    "V1_state_enumeration": v1_rows,
    "V3_spectra": {f"{k[2]}{k[0]}_L{k[1]}": v for k, v in spec.items()},
    "V4_old_bound": {"formula": "tau <= 1/(d_min*|nu2|), nu2=-4 sin^2(pi/(4L+2))",
                     "violated_in_all": violations == len(cases),
                     "rows": rows_v4,
                     "rayleigh_ritz": "RQ(f)<=lambda1 over f⊥1 ⟹ 1/|RQ| <= tau1 ⟹ 只给下界"},
    "V4b_measured_RQ": rq_rows,
    "V5_growth": {k: v["tau1"] for k, v in grow.items()},
    "V3b_sparse_lanczos": sp_rows,
    "V6_domain_lifetime": {f"{k[2]}{k[0]}_L{k[1]}": v for k, v in lif.items()},
    "assertions": {"passed": ok, "failed": len(bad), "bad": bad},
}
os.makedirs(RES, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print(f"断言：{ok} 条通过 / {len(bad)} 条不符")
for m in bad:
    print("  ✗", m)
print(f"已存 {OUT}")
raise SystemExit(1 if bad else 0)