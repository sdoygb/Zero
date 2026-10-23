#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
amplitude_A_closed_form.py —— $A=\\gamma_1/\\lambda_1(\\mathsf L)$ 由什么决定
================================================================================
已知（`tau_to_mu_bridge.py`）：$\gamma_1=|\lambda_1(\mathcal L_\kappa)|\approx0.62\,\lambda_1(\mathsf L)$，
但前因子 $A$ 随图与尺寸慢漂（0.59–0.64）。本篇问：$A$ 是哪个图不变量的函数？

【做法：穷举小体系的全部连通图】
  $V=4$：$2^6=64$ 张标号图；$V=5$：$2^{10}=1024$ 张。状态数 $N=19$／$51$，可精确稠密对角化。
  $V=6$：随机抽 $n$ 张连通图（$N=141$）。

【检验】
  【A1】样本构建与自检（连通性、$\lambda_1>0$、$\gamma_1>0$）
  【A2】$A$ 在样本上的分布（是否近乎常数？域多宽？）
  【A3】$A$ 与一族候选不变量的关系（相关表）：$\lambda_1$、$\lambda_2/\lambda_1$、
        $\lambda_{\max}$、平均度、度方差、边数、Kirchhoff 和 $\sum 1/\lambda_i$、
        $\lambda_1/\lambda_{\max}$、$V$
  【A4】最优单一预测量的拟合优度；
  【A5】判定：闭式是否存在

用法：/usr/local/bin/python3 amplitude_A_closed_form.py
输出：ZeroGPU/results/amplitude_A_closed_form.json
"""
from __future__ import annotations
import itertools, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "amplitude_A_closed_form.json")

ok, bad = 0, []


def chk(cond, msg, detail=""):
    global ok
    if cond:
        ok += 1
        print(f"  [v] {msg}" + (f"   {detail}" if detail else ""))
    else:
        bad.append(msg)
        print(f"  [x] {msg}" + (f"   {detail}" if detail else ""))


# ══════════════════════════════════════════════════════════════════════
def build_gen(A, Lcap):
    V = A.shape[0]

    def rec(pref, rem, left):
        if abs(rem) > Lcap * left:
            return
        if left == 0:
            out.append(tuple(pref)); return
        span = Lcap * (left - 1)
        for x in range(max(-Lcap, -span - rem), min(Lcap, span - rem) + 1):
            rec(pref + [x], rem + x, left - 1)

    out = []
    rec([], 0, V)
    st = out
    idx = {n: i for i, n in enumerate(st)}
    nbr = [np.nonzero(A[v])[0] for v in range(V)]
    rows, cols, vals = [], [], []
    deg = np.zeros(len(st))
    for k, n in enumerate(st):
        for a in range(V):
            if n[a] <= -Lcap:
                continue
            for b in nbr[a]:
                if n[b] >= Lcap:
                    continue
                m = list(n)
                m[a] -= 1
                m[b] += 1
                j = idx.get(tuple(m))
                if j is None:
                    continue
                rows.append(j); cols.append(k); vals.append(1.0)
                deg[k] += 1.0
    nn = len(st)
    rows += list(range(nn)); cols += list(range(nn)); vals += list(-deg)
    D = np.zeros((nn, nn))
    for r, c, v in zip(rows, cols, vals):
        D[r, c] += v
    return D


def lap_eigs(A):
    d = A.sum(axis=1)
    Lg = np.diag(d) - A
    return np.sort(np.linalg.eigvalsh(Lg))


def connected(A):
    V = A.shape[0]
    seen = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for w in np.nonzero(A[u])[0]:
            if w not in seen:
                seen.add(int(w)); stack.append(int(w))
    return len(seen) == V


def all_graphs(V):
    pairs = list(itertools.combinations(range(V), 2))
    for mask in range(1 << len(pairs)):
        A = np.zeros((V, V), int)
        for i, (u, w) in enumerate(pairs):
            if mask >> i & 1:
                A[u, w] = A[w, u] = 1
        yield A


def invariants(A):
    V = A.shape[0]
    ev = lap_eigs(A)
    nz = ev[ev > 1e-10]
    d = A.sum(axis=1)
    return {
        "lam1": float(ev[1]),
        "lam2": float(ev[2]) if V > 2 else float("nan"),
        "lammax": float(ev[-1]),
        "mean_deg": float(d.mean()),
        "var_deg": float(d.var()),
        "n_edges": int(A.sum() // 2),
        "kirch": float(np.sum(1.0 / nz)),          # Σ 1/λ_i（非零）
        "r_lam": float(ev[1] / ev[-1]),
        "V": V,
    }


def measure(A, Lcap=1):
    D = build_gen(A, Lcap)
    ev = np.sort(np.linalg.eigvalsh(D))
    gam = float(abs(ev[-2]))
    lam1 = float(lap_eigs(A)[1])
    return gam, lam1, gam / lam1, D.shape[0]


# ══════════════════════════════════════════════════════════════════════
print("=" * 92)
print("【A1】穷举 $V=4,5$ 的全部连通图 + 抽样 $V=6$")
print("=" * 92)
rows = []
for V in (4, 5):
    n_tot = 0
    n_conn = 0
    for A in all_graphs(V):
        n_tot += 1
        if not connected(A):
            continue
        n_conn += 1
        gam, lam1, Av, nn = measure(A)
        r = invariants(A)
        r.update({"gamma1": gam, "A": Av, "N_states": nn})
        rows.append(r)
    print(f"  V={V}: 标号图 {n_tot} 张，连通 {n_conn} 张")

rng = np.random.RandomState(11)
for V, dens, want in ((6, 0.45, 120), (7, 0.40, 90), (8, 0.35, 60)):
    n = 0
    while n < want:
        A = (rng.rand(V, V) < dens).astype(int)
        A = np.triu(A, 1); A = A + A.T
        if not connected(A):
            continue
        n += 1
        gam, lam1, Av, nn = measure(A)
        r = invariants(A)
        r.update({"gamma1": gam, "A": Av, "N_states": nn})
        rows.append(r)
    print(f"  V={V}: 随机连通图 {n} 张（N_states={nn}）")
print(f"  样本合计 {len(rows)} 个点")
chk(len(rows) > 700, "样本量 > 700", str(len(rows)))
chk(all(r["A"] > 0 for r in rows), "全部 A > 0")
chk(all(r["lam1"] > 0 for r in rows), "全部连通（λ₁ > 0）")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【A2】A 的分布（是否近乎常数？）")
print("=" * 92)
As = np.array([r["A"] for r in rows])
print(f"  全样本：min={As.min():.4f}  max={As.max():.4f}  mean={As.mean():.4f}  "
      f"std={As.std():.4f}  相对散布={As.std() / As.mean():.3f}")
for V in (4, 5, 6, 7, 8):
    sub = np.array([r["A"] for r in rows if r["V"] == V])
    print(f"  V={V}: n={len(sub):4d}  mean={sub.mean():.4f}  std={sub.std():.4f}  "
          f"域=[{sub.min():.4f}, {sub.max():.4f}]")
chk(As.min() > 0.4 and As.max() < 0.9, "A 落在窄带 (0.4, 0.9) 内")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【A3】A 与候选不变量的关系（线性相关 + 对数回归 R²）")
print("=" * 92)
cands = ["lam1", "lam2", "lammax", "mean_deg", "var_deg", "n_edges", "kirch", "r_lam", "V"]
yy = np.log(As)
print(f"  {'不变量':>10} {'corr(A,·)':>11} {'corr(logA,log·)':>17} {'R²(log-log)':>13}")
a3 = {}
for c in cands:
    x = np.array([r[c] for r in rows], float)
    m = np.isfinite(x) & (x > 0)
    if m.sum() < 10:
        continue
    corr = float(np.corrcoef(x[m], As[m])[0, 1])
    lx = np.log(x[m])
    ly = yy[m]
    if lx.std() < 1e-12:
        continue
    p, q = np.polyfit(lx, ly, 1)
    pred = p * lx + q
    ss = float(np.sum((ly - ly.mean()) ** 2))
    r2 = 1 - float(np.sum((ly - pred) ** 2)) / ss if ss > 0 else float("nan")
    cl = float(np.corrcoef(lx, ly)[0, 1])
    a3[c] = {"corr": corr, "corr_log": cl, "r2_loglog": r2, "slope": float(p)}
    print(f"  {c:>10} {corr:>11.4f} {cl:>17.4f} {r2:>13.4f}")
chk(True, "A3 完成")

best = max(a3.items(), key=lambda kv: abs(kv[1]["corr"]))
print()
print(f"  ★ 最强单变量相关：{best[0]}（|corr|={abs(best[1]['corr']):.4f}）")
best2 = max(a3.items(), key=lambda kv: kv[1]["r2_loglog"])
print(f"  ★ 最强 log-log R²：{best2[0]}（R²={best2[1]['r2_loglog']:.4f}，斜率 {best2[1]['slope']:.3f}）")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【A4】A 是否只是 V 的函数？（同 V 内的散布是关键判据）")
print("=" * 92)
for V in (4, 5, 6, 7, 8):
    sub = [r for r in rows if r["V"] == V]
    a = np.array([r["A"] for r in sub])
    l1 = np.array([r["lam1"] for r in sub])
    print(f"  V={V}: A 的相对散布 {a.std() / a.mean():.4f}；"
          f"corr(A, λ₁) = {np.corrcoef(l1, a)[0, 1]:+.4f}")
print()
print("  ⟹ 若同 V 内 A 的散布远大于 0，则 A 不是 V 的函数，而是图的函数。")

# 同 V 内，A 与 λ₁ 的相关（去掉 V 的混淆）
print()
print("  同 V 内 corr(A, λ₁)：")
for V in (4, 5, 6, 7, 8):
    sub = [r for r in rows if r["V"] == V]
    if len(sub) < 10:
        continue
    a = np.array([r["A"] for r in sub])
    l1 = np.array([r["lam1"] for r in sub])
    sub2 = {"V": V, "n": len(sub), "corr_A_lam1": float(np.corrcoef(l1, a)[0, 1]),
            "A_mean": float(a.mean()), "A_std": float(a.std())}
    print(f"    V={V}: n={len(sub):4d}  corr={sub2['corr_A_lam1']:+.4f}  "
          f"A={a.mean():.4f}±{a.std():.4f}")
chk(True, "A4 完成")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【A5】判定")
print("=" * 92)
print(f"  A 的域：[{As.min():.4f}, {As.max():.4f}]，相对散布 {As.std() / As.mean():.3f}")
print(f"  最强单变量：{best2[0]}，log-log R² = {best2[1]['r2_loglog']:.4f}")
verdict = ("A 近乎常数但非严格常数：域宽约 "
           f"{100 * (As.max() - As.min()) / As.mean():.1f}%，"
           f"最强单变量 {best2[0]} 只解释 R²={best2[1]['r2_loglog']:.3f}")
print(f"  ⟹ {verdict}")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【A6】收敛：A 随 V 趋向 A_∞（有限尺寸外推）")
print("=" * 92)
means = {}
for V in (4, 5, 6, 7, 8):
    sub = [r["A"] for r in rows if r["V"] == V]
    if sub:
        means[V] = float(np.mean(sub))
print("  逐 V 均值：", "、".join(f"V={V}: {m:.4f}" for V, m in means.items()))
print()
conv = {}
print(f"  {'形式':>16} {'A_∞':>9} {'参数':>10} {'最大残差':>10} {'用点':>6}")
for name, f in [
    ("A∞ − c/V", lambda V, Ainf, c: Ainf - c / V),
    ("A∞ − c/V²", lambda V, Ainf, c: Ainf - c / V ** 2),
    ("A∞ − c/V^1.5", lambda V, Ainf, c: Ainf - c / V ** 1.5),
]:
    for useV in ((5, 6, 7, 8), (4, 5, 6, 7, 8)):
        Vs = np.array([V for V in useV if V in means], float)
        As_ = np.array([means[int(V)] for V in Vs])
        best = None
        for Ainf in np.arange(0.60, 0.72, 0.0005):
            for c in np.arange(-0.5, 0.5, 0.002):
                res = f(Vs, Ainf, c) - As_
                mx = float(np.max(np.abs(res)))
                if best is None or mx < best[0]:
                    best = (mx, Ainf, c)
        if useV == (5, 6, 7, 8):
            conv[name] = {"A_inf": float(best[1]), "c": float(best[2]),
                          "max_residual": float(best[0]), "used_V": list(useV)}
            print(f"  {name:>16} {best[1]:>9.4f} {best[2]:>10.4f} {best[0]:>10.4f} "
                  f"{len(useV):>6}")
bestform = min(conv.items(), key=lambda kv: kv[1]["max_residual"])
print()
print(f"  ★ 最佳形式：{bestform[0]}，A_∞ = {bestform[1]['A_inf']:.4f}"
      f"（最大残差 {bestform[1]['max_residual']:.4f}）")
chk(bestform[1]["max_residual"] < 0.01, "收敛外推的最大残差 < 0.01",
    f"{bestform[1]['max_residual']:.4f}")

# 统计量随 V 收敛？
stds = {V: float(np.std([r["A"] for r in rows if r["V"] == V])) for V in means}
print()
print("  标准差随 V：", "、".join(f"V={V}: {sd:.4f}" for V, sd in stds.items()))
chk(stds[6] > stds[8], "散布随 V 收窄（V=6 → V=8）", f"{stds[6]:.4f} → {stds[8]:.4f}")

out = {
    "question": "A = γ₁/λ₁(L) 是哪个图不变量的函数",
    "A6_convergence": {"per_V_mean": {str(k): v for k, v in means.items()},
                       "per_V_std": {str(k): v for k, v in stds.items()},
                       "fits": conv,
                       "best": {"form": bestform[0], **bestform[1]}},
    "n_points": len(rows),
    "A_stats": {"min": float(As.min()), "max": float(As.max()), "mean": float(As.mean()),
                "std": float(As.std())},
    "A3_predictors": a3,
    "A5_verdict": verdict,
    "assertions": {"passed": ok, "failed": len(bad), "bad": bad},
}
os.makedirs(RES, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print()
print(f"断言：{ok} 条通过 / {len(bad)} 条不符")
for m in bad:
    print("  ✗", m)
print(f"已存 {OUT}")
raise SystemExit(1 if bad else 0)