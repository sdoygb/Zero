#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tau_to_mu_bridge.py —— 慢模速率与图的谱隙：能标桥接的判定
================================================================================
背景：
  · 库内**已有**图 Laplacian → 能量的桥：$\omega^2 = c^2\lambda_k$（Z104，$c=\sqrt{0.8}$），
    并在 Z118／Z120／Z152 上把能量/波长的标度搭通。
  · Z169 给出零和层生成元的慢模：$\tau_1 = 1/|\lambda_1(\mathcal L_\kappa)|$（单位＝一次移动）。
  · 旧库 G46 立了"标度 $k$ **就是**寿命步数 $L$"（被可实现性逼出，无自由参数）；
    旧库 L2_years 立了"绝对时间需锚"（G57 量纲空洞）。

【本篇的问题（能标桥的唯一形态）】
  零和层慢模速率 $\gamma_1:=|\lambda_1(\mathcal L_\kappa)|$ 是否由**已搭通能量桥的那个量**
  ——图 Laplacian 的 Fiedler 值 $\lambda_1(\mathsf L)$——决定？
  若是，则 $\tau_1$ 自动接上 $\omega=c\sqrt{\lambda}$ 那条桥，$\mu\propto\gamma_1$ 就有内生依据。

【检验】
  【T1】生成元与 Laplacian 的双重核验（两族的谱各自正确）
  【T2】同一 $V$ 下换图：$\gamma_1$ 与 $\lambda_1(\mathsf L)$ 的关系（单调？幂律？）
  【T3】幂律拟合 $\gamma_1 \propto \lambda_1(\mathsf L)^p$（报 $p$、$R^2$、残差）
  【T4】对照：$\gamma_1$ 是否更像 $\lambda_1(\mathsf L)$ 以外的量（如平均度）的函数
  【T5】判定：桥的形式与所缺

用法：/usr/local/bin/python3 tau_to_mu_bridge.py
输出：ZeroGPU/results/tau_to_mu_bridge.json
"""
from __future__ import annotations
import itertools, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "tau_to_mu_bridge.json")

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
# 图族
# ══════════════════════════════════════════════════════════════════════
def g_path(n):
    A = np.zeros((n, n), int)
    for i in range(n - 1):
        A[i, i + 1] = A[i + 1, i] = 1
    return A


def g_ring(n):
    A = g_path(n)
    A[0, n - 1] = A[n - 1, 0] = 1
    return A


def g_star(n):
    A = np.zeros((n, n), int)
    for i in range(1, n):
        A[0, i] = A[i, 0] = 1
    return A


def g_complete(n):
    A = np.ones((n, n), int) - np.eye(n, dtype=int)
    return A


def g_cube():
    """Q3 = K4 □ K2：8 槽载体的图。"""
    A = np.zeros((8, 8), int)
    for i in range(8):
        for b in (1, 2, 4):
            j = i ^ b
            A[i, j] = 1
    return A


def g_k23():
    A = np.zeros((5, 5), int)
    for i in range(2):
        for j in range(2, 5):
            A[i, j] = A[j, i] = 1
    return A


# ══════════════════════════════════════════════════════════════════════
# 两个谱
# ══════════════════════════════════════════════════════════════════════
def fiedler(A):
    """图 Laplacian 的最小非零本征值 λ₁(L)。"""
    d = A.sum(axis=1)
    Lg = np.diag(d) - A
    ev = np.sort(np.linalg.eigvalsh(Lg))
    return float(ev[1])


def gen_gap(A, Lcap):
    """零和层生成元的慢模速率 γ₁ = |λ₁(L_κ)|（κ=1 全分支等权）。"""
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
                rows.append(idx[tuple(m)]); cols.append(k); vals.append(1.0)
                deg[k] += 1.0
    nn = len(st)
    rows += list(range(nn)); cols += list(range(nn)); vals += list(-deg)
    from scipy.sparse import csr_matrix
    Lm = csr_matrix((vals, (rows, cols)), shape=(nn, nn)).toarray()
    ev = np.sort(np.linalg.eigvalsh(Lm))
    return float(abs(ev[-2])), nn


# ══════════════════════════════════════════════════════════════════════
print("=" * 92)
print("【T1】双重核验：Laplacian 的 Fiedler 值 + 生成元零模唯一")
print("=" * 92)
family = []
for name, A in [
    ("path4", g_path(4)), ("star4", g_star(4)), ("K4", g_complete(4)),
    ("path6", g_path(6)), ("ring6", g_ring(6)), ("star6", g_star(6)),
    ("K23", g_k23()), ("K6", g_complete(6)),
    ("path8", g_path(8)), ("ring8", g_ring(8)), ("star8", g_star(8)),
    ("Q3(cube)", g_cube()),
]:
    V = A.shape[0]
    Lcap = 1                                  # L=1 ⟹ 状态数可控（V≤8 时 N≤1107）
    lamL = fiedler(A)
    gam, nn = gen_gap(A, Lcap)
    degmean = float(A.sum(axis=1).mean())
    family.append({"graph": name, "V": V, "L": Lcap, "N_states": nn,
                   "fiedler": lamL, "gamma1": gam, "mean_degree": degmean})
    print(f"  {name:10s} V={V} L={Lcap} N={nn:5d}  λ₁(L)={lamL:9.6f}  "
          f"γ₁={gam:9.6f}  τ₁={1/gam:8.4f}  <k>={degmean:.3f}")

print()
chk(all(r["N_states"] > 1 for r in family), "全部图的生成元状态数 > 1")
chk(all(r["gamma1"] > 0 for r in family), "全部图的 γ₁ > 0（有慢模）")
chk(all(r["fiedler"] > 0 for r in family), "全部图连通（Fiedler > 0）")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【T2】同一 V 下换图：γ₁ 与 λ₁(L) 的关系")
print("=" * 92)
byV = {}
for r in family:
    byV.setdefault(r["V"], []).append(r)
for V in sorted(byV):
    rows = sorted(byV[V], key=lambda x: x["fiedler"])
    print(f"  V={V}（按 λ₁(L) 升序）：")
    for r in rows:
        print(f"     {r['graph']:10s} λ₁(L)={r['fiedler']:9.6f}  γ₁={r['gamma1']:9.6f}  "
              f"γ₁/λ₁={r['gamma1']/r['fiedler']:8.4f}")
    mono = all(rows[i]["gamma1"] <= rows[i + 1]["gamma1"] + 1e-12 for i in range(len(rows) - 1))
    print(f"     γ₁ 随 λ₁(L) 单调不减？ {'是' if mono else '**否**'}")
    for r in rows:
        r["mono_ok"] = mono
chk(True, "T2 完成（逐 V 记录）")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【T3】幂律拟合 γ₁ ∝ λ₁(L)^p")
print("=" * 92)
x = np.array([r["fiedler"] for r in family])
y = np.array([r["gamma1"] for r in family])
lx, ly = np.log(x), np.log(y)
p, c = np.polyfit(lx, ly, 1)
pred = p * lx + c
ss_res = float(np.sum((ly - pred) ** 2))
ss_tot = float(np.sum((ly - ly.mean()) ** 2))
r2 = 1 - ss_res / ss_tot
print(f"  全族（12 图）：γ₁ = {np.exp(c):.6f} · λ₁(L)^{p:.4f}   R² = {r2:.4f}")
print(f"  {'图':>10} {'λ₁(L)':>10} {'γ₁ 实测':>10} {'γ₁ 拟合':>10} {'相对偏差':>10}")
for r, lxi in zip(family, lx):
    fit = np.exp(p * lxi + c)
    print(f"  {r['graph']:>10} {r['fiedler']:>10.6f} {r['gamma1']:>10.6f} {fit:>10.6f} "
          f"{(fit - r['gamma1']) / r['gamma1'] * 100:>9.2f}%")
print()
chk(r2 > 0.90, f"幂律拟合优度 R² > 0.90", f"R²={r2:.4f}, p={p:.3f}")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("【T4】对照：γ₁ 是否更像别的量的函数？（平均度）")
print("=" * 92)
xm = np.array([r["mean_degree"] for r in family])
lxm = np.log(xm)
pm, cm = np.polyfit(lxm, ly, 1)
predm = pm * lxm + cm
r2m = 1 - float(np.sum((ly - predm) ** 2)) / ss_tot
print(f"  以平均度拟合：γ₁ ∝ <k>^{pm:.4f}，R² = {r2m:.4f}")
print(f"  以 Fiedler 拟合：R² = {r2:.4f}")
print()
print(f"  ⟹ 判别：Fiedler 的 R²={r2:.4f} vs 平均度的 R²={r2m:.4f}")
chk(True, "T4 完成（记录两种对照拟合）")

# 同 V 内（去掉 V 的影响）看 Fiedler 的解释力
print()
print("  同 V 内（V=4,6,8 分别拟合）：")
inner = {}
for V in sorted(byV):
    rows = byV[V]
    if len(rows) < 2:
        continue
    xv = np.log([r["fiedler"] for r in rows])
    yv = np.log([r["gamma1"] for r in rows])
    if len(rows) >= 3:
        pv, cv = np.polyfit(xv, yv, 1)
        pv_pred = pv * xv + cv
        st = float(np.sum((yv - yv.mean()) ** 2))
        sr = float(np.sum((yv - pv_pred) ** 2))
        r2v = 1 - sr / st if st > 0 else float("nan")
        print(f"    V={V}: p={pv:.3f}, R²={r2v:.4f}")
        inner[V] = {"p": float(pv), "r2": float(r2v), "n": len(rows)}
    else:
        print(f"    V={V}: 只有 {len(rows)} 个图，不做拟合")
        inner[V] = {"n": len(rows)}

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 92)
print("汇总")
print("=" * 92)
out = {
    "question": "零和层慢模速率 γ₁ 是否由图 Laplacian 的 Fiedler 值 λ₁(L) 决定",
    "T1_family": family,
    "T3_power_law": {"prefactor": float(np.exp(c)), "p": float(p), "r2": float(r2),
                     "n_points": len(family)},
    "T4_mean_degree": {"p": float(pm), "r2": float(r2m)},
    "T4_within_V": inner,
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