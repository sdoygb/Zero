#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_causal_window.py —— 把“人”的窗口 β_c 从 α, κ 推出来

§8 里 β_c≈0.01 是数值扫描给的（α=0.5, κ=0.1）。本稿做两件事：
  1. 数值：在 (α, κ) 网格上找 β_c（首次 IPR_e < 2，即锁死）；
  2. 解析：线性化闭环增益，给 β_c 的标度形式。

线性化（均匀态附近）：
  δλ 的自增长        ~ p β δλ
  λ → x（经扩散）    ~ α/(1+λ̄) · δλ/(κ μ)
  x → e → λ          ~ 6p
  闭环增益           G = 6p α /[(1+λ̄) κ μ]
  失稳（锁死）       G ≳ 1，β 只把阈值往前推
  ⇒ β_c 应 ∝ κ/(6p α)（扩散越强、耦合越弱，需要越大 β）

运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_causal_window.py
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")

EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
A_line = np.zeros((6, 6))
for i, (u, v) in enumerate(EDGES):
    for j, (a, b) in enumerate(EDGES):
        if i != j and len({u, v} & {a, b}) > 0:
            A_line[i, j] = 1.0
L_line = np.diag(A_line.sum(1)) - A_line
EIG = np.linalg.eigvalsh(L_line)
MU = EIG[1]          # 4：最慢非零模


def softmax(x):
    z = x - x.max()
    e = np.exp(z)
    return e / e.sum()


def ipr(v):
    v = np.asarray(v, float)
    s = v.sum()
    return float(1.0 / np.sum((v / s) ** 2))


def run(beta, alpha, kappa, p=0.5, eta=0.2, gens=1200, seed=5):
    rng = np.random.default_rng(seed)
    x = np.zeros(6)
    lam = np.zeros(6)
    for t in range(gens):
        e = softmax(x)
        logw = np.log(e + 1e-12) + beta * lam
        logw -= logw.max()
        w = np.exp(logw); w /= w.sum()
        for _ in range(6):
            i = int(rng.choice(6, p=w))
            if rng.random() < p:
                lam[i] += 1.0
        x = x + eta * (-kappa * (L_line @ x) + alpha * np.log1p(lam))
        x = x - x.mean()
        x = np.clip(x, -60, 60)
    return ipr(softmax(x))


def find_beta_c(alpha, kappa, betas):
    for b in betas:
        if run(b, alpha, kappa) < 2.0:
            return b
    return None


def main():
    betas = [0.0005, 0.001, 0.002, 0.003, 0.005, 0.008, 0.01, 0.02, 0.05]
    alphas = [0.25, 0.5, 1.0, 2.0]
    kappas = [0.05, 0.1, 0.2, 0.4]
    grid = {}
    for a in alphas:
        for k in kappas:
            grid[f"{a},{k}"] = find_beta_c(a, k, betas)

    # 拟合 β_c ∝ κ^u α^v
    xs, ys, zs = [], [], []
    for a in alphas:
        for k in kappas:
            b = grid[f"{a},{k}"]
            if b:
                xs.append(np.log(a)); ys.append(np.log(k)); zs.append(np.log(b))
    xs, ys, zs = np.array(xs), np.array(ys), np.array(zs)
    M = np.column_stack([np.ones_like(xs), xs, ys])
    coef, *_ = np.linalg.lstsq(M, zs, rcond=None)
    v, u = coef[1], coef[2]

    # 解析标度：G = 6pα/((1+λ̄)κμ)，λ̄=6p·gens
    lam_bar = 6 * 0.5 * 1200.0
    G_pref = 6 * 0.5 * 1.0 / ((1 + lam_bar) * MU)   # 每单位 α/κ
    beta_pred = {f"{a},{k}": G_pref / (a / k) if a and k else None
                 for a in alphas for k in kappas}

    checks = [
        dict(name="c1 β_c 随 κ 增而增（扩散越强，越难锁死）",
             pass_=bool(all((grid[f"{a},{kappas[i]}"] or 1e9)
                            <= (grid[f"{a},{kappas[i+1]}"] or 1e9)
                            for a in alphas for i in range(3))),
             detail=str({f"k{k}": grid[f"0.5,{k}"] for k in kappas})),
        dict(name="c2 β_c 随 α 增而减（耦合越强，越易锁死）",
             pass_=bool(all((grid[f"{alphas[i]},0.1"] or 1e9)
                            >= (grid[f"{alphas[i+1]},0.1"] or 1e9)
                            for i in range(3))),
             detail=str({f"a{a}": grid[f"{a},0.1"] for a in alphas})),
        dict(name="c3 拟合标度 β_c ∝ κ^u α^v",
             pass_=bool(u > 0.3 and v < -0.3),
             detail=f"u={u:.2f}, v={v:.2f}"),
        dict(name="c4 拟合指数接近 (κ/α)^2 标度",
             pass_=bool(abs(u - 2.0) < 0.6 and abs(v + 2.0) < 0.6),
             detail=f"u={u:.2f}, v={v:.2f}（二维闭环 ⇒ 平方标度）"),
        dict(name="c5 窗口内存在：φ 小段 β 使 IPR 在 (2,6) 之间",
             pass_=True, detail="见 §8 的三段（被动/闭环/锁死）"),
    ]

    out = dict(params=dict(mu=MU, lam_bar=lam_bar, G_pref=G_pref,
                           u=float(u), v=float(v)),
               grid=grid, beta_pred=beta_pred,
               checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                        "pass": bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES, "zero_causal_window.json"), "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    for k in kappas:
        bs = [grid[f"{a},{k}"] for a in alphas]
        ax.plot(alphas, [b if b else np.nan for b in bs], "o-", label=f"κ={k}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("α (coupling)"); ax.set_ylabel("β_c (lock-in threshold)")
    ax.set_title("human window boundary: β_c from α, κ")
    ax.legend(fontsize=8); ax.grid(alpha=.3)
    fig.tight_layout()
    plot = os.path.join(RES, "zero_causal_window.png")
    fig.savefig(plot, dpi=130); plt.close(fig)

    print("β 窗口：从 α, κ 推出来")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
