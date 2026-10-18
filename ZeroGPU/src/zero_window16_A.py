#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""A 判定：16/8/4 作为同一 W 的三分辨率读出。

只有一份底层 W（4x4，16 类）。
- 底 16 = W 的 16 项；
- 中 8 = 对 parity 求和；
- 顶 4 = 对 sign,parity 求和。
毁灭只作用底层 W：T=4 时 W 全灭；中/顶是读出，跟着底层走，没有“顶不动”。
record 从 W 的对角读；1:2:4 是读出权重（此时 x=W，权重自动一致）。
输出：results/zero_window16_A.json/.png
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zero_sim as zs  # noqa: E402

RES = os.path.join(HERE, "..", "results")
L, T, V = 4, 4, 4
TARGET = np.array([0.1724, 0.2759, 0.5517])
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


def main():
    adj = zs.k4(); A = np.zeros((V, V))
    for u, nb in enumerate(adj):
        for v in nb: A[u, v] = 1.0
    AL = np.linalg.matrix_power(A, L)
    seed = np.array([8.0, 4.0, 1.5, 1.0])
    W = np.zeros((V, V)); record = np.zeros(V); layers = []; rows = []
    for gen in range(1, 201):
        W[np.arange(V), np.arange(V)] += seed
        W = W @ AL
        clo = np.diag(W).copy()
        record += clo
        np.fill_diagonal(W, 0)
        mw = np.abs(W).sum()
        if mw > 1e200: W /= mw
        mr = record.sum()
        if mr > 1e200: record /= mr
        dest = gen % T == 0
        if dest:
            W[:] = 0.0
        layers.append(clo.copy())
        if len(layers) > 2: layers.pop(0)
        seed = sum(layers)
        if seed.sum() > 0: seed = seed / seed.sum()
        rec = record / record.sum() if record.sum() > 0 else record
        # 三层读出（同一 W 的聚合；总质量相同）
        W4 = W.reshape(V, 4)
        mid = W4.reshape(V, 2, 2).sum(axis=2).sum()  # 8 类总和 = W 总和
        top = W4.sum()
        bot = W4.sum()
        tot = bot + mid + top
        sh = [bot / tot, mid / tot, top / tot] if tot > 0 else [1/3, 1/3, 1/3]
        rows.append({"gen": gen, "destroyed": dest, "shares": sh,
                     "q": float((rec ** 2).sum()), "n_std": float(rec.std())})
        if gen <= 8 or dest:
            print(f"gen={gen:>3} D={int(dest)} shares={np.round(sh,4)} q={rows[-1]['q']:.6f}")

    sh = np.array([r["shares"] for r in rows]); qs = np.array([r["q"] for r in rows])
    print("=" * 70); print("A 判定"); print("=" * 70)
    ck("只有一份 W", W.shape == (V, V), "")
    ck("毁灭代 W 清零", all(np.abs(W).sum() == 0 for r in rows if r["destroyed"]), "post-destruction")
    ck("q 有周期 4 分量", abs(np.fft.rfft(qs - qs.mean())[50]) > 1e-9, f"bin50={abs(np.fft.rfft(qs-qs.mean())[50]):.3g}")
    ck("无 NaN", all(np.all(np.isfinite([r["q"], r["n_std"]])) for r in rows), "")
    print("q mean/std", qs.mean(), qs.std(), "shares tail", np.round(sh[-20:].mean(0), 4))
    dist = float(np.abs(sh[-20:].mean(0) - TARGET).sum())
    print("shares distance to target", dist)

    g = np.arange(1, 201)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].plot(g, qs, lw=1); ax[0].set_title("(a) q (single W)")
    ax[1].plot(g, sh[:, 0], label="bottom", lw=1)
    ax[1].plot(g, sh[:, 1], label="middle", lw=1)
    ax[1].plot(g, sh[:, 2], label="top", lw=1)
    ax[1].set_title("(b) readout shares"); ax[1].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); fp = os.path.join(RES, "zero_window16_A.png")
    fig.savefig(fp, dpi=130)
    out = {"rows": rows, "checks": CK, "n_fail": sum(1 for c in CK if not c["pass"]),
           "q_mean": float(qs.mean()), "shares_tail": sh[-20:].mean(0).tolist(),
           "dist": dist,
           "anchor": "single W; 16/8/4 as readout; no independent top; no dimensionful anchor"}
    with open(os.path.join(RES, "zero_window16_A.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_window16_A.json"), "和", fp)


if __name__ == "__main__":
    main()
