#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""窗口化 16/8/4：L=8 动力学 + L=4 读出窗口的 4 个 site。

类 c=(site j, sign, parity):
  j=0..3 是 L=4 读出窗口内的 4 个位置；
  sign=该位置字母 ±1；
  parity=窗口面积奇偶。
16 类；中=对 parity 求和（8）；顶=对 sign,parity 求和（4）。
毁灭：底 16 全灭；中留最高 site 的两个 sign；顶 4 个 site 不动。
reseed 用窗口类 rho 的 push-forward；未毁层继续收 seed/flow。
输出：results/zero_window_16.json/.png
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
TARGET = np.array([0.1724, 0.2759, 0.5517])
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


def main():
    W8 = list(zs.balanced_words(8))
    cnt = np.zeros((4, 2, 2))
    nwin = 0
    for w in W8:
        for i in range(5):
            sub = w[i:i + 4]
            if sum(sub) != 0:
                continue
            par = zs.word_area(sub) % 2
            for j in range(4):
                cnt[j, 0 if sub[j] > 0 else 1, par] += 1
            nwin += 1
    rho = cnt / cnt.sum()
    print(f"W8={len(W8)}, balanced 4-windows={nwin}, rho sum={rho.sum():.12f}")
    print("rho per site:")
    for j in range(4):
        print(f"  site {j}:", np.round(rho[j], 6).tolist())

    bottom = rho.reshape(-1).copy(); bottom /= bottom.sum()
    middle = rho.sum(axis=2).reshape(-1).copy(); middle /= middle.sum()
    top = rho.sum(axis=(1, 2)).copy(); top /= top.sum()
    k, seed = 0.2, 0.05
    rows = []
    for gen in range(1, 201):
        dest = gen % 4 == 0
        pre_top = top.copy()
        if dest:
            jstar = int(np.argmax(rho.sum(axis=(1, 2))))
            bottom[:] = 0.0
            mid = middle.reshape(4, 2)
            keep = np.zeros_like(mid); keep[jstar, :] = mid[jstar, :]
            middle = keep.reshape(-1)
        j32 = k * top.sum()
        top -= top * (j32 / max(top.sum(), 1e-300))
        rho_mid = rho.sum(axis=2).reshape(-1); rho_mid /= rho_mid.sum()
        middle += j32 * rho_mid
        j21 = k * middle.sum()
        middle -= middle * (j21 / max(middle.sum(), 1e-300))
        rho_bot = rho.reshape(-1); rho_bot /= rho_bot.sum()
        bottom += j21 * rho_bot
        top += seed * rho.sum(axis=(1, 2))
        tot = bottom.sum() + middle.sum() + top.sum()
        bottom /= tot; middle /= tot; top /= tot
        rows.append({"gen": gen, "destroyed": dest,
                     "shares": [float(bottom.sum()), float(middle.sum()), float(top.sum())],
                     "top_change": float(np.abs(top - pre_top).sum())})
        if gen <= 8 or dest:
            print(f"gen={gen:>3} D={int(dest)} shares={np.round(rows[-1]['shares'],4)}")

    sh = np.array([r["shares"] for r in rows])
    mean = sh.mean(axis=0)
    dist = float(np.abs(mean - TARGET).sum())
    phases = [sh[p::4].mean(axis=0).tolist() for p in range(4)]
    print("=" * 70); print("判据"); print("=" * 70)
    ck("4 site ×2×2 = 16 类", rho.shape == (4, 2, 2), str(rho.shape))
    ck("rho 归一", abs(rho.sum() - 1) < 1e-12, f"{rho.sum():.12f}")
    ck("180 个平衡 4-窗口", nwin == 180, str(nwin))
    ck("毁灭代底份额低", all(r["shares"][0] < r["shares"][1] and r["shares"][0] < r["shares"][2]
                         for r in rows if r["destroyed"]), "")
    ck("未毁顶在变", max(r["top_change"] for r in rows) > 1e-6,
       f"max={max(r['top_change'] for r in rows):.3e}")
    ns = sh[:, 0]
    fft = np.abs(np.fft.rfft(ns - ns.mean()))
    ck("底有周期 4 分量", fft[50] > np.median(fft[1:]), f"bin50={fft[50]:.3g}")
    print("mean", np.round(mean, 4), "dist", round(dist, 4), "phases", np.round(phases, 4).tolist())

    g = np.arange(1, 201)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].plot(g, sh[:, 0], label="bottom 16", lw=1)
    ax[0].plot(g, sh[:, 1], label="middle 8", lw=1)
    ax[0].plot(g, sh[:, 2], label="top 4", lw=1)
    ax[0].set_title("(a) windowed 16/8/4"); ax[0].legend(fontsize=8)
    ax[1].plot(g, [r["top_change"] for r in rows], lw=1)
    ax[1].set_title("(b) top reseed change")
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout()
    fp = os.path.join(RES, "zero_window_16.png")
    fig.savefig(fp, dpi=130)
    out = {"windows": nwin, "rho": rho.tolist(), "mean": mean.tolist(),
           "dist": dist, "phases": phases, "rows": rows,
           "checks": CK, "n_fail": sum(1 for c in CK if not c["pass"]),
           "anchor": "L=8 dynamics, L=4 readout window; 4 site x 2 x 2; no dimensionful anchor"}
    with open(os.path.join(RES, "zero_window_16.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_window_16.json"), "和", fp)


if __name__ == "__main__":
    main()
