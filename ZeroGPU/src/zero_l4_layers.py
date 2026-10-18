#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""L=4（D=4 锁定）下的 8/4/2 分级：16/8/4 的 L=4-native 替代。

L=4 -> H=2。L1 类 (h,sigma,pi): 2*2*2=8（底）；
中层对 pi 求和 -> 4；顶层对 sigma,pi 求和 -> 2。
毁灭：底 8 全灭；中只留最高 h* 的两个 sigma；顶 2 不动。
reseed 用 L1 词类 rho 的 push-forward；未毁层继续收 seed/flow。
输出：results/zero_l4_layers.json/.png
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
L = 4
H = L // 2
TARGET = np.array([0.1724, 0.2759, 0.5517])
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


def height(w):
    s = 0; h = 0
    for x in w:
        s += x; h = max(h, abs(s))
    return h


def main():
    words = zs.balanced_words(L)
    cnt = np.zeros((H, 2, 2))
    for w in words:
        cnt[height(w) - 1, 0 if w[0] > 0 else 1, zs.word_area(w) % 2] += 1
    rho = cnt / cnt.sum()
    print(f"L={L}, words={len(words)}, H={H}")
    print("rho:")
    for h in range(H):
        for s in range(2):
            for p in range(2):
                print(f"  h={h+1} sig={'+' if s==0 else '-'} par={p}: {rho[h,s,p]:.6f}")

    bottom = rho.reshape(-1).copy(); bottom /= bottom.sum()
    middle = rho.sum(axis=2).reshape(-1).copy(); middle /= middle.sum()
    top = rho.sum(axis=(1, 2)).copy(); top /= top.sum()
    k, seed = 0.2, 0.05
    rows = []
    for gen in range(1, 201):
        dest = gen % 4 == 0
        pre_top = top.copy()
        if dest:
            hstar = int(np.argmax(rho.sum(axis=(1, 2))))
            bottom[:] = 0.0
            mid = middle.reshape(H, 2)
            keep = np.zeros_like(mid); keep[hstar, :] = mid[hstar, :]
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
    phases = [sh[p::4].mean(axis=0).tolist() for p in range(4)]
    print("=" * 70); print("判据"); print("=" * 70)
    ck("L=4, H=2, 类数=8", rho.size == 8 and H == 2, str(rho.shape))
    ck("8/4/2 三层", len(bottom) == 8 and len(middle) == 4 and len(top) == 2, "")
    ck("rho 归一", abs(rho.sum() - 1) < 1e-12, f"{rho.sum():.12f}")
    ck("毁灭代底份额低", all(r["shares"][0] < r["shares"][1] and r["shares"][0] < r["shares"][2]
                         for r in rows if r["destroyed"]), "")
    ck("未毁顶在变", max(r["top_change"] for r in rows) > 1e-6,
       f"max={max(r['top_change'] for r in rows):.3e}")
    ns = sh[:, 0]
    fft = np.abs(np.fft.rfft(ns - ns.mean()))
    ck("底有周期 4 分量", fft[50] > np.median(fft[1:]), f"bin50={fft[50]:.3g}")

    g = np.arange(1, 201)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].plot(g, sh[:, 0], label="bottom 8", lw=1)
    ax[0].plot(g, sh[:, 1], label="middle 4", lw=1)
    ax[0].plot(g, sh[:, 2], label="top 2", lw=1)
    ax[0].set_title("(a) 8/4/2 shares"); ax[0].legend(fontsize=8)
    ax[1].plot(g, [r["top_change"] for r in rows], lw=1)
    ax[1].set_title("(b) top reseed change")
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout()
    fp = os.path.join(RES, "zero_l4_layers.png")
    fig.savefig(fp, dpi=130)
    out = {"L": L, "H": H, "classes": cnt.tolist(), "rho": rho.tolist(),
           "mean": mean.tolist(), "phases": phases,
           "rows": rows, "checks": CK,
           "n_fail": sum(1 for c in CK if not c["pass"]),
           "anchor": "L=4 fixed; classes from L1 words; no dimensionful anchor"}
    with open(os.path.join(RES, "zero_l4_layers.json"), "w") as f:
        json.dump(out, f, indent=1)
    print("mean", np.round(mean, 4), "phases", np.round(phases, 4).tolist())
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_l4_layers.json"), "和", fp)


if __name__ == "__main__":
    main()
