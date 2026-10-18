#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""由 L1 记录 push-forward 得到的 16/8/4 reseed（推导版 prototype）。

L1（L=8，H=4）的词类：c=(h,σ,π)
  h=1..4 高度；σ=首字母符号（Z2）；π=面积奇偶。
  16 类；中间层对 π 求和（8）；顶层对 σ,π 求和（4）。
reseed 权重 rho(c)=该类词数/总词数；不设自由 profile。

模型：三层持久通道；顶收 L1 seed，顶→中→底下行流；
每 T=4：底 16 全灭；中只留最高质量基 h* 的两个 σ；顶不动。
未毁灭层继续收 seed/流（reseed），被毁层从 rho 目标重新长。
输出：results/zero_reseed_derivation.json、.png
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
L = 8
H = L // 2
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


def height(w):
    s = 0
    h = 0
    for x in w:
        s += x
        h = max(h, abs(s))
    return h


def classes():
    words = zs.balanced_words(L)
    cnt = np.zeros((H, 2, 2), dtype=float)  # h, sigma(+-), parity
    for w in words:
        h = height(w) - 1
        sig = 0 if w[0] > 0 else 1
        par = zs.word_area(w) % 2
        cnt[h, sig, par] += 1
    rho = cnt / cnt.sum()
    return cnt, rho


def main():
    cnt, rho = classes()
    print(f"L={L}, H={H}, words={int(cnt.sum())}")
    print("rho(h,sigma,parity):")
    for h in range(H):
        for sig in range(2):
            for par in range(2):
                print(f"  h={h+1} sig={'-' if sig==0 else '+'} par={par}: {rho[h,sig,par]:.6f}")

    # 三层通道
    bottom = rho.reshape(-1).copy()                 # 16
    middle = rho.sum(axis=2).reshape(-1).copy()     # 8: h,sigma
    top = rho.sum(axis=(1, 2)).copy()               # 4: h
    # 每层各自归一到 1（初始）
    bottom /= bottom.sum(); middle /= middle.sum(); top /= top.sum()
    T = 4
    k = 0.2
    seed = 0.05
    rows = []
    for gen in range(1, 201):
        destroyed = (gen % T == 0)
        pre_top = top.copy()
        if destroyed:
            # 底全灭；中只留最高质量 h* 的两个 sigma；顶不动
            hstar = int(np.argmax(rho.sum(axis=(1, 2))))
            bottom[:] = 0.0
            mid = middle.reshape(H, 2)
            keep = np.zeros_like(mid); keep[hstar, :] = mid[hstar, :]
            middle = keep.reshape(-1)
        # 下行流：顶→中→底（按目标 rho 分布注入）
        j32 = k * top.sum()
        top -= top * (j32 / max(top.sum(), 1e-300))
        # push-forward：注入按 L1 的 rho 目标，不按当前（可能为零）分布
        rho_mid = rho.sum(axis=2).reshape(-1); rho_mid = rho_mid / rho_mid.sum()
        middle += j32 * rho_mid
        j21 = k * middle.sum()
        middle -= middle * (j21 / max(middle.sum(), 1e-300))
        rho_bot = rho.reshape(-1); rho_bot = rho_bot / rho_bot.sum()
        bottom += j21 * rho_bot
        # L1 seed 到顶（继续播种，未毁灭层也收）
        top += seed * rho.sum(axis=(1, 2))
        # 归一化总份额
        tot = bottom.sum() + middle.sum() + top.sum()
        if tot > 0:
            bottom /= tot; middle /= tot; top /= tot
        sh = np.array([bottom.sum(), middle.sum(), top.sum()])
        rows.append({"gen": gen, "destroyed": destroyed,
                     "shares": sh.tolist(),
                     "top_change": float(np.abs(top - pre_top).sum())})
        if gen <= 8 or destroyed:
            print(f"gen={gen:>3} D={int(destroyed)} shares={np.round(sh,4)} "
                  f"top_change={np.abs(top-pre_top).sum():.4f}")

    print("=" * 70)
    print("判据")
    print("=" * 70)
    ck("类数 = 16（4×2×2）", rho.size == 16, str(rho.shape))
    ck("rho 归一", abs(rho.sum() - 1) < 1e-12, f"{rho.sum():.12f}")
    dd = [r for r in rows if r["destroyed"]]
    ck("毁灭事件每 T=4", len(dd) == 200 // T, f"n={len(dd)}")
    # 毁灭后（下一步记录的是流后）底份额应小；看毁后瞬间需另存，这里用 event gen 的 shares
    ck("毁灭代底份额显著低于中/顶",
       all(r["shares"][0] < r["shares"][1] and r["shares"][0] < r["shares"][2] for r in dd),
       "bottom < middle, top")
    ck("未毁灭顶仍在变（被 reseed）",
       max(r["top_change"] for r in rows) > 1e-6,
       f"max={max(r['top_change'] for r in rows):.3e}")
    # 周期 4
    ns = np.array([r["shares"][0] for r in rows])
    fft = np.abs(np.fft.rfft(ns - ns.mean()))
    k4 = fft[50] if len(fft) > 50 else 0.0  # 200/4=50
    ck("底份额有周期 4 分量", k4 > np.median(fft[1:]), f"bin50={k4:.3g}")

    g = np.array([r["gen"] for r in rows])
    sh = np.array([r["shares"] for r in rows])
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].plot(g, sh[:, 0], label="bottom 16", lw=1)
    ax[0].plot(g, sh[:, 1], label="middle 8", lw=1)
    ax[0].plot(g, sh[:, 2], label="top 4", lw=1)
    ax[0].set_title("(a) shares (post-flow)"); ax[0].legend(fontsize=8)
    ax[1].plot(g, [r["top_change"] for r in rows], lw=1)
    ax[1].set_title("(b) top reseed change")
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout()
    fp = os.path.join(RES, "zero_reseed_derivation.png")
    fig.savefig(fp, dpi=130)

    out = {"L": L, "H": H, "classes": cnt.tolist(), "rho": rho.tolist(),
           "rows": rows, "checks": CK,
           "n_fail": sum(1 for c in CK if not c["pass"]),
           "anchor": "rho from L1 word counts; no free profile; no dimensionful anchor"}
    with open(os.path.join(RES, "zero_reseed_derivation.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_reseed_derivation.json"), "和", fp)


if __name__ == "__main__":
    main()
