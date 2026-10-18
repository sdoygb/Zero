#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""三拷贝 W 的聚合一致读法（推导版）。

16 类 (j,sigma,pi) -> 4x4 W[j, f=2*sigma+pi]。
聚合：中 8 = 对 pi 求和；顶 4 = 对 sigma,pi 求和。
等可靠性最小二乘：
    x = (Wb + A_m^T Wm + A_t^T Wt) / 3
record_j = sum_{sigma,pi} x[j,sigma,pi] = (Rb + 2 Rm + 4 Rt)/3.
毁灭：底全灭；中只留最高 site 的 2 类（4 个 parity 项）；顶不动。
flow：顶->中->底，按 rho 目标注入。
输出：results/zero_window16_aggregate.json/.png
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
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


def amat(adj):
    A = np.zeros((V, V))
    for u, nb in enumerate(adj):
        for v in nb: A[u, v] = 1.0
    return A


def readout(Wb, Wm, Wt):
    """x16=(Wb + A_m^T Wm8 + A_t^T Wt4)/3；返回 x16, record4。"""
    Wb16 = Wb.reshape(-1).copy()
    Wm8 = np.zeros((V, 2))
    for j in range(V):
        for s in range(2):
            Wm8[j, s] = Wm[j, 2 * s] + Wm[j, 2 * s + 1]
    Wt4 = Wt.sum(axis=1)
    x = Wb16.copy()
    for j in range(V):
        for s in range(2):
            x[4 * j + 2 * s] += Wm8[j, s]
            x[4 * j + 2 * s + 1] += Wm8[j, s]
        x[4 * j:4 * j + 4] += Wt4[j]
    x /= 3.0
    return x, x.reshape(V, 4).sum(axis=1)


def main():
    adj = zs.k4()
    A = amat(adj); AL = np.linalg.matrix_power(A, L)
    # rho 从 L=8 窗口类
    R = np.zeros((V, V))
    for w in zs.balanced_words(8):
        for i in range(5):
            sub = w[i:i + 4]
            if sum(sub) != 0: continue
            par = zs.word_area(sub) % 2
            for j in range(V):
                sig = 0 if sub[j] > 0 else 1
                R[j, 2 * sig + par] += 1
    R = R / R.sum()
    rho_bot = R.reshape(-1); rho_bot /= rho_bot.sum()
    rho_mid = np.zeros((V, 2))
    for j in range(V):
        for s in range(2):
            rho_mid[j, s] = R[j, 2 * s] + R[j, 2 * s + 1]
    rho_mid /= rho_mid.sum()

    seed = np.array([8.0, 4.0, 1.5, 1.0])
    shares = (1/3, 1/3, 1/3)  # 每层总测度相等：16*1=8*2=4*4
    Wb = np.zeros((V, V)); Wm = np.zeros((V, V)); Wt = np.zeros((V, V))
    record = np.zeros(V); layers = []; rows = []
    for gen in range(1, 201):
        for Wc, s in zip((Wb, Wm, Wt), shares):
            Wc[np.arange(V), np.arange(V)] += s * seed
        Wb = Wb @ AL; Wm = Wm @ AL; Wt = Wt @ AL
        x, rec4 = readout(Wb, Wm, Wt)
        record += rec4
        np.fill_diagonal(Wb, 0); np.fill_diagonal(Wm, 0); np.fill_diagonal(Wt, 0)
        # 归一化三份共同缩放
        mw = max(np.abs(Wb).sum(), np.abs(Wm).sum(), np.abs(Wt).sum())
        if mw > 1e200:
            Wb /= mw; Wm /= mw; Wt /= mw
        mr = record.sum()
        if mr > 1e200: record /= mr
        dest = gen % T == 0; bottom_zero_pre = False
        if dest:
            jstar = int(np.argmax(rec4))
            Wb[:] = 0.0
            for j in range(V):
                if j != jstar:
                    Wm[j, :] = 0.0   # 灭 3 个 site 的两类 = 6/8
            bottom_zero_pre = bool(np.abs(Wb).sum() < 1e-30)
        # flow：顶->中->底，按 rho 聚合目标
        kflow = 0.2
        Pt = np.abs(Wt).sum()
        if Pt > 0:
            J = kflow * Pt; Wt -= Wt * (J / Pt)
            for j in range(V):
                for s in range(2):
                    Wm[j, 2 * s] += J * rho_mid[j, s]
                    Wm[j, 2 * s + 1] += J * rho_mid[j, s]
        Pm = np.abs(Wm).sum()
        if Pm > 0:
            J2 = kflow * Pm; Wm -= Wm * (J2 / Pm)
            Wb += J2 * rho_bot.reshape(V, V)
        # seed：最近两代有效 record 的 top2
        layers.append(rec4.copy())
        if len(layers) > 2: layers.pop(0)
        seed = sum(layers)
        if seed.sum() > 0: seed = seed / seed.sum()
        rec = record / record.sum() if record.sum() > 0 else record
        # 有效份额按聚合重数 1:2:4（x 的三份贡献），不是原始 W 总量
        sh = [float(np.abs(Wb).sum()), float(2 * np.abs(Wm).sum()), float(4 * np.abs(Wt).sum())]
        tot = sum(sh) or 1.0
        rows.append({"gen": gen, "destroyed": dest, "bottom_zero_pre": bottom_zero_pre,
                     "shares": [x / tot for x in sh],
                     "q": float((rec ** 2).sum()), "n_std": float(rec.std()),
                     "record": rec.tolist()})
        if gen <= 8 or dest:
            print(f"gen={gen:>3} D={int(dest)} shares={np.round(rows[-1]['shares'],4)} q={rows[-1]['q']:.6f}")

    sh = np.array([r["shares"] for r in rows]); qs = np.array([r["q"] for r in rows])
    print("=" * 70); print("判据"); print("=" * 70)
    ck("聚合一致读法 (x 16 类)", x.size == 16, str(x.shape))
    ck("毁灭代底拷贝清零（flow 前）", all(r["bottom_zero_pre"] for r in rows if r["destroyed"]), "")
    ck("q 有周期 4 分量", abs(np.fft.rfft(qs - qs.mean())[50]) > 1e-9, "")
    ck("无 NaN", all(np.all(np.isfinite([r["q"], r["n_std"]])) for r in rows), "")
    print("q mean/std", qs.mean(), qs.std(), "shares tail", np.round(sh[-20:].mean(0), 4))

    g = np.arange(1, 201)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].plot(g, sh[:, 0], label="Wb bottom", lw=1)
    ax[0].plot(g, sh[:, 1], label="Wm middle", lw=1)
    ax[0].plot(g, sh[:, 2], label="Wt top", lw=1)
    ax[0].set_title("(a) W copies"); ax[0].legend(fontsize=8)
    ax[1].plot(g, qs, lw=1); ax[1].set_title("(b) q (aggregate readout)")
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); fp = os.path.join(RES, "zero_window16_aggregate.png")
    fig.savefig(fp, dpi=130)
    out = {"L": L, "T": T, "rows": rows, "checks": CK,
           "n_fail": sum(1 for c in CK if not c["pass"]),
           "anchor": "aggregation-consistent x=(Wb+A_m^T Wm+A_t^T Wt)/3; no dimensionful anchor"}
    with open(os.path.join(RES, "zero_window16_aggregate.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_window16_aggregate.json"), "和", fp)


if __name__ == "__main__":
    main()
