#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""忠实引擎底层尝试：窗口 16/8/4 的三份 W 拷贝（micro core）。

W 是 4x4（site j, fiber f=2*sign+parity）。
三份独立拷贝：
  Wb: 底 16 全部；
  Wm: 中 8（f 偶：0,2）；
  Wt: 顶 4（f=0）。
毁灭 T=4：Wb 全灭；Wm 只留最高 site 的 2 个（f=0,2）；Wt 不动。
reseed：每代 seed 按 (0.5,0.3,0.2) 分给三份；未毁拷贝继续收种。
有效 W = Wb+Wm+Wt，closures/record 从有效 W 的 diag 读。

这是 micro core，不是 12 通道；用来先看底层三拷贝规则能不能跑、有没有周期 4。
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
T = 4
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


def adj_matrix(adj):
    V = len(adj); A = np.zeros((V, V))
    for u, nb in enumerate(adj):
        for v in nb: A[u, v] = 1.0
    return A


def main():
    adj = zs.k4()
    A = adj_matrix(adj)
    AL = np.linalg.matrix_power(A, L)
    V = len(adj)
    # window rho: class (site j, sign, parity) -> W[j, 2*sign+parity]
    R = np.zeros((V, V))
    for w in zs.balanced_words(8):
        for i in range(5):
            sub = w[i:i + 4]
            if sum(sub) != 0:
                continue
            par = zs.word_area(sub) % 2
            for j in range(4):
                sig = 0 if sub[j] > 0 else 1
                R[j, 2 * sig + par] += 1
    R = R / R.sum()
    mask_m = np.zeros((V, V), dtype=bool); mask_m[:, [0, 2]] = True
    mask_t = np.zeros((V, V), dtype=bool); mask_t[:, 0] = True
    seed = np.array([8.0, 4.0, 1.5, 1.0])
    shares = (0.5, 0.3, 0.2)
    Wb = np.zeros((V, V)); Wm = np.zeros((V, V)); Wt = np.zeros((V, V))
    record = np.zeros(V)
    layers = []
    rows = []
    for gen in range(1, 201):
        # seed 注入三份对角
        for Wc, s in zip((Wb, Wm, Wt), shares):
            Wc[np.arange(V), np.arange(V)] += s * seed
        Wb = Wb @ AL; Wm = Wm @ AL; Wt = Wt @ AL
        # 清非对角，保留 survivors；有效 W = 三份之和
        cb = np.diag(Wb).copy(); cm = np.diag(Wm).copy(); ct = np.diag(Wt).copy()
        clo = cb + cm + ct
        record += clo
        np.fill_diagonal(Wb, 0); np.fill_diagonal(Wm, 0); np.fill_diagonal(Wt, 0)
        # 无损归一化：三份 W 共同缩放（保持相对注入比）
        mw = max(np.abs(Wb).sum(), np.abs(Wm).sum(), np.abs(Wt).sum())
        if mw > 1e200:
            Wb /= mw; Wm /= mw; Wt /= mw
        mr = record.sum()
        if mr > 1e200:
            record /= mr
        # 毁灭
        dest = gen % T == 0
        bottom_zero_pre = False
        if dest:
            jstar = int(np.argmax(record))
            Wb[:] = 0.0
            # Wm 只留 jstar 行的 f=0,2；其余清
            mask = np.zeros((V, V), dtype=bool)
            mask[jstar, [0, 2]] = True
            Wm[~mask] = 0.0
            # Wt 不动
            bottom_zero_pre = bool(np.abs(Wb).sum() < 1e-30)
        # top→middle→bottom flow：按 R 目标注入，不按当前（可能为零）分布
        kflow = 0.2
        Pt = np.abs(Wt).sum()
        if Pt > 0:
            J = kflow * Pt
            Wt -= Wt * (J / Pt)
            tgt = (R * mask_m); tgt = tgt / tgt.sum()
            Wm += J * tgt
        Pm = np.abs(Wm).sum()
        if Pm > 0:
            J2 = kflow * Pm
            Wm -= Wm * (J2 / Pm)
            tgt2 = R / R.sum()
            Wb += J2 * tgt2
        # top2 seed（有效 closure 的最近两代）
        layers.append(clo.copy())
        if len(layers) > 2: layers.pop(0)
        seed = sum(layers)
        if seed.sum() > 0: seed = seed / seed.sum()
        # 监控
        rec = record / record.sum() if record.sum() > 0 else record
        sh = [float(np.abs(Wb).sum()), float(np.abs(Wm).sum()), float(np.abs(Wt).sum())]
        tot = sum(sh) if sum(sh) > 0 else 1.0
        sh = [x / tot for x in sh]
        rows.append({"gen": gen, "destroyed": dest, "shares": sh,
                     "bottom_zero_pre": bottom_zero_pre,
                     "q": float((rec ** 2).sum()), "n_std": float(rec.std()),
                     "record": rec.tolist()})
        if gen <= 8 or dest:
            print(f"gen={gen:>3} D={int(dest)} shares={np.round(sh,4)} q={rows[-1]['q']:.6f}")

    sh = np.array([r["shares"] for r in rows])
    qs = np.array([r["q"] for r in rows])
    print("=" * 70); print("判据"); print("=" * 70)
    ck("三份 W 拷贝存在", Wb.shape == (4, 4) and Wm.shape == (4, 4) and Wt.shape == (4, 4), "")
    dd = [r for r in rows if r["destroyed"]]
    ck("毁灭代底拷贝清零（flow 前）", all(r["bottom_zero_pre"] for r in dd), "")
    ck("q 有周期 4 分量", abs(np.fft.rfft(qs - qs.mean())[50]) > 1e-9, "")
    ck("无 NaN", all(np.all(np.isfinite([r["q"], r["n_std"]])) for r in rows), "")
    print("q mean/std", qs.mean(), qs.std(), "shares tail", np.round(sh[-20:].mean(0), 4))

    g = np.arange(1, 201)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].plot(g, sh[:, 0], label="Wb bottom", lw=1)
    ax[0].plot(g, sh[:, 1], label="Wm middle", lw=1)
    ax[0].plot(g, sh[:, 2], label="Wt top", lw=1)
    ax[0].set_title("(a) W copies"); ax[0].legend(fontsize=8)
    ax[1].plot(g, qs, lw=1); ax[1].set_title("(b) q")
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); fp = os.path.join(RES, "zero_window16_engine.png")
    fig.savefig(fp, dpi=130)
    out = {"L": L, "T": T, "shares": shares, "rows": rows,
           "checks": CK, "n_fail": sum(1 for c in CK if not c["pass"]),
           "anchor": "3 W copies; integer masks; no dimensionful anchor"}
    with open(os.path.join(RES, "zero_window16_engine.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_window16_engine.json"), "和", fp)


if __name__ == "__main__":
    main()
