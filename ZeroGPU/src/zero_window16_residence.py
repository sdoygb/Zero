#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""B 判定：三份独立 W + 顶不动；record 权重从 flow 的驻留质量给。

三份拷贝 Wb/Wm/Wt 独立演化；每 T=4 底全灭、中留最高 site、顶不动；
flow 顶->中->底。record 不取固定 1:2:4，而取
    R = (m_b r_b + m_m r_m + m_t r_t) / (m_b+m_m+m_t)
其中 m_c 是该拷贝当代的驻留质量，r_c 是其归一化记录。
再对最近 40 代做时间平均，得到驻留权重 w_c。
判据：权重是否稳定；q 是否比 raw 更接近观测。
输出：results/zero_window16_residence.json/.png
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
TARGET_Q = 0.06830296233
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


def main():
    adj = zs.k4(); A = np.zeros((V, V))
    for u, nb in enumerate(adj):
        for v in nb: A[u, v] = 1.0
    AL = np.linalg.matrix_power(A, L)
    R = np.zeros((V, V))
    for w in zs.balanced_words(8):
        for i in range(5):
            sub = w[i:i+4]
            if sum(sub) != 0: continue
            par = zs.word_area(sub) % 2
            for j in range(V):
                R[j, 2*(0 if sub[j] > 0 else 1) + par] += 1
    R = R / R.sum()
    rho_bot = R.reshape(-1); rho_bot /= rho_bot.sum()
    rho_mid = np.zeros((V, 2))
    for j in range(V):
        for s in range(2): rho_mid[j, s] = R[j, 2*s] + R[j, 2*s+1]
    rho_mid /= rho_mid.sum()

    seed = np.array([8.0, 4.0, 1.5, 1.0])
    shares = (0.5, 0.3, 0.2)
    Wb = np.zeros((V, V)); Wm = np.zeros((V, V)); Wt = np.zeros((V, V))
    record_raw = np.zeros(V); rows = []
    for gen in range(1, 201):
        for Wc, s in zip((Wb, Wm, Wt), shares):
            Wc[np.arange(V), np.arange(V)] += s * seed
        Wb = Wb @ AL; Wm = Wm @ AL; Wt = Wt @ AL
        rb = np.diag(Wb).copy(); rm = np.diag(Wm).copy(); rt = np.diag(Wt).copy()
        record_raw += rb + rm + rt
        np.fill_diagonal(Wb, 0); np.fill_diagonal(Wm, 0); np.fill_diagonal(Wt, 0)
        mw = max(np.abs(Wb).sum(), np.abs(Wm).sum(), np.abs(Wt).sum())
        if mw > 1e200: Wb /= mw; Wm /= mw; Wt /= mw
        mr = record_raw.sum()
        if mr > 1e200: record_raw /= mr
        dest = gen % T == 0
        if dest:
            jstar = int(np.argmax(rb + rm + rt))
            Wb[:] = 0.0
            for j in range(V):
                if j != jstar: Wm[j, :] = 0.0
        k = 0.2
        Pt = np.abs(Wt).sum()
        if Pt > 0:
            J = k * Pt; Wt -= Wt * (J / Pt)
            for j in range(V):
                for s in range(2):
                    Wm[j, 2*s] += J * rho_mid[j, s]
                    Wm[j, 2*s+1] += J * rho_mid[j, s]
        Pm = np.abs(Wm).sum()
        if Pm > 0:
            J2 = k * Pm; Wm -= Wm * (J2 / Pm)
            Wb += J2 * rho_bot.reshape(V, V)
        # 记录（归一化）与驻留质量
        mb, mm, mt = np.abs(Wb).sum(), np.abs(Wm).sum(), np.abs(Wt).sum()
        msum = mb + mm + mt
        w = np.array([mb, mm, mt]) / msum if msum > 0 else np.zeros(3)
        def norm(x):
            s = x.sum()
            return x / s if s > 0 else np.full(V, 1/V)
        Rres = (mb*norm(rb) + mm*norm(rm) + mt*norm(rt)) / msum if msum > 0 else np.full(V, 1/V)
        qres = float((Rres ** 2).sum())
        raw = record_raw / record_raw.sum() if record_raw.sum() > 0 else record_raw
        qraw = float((raw ** 2).sum())
        rows.append({"gen": gen, "destroyed": dest, "w": w.tolist(),
                     "q_res": qres, "q_raw": qraw,
                     "shares": w.tolist()})
        if gen <= 8 or dest:
            print(f"gen={gen:>3} D={int(dest)} w={np.round(w,4)} q_res={qres:.6f} q_raw={qraw:.6f}")

    tail = rows[-40:]
    Wt_arr = np.array([r["w"] for r in tail])
    wmean = Wt_arr.mean(axis=0); wstd = Wt_arr.std(axis=0)
    qres_mean = np.mean([r["q_res"] for r in tail])
    qraw_mean = np.mean([r["q_raw"] for r in tail])
    print("=" * 70); print("B 判定"); print("=" * 70)
    ck("驻留权重归一", abs(wmean.sum() - 1) < 1e-12, f"{wmean.sum():.12f}")
    ck("权重稳定（std/mean < 0.1）", np.all(wstd / np.maximum(wmean, 1e-30) < 0.1),
       f"wmean={np.round(wmean,4)}, relstd={np.round(wstd/np.maximum(wmean,1e-30),4)}")
    ck("驻留 record 比 raw 更接近观测 q", abs(qres_mean - TARGET_Q) < abs(qraw_mean - TARGET_Q),
       f"q_res={qres_mean:.4f} q_raw={qraw_mean:.4f} target={TARGET_Q:.4f}")
    ck("驻留 record 的 q 有周期 4", abs(np.fft.rfft(np.array([r['q_res'] for r in rows]) -
                                      np.mean([r['q_res'] for r in rows]))[50]) > 1e-9, "")
    ck("无 NaN", all(np.all(np.isfinite([r["q_res"], r["q_raw"]])) for r in rows), "")

    g = np.arange(1, 201)
    qres = np.array([r["q_res"] for r in rows]); qraw = np.array([r["q_raw"] for r in rows])
    w = np.array([r["w"] for r in rows])
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].plot(g, w[:, 0], label="w_b", lw=1)
    ax[0].plot(g, w[:, 1], label="w_m", lw=1)
    ax[0].plot(g, w[:, 2], label="w_t", lw=1)
    ax[0].set_title("(a) residence weights"); ax[0].legend(fontsize=8)
    ax[1].plot(g, qres, label="q_res", lw=1)
    ax[1].plot(g, qraw, label="q_raw", lw=1)
    ax[1].axhline(TARGET_Q, color="k", ls=":", lw=1, label="obs")
    ax[1].set_title("(b) q"); ax[1].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); fp = os.path.join(RES, "zero_window16_residence.png")
    fig.savefig(fp, dpi=130)
    out = {"wmean": wmean.tolist(), "wstd": wstd.tolist(),
           "qres_mean": qres_mean, "qraw_mean": qraw_mean, "target_q": TARGET_Q,
           "rows": rows, "checks": CK, "n_fail": sum(1 for c in CK if not c["pass"]),
           "anchor": "residence weights from flow; no free 1:2:4; no dimensionful anchor"}
    with open(os.path.join(RES, "zero_window16_residence.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_window16_residence.json"), "和", fp)
    print("wmean", wmean, "q_res/q_raw", qres_mean, qraw_mean)


if __name__ == "__main__":
    main()
