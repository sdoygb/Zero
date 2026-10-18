#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""C′：顶偏置——只约束中、不约束顶；三层 seed 按等测度 1/3。

三份独立 W；每 T=4 底全灭、中留最高 site、顶不动；
flow 顶->中->底；软约束只作用在中：
    Wm <- (1-λ) Wm + λ * (底的 parity 聚合，每 parity 各 m/2)
顶完全不约束，保持独立；seed 份额按等测度 (1/3,1/3,1/3)。
扫 λ，看顶能否上到 0.55，且权重稳定。
输出：results/zero_window16_toppias.json/.png
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
TARGET_Q = 0.06830296233
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


def run(lam, steps=200, seed0=None):
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

    seed = np.array(seed0 if seed0 is not None else [8.0, 4.0, 1.5, 1.0])
    shares = (1/3, 1/3, 1/3)
    Wb = np.zeros((V, V)); Wm = np.zeros((V, V)); Wt = np.zeros((V, V))
    record = np.zeros(V); rows = []
    for gen in range(1, steps+1):
        for Wc, s in zip((Wb, Wm, Wt), shares):
            Wc[np.arange(V), np.arange(V)] += s * seed
        Wb = Wb @ AL; Wm = Wm @ AL; Wt = Wt @ AL
        rb = np.diag(Wb).copy(); rm = np.diag(Wm).copy(); rt = np.diag(Wt).copy()
        record += rb + rm + rt
        np.fill_diagonal(Wb, 0); np.fill_diagonal(Wm, 0); np.fill_diagonal(Wt, 0)
        mw = max(np.abs(Wb).sum(), np.abs(Wm).sum(), np.abs(Wt).sum())
        if mw > 1e200: Wb /= mw; Wm /= mw; Wt /= mw
        mr = record.sum()
        if mr > 1e200: record /= mr
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
        # C′：只拉中，不拉顶
        Wm_t = np.zeros((V, V))
        for j in range(V):
            for s in range(2):
                m = Wb[j, 2*s] + Wb[j, 2*s+1]
                Wm_t[j, 2*s] = Wm_t[j, 2*s+1] = 0.5 * m
        Wm = (1-lam) * Wm + lam * Wm_t
        mb, mm, mt = np.abs(Wb).sum(), np.abs(Wm).sum(), np.abs(Wt).sum()
        msum = mb + mm + mt
        w = np.array([mb, mm, mt]) / msum if msum > 0 else np.zeros(3)
        def norm(x):
            s = x.sum(); return x / s if s > 0 else np.full(V, 1/V)
        Rres = (mb*norm(rb) + mm*norm(rm) + mt*norm(rt)) / msum if msum > 0 else np.full(V, 1/V)
        raw = record / record.sum() if record.sum() > 0 else record
        p16 = np.abs(Wb).reshape(-1)
        q16 = float(((p16 / p16.sum()) ** 2).sum()) if p16.sum() > 0 else np.nan
        rows.append({"gen": gen, "destroyed": dest, "w": w.tolist(),
                     "q_res": float((Rres**2).sum()), "q_raw": float((raw**2).sum()),
                     "q16": q16})
    tail = rows[-40:]
    Wt_arr = np.array([r["w"] for r in tail])
    wmean = Wt_arr.mean(0); wstd = Wt_arr.std(0)
    # 周期内涨落 vs 周期间漂移：按 4 代一块取周期均值
    ncyc = len(Wt_arr) // 4
    cyc = Wt_arr[:4*ncyc].reshape(ncyc, 4, 3).mean(axis=1)
    cyc_mean = cyc.mean(0); cyc_std = cyc.std(0)
    relstd_cycle = (cyc_std / np.maximum(cyc_mean, 1e-30)).max()
    qres = np.mean([r["q_res"] for r in tail]); qraw = np.mean([r["q_raw"] for r in tail])
    q16 = np.mean([r["q16"] for r in tail])
    dist = float(np.abs(wmean - TARGET).sum())
    relstd = (wstd / np.maximum(wmean, 1e-30)).max()
    return {"lam": lam, "wmean": wmean.tolist(), "wstd": wstd.tolist(),
            "relstd": float(relstd), "relstd_cycle": float(relstd_cycle),
            "cyc_mean": cyc_mean.tolist(), "cyc_std": cyc_std.tolist(),
            "q_res": float(qres), "q_raw": float(qraw), "q16": float(q16),
            "dist": dist, "rows": rows}


def main():
    lams = [0.0, 0.25, 0.5, 0.75, 1.0]
    res = [run(l) for l in lams]
    print("lam   w_b     w_m     w_t     relstd  q_res   q_raw   dist")
    for r in res:
        print(f"{r['lam']:.2f}  " + "  ".join(f"{x:.4f}" for x in r["wmean"]) +
              f"   {r['relstd']:.4f} {r['relstd_cycle']:.4f}  {r['q_res']:.4f}  {r['q16']:.4f}  {r['dist']:.4f}")
    print("=" * 70)
    print("λ=1 + 均匀 seed [1,1,1,1]（K4 无优先顶点）")
    uni = run(1.0, seed0=[1.0, 1.0, 1.0, 1.0])
    print(f"uniform: w={np.round(uni['wmean'],4)} q16={uni['q16']:.4f} "
          f"dist={uni['dist']:.4f} relstd_cycle={uni['relstd_cycle']:.4f}")
    print("=" * 70); print("C′ 判据"); print("=" * 70)
    ck("顶不被压到 1/4（w_t > 0.25）", all(r["wmean"][2] > 0.25 for r in res[1:]),
       str([round(r["wmean"][2],3) for r in res]))
    ck("存在 λ 使份额接近观测（dist<0.1）", any(r["dist"] < 0.1 for r in res),
       str([round(r["dist"],3) for r in res]))
    ck("存在 λ 使周期间权重稳定（relstd_cycle<0.1）", any(r["relstd_cycle"] < 0.1 for r in res),
       str([round(r["relstd_cycle"],3) for r in res]))
    ck("q16 接近观测 0.068（|Δ|<0.02）", any(abs(r["q16"] - TARGET_Q) < 0.02 for r in res),
       str([round(r["q16"],4) for r in res]))
    ck("无 NaN", all(all(np.isfinite(r["wmean"])) for r in res), "")

    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    for i, lab in enumerate(["w_b", "w_m", "w_t"]):
        ax[0].plot(lams, [r["wmean"][i] for r in res], marker="o", label=lab)
    ax[0].set_title("(a) weights vs λ"); ax[0].legend(fontsize=8)
    ax[1].plot(lams, [r["dist"] for r in res], marker="o", label="dist to obs")
    ax[1].plot(lams, [r["relstd_cycle"] for r in res], marker="o", label="cycle relstd")
    ax[1].set_title("(b) dist / stability"); ax[1].legend(fontsize=8)
    ax[2].plot(lams, [r["q_res"] for r in res], marker="o", label="q_res")
    ax[2].plot(lams, [r["q_raw"] for r in res], marker="o", label="q_raw")
    ax[2].axhline(TARGET_Q, color="k", ls=":", lw=1, label="obs")
    ax[2].set_title("(c) q"); ax[2].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); fp = os.path.join(RES, "zero_window16_toppias.png")
    fig.savefig(fp, dpi=130)
    out = {"results": res, "uniform": uni, "checks": CK, "n_fail": sum(1 for c in CK if not c["pass"]),
           "anchor": "middle-only constraint; equal level seed; top independent; no dimensionful anchor"}
    with open(os.path.join(RES, "zero_window16_toppias.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_window16_toppias.json"), "和", fp)


if __name__ == "__main__":
    main()
