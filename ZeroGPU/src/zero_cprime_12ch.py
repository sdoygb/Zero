#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""C′ 接 12 通道：把 C′ 的底 W（4x4=16 类）映射成 cycle16 的 seed，测真实 content q。

C′ 在 K4 上独立演化（三份 W、flow、只约束中、顶独立）；
每代把 Wb 的 16 项归一化后作为 cycle16 引擎的 seed（4x4 -> 16 顶点）。
对照：同一 12 通道引擎、seed=np.ones(16)。
输出：results/zero_cprime_12ch.json/.png
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
from zero_evo import ZeroEvo, Params  # noqa: E402

RES = os.path.join(HERE, "..", "results")
L, T, V4 = 4, 4, 4
CK = []


def window_rho():
    R = np.zeros((V4, V4))
    for w in zs.balanced_words(8):
        for i in range(5):
            sub = w[i:i+4]
            if sum(sub) != 0:
                continue
            par = zs.word_area(sub) % 2
            for j in range(V4):
                R[j, 2*(0 if sub[j] > 0 else 1) + par] += 1
    R /= R.sum()
    rb = R.reshape(-1); rb /= rb.sum()
    rm = np.zeros((V4, 2))
    for j in range(V4):
        for s in range(2):
            rm[j, s] = R[j, 2*s] + R[j, 2*s+1]
    rm /= rm.sum()
    return R, rb, rm


RHO, RHO_B, RHO_M = window_rho()


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


def cprime_update(Wb, Wm, Wt, seed, gen):
    A = np.zeros((V4, V4))
    for u, nb in enumerate(zs.k4()):
        for v in nb: A[u, v] = 1.0
    AL = np.linalg.matrix_power(A, L)
    for Wc, s in zip((Wb, Wm, Wt), (1/3, 1/3, 1/3)):
        Wc[np.arange(V4), np.arange(V4)] += s * seed
    Wb = Wb @ AL; Wm = Wm @ AL; Wt = Wt @ AL
    mw = max(np.abs(Wb).sum(), np.abs(Wm).sum(), np.abs(Wt).sum())
    if np.isfinite(mw) and mw > 1e200:
        Wb /= mw; Wm /= mw; Wt /= mw
    if not (np.all(np.isfinite(Wb)) and np.all(np.isfinite(Wm)) and np.all(np.isfinite(Wt))):
        Wb[:] = 0.0; Wm[:] = 0.0; Wt[:] = 0.0
    clo = np.diag(Wb) + np.diag(Wm) + np.diag(Wt)
    np.fill_diagonal(Wb, 0); np.fill_diagonal(Wm, 0); np.fill_diagonal(Wt, 0)
    if gen % T == 0:
        jstar = int(np.argmax(clo))
        Wb[:] = 0.0
        for j in range(V4):
            if j != jstar: Wm[j, :] = 0.0
    # flow 顶->中->底，按窗口 rho 目标
    k = 0.2
    Pt = np.abs(Wt).sum()
    if np.isfinite(Pt) and Pt > 0:
        J = k * Pt; Wt -= Wt * (J / Pt)
        for j in range(V4):
            for s in range(2):
                Wm[j, 2*s] += J * RHO_M[j, s]
                Wm[j, 2*s+1] += J * RHO_M[j, s]
    Pm = np.abs(Wm).sum()
    if np.isfinite(Pm) and Pm > 0:
        J2 = k * Pm; Wm -= Wm * (J2 / Pm)
        Wb += J2 * RHO_B.reshape(V4, V4)
    # 只约束中：中 = 底的 parity 聚合
    for j in range(V4):
        for s in range(2):
            m = Wb[j, 2*s] + Wb[j, 2*s+1]
            Wm[j, 2*s] = Wm[j, 2*s+1] = 0.5 * m
    # seed = 最近两代闭合（top2）
    if not hasattr(cprime_update, 'layers'): cprime_update.layers = []
    cprime_update.layers.append(clo.copy())
    if len(cprime_update.layers) > 2: cprime_update.layers.pop(0)
    seed = sum(cprime_update.layers)
    if seed.sum() > 0 and np.all(np.isfinite(seed)):
        seed = seed / seed.sum()
    else:
        seed = np.ones(V4) / V4
    for W in (Wb, Wm, Wt):
        if not np.all(np.isfinite(W)):
            W[:] = 0.0
    return Wb, Wm, Wt, seed


def run(coupled, steps=500):
    adj = zs.cycle(16)
    e = ZeroEvo(adj, Params(L=4, T=4, micro=True, c_h=1.0, dt=0.05),
                seed=np.ones(16))
    if coupled:
        Wb = np.zeros((V4, V4)); Wm = np.zeros((V4, V4)); Wt = np.zeros((V4, V4))
        seed = np.ones(V4) / V4
    rows = []
    for g in range(1, steps+1):
        w = None
        if coupled:
            Wb, Wm, Wt, seed = cprime_update(Wb, Wm, Wt, seed, g)
            mb, mm, mt = np.abs(Wb).sum(), np.abs(Wm).sum(), np.abs(Wt).sum()
            ms = mb + mm + mt
            w = [float(mb/ms), float(mm/ms), float(mt/ms)] if ms > 0 else [1/3, 1/3, 1/3]
            p = np.abs(Wb).reshape(-1)
            if p.sum() > 0 and np.all(np.isfinite(p)):
                p = p / p.sum()
                e.seed = float(e.seed.sum()) * p
            else:
                e.seed = np.ones(16)
        e.step()
        s = e.snapshot()
        rows.append({"gen": g, "q": s["q"], "n_std": s["n_std"],
                     "q_live": s["q_live"], "q_hist": s["q_hist"],
                     "psi_coh": s["psi_coh"], "w": w})
    qs = np.array([r["q"] for r in rows])
    ql = np.array([r["q_live"] for r in rows])
    return rows, qs, ql


def main():
    rb, qb, qlb = run(False)
    rc, qc, qlc = run(True)
    print("baseline q_live tail", np.round(qlb[-100:].mean(), 6))
    print("coupled  q_live tail", np.round(qlc[-100:].mean(), 6))
    print("baseline q tail", np.round(qb[-100:].mean(), 6), "std", np.round(qb[-100:].std(), 6))
    print("coupled  q tail", np.round(qc[-100:].mean(), 6), "std", np.round(qc[-100:].std(), 6))
    ck("baseline q ≈ 1/16（真实 content）", abs(qb[-100:].mean() - 1/16) < 0.01,
       f"{qb[-100:].mean():.6f}")
    ck("coupled q 有定义", np.all(np.isfinite(qc)), "")
    ck("coupled q_live 接近观测 0.0683（|Δ|<0.01）",
       abs(qlc[-100:].mean() - 0.06830296233) < 0.01,
       f"coupled_q_live={qlc[-100:].mean():.6f} base={qlb[-100:].mean():.6f}")
    # 周期 4：FFT bin = steps/T
    def p4(x):
        x = np.asarray(x, dtype=float); x = x - x.mean()
        f = np.abs(np.fft.rfft(x)); return float(f[len(x)//T]) if len(f) > len(x)//T else 0.0
    qlh = np.array([r["q_live"] for r in rc])
    qhh = np.array([r["q_hist"] for r in rc])
    wc = np.array([r["w"] for r in rc if r["w"] is not None])
    p4_ql, p4_qh = p4(qlh), p4(qhh)
    p4_wt = p4(wc[:, 2]) if len(wc) else 0.0
    ck("q_live 有周期 4 分量", p4_ql > 1e-6, f"bin={p4_ql:.3g}")
    ck("q_hist 有周期 4 分量", p4_qh > 1e-6, f"bin={p4_qh:.3g}")
    ck("C′ 顶层份额有周期 4 分量", p4_wt > 1e-6, f"bin={p4_wt:.3g}")
    if len(wc):
        print("C′ w mean", np.round(wc.mean(0), 4), "dist", round(float(np.abs(wc.mean(0)-np.array([0.1724,0.2759,0.5517])).sum()),4))
    g = np.arange(1, 501)
    fig, ax = plt.subplots(figsize=(6, 3.4))
    ax.plot(g, qb, label="baseline", lw=1)
    ax.plot(g, qc, label="C′ seed q", lw=1)
    ax.plot(g, qlb, label="baseline q_live", lw=1)
    ax.plot(g, qlc, label="C′ q_live", lw=1)
    ax.axhline(0.06830296233, color="r", ls=":", lw=1, label="obs")
    ax.set_title("real content q (cycle16)"); ax.legend(fontsize=8); ax.grid(alpha=.3)
    fig.tight_layout(); fp = os.path.join(RES, "zero_cprime_12ch.png")
    fig.savefig(fp, dpi=130)
    out = {"baseline_q": float(qb[-100:].mean()), "coupled_q": float(qc[-100:].mean()),
           "p4_q_live": p4_ql, "p4_q_hist": p4_qh, "p4_w_top": p4_wt,
           "w_mean": wc.mean(0).tolist() if len(wc) else None,
           "baseline_q_live": float(qlb[-100:].mean()), "coupled_q_live": float(qlc[-100:].mean()),
           "rows_baseline": rb, "rows_coupled": rc,
           "checks": CK, "n_fail": sum(1 for c in CK if not c["pass"]),
           "anchor": "C′ Wb -> cycle16 seed; real content q; no dimensionful anchor"}
    with open(os.path.join(RES, "zero_cprime_12ch.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_cprime_12ch.json"), "和", fp)


if __name__ == "__main__":
    main()
