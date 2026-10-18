#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""C′ 接 12 通道（记录版）：三份 W 作为 16 槽记录接进 m_live/m_hist，不接 seed。

C′ 在 K4 上独立演化；每代把三份 W 按“求和聚合”的最小二乘读法
    x = (Wb + A_m^T Wm8 + A_t^T Wt4) / 7
（A_m^T A_m=2I, A_t^T A_t=4I，正规矩阵 7I）变成 16 槽记录，
直接写入引擎的 m_live/m_hist；seed 保持默认 1。
输出：results/zero_cprime_record.json/.png
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
TARGET = np.array([0.1724, 0.2759, 0.5517])
TARGET_Q = 0.06830296233
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


def window_rho():
    R = np.zeros((V4, V4))
    for w in zs.balanced_words(8):
        for i in range(5):
            sub = w[i:i+4]
            if sum(sub) != 0: continue
            par = zs.word_area(sub) % 2
            for j in range(V4):
                R[j, 2*(0 if sub[j] > 0 else 1) + par] += 1
    R /= R.sum()
    rb = R.reshape(-1); rb /= rb.sum()
    rm = np.zeros((V4, 2))
    for j in range(V4):
        for s in range(2): rm[j, s] = R[j, 2*s] + R[j, 2*s+1]
    rm /= rm.sum()
    return R, rb, rm


RHO, RHO_B, RHO_M = window_rho()


def cprime_step(Wb, Wm, Wt, seed, gen, lam=1.0):
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
        Wb[:] = 0; Wm[:] = 0; Wt[:] = 0
    clo = np.diag(Wb) + np.diag(Wm) + np.diag(Wt)
    np.fill_diagonal(Wb, 0); np.fill_diagonal(Wm, 0); np.fill_diagonal(Wt, 0)
    if gen % T == 0:
        jstar = int(np.argmax(clo))
        Wb[:] = 0.0
        for j in range(V4):
            if j != jstar: Wm[j, :] = 0.0
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
    Wm_t = np.zeros((V4, V4))
    for j in range(V4):
        for s in range(2):
            m = Wb[j, 2*s] + Wb[j, 2*s+1]
            Wm_t[j, 2*s] = Wm_t[j, 2*s+1] = 0.5 * m
    Wm = (1-lam) * Wm + lam * Wm_t
    seed = sum([clo])
    if seed.sum() > 0 and np.all(np.isfinite(seed)): seed = seed / seed.sum()
    else: seed = np.ones(V4) / V4
    return Wb, Wm, Wt, seed


def record16(Wb, Wm, Wt):
    """x=(Wb+A_m^T Wm8+A_t^T Wt4)/7，求和聚合的最小二乘。"""
    Wm8 = np.zeros((V4, 2)); Wt4 = Wt.sum(axis=1)
    for j in range(V4):
        for s in range(2): Wm8[j, s] = Wm[j, 2*s] + Wm[j, 2*s+1]
    x = Wb.reshape(-1).copy()
    for j in range(V4):
        for s in range(2):
            x[4*j+2*s] += Wm8[j, s]
            x[4*j+2*s+1] += Wm8[j, s]
        x[4*j:4*j+4] += Wt4[j]
    return x / 7.0


def run(lam, steps=500):
    e = ZeroEvo(zs.cycle(16), Params(L=4, T=4, micro=True, c_h=1.0, dt=0.05),
                seed=np.ones(16))
    Wb = np.zeros((V4, V4)); Wm = np.zeros((V4, V4)); Wt = np.zeros((V4, V4))
    seed = np.ones(V4) / V4
    rows = []
    for g in range(1, steps+1):
        Wb, Wm, Wt, seed = cprime_step(Wb, Wm, Wt, seed, g, lam)
        e.step()
        x = record16(Wb, Wm, Wt)
        if x.sum() > 0 and np.all(np.isfinite(x)):
            e.m_live = x / x.sum(); e.m_hist = x / x.sum()
        s = e.snapshot()
        mb, mm, mt = np.abs(Wb).sum(), np.abs(Wm).sum(), np.abs(Wt).sum()
        ms = mb + mm + mt
        w = [float(mb/ms), float(mm/ms), float(mt/ms)] if ms > 0 else [1/3]*3
        rows.append({"gen": g, "q_live": s["q_live"], "q_hist": s["q_hist"],
                     "q_content": s["q"], "w": w})
    ql = np.array([r["q_live"] for r in rows]); w = np.array([r["w"] for r in rows])
    def p4(a):
        a = np.asarray(a) - np.mean(a); f = np.abs(np.fft.rfft(a)); return float(f[len(a)//T])
    return {"lam": lam, "rows": rows, "q": float(ql[-100:].mean()),
            "w": w[-100:].mean(0).tolist(),
            "dist": float(np.abs(w[-100:].mean(0)-TARGET).sum()),
            "p4q": p4(ql), "p4w": p4(w[:, 2])}


def main():
    lams = [0.5, 0.75, 0.9, 1.0]
    res = [run(l) for l in lams]
    print("lam   q_live   w_b     w_m     w_t     dist    p4q     p4w")
    for r in res:
        print(f"{r['lam']:.2f}  {r['q']:.4f}  " + "  ".join(f"{x:.4f}" for x in r["w"]) +
              f"  {r['dist']:.4f}  {r['p4q']:.3f}  {r['p4w']:.1f}")
    best = min(res, key=lambda r: r["dist"])
    print("best dist", best["lam"], best["dist"])
    ck("存在 λ 使 q_live 接近观测", any(abs(r["q"]-TARGET_Q) < 0.01 for r in res),
       str([round(r["q"],4) for r in res]))
    ck("存在 λ 使份额距离 <0.1", any(r["dist"] < 0.1 for r in res),
       str([round(r["dist"],4) for r in res]))
    ck("存在 λ 使 q_live 有周期 4", any(r["p4q"] > 1e-6 for r in res),
       str([round(r["p4q"],3) for r in res]))
    ck("存在 λ 使 w_top 有周期 4", any(r["p4w"] > 1e-6 for r in res),
       str([round(r["p4w"],1) for r in res]))

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    ax[0].plot(lams, [r["dist"] for r in res], marker="o", label="share dist")
    ax[0].plot(lams, [abs(r["q"]-TARGET_Q) for r in res], marker="o", label="|q-obs|")
    ax[0].set_title("(a) trade-off vs λ"); ax[0].legend(fontsize=8)
    for i, lab in enumerate(["w_b", "w_m", "w_t"]):
        ax[1].plot(lams, [r["w"][i] for r in res], marker="o", label=lab)
    ax[1].set_title("(b) shares vs λ"); ax[1].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); fp = os.path.join(RES, "zero_cprime_record.png")
    fig.savefig(fp, dpi=130)
    out = {"results": res, "best": best, "checks": CK,
           "n_fail": sum(1 for c in CK if not c["pass"]),
           "anchor": "three W as record; λ scan; seed untouched"}
    with open(os.path.join(RES, "zero_cprime_record.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_cprime_record.json"), "和", fp)


if __name__ == "__main__":
    main()
