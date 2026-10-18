#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""存在性测试：顶拷贝的初始偏置能否把份额推到观测？

C′ 记录版（λ=1，三份 W 接 16 槽记录）；只改初始 Wt 的对角：
    Wt[j,j] = top0
扫 top0，跑 500 代，看尾段份额/q/周期 4。
输出：results/zero_cprime_top0.json/.png
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


def rho():
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


RHO, RHO_B, RHO_M = rho()


def cprime_step(Wb, Wm, Wt, seed, gen):
    A = np.zeros((V4, V4))
    for u, nb in enumerate(zs.k4()):
        for v in nb: A[u, v] = 1.0
    AL = np.linalg.matrix_power(A, L)
    diag = np.arange(V4)
    Wb[diag, diag] += (1/3) * seed
    Wm[diag, diag] += (1/3) * seed
    Wt[diag, diag] += (1/3) * seed
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
    Wm = Wm_t
    seed = sum([clo])
    if seed.sum() > 0 and np.all(np.isfinite(seed)): seed = seed / seed.sum()
    else: seed = np.ones(V4) / V4
    return Wb, Wm, Wt, seed


def record16(Wb, Wm, Wt):
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


def run(top0, steps=500):
    e = ZeroEvo(zs.cycle(16), Params(L=4, T=4, micro=True, c_h=1.0, dt=0.05),
                seed=np.ones(16))
    Wb = np.zeros((V4, V4)); Wm = np.zeros((V4, V4)); Wt = np.zeros((V4, V4))
    Wt[np.arange(V4), np.arange(V4)] = top0
    seed = np.ones(V4) / V4
    rows = []
    for g in range(1, steps+1):
        Wb, Wm, Wt, seed = cprime_step(Wb, Wm, Wt, seed, g)
        e.step()
        x = record16(Wb, Wm, Wt)
        if x.sum() > 0 and np.all(np.isfinite(x)):
            e.m_live = x / x.sum(); e.m_hist = x / x.sum()
        s = e.snapshot()
        mb, mm, mt = np.abs(Wb).sum(), np.abs(Wm).sum(), np.abs(Wt).sum()
        ms = mb + mm + mt
        w = [float(mb/ms), float(mm/ms), float(mt/ms)] if ms > 0 else [1/3]*3
        rows.append({"gen": g, "q_live": s["q_live"], "w": w})
    ql = np.array([r["q_live"] for r in rows]); w = np.array([r["w"] for r in rows])
    def p4(a):
        a = np.asarray(a) - np.mean(a); f = np.abs(np.fft.rfft(a)); return float(f[len(a)//T])
    return {"top0": top0, "q": float(ql[-100:].mean()),
            "w": w[-100:].mean(0).tolist(),
            "dist": float(np.abs(w[-100:].mean(0)-TARGET).sum()),
            "p4q": p4(ql), "p4w": p4(w[:, 2]), "rows": rows}


def main():
    top0s = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0]
    res = [run(t) for t in top0s]
    print("top0  q_live  w_b     w_m     w_t     dist    p4q   p4w")
    for r in res:
        print(f"{r['top0']:.2f}  {r['q']:.4f}  " + "  ".join(f"{x:.4f}" for x in r["w"]) +
              f"  {r['dist']:.4f}  {r['p4q']:.3f}  {r['p4w']:.1f}")
    ck("存在 top0 使份额距离 <0.1", any(r["dist"] < 0.1 for r in res),
       str([round(r["dist"],3) for r in res]))
    ck("存在 top0 使 q_live 接近观测", any(abs(r["q"]-TARGET_Q) < 0.01 for r in res),
       str([round(r["q"],4) for r in res]))
    ck("存在 top0 使 q_live 有周期 4", any(r["p4q"] > 1e-6 for r in res),
       str([round(r["p4q"],3) for r in res]))
    ck("无 NaN", all(all(np.isfinite(r["w"])) for r in res), "")

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    ax[0].plot(top0s, [r["dist"] for r in res], marker="o", label="share dist")
    ax[0].plot(top0s, [abs(r["q"]-TARGET_Q) for r in res], marker="o", label="|q-obs|")
    ax[0].set_title("(a) trade-off vs top0"); ax[0].legend(fontsize=8)
    for i, lab in enumerate(["w_b", "w_m", "w_t"]):
        ax[1].plot(top0s, [r["w"][i] for r in res], marker="o", label=lab)
    ax[1].set_title("(b) shares vs top0"); ax[1].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); fp = os.path.join(RES, "zero_cprime_top0.png")
    fig.savefig(fp, dpi=130)
    out = {"results": [{k: v for k, v in r.items() if k != "rows"} for r in res],
           "rows": {str(r["top0"]): r["rows"] for r in res},
           "checks": CK, "n_fail": sum(1 for c in CK if not c["pass"]),
           "anchor": "top0 initial bias; record coupling; no dimensionful anchor"}
    with open(os.path.join(RES, "zero_cprime_top0.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_cprime_top0.json"), "和", fp)


if __name__ == "__main__":
    main()
