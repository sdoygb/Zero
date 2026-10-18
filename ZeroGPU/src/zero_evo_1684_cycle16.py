#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""16/8/4 固定底层接 cycle16：class 空间 4x4(16) 映射到 16 顶点。

macro/geometry/wave 在 cycle16；micro 仍是三份 W（底16/中8/顶4），
x=(Wb+A_m^T Wm+A_t^T Wt)/7 作为 16 槽记录，alpha 回写 content。
只问能不能跑、周期 4 是否还在。
输出：results/zero_evo_1684_cycle16.json/.png
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
V4, T = 4, 4
TARGET_Q = 0.06830296233
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


class ZeroEvo1684Cycle16(ZeroEvo):
    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.A4 = zs.adj_matrix(zs.k4())
        self.A4_L = np.linalg.matrix_power(self.A4, self.p.L)
        self.Wb = np.zeros((V4, V4)); self.Wm = np.zeros((V4, V4)); self.Wt = np.zeros((V4, V4))
        self.sub_seed = np.ones(V4) / V4

    def _rec16(self):
        Wm8 = np.zeros((V4, 2)); Wt4 = self.Wt.sum(axis=1)
        for j in range(V4):
            for s in range(2): Wm8[j, s] = self.Wm[j, 2*s] + self.Wm[j, 2*s+1]
        x = self.Wb.reshape(-1).copy()
        for j in range(V4):
            for s in range(2):
                x[4*j+2*s] += Wm8[j, s]
                x[4*j+2*s+1] += Wm8[j, s]
            x[4*j:4*j+4] += Wt4[j]
        return x / 7.0

    def _micro_gen(self):
        diag = np.arange(V4)
        self.Wb[diag, diag] += (1/3) * self.sub_seed
        self.Wm[diag, diag] += (1/3) * self.sub_seed
        self.Wt[diag, diag] += (1/3) * self.sub_seed
        self.Wb = self.Wb @ self.A4_L
        self.Wm = self.Wm @ self.A4_L
        self.Wt = self.Wt @ self.A4_L
        for Wc in (self.Wb, self.Wm, self.Wt):
            if (not np.all(np.isfinite(Wc))) or np.abs(Wc).sum() > 1e200:
                Wc[:] = 0.0
        self.gen += 1
        if self.p.T and self.gen % self.p.T == 0:
            clo = np.diag(self.Wb) + np.diag(self.Wm) + np.diag(self.Wt)
            jstar = int(np.argmax(clo))
            self.Wb[:] = 0.0
            for j in range(V4):
                if j != jstar: self.Wm[j, :] = 0.0
            k = 0.2
            Pt = np.abs(self.Wt).sum()
            if np.isfinite(Pt) and Pt > 0:
                J = k * Pt; self.Wt -= self.Wt * (J / Pt)
                for j in range(V4):
                    for s in range(2):
                        self.Wm[j, 2*s] += J / (2*V4)
                        self.Wm[j, 2*s+1] += J / (2*V4)
            Pm = np.abs(self.Wm).sum()
            if np.isfinite(Pm) and Pm > 0:
                J2 = k * Pm; self.Wm -= self.Wm * (J2 / Pm)
                self.Wb += J2 / (V4*V4)
        mw = max(np.abs(self.Wb).sum(), np.abs(self.Wm).sum(), np.abs(self.Wt).sum())
        if np.isfinite(mw) and mw > 0:
            self.Wb /= mw; self.Wm /= mw; self.Wt /= mw
        x = self._rec16()
        if not np.all(np.isfinite(x)):
            x = np.ones(16) / 16
        self.gen_closures = x
        self.record = self.record + x
        clo = np.diag(self.Wb) + np.diag(self.Wm) + np.diag(self.Wt)
        if np.all(np.isfinite(clo)) and clo.sum() > 0:
            self.sub_seed = clo / clo.sum()
        else:
            self.sub_seed = np.ones(V4) / V4

    def q16(self):
        p = np.abs(self._rec16()).reshape(-1)
        s = p.sum()
        return float(((p / s) ** 2).sum()) if s > 0 else 0.0

    def shares(self):
        m = np.array([np.abs(self.Wb).sum(), np.abs(self.Wm).sum(), np.abs(self.Wt).sum()])
        s = m.sum()
        return (m / s).tolist() if s > 0 else [1/3]*3


def run(alpha, steps=1000):
    e = ZeroEvo1684Cycle16(zs.cycle(16), Params(L=4, T=4, micro=True, c_h=1.0, dt=0.05, alpha=alpha),
                           seed=np.ones(16))
    rows = []
    for g in range(1, steps+1):
        e.step()
        s = e.snapshot()
        rows.append({"gen": g, "q": s["q"], "q_live": s["q_live"], "n_std": s["n_std"],
                     "q16": e.q16(), "w": e.shares()})
    q = np.array([r["q"] for r in rows]); q16 = np.array([r["q16"] for r in rows])
    nstd = np.array([r["n_std"] for r in rows]); w = np.array([r["w"] for r in rows])
    def p4(a):
        a = np.asarray(a) - np.mean(a); f = np.abs(np.fft.rfft(a)); return float(f[len(a)//T])
    return {"alpha": alpha, "q": float(q[-100:].mean()), "q16": float(q16[-100:].mean()),
            "n_std": float(nstd[-100:].mean()), "n_std_max": float(np.nanmax(nstd)),
            "w": w[-100:].mean(0).tolist(), "p4q": p4(q), "p4q16": p4(q16), "rows": rows}


def main():
    alphas = [0.0, 0.1, 0.3]
    res = [run(a) for a in alphas]
    print("alpha  q       q16     n_std   shares                 p4q    p4q16")
    for r in res:
        print(f"{r['alpha']:.2f}  {r['q']:.4f}  {r['q16']:.4f}  {r['n_std']:.4f}  "
              + str(np.round(r["w"], 4)) + f"  {r['p4q']:.3f}  {r['p4q16']:.2f}")
    ck("无 NaN", all(np.isfinite(r["q"]) and np.isfinite(r["q16"]) for r in res), "")
    ck("content 有界（n_std<1）", all(r["n_std_max"] < 1.0 for r in res),
       str([round(r["n_std_max"],3) for r in res]))
    ck("q16 有定义", all(r["q16"] > 0 for r in res), str([round(r["q16"],4) for r in res]))
    ck("份额归一", all(abs(sum(r["w"])-1) < 1e-9 for r in res), "")
    ck("alpha>0 时 content 有周期 4", any(r["p4q"] > 1e-6 for r in res if r["alpha"] > 0),
       str([round(r["p4q"],3) for r in res]))
    ck("q16 有周期 4", any(r["p4q16"] > 1e-6 for r in res),
       str([round(r["p4q16"],2) for r in res]))

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    ax[0].plot(alphas, [r["q"] for r in res], marker="o", label="content q")
    ax[0].plot(alphas, [r["q16"] for r in res], marker="o", label="q16")
    ax[0].set_title("(a) q vs alpha"); ax[0].legend(fontsize=8)
    for i, lab in enumerate(["Wb", "Wm", "Wt"]):
        ax[1].plot(alphas, [r["w"][i] for r in res], marker="o", label=lab)
    ax[1].set_title("(b) shares vs alpha"); ax[1].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); fp = os.path.join(RES, "zero_evo_1684_cycle16.png")
    fig.savefig(fp, dpi=130)
    out = {"results": [{k: v for k, v in r.items() if k != "rows"} for r in res],
           "rows": {str(r["alpha"]): r["rows"] for r in res},
           "checks": CK, "n_fail": sum(1 for c in CK if not c["pass"]),
           "anchor": "16/8/4 fixed; cycle16 macro; alpha scan"}
    with open(os.path.join(RES, "zero_evo_1684_cycle16.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_evo_1684_cycle16.json"), "和", fp)


if __name__ == "__main__":
    main()
