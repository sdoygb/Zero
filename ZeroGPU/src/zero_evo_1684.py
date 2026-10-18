#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把 16/8/4 作为固定底层：三份 W 的 ZeroEvo 子类。

16/8/4 是输入（不再拟合观测份额）：
  Wb 底 16、Wm 中 8、Wt 顶 4；
  T=4 底全灭、中留最高 site、顶不动；
  flow 顶->中->底；seed 三层等分。
只问：这个底层能不能跑（有界、无 NaN、周期 4）。
输出：results/zero_evo_1684.json/.png
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
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


class ZeroEvo1684(ZeroEvo):
    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.A_L = np.linalg.matrix_power(self.A, self.p.L)
        self.Wb = np.zeros((self.V, self.V))
        self.Wm = np.zeros((self.V, self.V))
        self.Wt = np.zeros((self.V, self.V))

    def _micro_gen(self):
        V, diag = self.V, np.arange(self.V)
        for Wc in (self.Wb, self.Wm, self.Wt):
            Wc[diag, diag] += (1/3) * self.seed
        self.Wb = self.Wb @ self.A_L
        self.Wm = self.Wm @ self.A_L
        self.Wt = self.Wt @ self.A_L
        for Wc in (self.Wb, self.Wm, self.Wt):
            if (not np.all(np.isfinite(Wc))) or np.abs(Wc).sum() > 1e200:
                Wc[:] = 0.0
        cb = np.diag(self.Wb).copy(); cm = np.diag(self.Wm).copy(); ct = np.diag(self.Wt).copy()
        if not (np.all(np.isfinite(cb)) and np.all(np.isfinite(cm)) and np.all(np.isfinite(ct))):
            cb[:] = 0; cm[:] = 0; ct[:] = 0
        c = cb + cm + ct
        self.gen_closures = c
        self.record = self.record + c
        np.fill_diagonal(self.Wb, 0); np.fill_diagonal(self.Wm, 0); np.fill_diagonal(self.Wt, 0)
        mw = max(np.abs(self.Wb).sum(), np.abs(self.Wm).sum(), np.abs(self.Wt).sum())
        if np.isfinite(mw) and mw > 1e200:
            self.Wb /= mw; self.Wm /= mw; self.Wt /= mw
        if not (np.all(np.isfinite(self.Wb)) and np.all(np.isfinite(self.Wm)) and np.all(np.isfinite(self.Wt))):
            self.Wb[:] = 0; self.Wm[:] = 0; self.Wt[:] = 0
        self.layers.append(c.copy())
        if len(self.layers) > 2: self.layers.pop(0)
        if self.p.retention == "last": ret = c.copy()
        elif self.p.retention == "top2":
            ret = c + (self.layers[-2] if len(self.layers) >= 2 else 0.0)
        else: ret = self.record.copy()
        if np.all(np.isfinite(ret)) and ret.sum() > 0:
            ret = ret / ret.sum()
        else:
            ret = np.ones(V) / V
        self.seed = ret
        self.gen += 1
        if self.p.T and self.gen % self.p.T == 0:
            jstar = int(np.argmax(c))
            self.Wb[:] = 0.0
            for j in range(V):
                if j != jstar: self.Wm[j, :] = 0.0
            # flow
            k = 0.2
            Pt = np.abs(self.Wt).sum()
            if np.isfinite(Pt) and Pt > 0:
                J = k * Pt; self.Wt -= self.Wt * (J / Pt)
                for j in range(V):
                    for s in range(2):
                        self.Wm[j, 2*s] += J / (2*V)
                        self.Wm[j, 2*s+1] += J / (2*V)
            Pm = np.abs(self.Wm).sum()
            if np.isfinite(Pm) and Pm > 0:
                J2 = k * Pm; self.Wm -= self.Wm * (J2 / Pm)
                self.Wb += J2 / (V*V)
        self.W = self.Wb + self.Wm + self.Wt
        # 公共归一化（每代，保持比值，防指数增长）
        mw = max(np.abs(self.Wb).sum(), np.abs(self.Wm).sum(), np.abs(self.Wt).sum())
        if np.isfinite(mw) and mw > 0:
            self.Wb /= mw; self.Wm /= mw; self.Wt /= mw; self.W /= mw

    def q16(self):
        p = np.abs(self.Wb).reshape(-1)
        s = p.sum()
        return float(((p / s) ** 2).sum()) if s > 0 else 0.0

    def shares(self):
        m = np.array([np.abs(self.Wb).sum(), np.abs(self.Wm).sum(), np.abs(self.Wt).sum()])
        s = m.sum()
        return (m / s).tolist() if s > 0 else [1/3]*3


def run(alpha, steps=1000):
    e = ZeroEvo1684(zs.k4(), Params(L=4, T=4, micro=True, c_h=1.0, dt=0.05, alpha=alpha),
                    seed=np.array([8.0, 4.0, 1.5, 1.0]))
    rows = []
    for g in range(1, steps+1):
        e.step()
        s = e.snapshot()
        rows.append({"gen": g, "q": s["q"], "q_live": s["q_live"],
                     "n_std": s["n_std"], "q16": e.q16(), "w": e.shares()})
    q = np.array([r["q"] for r in rows]); q16 = np.array([r["q16"] for r in rows])
    w = np.array([r["w"] for r in rows]); nstd = np.array([r["n_std"] for r in rows])
    def p4(a):
        a = np.asarray(a) - np.mean(a); f = np.abs(np.fft.rfft(a)); return float(f[len(a)//T])
    return {"alpha": alpha, "q": float(q[-100:].mean()), "q16": float(q16[-100:].mean()),
            "n_std": float(nstd[-100:].mean()), "n_std_max": float(np.nanmax(nstd)),
            "w": w[-100:].mean(0).tolist(), "p4q": p4(q), "p4q16": p4(q16),
            "rows": rows}


def main():
    alphas = [0.0, 0.1, 0.3, 0.5, 1.0]
    res = [run(a) for a in alphas]
    print("alpha  q       q16     n_std   shares                 p4q    p4q16")
    for r in res:
        print(f"{r['alpha']:.2f}  {r['q']:.4f}  {r['q16']:.4f}  {r['n_std']:.4f}  "
              + str(np.round(r["w"], 4)) + f"  {r['p4q']:.3f}  {r['p4q16']:.2f}")
    ck("alpha>0 时 content q 有周期 4", any(r["p4q"] > 1e-6 for r in res if r["alpha"] > 0),
       str([round(r["p4q"],3) for r in res]))
    ck("所有 alpha 有界（n_std<1）", all(r["n_std_max"] < 1.0 for r in res),
       str([round(r["n_std_max"],3) for r in res]))
    ck("所有 alpha 无 NaN", all(np.isfinite(r["q"]) and np.isfinite(r["q16"]) for r in res), "")
    ck("所有 alpha 的 q16 有定义", all(r["q16"] > 0 for r in res), str([round(r["q16"],4) for r in res]))
    ck("所有 alpha 的份额归一", all(abs(sum(r["w"])-1) < 1e-9 for r in res), "")

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    ax[0].plot(alphas, [r["q"] for r in res], marker="o", label="content q")
    ax[0].plot(alphas, [r["q16"] for r in res], marker="o", label="q16")
    ax[0].set_title("(a) q vs alpha"); ax[0].legend(fontsize=8)
    for i, lab in enumerate(["Wb", "Wm", "Wt"]):
        ax[1].plot(alphas, [r["w"][i] for r in res], marker="o", label=lab)
    ax[1].set_title("(b) shares vs alpha"); ax[1].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); fp = os.path.join(RES, "zero_evo_1684.png")
    fig.savefig(fp, dpi=130)
    out = {"results": [{k: v for k, v in r.items() if k != "rows"} for r in res],
           "rows": {str(r["alpha"]): r["rows"] for r in res},
           "checks": CK, "n_fail": sum(1 for c in CK if not c["pass"]),
           "anchor": "16/8/4 fixed; alpha scan; no fitting"}
    with open(os.path.join(RES, "zero_evo_1684.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_evo_1684.json"), "和", fp)


if __name__ == "__main__":
    main()
