#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""16/8/4 + cycle16：模依赖共形红移相位（精确化 beta_h=omega/2）。

相位的红移项用算子 omega_op = V diag(omega_k) V^T：
    G = (omega_op H + H omega_op)/2,  H=diag(h-hbar)
    psi <- exp(-i dt G) psi
与标量 beta_h=omega_char/2 对比。
输出：results/zero_evo_1684_mode.json/.png
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
from zero_evo import Params  # noqa: E402
from zero_evo_1684_cycle16 import ZeroEvo1684Cycle16, MAPS  # noqa: E402

RES = os.path.join(HERE, "..", "results")
T = 4
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


class ZeroEvo1684Mode(ZeroEvo1684Cycle16):
    def __init__(self, *a, order="sym", **kw):
        super().__init__(*a, **kw)
        self.order = order
        self.omega = np.sqrt(np.maximum(self.wave_eigs, 0.0))
        self.omega_op = self.wave_vecs @ np.diag(self.omega) @ self.wave_vecs.T

    def _apply_mode_phase(self, dt):
        H = np.diag(self.h - self.h.mean())
        if self.order == "sym":
            G = 0.5 * (self.omega_op @ H + H @ self.omega_op)
        elif self.order == "left":
            A = self.omega_op @ H
            G = 0.5 * (A + A.conj().T)
        elif self.order == "right":
            A = H @ self.omega_op
            G = 0.5 * (A + A.conj().T)
        else:
            raise ValueError(self.order)
        vals, vecs = np.linalg.eigh(G)
        U = vecs @ np.diag(np.exp(-1j * dt * vals)) @ vecs.conj().T
        self.psi = U @ self.psi

    def _wave_step(self, dt):
        if not self.p.wave_on:
            return
        self._apply_mode_phase(0.5 * dt)
        self.psi = self.wave_vecs @ (np.exp(-1j * self.wave_eigs * dt) *
                                     (self.wave_vecs.T @ self.psi))
        self._apply_mode_phase(0.5 * dt)
        if self.p.Gamma_w:
            import math
            self.psi = self.psi * math.exp(-self.p.Gamma_w * dt)


def omega_char():
    e = ZeroEvo1684Cycle16(zs.cycle(16), Params(L=4, T=4, micro=True, c_h=1.0, dt=0.05,
                                                  alpha=0.1, g_src=1.0, beta=0.0, beta_h=0.0),
                            seed=np.ones(16))
    e.perm = MAPS["seq"]
    for _ in range(2000): e.step()
    c = e.wave_vecs.T @ e.psi; P = np.abs(c) ** 2
    w = np.sqrt(np.maximum(e.wave_eigs, 0))
    return float((P * w).sum() / P.sum())


def run(mode="scalar", bh=0.0, steps=2000, order="sym"):
    if mode == "scalar":
        e = ZeroEvo1684Cycle16(zs.cycle(16), Params(L=4, T=4, micro=True, c_h=1.0, dt=0.05,
                                                      alpha=0.1, g_src=1.0, beta=0.0, beta_h=bh),
                                seed=np.ones(16))
    else:
        e = ZeroEvo1684Mode(zs.cycle(16), Params(L=4, T=4, micro=True, c_h=1.0, dt=0.05,
                                                   alpha=0.1, g_src=1.0, beta=0.0, beta_h=0.0),
                             seed=np.ones(16), order=order)
    e.perm = MAPS["seq"]
    rows = []
    for g in range(1, steps+1):
        e.step(); s = e.snapshot()
        rows.append({"gen": g, "q": s["q"], "q16": e.q16(), "h_E": s["h_energy"],
                     "psi_coh": s["psi_coh"], "n_std": s["n_std"]})
    return rows


def p4(a):
    a = np.asarray(a) - np.mean(a); f = np.abs(np.fft.rfft(a)); return float(f[len(a)//T])


def main():
    wc = omega_char(); bh = wc / 2
    print("omega_char", wc, "scalar beta_h", bh)
    rs = run("scalar", bh, 2000)
    orders = ["sym", "left", "right"]
    ro = {o: run("mode", 0.0, 2000, o) for o in orders}
    def tail(r, k): return float(np.mean([x[k] for x in r[-200:]]))
    print("coupling   q       q16     h_E       psi_coh  p4psi")
    print(f"scalar     {tail(rs,'q'):.4f}  {tail(rs,'q16'):.4f}  {tail(rs,'h_E'):.3f}  "
          f"{tail(rs,'psi_coh'):.4f}  {p4([x['psi_coh'] for x in rs]):.4f}")
    for o in orders:
        r = ro[o]
        print(f"mode-{o:5s} {tail(r,'q'):.4f}  {tail(r,'q16'):.4f}  {tail(r,'h_E'):.3f}  "
              f"{tail(r,'psi_coh'):.4f}  {p4([x['psi_coh'] for x in r]):.4f}")
    ck("三种顺序都无 NaN", all(np.all(np.isfinite([x['psi_coh'] for x in ro[o]])) for o in orders), "")
    ck("三种顺序 psi_coh 有界", all(max(x['psi_coh'] for x in ro[o]) <= 1.0 + 1e-9 for o in orders), "")
    ck("三种顺序都有周期4", all(p4([x['psi_coh'] for x in ro[o]]) > 1e-6 for o in orders),
       str([round(p4([x['psi_coh'] for x in ro[o]]),3) for o in orders]))
    ck("三种顺序 h_E 相近", max(tail(ro[o],'h_E') for o in orders) - min(tail(ro[o],'h_E') for o in orders) < 1e-6, "")

    print("长跑 mode-sym, 10000 代 ...")
    rl = run("mode", 0.0, 10000, "sym")
    psi = np.array([x['psi_coh'] for x in rl])
    print("long psi_coh", round(float(psi[-200:].mean()),4), "p4psi", round(p4(psi),4),
          "p4h", round(p4([x['h_E'] for x in rl]),3), "q16", round(float(np.mean([x['q16'] for x in rl[-200:]])),4))
    ck("长跑 psi 有周期4", p4(psi) > 1e-6, f"{p4(psi):.4f}")
    ck("长跑 psi_coh 有界", float(np.max(psi)) <= 1.0 + 1e-9, f"max={np.max(psi):.4f}")
    ck("长跑无 NaN", np.all(np.isfinite(psi)), "")

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    ax[0].plot([x['psi_coh'] for x in rs], lw=1, label="scalar ω/2")
    for o in orders:
        ax[0].plot([x['psi_coh'] for x in ro[o]], lw=1, label=f"mode-{o}")
    ax[0].set_title("(a) psi_coh (2000)"); ax[0].legend(fontsize=7)
    ax[1].plot(psi, lw=1, label="mode-sym long")
    ax[1].set_title("(b) psi_coh (10000)"); ax[1].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); fp = os.path.join(RES, "zero_evo_1684_mode.png")
    fig.savefig(fp, dpi=130)
    out = {"omega_char": wc, "beta_h_scalar": bh,
           "scalar": {"psi_coh": tail(rs,'psi_coh'), "p4": p4([x['psi_coh'] for x in rs]),
                      "h_E": tail(rs,'h_E'), "rows": rs},
           "orders": {o: {"psi_coh": tail(ro[o],'psi_coh'),
                          "p4": p4([x['psi_coh'] for x in ro[o]]),
                          "h_E": tail(ro[o],'h_E'), "rows": ro[o]} for o in orders},
           "long": {"psi_coh": float(psi[-200:].mean()), "p4": p4(psi),
                    "p4_h": p4([x['h_E'] for x in rl]), "rows": rl},
           "checks": CK, "n_fail": sum(1 for c in CK if not c["pass"]),
           "anchor": "mode-dependent conformal redshift; ordering comparison; no fitting"}
    with open(os.path.join(RES, "zero_evo_1684_mode.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_evo_1684_mode.json"), "和", fp)


if __name__ == "__main__":
    main()
