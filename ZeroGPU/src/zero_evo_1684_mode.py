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
    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.omega = np.sqrt(np.maximum(self.wave_eigs, 0.0))
        self.omega_op = self.wave_vecs @ np.diag(self.omega) @ self.wave_vecs.T

    def _apply_mode_phase(self, dt):
        H = np.diag(self.h - self.h.mean())
        G = 0.5 * (self.omega_op @ H + H @ self.omega_op)
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


def run(mode="scalar", bh=0.0, steps=2000):
    if mode == "scalar":
        e = ZeroEvo1684Cycle16(zs.cycle(16), Params(L=4, T=4, micro=True, c_h=1.0, dt=0.05,
                                                      alpha=0.1, g_src=1.0, beta=0.0, beta_h=bh),
                                seed=np.ones(16))
    else:
        e = ZeroEvo1684Mode(zs.cycle(16), Params(L=4, T=4, micro=True, c_h=1.0, dt=0.05,
                                                   alpha=0.1, g_src=1.0, beta=0.0, beta_h=0.0),
                             seed=np.ones(16))
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
    rm = run("mode", 0.0, 2000)
    def tail(r, k): return float(np.mean([x[k] for x in r[-200:]]))
    print("mode     q       q16     h_E       psi_coh  p4psi")
    for name, r in [("scalar", rs), ("mode", rm)]:
        print(f"{name:8s} {tail(r,'q'):.4f}  {tail(r,'q16'):.4f}  {tail(r,'h_E'):.3f}  "
              f"{tail(r,'psi_coh'):.4f}  {p4([x['psi_coh'] for x in r]):.4f}")
    ck("scalar beta_h=omega/2 有周期4", p4([x['psi_coh'] for x in rs]) > 1e-6, "")
    ck("mode 相位有周期4", p4([x['psi_coh'] for x in rm]) > 1e-6,
       f"{p4([x['psi_coh'] for x in rm]):.4f}")
    ck("两者 psi_coh 都有界", max(x['psi_coh'] for x in rs) <= 1.0 + 1e-9 and
       max(x['psi_coh'] for x in rm) <= 1.0 + 1e-9, "")
    ck("两者 h_E 相近", abs(tail(rs,'h_E') - tail(rm,'h_E')) < 1e-6, "")

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    ax[0].plot([x['psi_coh'] for x in rs], lw=1, label="scalar ω/2")
    ax[0].plot([x['psi_coh'] for x in rm], lw=1, label="mode")
    ax[0].set_title("(a) psi_coh"); ax[0].legend(fontsize=8)
    ax[1].plot([x['h_E'] for x in rm], lw=1, label="mode h_E")
    ax[1].set_title("(b) h_E"); ax[1].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout(); fp = os.path.join(RES, "zero_evo_1684_mode.png")
    fig.savefig(fp, dpi=130)
    out = {"omega_char": wc, "beta_h_scalar": bh,
           "scalar": {"psi_coh": tail(rs,'psi_coh'), "p4": p4([x['psi_coh'] for x in rs]),
                      "h_E": tail(rs,'h_E'), "rows": rs},
           "mode": {"psi_coh": tail(rm,'psi_coh'), "p4": p4([x['psi_coh'] for x in rm]),
                    "h_E": tail(rm,'h_E'), "rows": rm},
           "checks": CK, "n_fail": sum(1 for c in CK if not c["pass"]),
           "anchor": "mode-dependent conformal redshift omega_op; no fitting"}
    with open(os.path.join(RES, "zero_evo_1684_mode.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_evo_1684_mode.json"), "和", fp)


if __name__ == "__main__":
    main()
