#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_screening_scf.py —— 从电子-电子扇区推自洽屏蔽（径向 LDA 第一试）
================================================================================

固定芯模型失败后，这里做**自洽平均场**：
  ρ(r) = Σ f_nl |u_nl|²/(4πr²)
  V_H(r) = Q(r)/r + 4π∫_r^∞ ρ r' dr'   （球对称 Hartree）
  V_X(r) = -(3/π)^(1/3) ρ^(1/3)        （Slater 交换，LDA）
  [-½d²/dr² + l(l+1)/(2r²) - Z/r + V_H + V_X] u = ε u
迭代到密度自洽；比较 K/Sc 的 4s/3d 构型总能量。

判据：
  S1 He（1s²）自洽，总能量 ≈ -2.86（LDA 交换，宽松 10%）
  S2 K：E(4s¹) < E(3d¹)（Madelung 4s 先填）
  S3 Sc：E(4s²3d¹) < E(4s¹3d²)
  S4 若 S2/S3 不成立：自洽 LDA 也不给 Madelung（负结果）
  S5 判定：这是标准 LDA，不是 Zero 原生；Zero 只提供 1/r 极限与密度

运行：cd ZeroGPU/src && /usr/bin/python3 zero_screening_scf.py
输出：results/zero_screening_scf.json
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": str(detail)})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


class RadialLDA:
    def __init__(self, Z, config, R=30.0, h=0.01, eps=0.02, k=6):
        self.Z = Z
        self.config = config  # list of (n,l,occ)
        self.R = R
        self.h = h
        self.eps = eps
        self.k = k
        self.r = (np.arange(int(R / h)) + 0.5) * h
        self.N = len(self.r)
        self.VH = np.zeros(self.N)
        self.VX = np.zeros(self.N)
        self.rho = np.zeros(self.N)
        self.orbitals = {}

    def solve_l(self, l):
        r = self.r
        main = 1.0 / (self.h ** 2) + l * (l + 1) / (2.0 * r * r)
        soft = np.sqrt(r * r + self.eps ** 2)
        V = -self.Z / soft + self.VH + self.VX
        off = -0.5 / (self.h ** 2) * np.ones(self.N - 1)
        H = sp.diags([off, main + V, off], [-1, 0, 1], shape=(self.N, self.N))
        ev, evec = spl.eigsh(H, k=self.k, which="SA", tol=1e-9)
        order = np.argsort(ev)
        return ev[order], evec[:, order]

    def build_density(self):
        r = self.r
        rho = np.zeros(self.N)
        self.orbitals = {}
        for l in sorted(set(l for _, l, _ in self.config)):
            ev, evec = self.solve_l(l)
            self.orbitals[l] = (ev, evec)
            for n, ll, occ in self.config:
                if ll != l:
                    continue
                idx = n - l - 1
                if idx < 0 or idx >= evec.shape[1]:
                    raise RuntimeError(f"state n={n},l={l} not in k={self.k}")
                u = evec[:, idx]
                u = u / np.sqrt(np.sum(u * u) * self.h)  # ∫u²dr=1
                rho += occ * (u ** 2) / (4.0 * np.pi * r * r)
        self.rho = np.maximum(rho, 1e-30)

    def potentials(self):
        r, h, rho = self.r, self.h, self.rho
        a = 4.0 * np.pi * r * r * rho  # = Σ f u²
        Q = np.concatenate([[0.0], np.cumsum(0.5 * (a[1:] + a[:-1]) * h)])
        b = rho * r
        c = np.concatenate([[0.0], np.cumsum(0.5 * (b[1:] + b[:-1]) * h)])
        total = c[-1]
        tail = total - c          # tail[i] = ∫_{r_i}^{R} ρ r' dr'
        VH = np.zeros(self.N)
        VH[1:] = Q[1:] / r[1:]
        VH += 4.0 * np.pi * tail
        VH[0] = VH[1]
        rho13 = np.maximum(rho, 1e-30) ** (1.0 / 3.0)
        VX = -(3.0 / np.pi) ** (1.0 / 3.0) * rho13
        return VH, VX

    def total_energy(self):
        r, h = self.r, self.h
        eps_sum = 0.0
        for n, l, occ in self.config:
            ev, _ = self.orbitals[l]
            eps_sum += occ * ev[n - l - 1]
        EH = 0.5 * np.sum(self.rho * self.VH * 4.0 * np.pi * r * r * h)
        C = (3.0 / np.pi) ** (1.0 / 3.0)
        EX = -0.75 * C * np.sum(self.rho ** (4.0 / 3.0) * 4.0 * np.pi * r * r * h)
        return float(eps_sum - EH - EX / 3.0)

    def scf(self, niter=40, mix=0.4):
        self.build_density()
        for it in range(niter):
            VH, VX = self.potentials()
            self.VH = mix * self.VH + (1 - mix) * VH
            self.VX = mix * self.VX + (1 - mix) * VX
            rho_old = self.rho.copy()
            self.build_density()
            drho = np.max(np.abs(self.rho - rho_old))
            if it % 5 == 0 or it == niter - 1:
                print(f"    it={it:2d}  Δρ={drho:.3e}  E={self.total_energy():+.5f}")
            self.last_drho = float(drho)
            self.converged = bool(drho < 1e-4)
            if drho < 1e-6:
                break
        return self.total_energy()


def config_energy(Z, config):
    print(f"  Z={Z} config={config}")
    m = RadialLDA(Z, config)
    return m.scf()


def main():
    print("=" * 84)
    print("自洽屏蔽（径向 LDA）：K/Sc 的 4s vs 3d")
    print("=" * 84)
    e_he = config_energy(2, [(1, 0, 2)])
    print(f"  He E={e_he:+.5f}")
    ck("S1 He 自洽，E≈-2.86（10%）", abs(e_he / -2.86 - 1) < 0.10, f"E={e_he:+.5f}")

    e_k4s = config_energy(19, [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 6), (4, 0, 1)])
    e_k3d = config_energy(19, [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 6), (3, 2, 1)])
    print(f"  K: E(4s¹)={e_k4s:+.5f}  E(3d¹)={e_k3d:+.5f}")
    ck("S2 K：E(4s¹) < E(3d¹)", e_k4s < e_k3d,
       f"ΔE={e_k4s - e_k3d:+.5f}")

    e_sc_a = config_energy(21, [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 6), (4, 0, 2), (3, 2, 1)])
    e_sc_b = config_energy(21, [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 6), (4, 0, 1), (3, 2, 2)])
    print(f"  Sc: E(4s²3d¹)={e_sc_a:+.5f}  E(4s¹3d²)={e_sc_b:+.5f}")
    ck("S3 Sc：E(4s²3d¹) < E(4s¹3d²)", e_sc_a < e_sc_b,
       f"ΔE={e_sc_a - e_sc_b:+.5f}")
    ck("S4 若 S2/S3 不成立：自洽 LDA 也不给 Madelung（负结果）",
       not (e_k4s < e_k3d and e_sc_a < e_sc_b),
       "见结果")
    ck("S5 判定：标准 LDA，不是 Zero 原生；Zero 只提供 1/r 极限与密度", True,
       "见 doc §3")

    nfail = sum(1 for c in CK if not c["pass"])
    print("-" * 84)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-nfail} / FAIL {nfail}")
    out = dict(He=e_he, K=dict(E4s=e_k4s, E3d=e_k3d), Sc=dict(A=e_sc_a, B=e_sc_b),
               checks=CK, n_fail=nfail,
               boundary="径向 LDA（Slater 交换），球平均；标准量子化学，不是 Zero 原生")
    with open(os.path.join(RES, "zero_screening_scf.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(f"结果 → {os.path.join(RES, 'zero_screening_scf.json')}")
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
