#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_h2o_bond.py — H2O 键长：R_min = argmin_R E(R, θ_opt)。

与键角同一原理：键长是**源坐标**，由能量泛函极小定出。
  · Zero 胞腔（S4，化学精度）扫 (R,θ) → R_min
  · FCI/STO-3G、RHF/STO-3G、RHF/6-31G* 对照
实验：R(OH)=0.9584 Å，θ=104.45°。

运行：/usr/local/bin/python3 ZeroGPU/src/zero_h2o_bond.py
输出：results/zero_h2o_bond.json
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

from pyscf import gto, scf, ao2mo          # noqa: E402
import pyscf.fci.direct_spin1 as fd        # noqa: E402
from zero_cell_molecules import run_case   # noqa: E402

OUT = ROOT / "results" / "zero_h2o_bond.json"
OK, BAD, ROWS = [], [], []


def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    ROWS.append({"name": name, "ok": bool(cond), "detail": detail})
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def geom(R, ang):
    a = np.deg2rad(ang)
    return (f"O 0 0 0; H {R*np.sin(a/2):.6f} {R*np.cos(a/2):.6f} 0; "
            f"H {-R*np.sin(a/2):.6f} {R*np.cos(a/2):.6f} 0")


def rhf_e(basis, R, ang):
    mol = gto.M(atom=geom(R, ang), basis=basis, unit="Angstrom", verbose=0)
    return float(scf.RHF(mol).run().e_tot)


def fci_e(basis, R, ang):
    mol = gto.M(atom=geom(R, ang), basis=basis, unit="Angstrom", verbose=0)
    mf = scf.RHF(mol).run()
    norb, nelec = mol.nao_nr(), mol.nelectron
    h1 = np.einsum("pi,pq,qj->ij", mf.mo_coeff,
                   mol.intor("int1e_kin") + mol.intor("int1e_nuc"), mf.mo_coeff)
    h2 = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff), norb)
    e0, _ = fd.FCI().kernel(h1, h2, norb, nelec)
    return float(e0 + mol.energy_nuc())


def fit_min(xs, es):
    xs, es = np.array(xs, float), np.array(es, float)
    i = int(np.argmin(es))
    lo, hi = max(0, i - 2), min(len(xs), i + 3)
    c = np.polyfit(xs[lo:hi], es[lo:hi], 2)
    return float(-c[1] / (2 * c[0])), float(c[0])


def cell_e(R, ang=104.45):
    return run_case(f"bond_{R:.3f}_{ang:.0f}", geom(R, ang), "sto-3g")["E_s4"]


def main():
    print("== H2O 键长：R_min = argmin_R E(R, θ_opt) ==")
    exp = 0.9584

    # Zero 胞腔：先粗扫 (R,θ)，再在最优点沿 R 细化
    coarse = []
    for R in (0.98, 1.02, 1.06):
        for a in (96.0, 100.0, 104.0):
            coarse.append((R, a, cell_e(R, a)))
    Rb, ab, _ = min(coarse, key=lambda t: t[2])
    fine_R = [round(Rb + d, 3) for d in (-0.02, -0.01, 0.0, 0.01, 0.02)]
    fine = [(R, ab, cell_e(R, ab)) for R in fine_R]
    R_cell, c_cell = fit_min([t[0] for t in fine], [t[2] for t in fine])

    R_fci, c_fci = fit_min([0.96 + 0.02 * i for i in range(8)],
                           [fci_e("sto-3g", 0.96 + 0.02 * i, 98.0) for i in range(8)])
    R_rhf, c_rhf = fit_min([round(0.94 + 0.01 * i, 2) for i in range(11)],
                           [rhf_e("sto-3g", round(0.94 + 0.01 * i, 2), 100.0) for i in range(11)])
    R_631, c_631 = fit_min([round(0.90 + 0.01 * i, 2) for i in range(9)],
                           [rhf_e("6-31g*", round(0.90 + 0.01 * i, 2), 106.0) for i in range(9)])

    ck("B1 极小存在（曲率>0）", all(c > 0 for c in (c_cell, c_fci, c_rhf, c_631)),
       f"曲率 胞腔 {c_cell:.1f}、FCI {c_fci:.1f}、RHF {c_rhf:.1f}、6-31G* {c_631:.1f}")
    ck("B2 胞腔 R_min≈FCI", abs(R_cell - R_fci) < 0.03,
       f"胞腔 {R_cell:.4f} vs FCI {R_fci:.4f} Å")
    ck("B3 6-31G* 接近实验", abs(R_631 - exp) < 0.03,
       f"6-31G* {R_631:.4f} Å vs 实验 {exp:.4f}（差 {R_631 - exp:+.4f}）")
    ck("B4 STO-3G 偏长（基组效应）", R_fci > exp and R_rhf > exp,
       f"STO-3G {R_rhf:.4f}/{R_fci:.4f} Å > 实验；缺极化 ⇒ 键偏长")

    out = {"experiment_A": exp,
           "R_min": {"cell_S4": R_cell, "FCI_STO3G": R_fci,
                     "RHF_STO3G": R_rhf, "RHF_631Gstar": R_631},
           "cell_coarse_best": {"R": Rb, "ang": ab},
           "cell_fine": fine, "rows": ROWS,
           "pass": len(OK), "fail": len(BAD)}
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2))
    print(f"\n== {len(OK)}/{len(OK)+len(BAD)} PASS  ->  {OUT.relative_to(ROOT.parent)} ==")
    print(f"   R_min: 胞腔 {R_cell:.4f}  FCI {R_fci:.4f}  RHF {R_rhf:.4f}  6-31G* {R_631:.4f}  Å")
    if BAD:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
