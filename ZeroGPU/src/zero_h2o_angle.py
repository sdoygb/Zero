#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_h2o_angle.py — 键角住在哪一层？怎样算？

结论（本脚本核验）：
  · 键角不是任何一层的变量：L0 是抽象图、L1 是离散词、L1′ 是记录计数、
    L2 的 3 槽是基底常数。键角是**源**（核的位置）的坐标。
  · 它的值由 **L2 能量泛函在源上取极小** 定出（变分原理）：E(R,θ) 的极小点。
  · 数值：STO-3G（Zero 胞腔管线的基组）给 θ≈99–100°，实验 104.45°；
    6-31G 无极化给更宽的角——差的是基组，不是方法。

运行：/usr/local/bin/python3 ZeroGPU/src/zero_h2o_angle.py
输出：results/zero_h2o_angle.json
"""
from __future__ import annotations

import json
import pathlib

import numpy as np
from pyscf import gto, scf, ao2mo
import pyscf.fci.direct_spin1 as fd

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "zero_h2o_angle.json"
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


def main():
    print("== H2O 键角：住在哪层、怎样算 ==")
    R0 = 0.9584
    xs = list(range(90, 117, 2))
    e_rhf = [rhf_e("sto-3g", R0, x) for x in xs]
    t_rhf, c_rhf = fit_min(xs, e_rhf)
    e_fci = [fci_e("sto-3g", R0, x) for x in xs]
    t_fci, c_fci = fit_min(xs, e_fci)
    e_631 = [rhf_e("6-31g", R0, x) for x in xs]
    t_631, c_631 = fit_min(xs, e_631)

    ck("A1 极小存在（STO-3G，曲率>0）", c_rhf > 0 and c_fci > 0,
       f"RHF 曲率 {c_rhf:.2e}、FCI 曲率 {c_fci:.2e}")
    ck("A2 FCI/STO-3G 键角", abs(t_fci - 104.45) < 6.0,
       f"θ_min={t_fci:.2f}°（R={R0} Å；实验 104.45°，差 {t_fci - 104.45:+.2f}°）")
    ck("A3 基组效应方向", t_631 > t_fci,
       f"6-31G 无极化 θ_min={t_631:.2f}° > STO-3G {t_fci:.2f}°（极化函数才会把角压回 105° 附近）")

    best_star = None
    for R in np.arange(0.92, 0.99, 0.02):
        for a in range(102, 111, 2):
            e = rhf_e("6-31g*", float(R), a)
            if best_star is None or e < best_star[2]:
                best_star = (float(R), float(a), e)
    ck("A3b 极化函数把角拉回", abs(best_star[1] - 104.45) < abs(t_631 - 104.45),
       f"6-31G* θ_min={best_star[1]:.0f}°（R={best_star[0]:.2f} Å）vs 6-31G {t_631:.1f}°；实验 104.45°")

    # 2D：RHF/STO-3G 的 (R,θ) 极小
    best = None
    for R in np.arange(0.93, 1.06, 0.03):
        for a in range(96, 109, 2):
            e = rhf_e("sto-3g", float(R), a)
            if best is None or e < best[2]:
                best = (float(R), float(a), e)
    ck("A4 2D 极小（RHF/STO-3G）", abs(best[0] - 0.99) < 0.06 and abs(best[1] - 100) < 5,
       f"(R,θ)=({best[0]:.3f} Å, {best[1]:.0f}°)")

    # Zero 胞腔方法的三点（已有数据）也给出同一极小
    cell = json.loads((ROOT / "results" / "zero_cell_molecules.json").read_text())
    tbl = {r["name"]: r for r in cell["table"]}
    pts = {"H2O_96": 96.0, "H2O_100": 100.0, "H2O_104": 104.45}
    cxs = [pts[k] for k in pts]
    ces = [tbl[k]["E_s4"] for k in pts]
    t_cell, c_cell = fit_min(cxs, ces)
    ck("A5 Zero 胞腔三点与 FCI 同向", c_cell > 0 and abs(t_cell - t_fci) < 4.0,
       f"胞腔 S4 三点拟合 θ_min={t_cell:.2f}°（FCI {t_fci:.2f}°）")

    out = {
        "method": "键角=源坐标；由 E(R,θ) 的极小定出（变分原理）",
        "scan": {"R": R0, "angles": xs,
                 "RHF_STO3G": e_rhf, "FCI_STO3G": e_fci, "RHF_631G": e_631},
        "minima": {"RHF_STO3G": t_rhf, "FCI_STO3G": t_fci, "RHF_631G": t_631,
                   "RHF_631Gstar": best_star, "cell_S4_3pt": t_cell,
                   "2D_RHF_STO3G": best},
        "rows": ROWS, "pass": len(OK), "fail": len(BAD),
    }
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2))
    print(f"\n== {len(OK)}/{len(OK)+len(BAD)} PASS  ->  {OUT.relative_to(ROOT.parent)} ==")
    if BAD:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
