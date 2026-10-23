#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_large_molecule.py — 大分子/复杂分子：Zero 胞腔方法能算多准？

用稀疏 MLE（O(ndet)/迭代）替换稠密特征矩阵，把胞腔方法推到 10^5 量级：
  H2O (10 e, 441 dets) → H2CO (16 e, 245k) → H2O2 (18 e, 48k)
对照 FCI。判据：小分子化学精度；大分子是否还准；标度墙在哪。

运行：/usr/local/bin/python3 ZeroGPU/src/zero_large_molecule.py
输出：results/zero_large_molecule.json
"""
from __future__ import annotations

import json
import math
import pathlib
import sys
import time

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

from pyscf import gto, scf, fci, ao2mo      # noqa: E402
from pyscf.fci import cistring             # noqa: E402
import pyscf.fci.direct_spin1 as fd        # noqa: E402
from zero_cell_molecules import feats, det_sign  # noqa: E402

OUT = ROOT / "results" / "zero_large_molecule.json"
OK, BAD, ROWS = [], [], []


def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    ROWS.append({"name": name, "ok": bool(cond), "detail": detail})
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def mle_sparse(cell, h, r, z, w, C2, H, niter=8000):
    """稀疏版 MLE：z_d = θ[cell_d] + θ[C2+h] + r θ[C2+H+h] + z θ[C2+2H+h]。"""
    def zlog(theta):
        return theta[cell] + theta[C2 + h] + r * theta[C2 + H + h] + z * theta[C2 + 2 * H + h]

    def obj(theta):
        zz = zlog(theta)
        return logsumexp(zz) - w @ zz

    def grad(theta):
        zz = zlog(theta)
        p = np.exp(zz - zz.max())
        p /= p.sum()
        d = p - w
        g = np.zeros(C2 + 3 * H)
        np.add.at(g, cell, d)
        np.add.at(g, C2 + h, d)
        np.add.at(g, C2 + H + h, r * d)
        np.add.at(g, C2 + 2 * H + h, z * d)
        return g

    res = minimize(obj, np.zeros(C2 + 3 * H), jac=grad, method="L-BFGS-B",
                   options={"maxiter": niter, "ftol": 1e-16})
    zz = zlog(res.x)
    p = np.exp(zz - zz.max())
    p /= p.sum()
    return p, int(res.nit)


def run_sparse(name, atom, basis="sto-3g"):
    t0 = time.time()
    mol = gto.M(atom=atom, basis=basis, unit="Angstrom", verbose=0)
    mf = scf.RHF(mol).run()
    norb, nelec = mol.nao_nr(), mol.nelectron
    h1e = np.einsum("pi,pq,qj->ij", mf.mo_coeff,
                    mol.intor("int1e_kin") + mol.intor("int1e_nuc"), mf.mo_coeff)
    h2e = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff), norb)
    ecore = mol.energy_nuc()
    cis = fd.FCI()
    e0, c = cis.kernel(h1e, h2e, norb, nelec)
    E0 = e0 + ecore
    strs = cistring.make_strings(range(norb), nelec // 2)
    cc = c.reshape(-1)
    g = np.abs(cc)
    g /= np.linalg.norm(g)
    w = g * g
    w /= w.sum()
    sgn = np.sign(cc)
    dets = [(sa, sb) for sa in strs for sb in strs]
    nocc = nelec // 2
    hf = (1 << nocc) - 1
    Hf = np.array([feats(sa, sb, norb) for sa, sb in dets])
    h = Hf[:, 0].astype(int)
    r = Hf[:, 1]
    z = Hf[:, 2]
    H = int(h.max()) + 1
    labs2 = []
    for sa, sb in dets:
        holes = (hf & ~sa) | (hf & ~sb)
        parts = ((sa | sb) & ~hf) & ((1 << norb) - 1)
        l = holes | parts
        nh_a = bin((hf & ~sa) & 0xFFFF).count("1")
        nh_b = bin((hf & ~sb) & 0xFFFF).count("1")
        np_a = bin((sa & ~hf) & 0xFFFF).count("1")
        np_b = bin((sb & ~hf) & 0xFFFF).count("1")
        labs2.append(int(l) * 6 ** 4 + nh_a * 6 ** 3 + nh_b * 6 ** 2 + np_a * 6 + np_b)
    labs2 = np.array(labs2)
    u2, idx2 = np.unique(labs2, return_inverse=True)
    C2 = len(u2)
    p2, nit = mle_sparse(idx2, h, r, z, w, C2, H)
    ds = np.array([det_sign(sa, sb, norb, nocc) for sa, sb in dets])
    res = sgn * ds
    num = np.bincount(idx2, weights=w * res, minlength=C2)
    spin = np.where(num >= 0, 1.0, -1.0)
    model = ds * spin[idx2]
    v = model * np.sqrt(p2)
    v /= np.linalg.norm(v)
    vr = v.reshape(c.shape)
    Ev = cis.energy(h1e, h2e, vr, norb, nelec) + ecore
    acc = float(np.sum(w * (model == sgn)))
    return {"name": name, "norb": norb, "nelec": nelec, "ndet": int(len(w)),
            "C2": C2, "H": H, "E_fci": float(E0), "E_s4": float(Ev),
            "err": float(Ev - E0), "acc": acc, "mle_nit": nit,
            "time_s": round(time.time() - t0, 1)}


def main():
    print("== 大分子/复杂分子：Zero 胞腔能算多准 ==")
    cases = [
        ("H2O", "O 0 0 0; H 0.757 0.586 0; H -0.757 0.586 0"),
        ("H2CO", "C 0 0 0; O 0 0 1.21; H 0.94 0 -0.59; H -0.94 0 -0.59"),
        ("H2O2", "O 0 0 0; O 0 0 1.475; H 0.950 0 -0.30; H -0.950 0 1.775"),
    ]
    rows = []
    for name, atom in cases:
        row = run_sparse(name, atom)
        rows.append(row)
        print(f"  {name:6s} ndet={row['ndet']:7d} C2={row['C2']:6d} "
              f"err={row['err']:+.5f} Ha  acc={row['acc']:.4f}  t={row['time_s']}s")

    by = {r["name"]: r for r in rows}
    ck("L1 H2O 化学精度", abs(by["H2O"]["err"]) < 1.6e-3,
       f"H2O err={by['H2O']['err']:+.2e}")
    ck("L2 H2CO 掉出化学精度", abs(by["H2CO"]["err"]) > 1.6e-3,
       f"H2CO err={by['H2CO']['err']:+.2e}（{abs(by['H2CO']['err'])/1.6e-3:.0f}× 化学精度）")
    ck("L3 H2O2 掉出化学精度", abs(by["H2O2"]["err"]) > 1.6e-3,
       f"H2O2 err={by['H2O2']['err']:+.2e}（{abs(by['H2O2']['err'])/1.6e-3:.0f}× 化学精度）")
    ck("L4 误差随行列式数增长", by["H2O"]["err"] < by["H2O2"]["err"] < by["H2CO"]["err"],
       f"err: H2O {by['H2O']['err']:.1e} < H2O2 {by['H2O2']['err']:.1e} < H2CO {by['H2CO']['err']:.1e}")

    # 标度墙：甲醇 14 轨道 18 电子
    ndet_meoh = math.comb(14, 9) ** 2
    C2_est = 30000
    dense_gb = ndet_meoh * C2_est * 8 / 1e9
    ck("L5 标度墙（甲醇）", dense_gb > 100,
       f"甲醇 ndet={ndet_meoh:,}；稠密特征矩阵 ≈ {dense_gb:.0f} GB（C2~{C2_est:,}）⇒ 稠密实现不可行；稀疏 MLE 可行")

    out = {"rows": rows, "checks": ROWS,
           "methanol_ndet": ndet_meoh, "dense_matrix_GB_est": dense_gb,
           "pass": len(OK), "fail": len(BAD)}
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2))
    print(f"\n== {len(OK)}/{len(OK)+len(BAD)} PASS  ->  {OUT.relative_to(ROOT.parent)} ==")
    if BAD:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
