#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""H2O 轨道对空间胞腔：把行列式按 (空穴集, 粒子集) 分组，胞腔+3H 的 MLE。

模型：p_d = softmax( cell one-hot + layer one-hot + 每层 runs/zc )
  3H：      层 one-hot + 每层 runs/zc
  cell：    空间激发胞腔（α/β 空穴/粒子的并集）+ 3H
符号：oracle（FCI 符号），只测幅度模型的能量。
判据
  H1 FCI 与 -75.01243743 一致
  H2 3H 基线 err ~0.0997 Ha
  H3 cell F > 3H F
  H4 |cell err| < |3H err|
  H5 |cell err| < 1.6e-3 Ha（化学精度；目标）
  H6 胞腔数 C<=128，参数 poly
输出：results/zero_h2o_cells.json
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

try:
    from pyscf import gto, scf, fci, ao2mo
    import pyscf.fci.direct_spin1 as fd
except Exception as exc:
    raise SystemExit("需要 pyscf（/usr/local/bin/python3）：%r" % (exc,))

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "zero_h2o_cells.json"
OK, BAD, ROWS = [], [], []


def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    ROWS.append({"name": name, "ok": bool(cond), "detail": detail})
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def build():
    atom = "O 0 0 0; H 0.757 0.586 0; H -0.757 0.586 0"
    mol = gto.M(atom=atom, basis="sto-3g", unit="Angstrom", verbose=0)
    mf = scf.RHF(mol).run()
    norb, nelec = mol.nao_nr(), mol.nelectron
    h1e = np.einsum("pi,pq,qj->ij", mf.mo_coeff,
                    mol.intor("int1e_kin") + mol.intor("int1e_nuc"), mf.mo_coeff)
    h2e = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff), norb)
    ecore = mol.energy_nuc()
    cis = fd.FCI()
    e0, c = cis.kernel(h1e, h2e, norb, nelec)
    strs = fci.cistring.make_strings(range(norb), nelec // 2)
    return mol, mf, h1e, h2e, ecore, cis, e0, c, strs, norb, nelec


def feats(sa, sb, norb):
    bits = []
    for i in range(norb):
        bits += [(sa >> i) & 1, (sb >> i) & 1]
    s = mx = zc = 0
    runs = 1
    for i, x in enumerate(bits):
        s += 2 * x - 1
        mx = max(mx, abs(s))
        if i > 0 and bits[i] != bits[i - 1]:
            runs += 1
        if 0 < i < len(bits) - 1 and s == 0:
            zc += 1
    return mx, runs, zc


def mle(f, w, niter=4000):
    def obj(theta):
        z = f @ theta
        return logsumexp(z) - np.dot(w, z)
    def grad(theta):
        z = f @ theta
        p = np.exp(z - z.max()); p /= p.sum()
        return f.T @ (p - w)
    res = minimize(obj, np.zeros(f.shape[1]), jac=grad, method="L-BFGS-B",
                   options={"maxiter": niter, "ftol": 1e-16})
    z = f @ res.x
    p = np.exp(z - z.max()); p /= p.sum()
    return p


def main():
    print("=" * 74)
    print("zero_h2o_cells.py — H2O 轨道对空间胞腔")
    print("=" * 74)
    mol, mf, h1e, h2e, ecore, cis, e0, c, strs, norb, nelec = build()
    E0 = e0 + ecore
    print(f"    FCI {E0:.8f}  norb {norb}  nelec {nelec}  c {c.shape}")
    ndet = c.size
    cc = c.reshape(-1)
    g = np.abs(cc); g /= np.linalg.norm(g)
    w = g ** 2; w /= w.sum()
    sgn = np.sign(cc)
    dets = []
    for ia, sa in enumerate(strs):
        for ib, sb in enumerate(strs):
            dets.append((sa, sb))
    nocc = nelec // 2
    hf = (1 << nocc) - 1
    Hf = np.array([feats(sa, sb, norb) for sa, sb in dets])
    h = Hf[:, 0].astype(int); r = Hf[:, 1].astype(int); z = Hf[:, 2].astype(int)
    H = int(h.max()) + 1
    # 3H
    f3 = np.zeros((ndet, 3 * H))
    for i in range(ndet):
        f3[i, h[i]] = 1.0
        f3[i, H + h[i]] = r[i]
        f3[i, 2 * H + h[i]] = z[i]
    p3 = mle(f3, w)
    # cell label: union of alpha/beta holes and particles
    labels = []
    for sa, sb in dets:
        holes = (hf & ~sa) | (hf & ~sb)
        parts = ((sa | sb) & ~hf) & 0b1100000
        labels.append(holes | parts)
    labels = np.array(labels)
    uniq = np.unique(labels)
    C = len(uniq)
    lab_idx = {u: i for i, u in enumerate(uniq)}
    fc = np.zeros((ndet, C + 3 * H))
    for i in range(ndet):
        fc[i, lab_idx[labels[i]]] = 1.0
        fc[i, C + h[i]] = 1.0
        fc[i, C + H + h[i]] = r[i]
        fc[i, C + 2 * H + h[i]] = z[i]
    pc = mle(fc, w)
    def F(p): return float((g @ np.sqrt(p)) ** 2)
    def E(p):
        v = sgn * np.sqrt(p); v = v / np.linalg.norm(v); v = v.reshape(c.shape)
        return float(cis.energy(h1e, h2e, v, norb, nelec) + ecore)
    F3, Fc = F(p3), F(pc)
    E3, Ec = E(p3), E(pc)
    print(f"    3H: F={F3:.5f} E={E3:.8f} err={E3-E0:+.6f} params={3*H}")
    print(f"    cell: F={Fc:.5f} E={Ec:.8f} err={Ec-E0:+.6f} params={C+3*H} C={C}")
    # cell table (top)
    cell_w = np.zeros(C)
    for i in range(ndet):
        cell_w[lab_idx[labels[i]]] += w[i]
    top = np.argsort(-cell_w)[:8]
    print("    top cells (label, weight, holes, parts):")
    for t in top:
        lab = int(uniq[t])
        print(f"      {lab:3d}  w={cell_w[t]:.5f}  holes={lab & 31:05b}  parts={(lab>>5)&3:02b}")
    ck("H1 FCI matches -75.012437", abs(E0 + 75.01243743249299) < 1e-6, f"{E0:.8f}")
    ck("H2 3H baseline err ~0.0997", abs((E3 - E0) - 0.0996906) < 2e-3, f"{E3-E0:+.6f}")
    ck("H3 cell F > 3H F", Fc > F3, f"{Fc:.5f} vs {F3:.5f}")
    ck("H4 |cell err| < |3H err|", abs(Ec - E0) < abs(E3 - E0), f"{Ec-E0:+.6f} vs {E3-E0:+.6f}")
    ck("H5 cell err < 1.6e-3 (chemical accuracy)", abs(Ec - E0) < 1.6e-3, f"{Ec-E0:+.2e}")
    ck("H6 cell count C<=128, params poly", C <= 128, f"C={C}, params={C+3*H}")
    out = {"script": "zero_h2o_cells.py", "checks": ROWS, "n_pass": len(OK), "n_fail": len(BAD),
           "e_fci": float(E0), "F3": F3, "E3": E3, "err3": float(E3 - E0),
           "Fc": Fc, "Ec": Ec, "errc": float(Ec - E0), "C": C, "H": H,
           "params3": 3 * H, "paramsc": C + 3 * H,
           "top_cells": [{"label": int(uniq[t]), "w": float(cell_w[t])} for t in top]}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n{len(OK)} PASS / {len(BAD)} FAIL  ->  {OUT}")
    return 0 if not BAD else 1


if __name__ == "__main__":
    raise SystemExit(main())
