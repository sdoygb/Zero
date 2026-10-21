#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Zero 胞腔多分子批跑：H2O 几何/LiH 键长，自旋胞腔 + S4 native 符号 + 偶极/1-RDM。

对每个 (name, atom)：
  FCI；自旋胞腔（空间胞腔+α/β计数）+ 3H 的 MLE；激发符号×自旋胞腔多数（S4）；
  能量误差；偶极；1-RDM 差；胞腔数/参数。
输出：results/zero_cell_molecules.json
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
OUT = ROOT / "results" / "zero_cell_molecules.json"
OK, BAD, ROWS, TBL = [], [], [], []


def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    ROWS.append({"name": name, "ok": bool(cond), "detail": detail})
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


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


def det_sign(sa, sb, norb, nocc):
    hfstr = (1 << nocc) - 1
    def sstr(hf_s, d_s):
        holes = hf_s & ~d_s
        parts = d_s & ~hf_s
        sign = 1
        for h in range(norb):
            if (holes >> h) & 1:
                for p in range(norb):
                    if (parts >> p) & 1 and h < p:
                        sign = -sign
        return sign
    return sstr(hfstr, sa) * sstr(hfstr, sb)


def majority(groups, res, w, ng):
    out = np.ones(ng)
    for g in range(ng):
        m = groups == g
        if m.any():
            out[g] = 1.0 if np.sum(w[m] * res[m]) >= 0 else -1.0
    return out


def run_case(name, atom, basis="sto-3g"):
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
    strs = fci.cistring.make_strings(range(norb), nelec // 2)
    ndet = c.size
    cc = c.reshape(-1)
    g = np.abs(cc); g /= np.linalg.norm(g)
    w = g ** 2; w /= w.sum()
    sgn = np.sign(cc)
    dets = [(sa, sb) for sa in strs for sb in strs]
    nocc = nelec // 2
    hf = (1 << nocc) - 1
    Hf = np.array([feats(sa, sb, norb) for sa, sb in dets])
    h = Hf[:, 0].astype(int); r = Hf[:, 1].astype(int); z = Hf[:, 2].astype(int)
    H = int(h.max()) + 1
    labels = []
    for sa, sb in dets:
        holes = (hf & ~sa) | (hf & ~sb)
        parts = ((sa | sb) & ~hf) & ((1 << norb) - 1)
        labels.append(holes | parts)
    labels = np.array(labels)
    uniq = np.unique(labels); C = len(uniq)
    lab_idx = {u: i for i, u in enumerate(uniq)}
    lab_arr = np.array([lab_idx[l] for l in labels])
    labels2 = []
    for i, (sa, sb) in enumerate(dets):
        nh_a = bin((hf & ~sa) & 0xFFFF).count("1"); nh_b = bin((hf & ~sb) & 0xFFFF).count("1")
        np_a = bin((sa & ~hf) & 0xFFFF).count("1"); np_b = bin((sb & ~hf) & 0xFFFF).count("1")
        labels2.append(int(labels[i]) * 6 ** 4 + nh_a * 6 ** 3 + nh_b * 6 ** 2 + np_a * 6 + np_b)
    labels2 = np.array(labels2)
    uniq2 = np.unique(labels2); C2 = len(uniq2)
    idx2 = {u: i for i, u in enumerate(uniq2)}
    idx2_arr = np.array([idx2[l] for l in labels2])
    f2 = np.zeros((ndet, C2 + 3 * H))
    for i in range(ndet):
        f2[i, idx2_arr[i]] = 1.0
        f2[i, C2 + h[i]] = 1.0
        f2[i, C2 + H + h[i]] = r[i]
        f2[i, C2 + 2 * H + h[i]] = z[i]
    p2 = mle(f2, w)
    ds = np.array([det_sign(sa, sb, norb, nocc) for sa, sb in dets])
    res = sgn * ds
    spin_sign = majority(idx2_arr, res, w, C2)
    model = ds * spin_sign[idx2_arr]
    v = model * np.sqrt(p2); v = v / np.linalg.norm(v); vr = v.reshape(c.shape)
    Ev = cis.energy(h1e, h2e, vr, norb, nelec) + ecore
    acc = float(np.sum(w * (model == sgn)))
    # observables
    dm_v = cis.make_rdm1(vr, norb, nelec)
    dm_c = cis.make_rdm1(c, norb, nelec)
    dip_ao = mol.intor("int1e_r", comp=3)
    dip_mo = np.einsum("pi,xpq,qj->xij", mf.mo_coeff, dip_ao, mf.mo_coeff)
    nuc = np.einsum("i,ix->x", mol.atom_charges(), mol.atom_coords())
    el_v = np.array([np.trace(dm_v @ dip_mo[x]) for x in range(3)])
    el_c = np.array([np.trace(dm_c @ dip_mo[x]) for x in range(3)])
    dip_v = (nuc - el_v) * 2.541746
    dip_c = (nuc - el_c) * 2.541746
    d_dip = float(np.linalg.norm(dip_v - dip_c))
    d_rdm = float(np.linalg.norm(dm_v - dm_c))
    row = {"name": name, "norb": norb, "nelec": nelec, "ndet": ndet,
           "E_fci": float(E0), "E_s4": float(Ev), "err": float(Ev - E0), "acc": acc,
           "F": float((cc / np.linalg.norm(cc) @ v) ** 2),
           "dip_fci": float(dip_c[1]) if len(dip_c) > 1 else float(np.linalg.norm(dip_c)),
           "dip_s4": float(dip_v[1]) if len(dip_v) > 1 else float(np.linalg.norm(dip_v)),
           "dip_diff": d_dip, "rdm_diff": d_rdm,
           "C": C, "C2": C2, "H": H, "params2": C2 + 3 * H}
    return row


def main():
    print("=" * 74)
    print("zero_cell_molecules.py — 多分子胞腔批跑（S4 native 符号）")
    print("=" * 74)
    def h2o(R=0.9584, ang=104.45):
        th = np.deg2rad(ang); x = R * np.sin(th / 2); y = R * np.cos(th / 2)
        return f"O 0 0 0; H {x:.6f} {y:.6f} 0; H {-x:.6f} {y:.6f} 0"
    cases = [
        ("H2O_104", h2o(), "sto-3g"),
        ("H2O_96", h2o(0.9584, 96), "sto-3g"),
        ("H2O_100", h2o(0.9584, 100), "sto-3g"),
        ("H2O_R0.96", h2o(0.96, 104.45), "sto-3g"),
        ("H2O_R1.02", h2o(1.02, 104.45), "sto-3g"),
        ("LiH_1.3", "Li 0 0 0; H 0 0 1.3", "sto-3g"),
        ("LiH_1.6", "Li 0 0 0; H 0 0 1.6", "sto-3g"),
        ("LiH_2.0", "Li 0 0 0; H 0 0 2.0", "sto-3g"),
        ("LiH_1.6_631g", "Li 0 0 0; H 0 0 1.6", "6-31g"),
    ]
    for name, atom, basis in cases:
        row = run_case(name, atom, basis)
        TBL.append(row)
        print(f"  {name:10s} n={row['ndet']:4d} E_fci={row['E_fci']:.6f} err={row['err']:+.6f} "
              f"acc={row['acc']:.4f} dip={row['dip_s4']:.4f}/{row['dip_fci']:.4f} "
              f"dRDM={row['rdm_diff']:.4f} C2={row['C2']} p={row['params2']}")
    h2o_rows = [x for x in TBL if x["name"].startswith("H2O")]
    lih_rows = [x for x in TBL if x["name"].startswith("LiH")]
    ck("M1 H2O all err <1.6e-3", all(abs(x["err"]) < 1.6e-3 for x in h2o_rows),
       f"max={max(abs(x['err']) for x in h2o_rows):.2e}")
    ck("M2 LiH all err <1.6e-3", all(abs(x["err"]) < 1.6e-3 for x in lih_rows),
       f"max={max(abs(x['err']) for x in lih_rows):.2e}")
    ck("M3 dipoles diff <0.05 D", all(x["dip_diff"] < 0.05 for x in TBL),
       f"max={max(x['dip_diff'] for x in TBL):.4f}")
    ck("M4 1-RDM diff <0.01", all(x["rdm_diff"] < 0.01 for x in TBL),
       f"max={max(x['rdm_diff'] for x in TBL):.4f}")
    ck("M5 params < ndet (poly reduction)", all(x["params2"] < x["ndet"] for x in TBL),
       f"max ratio ndet/params={max(x['ndet']/x['params2'] for x in TBL):.1f}")
    g631 = [x for x in TBL if x["name"].endswith("631g")]
    ck("M6 LiH 6-31G err <1.6e-3", all(abs(x["err"]) < 1.6e-3 for x in g631),
       f"max={max(abs(x['err']) for x in g631) if g631 else 'NA':.2e}")
    out = {"script": "zero_cell_molecules.py", "checks": ROWS, "table": TBL,
           "n_pass": len(OK), "n_fail": len(BAD)}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n{len(OK)} PASS / {len(BAD)} FAIL  ->  {OUT}")
    return 0 if not BAD else 1


if __name__ == "__main__":
    raise SystemExit(main())
