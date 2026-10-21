#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""H2O 观测量：S4 native 符号态（自旋胞腔幅度 + 激发符号×自旋胞腔多数）的偶极/1-RDM。

对照：FCI 的偶极/1-RDM。
输出：results/zero_h2o_observables.json
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
OUT = ROOT / "results" / "zero_h2o_observables.json"
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


def main():
    print("=" * 74)
    print("zero_h2o_observables.py — H2O S4 态的偶极/1-RDM")
    print("=" * 74)
    mol, mf, h1e, h2e, ecore, cis, e0, c, strs, norb, nelec = build()
    E0 = e0 + ecore
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
    labels2 = []
    hfbits = (1 << nocc) - 1
    for i, (sa, sb) in enumerate(dets):
        holes = (hfbits & ~sa) | (hfbits & ~sb)
        parts = ((sa | sb) & ~hfbits) & 0b1111111
        nh_a = bin((hfbits & ~sa) & 0xFF).count("1"); nh_b = bin((hfbits & ~sb) & 0xFF).count("1")
        np_a = bin((sa & ~hfbits) & 0xFF).count("1"); np_b = bin((sb & ~hfbits) & 0xFF).count("1")
        labels2.append((holes | parts) * 6 ** 4 + nh_a * 6 ** 3 + nh_b * 6 ** 2 + np_a * 6 + np_b)
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
    # 1-RDM and dipole
    dm_v = cis.make_rdm1(vr, norb, nelec)
    dm_c = cis.make_rdm1(c, norb, nelec)
    dip_ao = mol.intor("int1e_r", comp=3)
    mo = mf.mo_coeff
    dip_mo = np.einsum("pi,xpq,qj->xij", mo, dip_ao, mo)
    nuc = np.einsum("i,ix->x", mol.atom_charges(), mol.atom_coords())
    el_v = np.array([np.trace(dm_v @ dip_mo[x]) for x in range(3)])
    el_c = np.array([np.trace(dm_c @ dip_mo[x]) for x in range(3)])
    dip_v = (nuc - el_v) * 2.541746
    dip_c = (nuc - el_c) * 2.541746
    d_dip = float(np.linalg.norm(dip_v - dip_c))
    d_rdm = float(np.linalg.norm(dm_v - dm_c))
    print(f"    E(S4)={Ev:.8f} err={Ev-E0:+.6f}")
    print(f"    dip FCI {np.round(dip_c,4)} D  S4 {np.round(dip_v,4)} D  |diff| {d_dip:.4f} D")
    print(f"    |1-RDM diff|_F {d_rdm:.4f}  trace {np.trace(dm_v):.6f}")
    ck("O1 FCI dipole 1.6166 D", abs(dip_c[1] - 1.6166) < 5e-3, f"{dip_c[1]:.4f}")
    ck("O2 S4 energy 7.89e-4", abs((Ev - E0) - 7.89e-4) < 5e-5, f"{Ev-E0:+.2e}")
    ck("O3 dipole diff < 0.05 D", d_dip < 0.05, f"{d_dip:.4f} D")
    ck("O4 1-RDM diff < 0.01", d_rdm < 0.01, f"{d_rdm:.4f}")
    ck("O5 trace = 10", abs(np.trace(dm_v) - nelec) < 1e-8, f"{np.trace(dm_v):.6f}")
    out = {"script": "zero_h2o_observables.py", "checks": ROWS,
           "e_fci": float(E0), "e_s4": float(Ev), "err_s4": float(Ev - E0),
           "dip_fci": dip_c.tolist(), "dip_s4": dip_v.tolist(),
           "dip_diff": d_dip, "rdm_diff": d_rdm, "C2": C2, "H": H,
           "n_pass": len(OK), "n_fail": len(BAD)}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n{len(OK)} PASS / {len(BAD)} FAIL  ->  {OUT}")
    return 0 if not BAD else 1


if __name__ == "__main__":
    raise SystemExit(main())
