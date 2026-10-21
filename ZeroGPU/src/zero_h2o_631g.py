#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""H2O 6-31G：FCI(1.66M) + 自旋胞腔 + 稀疏 MLE + S4 native 符号 + 偶极/1-RDM。

输出：results/zero_h2o_631g.json
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.sparse import csr_matrix

try:
    from pyscf import gto, scf, fci, ao2mo
    import pyscf.fci.direct_spin1 as fd
except Exception as exc:
    raise SystemExit("需要 pyscf（/usr/local/bin/python3）：%r" % (exc,))

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "zero_h2o_631g.json"
OK, BAD, ROWS = [], [], []


def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    ROWS.append({"name": name, "ok": bool(cond), "detail": detail})
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def main():
    print("=" * 74)
    print("zero_h2o_631g.py — H2O 6-31G 全链")
    print("=" * 74)
    atom = "O 0 0 0; H 0.757 0.586 0; H -0.757 0.586 0"
    mol = gto.M(atom=atom, basis="6-31g", unit="Angstrom", verbose=0)
    mf = scf.RHF(mol).run()
    norb, nelec = mol.nao_nr(), mol.nelectron
    nocc = nelec // 2
    print(f"    nao {norb}  nelec {nelec}  ndet {len(fci.cistring.make_strings(range(norb), nocc)) ** 2}")
    h1e = np.einsum("pi,pq,qj->ij", mf.mo_coeff,
                    mol.intor("int1e_kin") + mol.intor("int1e_nuc"), mf.mo_coeff)
    h2e = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff), norb)
    ecore = mol.energy_nuc()
    cis = fd.FCI()
    e0, c = cis.kernel(h1e, h2e, norb, nelec)
    E0 = e0 + ecore
    strs = fci.cistring.make_strings(range(norb), nocc)
    nstr = len(strs)
    ndet = c.size
    print(f"    FCI {E0:.8f}  c {c.shape}")
    cc = c.reshape(-1)
    g = np.abs(cc); g /= np.linalg.norm(g)
    w = g ** 2; w /= w.sum()
    sgn = np.sign(cc)
    # vectorized features: interleaved word
    sa = np.repeat(strs, nstr); sb = np.tile(strs, nstr)
    bits = np.empty((ndet, 2 * norb), dtype=np.int8)
    for i in range(norb):
        bits[:, 2 * i] = (sa >> i) & 1
        bits[:, 2 * i + 1] = (sb >> i) & 1
    S = np.cumsum(2 * bits - 1, axis=1)
    h = np.max(np.abs(S), axis=1)
    r = 1 + np.sum(bits[:, 1:] != bits[:, :-1], axis=1)
    z = np.sum(S[:, 1:-1] == 0, axis=1)
    H = int(h.max()) + 1
    hf = (1 << nocc) - 1
    holes = (hf & ~sa) | (hf & ~sb)
    parts = (sa | sb) & ~hf
    labels = holes | parts
    uniq = np.unique(labels); C = len(uniq)
    lab_idx = np.searchsorted(uniq, labels)
    nh_a = np.array([bin((hf & ~x) & 0xFFFF).count("1") for x in sa])
    nh_b = np.array([bin((hf & ~x) & 0xFFFF).count("1") for x in sb])
    np_a = np.array([bin((x & ~hf) & 0xFFFF).count("1") for x in sa])
    np_b = np.array([bin((x & ~hf) & 0xFFFF).count("1") for x in sb])
    labels2 = labels.astype(np.int64) * 6 ** 4 + nh_a * 6 ** 3 + nh_b * 6 ** 2 + np_a * 6 + np_b
    uniq2 = np.unique(labels2); C2 = len(uniq2)
    idx2 = np.searchsorted(uniq2, labels2)
    print(f"    C={C} C2={C2} H={H}")
    # sparse features
    rows = np.repeat(np.arange(ndet), 4)
    cols = np.empty(4 * ndet, dtype=np.int64)
    cols[0::4] = idx2
    cols[1::4] = C2 + h
    cols[2::4] = C2 + H + h
    cols[3::4] = C2 + 2 * H + h
    f = csr_matrix((np.ones(4 * ndet), (rows, cols)), shape=(ndet, C2 + 3 * H))
    def obj(theta):
        zz = f @ theta
        return logsumexp(zz) - np.dot(w, zz)
    def grad(theta):
        zz = f @ theta
        p = np.exp(zz - zz.max()); p /= p.sum()
        return np.asarray(f.T @ (p - w)).ravel()
    res = minimize(obj, np.zeros(C2 + 3 * H), jac=grad, method="L-BFGS-B",
                   options={"maxiter": 4000, "ftol": 1e-14})
    zz = f @ res.x
    p = np.exp(zz - zz.max()); p /= p.sum()
    print(f"    MLE nll={res.fun:.6f} nit={res.nit}")
    # S4 sign: excitation sign x spin-cell majority
    def sstr(hf_s, d_s):
        holes_ = hf_s & ~d_s
        parts_ = d_s & ~hf_s
        sign = np.ones(len(d_s))
        for hh in range(norb):
            mh = ((holes_ >> hh) & 1).astype(bool)
            for pp in range(norb):
                if hh < pp:
                    mp = ((parts_ >> pp) & 1).astype(bool)
                    sign[mh & mp] *= -1.0
        return sign
    ds = sstr(hf, sa) * sstr(hf, sb)
    res_sign = sgn * ds
    sw = np.bincount(idx2, weights=w * res_sign, minlength=C2)
    spin_sign = np.where(sw >= 0, 1.0, -1.0)
    model = ds * spin_sign[idx2]
    v = model * np.sqrt(p); v = v / np.linalg.norm(v); vr = v.reshape(c.shape)
    Ev = cis.energy(h1e, h2e, vr, norb, nelec) + ecore
    acc = float(np.sum(w * (model == sgn)))
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
    params = C2 + 3 * H
    print(f"    S4: err={Ev-E0:+.6f} acc={acc:.4f} F={float((cc/np.linalg.norm(cc))@v)**2:.5f}")
    print(f"    dip FCI {np.round(dip_c,4)} / S4 {np.round(dip_v,4)}  |diff| {d_dip:.4f} D")
    print(f"    |1-RDM| {d_rdm:.4f}  params {params}  ratio {ndet/params:.1f}")
    ck("G1 FCI 6-31G ran", ndet > 1_000_000, f"ndet={ndet}")
    ck("G2 S4 err < 1.6e-3", abs(Ev - E0) < 1.6e-3, f"{Ev-E0:+.2e}")
    ck("G3 dipole diff < 0.05 D", d_dip < 0.05, f"{d_dip:.4f}")
    ck("G4 1-RDM diff < 0.01", d_rdm < 0.01, f"{d_rdm:.4f}")
    ck("G5 params 30440 (poly < ndet)", params < ndet, f"{params} < {ndet}")
    out = {"script": "zero_h2o_631g.py", "checks": ROWS, "n_pass": len(OK), "n_fail": len(BAD),
           "nao": norb, "nelec": nelec, "ndet": ndet, "E_fci": float(E0),
           "E_s4": float(Ev), "err": float(Ev - E0), "acc": acc,
           "C": C, "C2": C2, "H": H, "params": params,
           "dip_fci": dip_c.tolist(), "dip_s4": dip_v.tolist(), "dip_diff": d_dip,
           "rdm_diff": d_rdm}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n{len(OK)} PASS / {len(BAD)} FAIL  ->  {OUT}")
    return 0 if not BAD else 1


if __name__ == "__main__":
    raise SystemExit(main())
