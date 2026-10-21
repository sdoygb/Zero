#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""H2O 6-31G 状态级混合：SVD 低模（4|9 分划） × 自旋胞腔权重，S4 符号 + SVD 幅度。

基准：cell 18 mHa；目标 <1.6e-3。
输出：results/zero_h2o_631g_hybrid.json
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
OUT = ROOT / "results" / "zero_h2o_631g_hybrid.json"
OK, BAD, ROWS = [], [], []


def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    ROWS.append({"name": name, "ok": bool(cond), "detail": detail})
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def sstr(hf_s, d_s, norb):
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


def main():
    print("=" * 74)
    print("zero_h2o_631g_hybrid.py — H2O 6-31G 状态级混合")
    print("=" * 74)
    atom = "O 0 0 0; H 0.757 0.586 0; H -0.757 0.586 0"
    mol = gto.M(atom=atom, basis="6-31g", unit="Angstrom", verbose=0)
    mf = scf.RHF(mol).run()
    norb, nelec = mol.nao_nr(), mol.nelectron
    nocc = nelec // 2
    h1e = np.einsum("pi,pq,qj->ij", mf.mo_coeff,
                    mol.intor("int1e_kin") + mol.intor("int1e_nuc"), mf.mo_coeff)
    h2e = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff), norb)
    ecore = mol.energy_nuc()
    cis = fd.FCI()
    e0, c = cis.kernel(h1e, h2e, norb, nelec)
    E0 = e0 + ecore
    strs = fci.cistring.make_strings(range(norb), nocc)
    nstr = len(strs); ndet = c.size
    print(f"    FCI {E0:.8f}  ndet {ndet}")
    cc = c.reshape(-1)
    g = np.abs(cc); g /= np.linalg.norm(g)
    w = g ** 2; w /= w.sum()
    sgn = np.sign(cc)
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
    labels = ((hf & ~sa) | (hf & ~sb)) | ((sa | sb) & ~hf)
    uniq = np.unique(labels); C = len(uniq)
    nh_a = np.array([bin((hf & ~x) & 0xFFFF).count("1") for x in sa])
    nh_b = np.array([bin((hf & ~x) & 0xFFFF).count("1") for x in sb])
    np_a = np.array([bin((x & ~hf) & 0xFFFF).count("1") for x in sa])
    np_b = np.array([bin((x & ~hf) & 0xFFFF).count("1") for x in sb])
    labels2 = labels.astype(np.int64) * 6 ** 4 + nh_a * 6 ** 3 + nh_b * 6 ** 2 + np_a * 6 + np_b
    uniq2 = np.unique(labels2); C2 = len(uniq2)
    idx2 = np.searchsorted(uniq2, labels2)
    rows = np.repeat(np.arange(ndet), 4)
    cols = np.empty(4 * ndet, dtype=np.int64)
    cols[0::4] = idx2; cols[1::4] = C2 + h; cols[2::4] = C2 + H + h; cols[3::4] = C2 + 2 * H + h
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
    zz = f @ res.x; p = np.exp(zz - zz.max()); p /= p.sum()
    print(f"    MLE nit={res.nit}")
    # SVD 4|9 partition
    A = [(x & 0b1111, y & 0b1111) for x, y in zip(sa, sb)]
    B = [(x >> 4, y >> 4) for x, y in zip(sa, sb)]
    uA = sorted(set(A)); uB = sorted(set(B))
    iA = {x: i for i, x in enumerate(uA)}; iB = {x: i for i, x in enumerate(uB)}
    print(f"    A dim {len(uA)}  B dim {len(uB)}")
    M = np.zeros((len(uA), len(uB)))
    for i in range(ndet):
        M[iA[A[i]], iB[B[i]]] = cc[i]
    MMt = M @ M.T
    ev, U = np.linalg.eigh(MMt)
    order = np.argsort(-ev)[:min(200, len(ev))]
    Sv = np.sqrt(np.maximum(ev[order], 0.0))
    print(f"    top Schmidt {np.round(Sv[:6],4)}")
    def state_from_k(k):
        Uk = U[:, order[:k]]
        Mk = Uk @ (Uk.T @ M)
        ck = np.zeros(ndet)
        for i in range(ndet):
            ck[i] = Mk[iA[A[i]], iB[B[i]]]
        ck /= np.linalg.norm(ck)
        return ck
    def energy_vec(v):
        vv = v / np.linalg.norm(v); vv = vv.reshape(c.shape)
        return float(cis.energy(h1e, h2e, vv, norb, nelec) + ecore)
    # S4 sign model
    ds = sstr(hf, sa, norb) * sstr(hf, sb, norb)
    res_sign = sgn * ds
    sw = np.bincount(idx2, weights=w * res_sign, minlength=C2)
    spin_sign = np.where(sw >= 0, 1.0, -1.0)
    model = ds * spin_sign[idx2]
    Pcell = np.bincount(idx2, weights=p, minlength=C2)
    results = []
    for kk in [10, 20, 50, 100]:
        csvd = state_from_k(kk)
        Ek = energy_vec(csvd)
        Psvd = np.bincount(idx2, weights=csvd ** 2, minlength=C2)
        ratio = Pcell / np.maximum(Psvd, 1e-30)
        vh = csvd * np.sqrt(ratio[idx2])
        Eh = energy_vec(vh)
        vs4 = model * np.abs(csvd)
        Es4 = energy_vec(vs4)
        results.append({"k": kk, "SVD": Ek, "hybrid": Eh, "S4xSVD": Es4})
        print(f"    k={kk:3d}: SVD err={Ek-E0:+.6f}  hybrid err={Eh-E0:+.6f}  S4xSVD err={Es4-E0:+.6f}")
    allv = [(r["SVD"], f"SVD k={r['k']}") for r in results] + \
           [(r["hybrid"], f"hybrid k={r['k']}") for r in results] + \
           [(r["S4xSVD"], f"S4xSVD k={r['k']}") for r in results]
    best = min(allv, key=lambda x: abs(x[0] - E0))
    print(f"    best {best[1]} err={best[0]-E0:+.6f}")
    ck("K1 FCI 6-31G", abs(E0 + 76.12084479) < 1e-5, f"{E0:.8f}")
    ck("K2 SVD k=10 err < 1.6e-3", abs(results[0]["SVD"] - E0) < 1.6e-3, f"{results[0]['SVD']-E0:+.2e}")
    ck("K3 best SVD err < 1.6e-3", abs(min(r["SVD"] for r in results) - E0) < 1.6e-3,
       f"{min(r['SVD'] for r in results)-E0:+.2e}")
    ck("K4 best hybrid err < 1.6e-3", abs(min(r["hybrid"] for r in results) - E0) < 1.6e-3,
       f"{min(r['hybrid'] for r in results)-E0:+.2e}")
    ck("K5 best overall err < 1.6e-3", abs(best[0] - E0) < 1.6e-3, f"{best[1]} {best[0]-E0:+.2e}")
    out = {"script": "zero_h2o_631g_hybrid.py", "checks": ROWS, "n_pass": len(OK), "n_fail": len(BAD),
           "e_fci": float(E0), "results": results, "best": best[1],
           "best_err": float(best[0] - E0), "A": len(uA), "B": len(uB), "C2": C2}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n{len(OK)} PASS / {len(BAD)} FAIL  ->  {OUT}")
    return 0 if not BAD else 1


if __name__ == "__main__":
    raise SystemExit(main())
