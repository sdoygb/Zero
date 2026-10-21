#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""H2O 混合：轨道对胞腔 + SVD 低模 + 3H 的联合 MLE。

三档：
  3H      层 one-hot + 每层 runs/zc（Z168 特征）
  cell    (空穴集, 粒子集) 胞腔 + 3H
  hybrid  cell + SVD 低模（3|4 空间分划 top-k 的支配模）+ 3H
基准：
  SVD k=10 截断（标准 Schmidt；5.69e-4 Ha）
符号：oracle（FCI 符号），只测幅度模型。
输出：results/zero_h2o_hybrid.json
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
OUT = ROOT / "results" / "zero_h2o_hybrid.json"
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
    print("zero_h2o_hybrid.py — H2O：胞腔 + SVD 低模 + 3H")
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
    # 3H
    f3 = np.zeros((ndet, 3 * H))
    for i in range(ndet):
        f3[i, h[i]] = 1.0
        f3[i, H + h[i]] = r[i]
        f3[i, 2 * H + h[i]] = z[i]
    # cell
    labels = []
    for sa, sb in dets:
        holes = (hf & ~sa) | (hf & ~sb)
        parts = ((sa | sb) & ~hf) & 0b1100000
        labels.append(holes | parts)
    labels = np.array(labels)
    uniq = np.unique(labels); C = len(uniq)
    lab_idx = {u: i for i, u in enumerate(uniq)}
    idx_arr = np.array([lab_idx[l] for l in labels])
    fc = np.zeros((ndet, C + 3 * H))
    for i in range(ndet):
        fc[i, lab_idx[labels[i]]] = 1.0
        fc[i, C + h[i]] = 1.0
        fc[i, C + H + h[i]] = r[i]
        fc[i, C + 2 * H + h[i]] = z[i]
    # SVD 3|4 partition
    A = [(sa & 7, sb & 7) for sa, sb in dets]
    B = [(sa >> 3, sb >> 3) for sa, sb in dets]
    uA = sorted(set(A)); uB = sorted(set(B))
    iA = {x: i for i, x in enumerate(uA)}; iB = {x: i for i, x in enumerate(uB)}
    M = np.zeros((len(uA), len(uB)))
    for i in range(ndet):
        M[iA[A[i]], iB[B[i]]] = cc[i]
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    k = 10
    proj = np.zeros((ndet, k))
    for i in range(ndet):
        proj[i] = U[iA[A[i]], :k] * S[:k] * Vt[:k, iB[B[i]]]
    dom = np.argmax(np.abs(proj), axis=1)
    # hybrid: cell + mode + 3H
    fh = np.zeros((ndet, C + k + 3 * H))
    fh[:, :C] = fc[:, :C]
    for i in range(ndet):
        fh[i, C + dom[i]] = 1.0
        fh[i, C + k + h[i]] = 1.0
        fh[i, C + k + H + h[i]] = r[i]
        fh[i, C + k + 2 * H + h[i]] = z[i]
    p3 = mle(f3, w)
    pc = mle(fc, w)
    ph = mle(fh, w)
    def F(p): return float((g @ np.sqrt(p)) ** 2)
    def E(p):
        v = sgn * np.sqrt(p); v = v / np.linalg.norm(v); v = v.reshape(c.shape)
        return float(cis.energy(h1e, h2e, v, norb, nelec) + ecore)
    # SVD baseline
    Mk = (U[:, :k] * S[:k]) @ Vt[:k, :]
    csvd = np.zeros_like(c)
    nstr = len(strs)
    for i, (sa, sb) in enumerate(dets):
        ia = i // nstr; ib = i % nstr
        csvd[ia, ib] = Mk[iA[A[i]], iB[B[i]]]
    csvd /= np.linalg.norm(csvd)
    Esvd = cis.energy(h1e, h2e, csvd, norb, nelec) + ecore
    Fsvd = float((g @ (np.abs(csvd.reshape(-1)) / np.linalg.norm(csvd.reshape(-1)))) ** 2)
    print(f"    SVD k=10: E={Esvd:.8f} err={Esvd-E0:+.6f} F={Fsvd:.5f}")
    # --- state-level hybrid: SVD state x cell weight ratio ---
    Pcell = np.zeros(C); Psvd = np.zeros(C)
    cv = csvd.reshape(-1)
    for i in range(ndet):
        ci = lab_idx[labels[i]]
        Pcell[ci] += pc[i]
        Psvd[ci] += float(abs(cv[i]) ** 2)
    ratio = Pcell / np.maximum(Psvd, 1e-30)
    vh = cv * np.sqrt(ratio[idx_arr])
    vh = vh / np.linalg.norm(vh); vh = vh.reshape(c.shape)
    Esh = cis.energy(h1e, h2e, vh, norb, nelec) + ecore
    Fsh = float((g @ (np.abs(vh.reshape(-1)) / np.linalg.norm(vh.reshape(-1)))) ** 2)
    # --- spin-resolved cell labels (alpha/beta hole/particle counts) ---
    labels2 = []
    for i, (sa, sb) in enumerate(dets):
        nh_a = bin((hf & ~sa) & 0xFF).count("1"); nh_b = bin((hf & ~sb) & 0xFF).count("1")
        np_a = bin((sa & ~hf) & 0xFF).count("1"); np_b = bin((sb & ~hf) & 0xFF).count("1")
        labels2.append(int(labels[i]) * 6 ** 4 + nh_a * 6 ** 3 + nh_b * 6 ** 2 + np_a * 6 + np_b)
    labels2 = np.array(labels2)
    uniq2 = np.unique(labels2); C2 = len(uniq2)
    idx2 = {u: i for i, u in enumerate(uniq2)}
    f2 = np.zeros((ndet, C2 + 3 * H))
    for i in range(ndet):
        f2[i, idx2[labels2[i]]] = 1.0
        f2[i, C2 + h[i]] = 1.0
        f2[i, C2 + H + h[i]] = r[i]
        f2[i, C2 + 2 * H + h[i]] = z[i]
    p2 = mle(f2, w)
    # state-level hybrid with spin-resolved cells
    Pc2 = np.zeros(C2); Ps2 = np.zeros(C2)
    for i in range(ndet):
        ci = idx2[labels2[i]]
        Pc2[ci] += p2[i]
        Ps2[ci] += float(abs(cv[i]) ** 2)
    r2 = Pc2 / np.maximum(Ps2, 1e-30)
    idx2_arr = np.array([idx2[l] for l in labels2])
    vh2 = cv * np.sqrt(r2[idx2_arr])
    vh2 = vh2 / np.linalg.norm(vh2); vh2 = vh2.reshape(c.shape)
    Esh2 = cis.energy(h1e, h2e, vh2, norb, nelec) + ecore
    Fsh2 = float((g @ (np.abs(vh2.reshape(-1)) / np.linalg.norm(vh2.reshape(-1)))) ** 2)
    print(f"    state-hybrid(cell): err={Esh-E0:+.6f} F={Fsh:.5f}")
    print(f"    spin-cell: F={F(p2):.5f} err={E(p2)-E0:+.6f} params={C2+3*H} C2={C2}")
    print(f"    state-hybrid(spin-cell): err={Esh2-E0:+.6f} F={Fsh2:.5f}")
    print(f"    3H: F={F(p3):.5f} err={E(p3)-E0:+.6f} params={3*H}")
    print(f"    cell: F={F(pc):.5f} err={E(pc)-E0:+.6f} params={C+3*H}")
    print(f"    hybrid: F={F(ph):.5f} err={E(ph)-E0:+.6f} params={C+k+3*H}")
    ck("H1 FCI matches", abs(E0 + 75.01243743249299) < 1e-6, f"{E0:.8f}")
    ck("H2 SVD k=10 base 5.69e-4", abs((Esvd - E0) - 5.6909e-4) < 5e-5, f"{Esvd-E0:+.2e}")
    ck("H3 hybrid F > cell F", F(ph) > F(pc), f"{F(ph):.5f} vs {F(pc):.5f}")
    ck("H4 hybrid err < cell err", abs(E(ph) - E0) < abs(E(pc) - E0),
       f"{E(ph)-E0:+.6f} vs {E(pc)-E0:+.6f}")
    ck("H5 hybrid err < 1.6e-3", abs(E(ph) - E0) < 1.6e-3, f"{E(ph)-E0:+.2e}")
    ck("H6 params poly", C + k + 3 * H < 200, f"{C+k+3*H}")
    ck("H7 spin-cell F > cell F", F(p2) > F(pc), f"{F(p2):.5f} vs {F(pc):.5f}")
    ck("H8 spin-cell err < cell err", abs(E(p2) - E0) < abs(E(pc) - E0),
       f"{E(p2)-E0:+.6f} vs {E(pc)-E0:+.6f}")
    ck("H9 state-hybrid(cell) err <= SVD+1e-5", abs(Esh - E0) <= abs(Esvd - E0) + 1e-5,
       f"{Esh-E0:+.6f} vs SVD {Esvd-E0:+.6f}")
    ck("H10 state-hybrid(spin) err <= SVD+1e-5", abs(Esh2 - E0) <= abs(Esvd - E0) + 1e-5,
       f"{Esh2-E0:+.6f} vs SVD {Esvd-E0:+.6f}")
    out = {"script": "zero_h2o_hybrid.py", "checks": ROWS, "n_pass": len(OK), "n_fail": len(BAD),
           "e_fci": float(E0), "e_svd": float(Esvd), "err_svd": float(Esvd - E0),
           "F_svd": Fsvd, "F3": F(p3), "err3": float(E(p3) - E0),
           "Fc": F(pc), "errc": float(E(pc) - E0),
           "Fh": F(ph), "errh": float(E(ph) - E0),
           "C": C, "k": k, "H": H, "params3": 3 * H, "paramsc": C + 3 * H,
           "paramsh": C + k + 3 * H,
           "Esh": float(Esh), "Fsh": Fsh, "Esh2": float(Esh2), "Fsh2": Fsh2,
           "F2": F(p2), "err2": float(E(p2) - E0), "C2": C2, "params2": C2 + 3 * H}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n{len(OK)} PASS / {len(BAD)} FAIL  ->  {OUT}")
    return 0 if not BAD else 1


if __name__ == "__main__":
    raise SystemExit(main())
