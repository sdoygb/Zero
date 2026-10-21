#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""H2O 首算：STO-3G FCI 基线、原生 3H 振幅拟合、纠缠/模隙与模截断。

运行（本机 pyscf 在 python.org 3.12 下）：
    /usr/local/bin/python3 ZeroGPU/src/zero_h2o.py

内容
  B. 基线：RHF/FCI 能量、关联能、HOMO-LUMO、偶极（RHF/FCI）。
  N. 原生 3H 振幅（层 h、段数 r、零点穿越 z）对 FCI 行列式权重的 MLE 拟合；
     用 FCI 符号（oracle）算能量误差。→ 粗粒化特征不够。
  E. 空间 3|4 分划的 Schmidt 谱、模隙、k 模截断能量。→ 10 模到化学精度。
输出：results/zero_h2o.json
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
from scipy.optimize import minimize

try:
    from pyscf import gto, scf, fci, ao2mo
    import pyscf.fci.direct_spin1 as fd
except Exception as exc:  # pragma: no cover
    raise SystemExit("需要 pyscf（/usr/local/bin/python3）：%r" % (exc,))

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "zero_h2o.json"
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


def feats(sa, sb, norb, interleaved=True):
    if interleaved:
        bits = []
        for i in range(norb):
            bits += [(sa >> i) & 1, (sb >> i) & 1]
    else:
        bits = [(sa >> i) & 1 for i in range(norb)] + [(sb >> i) & 1 for i in range(norb)]
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


def fit_3h(dets, w, norb, interleaved=True, with_rank=False):
    F = np.array([feats(d[0], d[1], norb, interleaved) for d in dets])
    h, r, z = F[:, 0].astype(int), F[:, 1].astype(int), F[:, 2].astype(int)
    H = h.max() + 1
    if with_rank:
        hf = (1 << (max(1, norb // 2))) - 1  # not used; rank feature below
        e = np.array([bin(d[0] & ~((1 << 5) - 1)).count("1") +
                      bin(d[1] & ~((1 << 5) - 1)).count("1") for d in dets])
        E = e.max() + 1
    else:
        e = np.zeros(len(dets), dtype=int)
        E = 0
    npar = 3 * H + E

    def unpack(x):
        u = x[:H]; mu = x[H:2 * H]; nu = x[2 * H:3 * H]
        u = u - u.max(); ph = np.exp(u); ph /= ph.sum()
        lw = mu[h] * r + nu[h] * z
        if E:
            lw = lw + x[3 * H:3 * H + E][e]
        lz = np.zeros(H)
        for k in range(H):
            m = h == k
            if m.any():
                lz[k] = np.log(np.exp(lw[m] - lw[m].max()).sum()) + lw[m].max()
        return np.log(ph[h] + 1e-300) + lw - lz[h], ph

    def nll(x):
        return -np.dot(w, unpack(x)[0])

    def grad(x):
        lp, ph = unpack(x)
        p = np.exp(lp); g = np.zeros_like(x)
        for k in range(H):
            m = h == k
            if not m.any():
                continue
            W, P = w[m].sum(), p[m].sum()
            g[k] = P - W
            g[H + k] = -(np.dot(w[m], r[m]) - W * np.dot(p[m] / max(P, 1e-300), r[m]))
            g[2 * H + k] = -(np.dot(w[m], z[m]) - W * np.dot(p[m] / max(P, 1e-300), z[m]))
        if E:
            for a in range(E):
                m = e == a
                g[3 * H + a] = p[m].sum() - w[m].sum()
        return g

    best = None
    for x0 in [np.zeros(npar), np.random.default_rng(0).normal(0, 0.2, npar)]:
        res = minimize(nll, x0, jac=grad, method="L-BFGS-B",
                       options={"maxiter": 4000, "ftol": 1e-15})
        if best is None or res.fun < best.fun:
            best = res
    lp, ph = unpack(best.x)
    return np.exp(lp), float(np.sum(np.sqrt(w * np.exp(lp))) ** 2), npar, H


def main():
    print("=" * 72)
    print("zero_h2o.py — H2O STO-3G：FCI、原生 3H、纠缠/模截断")
    print("=" * 72)
    mol, mf, h1e, h2e, ecore, cis, e0, c, strs, norb, nelec = build()
    E0 = e0 + ecore
    print(f"    RHF {mf.e_tot:.8f}  FCI {E0:.8f}  corr {E0-mf.e_tot:.6f}")
    homo, lumo = mf.mo_energy[nelec // 2 - 1], mf.mo_energy[nelec // 2]
    print(f"    HOMO {homo:.4f}  LUMO {lumo:.4f}  gap {lumo-homo:.4f} Ha")
    dm1 = cis.make_rdm1(c, norb, nelec)
    dip_ao = mol.intor("int1e_r", comp=3)
    mo = mf.mo_coeff
    dip_mo = np.einsum("pi,xpq,qj->xij", mo, dip_ao, mo)
    el = np.array([np.trace(dm1 @ dip_mo[x]) for x in range(3)])
    nuc = np.einsum("i,ix->x", mol.atom_charges(), mol.atom_coords())
    total = (nuc - el) * 2.541746  # Debye
    print(f"    dipole FCI total {np.round(total,4)} D  (RHF {np.round(mf.dip_moment(verbose=0),4)})")

    # native 3H fits
    dets = [(sa, sb, ia, ib, c[ia, ib]) for ia, sa in enumerate(strs) for ib, sb in enumerate(strs)]
    w = np.array([abs(d[4]) ** 2 for d in dets]); w /= w.sum()
    signs = np.sign(np.array([d[4] for d in dets]))
    native = []
    for inter in [True, False]:
        for rank in [False, True]:
            p, fmag, npar, H = fit_3h(dets, w, norb, inter, rank)
            v = signs * np.sqrt(p); v = v.reshape(c.shape); v /= np.linalg.norm(v)
            E = cis.energy(h1e, h2e, v, norb, nelec) + ecore
            native.append({"interleaved": inter, "rank": rank, "H": H, "params": npar,
                           "Fmag": fmag, "E": float(E), "err": float(E - E0)})
            print(f"    native inter={inter} rank={rank}: H={H} params={npar} Fmag={fmag:.4f}"
                  f"  E={E:.8f}  err={E-E0:.6f}")

    # entanglement / Schmidt across spatial 3|4
    A = [(sa & 7, sb & 7) for sa in strs for sb in strs]
    B = [(sa >> 3, sb >> 3) for sa in strs for sb in strs]
    uA = sorted(set(A)); uB = sorted(set(B))
    iA = {x: i for i, x in enumerate(uA)}; iB = {x: i for i, x in enumerate(uB)}
    M = np.zeros((len(uA), len(uB)))
    for ia, sa in enumerate(strs):
        for ib, sb in enumerate(strs):
            M[iA[(sa & 7, sb & 7)], iB[(sa >> 3, sb >> 3)]] = c[ia, ib]
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    p = S ** 2
    gap = float(np.log(S[0] / S[1]))
    print(f"    Schmidt top {np.round(S[:6],4)}  entropy {-np.sum(p[p>0]*np.log(p[p>0])):.4f}"
          f"  modular gap {2*gap:.4f}")
    trunc = []
    for k in [1, 2, 3, 4, 5, 7, 10, 15, 20]:
        Mk = (U[:, :k] * S[:k]) @ Vt[:k, :]
        ckvec = np.zeros_like(c)
        for ia, sa in enumerate(strs):
            for ib, sb in enumerate(strs):
                ckvec[ia, ib] = Mk[iA[(sa & 7, sb & 7)], iB[(sa >> 3, sb >> 3)]]
        ckvec /= np.linalg.norm(ckvec)
        Ek = cis.energy(h1e, h2e, ckvec, norb, nelec) + ecore
        trunc.append({"k": k, "fid": float(np.sum(p[:k])), "E": float(Ek), "err": float(Ek - E0)})
        print(f"    k={k:2d}  fid={np.sum(p[:k]):.6f}  E={Ek:.8f}  err={Ek-E0:.6f}")

    ck("H1 FCI matches -75.012437", abs(E0 + 75.01243743249299) < 1e-6, f"E0={E0:.8f}")
    ck("H2 correlation", abs((E0 - mf.e_tot) + 0.049491) < 1e-4, f"corr={E0-mf.e_tot:.6f}")
    best_nat = min(native, key=lambda r: abs(r["err"]))
    ck("H3 native 3H coarse (Fmag>0.95, err>1.6e-3)",
       best_nat["Fmag"] > 0.95 and abs(best_nat["err"]) > 1.6e-3,
       f"Fmag={best_nat['Fmag']:.4f} err={best_nat['err']:.4f}")
    ck("H4 modular gap > 1", 2 * gap > 1.0, f"{2*gap:.3f}")
    k10 = [x for x in trunc if x["k"] == 10][0]
    ck("H5 10 modes chemical accuracy", abs(k10["err"]) < 1.6e-3, f"err={k10['err']:.2e}")
    ck("H6 3 modes fidelity > 0.99", trunc[2]["fid"] > 0.99, f"{trunc[2]['fid']:.4f}")

    out = {"script": "zero_h2o.py", "checks": ROWS, "n_pass": len(OK), "n_fail": len(BAD),
           "geometry": "O 0 0 0; H 0.757 0.586 0; H -0.757 0.586 0 (Angstrom)",
           "basis": "STO-3G", "norb": norb, "nelec": nelec,
           "e_rhf": float(mf.e_tot), "e_fci": float(E0), "corr": float(E0 - mf.e_tot),
           "homo": float(homo), "lumo": float(lumo), "gap_homo_lumo": float(lumo - homo),
           "dipole_fci_D": [float(x) for x in total],
           "native": [{**r, "H": int(r["H"]), "params": int(r["params"])} for r in native],
           "schmidt": S[:20].tolist(),
           "entropy": float(-np.sum(p[p > 0] * np.log(p[p > 0]))),
           "modular_gap": float(2 * gap),
           "truncation": [{**x, "k": int(x["k"]), "fid": float(x["fid"])} for x in trunc]}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n{len(OK)} PASS / {len(BAD)} FAIL  ->  {OUT}")
    return 0 if not BAD else 1


if __name__ == "__main__":
    raise SystemExit(main())
