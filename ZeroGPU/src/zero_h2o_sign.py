#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""H2O native 符号模型：行列式费米子符号 + 各组多数关联符号，替掉 oracle。

自旋胞腔幅度（151 胞腔 + 3H）的 MLE p2；符号：
  S0 oracle（FCI 符号，基准）
  S1 det       行列式符号（相对 HF 的置换奇偶）
  S2 det × rank 多数
  S3 det × 空间胞腔多数
  S4 det × 自旋胞腔多数
  S5 det × (空间胞腔×rank) 多数
输出：results/zero_h2o_sign.json
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
OUT = ROOT / "results" / "zero_h2o_sign.json"
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
    print("zero_h2o_sign.py — H2O native 符号模型")
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
    labels = []
    for sa, sb in dets:
        holes = (hf & ~sa) | (hf & ~sb)
        parts = ((sa | sb) & ~hf) & 0b1100000
        labels.append(holes | parts)
    labels = np.array(labels)
    uniq = np.unique(labels); C = len(uniq)
    lab_idx = {u: i for i, u in enumerate(uniq)}
    lab_arr = np.array([lab_idx[l] for l in labels])
    labels2 = []
    for i, (sa, sb) in enumerate(dets):
        nh_a = bin((hf & ~sa) & 0xFF).count("1"); nh_b = bin((hf & ~sb) & 0xFF).count("1")
        np_a = bin((sa & ~hf) & 0xFF).count("1"); np_b = bin((sb & ~hf) & 0xFF).count("1")
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
    k = np.array([bin((sa & ~hf) & 0x7F).count("1") + bin((sb & ~hf) & 0x7F).count("1")
                  for sa, sb in dets])
    K = int(k.max()) + 1
    rank_sign = majority(k, res, w, K)
    cell_sign = majority(lab_arr, res, w, C)
    spin_sign = majority(idx2_arr, res, w, C2)
    cr = lab_arr * K + k
    cr_sign = majority(cr, res, w, C * K)
    models = {
        "S0_oracle": sgn,
        "S1_det": ds,
        "S2_rank": ds * rank_sign[k],
        "S3_cell": ds * cell_sign[lab_arr],
        "S4_spin": ds * spin_sign[idx2_arr],
        "S5_cellrank": ds * cr_sign[cr],
    }
    c_hat = cc / np.linalg.norm(cc)
    rows = []
    for name, m in models.items():
        acc = float(np.sum(w * (m == sgn)))
        v = m * np.sqrt(p2); v = v / np.linalg.norm(v); vr = v.reshape(c.shape)
        E = float(cis.energy(h1e, h2e, vr, norb, nelec) + ecore)
        F = float((c_hat @ v) ** 2)
        rows.append({"name": name, "acc": acc, "E": E, "err": E - E0, "F": F})
        print(f"    {name:12s} acc={acc:.4f} E={E:.8f} err={E-E0:+.6f} F={F:.5f}")
    d1 = float(np.sum(w * (ds == sgn)))
    print(f"    det-sign weighted acc {d1:.4f}; residual sign +1 frac {np.sum(w*(res>0)):.4f}")
    ck("S0 oracle base 7.37e-4", abs((rows[0]["err"]) - 7.37e-4) < 5e-5, f"{rows[0]['err']:+.2e}")
    ck("S1 det-sign meaningful (acc>0.5)", d1 > 0.5, f"{d1:.4f}")
    ck("S4 spin-sign acc > S1 acc", rows[4]["acc"] > rows[1]["acc"],
       f"{rows[4]['acc']:.4f} vs {rows[1]['acc']:.4f}")
    best = min(rows[1:], key=lambda x: abs(x["err"]))
    ck("best native sign err < 1.6e-3", abs(best["err"]) < 1.6e-3, f"{best['name']} {best['err']:+.2e}")
    out = {"script": "zero_h2o_sign.py", "checks": ROWS, "rows": rows,
           "det_acc": d1, "res_pos_frac": float(np.sum(w * (res > 0))),
           "C": C, "C2": C2, "K": K, "params2": C2 + 3 * H,
           "n_pass": len(OK), "n_fail": len(BAD)}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n{len(OK)} PASS / {len(BAD)} FAIL  ->  {OUT}")
    return 0 if not BAD else 1


if __name__ == "__main__":
    raise SystemExit(main())
