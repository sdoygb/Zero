#!/usr/bin/env python3
"""χ 与镜面判据的一致性检查（记录层用）

命题（本脚本要核的那条）：若分子图 Γ 上有**非常规**自同构 π（镜面/反演，判据判非手性），
且相位场在 π 下按镜像共轭变换 z_{πv} = conj(z_v)，则 ∏_v z_v 为实数，故 χ = sgn Im ∏ z_v = 0。
反过来，判据判手性（无非真自同构）时，把随机相位场对**常规**自同构群做平均后，
Im ∏ z_v 一般不为零——χ ≠ 0 允许存在。

做法（纯数值，不依赖 {z_v} 的动力学）：
  1. 随机取零和相位场 z（Σ z_v = 0）；
  2. 对自同构群平均：常规 π 按 z 的原像搬，非常规 π 按共轭搬；
  3. 算 Im ∏ z，与判据（zb092_chirality_topology.py 的 rms 判据）对照。
"""
from __future__ import annotations

import json
import sys

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

sys.path.insert(0, "/Users/oygb/Downloads/Zero/ZeroGPU/src")
from zb092_chirality_topology import (  # noqa: E402
    TOL, automorphisms, conf_improper_rms, has_cumulated,
)

RDLogger.DisableLog("rdApp.*")

OUT = "/Users/oygb/Downloads/Zero/ZeroGPU/results/zb092_chi_consistency.json"
TRIALS = 64


def phase_product(mol: Chem.Mol, proper: list, improper: list, rng, trials=TRIALS):
    """对给定对称群（常规 proper + 非常规 improper）平均相位场，返回 (均|Im∏z|, 最大|Im∏z|)。

    只对**同一构象**里真正被实现的置换做平均：只有这样才能保证 ε(π) 是群同态
    （可挠分子的「非常规」自同构只在特定构象里成立，跨构象混用会破坏同态性）。
    """
    n = mol.GetNumAtoms()
    grp = [(pi, False) for pi in proper] + [(pi, True) for pi in improper]
    vals = []
    for _ in range(trials):
        z = rng.normal(size=n) + 1j * rng.normal(size=n)
        z -= z.mean()                              # 零和
        acc = np.zeros(n, dtype=complex)
        for pi, imp in grp:
            acc += np.conj(z)[np.asarray(pi)] if imp else z[np.asarray(pi)]
        acc /= len(grp)
        vals.append(abs(float(np.imag(np.prod(acc)))))
    return float(np.mean(vals)), float(np.max(vals))


def _proper_rms(X, pi):
    """把 X[π] 拟合到 X，只允许真转动时的残差。"""
    Xq = X[pi]
    Xqc = Xq - Xq.mean(0)
    Xc = X - X.mean(0)
    U, S, Vt = np.linalg.svd(Xc.T @ Xqc)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt2 = Vt.copy()
        Vt2[-1] *= -1
        R = Vt2.T @ U.T
    return float(np.sqrt(((Xc @ R.T - Xqc) ** 2).sum(1).mean())), float(np.linalg.det(R))


CASES = [
    ("甲烷", "C"), ("乙醇", "CCO"), ("甘氨酸", "NCC(=O)O"),
    ("丙氨酸（L）", "C[C@@H](N)C(=O)O"), ("乳酸（L）", "C[C@@H](O)C(=O)O"),
    ("苯", "c1ccccc1"), ("内消旋酒石酸", "O=C(O)[C@@H](O)[C@@H](O)C(=O)O"),
    ("(S,S)-酒石酸", "O=C(O)[C@@H](O)[C@H](O)C(=O)O"),
    ("顺-1,2-二甲基环丙烷（内消旋）", "C[C@@H]1C[C@@H]1C"),
    ("反-1,2-二甲基环丙烷（R,R）", "C[C@@H]1C[C@H]1C"),
    ("柠檬烯（R）", "CC(=C)[C@H]1CCC(=CC1)C"),
    ("樟脑", "CC1(C)C2CC[C@@]1(C)C(=O)C2"),
    ("环己烷", "C1CCCCC1"),
    ("2,3-戊二烯（轴手性）", "CC=C=CC"),
    ("取代丙二烯 ClC=C=C(F)Br", "ClC=C=C(F)Br"),
    ("3-甲基戊烷", "CCC(C)CC"),
    ("1,3-二氯环戊烷（顺式）", "Cl[C@H]1C[C@@H](Cl)CC1"),
    ("葡萄糖（开链 D）", "OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C=O"),
]


def main() -> int:
    rng = np.random.default_rng(20261109)
    rows = []
    for name, smi in CASES:
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        auts = automorphisms(mol)
        ps = AllChem.ETKDGv3()
        ps.randomSeed = 0xF00D
        ps.pruneRmsThresh = -1.0
        cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=12, params=ps))
        cum = has_cumulated(Chem.RemoveHs(mol))
        if not cum:
            for cid in cids:
                try:
                    AllChem.MMFFOptimizeMolecule(mol, confId=cid, maxIters=300)
                except Exception:
                    pass
        # 逐构象求「真正被实现的对称群」：常规集合 + 非常规集合
        best_key, best_conf, best_imp, best_prop = (-1, -1), None, [], []
        for cid in cids:
            X = np.asarray(mol.GetConformer(cid).GetPositions(), dtype=np.float64)
            imp = [pi for pi in auts if conf_improper_rms(X, [pi])[0] < TOL]
            prop = [pi for pi in auts if _proper_rms(X, np.asarray(pi))[0] < TOL]
            key = (len(imp), len(prop))
            if key > best_key:
                best_key, best_conf, best_imp, best_prop = key, cid, imp, prop
        chiral = not bool(best_imp)
        mean_im, max_im = phase_product(mol, best_prop, best_imp, rng)
        rows.append({"name": name, "smiles": smi, "aut": len(auts),
                     "n_proper": len(best_prop), "n_improper": len(best_imp),
                     "chiral": chiral,
                     "mean_abs_im_prod": mean_im, "max_abs_im_prod": max_im})
        flag = "✓" if (chiral and max_im > 1e-6) or ((not chiral) and max_im < 1e-9) else "✗"
        print(f"{name:<26} 自同构{len(auts):>5} 群{len(best_prop)}+{len(best_imp):<4} "
              f"判定 {'手性' if chiral else '非手性':<4}  |Im∏z| 均值 {mean_im:.2e} 最大 {max_im:.2e}  {flag}")

    ok = sum(1 for r in rows if (r["chiral"] and r["max_abs_im_prod"] > 1e-6)
             or ((not r["chiral"]) and r["max_abs_im_prod"] < 1e-9))
    print(f"\n一致 {ok}/{len(rows)}：判非手性 ⟹ 对称化后 |Im∏z| ≈ 0；判手性 ⟹ 不为零")
    with open(OUT, "w") as f:
        json.dump({"tol": TOL, "trials": TRIALS, "n": len(rows), "consistent": ok, "rows": rows},
                  f, indent=2, ensure_ascii=False)
    print("落盘", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
