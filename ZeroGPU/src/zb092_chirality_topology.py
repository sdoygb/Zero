#!/usr/bin/env python3
"""手性的拓扑判据：自同构群里的镜面元素（v2，多构象）

问题（`_staging/推导/零宇宙分子_电子手征.md` §6 第一步）：
「骨架无镜面对称 ⟺ 允许 χ≠0」能不能写成**可枚举判据**，并在真实分子上验？

判据：分子（带元素标签的图 Γ，含显式氢）**非手性** ⟺ 存在一个构象 ρ 和一个自同构 π，
使得 π 被一个**非常规**等距（镜面/反演，det=−1）实现：ρ(π v) = O ρ(v) + t，det O < 0。

两条独立实现：
  (A) 自同构回溯枚举 + 自写 Kabsch（看最优正交阵的 det）；
  (B) RDKit 子结构匹配枚举原子对应 + `rdMolAlign.AlignMol`（只做真转动）叠镜像构象。
手性只取决于分子、不取决于构象 ⟹ 判定取「任一构象非手性即非手性」
（可挠的非手性分子，如乙醇、3-甲基戊烷，只有某个构象带镜面）。
"""
from __future__ import annotations

import json
import sys

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign

RDLogger.DisableLog("rdApp.*")

TOL = 0.20          # Å；坐标噪声上限
N_CONF = 40         # 每分子构象数（不剪枝：非手性分子往往只有个别构象带镜面）
MAX_AUT = 200_000


def has_cumulated(mol: Chem.Mol) -> bool:
    """是否含累积双键 C=C=C（丙二烯型）。"""
    pat = Chem.MolFromSmarts("[CX3]=[CX2]=[CX3]")
    return mol.HasSubstructMatch(pat)


def automorphisms(mol: Chem.Mol):
    """按元素标签+度枚举全部自同构（回溯）。"""
    n = mol.GetNumAtoms()
    labels = [a.GetAtomicNum() * 10 + int(a.GetFormalCharge()) for a in mol.GetAtoms()]
    adj = np.zeros((n, n), dtype=bool)
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        adj[i, j] = adj[j, i] = True
    deg = adj.sum(axis=1)
    cls: dict[tuple[int, int], list[int]] = {}
    for v in range(n):
        cls.setdefault((labels[v], int(deg[v])), []).append(v)
    order = sorted(range(n), key=lambda v: len(cls[(labels[v], int(deg[v]))]))
    pi = [-1] * n
    used = [False] * n
    out = []

    def rec(k):
        if len(out) >= MAX_AUT:
            return
        if k == n:
            out.append(tuple(pi))
            return
        v = order[k]
        for w in cls[(labels[v], int(deg[v]))]:
            if used[w]:
                continue
            ok = True
            for u in range(n):
                if pi[u] >= 0 and adj[v, u] != adj[w, pi[u]]:
                    ok = False
                    break
            if not ok:
                continue
            pi[v] = w
            used[w] = True
            rec(k + 1)
            used[w] = False
            pi[v] = -1

    rec(0)
    return out


def improper_rms(Xq, X):
    """把 X 拟合到 Xq，**强制**正交阵 det=−1：返回最小残差。

    平面构象（如苯、丁二烯）的最优拟合 det=+1、残差也为 0，单看「最优解 det」
    会漏判；直接求 det=−1 约束下的最小残差才是「是否存在非常规等距」。
    """
    Xqc = Xq - Xq.mean(0)
    Xc = X - X.mean(0)
    U, S, Vt = np.linalg.svd(Xc.T @ Xqc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, -d]) @ U.T      # 强制 det = −1
    return float(np.sqrt(((Xc @ R.T - Xqc) ** 2).sum(1).mean()))


def conf_improper_rms(X, auts):
    """对给定构象，返回 (最小非常规残差, 达到它的自同构数<tol)。"""
    best, n_hit = np.inf, 0
    for pi in auts:
        rms = improper_rms(X[np.asarray(pi)], X)
        if rms < best:
            best = rms
        if rms < TOL:
            n_hit += 1
    return best, n_hit


def mirror_rms_rdkit(mol, conf_id, auts):
    """(B) 镜像构象 + 只真转动叠合，枚举全部原子对应，取最小。"""
    n = mol.GetNumAtoms()
    conf = mol.GetConformer(conf_id)
    mm = Chem.Mol(mol)
    c2 = Chem.Conformer(conf)
    for i in range(n):
        p = conf.GetAtomPosition(i)
        c2.SetAtomPosition(i, (-p.x, p.y, p.z))
    mm.RemoveAllConformers()
    mm.AddConformer(c2, assignId=True)
    best = np.inf
    for pi in auts:
        amap = [(i, int(pi[i])) for i in range(n)]
        try:
            r = float(rdMolAlign.AlignMol(mm, mol, prbCid=0, refCid=conf_id, atomMap=amap))
        except Exception:
            continue
        if r < best:
            best = r
    return best


CASES = [
    ("甲烷", "C", "非手性"),
    ("甲醇", "CO", "非手性"),
    ("乙醇", "CCO", "非手性"),
    ("2-丙醇", "CC(C)O", "非手性"),
    ("甘氨酸", "NCC(=O)O", "非手性"),
    ("丙氨酸（L）", "C[C@@H](N)C(=O)O", "手性"),
    ("乳酸（L）", "C[C@@H](O)C(=O)O", "手性"),
    ("2-丁醇（R）", "CC[C@@H](C)O", "手性"),
    ("CHFClBr", "[C@@H](F)(Cl)Br", "手性"),
    ("CH2ClBr", "C(Cl)Br", "非手性"),
    ("苯", "c1ccccc1", "非手性"),
    ("甲苯", "Cc1ccccc1", "非手性"),
    ("内消旋酒石酸（R,S）", "O=C(O)[C@@H](O)[C@@H](O)C(=O)O", "非手性"),
    ("(S,S)-酒石酸", "O=C(O)[C@@H](O)[C@H](O)C(=O)O", "手性"),
    ("顺-1,2-二甲基环丙烷（R,S，内消旋）", "C[C@@H]1C[C@@H]1C", "非手性"),
    ("反-1,2-二甲基环丙烷（R,R）", "C[C@@H]1C[C@H]1C", "手性"),
    ("柠檬烯（R）", "CC(=C)[C@H]1CCC(=CC1)C", "手性"),
    ("α-蒎烯", "CC1=CC[C@@H](C)C1(C)C", "手性"),
    ("环己烷", "C1CCCCC1", "非手性"),
    ("樟脑", "CC1(C)C2CC[C@@]1(C)C(=O)C2", "手性"),
    ("葡萄糖（开链 D）", "OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C=O", "手性"),
    ("丝氨酸（L）", "N[C@@H](CO)C(=O)O", "手性"),
    ("半胱氨酸（L）", "N[C@@H](CS)C(=O)O", "手性"),
    ("异丙基氯", "CC(C)Cl", "非手性"),
    ("3-甲基戊烷", "CCC(C)CC", "非手性"),
    ("2-氯丁烷（R）", "CC[C@@H](C)Cl", "手性"),
    ("2,3-戊二烯（轴手性）", "CC=C=CC", "手性"),
    ("取代丙二烯 ClC=C=C(F)Br", "ClC=C=C(F)Br", "手性"),
    ("丁二烯（s-反式）", "C=CC=C", "非手性"),
    ("1,3-二氯环戊烷（顺式）", "Cl[C@H]1C[C@@H](Cl)CC1", "非手性"),
]


def main() -> int:
    rows = []
    for name, smi, expect in CASES:
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        n = mol.GetNumAtoms()
        auts = automorphisms(mol)
        ps = AllChem.ETKDGv3()
        ps.randomSeed = 0xF00D
        ps.pruneRmsThresh = -1.0            # 不剪枝
        cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=N_CONF, params=ps))
        if not cids:
            print(f"[skip] {name}")
            continue
        # 力场缺陷：MMFF/UFF 都把 C=C=C 的轴向二面角拉平到 180°（实测 ETKDG 嵌入时是 93.5°）。
        # 含累积双键的分子跳过力场优化，直接用 ETKDG 几何。
        cumulated = has_cumulated(Chem.RemoveHs(mol))
        if not cumulated:
            for cid in cids:
                try:
                    AllChem.MMFFOptimizeMolecule(mol, confId=cid, maxIters=400)
                except Exception:
                    pass
        best_a, hit_a, best_b = np.inf, 0, np.inf
        for cid in cids:
            X = np.asarray(mol.GetConformer(cid).GetPositions(), dtype=np.float64)
            ba, ha = conf_improper_rms(X, auts)
            bb = mirror_rms_rdkit(mol, cid, auts)
            best_a = min(best_a, ba)
            best_b = min(best_b, bb)
            hit_a += ha
        chiral_a = not (best_a < TOL)
        chiral_b = not (best_b < TOL)
        got = "手性" if chiral_a else "非手性"
        base = Chem.MolFromSmiles(smi)
        Chem.AssignStereochemistry(base, cleanIt=True, force=True)
        centers = Chem.FindMolChiralCenters(base, includeUnassigned=True, useLegacyImplementation=False)
        rows.append({
            "name": name, "smiles": smi, "expect": expect, "n_conf": len(cids),
            "n_atoms": n, "aut": len(auts), "geometry": "etkdg" if cumulated else "mmff",
            "improper_hits": int(hit_a),
            "rms_improper_best": float(best_a), "rms_mirror_best": float(best_b),
            "chiral_graph": chiral_a, "chiral_mirror": chiral_b,
            "rdkit_centers": len(centers),
            "agree_AB": chiral_a == chiral_b, "agree_expect": got == expect,
        })
        flag = "" if got == expect else "  ← 与期望不符"
        print(f"{name:<28} 构象{len(cids):>3} 自同构{len(auts):>5} 非真命中{hit_a:>3} "
              f"rms_a {best_a:7.3f} rms_b {best_b:7.3f} 判定 {got:<4} 期望 {expect:<4}{flag}")

    n = len(rows)
    ab = sum(r["agree_AB"] for r in rows)
    ex = sum(r["agree_expect"] for r in rows)
    print(f"\n共 {n} 个分子：两法一致 {ab}/{n}；与化学期望一致 {ex}/{n}")
    print("分歧：", [r["name"] for r in rows if not r["agree_AB"]])
    print("与期望不符：", [r["name"] for r in rows if not r["agree_expect"]])
    with open("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb092_chirality_topology.json", "w") as f:
        json.dump({"tol_angstrom": TOL, "n_conf": N_CONF, "n": n,
                   "agree_ab": ab, "agree_expect": ex, "cases": rows},
                  f, indent=2, ensure_ascii=False)
    print("落盘 ZeroGPU/results/zb092_chirality_topology.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
