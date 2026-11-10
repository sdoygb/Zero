#!/usr/bin/env python3
"""Zero 分子指纹器：SMILES → (β1, ω, Δ, k, χ, 计数)。

用法:
  python3 zerofp.py "c1ccccc1" "CCO" ...
  python3 zerofp.py --json "c1ccccc1"

输出字段:
  formula   分子式（含 H）
  V,E,c     重原子图（H 抑制）的顶点/边/连通分量
  beta1     β1 = E - V + c（独立环数）
  delta     共轭缺陷 = ||spec(Γ) - spec(L(Γ))||_2
  omega     L1′ 记录向量（全图，含 H；ω_v = deg_v / Σdeg）
  rings     最小环长度与绕数 k = Σb/3（无环则为空）
  chiral    手性中心数（RDKit FindMolChiralCenters）
  chi_flag  是否允许 χ≠0（手性中心>0 或 无镜面对称的近似判据）
"""
import sys, json
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def heavy(n):
    return [a.GetIdx() for a in n.GetAtoms() if a.GetAtomicNum() > 1]

def adjacency_heavy(mol):
    idx = heavy(mol); m = {v: i for i, v in enumerate(idx)}
    E = []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in m and j in m: E.append((m[i], m[j]))
    return len(idx), E

def spec(n, E):
    A = np.zeros((n, n))
    for i, j in E: A[i, j] = A[j, i] = 1
    return np.sort(np.linalg.eigvalsh(A)) if n else np.zeros(0)

def line_edges(E):
    out = []
    for a in range(len(E)):
        for b in range(a + 1, len(E)):
            if set(E[a]) & set(E[b]): out.append((a, b))
    return out

def delta(mol):
    n, E = adjacency_heavy(mol)
    s1, s2 = spec(n, E), spec(len(E), line_edges(E))
    L = max(len(s1), len(s2)); a = np.zeros(L); b = np.zeros(L)
    a[:len(s1)] = s1; b[:len(s2)] = s2
    return float(np.linalg.norm(np.sort(a) - np.sort(b)))

def omega(mol):
    m = Chem.AddHs(mol)
    deg = np.array([a.GetDegree() for a in m.GetAtoms()], float)
    s = deg.sum()
    return [round(float(x / s), 6) for x in deg]

def rings_info(mol):
    """SSSR 逐环：长度 n 与绕数 k = Σb/3。"""
    ri = mol.GetRingInfo()
    out = []
    for r in ri.AtomRings():
        s = 0.0
        for i in range(len(r)):
            b = mol.GetBondBetweenAtoms(r[i], r[(i + 1) % len(r)])
            if b is None: continue
            s += 1.5 if b.GetIsAromatic() else b.GetBondTypeAsDouble()
        out.append({"n": len(r), "sum_b": round(s, 3), "k": round(s / 3.0, 4)})
    out.sort(key=lambda x: x["n"])
    return out

def omega_entropy(om):
    """ω 的账本熵 H = -Σ ω ln ω（nat）与归一化 H/ln N。"""
    import math
    H = -sum(x * math.log(x) for x in om if x > 0)
    return round(H, 6), round(H / math.log(len(om)), 6) if len(om) > 1 else 0.0

def fingerprint(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None: return {"smiles": smi, "error": "parse failed"}
    n, E = adjacency_heavy(mol)
    frags = Chem.GetMolFrags(mol)
    c = len(frags)
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    return {
        "smiles": smi,
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "V": n, "E": len(E), "c": c,
        "beta1": len(E) - n + c,
        "delta": round(delta(mol), 6),
        "omega": omega(mol),
        "omega_H": omega_entropy(omega(mol))[0],
        "omega_H_norm": omega_entropy(omega(mol))[1],
        "rings": rings_info(mol),
        "chiral_centers": len(centers),
        "chi_allowed": len(centers) > 0,
    }

def main():
    args = sys.argv[1:]
    js = "--json" in args
    args = [a for a in args if a != "--json"]
    if not args: args = ["c1ccccc1"]
    out = [fingerprint(s) for s in args]
    if js:
        print(json.dumps(out, ensure_ascii=False, indent=1))
    else:
        for d in out:
            if "error" in d: print(d["smiles"], "ERROR", d["error"]); continue
            r = ",".join(f"n{x['n']}/k{x['k']}" for x in d["rings"]) or "-"
            print(f"{d['smiles']:22s} {d['formula']:12s} β1={d['beta1']}  Δ={d['delta']:.4f}  "
                  f"rings=[{r}]  H(ω)={d['omega_H']:.4f} (H/lnN={d['omega_H_norm']:.4f})  χ={d['chiral_centers']}")

if __name__ == "__main__":
    main()
