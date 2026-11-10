#!/usr/bin/env python3
"""Zb084 §6.4 的机器核验：分子层面的 β₁（环数）与因果可导出性的记录门。

输出 ZeroGPU/results/zb084_molecule_beta1.json：
  mols   : 一批分子（含路线中间体）的 β₁
  p0     : 由 β₁=0（无环、无记录）导出的前生命集
  routes : 三条合成路线的 β₁ 序列与单调性
"""
import json, os
from rdkit import Chem

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results", "zb084_molecule_beta1.json")

MOLS = {
    "葡萄糖(吡喃)": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    "山梨醇": "OC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO",
    "山梨糖(呋喃)": "OC[C@H]1O[C@](O)(CO)[C@@H](O)[C@@H]1O",
    "2-酮基古洛糖酸": "O=C(O)C(=O)[C@H](O)[C@@H](O)[C@H](O)CO",
    "抗坏血酸": "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O",
    "FPP": "CC(C)=CCC/C(C)=C/CC/C(C)=C/COP(=O)(O)OP(=O)(O)O",
    "紫穗槐二烯": "C=C1CC[C@@H]2[C@@H]1CC=C(C)C2",
    "青蒿素": "CC1CCC2[C@@H](C)C(=O)[C@]3(C)OO[C@@]12[C@@H](OC(=O)C3)C",
    "萘": "c1ccc2ccccc2c1",
    "邻苯二甲酸酐": "O=C1OC(=O)c2ccccc21",
    "沙利度胺": "O=C1CCC(N2C(=O)c3ccccc3C2=O)C(=O)N1",
}

P0 = {  # β₁=0（无环、无记录）的候选：名称 -> SMILES
    "H2": "[H][H]", "O2": "O=O", "N2": "N#N",
    "H2O": "O", "H2S": "S", "NH3": "N", "PH3": "P",
    "CH4": "C", "CO": "[C-]#[O+]", "CO2": "O=C=O", "HCN": "C#N",
}

ROUTES = {
    "维C植物": ["葡萄糖(吡喃)", "山梨糖(呋喃)", "抗坏血酸"],
    "维C_Reichstein": ["葡萄糖(吡喃)", "山梨醇", "山梨糖(呋喃)", "2-酮基古洛糖酸", "抗坏血酸"],
    "青蒿素生物合成": ["FPP", "紫穗槐二烯", "青蒿素"],
    "沙利度胺工业": ["萘", "邻苯二甲酸酐", "沙利度胺"],
}

def beta1(smi):
    m = Chem.MolFromSmiles(smi)
    return None if m is None else m.GetRingInfo().NumRings()

def main():
    mols = {k: beta1(v) for k, v in MOLS.items()}
    p0 = {k: beta1(v) for k, v in P0.items()}
    routes = {}
    for name, seq in ROUTES.items():
        vals = [mols[s] for s in seq]
        mono = all(b >= a for a, b in zip(vals, vals[1:]))
        routes[name] = {"seq": seq, "beta1": vals, "monotone": mono}
    res = {
        "mols": mols,
        "p0": p0,
        "p0_all_acyclic": all(v == 0 for v in p0.values()),
        "routes": routes,
        "note": "β1 = RDKit 环数（SSSR）；记录门：β1=0 ⟹ 记录内容为零（Zb084 命题 4.3）",
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(res, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(json.dumps(res, ensure_ascii=False, indent=1))

if __name__ == "__main__":
    main()
