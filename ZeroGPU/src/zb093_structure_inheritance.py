#!/usr/bin/env python3
"""结构可传承性检验：折叠记账量对**环境扰动**（pH/配体/晶体接触）的敏感性

问题（Zb085 §5 的"跨代"落到结构上）：一份结构记录要能传承，重建算子 g 对环境变化必须**不敏感**；
若敏感，则"重建"只是把扰动一起抄下来，不是传承。
把"环境扰动"换成"同一蛋白在不同条件下的两份实测结构"：
    open/closed（加/不加配体与 AlF3）   ：1AKE vs 4AKE
    NAC 态 vs 折叠态（同源片段的两种构象）：2K39 vs 1G6J
    pH 6.5/4.5（同蛋白、两张晶体）      ：1D3Z vs 2NR2 / 2MJB
    NMR 系综内（同一蛋白同一条件）       ：每个系综内部
记账量（Zb074 的 R 在同一舞台上的代理，全部只依赖 CA 坐标，接触=序列间隔≥2 且 CA≤8 Å）：
    n_contact  接触（闭合）总数
    const      全部模型/两张结构都保持的接触数（稳定骨架）
    modul      只在一部分结构里出现的接触数（可变部分）
    switch/T   相邻结构之间的接触开关率 = |R_i Δ R_{i+1}| / n_contact / 结构数
"""
import json
import math
import sys
from pathlib import Path

PDB = Path("/Users/oygb/Downloads/Zero/_staging/推导/pdb")
OUT = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb093_structure_inheritance.json")
CUT = 8.0
MINSEQ = 2


def parse_ca(path):
    """返回每个模型一个 dict{(chain, resSeq): (x,y,z)}，以及坐标空间差。"""
    models, cur = [], {}
    with open(path, errors="ignore") as f:
        for line in f:
            if line.startswith("MODEL"):
                cur = {}
            elif line.startswith("ATOM") and line[12:16].strip() == "CA":
                alt = line[16]
                if alt not in (" ", "A"):
                    continue
                key = (line[21], int(line[22:26]), line[26])
                cur[key] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            elif line.startswith("ENDMDL"):
                if cur:
                    models.append(cur)
                    cur = {}
    if cur:
        models.append(cur)
    return models


def contacts(model):
    """返回接触键集合 {(k1, k2)}，k=(chain, resSeq, icode)，要求序列间隔≥2。"""
    keys = sorted(model)
    pos = {k: i for i, k in enumerate(keys)}
    out = set()
    n = len(keys)
    for i in range(n):
        xi, yi, zi = model[keys[i]]
        for j in range(i + 1, n):
            if abs(pos[keys[j]] - pos[keys[i]]) < MINSEQ:
                continue
            xj, yj, zj = model[keys[j]]
            d2 = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
            if d2 <= CUT * CUT:
                out.add((keys[i], keys[j]))
    return out


def analyse(models):
    cons = [contacts(m) for m in models]
    cats = [frozenset(c) for c in cons]
    n = len(cats)
    if n == 0:
        return None
    inter = set(cats[0])
    union = set(cats[0])
    for c in cats[1:]:
        inter &= c
        union |= c
    switches = []
    for i in range(n - 1):
        switches.append(len(cats[i] ^ cats[i + 1]))
    denom = len(union) if union else 1
    return {
        "n_struct": n,
        "n_contact_mean": sum(len(c) for c in cats) / n,
        "const": len(inter),
        "modul": len(union) - len(inter),
        "const_frac": len(inter) / len(union) if union else 0.0,
        "switch_per_contact": (sum(switches) / (n - 1) / denom) if n > 1 else 0.0,
        "residues": len(models[0]),
    }


def main():
    sets = {
        # 名字: (文件列表, 类别, 说明)
        "AK_open": ([PDB / "4AKE.pdb"], "单结构", "腺苷酸激酶 无配体"),
        "AK_closed": ([PDB / "1AKE.pdb"], "单结构", "腺苷酸激酶 +AP5A"),
        "GB1": ([PDB / "2GB1.pdb"], "单结构", "蛋白 G B1 结构域"),
        "NAC_2K39": ([PDB / "2K39.pdb"], "系综（同一条件）", "折叠态系综"),
        "NAC_1G6J": ([PDB / "1G6J.pdb"], "系综（同一条件）", "NAC 系综"),
        "UBQ_1D3Z": ([PDB / "1D3Z.pdb"], "系综", "泛素 pH 6.5 晶体模型集"),
        "UBQ_2MJB": ([PDB / "2MJB.pdb"], "系综", "泛素 pH 4.5"),
        "UBQ_2NR2": ([PDB / "2NR2.pdb"], "系综", "泛素 pH 6.5 (NMR)"),
        "UBQ_env": ([PDB / "1D3Z.pdb", PDB / "2MJB.pdb", PDB / "2NR2.pdb"], "跨条件", "泛素三个条件合看"),
        "AK_env": ([PDB / "4AKE.pdb", PDB / "1AKE.pdb"], "跨条件", "腺苷酸激酶 open vs closed"),
        "NAC_env": ([PDB / "2K39.pdb", PDB / "1G6J.pdb"], "跨条件", "NAC 折叠态 vs 未折叠态"),
    }
    rows = {}
    for name, (files, cat, desc) in sets.items():
        models = []
        for f in files:
            models.extend(parse_ca(f))
        r = analyse(models)
        if r is None:
            print(f"[skip] {name}")
            continue
        r["category"] = cat
        r["desc"] = desc
        rows[name] = r
        print(f"{name:<10} {cat:<12} 结构{r['n_struct']:>4} 残基{r['residues']:>4} "
              f"接触均值{r['n_contact_mean']:7.0f} 稳定{r['const']:6d} 可变{r['modul']:6d} "
              f"稳定占比{r['const_frac']:6.1%} 开关/接触{r['switch_per_contact']:8.4f}")

    # 两类比较：同条件系综内 vs 跨条件
    within = [r for r in rows.values() if r["category"].startswith("系综")]
    across = [r for r in rows.values() if r["category"] == "跨条件"]
    if within and across:
        w = sum(r["const_frac"] for r in within) / len(within)
        a = sum(r["const_frac"] for r in across) / len(across)
        print(f"\n稳定占比：同条件系综均值 {w:.3f}；跨条件均值 {a:.3f}；差 {w - a:+.3f}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump({"cut_angstrom": CUT, "minseq": MINSEQ, "rows": rows}, f, indent=2, ensure_ascii=False)
    print("落盘", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())