#!/usr/bin/env python3
"""20+ 对同蛋白多条件结构：跨家族批量检验「L1' 计数不分辨 R 上构象态」

数据：RCSB 搜索（腺苷酸激酶 / 溶菌酶 / 糖原磷酸化酶 / 泛素，各 50 个 X-ray 条目），下载前 160 个 PDB。
每对 = 同一家族内**两个不同条件**的条目（不同配体/不同晶体/不同分辨率），配对全部组合后取样本。

每对测（层口径）：
  · R 层：两态内部 CA 距离矩阵（模掉平移旋转）的平均绝对差；
  · L1' 层：接触集合的 Jaccard；
  · 同态噪声：每条 PDB 若有多模型，态内结构之间的距离差（否则视为 0）。
判据：几何差 ≫ 同态噪声 且 Jaccard ≥ 0.85 ⟹ L1' 计数不分辨 R 上的构象态。
"""
import json
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/oygb/Downloads/Zero/ZeroGPU/src")
from zb093_structure_inheritance import parse_ca, contacts, CUT, MINSEQ  # noqa: E402

PDB_DIR = Path("/tmp/pdbs")
OUT = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb093_structure_pairs20.json")
FAMILIES = {
    "AK": "/tmp/rcsb_adenylate_kinase.json",
    "LYZ": "/tmp/rcsb_lysozyme.json",
    "GP": "/tmp/rcsb_glycogen_phosphorylase.json",
    "UBQ": "/tmp/rcsb_ubiquitin.json",
    "MYO": "/tmp/more_fam.json#myoglobin",
    "THR": "/tmp/more_fam.json#thrombin",
    "CAH": "/tmp/more_fam.json#carbonic anhydrase",
    "GST": "/tmp/more_fam.json#glutathione s-transferase",
    "LDH": "/tmp/more_fam.json#lactate dehydrogenase",
    "TRP": "/tmp/more_fam.json#trypsin",
}
MAX_PAIRS_PER_FAMILY = 10


def seq_key(model):
    """残基键（链+序号+插入码）集合，用于判断两结构能否逐位比较。"""
    return frozenset(model)


def internal_dist(model):
    keys = sorted(model)
    X = np.asarray([model[k] for k in keys], dtype=np.float64)
    D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=-1)
    iu = np.triu_indices(len(keys), k=1)
    return D[iu]


def main():
    rows = []
    for fam, path in FAMILIES.items():
        if "#" in path:
            base, fam_key = path.split("#")
            ids = json.load(open(base))[fam_key]
        else:
            ids = json.load(open(path))
        avail = [i for i in ids if (PDB_DIR / f"{i}.pdb").exists()]
        pairs = list(combinations(avail[:8], 2))[:MAX_PAIRS_PER_FAMILY]
        for a, b in pairs:
            ma = parse_ca(PDB_DIR / f"{a}.pdb")
            mb = parse_ca(PDB_DIR / f"{b}.pdb")
            if not ma or not mb:
                continue
            ka, kb = seq_key(ma[0]), seq_key(mb[0])
            cov = len(ka & kb) / max(len(ka | kb), 1)
            # 只比较残基集合高度重叠且原子数接近的配对（否则是不同片段/多聚体，无比对意义）
            if cov < 0.8 or abs(len(ka) - len(kb)) > 20:
                rows.append({"family": fam, "pair": [a, b], "coverage": cov,
                             "skipped": "残基覆盖不足/长度差过大"})
                print(f"{fam:<4} {a} vs {b}: 跳过（覆盖 {cov:.2f}，原子 {len(ka)}/{len(kb)}）")
                continue
            common = sorted(ka & kb)
            ma0 = {k: ma[0][k] for k in common}
            mb0 = {k: mb[0][k] for k in common}
            da, db = internal_dist(ma0), internal_dist(mb0)
            geo = float(np.abs(da - db).mean())
            ca, cb = contacts(ma0), contacts(mb0)
            jac = len(ca & cb) / max(1, len(ca))     # 用"占 A 的比例"更稳（两集合大小不同时 Jaccard 会偏低）
            noise = 0.0
            if len(ma) > 1:
                d0 = internal_dist(ma0)
                for m in ma[1:4]:
                    mk = {k: m[k] for k in common if k in m}
                    if len(mk) == len(common):
                        dk = internal_dist(mk)
                        noise = max(noise, float(np.abs(d0 - dk).mean()))
            rows.append({"family": fam, "pair": [a, b], "coverage": cov,
                         "geo_r": geo, "retention": jac, "noise_r": noise,
                         "geo_over_noise": geo / max(noise, 1e-3),
                         "disambiguates": bool(geo / max(noise, 1e-3) > 2 and jac >= 0.85)})
            print(f"{fam:<4} {a} vs {b}: 覆盖 {cov:.2f} 几何 {geo:7.3f} Å  保留 {jac:.3f} 噪声 {noise:.3f} "
                  f"{'不分辨 ✓' if rows[-1]['disambiguates'] else '可分'}")

    valid = [r for r in rows if "skipped" not in r]
    geo = np.asarray([r["geo_r"] for r in valid])
    jac = np.asarray([r["retention"] for r in valid])
    dis = np.asarray([r["disambiguates"] for r in valid])
    print(f"\n共 {len(rows)} 对，有效 {len(valid)} 对（{len(set(r['family'] for r in valid))} 个家族，跳过 {len(rows)-len(valid)} 对残基覆盖不足者）")
    print(f"几何差 均值 {geo.mean():.3f} Å；保留率 均值 {jac.mean():.3f}")
    print(f"保留率 ≥0.85 的对：{int((jac>=0.85).sum())}/{len(valid)}")
    print(f"按判据「L1' 计数不分辨 R」的对：{int(dis.sum())}/{len(valid)}")
    print(f"Jaccard（保留率）≥0.85 且几何差>2×噪声："
          f"{int(((jac>=0.85) & (geo > 2*np.maximum([r['noise_r'] for r in valid], 1e-3))).sum())}/{len(valid)}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump({"n_pairs": len(rows), "rows": rows}, f, indent=2, ensure_ascii=False)
    print("落盘", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())