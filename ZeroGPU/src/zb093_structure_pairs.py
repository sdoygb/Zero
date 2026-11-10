#!/usr/bin/env python3
"""A5 扩展：20+ 对同蛋白多条件结构的层口径核验

目标（用户指示）：把「L1' 计数不分辨 R 上构象态」从 2 个蛋白扩到 20+ 对。

每一对 = 同一蛋白在两种条件下（配体/构象/晶体 vs NMR/pH 等）的结构。对每对测：
  · R 层：内部 CA 距离矩阵的平均绝对差（模掉平移旋转）——两态的几何差；
  · L1' 层：接触集合的交集/并集（Jaccard）——两态的计数差；
  · 构象引擎（L2）的涨落：每态内部各结构之间的内部距离差的均值——同态噪声尺度。
判据：几何差 >> 同态噪声 且 Jaccard ≥ 0.85 ⟹ L1' 计数不分辨 R 上的构象态（该对成立）。
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/oygb/Downloads/Zero/ZeroGPU/src")
from zb093_structure_inheritance import parse_ca, contacts  # noqa: E402

PDB = Path("/Users/oygb/Downloads/Zero/_staging/推导/pdb")
OUT = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb093_structure_pairs.json")


def internal_dist(model, chain=None):
    keys = sorted(k for k in model if chain is None or k[0] == chain)
    X = np.asarray([model[k] for k in keys], dtype=np.float64)
    D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=-1)
    iu = np.triu_indices(len(keys), k=1)
    return D[iu], keys


def mean_abs_diff(a, b):
    n = min(len(a), len(b))
    return float(np.abs(a[:n] - b[:n]).mean())


def pair_metrics(fa, fb, chain=None):
    ma = parse_ca(PDB / fa)
    mb = parse_ca(PDB / fb)
    # 每态取第一个模型的代表，用于跨态比较
    da, ka = internal_dist(ma[0], chain)
    db, kb = internal_dist(mb[0], chain)
    geo = mean_abs_diff(da, db)
    ca = contacts({k: ma[0][k] for k in ka})
    cb = contacts({k: mb[0][k] for k in kb})
    jac = len(ca & cb) / max(1, len(ca | cb))
    # 同态噪声：态 A 内部各结构两两的内部距离差
    noise = []
    if len(ma) > 1:
        ref, kk = internal_dist(ma[0], chain)
        for m in ma[1:6]:
            d, _ = internal_dist(m, chain)
            noise.append(mean_abs_diff(ref, d))
    res = {"pair": [fa, fb], "geo_r": geo, "contact_jaccard": jac,
           "within_noise_r": float(np.mean(noise)) if noise else 0.0,
           "n_a": len(ma), "n_b": len(mb)}
    res["geo_over_noise"] = geo / max(res["within_noise_r"], 1e-3)
    res["disambiguates"] = bool(res["geo_over_noise"] > 2 and jac >= 0.85)
    res["differs_from_noise"] = bool(res["geo_over_noise"] > 2)
    return res


def build_structure_pairs():
    """先用本地 8 个结构做 6 对；再按需补 UniProt/PDB 网络对。"""
    structures = {
        "AK_open": ["4AKE.pdb"],
        "AK_closed": ["1AKE.pdb"],
        "GB1": ["2GB1.pdb"],
        "NAC_fold": ["2K39.pdb"],
        "NAC_unfold": ["1G6J.pdb"],
        "UBQ_65": ["1D3Z.pdb"],
        "UBQ_45": ["2MJB.pdb"],
        "UBQ_nmr": ["2NR2.pdb"],
    }
    pairs = [
        ("AK_open", "AK_closed"),
        ("NAC_fold", "NAC_unfold"),
        ("UBQ_65", "UBQ_45"),
        ("UBQ_65", "UBQ_nmr"),
        ("UBQ_45", "UBQ_nmr"),
    ]
    rows = []
    for a, b in pairs:
        r = pair_metrics(structures[a][0], structures[b][0])
        rows.append(r)
        print(f"{a:<12} vs {b:<12} 几何差 {r['geo_r']:6.3f} Å  Jaccard {r['contact_jaccard']:.3f} "
              f"同态噪声 {r['within_noise_r']:.3f}  {'不分辨 ✓' if r['disambiguates'] else '可分'}")
    return {"rows": rows, "note": "本地 8 个结构可供 5 对；补到 20+ 对需要更多同蛋白多条件 PDB"}


def main():
    res = build_structure_pairs()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=2, ensure_ascii=False)
    print(f"\n共 {len(res['rows'])} 对；落盘 {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())