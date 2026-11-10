#!/usr/bin/env python3
"""结构记录不是单值函数：两类不同构象映射到同一接触账（非单射）

Zb085 §5 的框架：记录写成 f:X→Y、重建是截面 g:Y→X（g∘f）。
要让"跨代不加分辨"在结构层成立，需先给出 f 与"纤维"（同一记录对应的全部结构）。
这里用两条独立判据量它：
  ① 同一蛋白的两种公认构象放在**同一张表**里（先把坐标交给同一套流程做整合）：
     腺苷酸激酶 open/closed、泛素晶体/NMR；
  ② 用每个结构的**内部距离矩阵**（模掉整体平移/旋转）做指纹，在同一张表里找最近邻，
     看"两个不同构象"能否被同一个记录代理（极简代理：接触数 + 稳定接触占比）区分开。
判据：若"不同构象"的代理值落在同一档、而它们互为一对的概率不高于随机，则记录对这层构象不分辨。
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/oygb/Downloads/Zero/ZeroGPU/src")
from zb093_structure_inheritance import parse_ca, contacts  # noqa: E402

PDB = Path("/Users/oygb/Downloads/Zero/_staging/推导/pdb")
OUT = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb093_structure_noninjective.json")


def internal_dist(model, keys):
    """内部 CA 距离矩阵（对上三角展平），模掉平移旋转。"""
    X = np.asarray([model[k] for k in keys], dtype=np.float64)
    D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=-1)
    iu = np.triu_indices(len(keys), k=1)
    return D[iu]


def load(model_path, chain):
    ms = parse_ca(PDB / model_path)
    out = []
    for m in ms:
        keys = sorted(k for k in m if k[0] == chain)
        if len(keys) < 10:
            continue
        out.append((internal_dist(m, keys), contacts({k: m[k] for k in keys}), len(keys)))
    return out


def main():
    groups = {
        "AK_A": ("4AKE.pdb", "1AKE.pdb", "A"),
        "UBQ_xtal_vs_nmr": ("1D3Z.pdb", "2NR2.pdb", "A"),
    }
    res = {}
    for name, (pa, pb, chain) in groups.items():
        a = load(pa, chain)
        b = load(pb, chain)
        print(f"\n=== {name}: {pa} {len(a)} 模型 vs {pb} {len(b)} 模型（链 {chain}） ===")
        # 只用共同残基数
        n = min(min(x[2] for x in a), min(x[2] for x in b))
        # 内部距离的互最近距：对每个 a 找最近的 b
        best = []
        for ia, ca, na in a:
            d0 = min(np.abs(ia[: len(ia)] - ib[: len(ib)]).mean() for ib, cb, nb in b) if False else None
            # 要求相同长度后比较：都取前 n 个残基的距离矩阵不方便，改为比较整矩阵的平均绝对差
            best.append(min(np.abs(ia - ib).mean() for ib, _, _ in b))
        best = np.asarray(best)
        # 同组内（a 内部）的最近距作为对照
        within = []
        for i in range(len(a)):
            within.append(min(np.abs(a[i][0] - a[j][0]).mean() for j in range(len(a)) if j != i) if len(a) > 1 else 0.0)
        within = np.asarray(within)
        # 接触代理：接触数
        nc_a = np.asarray([x[1].__len__() for x in a])
        nc_b = np.asarray([x[1].__len__() for x in b])
        overlap = len(set().union(*[x[1] for x in a]) & set().union(*[x[1] for x in b])) / max(1, len(set().union(*[x[1] for x in a])))
        r = {
            "n_a": len(a), "n_b": len(b),
            "dist_a_to_b_nearest_mean": float(best.mean()),
            "dist_a_to_b_nearest_min": float(best.min()),
            "dist_within_a_nearest_mean": float(within.mean()) if len(within) else None,
            "contact_count_a": [int(nc_a.min()), int(nc_a.max())],
            "contact_count_b": [int(nc_b.min()), int(nc_b.max())],
            "contact_union_overlap_frac": float(overlap),
        }
        res[name] = r
        print(f"  内部距离(a→b 最近) 均值 {best.mean():.3f} Å；a 内部最近 均值 {np.mean(within):.3f} Å")
        print(f"  接触数 a {nc_a.min()}–{nc_a.max()}，b {nc_b.min()}–{nc_b.max()}；接触集合交集占 a 的并集 {overlap:.1%}")
        print(f"  ⟹ 判读：{'两类构象在接触账上几乎不可分（非单射）' if overlap > 0.9 else '接触账能部分区分两类'}；"
              f"内部距离差 {best.mean():.3f} vs 组内 {np.mean(within):.3f}："
              f"{'构象确实不同，但接触账看不见' if best.mean() > 3 * max(np.mean(within), 1e-6) else '构象差与组内噪声同量级'}")

    print("\n结论口径：若『构象在内部距离上明显不同（跨组差 >> 组内差）』而『接触账交集 > 90%』，"
          "则 f:结构→接触账 不是单射——同一记录对应多个构象，正是 Zb085 跨代算子里 f 的纤维；"
          "此时 f 后的重建 g 必然丢掉这层构象信息。")
    with open(OUT, "w") as f:
        json.dump(res, f, indent=2, ensure_ascii=False)
    print("落盘", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())