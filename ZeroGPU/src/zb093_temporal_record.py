#!/usr/bin/env python3
"""时序接触账（fold-order record）能否把构象态分开？

背景：`_staging/推导/结构可传承性_接触账分辨.md` 命题 1 说无时序的接触集合对构象不是单射。
若把记录换成**带时序的接触账**——按"先形成的接触"排序写出接触序列——能否把 open/closed 分开？

时序怎么定（两套独立口径，都只用几何）：
  A. 接触量级排序：接触按"何时出现"排（序列间隔 = 形成顺序），形成顺序由 Zb074 的"局部闭合先出现"给，
     代理取 |i−j|（局部接触先闭合），同量级按 CA 距离排；
  B. 折叠核顺序：先排 Cα 根序列（i<j，按 j 递增），即"从 N 端向 C 端生长"的闭合顺序。

判据：两份结构（open vs closed）的"有序接触签名"的前缀一致度——记录只被压成前缀时还能否区分两态。
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/oygb/Downloads/Zero/ZeroGPU/src")
from zb093_structure_inheritance import parse_ca, contacts, CUT, MINSEQ  # noqa: E402

PDB = Path("/Users/oygb/Downloads/Zero/_staging/推导/pdb")
OUT = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb093_temporal_record.json")


def ordered_signature(model, order="local"):
    """返回有序接触列表。order='local'：按 |i−j| 升序（局部先闭合）；order='growth'：按 j 升序（N→C 生长）。"""
    keys = sorted(model)
    pos = {k: i for i, k in enumerate(keys)}
    items = []
    n = len(keys)
    for i in range(n):
        xi, yi, zi = model[keys[i]]
        for j in range(i + 1, n):
            if pos[keys[j]] - pos[keys[i]] < MINSEQ:
                continue
            xj, yj, zj = model[keys[j]]
            d2 = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
            if d2 <= CUT * CUT:
                d = np.sqrt(d2)
                if order == "local":
                    items.append((pos[keys[j]] - pos[keys[i]], d, keys[i], keys[j]))
                else:
                    items.append((pos[keys[j]], pos[keys[i]], keys[i], keys[j]))
    items.sort(key=lambda t: (t[0], t[1]))
    return [frozenset(t[2:]) for t in items]


def prefix_overlap(sig_a, sig_b):
    """按"第 k 个接触"对齐：前 k 个接触集合的交集大小 / k。"""
    k_max = min(len(sig_a), len(sig_b))
    set_a = set(sig_a)
    set_b = set(sig_b)
    # 集合层面与"每个位置上的接触"
    return [len(set_a & set_b) / max(1, len(set_a))]


def main():
    res = {}
    for label, order in (("local", "local"), ("growth", "growth")):
        o = ordered_signature(parse_ca(PDB / "4AKE.pdb")[0], order)
        c = ordered_signature(parse_ca(PDB / "1AKE.pdb")[0], order)
        # 前缀一致度：取前 k 个接触，算 Jaccard
        ks = [10, 25, 50, 100, 200, 500, 1000]
        rows = []
        for k in ks:
            if k > min(len(o), len(c)):
                continue
            a, b = set(o[:k]), set(c[:k])
            jac = len(a & b) / len(a | b)
            rows.append({"k": k, "jaccard": jac})
        res[order] = {
            "n_open": len(o), "n_closed": len(c),
            "unordered_intersection_over_union": len(set(o) & set(c)) / len(set(o) | set(c)),
            "prefix": rows,
        }
        print(f"\n[{order}] open {len(o)} 接触 / closed {len(c)}；无序并集交集 {len(set(o) & set(c)) / len(set(o) | set(c)):.3f}")
        for r in rows:
            print(f"   前 {r['k']:>4} 个接触的 Jaccard：{r['jaccard']:.3f}")
        # 判读：前缀 Jaccard 明显低于无序版 ⟹ 时序记录能把两态分开
        if rows:
            low = min(r["jaccard"] for r in rows)
            print(f"   判读：前缀最差 Jaccard {low:.3f} vs 无序 {res[order]['unordered_intersection_over_union']:.3f} ⟹ "
                  f"{'时序记录能部分分开两态' if low < res[order]['unordered_intersection_over_union'] - 0.05 else '时序记录也分不开'}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=2, ensure_ascii=False)
    print("\n落盘", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())