#!/usr/bin/env python3
"""③ 的 (a) 层落地：用真核基因组 ORF 做"闭合碰撞"统计

判据式 7.1 的 (a)：每一步是**局部闭合**（成键/成环）。落到分子式层面就是"闭合计数 β₁ 只按整数步增长"。
现实里"分子式"允许同分异构/同质量替代，于是同一个 $P_0$（元素计数）可对应多个不同的 β₁ 结构——
**闭合碰撞** c = 同一元素计数下的不同 β₁ 结构数。c 越大，判据 (a) 的"一步一闭合"越粗，可分辨性越低。

可测的代理：真核基因组里的 ORF（蛋白编码序列）—— 同一长度上的不同序列就是一个"同计数、不同结构"的例子。
测两件事：
  · 每个长度的 ORF 数（= 同计数、不同序列的计数）；
  · 同一长度内的序列两两重数（同计数下的"不同结构"）。
"""
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

OUT = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb093_closure_collision.json")
GENOME = Path("/tmp/yeast6.fa")     # S. cerevisiae chr VI，270 kb


def read_seq(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def orfs(seq, min_len=90, max_len=3000):
    """六个读框里找 ATG…TAA/TAG/TGA 的 ORF（不考虑内含子，作为代理）。"""
    stops = {"TAA", "TAG", "TGA"}
    out = []
    for frame in range(3):
        cur = []
        for i in range(frame, len(seq) - 2, 3):
            c = seq[i:i + 3]
            if c == "ATG" and not cur:
                cur = [i // 3]
            if cur and c in stops:
                ln = i + 3 - cur[0] * 3
                if min_len <= ln <= max_len:
                    out.append((cur[0] * 3, ln, seq[cur[0] * 3:i + 3]))
                cur = []
    return out


def main():
    seq = read_seq(GENOME)
    found = orfs(seq)
    print(f"基因组 {len(seq)} bp；ORF（≥90 nt，含终止子）{len(found)} 个")
    by_len = defaultdict(list)
    for pos, ln, s in found:
        by_len[ln].append(s)
    lens = np.asarray(sorted(by_len))
    counts = np.asarray([len(by_len[L]) for L in lens])
    res = {
        "genome_bp": len(seq),
        "n_orfs": len(found),
        "n_distinct_lengths": int(len(lens)),
        "length_range": [int(lens.min()), int(lens.max())] if len(lens) else None,
        "orfs_per_length_median": float(np.median(counts)) if len(counts) else None,
        "orfs_per_length_max": int(counts.max()) if len(counts) else None,
        "collision_frac_lengths_with_gt1": float((counts > 1).mean()) if len(counts) else None,
    }
    print(f"  不同长度 {res['n_distinct_lengths']}；长度区间 {res['length_range']}")
    print(f"  同长度 ORF 数：中位 {res['orfs_per_length_median']:.1f}，最大 {res['orfs_per_length_max']}")
    print(f"  出现「同计数多个结构」的长度占比 {res['collision_frac_lengths_with_gt1']:.1%}")

    # 闭合增量代理：每步成肽键 = 一次闭合；β₁ 在链状分子上恒为 0，环化才是闭合
    # 故这里量的是"同长度下不同序列结构"的数目——(a) 的分辨下限
    if len(counts):
        print("\n判读：\n"
              f"  同一元素计数（同长度）下平均 {np.median(counts):.1f} 个不同结构；\n"
              f"  若判据 (a) 只按'闭合计数'分辨，则这些结构必须共用同一 β₁，\n"
              f"  (a) 的分辨下限 = 同计数结构数 ≈ {np.median(counts):.0f}。\n"
              "  即：(a) 单独不能把'同分子式的不同结构'分开，必须叠加 (d) 的锥界/机制才够。")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=2, ensure_ascii=False)
    print("\n落盘", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())