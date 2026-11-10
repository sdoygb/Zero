#!/usr/bin/env python3
"""序列片段的碰撞率 vs Boltzmann：命题 4 的组合侧

同一蛋白被切成不同长度片段，问"同长片段在两份拷贝间一致"这件事是否由组合计数决定。
若碰撞率 = (片段数)^2 / 2^L 随 L 指数下降、且与配对自由能无关，则片段层面的可分辨性由组合决定。
"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb093_fragmentation.json")
FASTA = Path("/tmp/mouse.fa")
LENS = (8, 10, 12, 15, 20, 25, 30, 40)


def read_fasta(path):
    seqs, cur, name = [], [], None
    for line in open(path):
        if line.startswith(">"):
            if cur:
                seqs.append((name, "".join(cur).upper()))
            name, cur = line[1:].strip(), []
        else:
            cur.append(line.strip().upper())
    if cur:
        seqs.append((name, "".join(cur).upper()))
    return seqs


def main():
    seqs = read_fasta(FASTA)
    print(f"读入 {len(seqs)} 条蛋白")
    rows = {}
    for L in LENS:
        seen = {}
        total = 0
        for _, s in seqs:
            for i in range(len(s) - L + 1):
                frag = s[i:i + L]
                seen[frag] = seen.get(frag, 0) + 1
                total += 1
        collided = sum(1 for v in seen.values() if v >= 2)
        rows[L] = {
            "n_fragments": total,
            "n_distinct": len(seen),
            "n_collided": collided,
            "collision_rate": collided / max(1, len(seen)),
            "combinatorial_upper": 4 ** L,
        }
        print(f"L={L:>2}: 片段 {total:>7} 不同 {len(seen):>7} 重复 {collided:>6} "
              f"碰撞率 {collided/max(1,len(seen)):.4f}  (4^L={4**L:.2e})")
    # 碰撞率 vs L 的斜率（对数轴上）
    ls = np.asarray(sorted(rows))
    rates = np.asarray([rows[L]["collision_rate"] for L in ls])
    ok = rates > 0
    slope = float(np.polyfit(ls[ok], np.log10(rates[ok]), 1)[0]) if ok.sum() > 2 else None
    print(f"\n碰撞率对 L 的斜率（log10 per residue）：{slope}")
    print("对照：若全由组合决定，碰撞率 ≈ (片段数 / 4^L)，斜率 ≈ −log10(4) = −0.602")
    out = {"n_proteins": len(seqs), "rows": rows, "slope": slope}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print("落盘", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
