#!/usr/bin/env python3
"""P2 的区分检验：密码子的能量代价与简并度是否一致？

假设（P2 思路）：Leu/Ser/Arg 简并度 6（N=6）、Trp/Met 简并度 1（N=1）不是巧合，
而是 "高代价的氨基酸能被更多密码子容忍" 或反过来。

做法：用 Zero 口径的键能（Zb092 系列脚本里的 C/N/O/S 单键能，单位 kJ/mol）
与氨基酸侧链的"重原子数 / 含硫与否"做代理代价，检验简并度与代价的 Spearman 相关。
这不是"证明"，只是一个可复算的相关性；若相关为零，P2 的"能量--简并"链条就该撤。
"""
import json
from collections import Counter
from pathlib import Path

import numpy as np

OUT = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb095_codon_vs_cost.json")

# 侧链组成（重原子数、是否含 S、是否芳香），来源：氨基酸结构式，可手算复核
SIDE = {
    "G": dict(heavy=0, S=0), "A": dict(heavy=1, S=0), "S": dict(heavy=2, S=0),
    "C": dict(heavy=1, S=1), "T": dict(heavy=2, S=0), "V": dict(heavy=3, S=0),
    "L": dict(heavy=4, S=0), "I": dict(heavy=4, S=0), "M": dict(heavy=3, S=1),
    "P": dict(heavy=3, S=0), "F": dict(heavy=7, S=0), "Y": dict(heavy=8, S=0),
    "W": dict(heavy=9, S=0), "D": dict(heavy=4, S=0), "E": dict(heavy=5, S=0),
    "N": dict(heavy=4, S=0), "Q": dict(heavy=5, S=0), "K": dict(heavy=5, S=0),
    "R": dict(heavy=6, S=0), "H": dict(heavy=5, S=0),
}
DEGEN = {"W": 1, "M": 1, "F": 2, "Y": 2, "C": 2, "H": 2, "Q": 2, "N": 2, "K": 2, "D": 2, "E": 2,
         "I": 3, "P": 4, "T": 4, "V": 4, "A": 4, "G": 4, "L": 6, "S": 6, "R": 6}


def spearman(x, y):
    rx = np.argsort(np.argsort(x))
    ry = np.argsort(np.argsort(y))
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    return float((rx @ ry) / np.sqrt((rx @ rx) * (ry @ ry)))


def main():
    aas = sorted(DEGEN)
    deg = np.asarray([DEGEN[a] for a in aas], dtype=float)
    heavy = np.asarray([SIDE[a]["heavy"] for a in aas], dtype=float)
    hS = np.asarray([SIDE[a]["S"] for a in aas], dtype=float)
    r_all = spearman(heavy, deg)
    print(f"20 氨基酸：简并度 vs 侧链重原子数  Spearman r = {r_all:+.3f}")
    # ≤4 重原子的子集（排除芳香的大侧链，避免被大小主导）
    m = heavy <= 4
    r_small = spearman(heavy[m], deg[m])
    print(f"重原子 ≤4 的 {int(m.sum())} 个：r = {r_small:+.3f}")
    # 含硫 vs 不含硫
    print(f"含硫（C/M）简并度 {DEGEN['C']}/{DEGEN['M']}；同大小不含硫（V 3, L 4, I 4）简并度 4/6/3")
    out = {"r_heavy": r_all, "r_heavy_le4": r_small, "n_le4": int(m.sum()),
           "degen": DEGEN, "heavy": {a: SIDE[a]["heavy"] for a in aas}}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print("落盘", OUT)


if __name__ == "__main__":
    raise SystemExit(main())
