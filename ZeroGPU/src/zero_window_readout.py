#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""L=8 词里的 L=4 读出窗口：q_rot=5/9 的核验。

问题：Z062 用整个 W_L 算 q_rot(L)；L=8 给 19/175，掉出四维窗。
但如果物理读出是 L=4 窗口（滑动窗），则 L=8 词里所有平衡 4-窗口的
q_rot 应回到 5/9。本脚本核验这一点。
"""
from __future__ import annotations
import json, os, sys
from collections import Counter
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zero_sim as zs  # noqa: E402

RES = os.path.join(HERE, "..", "results")
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


def rot_class(w):
    L = len(w)
    return min(tuple(w[k:] + w[:k]) for k in range(L))


def main():
    W8 = list(zs.balanced_words(8))
    wins = []
    positions = {i: 0 for i in range(5)}
    words_with = 0
    for w in W8:
        got = False
        for i in range(5):
            sub = w[i:i + 4]
            if sum(sub) == 0:
                wins.append(sub)
                positions[i] += 1
                got = True
        words_with += int(got)
    c = Counter(rot_class(w) for w in wins)
    N = sum(c.values())
    q = sum((v / N) ** 2 for v in c.values())
    fh = sum(1 for w in W8 if sum(w[:4]) == 0)
    lh = sum(1 for w in W8 if sum(w[4:]) == 0)
    print(f"W8={len(W8)}; words with >=1 balanced 4-window={words_with}; "
          f"windows={N}; positions={positions}; q_rot={q:.12f}")
    print(f"first half balanced={fh}; last half balanced={lh}")
    ck("每个 W8 词至少一个平衡 4-窗口", words_with == len(W8), f"{words_with}/{len(W8)}")
    ck("窗口数 = 180", N == 180, str(N))
    ck("q_rot(窗口) = 5/9", abs(q - 5 / 9) < 1e-12, f"{q:.12f}")
    ck("五个位置各 36 个窗口", all(v == 36 for v in positions.values()), str(positions))
    ck("首/尾半词平衡各 36", fh == 36 and lh == 36, f"{fh}/{lh}")
    out = {"W8": len(W8), "words_with_window": words_with,
           "windows": N, "positions": positions, "q_rot": q,
           "first_half_balanced": fh, "last_half_balanced": lh,
           "checks": CK, "n_fail": sum(1 for c in CK if not c["pass"])}
    with open(os.path.join(RES, "zero_window_readout.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_window_readout.json"))


if __name__ == "__main__":
    main()
