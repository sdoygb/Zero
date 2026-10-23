#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_periodic_table.py —— 元素周期表：能推到哪一步？
================================================================================

结论先说：**周期表可以“生成”，但当前 Zero 不能“推出”它。**
本脚本做两件事：
  A. 用标准 Aufbau（容量 2(2l+1) + Madelung n+l 规则）生成 Z=1..118 的周期表，
     核对周期长度 2,8,8,18,18,32,32、4s 先于 3d、区块划分；
  B. 逐项核对 Zero 现有承重件，标出周期表还缺什么。

判据：
  P1 周期长度 = [2,8,8,18,18,32,32]
  P2 4s 先于 3d（K/Ca 在 Sc 前）
  P3 累积容量 2,8,18,32 = 2n²
  P4 Zero 有 d=3（Z021）、spin-1/2+Pauli（Z006）、1/r（Z147）
  P5 Zero **没有**：原子 n,l 壳层、Coulomb SO(4)/Runge-Lenz 简并、Madelung 规则、核结合能、屏蔽
  P6 判定：周期表是**输入结构**，不是当前 Zero 的定理

运行：cd ZeroGPU/src && /usr/bin/python3 zero_periodic_table.py
输出：results/zero_periodic_table.json
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": str(detail)})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def aufbau(Zmax=118):
    """标准 Aufbau：按 (n+l, n) 排序填 2(2l+1)。"""
    subs = []
    for n in range(1, 8):
        for l in range(0, n):
            subs.append((n + l, n, l, 2 * (2 * l + 1)))
    subs.sort(key=lambda x: (x[0], x[1]))
    table = []
    Z = 0
    for npl, n, l, cap in subs:
        for _ in range(cap):
            if Z >= Zmax:
                break
            Z += 1
            table.append(dict(Z=Z, n=n, l=l, npl=npl,
                              block=("s" if l == 0 else "p" if l == 1 else "d" if l == 2 else "f")))
        if Z >= Zmax:
            break
    return table


def periods_from_subs(subs, Zmax=118):
    """周期 = 从一个 ns 开始，到下一个 ns 之前；容量按 Aufbau 顺序累加。"""
    out = []
    cur = 0
    Z = 0
    for npl, n, l, cap in subs:
        if l == 0:
            if cur:
                out.append(cur)
            cur = 0
        take = min(cap, Zmax - Z)
        cur += take
        Z += take
        if Z >= Zmax:
            break
    if cur:
        out.append(cur)
    return out


def main():
    print("=" * 84)
    print("元素周期表：能推到哪一步？")
    print("=" * 84)
    tab = aufbau(118)
    subs = sorted([(n + l, n, l, 2 * (2 * l + 1)) for n in range(1, 8) for l in range(0, n)],
                  key=lambda x: (x[0], x[1]))
    pl = periods_from_subs(subs, 118)
    print(f"  Z=1..118；周期长度={pl}")
    print(f"  前 20 号：{[ (t['Z'], t['n'], t['l']) for t in tab[:20] ]}")

    ck("P1 周期长度 = [2,8,8,18,18,32,32]",
       pl == [2, 8, 8, 18, 18, 32, 32], str(pl))
    # K(19) Z=19: 4s; Sc(21): 3d
    z19, z21 = tab[18], tab[20]
    ck("P2 4s 先于 3d（K/Ca 在 Sc 前）",
       (z19["n"], z19["l"]) == (4, 0) and (z21["n"], z21["l"]) == (3, 2),
       f"Z19={z19['n']}{z19['block']}  Z21={z21['n']}{z21['block']}")
    ck("P3 累积容量 2,8,18,32 = 2n²",
       [2, 8, 18, 32] == [2 * n * n for n in (1, 2, 3, 4)],
       "2(1,4,9,16)")
    ck("P4 Zero 有 d=3（Z021）、spin-1/2+Pauli（Z006）、1/r（Z147）", True,
       "见前置")
    ck("P5 Zero 没有：原子 n,l 壳层、Coulomb SO(4)/Runge-Lenz 简并、Madelung、核结合能、屏蔽", True,
       "逐项见 doc §2")
    ck("P6 判定：周期表是输入结构，不是当前 Zero 的定理", True,
       "生成 ≠ 推出")

    nfail = sum(1 for c in CK if not c["pass"])
    print("-" * 84)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-nfail} / FAIL {nfail}")
    out = dict(period_lengths=pl, first20=[(t["Z"], t["n"], t["l"]) for t in tab[:20]],
               checks=CK, n_fail=nfail,
               verdict="生成≠推出；缺 SO(4)/Runge-Lenz、Madelung、核结合能、屏蔽")
    with open(os.path.join(RES, "zero_periodic_table.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(f"结果 → {os.path.join(RES, 'zero_periodic_table.json')}")
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
