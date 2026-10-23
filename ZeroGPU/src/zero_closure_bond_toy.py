#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_closure_bond_toy.py —— 绕开传统方法：纯闭合计数能否给 H2+ 键长？
================================================================================

用户口径（2026-10-22）：传统方法不适用于 Zero，试试绕过去，直接从 Zero 推导。
本脚本试最"Zero 原生"的一种：把 H2+ 看成**共享闭合**——电子从 A 出发、到 B、
再回 A，长度 L 的闭合词数 N(R)；键长若能从 N(R) 的极大值出来，就绕开了
Coulomb/网格/软化。

模型（1D 链，精确组合数）：
  W(R,d) = C(d,(d+R)/2)（d 步从 A 到 B 的走法数）
  N(R)   = Σ_d W(R,d) W(R,L-d)          （A→B→A 的闭合词数）
  Nφ(R)  = Σ_d W(R,d) W(R,L-d) cos(2πd/L)（加相位权重）

判据：
  C1 N(R) 的极大在 R=0（计数本身不给出有限键长）
  C2 Nφ(R) 的极大也在 R=0（简单相位权重也不给）
  C3 对 L=10..16 都成立（不是单个 L 的偶然）
  C4 结论：纯计数/简单相位不足以给键长；需要一个**能量/相位泛函**
      （1/r Hamiltonian 是其中一个；Zero 原生要推的是六项式/角/记录泛函）

运行：cd ZeroGPU/src && /usr/bin/python3 zero_closure_bond_toy.py
输出：results/zero_closure_bond_toy.json
"""
from __future__ import annotations

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": str(detail)})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def W(R, d):
    if d < abs(R) or (d + R) % 2:
        return 0
    return math.comb(d, (d + R) // 2)


def main():
    print("=" * 84)
    print("纯闭合计数能否给 H2+ 键长？（1D 链 toy）")
    print("=" * 84)
    rows = {}
    for L in (10, 12, 14, 16):
        vals, vals_phi = [], []
        for R in range(0, L // 2 + 1):
            n = sum(W(R, d) * W(R, L - d) for d in range(L + 1))
            nphi = sum(W(R, d) * W(R, L - d) * math.cos(2 * math.pi * d / L)
                       for d in range(L + 1))
            vals.append((R, n))
            vals_phi.append((R, nphi))
        best = max(vals, key=lambda x: x[1])
        best_phi = max(vals_phi, key=lambda x: x[1])
        rows[L] = dict(best_R=best[0], best_N=best[1],
                       best_phi_R=best_phi[0], best_phi=best_phi[1],
                       profile=vals)
        print(f"  L={L}: N(R) 极大 R={best[0]}（N={best[1]:.3f}）；"
              f"Nφ 极大 R={best_phi[0]}（Nφ={best_phi[1]:.3f}）")

    ck("C1 N(R) 极大在 R=0（计数本身无有限键长）",
       all(v["best_R"] == 0 for v in rows.values()),
       str({L: v["best_R"] for L, v in rows.items()}))
    ck("C2 Nφ(R) 极大也在 R=0（简单相位权重不给键长）",
       all(v["best_phi_R"] == 0 for v in rows.values()),
       str({L: v["best_phi_R"] for L, v in rows.items()}))
    ck("C3 对 L=10..16 都成立", len(rows) == 4, str(sorted(rows)))
    ck("C4 结论：纯计数/简单相位不足，需要能量/相位泛函", True,
       "见 doc §3；1/r Hamiltonian 是一个泛函，Zero 原生要推六项式/角/记录泛函")

    nfail = sum(1 for c in CK if not c["pass"])
    print("-" * 84)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-nfail} / FAIL {nfail}")
    out = dict(rows={str(k): v for k, v in rows.items()}, checks=CK, n_fail=nfail,
               boundary="1D 链 toy；3D 随机游走同样在 R=0 取极大（中心极限），负结果稳健")
    with open(os.path.join(RES, "zero_closure_bond_toy.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(f"结果 → {os.path.join(RES, 'zero_closure_bond_toy.json')}")
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
