#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_hydrogen_like.py —— 类氢离子光谱：把 Z166 的氢结果按 1/Z² 扩到 He⁺、Li²⁺、Be³⁺…
================================================================================

这是元素审计里最容易补的缺口：
  · 一个电子 + 核电荷 Z（类氢离子）不需要核结合能、不需要多电子壳层；
  · 能级 E_n = -13.5983 eV · Z²/n²（含约化质量 μ），波长 1/λ = R_μ Z² (1/n_lo² - 1/n_hi²)；
  · 与 Z166 的氢链同一套 α/m_e/m_p，只是把 m_N 换成对应核质量。

核验：
  H1  H Lyman-α 121.567 nm（Z=1）
  H2  He⁺ 30.375 nm；Li²⁺ 13.500；Be³⁺ 7.594（1% 内）
  H3  固定跃迁 Z²λ 近常数（约化质量造成的偏离 <0.1%）
  H4  约化质量修正：H 的 λ 比无限质量值长 m_e/m_N ≈ 0.0545%
  H5  Balmer-α：H 656.28；He⁺ 164.0；Li²⁺ 72.9；Be³⁺ 41.0 nm（1% 内）

运行：cd ZeroGPU/src && /usr/bin/python3 zero_hydrogen_like.py
输出：results/zero_hydrogen_like.json
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


# Zero 基础：α 用 Z164 闭式 137.036；m_e 用 Z126/Z164 的电子质量；
# c、h 是 SI 单位换算（Z166 同款），不是物理输入。
ALPHA_INV = 137.036          # Z164
ALPHA = 1.0 / ALPHA_INV
C = 299792458.0              # SI 换算
H_PLANCK = 6.62607015e-34    # SI 换算
ME = 9.1093837015e-31        # Z126/Z164 电子质量（SI 实现）
RINF = ALPHA * ALPHA * ME * C / (2.0 * H_PLANCK)  # m^-1，Zero 侧尺度
U = 1.66053906660e-27        # kg，核质量单位（SI 换算）

# (name, Z, atomic mass number A)
IONS = [("H", 1, 1.00784), ("He+", 2, 4.002602), ("Li2+", 3, 6.941),
        ("Be3+", 4, 9.0121831), ("B4+", 5, 10.81), ("C5+", 6, 12.011)]

# 已知真空波长（nm），只作对照
KNOWN = {
    ("H", 1, 2): 121.567, ("He+", 1, 2): 30.375, ("Li2+", 1, 2): 13.500,
    ("Be3+", 1, 2): 7.594,
    ("H", 2, 3): 656.28, ("He+", 2, 3): 164.0, ("Li2+", 2, 3): 72.9,
    ("Be3+", 2, 3): 41.0,
}


def wavelengths(Z, A):
    mN = A * U
    mu = ME * mN / (ME + mN)
    RM = RINF * mu / ME
    out = {}
    for n_lo in (1, 2):
        for n_hi in range(n_lo + 1, n_lo + 4):
            inv = RM * Z * Z * (1.0 / n_lo ** 2 - 1.0 / n_hi ** 2)
            out[(n_lo, n_hi)] = 1e9 / inv
    return out, mu, RM


def main():
    print("=" * 84)
    print("类氢离子光谱：1/Z² 标度 + 约化质量")
    print("=" * 84)
    rows = {}
    for name, Z, A in IONS:
        w, mu, RM = wavelengths(Z, A)
        rows[name] = dict(Z=Z, A=A, mu=mu, R=RM,
                          lam={f"{a}->{b}": v for (a, b), v in w.items()})
        print(f"  {name:6s} Z={Z}  μ={mu:.6e} kg  "
              f"Lyman-α={w[(1,2)]:.4f} nm  Balmer-α={w[(2,3)]:.4f} nm")

    # H1/H2
    for key in [("H", 1, 2), ("He+", 1, 2), ("Li2+", 1, 2), ("Be3+", 1, 2)]:
        name = key[0]
        got = rows[name]["lam"][f"{key[1]}->{key[2]}"]
        ck(f"H1 {name} Lyman-α 与已知一致（1%）", abs(got / KNOWN[key] - 1) < 0.01,
           f"{got:.4f} vs {KNOWN[key]:.3f} nm")

    # H3 Z²λ 常数（同跃迁）
    for (a, b) in [(1, 2), (2, 3)]:
        vals = []
        for name, Z, A in IONS:
            vals.append(rows[name]["lam"][f"{a}->{b}"] * Z * Z)
        ck(f"H3 固定跃迁 {a}->{b}：Z²λ 近常数（相对散布 <0.5%）",
           (max(vals) / min(vals) - 1) < 0.005,
           f"min={min(vals):.3f}, max={max(vals):.3f} nm")

    # H4 约化质量修正（H）
    lam_inf = 1e9 / (RINF * (1.0 - 0.25))
    lam_H = rows["H"]["lam"]["1->2"]
    corr = lam_H / lam_inf - 1.0
    ck("H4 H 的约化质量修正 ≈ m_e/m_N = 0.0545%（0.05%–0.06%）",
       0.0005 < corr < 0.0006, f"{corr*100:.4f}%")

    # H5 Balmer
    for key in [("H", 2, 3), ("He+", 2, 3), ("Li2+", 2, 3), ("Be3+", 2, 3)]:
        name = key[0]
        got = rows[name]["lam"][f"{key[1]}->{key[2]}"]
        ck(f"H5 {name} Balmer-α 与已知一致（1%）", abs(got / KNOWN[key] - 1) < 0.01,
           f"{got:.3f} vs {KNOWN[key]:.2f} nm")

    nfail = sum(1 for c in CK if not c["pass"])
    print("-" * 84)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-nfail} / FAIL {nfail}")
    out = dict(ions=rows, checks=CK, n_fail=nfail,
               boundary="类氢离子=一个电子+核电荷 Z；中性多电子原子仍需核结合能与原子 n,l 壳层")
    with open(os.path.join(RES, "zero_hydrogen_like.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(f"结果 → {os.path.join(RES, 'zero_hydrogen_like.json')}")
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
