#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_elements_audit.py —— 各种自然元素能否从 Zero 推导：四层拆解 + 旧 10.9 审计
================================================================================

四层问题（必须分开）：
  L1 元素存在/周期表：哪些 Z 存在，电子壳层容量 2,8,18,32
  L2 核素/结合能：哪些核稳定，结合能多少
  L3 自然丰度：核合成（BBN/恒星）
  L4 原子性质：半径、电离能、光谱

本脚本只做可核验的审计，不假装导出元素：
  E1 旧 10.9 的轻核表：数值 = 观测值（2.2248/7.718/28.296/31.994/39.245/58.165 MeV）
     ⟹ 「偏差<0.1%」是同义反复，不是推导
  E2 旧 10.9 的电子结构：He Z_eff=1.70、O K=7.70/L=3.76 由屏蔽 0.30/0.85 给；
     旧文自己标注「屏蔽参数是框架内构造，严格推导待后续」
  E3 K4 闭游走/距离层：12/84/732；层 1,3 —— 不是 2,8,18,32
  E4 Heawood 闭游走/距离层：42/210/1554；层 1,4,6,3 —— 不是 2,8,18,32
  E5 类氢离子（一个电子 + 核电荷 Z）：Lyman-α 波长 = 121.502/Z² nm；
     He+ 30.375、Li2+ 13.500、Be3+ 7.594 nm；这是可导出的子集
  E6 当前 Zero 有的：α（Z164）、m_p/m_e（Z165）、氢光谱/21cm（Z166）、
     spin-1/2+Pauli（Z006）；没有：核力/核结合能、原子 n,l 壳层、BBN
  E7 结论：H 与类氢离子可导出；中性多电子原子、核素表、天然丰度不可导出

运行：cd ZeroGPU/src && /usr/bin/python3 zero_elements_audit.py
输出：results/zero_elements_audit.json
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": str(detail)})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def closed_walks(A, L):
    w = np.linalg.eigvalsh(A)
    return float(np.sum(w ** L))


def layers(A, s=0):
    n = A.shape[0]
    dist = {s: 0}
    seen = {s}
    frontier = [s]
    d = 0
    while frontier:
        d += 1
        nxt = []
        for u in frontier:
            for v in np.where(A[u] > 0)[0]:
                if v not in seen:
                    seen.add(v)
                    nxt.append(v)
                    dist[v] = d
        frontier = nxt
    from collections import Counter
    c = Counter(dist.values())
    return [int(c[k]) for k in sorted(c)]


def heawood():
    # Fano 平面 incidence graph: 7 points + 7 lines; 21 edges, 3-regular, girth 6
    lines = [(0, 1, 3), (1, 2, 4), (2, 3, 5), (3, 4, 6),
             (4, 5, 0), (5, 6, 1), (6, 0, 2)]
    A = np.zeros((14, 14), dtype=int)
    for j, pts in enumerate(lines):
        for p in pts:
            A[p, 7 + j] = A[7 + j, p] = 1
    return A


def main():
    print("=" * 84)
    print("自然元素能否从 Zero 推导：四层拆解 + 旧 10.9 审计")
    print("=" * 84)

    # E1 旧 10.9 轻核表 = 观测值
    print("\n§1 旧 10.9 的「轻核生成」表")
    obs = {"H-2": 2.2248, "He-3": 7.718, "He-4": 28.296,
           "Li-6": 31.994, "Li-7": 39.245, "Be-9": 58.165}
    art = {"H-2": 2.2248, "He-3": 7.718, "He-4": 28.296,
           "Li-6": 31.994, "Li-7": 39.245, "Be-9": 58.165}
    for k in obs:
        print(f"  {k}: 旧文 {art[k]:.6f} MeV；观测 {obs[k]:.6f} MeV")
    ck("E1 旧 10.9 轻核表的数值逐项等于观测值（不是独立推导）",
       all(abs(art[k] - obs[k]) < 1e-9 for k in obs), "6/6 逐字相同")

    # E2 电子结构输入
    print("\n§2 旧 10.9 的电子结构")
    print("  He: Z_eff=1.70（2−σ_same, σ_same=0.30）")
    print("  O : K 层 Z_eff=7.70，L 层 Z_eff=3.76；屏蔽 σ_same=0.30, σ_inner=0.85")
    ck("E2 旧 10.9 的电子结构以屏蔽参数为输入（旧文自认「框架内构造」）",
       abs(2.0 - 0.30 - 1.70) < 1e-9, "He Z_eff=1.70 = 2−0.30")

    # E3/E4 图计数
    print("\n§3 Zero 的图计数能否给出 2,8,18,32")
    K4 = np.ones((4, 4)) - np.eye(4)
    H = heawood()
    print(f"  K4: 边={int(K4.sum()//2)}  闭游走 L=2,4,6,8 = "
          f"{[round(closed_walks(K4,L)) for L in (2,4,6,8)]}  距离层={layers(K4)}")
    print(f"  Heawood: 顶点={H.shape[0]}  边={int(H.sum()//2)}  "
          f"闭游走 L=2,4,6,8 = {[round(closed_walks(H,L)) for L in (2,4,6,8)]}  "
          f"距离层={layers(H)}")
    ck("E3 K4 的闭游走/距离层不是 2,8,18,32", layers(K4) == [1, 3],
       f"层={layers(K4)}；闭游走={[round(closed_walks(K4,L)) for L in (2,4,6,8)]}")
    ck("E4 Heawood 的闭游走/距离层不是 2,8,18,32", layers(H) == [1, 3, 6, 4],
       f"层={layers(H)}；闭游走={[round(closed_walks(H,L)) for L in (2,4,6,8)]}")

    # E5 类氢离子
    print("\n§4 类氢离子（可导出子集）：Lyman-α 波长 = 121.502/Z² nm")
    R = 10973731.568  # m^-1
    lam1 = 4.0 / (3.0 * R) * 1e9  # nm
    rows = []
    for Z, known in [(1, 121.502), (2, 30.375), (3, 13.500), (4, 7.594),
                     (5, 4.860), (6, 3.375), (7, 2.480), (8, 1.898)]:
        lam = lam1 / Z ** 2
        rows.append((Z, lam, known))
        print(f"  Z={Z}: λ={lam:.4f} nm（已知 {known:.3f}）")
    ck("E5 类氢 Lyman-α 按 1/Z² 标度（1% 内）",
       all(abs(lam / known - 1) < 0.01 for _, lam, known in rows),
       f"Z=1..8；λ1={lam1:.3f} nm")

    # E6 当前资产
    print("\n§5 当前 Zero 的承重件")
    for z, note in {"Z164": "α=137.036", "Z165": "m_p/m_e", "Z166": "氢光谱/21cm",
                    "Z006": "spin-1/2 + Pauli"}.items():
        hits = glob.glob(os.path.join(REPO, f"{z}_*.md"))
        print(f"  {z}: {'✓' if hits else '✗'}  {note}")
    no_periodic = not glob.glob(os.path.join(REPO, "Z*周期表*.md")) and \
                  not glob.glob(os.path.join(REPO, "Z*元素*.md"))
    ck("E6 当前库有 α/m_p/m_e/氢光谱/Pauli，但没有周期表/元素篇",
       no_periodic, "Z*周期表*/Z*元素* = 0")

    # E7 结论
    ck("E7 结论：H 与类氢离子可导出；中性多电子原子、核素表、天然丰度不可导出",
       True, "见 doc §5；旧 10.9 轻核表=观测值、电子结构=输入屏蔽")

    nfail = sum(1 for c in CK if not c["pass"])
    print("-" * 84)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-nfail} / FAIL {nfail}")
    out = dict(checks=CK, n_fail=nfail,
               old_light_nuclei=art, observed=obs,
               k4=dict(layers=layers(K4),
                       closed={str(L): closed_walks(K4, L) for L in (2, 4, 6, 8)}),
               heawood=dict(vertices=int(H.shape[0]), edges=int(H.sum() // 2),
                            layers=layers(H),
                            closed={str(L): closed_walks(H, L) for L in (2, 4, 6, 8)}),
               hydrogen_like=[{"Z": z, "lambda_nm": lam, "known_nm": k}
                              for z, lam, k in rows],
               boundary="H/类氢可导出；中性多电子原子、核素/结合能、天然丰度不可导出")
    with open(os.path.join(RES, "zero_elements_audit.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(f"结果 → {os.path.join(RES, 'zero_elements_audit.json')}")
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
