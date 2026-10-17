#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""chain_systematics.py —— 半槽链对剩余建模选择的稳定性（推导稿 §4）。

剩余的选择：
  D1 的精确形式：1/(N+1/2)  vs  (1/N)(1−1/(2N))；
  D3 的形式：N/(N+1/2)=32/33  vs  (1−1/(2N))=31/32。
本脚本列出四种组合对 m_ν3、Σ、mββ 与 Δm²31 的 pull 的影响——检验链的预言是否稳定。

输出：ZeroGPU/results/chain_systematics.json
用法：python3 ZeroGPU/src/chain_systematics.py
"""
from __future__ import annotations

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


V, M_Z = 246.22, 91.1876
M_I = M_Z * math.exp(4 * math.pi)
M_X = 1.009e13
M_R0 = math.sqrt(M_I * M_X)
A = math.sqrt(2)
d0 = math.pi / 12 - 2 / 9
de = 3 * math.pi / 4 - d0
f = [1 + A * math.cos(de + 2 * math.pi * k / 3) for k in range(3)]
R = (73 - 28 * math.sqrt(6)) / 25

N1 = 784          # 质量单元的通道数
N3 = 16           # M_R 的通道数

d1_variants = {
    "1/(N+1/2)": 1 / (N1 + 0.5),
    "(1/N)(1−1/(2N))": (1 / N1) * (1 - 1 / (2 * N1)),
}
d3_variants = {
    "N/(N+1/2)=32/33": N3 / (N3 + 0.5),
    "1−1/(2N)=31/32": 1 - 1 / (2 * N3),
}
rows = []
for d1name, c in d1_variants.items():
    mtau = V * c * f[2] ** 2
    for d3name, alpha in d3_variants.items():
        MR = M_R0 * alpha
        m3 = (mtau / 2) ** 2 / MR * 1e12
        m2 = R * m3
        Sigma = m2 + m3
        mbb = 0.5443**2 * m2 + 0.1492**2 * m3
        d31 = (m3 * 1e-3) ** 2
        pull = (d31 - 2.510e-3) / 0.027e-3
        rows.append({"D1": d1name, "D3": d3name, "m_tau": mtau, "m_nu3": m3,
                     "Sigma": Sigma, "m_bb": mbb, "d31_pull": pull})
        print(f"   {d1name:18s} × {d3name:16s}: m_τ={mtau:.6f}  m_ν3={m3:.3f} meV  Σ={Sigma:.2f}  pull={pull:+.2f}σ")

spread = (max(r["m_nu3"] for r in rows) - min(r["m_nu3"] for r in rows)) / sum(r["m_nu3"] for r in rows) * len(rows)
print(f"\n   m_ν3 的跨选项散布: {100*spread:.3f}%")
check("四种组合的 m_ν3 散布 <0.5%", spread < 0.005, f"{100*spread:.3f}%")
check("四种组合的 Δm²31 pull 都在 ±1σ 内", all(abs(r["d31_pull"]) < 1 for r in rows),
      f"{[round(r['d31_pull'],2) for r in rows]}")
check("Σ 都 <64 meV", all(r["Sigma"] < 64 for r in rows),
      f"{[round(r['Sigma'],1) for r in rows]}")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "rows": rows, "spread_mnu3": spread}
with open(os.path.join(OUT, "chain_systematics.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/chain_systematics.json")
sys.exit(0 if not FAIL else 1)
