#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_mode_energy_check.py —— 相位零模能量路线对 c_ℓ 的检验（Z166 §2.2）。

路线：质量单元 A² 是相位偏离振动在零点（3π/4）的零模能量。
取 V_eff(θ)=(1+√2cosθ)²（质量²的势），零点曲率 V''=2。若 A²=v·V''/N（N 为结构惯性），
则 c_ℓ = A²/v = V''/N。反解：N = V''/c_ℓ = 2/c_ℓ。

检查：所需 N 与干净计数 2·C(8,2)²=2·784=1568 的关系（是否独立给出新结构）。

输出：ZeroGPU/results/zero_mode_energy_check.json
用法：python3 ZeroGPU/src/zero_mode_energy_check.py
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


V, M_TAU = 246.22, 1.77686
A_L = math.sqrt(2)
delta0 = math.pi / 12 - 2 / 9
delta_e = 3 * math.pi / 4 - delta0
f_tau = max(1 + A_L * math.cos(delta_e + 2 * math.pi * k / 3) for k in range(3))
c_ell = M_TAU / (V * f_tau**2)

# 零点曲率：V_eff=(1+√2cosθ)²，在 M=0 处 V''=2
theta0 = 3 * math.pi / 4
curv = 2.0 * (math.sqrt(2) * math.sin(theta0)) ** 2 + 0.0
print(f"   c_ℓ={c_ell:.9e}；零点曲率 V''={curv:.6f}（解析）")
N_req = curv / c_ell
print(f"   所需惯性 N=V''/c_ℓ = {N_req:.3f}")
for name, val in [("2·C(8,2)²=1568", 2 * math.comb(8, 2) ** 2),
                  ("C(8,2)²=784", math.comb(8, 2) ** 2),
                  ("16·7²=784", 16 * 49),
                  ("2·16·7²=1568", 2 * 16 * 49)]:
    print(f"     {name:16s} = {val:8.0f}   与 N_req 差 {100*(val/N_req-1):+.3f}%")
check("所需 N 与 2·784 差 <0.1%", abs(2 * 784 / N_req - 1) < 1e-3,
      f"{100*(2*784/N_req-1):+.3f}%")
check("零模路线不独立：所需 N 就是 2·C(8,2)²", abs(N_req - 2 * 784) / N_req < 1e-3)

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "c_ell": c_ell, "curvature": curv,
       "N_required": N_req, "N_count": 2 * 784}
with open(os.path.join(OUT, "zero_mode_energy_check.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/zero_mode_energy_check.json")
sys.exit(0 if not FAIL else 1)
