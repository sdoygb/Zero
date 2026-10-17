#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""quark_phase_fit.py —— 夸克质量谱的相位律拟合（推导稿：跨扇区统一）。

质量律（Z121）：m = A_s² · M(θ)²，M(θ)=1+A_s cosθ；三代按 2π/3 间隔。
对每个夸克扇区（up: u,c,t；down: d,s,b），用两个 √m 比值定 (θ, A_s)，
再检查绝对标度相对轻子单元 A²=v/784.5 的因子。

输出：ZeroGPU/results/quark_phase_fit.json
用法：python3 ZeroGPU/src/quark_phase_fit.py
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


# PDG MS-bar（MeV）：u=2.16, d=4.67, s=93.4, c=1270, b=4180, t=172690
Q = {"u": 2.16, "d": 4.67, "s": 93.4, "c": 1270.0, "b": 4180.0, "t": 172690.0}
A2 = 246.22 / 784.5   # 轻子质量单元（GeV）

sectors = [("up", ["u", "c", "t"], 2.10, 2.00), ("down", ["d", "s", "b"], 2.00, 2.70)]


def solve(r2, r3, th, A):
    """牛顿解：log(M2/M1)=log r2, log(M3/M1)=log r3。"""
    def M(th, A):
        return [abs(1 + A * math.cos(th + 2 * math.pi * k / 3)) for k in range(3)]
    for _ in range(60):
        f0 = M(th, A)
        F = [math.log(f0[1] / f0[0] / r2), math.log(f0[2] / f0[0] / r3)]
        h = 1e-7
        fth, fA = M(th + h, A), M(th, A + h)
        J = [[(fth[1] / f0[1] - fth[0] / f0[0]) / h, (fA[1] / f0[1] - fA[0] / f0[0]) / h],
             [(fth[2] / f0[2] - fth[0] / f0[0]) / h, (fA[2] / f0[2] - fA[0] / f0[0]) / h]]
        det = J[0][0] * J[1][1] - J[0][1] * J[1][0]
        if abs(det) < 1e-16:
            break
        dth = (-F[0] * J[1][1] + F[1] * J[0][1]) / det
        dA = (-J[0][0] * F[1] + J[1][0] * F[0]) / det
        th += dth
        A += dA
    return th, A


rows = []
for name, qs, th0, A0 in sectors:
    m = [Q[q] for q in qs]
    r2, r3 = math.sqrt(m[1] / m[0]), math.sqrt(m[2] / m[0])
    th, A = solve(r2, r3, th0, A0)
    M = [abs(1 + A * math.cos(th + 2 * math.pi * k / 3)) for k in range(3)]
    pred = [A2 * M[k] ** 2 * 1e3 for k in range(3)]   # MeV，轻子单元
    factors = [m[k] / pred[k] for k in range(3)]
    print(f"   {name}: θ={math.degrees(th):.4f}°  A={A:.6f}  A²={A*A:.6f}")
    print(f"      预测(轻子单元)={[round(p,3) for p in pred]} MeV  实际={m}")
    print(f"      单元因子={[round(f,3) for f in factors]}")
    rows.append({"sector": name, "theta_deg": math.degrees(th), "A": A, "A2": A * A,
                 "pred_unit_MeV": pred, "measured_MeV": m, "unit_factors": factors})

up = rows[0]["unit_factors"]
dn = rows[1]["unit_factors"]
check("上型单元因子一致（散差 <2%）", (max(up) - min(up)) / sum(up) * 3 < 0.02,
      f"{[round(f,2) for f in up]}")
check("下型单元因子≈1（轻子单元，<3%）", all(abs(f - 1) < 0.03 for f in dn),
      f"{[round(f,3) for f in dn]}")
check("上型因子≈72（<1%）", all(abs(f / 72 - 1) < 0.01 for f in up),
      f"{[round(f/72,4) for f in up]}")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "rows": rows, "A2_GeV": A2}
with open(os.path.join(OUT, "quark_phase_fit.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/quark_phase_fit.json")
sys.exit(0 if not FAIL else 1)
