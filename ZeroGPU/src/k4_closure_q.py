#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""k4_closure_q.py —— L0–L1–L2 三层闭环动力学在 K4 上的稳态纯度 q。

Z088 的 q = Σ_v ω_v²（稳态纯度）；路径图给 q=5/9（两支 2:1）。
本脚本把同一动力学（`layer_closure.LayerClosure`）搭在 K4 上，测稳态 q。

预期：K4 上稳态为均匀四态（ω=1/4），故 q=4·(1/4)²=1/4。
若成立，"四态等权 1/4"就是 L0–L1–L2 动力学的稳态纯度（原生），
而不是统计读法。

输出：ZeroGPU/results/k4_closure_q.json
用法：python3 ZeroGPU/src/k4_closure_q.py
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from layer_closure import LayerClosure  # noqa: E402

OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


adj = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]   # K4
lc = LayerClosure(adj, L=4)
lc.run(600, log_every=200)
ql, qh = lc.q_live(), lc.q_hist()
print(f"   K4 稳态：q_live={ql:.12f}  q_hist={qh:.12f}  （1/4 = {1/4:.12f}）")

check("K4 稳态 q = 1/4（q_live）", abs(ql - 0.25) < 1e-9, f"{ql:.12f}")
check("K4 稳态 q = 1/4（q_hist）", abs(qh - 0.25) < 1e-9, f"{qh:.12f}")
check("1−q = 3/4", abs((1 - ql) - 0.75) < 1e-9, f"{1-ql:.12f}")

print("\n== 汇总 ==")
print(f"   PASS {len(PASS)} / FAIL {len(FAIL)}")
res = {"pass": PASS, "fail": FAIL, "q_live": ql, "q_hist": qh}
with open(os.path.join(OUT, "k4_closure_q.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/k4_closure_q.json")
sys.exit(0 if not FAIL else 1)
