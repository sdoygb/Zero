#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""neutrino_suite.py —— 中微子程序的「一条命令」全套核验。

依次运行 11 个中微子脚本，汇总退出码与 PASS/FAIL。任一失败则整体失败。
用途：goal 的「可复跑核验」——一个入口跑完全部检查（约数秒）。

输出：ZeroGPU/results/neutrino_suite.json
用法：python3 ZeroGPU/src/neutrino_suite.py
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)

SCRIPTS = [
    "neutrino_lambda_amplitude.py",
    "neutrino_amplitude_geometry.py",
    "neutrino_amplitude_spectrum.py",
    "neutrino_scale_conditional.py",
    "neutrino_branch_test.py",
    "neutrino_predictions.py",
    "neutrino_shape_residual.py",
    "neutrino_tension_decomposition.py",
    "neutrino_scale_sensitivity.py",
    "neutrino_prediction_band.py",
    "neutrino_half_factor_rge.py",
    "neutrino_independent_check.py",
]

rows = []
for name in SCRIPTS:
    path = os.path.join(HERE, name)
    if not os.path.exists(path):
        rows.append({"script": name, "ok": False, "detail": "missing"})
        print(f"  [x] {name}: 缺文件")
        continue
    r = subprocess.run([sys.executable, path], capture_output=True, text=True)
    tail = [ln for ln in r.stdout.strip().splitlines() if ln.strip()]
    summary = next((ln.strip() for ln in tail if "PASS" in ln and "FAIL" in ln), "")
    ok = r.returncode == 0
    rows.append({"script": name, "ok": ok, "returncode": r.returncode,
                 "summary": summary})
    print(f"  [{'v' if ok else 'x'}] {name}: {summary or ('rc=' + str(r.returncode))}")

n_ok = sum(1 for r in rows if r["ok"])
print(f"\n== 汇总：{n_ok}/{len(rows)} 脚本通过 ==")

res = {"scripts": rows, "n_ok": n_ok, "n_total": len(rows)}
with open(os.path.join(OUT, "neutrino_suite.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print("   → results/neutrino_suite.json")
sys.exit(0 if n_ok == len(rows) else 1)
