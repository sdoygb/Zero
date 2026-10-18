#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""amplitude_closed_form_test.py —— 上/下型振幅闭式候选：同一把尺上的取舍。

尺：质量律 m = u_s A² M(θ)² 里进入的是 **A²**，所以候选闭式以 A² 为对象、以 A² 为尺比。
参考值取库内谱拟合（`quark_phase_fit.py` 的 A），本脚本只做**候选比较**，不重拟合。

背景：每层 3 个质量由 (A, θ, u_s) 三个量定，拟合是**恰好定解**（无余量），
所以振幅没有"预言残差"可作判据——唯一可用的判据是闭式对拟合值的残差 + 结构可读性。

判据语义：FAIL = 表中某一步不成立（要报的判定）。
输出：ZeroGPU/results/amplitude_closed_form_test.json
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


# 库内谱拟合值（quark_phase_fit.py → results/quark_phase_fit.json）
A_FIT = {"up": 1.7778326303083105, "down": 2.682849902237757}

CANDS = {
    "up": {
        "(4/3)^4 = 256/81": (4 / 3) ** 4,
        "79/25": 79 / 25,
        "√10": math.sqrt(10),
    },
    "down": {
        "36/5": 36 / 5,
        "(8/3)² = 64/9": (8 / 3) ** 2,
        "2π+1": 2 * math.pi + 1,
        "(4/3)^4 · (8/3)": (4 / 3) ** 4 * (8 / 3),
    },
}

rows = []
print("=" * 84)
print("上/下型振幅闭式候选：以 A² 为尺")
print("=" * 84)
for sector in ("up", "down"):
    A2 = A_FIT[sector] ** 2
    print(f"\n  {sector} 层：拟合 A = {A_FIT[sector]:.9f}，A² = {A2:.9f}")
    for name, c2 in sorted(CANDS[sector].items(), key=lambda kv: abs(kv[1] / A2 - 1)):
        rel = 100 * (c2 / A2 - 1)
        rows.append(dict(sector=sector, name=name, A2_cand=c2, A2_fit=A2, rel_pct=rel))
        print(f"    A² = {c2:.9f}（{name:22s}）  残差 {rel:+.4f}%")

up = [r for r in rows if r["sector"] == "up"]
dn = [r for r in rows if r["sector"] == "down"]
bu = min(up, key=lambda r: abs(r["rel_pct"]))
bd = min(dn, key=lambda r: abs(r["rel_pct"]))
d83 = [r for r in dn if r["name"].startswith("(8/3)")][0]
check("上型最佳候选 = (4/3)^4 = 256/81，残差 <0.01%",
      bu["name"].startswith("(4/3)^4") and abs(bu["rel_pct"]) < 0.01,
      f"{bu['name']} {bu['rel_pct']:+.4f}%")
check("下型最佳候选 = 36/5，残差 <0.05%；(8/3)² 差 >1%",
      bd["name"] == "36/5" and abs(bd["rel_pct"]) < 0.05 and abs(d83["rel_pct"]) > 1.0,
      f"36/5 {bd['rel_pct']:+.4f}%，(8/3)² {d83['rel_pct']:+.3f}%")
check("A 尺上的 8/3（0.60%）与 A² 尺上的 (8/3)²（1.2%）是同一个候选的两种写法——尺换一次残差翻倍",
      abs(100 * ((8 / 3) / A_FIT["down"] - 1)) < 0.7
      and abs(100 * ((8 / 3) ** 2 / A_FIT["down"] ** 2 - 1)) > 1.0,
      "0.60% → 1.20%")
print()
print("  判定：")
print("  · 下型：8/3 在 A 尺上 0.60%、在 A² 尺上 1.20%——换尺翻倍，说明它只是「看着近」；退役，取 36/5。")
print("  · 上型：256/81 = (4/3)^4 残差 −0.006%，且 4/3 是亚层比的公比（有结构来源）；79/25 无读法。")
print("  · 于是振幅定在 A² 尺上：A_up² = (4/3)^4、A_down² = 36/5。")

res = {"A_fit": A_FIT, "rows": rows, "pass": PASS, "fail": FAIL}
with open(os.path.join(OUT, "amplitude_closed_form_test.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print()
print("=" * 84)
print(f"汇总：PASS {len(PASS)} / FAIL {len(FAIL)}")
if FAIL:
    print("  FAIL:", FAIL)
print("  → results/amplitude_closed_form_test.json")
sys.exit(0 if not FAIL else 1)
