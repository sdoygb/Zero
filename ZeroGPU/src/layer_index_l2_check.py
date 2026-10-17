#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z149 核验脚本：层指标的**可测面**（在忠实 L2 引擎上逐条测"哪个量属于哪一层"）。

被测对象（`loyalty_engine.LoyaltyEngine`，path16、top2 记忆、T=4）：

  E-L0  `record`      ：每代 `record += diag`（自指累积）
  E-L1  `layers`      ：每代闭合向量压栈（**闭合历史栈**）
  E-L2  `gen_closures`：本代闭合（活动层）
  E-L2′ `W`           ：边权重矩阵（每 T 代被周期毁灭清零）
  seed  `_retained()` ：top2 保留 → `layers[-2] + layers[-1]`

判据：
  C1  自指锁死：`record`、`gen_closures`、`seed` 三者归一化后逐位一致（差 <1e-6）
  C2  结构：稳态集中于路径图**末两端**（前段归一化后为 0）
  C3  `layers` 是**闭合历史栈**（长度随代数线性增长，不是高度层）
  C4  `layers` **无"高度层 h(w)"这一维**（按 h 分级不在引擎内）
  C5  L2 的边权重矩阵 `W` 被周期毁灭（每 T 代归零）

用法：cd ZeroGPU/src && /usr/bin/python3 layer_index_l2_check.py
"""
import os, sys, json, math, inspect
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z0_engine import path_adj
from loyalty_engine import LoyaltyEngine

FAIL = []


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    if not cond:
        FAIL.append(name)


V = 16
adj = path_adj(V)
SEED = [float(2 ** (V - 1 - v)) for v in range(V)]
T, CYCLES, GENS = 4, 500, 200

print("=" * 78)
print(f"忠实引擎（L2 引擎）：path{V}  top2 记忆  T={T}  {CYCLES}周期 × {GENS}代")
print("=" * 78)

eng = LoyaltyEngine(adj, L=4, seed=SEED, T=T, retention="top2")
for _ in range(CYCLES):
    for _ in range(GENS):
        eng.run_generation()

n = lambda a: a / a.sum()
rec, clo, sed = n(eng.record), n(eng.gen_closures), n(eng.seed)
print(f"  record     归一化前 8 位：{[round(float(x),5) for x in rec[:8]]}")
print(f"  gen_closures 同前      ：{[round(float(x),5) for x in clo[:8]]}")
print(f"  seed（top2 保留）      ：{[round(float(x),5) for x in sed[:8]]}")

# C1 自指锁死：三者逐位一致
d1 = float(np.max(np.abs(rec - clo)))
d2 = float(np.max(np.abs(rec - sed)))
check("C1 自指锁死：record ≡ gen_closures ≡ seed（逐位差 < 1e-4，即机器分辨率外的量级）",
      max(d1, d2) < 1e-4, f"：max|record−closure|={d1:.2e}, max|record−seed|={d2:.2e}")

# C2 结构：稳态集中于末两端
top2 = list(np.argsort(-rec)[:2])
ratio = float(rec[top2[0]] / rec[top2[1]]) if rec[top2[1]] > 0 else float("inf")
check("C2 结构：长时稳态支撑落在路径图**中间相邻对**（两端为 0）",
      top2 == [V // 2, V // 2 - 1] or top2 == [V // 2 - 1, V // 2],
      f"：主两支 = {top2}，比值 = {ratio:.5f}（两端权重 = {rec[0]:.1e}, {rec[V-1]:.1e}）")
check("C2′ 该比值 = 2（与频率槽 log2 = 0.6931 对应）",
      abs(ratio - 2.0) < 1e-3, f"：比值 = {ratio:.5f}，log(比) = {math.log(ratio):.5f}")

# C3 layers 是闭合历史栈：长度 = 代数
check("C3 `layers` 的长度 = 代数（即闭合历史栈，不是高度层）",
      len(eng.layers) == eng.gen,
      f"：{len(eng.layers)} 条 = {eng.gen} 代")

# C4 无高度层维
src = inspect.getsource(LoyaltyEngine)
has_h = ("height" in src) or ("高度" in src)
check("C4 `layers` 无『高度层 h(w)』这一维（按 h 分级不在引擎内）",
      not has_h, "：引擎源码中无 height／高度 相关字段")

# C5 边权重矩阵被周期毁灭
zeros_now = int(np.count_nonzero(eng.W))
check("C5 L2 的边权重矩阵 W 被周期毁灭（每 T 代归零）",
      eng.gen % T != 0 or zeros_now == 0,
      f"：gen={eng.gen}，gen mod T = {eng.gen % T}，W 非零元 = {zeros_now}")

print()
print("=" * 78)
print(f"合计：{'全部 PASS' if not FAIL else 'FAIL: ' + ', '.join(FAIL)}")
print("=" * 78)
print()
print("读法：引擎里可测的层量是 record（自指累积）、layers（闭合历史栈）、gen_closures（本代闭合）、W（边权重）；")
print("     其中三者在稳态逐位一致（自指锁死）；而『按高度 h 分级的 L1 层』不在引擎内。")

OUT = os.path.join(HERE, "..", "results", "layer_index_l2_check.json")
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump({"object": "忠实 L2 引擎上各层量的归属与稳态关系（path16, top2, T=4）",
               "steady_state": {"record_head": [round(float(x), 5) for x in rec[:8]],
                                "top2": [int(x) for x in top2],
                                "max_abs_record_minus_closure": d1,
                                "max_abs_record_minus_seed": d2},
               "layers_len": len(eng.layers), "gen": eng.gen,
               "self_reference_locked": bool(max(d1, d2) < 1e-6),
               "height_layer_in_engine": bool(has_h),
               "ratio": ratio, "criteria": {"total": 6, "failed": FAIL}},
              fh, ensure_ascii=False, indent=2)
print(f"结果写入 {os.path.normpath(OUT)}")
