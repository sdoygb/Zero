#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""l2_sublayer_mass_check.py —— 「质量与 L2 三个亚层」的检验（用户猜想）。

L2 的亚层结构（Z129／Z089／Z091）：
  · L2 分 3 层（K₄ 的 3 频率槽 = 3 代）；
  · 忠实引擎稳态 = seed 归一化，seed = [8, 4, 1.5, 1]（和 14.5）；
  · 三个相邻比（升序）= (3/2, 2, 8/3)：等比，公比 4/3；
  · 频率槽 = 三个比的 log。

本脚本检验三件事：
  A  seed 与三个比的结构（等比、公比、频率槽）
  B  库内声称的「3 代质量（归一）」是否独立（Z091 / k4_nuR_features.py）
  C  四个扇区振幅与三个比的关系

判据语义：FAIL = 表中某一步不成立（要报的判定）。
输出：ZeroGPU/results/l2_sublayer_mass_check.json
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


# ───────────────── A  亚层结构 ─────────────────
print("=" * 78)
print("A  L2 亚层：seed、三个比、频率槽")
print("=" * 78)
SEED = [8.0, 4.0, 1.5, 1.0]
S = sum(SEED)
wn = [x / S for x in SEED]
print(f"  seed = {SEED}   和 = {S}   归一 = {[round(x, 5) for x in wn]}")
raw = [SEED[i] / SEED[i + 1] for i in range(3)]
r = sorted(raw)
print(f"  相邻比（按 seed 次序）= {[round(x, 6) for x in raw]}")
print(f"  相邻比（升序）        = {[round(x, 6) for x in r]}")
print(f"  频率槽 = ln r         = {[round(math.log(x), 6) for x in r]}")
check("三个比升序后成等比、公比 4/3",
      abs(r[1] / r[0] - 4 / 3) < 1e-12 and abs(r[2] / r[1] - 4 / 3) < 1e-12,
      f"公比 {r[1]/r[0]:.6f}")
check("三个比 = (3/2, 2, 8/3)", all(abs(r[i] - v) < 1e-12 for i, v in enumerate([1.5, 2, 8 / 3])))
check("频率槽与 Z091 的 [0.4055, 0.6931, 0.9808] 一致（舍入）",
      all(abs(math.log(r[i]) - v) < 5e-5 for i, v in enumerate([0.4055, 0.6931, 0.9808])))
check("seed 的头三项是 2 的幂（8=2³、4=2²、1=2⁰），第四项 1.5 不是",
      SEED[0] == 2 ** 3 and SEED[1] == 2 ** 2 and SEED[3] == 2 ** 0 and SEED[2] == 1.5,
      "8, 4, 1.5, 1")

# ───────────────── B  库内那条「3 代质量」是否独立 ─────────────────
print()
print("=" * 78)
print("B  库内的『3 代质量（归一）』是不是频率槽本身")
print("=" * 78)
Z091 = [0.41339, 0.70670, 1.0]            # Z091 正文 / k4_nuR_features.py 硬编码
fslots = [math.log(x) for x in r]          # 升序频率槽
fmax = fslots[-1]
fs_norm = [x / fmax for x in fslots]
print(f"  Z091 的数        = {Z091}")
print(f"  频率槽 / 最大者   = {[round(x, 5) for x in fs_norm]}")
print(f"  逐项相对差        = {[round(100*(Z091[i]/fs_norm[i]-1), 4) for i in range(3)]} (%)")
check("Z091 的『3 代质量』逐项 = 频率槽把最大者归一到 1（差 <0.02%）",
      all(abs(Z091[i] / fs_norm[i] - 1) < 2e-4 for i in range(3)),
      "同一个因子")
print("  ⟹ 这不是独立的『质量↔频率』比较：它就是频率槽本身换了个标。")

# 该路线的数值后果（k4_nuR_features.py 的做法：Δm²_atm 定标后预测 Δm²_solar）
dm2atm, dm2sol_obs = 2.45e-3, 7.39e-5
A2 = dm2atm / (Z091[2] ** 2 - Z091[1] ** 2)
A = math.sqrt(A2)
m = [A * x for x in Z091]
pred = m[1] ** 2 - m[0] ** 2
print(f"  该路线：A = {A:.5f} eV，m = ({m[0]:.4f}, {m[1]:.4f}, {m[2]:.4f}) eV")
print(f"  Δm²_solar 预测 = {pred:.4e} eV²，观测 {dm2sol_obs:.3e} ⟹ 差 {pred/dm2sol_obs:.1f} 倍")
check("该路线在中微子上差 21 倍量级（库内记为 21.8 倍）", 15 < pred / dm2sol_obs < 30,
      f"{pred/dm2sol_obs:.1f} 倍")

# ───────────────── C  四个扇区振幅与三个亚层比 ─────────────────
print()
print("=" * 78)
print("C  四个扇区振幅 vs 三个亚层比（(3/2, 2, 8/3)，公比 4/3）")
print("=" * 78)
AMP = {"轻子": math.sqrt(2), "中微子": math.sqrt(1.5),
       "上型": 1.7778326303083105, "下型": 2.682849902237757}
cands = [
    ("A_ν²", AMP["中微子"] ** 2, r[0], "r₁ = 3/2"),
    ("A_ℓ²", AMP["轻子"] ** 2, r[1], "r₂ = 2"),
    ("A_down", AMP["下型"], r[2], "r₃ = 8/3"),
    ("A_up", AMP["上型"], r[2] / r[0], "r₃/r₁ = 16/9"),
]
print(f"  {'量':<8}{'值':>12}{'亚层比组合':>16}{'值':>12}{'偏差':>10}")
rows = []
for nm, v, w, tag in cands:
    dev = 100 * (v / w - 1)
    rows.append(dict(name=nm, value=v, combo=tag, combo_value=w, dev_pct=dev))
    print(f"  {nm:<8}{v:>12.6f}{tag:>16}{w:>12.6f}{dev:>9.3f}%")
check("两个精确命中（A_ν²=r₁、A_ℓ²=r₂，偏差 <1e-12）",
      abs(rows[0]["dev_pct"]) < 1e-9 and abs(rows[1]["dev_pct"]) < 1e-9)
check("下型 = 8/3 落在 1% 内", abs(rows[2]["dev_pct"]) < 1.0, f"{rows[2]['dev_pct']:.3f}%")
check("上型 = r₃/r₁ = 16/9 落在 0.01% 内", abs(rows[3]["dev_pct"]) < 1e-2, f"{rows[3]['dev_pct']:.4f}%")
print("  注意：3 个比 ↔ 4 个振幅，映射不唯一（上型需用『比之比』）⟹ 这是候选，不是推导。")

# ───────────────── D  「三个不同大小的圆」的判别 ─────────────────
print()
print("=" * 78)
print("D  三代是一个圆还是三个圆（大小不同）")
print("=" * 78)
Al = math.sqrt(2)
A2_MeV = 246.2205 / 784.5 * 1e3
PDG = {"e": 0.51099895, "mu": 105.6584, "tau": 1776.86}
PH = {"e": 2 * math.pi / 3 + 2 / 9, "mu": 2 * math.pi / 3 + 2 / 9 + 2 * math.pi / 3,
      "tau": 2 / 9}
Rk = {}
print("  每一代反解出的单元 R_k = m_k / M_k^2（三代共用一个圆 ⟹ R_k 相同）")
for nm in ("e", "mu", "tau"):
    M = 1 + Al * math.cos(PH[nm])
    Rk[nm] = PDG[nm] / M ** 2
    print(f"    {nm:<3} θ={math.degrees(PH[nm]) % 360:9.4f}°  M={M:+.6f}  R_k={Rk[nm]:11.4f} MeV")
disp = max(Rk.values()) / min(Rk.values())
check("三代反解出的单元相同（离散度 <0.1%）——即三代共用同一个圆",
      disp - 1 < 1e-3, f"max/min−1 = {100*(disp-1):.4f}%")
print("  ⟹ 『三个亚层＝三个大小不同的圆（各有自己的振幅）』在 0.01% 水平上不成立。")

print()
print("  『大小不同的圆』在库里的位置是**扇区**：")
SEC = {"轻子": (math.sqrt(2), 1.0), "中微子": (math.sqrt(1.5), 1.0),
       "上型": (1.7778326303083105, 71.6020533976),
       "下型": (2.682849902237757, 0.9879791572463)}
for nm, (amp, u) in SEC.items():
    print(f"    {nm:<4} A_s={amp:.6f}  u_s={u:.5f}  有效振幅 A·√u = {amp*math.sqrt(u):.4f}")
check("上型扇区标度 ≈72（不同大小的圆在扇区层）",
      abs(71.6020533976 / 72 - 1) < 0.01, f"{100*(71.6020533976/72-1):+.2f}%")

res = {"seed": SEED, "sum": S, "normalized": wn, "ratios_sorted": r,
       "frequencies": [math.log(x) for x in r], "z091": Z091, "z091_ratio_checks": fs_norm,
       "neutrino_route_pred": pred, "amplitude_table": rows, "generation_units_MeV": Rk, "generation_unit_spread": disp,
       "pass": PASS, "fail": FAIL}
with open(os.path.join(OUT, "l2_sublayer_mass_check.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print()
print("=" * 78)
print(f"汇总：PASS {len(PASS)} / FAIL {len(FAIL)}")
if FAIL:
    print("  FAIL:", FAIL)
print("  → results/l2_sublayer_mass_check.json")
sys.exit(0 if not FAIL else 1)
