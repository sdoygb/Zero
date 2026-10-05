#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R23_check.py -- 核验 DIM-DESC 条件选维模型。

对应 R23_dimension_descendant_selection.md。

  F1  文档与状态边界
  F2  乘积扇区 no-go 与复制规则不选维
  F3  F_D=B*D*q^D 的四维唯一窗口
  F4  q=exp(-1/L) 时唯一峰 D=L
  F5  数值探针结果可独立复算
  F6  与 R3 / R0 / STATUS / INDEX 的当前状态一致
  F7  反过度主张：不把条件模型写成 Zero 无条件推出四维 GR
"""

from __future__ import annotations

import io
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


def head(title: str) -> None:
    print("")
    print("=" * 72)
    print(title)
    print("=" * 72)


def check(name: str, cond: bool, detail: str = "") -> None:
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def read(name: str) -> str:
    with io.open(os.path.join(HERE, name), encoding="utf-8") as handle:
        return handle.read()


def weights(D: int, q: float) -> float:
    return D * q**D


def argmax(values: dict[int, float], rel_tol: float = 1.0e-12) -> tuple[list[int], float]:
    maximum = max(values.values())
    winners = [
        D
        for D, value in values.items()
        if math.isclose(value, maximum, rel_tol=rel_tol, abs_tol=0.0)
    ]
    return winners, maximum


DIM = list(range(1, 25))
PY = read("R23_dim_desc_probe.py")
R23 = read("R23_dimension_descendant_selection.md")
R3 = read("R3_dimension_selection.md")
R0 = read("R0_publication_theorem.md")
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")

# ======================================================================
head("F1  文档与状态边界")

check("标题与对象在位", "后代优势选维" in R23 and "DIM-DESC" in R23)
check("不声称无条件从 Zero 推出四维",
      "不声称已从 Zero 无条件推出四维时空" in R23)
check("三个具名输入写清",
      all(name in R23 for name in ["DIM-SECTOR", "DIM-COST", "L=4"]))
check("结论区分已证、条件证成、开放、排除",
      all(word in R23 for word in ["已证", "条件证成", "开放", "排除"]))
check("明确 D=4 不等于四维 GR",
      "有效维数扇区标签" in R23 and "四维 GR" in R23)

# ======================================================================
head("F2  乘积扇区 no-go 与复制不选维")

for a in (0.8, 1.0, 1.2):
    ratios = [((a ** (D + 1)) / (a**D)) for D in DIM[:-1]]
    check("乘积扇区 a=%.1f 的相邻增长率相同" % a,
          max(ratios) - min(ratios) < 1.0e-14,
          "ratio=%.12g" % ratios[0])

copy_values = {D: 2.0 for D in DIM}
without_4 = {D: value for D, value in copy_values.items() if D != 4}
check("copy 规则在所有维数给相同 F_D=2",
      all(abs(value - 2.0) < 1.0e-15 for value in copy_values.values()))
check("去掉 D=4 不改变 copy 的权重集合，故不选四维",
      set(without_4.values()) == {2.0})
check("文档明写只靠复制不能选出四维",
      "仅靠“复制规则存在”不能选出四维" in R23)
check("探针代码实现乘积与 descendant 两组权重",
      "def product_weight" in PY and "def descendant_weight" in PY)

# ======================================================================
head("F3  F_D=B*D*q^D 的四维唯一窗口")

q_low = 0.7
w_low, _ = argmax({D: weights(D, q_low) for D in DIM})
check("q=0.70 时峰在 D=3", w_low == [3], str(w_low))

w_boundary_low, _ = argmax({D: weights(D, 0.75) for D in DIM})
check("q=3/4 时 D=3 与 D=4 并列",
      w_boundary_low == [3, 4], str(w_boundary_low))

w_inside, _ = argmax({D: weights(D, math.exp(-0.25)) for D in DIM})
check("q=exp(-1/4) 时唯一峰在 D=4", w_inside == [4], str(w_inside))

w_boundary_high, _ = argmax({D: weights(D, 0.8) for D in DIM})
check("q=4/5 时 D=4 与 D=5 并列",
      w_boundary_high == [4, 5], str(w_boundary_high))

w_high, _ = argmax({D: weights(D, 0.82) for D in DIM})
check("q=0.82 时峰移到 D=5", w_high == [5], str(w_high))

four_regions = []
previous = None
start = None
for index in range(40001):
    q = 0.50 + 0.40 * index / 40000
    winners, _ = argmax({D: weights(D, q) for D in DIM})
    inside = winners == [4]
    if inside and start is None:
        start = q
    if not inside and start is not None:
        four_regions.append((start, previous))
        start = None
    previous = q
if start is not None:
    four_regions.append((start, previous))

check("扫描只发现一个四维区间，且不越过 3/4 与 4/5",
      len(four_regions) == 1
      and four_regions[0][0] > 0.75 - 2.0e-5
      and four_regions[0][1] < 0.80 + 2.0e-5,
      str(four_regions))
check("文档写出四维窗口 3/4<q<4/5",
      "\\frac34<q<\\frac45" in R23 and "(0.75,0.80)" in R23)

# ======================================================================
head("F4  q=exp(-1/L) 时唯一峰 D=L")

ok = True
details = []
for L in range(1, 17):
    winners, _ = argmax({D: D * math.exp(-D / L) for D in DIM})
    ok &= winners == [L]
    details.append("%d:%s" % (L, winners))
check("L=1..16 全部唯一峰 D=L", ok, "; ".join(details[:6]) + " ...")

ratio_identity_ok = True
for L in range(2, 17):
    rising = (L / (L - 1)) * math.exp(-1.0 / L) > 1.0
    falling = (1.0 + 1.0 / L) * math.exp(-1.0 / L) < 1.0
    ratio_identity_ok &= rising and falling
check("定理 R23.4 的上升／下降比值成立", ratio_identity_ok)
check("q=exp(-1/L) 与寿命定理在位",
      "q_L=e^{-1/L}" in R23 and "定理 R23.4" in R23)

even_lifetimes = list(range(4, 25, 2))
per_cycle_peaks = [L * math.exp(-1.0) for L in even_lifetimes]
per_step_rates = [(math.log(L) - 1.0) / L for L in even_lifetimes]
per_step_winner = even_lifetimes[
    max(range(len(even_lifetimes)), key=lambda index: per_step_rates[index])
]
check("自由寿命的单周期峰值 L/e 严格递增，无有限极大点",
      all(
          per_cycle_peaks[index] < per_cycle_peaks[index + 1]
          for index in range(len(per_cycle_peaks) - 1)
      ))
check("自由寿命的每步增长率在偶数闭圈中唯一选中 L=8",
      per_step_winner == 8
      and max(per_step_rates) == per_step_rates[even_lifetimes.index(8)])
check("自由寿命 no-go 与 R23.6 在位",
      "命题 R23.6（自由寿命的归一化 no-go）" in R23
      and "L_*=e^2\\approx7.389" in R23
      and "\\boxed{L=8}" in R23
      and "“最容易留下后代”并不自动选出四维" in R23)

# ======================================================================
head("F4.5  GR 兼容扇区的生存优势证明")


def product_survival(D: int, a: float) -> float:
    return a**D

GR_DIM = [D for D in DIM if D >= 4]
for a in (0.50, 0.90, 0.99):
    surviving = {D: product_survival(D, a) for D in GR_DIM}
    winners, _ = argmax(surviving)
    check("a=%.2f 时 GR 兼容扇区生存峰唯一在 D=4" % a,
          winners == [4])
check("四维优势来自 a<1 的乘积存活，而非四维生存奖励",
      "命题 R23.7（GR 兼容扇区的条件存活率排序）" in R23
      and "S_D^{\\rm evo}(L)" in R23
      and "a_L^D" in R23
      and "0<a_L<1" in R23
      and "\\binom{L}{L/2}" in R23
      and "不再预设“四维生存率更高”" in R23
      and "GR-LB" in R23
      and "GR 兼容的演化层扇区" in R23
      and "条件存活率" in R23
      and "只比较演化层" in R23
      and "不声明整个宇宙或其它层" in R23
      and "EVO-NORM" in R23
      and "绝对存活分支数" in R23
      and "绝对后代数或长期演化层占比" in R23
      and "开放：}DIM\\text{-}SECTOR\\text{ 的 Zero 原生构造" in R23
      and "PROD\\text{-}SURV\\text{、}DIM\\text{-}COST\\text{ 的 Zero 构造" not in R23)
check("绝对存活数与条件存活率没有被混为一谈",
      math.comb(4, 2) ** 4 == 1296
      and math.comb(4, 2) ** 5 == 7776
      and (3 / 8) ** 4 > (3 / 8) ** 5
      and "S_4^{\\rm evo}" in R23
      and "S_5^{\\rm evo}" in R23)


def survival_thresholds(q: float) -> tuple[float, int]:
    base = {D: weights(D, q) for D in DIM}
    ratios = {D: base[D] / base[4] for D in DIM if D != 4}
    threshold_dimension = max(ratios, key=ratios.get)
    return ratios[threshold_dimension], threshold_dimension


def survival_winners(q: float, s_four: float) -> list[int]:
    values = {
        D: weights(D, q) * (s_four if D == 4 else 1.0)
        for D in DIM
    }
    winners, _ = argmax(values)
    return winners

threshold_070, threshold_dim_070 = survival_thresholds(0.70)
threshold_082, threshold_dim_082 = survival_thresholds(0.82)
check("q=0.70 的生存阈值由 D=3 给出且约 1.071",
      threshold_dim_070 == 3 and abs(threshold_070 - 3.0 / 2.8) < 1.0e-12)
check("q=0.82 的生存阈值由 D=5 给出且等于 1.025",
      threshold_dim_082 == 5 and abs(threshold_082 - 1.025) < 1.0e-12)
check("仅四维生存率略高仍输给五维",
      survival_winners(0.82, 1.01) == [5])
check("四维生存优势超过阈值后唯一胜出",
      survival_winners(0.82, 1.03) == [4])
check("任意外部生存乘子的阈值判据在位",
      "命题 R23.8（外部生存乘子的阈值判据）" in R23
      and "\\frac{s_4}{s_D}" in R23
      and "这不是本轮的主证明" in R23)

# ======================================================================
head("F5  数值探针实现可复算")

check("探针给出 q 扫描与 lifetime 扫描",
      "def q_scan" in PY
      and "def lifetime_case" in PY
      and "def lifetime_family_case" in PY
      and "def gr_survival_case" in PY
      and "def product_survival_case" in PY)
check("探针状态明写 DIM-SECTOR/DIM-COST 未从 Z0 导出",
      "not derived from Z0" in PY)
check("探针把 EVO-NORM 保持为开放归一化桥",
      "EVO-NORM is open" in PY and '"evo_norm_open": True' in PY)
check("探针输出无无条件 Zero→4 的字段",
      "no_unconditional_zero_to_four" in PY)
check("探针默认不写盘；只有 --json 才输出",
      "parser.add_argument(\"--json\"" in PY)

# ======================================================================
head("F6  与 R3 / R0 / STATUS / INDEX 一致")

check("R3 登记 DIM-DESC 候选",
      "DIM-DESC" in R3
      and "R23" in R3
      and "PROD-SURV" in R3
      and "GR-LB" in R3)
check("R0 把 R23 纳入 O3 条件候选",
      "R23" in R0
      and "DIM-DESC" in R0
      and "PROD-SURV" in R0
      and "GR-LB" in R0)
check("STATUS 登记 R23 且保持 L1/E2 边界",
      "### 2.23" in STATUS
      and "R23" in STATUS
      and "DIM-DESC" in STATUS
      and "自由寿命归一化" in STATUS
      and "PROD-SURV" in STATUS
      and "GR-LB" in STATUS)
check("INDEX 已由生成器纳入 R23", "R23_dimension_descendant_selection.md" in INDEX)
check("R23 明确 GR 只是临时脚手架并指向 R24",
      "R24_global_four_survival_gate.md" in R23
      and "SURV4-GR-SCAFFOLD" in R23
      and "DIM-COST-Q" in R23
      and "临时" in R23)
check("R23 已同步 R25 单方向 no-go 与成对窗口",
      "R25_native_pair_cost_and_four_dim_peak.md" in R23
      and "LEDGER-ROT" in R23
      and "PAIR-CARRIER" in R23
      and "1/2<q<3/5" in R23
      and "q=5/9" in R23)
check("R23 已同步 R26 的成对载体缺口分解",
      "R26_pair_carrier_reduction_no_go.md" in R23
      and "PAIR-CARRIER-DER" in R23
      and "D=m-1" in R23
      and "D=3" in R23)

# ======================================================================
head("F7  反过度主张")

for phrase in [
    "不能写成“Zero 无条件推出四维”",
    "仍然需要的是",
    "当前还缺什么",
    "即使选出 `D=4`，也还没有四维 GR",
]:
    check("反过度主张句在位：%s" % phrase, phrase in R23)

check("文档保留 E1–E4 / Lorentz / Lovelock 分离条件",
      all(mark in R23 for mark in ["Lorentz", "反射正真空", "Lovelock", "R1", "R2"]))
check("没有把 R23 写成 O3 已关闭",
      "本轮没有关闭 O3" in R23 and "L=4` 仍只是条件" in R23)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
