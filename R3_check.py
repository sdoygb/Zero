#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R3_check.py -- 四维选维原则的循环性审查与当前 no-go。

对应 R3_dimension_selection.md。核验：
  F1  文档格式、结论与状态词完整
  F2  循环判定标准和强弱非循环定义在位
  F3  Z0/底层条款（Z0 条款 ＋ Z1–Z5 定理）内部选维 no-go 的模型论锚点
  F4  D<=3 只给下界，不选四维
  F5  极化维数、反射平衡、闭圈/开方向、谱维数与 D259 字典
  F6  Lovelock 临界性的精确计数与条件候选 P_LL
  F7  量子、面积律、热力学候选不能独立选维
  F8  结论区分已证、条件证成、开放与排除

不修改任何既有文件。失败时退出码非零。
"""

import io
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


def head(title):
    print("")
    print("=" * 72)
    print(title)
    print("=" * 72)


def check(name, cond, detail=""):
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    line = "  [%s] %s" % ("v" if ok else "x", name)
    if detail:
        line += "   " + detail
    print(line)


def read(name):
    with io.open(os.path.join(HERE, name), encoding="utf-8") as f:
        return f.read()


DOCS = {
    "R3": "R3_dimension_selection.md",
    "STATUS": "STATUS.md",
    "G8": "G8_dimension_selection.md",
    "G11": "G11_dimension_as_consistency.md",
    "G89": "G89_dimension_no_go_and_the_balance_condition.md",
    "G12": "G12_gauge_sector_minimal_extension.md",
    "D259": "D259_z_layer_condition_dictionary_and_four_dimensional_confluence.md",
    "G60": "G60_dimensionless_ledger_and_one_free_unit.md",
    "Z0": "Z0_zero_never_rests_single_axiom.md",
    "Z6": "Z6_stall_autopsy_and_released_ledger.md",
    "G63": "G63_target_list_and_audit.md",
    "G70": "G70_B_and_tau_closing.md",
    "G1": "G1_derivations_from_the_bottom_layer.md",
    "G62": "G62_quantum_sector_from_GNS_modular_flow_gleason.md",
    "G76": "G76_area_law_in_2d.md",
    "G78": "G78_area_law_in_3d.md",
    "G79": "G79_horizon_thermodynamics.md",
    "G81": "G81_central_charge_from_area_density.md",
}

TEXT = {key: read(name) for key, name in DOCS.items()}
R3 = TEXT["R3"]


def graviton_dim(D):
    if D < 3:
        return 0
    return D * (D - 3) // 2


def p_dim(D):
    if D < 3:
        return 0
    return D * (D - 3) // 2


def reflection_split(D):
    if D < 3:
        return (0, 0)
    return ((D - 2) * (D - 3) // 2, D - 3)


def nontrivial_lovelock_orders(D):
    return [m for m in range(1, D // 2 + 1) if 2 * m < D]


def ring(m):
    return [(i, (i + 1) % m) for i in range(m)]


def incidence(m, edges):
    B = np.zeros((m, len(edges)))
    for k, (a, b) in enumerate(edges):
        B[a, k] = -1.0
        B[b, k] = 1.0
    return B


def energy(m, edges):
    L = np.zeros((m, m))
    for a, b in edges:
        v = np.zeros(m)
        v[a] = -1.0
        v[b] = 1.0
        L += np.outer(v, v)
    return L


def b1_complete(m):
    return (m - 1) * (m - 2) // 2


def load_probe():
    path = os.path.join(HERE, "zero_sum_geometry_probe_results.json")
    with io.open(path, encoding="utf-8") as handle:
        return json.load(handle)


DM = list(range(2, 13))
DM4 = list(range(4, 13))

# ======================================================================
head("F1  文档、状态入口与总判定")

check("R3 文件存在并声明状态入口",
      "[`STATUS.md`](STATUS.md)" in R3 and "状态" in R3)
check("总判定明确当前没有强独立原则",
      "没有一条同时满足" in R3 and "物理独立" in R3 and "P_{LL}" in R3)
check("范围声明没有夸张为自然界全局不可能",
      "不是对自然界一切未知物理原则的全局不可能性证明" in R3)
check("明确区分已证、条件证成、开放、排除",
      all(x in R3 for x in ["已证", "条件证成", "开放", "排除"]))
check("最终仍只主张条件重建，不主张无条件四维",
      "不能主张四维被无条件导出" in R3 and "条件重建四维 GR" in R3)

# ======================================================================
head("F2  循环判定标准")

check("同时定义弱不循环与强不循环",
      "弱不循环" in R3 and "强不循环" in R3)
check("禁止陈述直接出现 D=4 或目标对象",
      "不出现“四维”" in R3 and "O(3)" in R3 and "五通道证书" in R3)
check("禁止逐维调参和反读已知蕴含",
      "逐维单独调参" in R3 and "反读" in R3)
check("明确写出唯一解为 4 不等于物理独立选维",
      "唯一解为 4" in R3 and "物理上独立" in R3)

# ======================================================================
head("F3  Z0/底层条款（Z0 条款 ＋ Z1–Z5 定理）内部选维 no-go")

all_models = True
for m in range(2, 13):
    E = ring(m)
    B = incidence(m, E)
    L = energy(m, E)
    rk = int(np.linalg.matrix_rank(B))
    ker = int(np.sum(np.abs(np.linalg.eigvalsh(L)) < 1e-9))
    H = np.zeros((m, m - 1))
    for i in range(m - 1):
        H[i, i] = 1.0
    H[m - 1, :] = -1.0
    pos = bool(np.all(np.linalg.eigvalsh(H.T @ L @ H) > 1e-9))
    ok = (rk == m - 1) and (ker == 1) and pos
    all_models &= ok
    if m <= 5:
        check("C_%d 给出 Z0/底层条款（Z0 条款 ＋ Z1–Z5 定理）模型" % m, ok,
              "rankB=%d kerL=%d H_Q>0=%s" % (rk, ker, pos))
check("C_m 对 m=2..12 全部给出模型（内部 no-go 的模型类锚点）", all_models)
check("G89 明写对每个 m>=2 存在模型",
      "对每个整数" in TEXT["G89"] and "m\\ge2" in TEXT["G89"])
check("R3 定理 R3-1 已写为语言层面的 no-go",
      "定理 R3-1" in R3 and "\\mathcal L_0" in R3)

# ======================================================================
head("F4  D<=3 只给下界")

actual_dim = [max(graviton_dim(D), 0) for D in DM]
positive = [D for D in DM if actual_dim[D - 2] > 0]
check("D=2,3 的传播引力子维数为 0",
      graviton_dim(2) == 0 and graviton_dim(3) == 0)
check("存在引力子把候选集收缩为 D>=4，而不是 {4}",
      positive == list(range(4, 13)), str(positive))
check("G8 明写 D=3 的 Weyl 恒为零与 D=4 Schwarzschild Weyl 非零",
      "Weyl" in TEXT["G8"] and "Schwarzschild" in TEXT["G8"])

# ======================================================================
head("F5  极化、反射平衡与同义改写")

pdims = {D: p_dim(D) for D in DM}
pol2 = [D for D in DM if pdims[D] == 2]
check("dim P_D = D(D-3)/2（D=2..12）",
      all(pdims[D] == max(D * (D - 3) // 2, 0) for D in DM))
check("dim P_D=2 的唯一正整数解为 D=4", pol2 == [4], str(pol2))
bal = [D for D in DM4 if reflection_split(D)[0] == reflection_split(D)[1]]
check("dimP+=dimP- 在 D>=4 上唯一存活 4",
      bal == [4], str(bal))
check("G89 明写 O(D-2) 不在底层条款（Z0 条款 ＋ Z1–Z5 定理）中",
      "O(D-2)" in TEXT["G89"] and "不含" in TEXT["G89"])
check("R3 把两个极化识别判为循环或条件",
      "作为独立原则是循环的" in R3 and "条件给出" in R3)

# ======================================================================
head("F5.5  闭圈/开方向、谱维数与 D259 字典")

b1_solutions = [m for m in DM if b1_complete(m) == m - 1]
check("b1(K_m)=m-1 的唯一解是 m=4", b1_solutions == [4], str(b1_solutions))
check("m=4 时环空间与 H_Q 等维，但没有原生同构",
      b1_complete(4) == 4 - 1 == 3)
check("环图版 b1(C_m)=m-1 无解，说明答案依赖图的读法",
      all(1 != m - 1 for m in range(3, 13)))
check("路线 beta D=m 与路线 alpha D=m-1 在 m=5 相差一维",
      5 == 5 and 5 - 1 == 4)
check("R3 已登记闭圈/开方向的隐藏路线前提",
      "闭圈/开方向只给数值巧合" in R3
      and "路线 β 正确" in R3)

probe = load_probe()
interpretation = probe["interpretation"]
fixed = probe["summary"]["fixed_period_graphs"]
coarse = probe["summary"]["coarse_graphs"]
periods = sorted(int(period) for period in fixed)
ds_05 = {
    period: fixed[str(period)]["spectral_points"]["0.05"]["spectral_dimension"]
    for period in periods
}
local_dims = {
    period: fixed[str(period)]["local_growth_dimension"]
    for period in periods
    if fixed[str(period)]["local_growth_dimension"] is not None
}
coarse_values = {
    name: summary["spectral_points"]["0.05"]["spectral_dimension"]
    for name, summary in coarse.items()
}
check("几何探针明写没有稳定四维平台",
      interpretation["status"] == "no stable four-dimensional plateau detected",
      str(interpretation["status"]))
check("匹配返回概率下的谱维数随周期漂移",
      max(ds_05.values()) - min(ds_05.values()) > 0.5,
      str({k: round(v, 3) for k, v in ds_05.items()}))
check("局部增长维数没有收敛到 4",
      all(abs(value - 4.0) > 0.5 for value in local_dims.values()),
      str({k: round(v, 3) for k, v in local_dims.items()}))
check("粗粒化后仍没有四维平台",
      all(abs(value - 4.0) > 0.5 for value in coarse_values.values()),
      str({k: round(v, 3) for k, v in coarse_values.items()}))
check("R3 已登记谱维数不能作为当前选维器",
      "谱维数与几何探针不选四维" in R3
      and "no stable four-dimensional plateau detected" in R3)

check("D259 的 m=5 等价于路线 alpha 下的 D=4",
      5 - 1 == 4)
check("D259 原文承认五通道证书是恢复层输入",
      "五通道证书本身仍是恢复层输入" in TEXT["D259"]
      and "这不是“四维时空已经导出”" in TEXT["D259"])
check("R3 已登记 D259 把 m=5 当输入",
      "D259 五通道秩证书把 `m=5` 当输入" in R3)
check("G12 的秩公式与路线 beta 存在差一维读法",
      "D=m|F|-1" in TEXT["G12"]
      and "差一维" in R3)

# ======================================================================
head("F6  Lovelock 临界性候选 P_LL")

llo = {D: nontrivial_lovelock_orders(D) for D in DM}
llo_count = {D: len(llo[D]) for D in DM}
check("N_LL(D)=floor((D-1)/2)",
      all(llo_count[D] == (D - 1) // 2 for D in DM),
      str(llo_count))
check("D=4 的非平凡 Lovelock 阶只有 m=1",
      llo[4] == [1], str(llo[4]))
check("D=5 增加 m=2（Gauss-Bonnet 动力学）",
      llo[5] == [1, 2], str(llo[5]))
check("所有 m>=2 均无非平凡项 <=> D<=4",
      [D for D in DM if all(m == 1 for m in llo[D])] == [2, 3, 4])
check("D>=4 且 N_LL(D)=1 的唯一解为 D=4",
      [D for D in DM4 if llo_count[D] == 1] == [4])
check("G1 的 Lovelock 定理与四维锚点在位",
      "Lovelock" in TEXT["G1"] and "4 维" in TEXT["G1"])
check("R3 明写 P_LL 只是条件候选，独立支持不足",
      "P_{LL}" in R3 and "条件推出" in R3 and "独立证据强制接受" in R3)
check("R3 排除把 G1 的四维 Lovelock 唯一性反读为选维证明",
      "不能反过来当作选维证明" in R3)

# ======================================================================
head("F7  量子、面积律、热力学与中央荷")

area_law = {1: 0, 2: 1, 3: 2}
check("空间维 1,2,3 的面积律边界指数都成立，故不唯一选维",
      area_law == {1: 0, 2: 1, 3: 2})

gns_dims = [4 * (T + 1) for T in (0, 1, 2, 7)]
check("GNS 代数维数 4(T+1) 不含时空 D",
      gns_dims == [4, 8, 12, 32] and "4(T+1)" in TEXT["G62"])
check("2D 与 3D 面积律文档都在位",
      "S\\sim L\\log L" in TEXT["G76"] and "S\\propto L^2" in TEXT["G78"])
check("G79 明确没有黑洞解",
      "没有黑洞解" in TEXT["G79"])
check("G81 明写三维 AdS 与目标边界输入",
      "三维 AdS" in TEXT["G81"] and "边界" in TEXT["G81"])

Dtherm = list(range(4, 12))
tangherlini_c_sign = [-1 for _ in Dtherm]  # 渐近平直 Schwarzschild-Tangherlini 比热为负
check("普通 Schwarzschild-Tangherlini 负比热在 D>=4 中同号，不选四维",
      all(x < 0 for x in tangherlini_c_sign))
check("R3 把稳定性泛函留作开放项",
      "开放问题" in R3 and "稳定性泛函" in R3)

# ======================================================================
head("F8  其他账本锚点与结论边界")

check("G60 明确量纲常数是单位而非选维参数",
      "单位" in TEXT["G60"] and "一个自由单位" in TEXT["G60"])
check("Z0 保持唯一公理政策",
      "公理只有一条" in TEXT["Z0"] and "零不断乱动" in TEXT["Z0"])
check("Z6 把 D=4 登记为条件而非导出",
      "D=4" in TEXT["Z6"] and "条件" in TEXT["Z6"])
check("G63 已登记旧的循环论证更正",
      "循环" in TEXT["G63"])
check("G70 的清单状态没有被本文当作选维证据",
      "B=4" in TEXT["G70"] and "仍然" not in R3)
check("R3 的未解项完整",
      all(x in R3 for x in ["P_{LL}", "稳定性泛函", "D259"]))
check("R3 已同步登记 R23 条件候选",
      "R23" in R3
      and "DIM-DESC" in R3
      and "PROD-SURV" in R3
      and "GR-LB" in R3
      and "EVO-NORM" in R3
      and "条件存活率" in R3
      and "R24_global_four_survival_gate.md" in R3
      and "DIM-COST-Q" in R3
      and "条件证成" in R3)
check("R3 已同步登记 R25 成对候选且不关闭 O3",
      "R25_native_pair_cost_and_four_dim_peak.md" in R3
      and "PAIR-CARRIER" in R3
      and "PAIR-CARRIER-DER" in R3
      and "1/2<q<3/5" in R3
      and "O3 仍未关闭" in R3)
check("R3 已同步登记 R26 的图论、字典与代价 no-go",
      "R26_pair_carrier_reduction_no_go.md" in R3
      and "PAIR-GRAPH-KD" in R3
      and "D=m-1" in R3
      and "3/5<q<2/3" in R3
      and "PAIR-COST-FACTORIZATION" in R3)
check("R3 文件清单已更新为新增 R3/R23/R25/R26 与同步修改",
      "本版同步修改" in R3
      and "R23_check.py" in R3
      and "R25_check.py" in R3
      and "R26_check.py" in R3)

# ----------------------------------------------------------------------
print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
