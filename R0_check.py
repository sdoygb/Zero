#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R0_check.py -- 核验可发表主定理接口与五项推进门槛。

对应文档 R0_publication_theorem.md。

  F1  主定理接口、目标对象与结论在位
  F2  O1-O5 五条交付完整
  F3  M0-M4 门槛在位，且当前状态未被夸大为完成
  F4  反过度主张规则在位
  F5  STATUS.md 仍是唯一当前状态源
  F6  R1-R29 存在时，每项都声明自身状态与证据
  F7  R8-R22 只作外部基准，R23-R29 只作 O3 条件候选，不抬高主状态
"""

import io
import os
import sys


HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def read(name):
    with io.open(os.path.join(HERE, name), encoding="utf-8") as f:
        return f.read()


R0 = read("R0_publication_theorem.md")
STATUS = read("STATUS.md")

head("F1  主定理接口")

check("只主张条件重建，不主张无条件导出",
      "条件重建定理" in R0 and "不能写成" in R0 and "无条件导出四维广义相对论" in R0)
check("离散输入对象完整",
      all(x in R0 for x in ["\\Gamma_N", "\\mathcal W_N", "\\iota_N", "k_N=L_N", "d=4"]))
check("目标连续对象完整",
      all(x in R0 for x in ["(M,g,\\varphi,T_{ab},\\mu_*)", "Lorentz 度规", "守恒应力张量", "pulled-front"]))
check("主结论包含度规、Lovelock、源与前沿",
      all(x in R0 for x in ["Lovelock 前提", "G_{ab}+\\Lambda g_{ab}=8\\pi G\\,T_{ab}", "源守恒", "有限有效前沿"]))

head("F2  O1-O5 五项交付")

for tag, name, markers in [
    ("O1", "Γ-收敛", ["R1_gamma_convergence_theorem.md", "紧性", "Γ-liminf", "恢复列"]),
    ("O2", "类到站点", ["R2_site_identification.md", "良定义", "与寿命", "准均匀"]),
    ("O3", "独立选维", ["R3_dimension_selection.md", "不引用“四维”", "若不可得", "最小性"]),
    ("O4", "可证伪预测", ["R4_falsifiable_prediction.md", "预测式", "排除阈值", "有效锥指数尾"]),
    ("O5", "路线汇流", ["R5_route_equivalence.md", "证明关键结构保持", "明确选出唯一路线", "五通道证书"]),
]:
    check("%s %s 接口完整" % (tag, name), all(m in R0 for m in markers))

head("F3  门槛与当前状态")

check("M0-M4 门槛在位",
      all(("**%s" % x) in R0 for x in ["M0", "M1", "M2", "M3", "M4"]))
check("当前只声明 M0，未夸大为完成",
      "当前状态：**M0；正在推进 M1–M3**" in R0)

head("F4  反过度主张")

check("六条禁止句与允许句均在位",
      all(x in R0 for x in ["从零无条件导出 GR", "证明四维必然", "连续极限已证明", "有效锥就是光锥"]))

head("F5  唯一状态源")

check("R0 指向 STATUS.md", "[`STATUS.md`](STATUS.md)" in R0)
check("STATUS 仍是唯一当前状态源", "唯一状态源" in STATUS and "当前状态" in STATUS)

head("F6  R1-R29 状态登记")

for name in [
    "R1_gamma_convergence_theorem.md",
    "R2_site_identification.md",
    "R3_dimension_selection.md",
    "R4_falsifiable_prediction.md",
    "R5_route_equivalence.md",
    "R6_h5_dictionary_error_bound.md",
    "R7_h3_h7_regularity.md",
    "R8_jacobson_entanglement_equilibrium_completion.md",
    "R9_external_GR_derivations_landscape.md",
    "R10_cao_carroll_bulk_entanglement_completion.md",
    "R12_zero_native_gap_filling.md",
    "R13_L1_strong_resolvent_attempt.md",
    "R14_L1_from_zero_assembly.md",
    "R15_zcar_double_cover_and_zstress_scale.md",
    "R16_direction_audit_reduction_tree.md",
    "R17_L1_critical_path_and_L5_gate.md",
    "R18_L5_cert_dimension_normalization_gate.md",
    "R19_L1_upstream_probability_phase_and_missing_boost.md",
    "R20_A1_verdict_shape_holds_constant_fails.md",
    "R21_vf_normalization_resolves_the_r20_factor.md",
    "R22_principal_symbol_vs_r20_estimator.md",
    "R23_dimension_descendant_selection.md",
    "R24_global_four_survival_gate.md",
    "R25_native_pair_cost_and_four_dim_peak.md",
    "R26_pair_carrier_reduction_no_go.md",
    "R27_phase_cochain_quotient_and_dimension_dictionary.md",
    "R28_phase_identity_cluster_reduction.md",
    "R29_full_support_ledger_factorization_no_go.md",
]:
    p = os.path.join(HERE, name)
    if os.path.exists(p):
        text = read(name)
        check("%s 已登记状态与证据" % name,
              any(x in text for x in ["状态", "结论", "边界", "未证", "证明"]))
    else:
        print("  [i] %s 尚未落地（并行分支运行中）" % name)

head("F7  外部基准不冒充主证明义务")

check("R8-R22 明列为外部基准，不作为 O6",
      any(("`R8`–`R%d` 不属于 O1–O5" % n) in R0 for n in range(22, 40))
      and "不计入五条主证明义务" in R0)
check("R23 登记为 O3 条件候选且不关闭 O3",
      "R23_dimension_descendant_selection.md" in R0
      and "O3 仍未关闭" in R0
      and "DIM-SECTOR" in R0
      and "DIM-COST" in R0
      and "GR-LB" in R0
      and "EVO-NORM" in R0
      and "S_D^{\\rm evo}(L)" in R0
      and "只关于演化层" in R0)
check("R24 固定全局目标并声明 GR 最终撤掉",
      "R24_global_four_survival_gate.md" in R0
      and "SURV4-GLOBAL" in R0
      and "GR-LB" in R0
      and "最终要撤掉" in R0
      and "DIM-COST-Q" in R0)
check("R25 收窄 O3 为成对候选且不关闭 O3",
      "R25_native_pair_cost_and_four_dim_peak.md" in R0
      and "PAIR-CARRIER" in R0
      and "PAIR-CARRIER-DER" in R0
      and "1/2<q<3/5" in R0
      and "O3 仍未关闭" in R0)
check("R26 把成对载体拆为可审计接口且不关闭 O3",
      "R26_pair_carrier_reduction_no_go.md" in R0
      and "PAIR-GRAPH-KD" in R0
      and "D=m-1" in R0
      and "PAIR-COST-FACTORIZATION" in R0
      and "O3 仍未关闭" in R0)
check("R27 汇流相位商与 D=m-1 字典且不关闭 O3",
      "R27_phase_cochain_quotient_and_dimension_dictionary.md" in R0
      and "PHASE-1-COCHAIN" in R0
      and "PAIR-ID-QUOTIENT" in R0
      and "FULL-SUPPORT-LEDGER" in R0
      and "WIPE-RESET-LEDGER" in R0
      and "PARALLEL-PERIOD" not in R0
      and "共同毁灭代际" in R0
      and "活动层全灭" in R0
      and "历史层只保留最高两层亚层" in R0
      and "零层 `\\mathcal Z_\\ast` 全部保留" in R0
      and "EVO-NORM" in R0
      and "R27 不关闭 O3" in R0)
check("R28 把相位身份簇归约并证明两条输入不足",
      "R28_phase_identity_cluster_reduction.md" in R0
      and "PHASE-IDENTITY-DER" in R0
      and "EDGE-CONNECTION" in R0
      and "HOLONOMY-FULL-SPAN" in R0
      and "INHERITANCE-IDENTITY" in R0
      and "O3 仍关闭失败" in R0)
check("R29 把全支撑账本拆为四项并证明乘积记录缺口",
      "R29_full_support_ledger_factorization_no_go.md" in R0
      and "LEDGER-FACTORIZATION" in R0
      and "DIR-SUPPORT-D" in R0
      and "RECORD-FAMILY-D" in R0
      and "PRODUCT-LEDGER" in R0
      and "SAME-Q" in R0
      and "O3 仍未关闭" in R0)
check("R8 仍为条件恢复且 L1 是决定性缺口",
      "精确熵差恒等式" in R0 and "L1" in R0 and "条件恢复" in R0)

head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
