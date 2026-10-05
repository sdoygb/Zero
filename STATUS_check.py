#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
STATUS_check.py -- 核验全项目唯一当前状态源。

对应文档 STATUS.md。

  F1  唯一状态源存在，并包含统一状态词与 I1-I12 当前表
  F2  INDEX / G10 / G0 / Z6 / Z1 / Z2 / D259 / verify 均指向 STATUS.md
  F3  关键旧状态页已标为历史，且旧“唯一账本”措辞清零
  F4  主 Z/G 路线与借入 D259 路线分开记账
  F5  参数与目标清单只采用 G70 的当前计数
  F6  ledger_sync.py 会运行本核验
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


def has(name, *needles):
    text = read(name)
    return all(n in text for n in needles)


STATUS = read("STATUS.md")
INDEX = read("INDEX.md")

# ======================================================================
head("F1  唯一状态源与状态表")

check("STATUS.md 声明自己是唯一当前状态源",
      "唯一状态源" in STATUS and "只有本文件给出的状态可以当作当前项目状态" in STATUS)
check("统一状态词在位",
      all(x in STATUS for x in ["已证", "条件证成", "具名输入", "条件", "开放", "排除", "不可导出", "已撤回", "历史"]))
check("I1-I12 当前表在位",
      all(("**%s**" % x) in STATUS for x in
          ["I1", "I2a", "I2b", "I2c", "I3a", "I3b", "I3c", "I4", "I5",
           "I6", "I7", "I8a", "I8b", "I9a", "I9b", "I10", "I11", "I12"]))
check("I2a 写明条件定理（R1-A/R1-B）、H5 已由 R6 关闭、H3/H7 已由 R7 条件闭合",
      "条件定理（R1-A" in STATUS and "H5 已由 R6" in STATUS
      and "H3/H7 已由 R7 条件闭合" in STATUS)
check("I4 已写为不可导出定理而非待解缺口",
      "**不可导出**；作为单位／约定 E4" in STATUS and "旧“未建立”撤回" in STATUS)
check("I8 条款与几何剖面分开记账",
      "**I8a**" in STATUS and "**I8b**" in STATUS and "旧“输入 I8”与“公理 Z5”均历史化" in STATUS)
check("R8 外部基准以条件恢复登记，L1 为决定性开放项",
      "### 2.5 外部强论文补前提（`R8`）当前状态" in STATUS
      and "Jacobson 2016" in STATUS
      and "引理 R8.1" in STATUS
      and "L1 几何模极限" in STATUS
      and "未关闭时 R8 只能称条件恢复" in STATUS
      and "不能说本项目已经补上 Jacobson 2016" in STATUS)
check("R8 没有被写成主路线已关闭或无条件成功",
      "外部基准，不替代主 `Z/G` 路线" in STATUS
      and "R8 已无条件导出四维 GR" not in STATUS
      and "已补上 Jacobson 2016 了" not in STATUS)
check("R9 外部排序与 L5 已登记",
      "### 2.6 外部路线排序（`R9`）当前状态" in STATUS
      and "Cao-Carroll 2018 取代 Jacobson 1995" in STATUS
      and "L5 低维相关算符污染" in STATUS
      and "R^{2\\Delta}" in STATUS)
check("R10 Cao-Carroll 次选桥只到弱场且决定性缺口已登记",
      "Cao-Carroll 条件桥" in STATUS
      and "不推出完整非线性 GR" in STATUS
      and "RC 割函数、跨切割面积–互信息比例、Radon 反演、Lorentzian 组装" in STATUS
      and "R10_check.py" in STATUS)
check("R11 旧理论审计只作历史与外部基准",
      "### 2.7 旧理论与历史 GR 推导审计（`R11`）当前状态" in STATUS
      and "不计入 O1–O5" in STATUS
      and "旧材料不能关闭 J1 几何 boost" in STATUS
      and "旧材料不能关闭 L5 低维算符污染" in STATUS
      and "R11_check.py" in STATUS)
check("R23 后代选维只作 O3 条件候选",
      "### 2.23 后代优势选维与 `DIM-DESC` 条件模型（`R23`）当前状态" in STATUS
      and "乘积扇区 no-go" in STATUS
      and "3/4<q<4/5" in STATUS
      and "不关闭 O3" in STATUS
      and "不含 Lorentz、反射正性、boost、Lovelock" in STATUS
      and "PROD-SURV" in STATUS
      and "GR-LB" in STATUS
      and "演化层" in STATUS
      and "EVO-NORM" in STATUS
      and "条件存活率" in STATUS
      and "s_4/s_D>F_D^0/F_4^0" in STATUS
      and "R23_check.py" in STATUS)
check("R24 登记全局四维生存峰门槛且保持未证",
      "### 2.24 全局四维生存峰门槛与移除 GR（`R24`）当前状态" in STATUS
      and "TEMP-GR" in STATUS
      and "SURV4-GLOBAL" in STATUS
      and "未证" in STATUS
      and "DIM-COST-Q" in STATUS
      and "DIM-INTERACT" in STATUS
      and "R24_check.py" in STATUS)
check("R25 登记成对候选、单方向 no-go 与开放原生桥",
      "### 2.25 成对连接与旋转类读出账本四维峰（`R25`）当前状态" in STATUS
      and "LEDGER-ROT" in STATUS
      and "PAIR-CARRIER" in STATUS
      and "PAIR-CARRIER-DER" in STATUS
      and "1/2<q<3/5" in STATUS
      and "q_4=5/9" in STATUS
      and "条件证成" in STATUS
      and "R25_check.py" in STATUS)
check("R26 登记成对载体约化、字典峰移动与代价缺口",
      "### 2.26 成对载体约化与两条 no-go（`R26`）当前状态" in STATUS
      and "PAIR-GRAPH-KD" in STATUS
      and "PAIR-ID-EDGE" in STATUS
      and "PAIR-COST-FACTORIZATION" in STATUS
      and "PAIR-NO-EXTRA-MULT" in STATUS
      and "D=m-1" in STATUS
      and "q=5/9" in STATUS
      and "D=3" in STATUS
      and "R26_check.py" in STATUS)
check("R27 登记相位上链商、字典汇流与共同毁灭代际口径",
      "### 2.27 相位上链商与维数字典汇流（`R27`）当前状态" in STATUS
      and "PHASE-1-COCHAIN" in STATUS
      and "PAIR-ID-QUOTIENT" in STATUS
      and "FULL-SUPPORT-LEDGER" in STATUS
      and "WIPE-RESET-LEDGER" in STATUS
      and "PARALLEL-PERIOD" not in STATUS
      and "活动层全清" in STATUS
      and "历史层只保留最高两层亚层" in STATUS
      and "`\\mathcal Z_\\ast` 全保留" in STATUS
      and "D222 允许不同区域不同寿命" in STATUS
      and "绝对层占比" in STATUS
      and "C(D,2)" in STATUS
      and "R27_check.py" in STATUS)
check("R28 登记纯规范 no-go 与身份簇三重归约",
      "### 2.28 相位身份簇归约与纯规范 no-go（`R28`）当前状态" in STATUS
      and "PHASE-IDENTITY-DER" in STATUS
      and "EDGE-CONNECTION" in STATUS
      and "HOLONOMY-FULL-SPAN" in STATUS
      and "INHERITANCE-IDENTITY" in STATUS
      and "纯规范记录" in STATUS
      and "身份秩与代价" in STATUS
      and "R28_check.py" in STATUS)
check("R29 登记全支撑账本四项分解与乘积记录 no-go",
      "### 2.29 全支撑账本四项分解与乘积记录 no-go（`R29`）当前状态" in STATUS
      and "LEDGER-FACTORIZATION" in STATUS
      and "DIR-SUPPORT-D" in STATUS
      and "RECORD-FAMILY-D" in STATUS
      and "PRODUCT-LEDGER" in STATUS
      and "SAME-Q" in STATUS
      and "全支撑自动给乘积记录" in STATUS
      and "R29_check.py" in STATUS)
check("R30 登记小群—相位路线的对抗审计与撤回",
      "### 2.30 小群—相位路线的对抗审计（`R30`）当前状态" in STATUS
      and "LG-FUNCTOR" in STATUS
      and "空解" in STATUS
      and "同义改写" in STATUS
      and "Z4-EVIDENCE-SCOPE" in STATUS
      and "R30_little_group_phase_route_audit.md" in STATUS
      and "R30_check.py" in STATUS)
check("R31 登记相位接通账本与 L=4 唯一选择",
      "### 2.31 相位接通账本：`L=4` 的唯一选择（`R31`）当前状态" in STATUS
      and "PEAK-IN-GRAVITON-DOMAIN" in STATUS
      and "定理 R31.1" in STATUS
      and "q_L" in STATUS
      and "相位**做不了小群**，但**能做账本**" in STATUS
      and "R31_phase_ledger_and_lifetime_selection.md" in STATUS
      and "R31_check.py" in STATUS)
check("R32 登记账本读出选择与 L=8 张力消解",
      "### 2.32 账本读出的选择与 `L=8` 张力的消解（`R32`）当前状态" in STATUS
      and "定理 R32.1" in STATUS
      and "定理 R32.2" in STATUS
      and "LEDGER-ROT" in STATUS
      and "条件选择" in STATUS
      and "R32_ledger_readout_selection_and_L8_resolution.md" in STATUS
      and "R32_check.py" in STATUS)
check("R33 登记作用量相位立项",
      "### 2.33 作用量相位立项（`R33`）当前状态" in STATUS
      and "ACTION-PHASE-MATCH" in STATUS
      and "子目标 **S1**" in STATUS and "子目标 **S4**" in STATUS
      and "R33_action_phase_match_project.md" in STATUS
      and "R33_check.py" in STATUS)
check("R34 登记有限维 boost 不可能定理",
      "### 2.34 P2 探针：有限维不可能承载 boost（`R34`）当前状态" in STATUS
      and "定理 R34.1" in STATUS
      and "type III" in STATUS
      and "R34_finite_dimensional_boost_obstruction.md" in STATUS
      and "R34_check.py" in STATUS)
check("R35 登记 P2″ 类型判据与对 π 的新约束",
      "### 2.35 P2″：极限类型由粗粒化轮廓决定，并反向约束 `π`（`R35`）当前状态" in STATUS
      and "非等差" in STATUS
      and "III$_{1/2}$" in STATUS
      and "R35_type_iii_classification.md" in STATUS
      and "R35_check.py" in STATUS)
check("R36 登记 CHSH 探针与 Bell 局域判定",
      "### 2.36 CHSH 探针：现有二分割上 Bell 局域（`R36`）当前状态" in STATUS
      and "Bell 局域" in STATUS
      and "违反门槛" in STATUS
      and "R36_chsh_bell_locality.md" in STATUS
      and "R36_check.py" in STATUS)
check("R37 登记 KCBS 语境性检验",
      "### 2.37 KCBS 探针：单体统计**是语境的**（`R37`）当前状态" in STATUS
      and "语境" in STATUS
      and "0.7236" in STATUS
      and "R37_kcbs_contextuality.md" in STATUS
      and "R37_check.py" in STATUS)
check("R38 登记纠缠涌现机制",
      "### 2.38 纠缠的涌现机制：共同起因＋未记录自由度（`R38`）当前状态" in STATUS
      and "2.828427" in STATUS
      and "共同起因" in STATUS
      and "R38_entanglement_from_shared_closure_origin.md" in STATUS
      and "R38_check.py" in STATUS)
check("R39 登记历史层位点失明裁决",
      "### 2.39 历史层对位点的失明裁决（`R39`）当前状态" in STATUS
      and "失明" in STATUS
      and "播种点" in STATUS
      and "R39_history_layer_site_blindness.md" in STATUS
      and "R39_check.py" in STATUS)
check("R40 登记播种不注入位点信息",
      "### 2.40 播种不注入位点信息（`R40`）当前状态" in STATUS
      and "失明是结构性的" in STATUS
      and "R40_seeding_does_not_inject_site.md" in STATUS
      and "R40_check.py" in STATUS)
check("R41 登记维数无关性与边界更正",
      "### 2.41 维数无关性（`R41`）当前状态" in STATUS
      and "维数无关" in STATUS
      and "2.828427" in STATUS
      and "R41_dimension_independence.md" in STATUS
      and "R41_check.py" in STATUS)
check("R42 登记显式 pi 与双侧约束",
      "### 2.42 显式 `π`（旋转类）与两项二值判据的实际裁决（`R42`）当前状态" in STATUS
      and "双侧" in STATUS
      and "R42_explicit_pi_rotation_class.md" in STATUS
      and "R42_check.py" in STATUS)
check("R43 登记两层族构造",
      "### 2.43 `π` 的构造：两层族同时满足双侧约束（`R43`）当前状态" in STATUS
      and "两层族" in STATUS
      and "幂律尾" in STATUS
      and "R43_pi_two_layer_construction.md" in STATUS
      and "R43_check.py" in STATUS)
check("R44 登记生存与语境性互斥 no-go",
      "### 2.44 no-go：选维（生存）与单体语境性互斥（`R44`）当前状态" in STATUS
      and "互斥" in STATUS
      and "3/5" in STATUS
      and "R44_survival_vs_contextuality_no_go.md" in STATUS
      and "R44_check.py" in STATUS)
check("R45 登记账本形式扫描与逃生口",
      "### 2.45 账本形式扫描：`R44` 的 no-go 可逃，代价是 `R32` 唯一性被削弱（`R45`）当前状态" in STATUS
      and "欧氏对" in STATUS
      and "39" in STATUS
      and "R45_ledger_form_scan.md" in STATUS
      and "R45_check.py" in STATUS)
check("R46 登记单纯形导出与价格",
      "### 2.46 对账本的对象与单纯形导出（`R46`）当前状态" in STATUS
      and "单纯形" in STATUS
      and "签名" in STATUS
      and "R46_pair_ledger_objects_and_simplex.md" in STATUS
      and "R46_check.py" in STATUS)
check("Z17 落盘登记（定理 Z17.1 与 §9 更新）",
      "### 2.0z Z17 落盘" in STATUS
      and "定理 Z17.1" in STATUS
      and "Z17_A0_A5_retirement_vacancy_ledger.md" in STATUS
      and "Z17_check.py" in STATUS)

check("R47 登记因果锥供出签名",
      "### 2.47 因果锥供出洛伦兹签名：`R46` 的价格可以付（`R47`）当前状态" in STATUS
      and "洛伦兹签名" in STATUS
      and "离散不变量" in STATUS
      and "R47_signature_from_causal_cone.md" in STATUS
      and "R47_check.py" in STATUS)

# ======================================================================
head("F2  关键入口都指向 STATUS.md")

check("INDEX 有当前状态源标题", "## 0 当前状态源（先读）" in INDEX and "](STATUS.md)" in INDEX)
check("G10 指向 STATUS 且旧 §4/§7 标为历史",
      has("G10_final_derivation_and_input_ledger.md",
          "[`STATUS.md`](STATUS.md)", "§4 的旧终态表", "§7 的未建立清单", "都不是当前状态"))
check("G0 指向 STATUS 且 §4 标为历史",
      has("G0_bottom_layer_and_derivation_route.md",
          "[`STATUS.md`](STATUS.md)", "本节只保留早期输入表", "历史推导记录"))
check("Z6 指向 STATUS 且不再自称唯一账本",
      has("Z6_stall_autopsy_and_released_ledger.md",
          "[`STATUS.md`](STATUS.md)", "不再是项目唯一账本"))
check("Z1 标出 Z0 取代五条公理的历史",
      has("Z1_zero_layer_as_the_foundation.md",
          "[`STATUS.md`](STATUS.md)", "随后 [`Z0`](Z0_zero_never_rests_single_axiom.md) 已把基础压成**唯一公理**"))
check("Z2 标出 I2a 当前为 E1 条件链",
      has("Z2_zero_to_gr_direct_route.md",
          "[`STATUS.md`](STATUS.md)", "当前归入 E1 条件链"))
check("D259 标出路线本地状态",
      has("D259_z_layer_condition_dictionary_and_four_dimensional_confluence.md",
          "[`STATUS.md`](STATUS.md)", "D259 借入路线本地"))
check("借入层说明指向 STATUS",
      has("verify/README.md", "[`STATUS.md`](../STATUS.md)", "D259 路线本地状态"))

# ======================================================================
head("F3  关键旧状态页已历史化，旧唯一账本措辞清零")

for name in [
    "G58_I2a_resolved_as_embedding_input.md",
    "G59_I7_settled_native_cone_and_its_residue.md",
    "G63_target_list_and_audit.md",
    "G70_B_and_tau_closing.md",
]:
    check("%s 已链接唯一状态源" % name, "[`STATUS.md`](STATUS.md)" in read(name))

md = [f for f in os.listdir(HERE) if f.endswith(".md")]
bad = []
for name in md:
    text = read(name)
    if "读者应把本文 §2 当作唯一账本" in text:
        bad.append(name)
check("旧“本文 §2 是唯一账本”措辞为 0", not bad, "遗留: %s" % bad if bad else "无")

check("G10 §7 的旧未建立列表只作历史引用",
      "历史“未建立”清单（已失效）" in read("G10_final_derivation_and_input_ledger.md"))
check("G63 旧计数被 G70 取代",
      "均已被 [`G70`](G70_B_and_tau_closing.md) 取代" in read("G63_target_list_and_audit.md"))

# ======================================================================
head("F4  主路线与 D259 路线分开记账")

check("STATUS 明确两条路线未证明等价",
      "两条路线尚未打通" in STATUS and "进度不能相加" in STATUS and "没有等价性证明" in STATUS)
check("D259 的未解输入没有被写成主路线当前缺口",
      "中的“未解输入”应继续保留" in STATUS and "D259 路线本地未解" in STATUS)
check("旧年龄路线被标为非当前承重链",
      "旧年龄／几何路线" in STATUS and "不再是当前主线的承重链" in STATUS)

# ======================================================================
head("F5  参数与目标清单只采用当前计数")

check("当前目标清单采用 21 / 2 / 3 / 0",
      "| 已关闭／已导出条目 | **21** |" in STATUS and
      "| 在册开放 | **2**" in STATUS and
      "| 约定／口径 | **3**" in STATUS and
      "| 无记录 | **0** |" in STATUS)
check("G63 的两套旧计数被明确列为历史",
      "17 已导出、4 开放" in STATUS and "18 已导出、3 开放" in STATUS)
check("B=4 仍标为候选识别",
      "`B=4` 仍是候选识别，因果论证不成立" in STATUS and "仍是【候选识别】" in read("G70_B_and_tau_closing.md"))

# ======================================================================
head("F6  同步入口运行本核验")

ledger = read("ledger_sync.py")
check("ledger_sync.py 调用 STATUS_check.py", "STATUS_check.py" in ledger)
check("INDEX 登记 STATUS_check.py", "STATUS_check.py" in INDEX)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
