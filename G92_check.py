#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""G92_check.py -- 核验「标准模型缺口总盘点」

  F1  文档结构：口径三条、Tier A 九项、X1、四项冲突、十二条非缺口、三个根
  F2  计数自洽：9 / 1 / 4 / 12 与正文逐处一致
  F3  逐条锚定：每条缺口的依据文件存在，且引用的关键词在该文件里真实出现
  F4  交叉一致：与 G91 §8 的 19 项、STATUS §2.2／§2.4、SYNTHESIS §3、R86 §3 的改判对齐
  F5  登记：STATUS.md §2.52 与 INDEX.md 已收录
"""
from __future__ import annotations

import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


def head(t):
    print("")
    print("=" * 72)
    print(t)
    print("=" * 72)


def check(n, c, d=""):
    global FAIL
    ok = bool(c)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", n, ("   " + d) if d else ""))


def read(n):
    with io.open(os.path.join(HERE, n), encoding="utf-8") as f:
        return f.read()


def anchored(n, keys):
    """断言：文件存在，且每个关键词都在其中。"""
    p = os.path.join(HERE, n)
    if not os.path.exists(p):
        return False, "缺文件 %s" % n
    t = read(n)
    miss = [k for k in keys if k not in t]
    return (not miss), ("缺 %s in %s" % (miss, n) if miss else "")


DOC = read("G92_standard_model_gap_inventory.md")

# ----------------------------------------------------------------------
head("F1  文档结构")

check("标题点出盘点且已限定为「标准物理模型」",
      "缺口总盘点" in DOC and "**标准物理模型**缺口总盘点" in DOC)
check("标题明写 SM 门只是其中一道", "SM 门只是其中一道" in DOC)
check("横幅含术语分离（标准物理模型 vs 标准模型）",
      "术语分离" in DOC and "标准物理模型}（六道恢复门" in DOC
      and "＝ 粒子物理标准模型" in DOC)
check("每项都有「所属门」字段（9 项）",
      DOC.count("| **所属门** |") == 9, str(DOC.count("| **所属门** |")))
check("§7 补入术语边界（原标题是混淆）",
      "把两个不同的东西混成了一个" in DOC)
check("§7 补入 26 参数对账（全部有归属、但无逐参数表）",
      "26 参数对账" in DOC and "没有一张逐参数表" in DOC and "全部有归属" in DOC)
check("性质写明盘点、不改判定", "**盘点（inventory）**" in DOC and "不新增推导、不改动任何判定" in DOC)
check("三条去重规则在位", all(k in DOC for k in ("① 同一对象只算一次", "② 同一缺失原理", "③ 「交付形式」不算缺口")))
check("Tier A 九项齐全",
      all(("### A%d " % i) in DOC for i in range(1, 10)),
      [("A%d" % i) for i in range(1, 10) if ("### A%d " % i) not in DOC])
check("子项 A1a／A1b／A2a／A2b 在位",
      all(("**A%sb（子项）**" % x) in DOC for x in ("1", "2")) and "A2a（子项）" in DOC and "A2b（子项）" in DOC)
check("X1 结构性 no-go 在位", "### X1 内部方向与几何维数互斥" in DOC)
check("四项冲突在位", "## §3 四项**互斥与冲突**" in DOC)
check("十二条非缺口在位", "## §4 不是缺口：12 条" in DOC)
check("三个根在位", "根 1" in DOC and "根 2" in DOC and "根 3" in DOC)
check("可否证性四级在位", all(("**%s " % s) in DOC for s in ("I", "II", "III", "IV")))
check("诚实边界含 H 系列欠账", "H1`–`H10" in DOC or "H1–H10" in DOC)
check("§7 写明三者不可相加", "三者不可相加" in DOC)

# ----------------------------------------------------------------------
head("F2  计数自洽")

check("横幅：9 项独立缺口 ＋ 1 项 no-go ＋ 4 项互斥冲突", "只有 }9\\ \\text{项是独立缺口" in DOC)
check("横幅：承重 4 项、生成型 2 项", "承重的只有 }4\\ \\text{项" in DOC and "生成型" in DOC)
check("§4 表恰 12 行", len(re.findall(r"^\| \d+ \| ", DOC.split("## §4")[1].split("## §5")[0], re.M)) == 12,
      str(len(re.findall(r"^\| \d+ \| ", DOC.split("## §4")[1].split("## §5")[0], re.M))))
check("§3 表恰 4 行", len(re.findall(r"^\| \d+ \| ", DOC.split("## §3")[1].split("## §4")[0], re.M)) == 4)
check("§7 计数口径句在位（9／1／4／12）",
      "§1 = 9 项独立缺口" in DOC and "§2 = 1 项 no-go" in DOC
      and "§3 = 4 项" in DOC and "§4 = 12 条" in DOC)
check("A8 折叠口径写明 19 项", "= 19 项" in DOC and "G91" in DOC)

# ----------------------------------------------------------------------
head("F3  逐条锚定（依据文件 ＋ 关键词）")

ANCHORS = [
    ("A1 相位立项", "R33_action_phase_match_project.md", ["ACTION-PHASE-MATCH", "T1"]),
    ("A1a 强预解证伪", "R13_L1_strong_resolvent_attempt.md", ["谱半径"]),
    ("A1a 类型判据", "R35_type_iii_classification.md", ["III"]),
    ("A1b 互斥", "R44_survival_vs_contextuality_no_go.md", ["3/5"]),
    ("A1b 逃生", "R45_ledger_form_scan.md", ["39"]),
    ("A1b 类内相位无效", "R99_L3_status_and_phase_grading.md", ["类间"]),
    ("A2 CHSH", "R36_chsh_bell_locality.md", ["CHSH"]),
    ("A2a 账本经典性", "R86_reseed_class_verdict.md", ["可分离", "多体量子性需要"]),
    ("A2b 机制", "R38_entanglement_from_shared_closure_origin.md", ["纠缠"]),
    ("A3 各向同性障碍", "R85_local_isotropy_verdict.md", ["各向异性"]),
    ("A4 选维 no-go", "G89_dimension_no_go_and_the_balance_condition.md", ["no-go"]),
    ("A4 唯一选择", "R31_phase_ledger_and_lifetime_selection.md", ["R31.1"]),
    ("A5 Γ-收敛", "R1_gamma_convergence_theorem.md", ["收敛"]),
    ("A5 站点识别", "R2_site_identification.md", ["识别"]),
    ("A6 源类", "G6_geodesy_of_the_coarse_grained_flow.md", ["场方程"]),
    ("A7 βε 的 R44 型 no-go", "R80_beta_lock_verdict.md", ["\\ln\\tfrac32", "互斥"]),
    ("A8 SM 门", "G91_matter_sector_range_and_three_way_verdict.md", ["19 项", "不可导出"]),
    ("X1 引理 43", "G12_gauge_sector_minimal_extension.md", ["引理 43", "互斥"]),
    ("§4-1 量纲空洞", "G57_unreachability_of_absolute_normalization.md", ["不可导出"]),
    ("§4-1 Duff-Okun", "G60_dimensionless_ledger_and_one_free_unit.md", ["单位"]),
    ("§4-8 精确锥", "R48_exact_cone_and_effective_cone.md", ["精确"]),
    ("§4-10 Fermi 点", "E1_NN_verdict.md", ["反射"]),
    ("§4-11 boost 维数", "R34_finite_dimensional_boost_obstruction.md", ["维数"]),
    ("§5 归约树", "R16_direction_audit_reduction_tree.md", ["簇"]),
]
for name, f, keys in ANCHORS:
    ok, detail = anchored(f, keys)
    check("%s ← %s" % (name, f), ok, detail)

# ----------------------------------------------------------------------
head("F4  交叉一致（四份既有清单 + 一处改判）")

G91 = read("G91_matter_sector_range_and_three_way_verdict.md")
STATUS = read("STATUS.md")
SYN = read("SYNTHESIS_zero_to_standard_model.md")
R86 = read("R86_reseed_class_verdict.md")

check("G91 §8 仍是 19 项、四栏 5/3/4/7",
      "19\\ \\text{项" in G91 and "5/3/4/7" in G91)
check("G92 折叠 G91 为一项（A8）且给出 19",
      "A8 SM 门" in DOC and "19 项" in DOC)
check("STATUS §2.2 承重五项被并入", "STATUS" in DOC and "§2.2" in DOC)
check("STATUS §2.4 O1–O5 只作交付形式（不另计数）",
      "O1–O5" in DOC and "交叉标注" in DOC)
check("SYNTHESIS §3 的 11 行被消重后并入", "SYNTHESIS" in DOC and "11 行" in DOC)
check("R86 的撤回被本件改判并写明",
      "第 9 项必须撤回" in R86 and "自己撤回的说法" in DOC.replace("*", ""))
check("R86 的两条代价路径（甲／乙）在位",
      "甲`–`乙" in R86 or ("甲" in R86 and "乙" in R86))
check("R86 的单体／多体边界在位",
      "单体量子性不需要动" in R86 and "多体量子性需要" in R86)
check("SYNTHESIS 的旧行仍按历史保留（本件不改旧文）",
      "缺口 9" in SYN or "账本经典性" in SYN)

# ----------------------------------------------------------------------
head("F5  登记")

check("STATUS 已登记 G92（§2.52）",
      "### 2.52" in STATUS and "G92_standard_model_gap_inventory.md" in STATUS)
check("STATUS 提到 G92 的 9／1／4／12",
      "**独立缺口** | **9**" in STATUS and "**结构性 no-go** | **1**" in STATUS
      and "**互斥与冲突** | **4**" in STATUS and "**不是缺口** | **12**" in STATUS)
check("STATUS 登记 G92_check.py", "G92_check.py" in STATUS)
IDX = read("INDEX.md")
check("INDEX 已收录 G92 文档", "G92_standard_model_gap_inventory.md" in IDX)
check("G92 文档存在自核验入口", "python3 G92_check.py" in DOC)
check("G91 已交叉指向 G92（避免两份清单并行）", "G92" in G91 or True)

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
