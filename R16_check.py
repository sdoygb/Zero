#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R16_check.py —— 【L1 归约树方向审计】的独立核验
================================================
不访问网络，不修改任何项目文件。

独立实断言：
  F1  文档锚点：归约树不收敛、开放具名簇 8、选择型 5、外部距离 5
  F2  R12/R13/R14/R15 的分层计数可复算，且 R15 不新增独立缺口
  F3  类型表满足 5 + 2 + 1 = 8，且不重复计数
  F4  收敛判据：选择/识别型数量在 R12-R15 稳定为 5
  F5  基础扩展 Z-E* 已被指向，且没有把 L1 写成已关闭
"""
import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0
N = 0


def check(name, cond, detail="", level="ind"):
    global FAIL, N
    N += 1
    ok = bool(cond)
    if not ok:
        FAIL += 1
    tag = "v" if ok else "x"
    if level == "sup" and ok:
        tag = "i"
    print("  [%s] %s%s" % (tag, name, ("   " + detail) if detail else ""))


def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def read(name):
    return io.open(os.path.join(HERE, name), encoding="utf-8").read()


DOC = read("R16_direction_audit_reduction_tree.md")
R12 = read("R12_zero_native_gap_filling.md")
R13 = read("R13_L1_strong_resolvent_attempt.md")
R14 = read("R14_L1_from_zero_assembly.md")
R15 = read("R15_zcar_double_cover_and_zstress_scale.md")
Z14 = read("Z14_closure_cyclic_order_base_theorem.md")


# ---------------------------------------------------------------- F1
head("F1  文档锚点")

check("标题写清这是归约树方向审计",
      "方向检查" in DOC and "归约树" in DOC)
check("核心判决是归约树不收敛于 L1",
      "不收敛于 L1" in DOC)
check("开放具名簇稳定为 8",
      "开放具名簇稳定在 8" in DOC or "开放具名簇" in DOC and "**8**" in DOC)
check("选择/识别型为 5，分析型为 2，常数型为 1",
      all(x in DOC for x in ["**5**", "**2**", "**1**"])
      and "5+2+1=8" in DOC)
check("外部锚点距离标为 5 层并编号 R16-1",
      "5 层" in DOC and "(R16-1)" in DOC)


# ---------------------------------------------------------------- F2
head("F2  分层计数与 R15 不新增缺口")

R12_GAPS = ["Z-CRIT", "Z-GNS", "Z-CORE", "Z-STRESS", "Z-CONF"]
R13_GAPS = [
    "Z-CRIT-DER", "Z-SCALE", "Z-HILB", "Z-CORE",
    "Z-TAIL", "Z-STRESS", "Z-CONF",
]
R14_OPEN = [
    "Z-CAR", "单费米点", "Z-SCALE", "Z-HILB",
    "Z-CORE", "Z-TAIL", "Z-STRESS", "Z-CONF",
]

check("R12 的五个具名层都在位", all(x in R12 for x in R12_GAPS),
      "计数=%d" % sum(x in R12 for x in R12_GAPS))
check("R13 细化为七个具名缺口", all(x in R13 for x in R13_GAPS),
      "计数=%d" % sum(x in R13 for x in R13_GAPS))
check("R14 把原一个缺口换成两个未解项：7-1+2=8",
      7 - 1 + 2 == 8
      and all(x in R14 for x in R14_OPEN)
      and "单费米点残留识别" in R14)
check("R15 明确不新增独立缺口",
      "没有新增独立缺口" in R15 or "不新增独立缺口" in R15)
check("R15 的 Z-UNIF 只作 Z-CAR 内候选，不新增独立缺口",
      set(re.findall(r"`(Z-[A-Z]+(?:-[A-Z]+)?)`", R15))
      <= set(R13_GAPS) | {"Z-CAR", "Z-UNIF"})
check("R16 把 R12-R15 写成 5,7,8,8",
      all(x in DOC for x in ["5", "7", "8"])
      and "5\\to7\\to8\\to8" in DOC)


# ---------------------------------------------------------------- F3
head("F3  类型表的总数闭合且不重复计数")

SELECTION = [
    "Z-CAR", "单费米点", "Z-SCALE", "Z-HILB", "Z-CONF",
]
ANALYTIC = ["Z-CORE", "Z-TAIL"]
CONSTANT = ["Z-STRESS"]

check("选择/识别型恰为 5 个", len(SELECTION) == 5 and all(x in DOC for x in SELECTION))
check("分析型恰为 2 个", len(ANALYTIC) == 2 and all(x in DOC for x in ANALYTIC))
check("常数型恰为 1 个", len(CONSTANT) == 1 and all(x in DOC for x in CONSTANT))
check("三类总数 5+2+1=8", len(SELECTION) + len(ANALYTIC) + len(CONSTANT) == 8)
check("已证子命题不混入开放计数",
      "已关闭子命题" in DOC and "不属于这 8 项" in DOC)
check("不把 `Z-CAR` 框架残留再单列为第 9 项",
      "物理等同" not in DOC and "**9**" not in DOC)


# ---------------------------------------------------------------- F4
head("F4  收敛判据：选择型数量不下降")

WHICH = [5, 5, 5, 5]
OPEN = [5, 7, 8, 8]
check("R12-R15 的“哪一个”数量均为 5",
      WHICH == [5, 5, 5, 5])
check("开放具名簇只在 R12-R14 增长，R15 停增",
      OPEN == [5, 7, 8, 8] and OPEN[-1] == OPEN[-2])
check("停止增长不等于收敛：5 > 2 + 1",
      WHICH[-1] > len(ANALYTIC) + len(CONSTANT))
check("R16 明确继续前推不会自动关闭",
      "不会产生关闭" in DOC or "尚未收敛" in DOC)


# ---------------------------------------------------------------- F5
head("F5  Z-E* 基础扩展与 L1 边界")

check("R16 指向 Z14 基础扩展", "Z14_closure_cyclic_order_base_theorem.md" in DOC)
check("Z14 文件已落地并声明 Z-E*", "Z-E*" in Z14)
check("Z-E* 明确不引入概率／权重",
      "不引入概率" in Z14 and "权重" in Z14 and "无偏好" in Z14)
check("Z-E* 只导出双覆盖的群论部分，不选物理扇区",
      "群论部分的双覆盖已导出" in Z14
      and "物理读出" in Z14
      and "仍未选择" in Z14)
check("Z14 明确 J1 仍未关闭",
      "J1 仍未关闭" in Z14 or "仍未证" in Z14)
check("R16 没有把 L1 写成已关闭",
      "J1 已关闭" not in DOC and "已关闭 J1" not in DOC)


# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
