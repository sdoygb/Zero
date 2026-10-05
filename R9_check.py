#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R9_check.py -- 外部 GR/EFE 恢复路线再排序的独立核验。

对应文档 R9_external_GR_derivations_landscape.md。

  F1  R8 主目标保留 Jacobson 2016，次选改为 Cao-Carroll 2018
  F2  2015-2026 代表路线与全部可核验链接在位
  F3  六个比较维度与路线矩阵在位
  F4  Jacobson 2016 的外部批评被登记
  F5  Oh-Park-Sin、Faulkner、Gorard/Wolfram、Bianconi 的强弱被准确区分
  F6  Zero 的相对强项、弱项和继续攻击顺序在位
  F7  文档没有把外部结果或 Zero 条件恢复写成无条件成功
"""

import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0
NCHECK = 0


def check(name, cond, detail=""):
    global FAIL, NCHECK
    NCHECK += 1
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


DOC = read("R9_external_GR_derivations_landscape.md")

# ======================================================================
head("F1  R8 主目标与次选排序")

for token in [
    "R8 把 Jacobson 2016 作为第一接口的选择仍然成立",
    "次选改为 Cao-Carroll 2018",
    " Jacobson 1995 仍重要",
    "不再作为 Zero 的次选目标",
    "保留第一目标",
    "次选。它不要求先搬入 AdS/CFT",
    "R8 选择 Jacobson 2016 仍是最强接口",
]:
    check("排序锚点：%s" % token, token in DOC)

check("最终选择明确是 Jacobson 2016 与 Cao-Carroll 2018",
      "最适合补的一到两篇是 Jacobson 2016 与 Cao-Carroll 2018" in DOC)

check("Oh-Park-Sin 没有被选为第一或第二目标",
      "不能当 Zero 的第一或第二目标" in DOC
      and "不能当 Zero 的首选接口" in DOC)

# ======================================================================
head("F2  代表路线与可核验链接")

SOURCES = {
    "Jacobson 1995": [
        "10.1103/PhysRevLett.75.1260",
        "gr-qc/9504004",
    ],
    "Jacobson 2016": [
        "10.1103/PhysRevLett.116.201101",
        "1505.04753",
    ],
    "Casini-Galante-Myers 2016": [
        "10.1007/JHEP03(2016)194",
        "1601.00528",
    ],
    "Speranza 2016": [
        "10.1007/JHEP04(2016)105",
        "1602.01380",
    ],
    "Faulkner et al. 2017": [
        "10.1007/JHEP08(2017)057",
        "1705.03026",
    ],
    "Bueno et al. 2017": [
        "10.1103/PhysRevD.95.046003",
        "1612.04374",
    ],
    "Oh-Park-Sin 2017/2018": [
        "10.1103/PhysRevD.98.026020",
        "1709.05752",
    ],
    "Cao-Carroll 2018": [
        "10.1103/PhysRevD.97.086003",
        "1712.02803",
    ],
    "Leichenauer et al. 2018": [
        "10.1103/PhysRevD.98.086013",
    ],
    "Dong-Lewkowycz 2018": [
        "10.1007/JHEP01(2018)081",
        "1705.08453",
    ],
    "Alonso-Serrano-Liska 2020": [
        "10.1103/PhysRevD.102.104056",
        "2008.04805",
    ],
    "Gorard 2020": [
        "10.25088/ComplexSystems.29.2.599",
        "2004.14810",
    ],
    "Wolfram 2020": [
        "10.25088/ComplexSystems.29.2.107",
        "2004.08210",
    ],
    "Carrasco et al. 2023": [
        "10.1007/JHEP09(2023)167",
        "2306.08503",
    ],
    "Kumar 2024/2025": [
        "10.1007/s10714-023-03172-x",
        "2404.16912",
    ],
    "Bianconi 2025": [
        "10.1103/PhysRevD.111.066001",
        "2408.14391",
    ],
}

for name, tokens in SOURCES.items():
    for token in tokens:
        check("%s：%s" % (name, token), token in DOC)

doi_count = len(re.findall(r"https://doi\.org/", DOC))
arxiv_count = len(re.findall(r"https://arxiv\.org/abs/", DOC))
check("至少 15 个 DOI 链接", doi_count >= 15, "count=%d" % doi_count)
check("至少 14 个 arXiv 链接", arxiv_count >= 14, "count=%d" % arxiv_count)
check("明确记录 arXiv API 当时限流及替代核查",
      "Rate exceeded" in DOC and "arXiv 页面 metadata" in DOC)
check("明确声明不采用二手摘要作为结论来源",
      "不把评论文章、博客或二手摘要当作结论来源" in DOC)

# ======================================================================
head("F3  六维比较与路线矩阵")

for token in [
    "连续流形与区域代数",
    "modular / boost 桥",
    "面积密度",
    "固定体积平衡",
    "完整非线性 EFE",
    "可检验预测",
]:
    check("比较维度：%s" % token, token in DOC)

for token in [
    "| Jacobson 1995 |",
    "| Jacobson 2016 |",
    "| Casini-Galante-Myers 2016 |",
    "| Speranza 2016 |",
    "| Faulkner et al. 2017 |",
    "| Bueno-Min-Speranza-Visser 2017 |",
    "| Oh-Park-Sin 2017/2018 |",
    "| Cao-Carroll 2018 |",
    "| Leichenauer-Levine-Shahbazi-Moghaddam 2018 |",
    "| Dong-Lewkowycz 2018 |",
    "| Alonso-Serrano-Liska 2020 |",
    "| Gorard 2020 + Wolfram 2020 |",
    "| Carrasco-Pedraza-Svesko-Weller-Davies 2023 |",
    "| Kumar 2024/2025 |",
    "| Bianconi 2025 |",
]:
    check("矩阵行：%s" % token, token in DOC)

check("论文强弱界定包含完整非线性、微扰、弱场、线性化与修改重力",
      all(x in DOC for x in ["完整", "微扰", "弱场", "线性化", "修改重力"]))

# ======================================================================
head("F4  Jacobson 2016 的关键批评")

for token in [
    "Casini, Galante, Myers",
    "R^{2\\Delta}\\delta\\langle O_\\Delta\\rangle^2",
    "\\Delta\\le \\frac d2",
    "Speranza",
    "低维污染不是技术小项",
    "若 Zero 的局部相关算符有效维数接近二",
    "在补几何 boost 之前，先证明低维相关算符不会污染固定体积小球首阶熵变。",
]:
    check("Jacobson 批评锚点：%s" % token, token in DOC)

check("批评没有被淡化成一阶技术项",
      "可以改变小球熵的首阶标度" in DOC and "R9 的结论更强" in DOC)

# ======================================================================
head("F5  强路线强弱区分")

for token in [
    "微扰到二阶",
    "明确说明，对高阶曲率引力，小球线性化方程不能像 Einstein 情形那样推出完整非线性场方程",
    "只是把“怎样从量子信息得到完整 EFE”改名成“怎样从底层得到全阶 GDERE”",
    "Gorard 论文的摘要明确使用“更新规则在极限中保持因果图维数”这一假设",
    "低耦合时回到零宇宙学常数的 Einstein 方程",
    "热力学路线自然给出 unimodular gravity，而不是完整 GR",
    "大部分证明在二维膨胀引力中完成",
    "现阶段不能据此宣称 Jacobson 路线已被无假设替代",
]:
    check("强弱锚点：%s" % token, token in DOC)

check("Faulkner 被定为外部强结果而不是 Zero 首选",
      "作为“完整非线性恢复的数学上限”阅读，不作为 Zero 下一条主链" in DOC)
check("Gorard/Wolfram 被定为离散对照",
      "离散对照，不作为第一或第二接口" in DOC)
check("Bianconi 被定为最新强路线但不是当前直接补前提目标",
      "最新强路线之一" in DOC and "当前不适合作为 Zero 的直接补前提目标" in DOC)

# ======================================================================
head("F6  Zero 的强项、弱项与下一步")

for token in [
    "原语更少",
    "账本更细",
    "有限维第一定律已闭合",
    "有可证伪残余",
    "对抗审计习惯",
    "连续区域代数网未构造",
    "几何 boost 未证",
    "固定体积平衡未连续化",
    "面积密度尚未普适化",
    "低维相关算符控制缺失",
    "Lorentzian 四维组装不完整",
    "绝对单位映射仍是 E4",
    "先攻低维相关算符污染",
    "再攻算子级几何 boost",
    "同时补面积密度",
    "最后才做固定体积连续变分",
    "把 Zero 的零和图与互信息图对应起来",
    "只在弱场范围内声称 Einstein 方程",
]:
    check("Zero 强弱/工作顺序：%s" % token, token in DOC)

check("明确列出不应先做的全息、离散与无条件升格错误",
      all(x in DOC for x in [
          "不先搬入 AdS/CFT 或 Ryu-Takayanagi",
          "不把 Oh-Park-Sin 的全阶 GDERE 当作已给输入",
          "不把外部论文的“条件恢复”写成项目自己的“无条件导出”",
      ]))

# ======================================================================
head("F7  诚实边界")

for token in [
    "外部路线审计，不修改 `STATUS.md`",
    "Zero 目前仍只支持条件恢复，不支持无条件导出四维 GR",
    "还不能拿出手的是",
    "“Zero 已无条件导出四维 GR”",
    "“Zero 已补完 Jacobson 2016”",
    "“R8 与 Oh-Park-Sin 或 Faulkner 路线等价”",
    "“主 Z/G 路线与 D259 图册路线已经汇流”",
]:
    check("诚实边界锚点：%s" % token, token in DOC)

forbidden = [
    "Zero 已经无条件导出四维 GR。",
    "Zero 已经补完 Jacobson 2016。",
    "主 Z/G 路线与 D259 图册路线已经汇流。",
]
check("没有出现项目级无条件成功宣告", not any(x in DOC for x in forbidden))

# ======================================================================
head("汇总")
print("  断言 %d 项，不符项：%d" % (NCHECK, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
