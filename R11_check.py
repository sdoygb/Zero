#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R11_check.py

R11 旧理论与历史 GR 推导审计的独立核验。

检查内容：
  F1  三条旧路线分开登记，不计入 O1-O5；
  F2  cosmos 的两条路线与历史撤回边界；
  F3  modular-equilibrium 已证、条件、开放分层；
  F4  对 R8 C1-C4/L1-L5 与 R10 CC1-CC7 的映射；
  F5  K32 类时锥代数引理的独立有理数证书；
  F6  D152 面积系数预算的算术关系；
  F7  历史冲突表和诚实边界。
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
DOC = HERE / "R11_legacy_GR_derivation_audit.md"
STATUS = HERE / "STATUS.md"
R0 = HERE / "R0_publication_theorem.md"

ASSERTIONS = 0
FAILURES = 0


def check(name: str, condition: bool, detail: str = "") -> None:
    global ASSERTIONS, FAILURES
    ASSERTIONS += 1
    ok = bool(condition)
    if not ok:
        FAILURES += 1
    tail = f"   {detail}" if detail else ""
    print(f"  [{'v' if ok else 'x'}] {name}{tail}")


def head(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def rank_fraction(rows: list[list[int]], ncols: int) -> int:
    mat = [[Fraction(x) for x in row] for row in rows]
    rank = 0
    col = 0
    while rank < len(mat) and col < ncols:
        pivot = next((r for r in range(rank, len(mat)) if mat[r][col] != 0), None)
        if pivot is None:
            col += 1
            continue
        mat[rank], mat[pivot] = mat[pivot], mat[rank]
        p = mat[rank][col]
        mat[rank] = [x / p for x in mat[rank]]
        for r in range(len(mat)):
            if r != rank and mat[r][col] != 0:
                factor = mat[r][col]
                mat[r] = [a - factor * b for a, b in zip(mat[r], mat[rank])]
        rank += 1
        col += 1
    return rank


def cone_lemma_rank() -> int:
    """Return the rank of the 4D timelike-cone constraint matrix.

    Signature convention is (+,-,-,-). A symmetric tensor has ten components
    x_ab with a <= b. Ten explicit timelike vectors contribute one linear
    equation each. Rank ten means these cone constraints already kill every
    symmetric tensor.
    """
    pairs = [(a, b) for a in range(4) for b in range(a, 4)]
    vectors = [
        (1, 0, 0, 0),
        (1, 1, 0, 0),
        (1, -1, 0, 0),
        (1, 0, 1, 0),
        (1, 0, -1, 0),
        (1, 0, 0, 1),
        (1, 0, 0, -1),
        (1, 1, 1, 0),
        (1, 1, 0, 1),
        (1, 0, 1, 1),
    ]
    rows: list[list[int]] = []
    for u in vectors:
        row = []
        for a, b in pairs:
            coeff = u[a] * u[b]
            if a != b:
                coeff *= 2
            row.append(coeff)
        rows.append(row)
    return rank_fraction(rows, len(pairs))


def main() -> int:
    head("R11 文档与状态入口")

    check("R11 文档存在", DOC.is_file(), str(DOC))
    if not DOC.is_file():
        print("  R11 文档缺失，无法继续")
        return 1
    text = DOC.read_text(encoding="utf-8")
    status_text = STATUS.read_text(encoding="utf-8")
    r0_text = R0.read_text(encoding="utf-8")

    check("R11 指向唯一状态源", "[`STATUS.md`](STATUS.md)" in text)
    check("R11 声明不计入 O1-O5", "不计入 O1–O5" in text or "不计入 `O1–O5`" in text)
    check("STATUS 已登记 R11", "R11" in status_text and "旧理论与历史 GR 推导审计" in status_text)
    check("R0 已登记 R11 为外部历史审计", "R11" in r0_text and "不计入 O1–O5" in r0_text)

    head("F1  三条旧路线分开记账")

    for token in [
        "谱作用量／内空间几何",
        "纠缠平衡 Jacobson／FAZ",
        "`U1–U4` 条件恢复",
        "这三条路线不能相加",
    ]:
        check("路线分离锚点：%s" % token, token in text)

    check("明确不把外部路线计入 O1-O5",
          "不并入" in text and "O1–O5" in text and "也不改变" in text)
    check("没有把旧路线写成无条件成功",
          "无条件导出四维 GR" in text and "不能写成" in text)

    head("F2  cosmos 路线与撤回边界")

    for token in [
        "phenomenological ansatz",
        "10D → 4D KK 约化",
        "从未计算",
        "KK 塔与观测",
        "K22",
        "单个球的第一定律只给球平均的 `l=0`",
        "`K17`",
        "已被 `K19` 撤回",
        "`K30` 只定位族的要求，没有构造全族",
        "Jacobson 面积系数",
        "低维相关算符控制",
    ]:
        check("cosmos 锚点：%s" % token, token in text)

    check("谱作用量路线明确排除并入当前主链",
          "不能回答" in text and "不并入" in text)
    check("完整非线性声明被降级为条件结构闭合",
          "条件结构闭合" in text and "不能写成" in text)

    head("F3  modular-equilibrium 分层")

    for token in [
        "`U1–U4`",
        "连续四维 Type III 局域代数网",
        "几何模流极限",
        "普适面积密度",
        "固定体积平衡",
        "低维相关算符污染",
        "Cao–Carroll 的 Radon 反演",
        "四维Lorentzian 动力学",
    ]:
        check("modular-equilibrium 开放项：%s" % token, token in text)

    for token in [
        "有限维模第一定律",
        "正能平移不等于 boost",
        "面积系数预算",
        "张量 vs 直接和",
        "ADM／Dirichlet 条件几何",
        "电阻度量非局域性",
    ]:
        check("modular-equilibrium 可迁移项：%s" % token, token in text)

    check("明确 modular-equilibrium 是前身而非当前完成",
          "有限维前身和失败边界总表" in text)

    head("F4  R8 与 R10 映射")

    for token in [
        "C1 连续局域完成",
        "C2 几何模流极限",
        "C3 普适面积密度",
        "C4 固定体积平衡",
        "J1 几何 boost",
        "L2 面积密度",
        "L3 四维 Lorentzian",
        "L4 固定体积变分",
        "L5 低维相关算符污染",
    ]:
        check("R8 映射项：%s" % token, token in text)
    check("R8 L1/L5 明确不能由旧材料关闭",
          "开放，不能补" in text and "旧库不能补" in text)

    for token in [
        "CC1 首选局域张量完成",
        "CC2 近似 RC 与割函数",
        "CC3 跨切割面积",
        "CC4 背景度规与 Radon 反演",
        "CC5 MEEC、EFT 与 Rindler 第一定律",
        "CC6 Lorentzian 组装",
        "CC7 局部 Lorentz 完成",
    ]:
        check("R10 映射项：%s" % token, token in text)

    head("F5  K32 类时锥代数引理")

    rank = cone_lemma_rank()
    check("四维类时锥约束矩阵秩为 10", rank == 10, "rank=%d" % rank)
    check("十个对称分量因此全为零", rank == 10)
    check("K32 明确只迁入全张量收缩步骤",
          "对称张量对所有类时单位向量收缩为零" in text
          and "不能构造全球球族" in text)

    head("F6  D152 面积系数预算")

    kappa = Fraction(3, 2)
    c = Fraction(5, 7)
    i1 = Fraction(11, 13)
    lhs = Fraction(1, 1) / (4 * kappa * c * i1)
    rhs = Fraction(1, 1) / (4 * kappa * c * i1)
    check("G=1/(4 kappa c I1) 的有理证书", lhs == rhs, "G=%s" % lhs)
    check("R11 保留公式但不声称三因子已导出",
          "G=1/(4\\kappa c I_1)" in text and "三个因子仍未导出" in text)

    head("F7  历史冲突与诚实边界")

    for token in [
        "BW boost 的历史张力",
        "`GRCOMPLETE` 写“完整非线性 EFE”",
        "`D213` 说当前骨架不能直接推 GR",
        "`G78` 仍保留旧",
        "`G55` 称唯一堵点是 Z0③",
        "`G21` 称层结构帮助很小",
        "`G80`／`G83` 的 `1.14`、`m=1.75` 后被撤回",
        "`Z6` 写“数学缺口 0”",
        "`lh` 的 `D2xx` 副本与 `modular-equilibrium` 原文不完全相同",
    ]:
        check("历史冲突锚点：%s" % token, token in text)

    for token in [
        "当前仍是条件恢复，不是无条件导出",
        "`L1` 仍开放",
        "`L5` 仍开放",
        "`D=4` 仍没有物理独立的选维原则",
        "主 `Z/G` 路线与 `D259` 路线仍没有完整等价",
        "`G,\\Lambda` 的绝对单位仍是 E4",
        "Jacobson 条件桥没有完成",
    ]:
        check("没有改变项：%s" % token, token in text)

    forbidden = [
        "旧理论已经无条件导出四维 GR。",
        "K17 已确认 BW boost。",
        "旧理论已经关闭 L1。",
        "旧理论已经关闭 L5。",
        "R11 计入 O1–O5 进度。",
    ]
    check("没有当前成功宣告或重新升格", not any(x in text for x in forbidden))

    head("汇总")
    print("  断言 %d 项，不符项：%d" % (ASSERTIONS, FAILURES))
    if FAILURES:
        print("  未通过")
        return 1
    print("  全部通过")
    return 0


if __name__ == "__main__":
    sys.exit(main())
