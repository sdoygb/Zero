#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R8_L1_check.py -- R8 C2/R8-GMA L1 对抗审计的离线锚点与最小代数核验。

检查对象:
  - R8_L1_refutation_attempt.md 的判决、充分条件和失败机制
  - R8/G62/G75/G79/D231/D234 的原始锚点
  - 抛物型几何权重的边界退化与表面重力归一
  - "同一支持/剖面相关 => 算子相等" 的最小反例

不访问网络，不修改任何项目文件。
"""

from __future__ import annotations

import io
import os
import sys

import numpy as np



HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = "R8_L1_refutation_attempt.md"

LOCAL_FILES = [
    REPORT,
    "R8_jacobson_entanglement_equilibrium_completion.md",
    "G62_quantum_sector_from_GNS_modular_flow_gleason.md",
    "G75_quantum_geometry_modular_readout.md",
    "G79_horizon_thermodynamics.md",
    "D231_modular_density_profile_gap.md",
    "D234_geometric_ball_profile_candidate.md",
]

checks: list[tuple[str, bool, str]] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    checks.append((name, bool(condition), detail))


def read_local(name: str) -> str:
    path = os.path.join(HERE, name)
    with io.open(path, encoding="utf-8") as handle:
        return handle.read()


def all_tokens(text: str, tokens: list[str]) -> bool:
    return all(token in text for token in tokens)


def main() -> int:
    for name in LOCAL_FILES:
        check("存在文件 %s" % name, os.path.isfile(os.path.join(HERE, name)))

    if not os.path.isfile(os.path.join(HERE, REPORT)):
        for name, ok, detail in checks:
            print("[%s] %s%s" % ("v" if ok else "x", name, ("  " + detail) if detail else ""))
        print("\nASSERTIONS=%d" % len(checks))
        print("FAILURES=%d" % sum(not ok for _, ok, _ in checks))
        print("RESULT=FAIL")
        return 1

    report = read_local(REPORT)
    r8 = read_local("R8_jacobson_entanglement_equilibrium_completion.md")
    g62 = read_local("G62_quantum_sector_from_GNS_modular_flow_gleason.md")
    g75 = read_local("G75_quantum_geometry_modular_readout.md")
    g79 = read_local("G79_horizon_thermodynamics.md")
    d231 = read_local("D231_modular_density_profile_gap.md")
    d234 = read_local("D234_geometric_ball_profile_candidate.md")

    # ------------------------------------------------------------------
    # A. 原始 R8 目标与张力
    # ------------------------------------------------------------------
    check(
        "R8 原始 C2 锚点",
        all_tokens(
            r8,
            [
                "C2｜强图／预解意义的几何模流极限",
                "K_{B,a}",
                "2\\pi[B_B,P]",
                "J1｜BW/几何 boost 极限",
                "未证",
            ],
        ),
    )
    check(
        "R8 明列几何 boost 未证",
        "小球的几何模流极限 $K_B\\to 2\\pi B_B$" in r8
        and "**开放**" in r8,
    )

    # ------------------------------------------------------------------
    # B. 本地材料进度锚点
    # ------------------------------------------------------------------
    check(
        "G62 有限维 GNS/模流锚点",
        all_tokens(g62, ["有限维忠实态", "模流", "Tomita–Takesaki"]),
    )
    check(
        "G62 未冒充连续区域代数",
        "有限维年龄代数" in r8 and "空间区域代数未构造" in r8,
    )
    check(
        "G75 Gibbs 态条件锚点",
        all_tokens(g75, ["取 Gibbs 态", "\\beta L_W", "模 Hamiltonian"]),
    )
    check(
        "G75 不把 Gibbs 读数写成区域真空定理",
        "取 Gibbs 态" in g75 and "几何模条件" not in g75[:180],
    )
    check(
        "G79 边界退化和 0.865 如实登记",
        all_tokens(g79, ["边界", "0.865", "CFT 的精确抛物线", "**不是**"]),
    )
    check(
        "G79 把几何识别列为输入",
        "几何模条件" in g79
        and "boost 生成元的识别仍是【输入】" in g79,
    )
    check(
        "D231 同一支持不决定剖面/生成元",
        all_tokens(
            d231,
            [
                "不能唯一选出 boost 剖面",
                "支持数据不能区分不同 boost 候选剖面",
                "本文没有导出上述任一项",
            ],
        ),
    )
    check(
        "D234 抛物核与接触项均为候选/缺口",
        all_tokens(
            d234,
            [
                "\\frac{R^2-r^2}{2R}",
                "接触项",
                "未由零层导出",
                "格子上的 Bisognano-Wichmann 形式只是这些连续核的离散近似",
            ],
        ),
    )

    # ------------------------------------------------------------------
    # C. 报告必须给出的判决和条件
    # ------------------------------------------------------------------
    for verdict in ["已证", "条件证成", "开放", "排除"]:
        check("报告含判决类别：%s" % verdict, verdict in report)

    for condition in ["C-L1a", "C-L1b", "C-L1c", "C-L1d",
                      "C-L1e", "C-L1f", "C-L1g", "C-L1h"]:
        check("报告含充分条件 %s" % condition, condition in report)

    check(
        "报告明确区分流、生成元和范数三种目标",
        all_tokens(report, ["流的相等", "生成元在公共核心上的相等", "算子范数收敛"]),
    )
    check(
        "报告指出字面无界算子范数不良定义",
        "字面无界算子范数" in report and "不是待证明定理，而是待定义的命题" in report,
    )
    check(
        "报告修正局部 boost 为共形几何生成元",
        "conformal Killing field" in report and "区域保持共形 boost/几何生成元" in report,
    )

    # ------------------------------------------------------------------
    # D. 外部定理与反例桥
    # ------------------------------------------------------------------
    check(
        "Jacobson 2016 外部锚点",
        all_tokens(report, ["1505.04753", "10.1103/PhysRevLett.116.201101", "dKC2"]),
    )
    check(
        "Bisognano-Wichmann 外部锚点",
        all_tokens(report, ["10.1063/1.522605", "10.1063/1.522898", "Rindler wedge"]),
    )
    check(
        "Hislop-Longo 外部锚点",
        all_tokens(report, ["10.1007/BF01208372", "自由无质量标量", "任意维"]),
    )
    check(
        "BGL 外部定理与全部关键前提",
        all_tokens(
            report,
            [
                "10.1007/BF02096738",
                "funct-an/9302008",
                "isotony",
                "causality",
                "additivity",
                "Delta_O^{it}=U_O(t)",
                "正能",
            ],
        ),
    )
    check(
        "Fredenhagen 只给角落/渐近结论",
        "10.1007/BF01206179" in report and "角落极限" in report,
    )
    check(
        "Brunetti-Moretti 非共形失败机制",
        all_tokens(
            report,
            [
                "1009.4990",
                "L^0_{1,1}",
                "一般不为零",
                "非共形质量项",
            ],
        ),
    )
    check(
        "报告不把候选剖面当作几何 boost",
        "候选剖面相关不得当作几何 boost" in report,
    )
    check(
        "报告没有越过现有 Zero 基础",
        "现有 Zero 基础不能补上 C2" in report,
    )

    # 直接禁止的越界表述。报告可在否定句中讨论这些概念，但不能用这些短句下结论。
    forbidden = [
        "L1 已证",
        "R8-GMA 已证明",
        "候选剖面就是几何 boost",
        "已从 Zero 无条件推出 R8-GMA",
        "已无条件导出四维 GR",
    ]
    for phrase in forbidden:
        check("报告不含越界短语：%s" % phrase, phrase not in report)

    # ------------------------------------------------------------------
    # E. 抛物型候选核的基本代数性质
    # ------------------------------------------------------------------
    radius = 2.0
    radii = np.linspace(0.0, radius, 101)
    profile = (radius**2 - radii**2) / (2.0 * radius)
    derivative = -radii / radius
    second_derivative = np.full_like(radii, -1.0 / radius)

    check("抛物型核中心值 R/2", np.isclose(profile[0], radius / 2.0))
    check("抛物型核边界为零", np.isclose(profile[-1], 0.0))
    check("抛物型核严格凹", bool(np.all(second_derivative < 0.0)))
    check(
        "抛物型核在内部严格下降",
        bool(np.all(np.diff(profile) < 0.0)),
    )
    check(
        "边界斜率大小为一（表面重力归一）",
        np.isclose(abs(derivative[-1]), 1.0),
    )

    # ------------------------------------------------------------------
    # F. 最小矩阵反例
    #   同一支持上的不同剖面给出不同交换子。
    # ------------------------------------------------------------------
    sites = np.array([-1.5, -0.5, 0.5, 1.5])
    local_radius = 2.0
    geometric_diag = (local_radius**2 - sites**2) / (2.0 * local_radius)
    constant_diag = np.ones_like(sites)

    phi = np.array(
        [
            [0.0, 1.0, 0.0, 0.0],
            [1.0, 0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 0.0],
        ],
        dtype=float,
    )
    k_geom = np.diag(geometric_diag)
    k_flat = np.diag(constant_diag)

    comm_geom = k_geom @ phi - phi @ k_geom
    comm_flat = k_flat @ phi - phi @ k_flat
    mismatch = comm_geom - comm_flat

    check("最小反例中 phi 自伴", np.allclose(phi, phi.T))
    check(
        "同一支持常量剖面是单位矩阵",
        np.allclose(k_flat, np.eye(4)),
    )
    check(
        "几何与常量剖面支持相同",
        np.all(np.abs(geometric_diag) > 0.0)
        and np.all(np.abs(constant_diag) > 0.0),
    )
    check(
        "最小反例给出非零交换子差异",
        np.linalg.norm(mismatch, ord="fro") > 1e-12,
        "mismatch_fro=%.6f" % np.linalg.norm(mismatch, ord="fro"),
    )
    check(
        "命题不是流相等，而是非唯一性反例",
        "这个反例不否定真实 R8-GMA" in report
        and "算子级几何 boost" in report,
    )

    failures = sum(not ok for _, ok, _ in checks)
    for name, ok, detail in checks:
        suffix = ("  " + detail) if detail else ""
        print("[%s] %s%s" % ("v" if ok else "x", name, suffix))
    print()
    print("ASSERTIONS=%d" % len(checks))
    print("FAILURES=%d" % failures)
    print("RESULT=%s" % ("PASS" if failures == 0 else "FAIL"))
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
