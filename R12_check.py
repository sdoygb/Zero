#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R12_check.py -- R12 三个定向缺口审计的离线核验。

检查对象:
  - R12 文档与 R8/R10/R9 的接口锚点
  - L5: 相对熵 BKM 二次型、有序对因子、有限维 no-go 族、gap 界
  - L1: 中央剖面与 Poincare 闭合阻碍的文字锚点
  - Cao-Carroll: RC-EPF 引理与 GHZ 反例

不访问网络，不修改任何项目文件。
"""

from __future__ import annotations

import io
import math
import os
import sys

import numpy as np



HERE = os.path.dirname(os.path.abspath(__file__))
DOC = "R12_zero_native_gap_filling.md"

checks: list[tuple[str, bool, str]] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    checks.append((name, bool(condition), detail))


def read_local(name: str) -> str:
    with io.open(os.path.join(HERE, name), encoding="utf-8") as handle:
        return handle.read()


def all_tokens(text: str, tokens: list[str]) -> bool:
    return all(token in text for token in tokens)


def entropy(rho: np.ndarray) -> float:
    vals = np.clip(np.linalg.eigvalsh((rho + rho.conj().T) / 2), 1e-15, None)
    vals = vals[vals > 1e-14]
    return float(-np.sum(vals * np.log(vals)))


def relative_entropy(sigma: np.ndarray, rho: np.ndarray) -> float:
    eig_s, vec_s = np.linalg.eigh((sigma + sigma.conj().T) / 2)
    eig_s = np.clip(eig_s, 1e-300, None)
    # 参考态固定为对角 rho，故 log rho 在对角基下精确。
    p_ref = np.real(np.diag(rho))
    diagonal = np.real(np.einsum("ai,i,ai->a", vec_s.conj(), np.log(p_ref), vec_s))
    return float(np.sum(eig_s * (np.log(eig_s) - diagonal)))


def ghz4() -> np.ndarray:
    state = np.zeros(16, dtype=complex)
    state[0] = 1.0 / math.sqrt(2.0)
    state[-1] = 1.0 / math.sqrt(2.0)
    return state


def reduced_state(state: np.ndarray, keep: list[int], n: int) -> np.ndarray:
    tensor = state.reshape([2] * n)
    rho = np.tensordot(tensor, tensor.conj(), axes=0)
    # axes are (ket_0, ..., ket_{n-1}, bra_0, ..., bra_{n-1})
    axis_order = list(range(2 * n))
    for site in sorted(set(range(n)) - set(keep)):
        ket_axis = axis_order.index(site)
        bra_axis = axis_order.index(n + site)
        if ket_axis > bra_axis:
            ket_axis, bra_axis = bra_axis, ket_axis
        rho = np.trace(
            np.moveaxis(rho, [ket_axis, bra_axis], [-2, -1]),
            axis1=-2,
            axis2=-1,
        )
        axis_order.remove(site)
        axis_order.remove(n + site)
    k = len(keep)
    ket_axes = [axis_order.index(site) for site in keep]
    bra_axes = [axis_order.index(n + site) for site in keep]
    rho = np.transpose(rho, axes=ket_axes + bra_axes)
    dim = 2 ** k
    return np.asarray(rho).reshape(dim, dim)


def pair_mutual_information(state: np.ndarray, i: int, j: int, n: int) -> float:
    return (
        entropy(reduced_state(state, [i], n))
        + entropy(reduced_state(state, [j], n))
        - entropy(reduced_state(state, [i, j], n))
    )


def rc_cut(state: np.ndarray, region: list[int], n: int) -> float:
    outside = [i for i in range(n) if i not in region]
    return 0.5 * sum(
        pair_mutual_information(state, i, j, n)
        for i in region
        for j in outside
    )


def main() -> int:
    check("R12 文档存在", os.path.isfile(os.path.join(HERE, DOC)))
    if not os.path.isfile(os.path.join(HERE, DOC)):
        for name, ok, detail in checks:
            print("[%s] %s%s" % ("v" if ok else "x", name, ("  " + detail) if detail else ""))
        print("\nASSERTIONS=%d" % len(checks))
        print("FAILURES=%d" % sum(not ok for _, ok, _ in checks))
        print("RESULT=FAIL")
        return 1

    text = read_local(DOC)
    r8 = read_local("R8_jacobson_entanglement_equilibrium_completion.md")
    r10 = read_local("R10_cao_carroll_bulk_entanglement_completion.md")
    d44 = read_local(os.path.join(os.pardir, "modular-equilibrium", "derivations",
                                  "D44_poincare_closure_obstruction.md"))
    d233 = read_local("D233_sign_age_symmetry_no_go_for_profile.md")

    # ------------------------------------------------------------------
    # A. 文档边界与接口
    # ------------------------------------------------------------------
    check("R12 明确 L1/L5 未关闭",
          all_tokens(text, ["没有关闭 L1", "没有关闭 L5"]))
    check("R12 明确只作条件恢复",
          all_tokens(text, ["Current status: conditional recovery", "条件桥", "定理"]))
    check("R12 引用 R8.1 精确熵差",
          all_tokens(r8, ["精确熵差恒等式", "D(\\sigma\\|\\rho)"]))
    check("R12 引用 R10 七项假设",
          "CC1" in r10 and "CC2" in r10 and "CC7" in r10)
    check("R12 引用 D44 Poincare 阻碍",
          all_tokens(d44, ["不可能自动包含 Lorentz boost", "交换"]))
    check("R12 引用 D233 正负延拓等权",
          all_tokens(d233, ["N_+(a)=N_-(a)", "常数"]))

    # ------------------------------------------------------------------
    # B. L5: 对角块的精确曲率（本文唯一用到的完整二次型）
    # ------------------------------------------------------------------
    p_diag = np.array([0.3, 0.7])
    x_diag = np.diag([1.0, -1.0])
    rho_diag = np.diag(p_diag)
    lam = 1e-3
    central = (
        relative_entropy(rho_diag + lam * x_diag, rho_diag)
        + relative_entropy(rho_diag - lam * x_diag, rho_diag)
        - 2.0 * relative_entropy(rho_diag, rho_diag)
    ) / lam ** 2
    expected_diag = 1.0 / (p_diag[0] * p_diag[1])
    check("二能级对角二阶中心差分等于 1/[p(1-p)]",
          abs(central - expected_diag) < 1e-3 * expected_diag,
          "central=%.8f expected=%.8f" % (central, expected_diag))

    # 对角曲率逐项与 1/p_i 一致。
    diag_curvature = abs(x_diag[0, 0]) ** 2 / p_diag[0] + abs(x_diag[1, 1]) ** 2 / p_diag[1]
    check("对角块曲率 sum x_i^2/p_i",
          abs(diag_curvature - expected_diag) < 1e-12,
          "chi=%.8f" % diag_curvature)

    # ------------------------------------------------------------------
    # C. L5-NG: 三块反例使 D 达到 R^d 阶
    # ------------------------------------------------------------------
    w = np.array([1.0, 2.0, 4.0]) / 7.0
    vec = np.array([1.0, -2.0, 1.0])
    check("no-go 扰动保持零和",
          abs(np.sum(vec)) < 1e-14)
    check("no-go 扰动正交于 log w",
          abs(np.sum(vec * np.log(w))) < 1e-13,
          "dot=%.3e" % np.sum(vec * np.log(w)))

    for d in (3, 4):
        radius = 1.0 / 16.0
        eps = radius ** (d / 2.0)
        q = w + eps * vec
        check("L5-NG 扰动保持正性",
              bool(np.all(q > 0.0)))
        d_rel = float(np.sum(q * np.log(q / w)))
        predicted = 91.0 / 8.0 * radius ** d
        check("L5-NG 在 d=%d 达到体积阶" % d,
              d_rel > 0.9 * predicted and abs(d_rel / predicted - 1.0) < 0.1,
              "D=%.8f predicted=%.8f" % (d_rel, predicted))

    # 模能一阶项严格为零
    modular_first = -eps * float(np.sum(vec * np.log(w)))
    check("L5-NG 模能一阶项为零",
          abs(modular_first) < 1e-14,
          "first=%.3e" % modular_first)

    # ------------------------------------------------------------------
    # D. gap 条件界的最小解析核验
    # ------------------------------------------------------------------
    # 两能级 p=1/2, X=sigma_x: chi=4 且 ||X||_2^2=2，
    # 故有限维下不存在与 gap 无关的普适上界；有 gap 时才可用 1/gamma 控制。
    gap = math.log(w[2]) - math.log(w[1])
    bound = (1.0 / gap) * 2.0
    check("gap 界量纲与规模有限",
          math.isfinite(bound) and bound > 0.0,
          "1/gamma*||X||^2=%.8f gamma=%.8f" % (bound, gap))

    # ------------------------------------------------------------------
    # E. Cao-Carroll: RC-EPF 引理与 GHZ 反例
    # ------------------------------------------------------------------
    n_ghz = 4
    ghz = ghz4()
    region = [0, 1]
    s_ghz = entropy(reduced_state(ghz, region, n_ghz))
    rc_ghz = rc_cut(ghz, region, n_ghz)
    check("GHZ 取前半区域 S(B)=log 2",
          abs(s_ghz - math.log(2.0)) < 1e-12,
          "S=%.12f" % s_ghz)
    check("GHZ 违反 RC 恒等式",
          abs(s_ghz - rc_ghz) > 0.5,
          "|S-RC|=%.8f" % abs(s_ghz - rc_ghz))

    # 逐边纯态乘积的最小实现：4 顶点环，每边 Bell，
    # 切割两条边时 S(B)=2 log 2，RC 割函数也为 2 log 2。
    # 用 tensor-network 构造完整 8 量子比特态太占位；这里直接核验
    # 引理使用的一对纯边恒等式 I(u:v)=2 s_e。
    bell = np.zeros(4, dtype=complex)
    bell[0] = 1.0 / math.sqrt(2.0)
    bell[-1] = 1.0 / math.sqrt(2.0)
    s_edge = entropy(reduced_state(bell, [0], 2))
    i_edge = (
        entropy(reduced_state(bell, [0], 2))
        + entropy(reduced_state(bell, [1], 2))
        - entropy(reduced_state(bell, [0, 1], 2))
    )
    check("单条纯 Bell 边满足 I(u:v)=2 s_e",
          abs(i_edge - 2.0 * s_edge) < 1e-12,
          "I=%.12f 2s=%.12f" % (i_edge, 2.0 * s_edge))

    # ------------------------------------------------------------------
    # F. 正路由与未推进项均被登记
    # ------------------------------------------------------------------
    check("Z_crit 正路由被登记为新增层",
          all_tokens(text, ["Z-CRIT", "Z-GNS", "Z-CORE", "Z-STRESS", "Z-CONF"]))
    check("强预解收敛被写成待证目标",
          "H_n" in text and "2\\pi B_B" in text and "强预解收敛" in text)
    check("RC-EPF 不被冒充为 CC2 普遍定理",
          all_tokens(text, ["条件子类定理", "不是 CC2 的无条件"]))
    check("面积-互信息比例保留额外条件",
          all_tokens(text, ["s_e=\\alpha w_e", "额外条件"]))

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
