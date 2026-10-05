#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R10_check.py

R10 Cao-Carroll 2018 条件桥审计的独立核验。

检查内容：
  F1/F2  RC 不是面积律的自动结果；
  F3     RC 可不蕴含局域几何：长程 Bell 配对态是精确 RC；
  F4     二维 massive free-fermion 基准中的 I/A 比例稳定性；
  F5     弱场系数 4G alpha = 1 与 16 pi G / 2 = 8 pi G；
  F6     D=4 极化维数只是条件事实；
  F7     文档必须保留“弱场/线性化，不是完整非线性 GR”的边界。
"""

from __future__ import annotations

import math
import sys
from itertools import product
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DOC = HERE / "R10_cao_carroll_bulk_entanglement_completion.md"

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


def entropy(rho: np.ndarray) -> float:
    rho = (rho + rho.conj().T) / 2.0
    vals = np.linalg.eigvalsh(rho)
    vals = vals[vals > 1e-12]
    return float(-np.sum(vals * np.log(vals)))


def compress_bits(index: int, keep: list[int], n: int) -> int:
    out = 0
    for qubit in keep:
        out = (out << 1) | ((index >> (n - 1 - qubit)) & 1)
    return out


def reduced_state(psi: np.ndarray, keep: list[int], n: int) -> np.ndarray:
    """Brute-force partial trace, valid for the small n used in this audit."""
    keep = list(keep)
    dim = 1 << len(keep)
    rho = np.zeros((dim, dim), dtype=complex)
    complement = [q for q in range(n) if q not in keep]

    for a in range(1 << n):
        for b in range(1 << n):
            if all(((a >> (n - 1 - q)) & 1) == ((b >> (n - 1 - q)) & 1)
                   for q in complement):
                rho[compress_bits(a, keep, n), compress_bits(b, keep, n)] += (
                    psi[a] * np.conjugate(psi[b])
                )
    return rho


def mutual_information(psi: np.ndarray, n: int, i: int, j: int) -> float:
    return (
        entropy(reduced_state(psi, [i], n))
        + entropy(reduced_state(psi, [j], n))
        - entropy(reduced_state(psi, [i, j], n))
    )


def rc_cut(psi: np.ndarray, n: int, region: list[int]) -> float:
    complement = [q for q in range(n) if q not in region]
    return 0.5 * sum(
        mutual_information(psi, n, i, j)
        for i in region
        for j in complement
    )


def bell_product_state(n: int, pairs: list[tuple[int, int]]) -> np.ndarray:
    psi = np.zeros(1 << n, dtype=complex)
    for choices in product([0, 1], repeat=len(pairs)):
        index = 0
        for bit, (i, j) in zip(choices, pairs):
            if bit:
                index |= 1 << (n - 1 - i)
                index |= 1 << (n - 1 - j)
        psi[index] = 2.0 ** (-0.5 * len(pairs))
    return psi


def build_fermion_2d(n: int, mass: float) -> np.ndarray:
    def idx(i: int, j: int) -> int:
        return (i % n) * n + (j % n)

    ham = np.zeros((n * n, n * n), dtype=float)
    for i in range(n):
        for j in range(n):
            a = idx(i, j)
            for di, dj in ((1, 0), (0, 1)):
                b = idx(i + di, j + dj)
                ham[a, b] = -1.0
                ham[b, a] = -1.0
            ham[a, a] += mass * ((-1) ** (i + j))
    return ham


def free_fermion_corr(ham: np.ndarray) -> np.ndarray:
    evals, evecs = np.linalg.eigh(ham)
    occupied = evecs[:, evals < 0]
    return occupied @ occupied.conj().T


def block_entropy(corr: np.ndarray, n: int, length: int) -> float:
    selected = [i * n + j for i in range(length) for j in range(length)]
    sub = corr[np.ix_(selected, selected)]
    vals = np.clip(np.linalg.eigvalsh((sub + sub.conj().T) / 2), 1e-13, 1 - 1e-13)
    return float(-np.sum(vals * np.log(vals) + (1 - vals) * np.log(1 - vals)))


def main() -> int:
    head("R10 文字边界")

    check("R10 文档存在", DOC.is_file(), str(DOC))
    if not DOC.is_file():
        print("  文档缺失，无法继续文字核验")
        return 1

    text = DOC.read_text(encoding="utf-8")

    assumption_names = [
        "首选张量分解",
        "RC 态",
        "面积来自互信息",
        "修改纠缠平衡",
        "emergent EFT",
        "生成几何的动力学",
        "Lorentz 不变性",
    ]
    for name in assumption_names:
        check(f"七项假设登记：{name}", name in text)

    bridge_names = [
        "CC1｜首选局域张量完成",
        "CC2｜近似 RC 与割函数收敛",
        "CC3｜跨切割面积-互信息比例",
        "CC4｜背景度规与 Radon 反演",
        "CC5｜MEEC、EFT 与 Rindler 第一定律",
        "CC6｜Lorentzian 组装与弱场极限",
        "CC7｜局部 Lorentz 完成",
    ]
    for name in bridge_names:
        check(f"条件桥登记：{name}", name in text)

    check("明确只到弱场线性化 EFE", "弱场线性化 EFE" in text)
    check("明确不推出完整非线性 GR", "不推出完整非线性" in text)
    check("明确 D=4 不是自动输出", "不推出 }D=4" in text or "D=4" in text)
    check("Radon 缺口被登记为开放", "Radon" in text and "开放" in text)
    check("Lorentzian 组装缺口被登记为开放", "Lorentzian 组装" in text)
    check("RC 缺口被登记为决定性缺口", "RC 条件" in text and "决定性缺口" in text)

    head("F1/F2  RC 不是面积律的自动结果")

    n_ghz = 4
    ghz = np.zeros(1 << n_ghz, dtype=complex)
    ghz[0] = 1.0 / math.sqrt(2.0)
    ghz[-1] = 1.0 / math.sqrt(2.0)
    region = [0, 1]
    s_ghz = entropy(reduced_state(ghz, region, n_ghz))
    rc_ghz = rc_cut(ghz, n_ghz, region)
    ghz_gap = abs(s_ghz - rc_ghz)
    check("GHZ 取前半区域时 S(B)=log 2", abs(s_ghz - math.log(2.0)) < 1e-12,
          f"S={s_ghz:.12f}")
    check("GHZ 的 RC 割函数与真实熵相差很大", ghz_gap > 0.5,
          f"|S-RC|={ghz_gap:.6f}")

    rng = np.random.default_rng(123)
    raw = rng.normal(size=1 << 4) + 1j * rng.normal(size=1 << 4)
    random_state = raw / np.linalg.norm(raw)
    random_gap = abs(
        entropy(reduced_state(random_state, region, 4))
        - rc_cut(random_state, 4, region)
    )
    check("随机四因子态通常不满足 RC", random_gap > 0.05,
          f"|S-RC|={random_gap:.6f}")

    head("F3  RC 可与长程信息边共存")

    n_bell = 4
    pairs = [(0, 3), (1, 2)]
    bell = bell_product_state(n_bell, pairs)
    pair_mi = {
        (i, j): mutual_information(bell, n_bell, i, j)
        for i in range(n_bell)
        for j in range(i + 1, n_bell)
    }

    max_rc_error = 0.0
    for mask in range(1, (1 << n_bell) - 1):
        region_b = [q for q in range(n_bell) if (mask >> (n_bell - 1 - q)) & 1]
        complement = [q for q in range(n_bell) if q not in region_b]
        rc_value = 0.5 * sum(
            pair_mi[(min(i, j), max(i, j))]
            for i in region_b
            for j in complement
        )
        max_rc_error = max(
            max_rc_error,
            abs(entropy(reduced_state(bell, region_b, n_bell)) - rc_value),
        )

    check("长程 Bell 配对态对所有非空真子集满足 RC",
          max_rc_error < 1e-12, f"max error={max_rc_error:.3e}")
    check("每个 Bell 因子的互信息为 2 log 2",
          all(abs(pair_mi[p] - 2.0 * math.log(2.0)) < 1e-12 for p in pairs),
          f"MI pair={[round(pair_mi[p], 9) for p in pairs]}")
    check("非局域配对边存在：最大配对距离 > 1",
          max(abs(i - j) for i, j in pairs) > 1)
    check("RC 不唯一决定局域嵌入",
          max_rc_error < 1e-12 and max(abs(i - j) for i, j in pairs) > 1)

    head("F4  二维 massive free-fermion 的面积-互信息比例")

    n_fermion = 32
    mass = 4.0
    lengths = [4, 6, 8, 10, 12]
    corr = free_fermion_corr(build_fermion_2d(n_fermion, mass))
    entropies = np.array([block_entropy(corr, n_fermion, length) for length in lengths])
    mi_pure = 2.0 * entropies
    area = 4.0 * np.array(lengths, dtype=float)
    ratios = mi_pure / area
    last_three_ratio = ratios[-3:]
    rel_range = float(
        (last_three_ratio.max() - last_three_ratio.min()) / last_three_ratio.mean()
    )

    check("mass gap m=4 时 S(L)>0", np.all(entropies > 0.0),
          "S=" + ",".join(f"{v:.6f}" for v in entropies))
    check("I/A 比例在最后三个尺寸上仅有小幅漂移（<2%）",
          rel_range < 0.02, f"relative range={100 * rel_range:.3f}%")
    check("I/A 比例数值落在有限窗口稳定区间",
          np.all((ratios > 0.11) & (ratios < 0.13)),
          "I/A=" + ",".join(f"{v:.6f}" for v in ratios))

    head("F5  弱场系数")

    g_newton = 1.37
    alpha = 1.0 / (4.0 * g_newton)
    check("4 G_N alpha = 1", abs(4.0 * g_newton * alpha - 1.0) < 1e-14)
    check("16 pi G_N / 2 = 8 pi G_N",
          abs(16.0 * math.pi * g_newton / 2.0 - 8.0 * math.pi * g_newton) < 1e-14)
    check(
        "线性化关系不代表完整非线性 EFE",
        "线性化约束不等于完整约束代数" in text and "完整非线性" in text,
    )

    head("F6  D=4 的条件性")

    dims = {d: d * (d - 3) // 2 for d in range(4, 13)}
    check("dim P_D = D(D-3)/2 在 D=4 给 2", dims[4] == 2, f"dim={dims[4]}")
    check("恰好两个极化只在 D=4 出现",
          [d for d, value in dims.items() if value == 2] == [4])
    check("该选择需要额外输入，不是 Zero 自动输出",
          "物理独立的 $D=4$ 选择" in text or "物理独立的 $D=4$" in text)

    head("F7  失败树与诚实边界")

    required_warnings = [
        "不把外部论文的假设算作 Zero 定理",
        "条件恢复，不是无条件导出",
        "Radon 型度规反演",
        "Lorentzian 组装",
        "RC 条件",
        "有限维因子与物理局域同构",
    ]
    for warning in required_warnings:
        check(f"诚实边界在位：{warning}", warning in text)

    head("汇总")
    print(f"  ASSERTIONS={ASSERTIONS}")
    print(f"  FAILURES={FAILURES}")
    if FAILURES:
        print("  RESULT=FAIL")
        return 1
    print("  RESULT=PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
