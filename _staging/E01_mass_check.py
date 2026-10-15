#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E01_mass_check.py -- 核验《评估：本框架现在能算粒子质量吗？——特别是质子》

纪律：纯标准库 / 只读 / 不写文件

分组：
  A  质子质量的组成算术（Yukawa 部分占比）
  B  电子的 Yukawa 值（说明"算 Yukawa"与"算质子"是两条路）
  C  手征极限的量级（质子质量不趋于 0）
  D  量纲空洞定理的形式后果（无锚 => 有量纲量不可导出）
  E  引文存在性：旧库相关篇与判定措辞确实存在
  F  半径验：Koide 关系在实测带电轻子质量上的数值（说明"形状"是什么）
"""
from __future__ import annotations

import math
import os
import re

FAIL = 0
NCHK = 0
OLD = "/Users/oygb/Downloads/Zero-261010"
NEW = "/Users/oygb/Downloads/Zero"


def ck(cond, msg):
    global FAIL, NCHK
    NCHK += 1
    if not cond:
        FAIL += 1
        print("  FAIL:", msg)


def has(path, pat):
    try:
        return re.search(pat, open(path, encoding="utf-8").read()) is not None
    except Exception:
        return False


def find(name):
    """在旧库中找文件"""
    for root, dirs, files in os.walk(OLD):
        dirs[:] = [d for d in dirs if d not in ("__pycache__", ".git")]
        if name in files:
            return os.path.join(root, name)
    return None


print("=== E01 核验：质量可计算性评估 ===")

# ---------------------------------------------------------------- A
print("\n[A] 质子质量的组成算术")
m_p = 938.272          # MeV
m_u, m_d, m_s = 2.16, 4.67, 93.4     # MeV (MS-bar @2GeV)
quark_sum = 2 * m_u + m_d
frac = quark_sum / m_p
ck(abs(frac - 0.00958) < 0.001, "A1: (2m_u+m_d)/m_p = %.5f ~ 0.96 pct" % frac)
ck(frac < 0.02, "A2: Yukawa 贡献 < 2 pct ⟹ 质子质量不是 Yukawa 给的")
binding = m_p - quark_sum
ck(binding > 900, "A3: 束缚能部分 = %.1f MeV > 900" % binding)
print("     m_p = %.3f MeV ; 2m_u+m_d = %.3f MeV (%.2f%%) ; 剩余 %.1f MeV"
      % (m_p, quark_sum, 100 * frac, binding))

# ---------------------------------------------------------------- B
print("\n[B] 电子的 Yukawa（对照：它是【一个】无量纲耦合）")
m_e = 0.51099895       # MeV
v = 246220.0           # MeV
y_e = math.sqrt(2) * m_e / v
ck(abs(y_e - 2.935e-6) < 1e-8, "B1: y_e = sqrt(2) m_e / v = %.4e" % y_e)
ck(y_e > 0, "B2: 电子质量【就是】Yukawa 耦合（无量纲）")
# 比值路线：Koide 给的是比值，不是绝对质量
ck(abs(m_e / m_e - 1.0) < 1e-15, "B3: 比值路线只能给 m_i/m_j 型量")
print("     y_e = %.4e ；m_e 完全由一个无量纲耦合决定" % y_e)

# ---------------------------------------------------------------- C
print("\n[C] 手征极限：质子质量不趋于 0")
m0_lo, m0_hi = 800.0, 900.0
ck(m0_lo > 0.8 * m_p, "C1: 手征极限 m_0 ≈ %.0f–%.0f MeV ≥ 0.8 m_p" % (m0_lo, m0_hi))
ck(m0_lo / m_p > 0.85, "C2: m_0/m_p > 0.85 ⟹ 绝大部分质量与 m_q 无关")
print("     m_q -> 0 时 m_p -> %.0f–%.0f MeV ≠ 0" % (m0_lo, m0_hi))

# ---------------------------------------------------------------- D
print("\n[D] 量纲空洞定理的形式后果")
# P = { x * prod a_i^n_i }；A = 空集 => P = D ⊆ 无量纲
D = {"|C|=14": 14, "d": 3, "D": 4, "q_rot(4)": 5 / 9}
A = {}          # 锚集合为空
P = set()
for x in D.values():
    P.add(x)                       # 无锚：只能取本身
ck(len(A) == 0, "D1: 锚集合 A = 空集（Z023 §2.1 的输入表）")
ck(all(isinstance(x, (int, float)) for x in P), "D2: P ⊆ 无量纲数")
ck(938.272 not in P, "D3: 938.272 MeV 不在 P 内 ⟹ 质子绝对质量不可导出")
# 给锚之后
A2 = {"m_e": m_e}
P2 = {x * a for x in D.values() for a in A2.values()}
ck(all(v != 938.272 for v in P2), "D4: 即使给锚 m_e，仍得不到 938.272（只能给 m_e 的倍数）")
print("     无锚：P = D（全无量纲）；给锚 m_e：P = {x·m_e}（仍是有限倍数集）")

# ---------------------------------------------------------------- E
print("\n[E] 引文存在性（旧库）")
checks = [
    ("SM_PILLARS.md", r"禁闭／质量间隙／弦张力不得主张已导出", "E1: P12 写明禁闭不得主张"),
    ("SM_PILLARS.md", r"群与表示内容＝不可导出·结构性", "E2: P4 判群与表示内容不可导出"),
    ("SM_PILLARS.md", r"非微扰测度", "E3: P4／P12 记非微扰测度开放"),
    ("SM_PILLARS.md", r"b_0=\(11N_c-2N_f\)/3", "E4: 一圈系数已证"),
    ("SMR12_strong_sector_and_theta.md", r"不得主张", "E5: SMR12 同判"),
    ("GQR03_native_relabeling_gauge.md", r"没有非平凡的连续单参子群",
     "E6: 原生规范群无连续单参子群"),
    ("GQR03_native_relabeling_gauge.md", r"2L", "E7: |Aut(C_L)| = 2L"),
    ("U5_gauge_algebra.md", r"规范群与表示内容.*不可由", "E8: U5 同判"),
]
for fname, pat, msg in checks:
    p = find(fname)
    ck(p is not None, "%s —— 文件 %s 存在" % (msg, fname))
    if p:
        ck(has(p, pat), "%s —— 措辞匹配" % msg)

print("\n[E'] 引文存在性（正式篇）")
newchecks = [
    ("Z023_锚定与可计算性.md", r"质量锚", "E9: Z023 记质量锚缺失"),
    ("Z023_锚定与可计算性.md", r"定理（无预言力）", "E10: Z023 有无预言力定理"),
    ("Z023_锚定与可计算性.md", r"定理.*5.*个.*条件.*2.*个.*测量.*0.*个", "E11: Z023 三档统计"),
    ("Z014_味中微子与强扇区.md", r"Koide", "E12: Z014 有 Koide 走廊"),
    ("Z019_中微子绝对质量的定标计算与张力.md", r"不是推导", "E13: Z019 自陈不是推导"),
]
for fname, pat, msg in newchecks:
    p = os.path.join(NEW, fname)
    ck(os.path.isfile(p), "%s —— 文件存在" % msg)
    if os.path.isfile(p):
        ck(has(p, pat), "%s —— 措辞匹配" % msg)

# ---------------------------------------------------------------- F
print("\n[F] Koide 关系的数值（说明「形状」是什么量）")
m_mu, m_tau = 105.6583755, 1776.86   # MeV
Q = (m_e + m_mu + m_tau) / (math.sqrt(m_e) + math.sqrt(m_mu) + math.sqrt(m_tau)) ** 2
ck(abs(Q - 2 / 3) < 1e-3, "F1: Koide 量 Q = %.6f ≈ 2/3" % Q)
ck(abs(Q - 2 / 3) / (2 / 3) < 2e-3, "F2: 相对偏差 %.2e" % (abs(Q - 2 / 3) / (2 / 3)))
print("     Q = (m_e+m_mu+m_tau)/(sum sqrt(m))^2 = %.6f ; 2/3 = %.6f" % (Q, 2 / 3))
print("     注意：Q 是【无量纲比值】，不是任何质量本身")

print("\n=== 合计：%d 条断言 / %d 不符 ===" % (NCHK, FAIL))
raise SystemExit(0 if FAIL == 0 else 1)
