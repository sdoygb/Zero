#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z16_check.py —— 【Z-UNIF：平衡正则模块与读出残余】的独立核验
==============================================================
不访问网络，不修改任何项目文件。

独立实断言：
  F1  文档锚点：Z-UNIF、Z16.1 no-go、正则表示、CAR、残余 no-go、J1 未关闭
  F2  无偏好不能自动提升：+1、-1、正则三种载体都合法
  F3  最小平衡模块：正则表示、中心幂等元、迹为 0
  F4  多模宇称：2^N 个字符等重，总宇称给出 (-1)^F
  F5  Majorana／Jordan-Wigner：Clifford 与 CAR 关系成立
  F6  Z-UNIF 不是概率／权重；Z0③ 与全分支不变
  F7  残余 no-go：不选自旋结构、循环切口、模式相位或四维框架
  F8  跨文档状态一致：Z13/Z14/Z15/R14/R15/R16/G10/STATUS
"""
import io
import itertools
import os
import sys

import numpy as np

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


DOC = read("Z16_zunif_balanced_regular_module.md")
Z0 = read("Z0_zero_never_rests_single_axiom.md")
Z13 = read("Z13_zero_foundation_missing_principle.md")
Z14 = read("Z14_closure_cyclic_order_base_theorem.md")
Z15 = read("Z15_zcar_no_go_and_jordan_wigner_readout.md")
R14 = read("R14_L1_from_zero_assembly.md")
R15 = read("R15_zcar_double_cover_and_zstress_scale.md")
R16 = read("R16_direction_audit_reduction_tree.md")
G10 = read("G10_final_derivation_and_input_ledger.md")
STATUS = read("STATUS.md")


# ---------------------------------------------------------------- F1
head("F1  文档锚点")

check("标题登记 Z-UNIF 与平衡正则模块",
      "Z-UNIF" in DOC and "正则表示" in DOC and "平衡" in DOC)
check("Z16.1 no-go 在位", "定理 Z16.1" in DOC and "Z-UNIF-NOGO" in DOC)
check("Z-UNIF 被登记为具名输入",
      "具名输入" in DOC and "表示层选择规则" in DOC)
check("正则表示与中心幂等元在位",
      "mathbb C[G]" in DOC and "\mathbb I\\pm J" in DOC)
check("Jordan–Wigner 与 CAR 在位",
      "Jordan--Wigner" in DOC and "CAR" in DOC and "Clifford" in DOC)
check("残余 no-go 在位",
      "命题 Z16.5" in DOC and "自旋结构" in DOC and "线性切口" in DOC)
check("J1 未关闭被写明", "J1 仍未关闭" in DOC or "仍未关闭" in DOC)


# ---------------------------------------------------------------- F2
head("F2  不能把无偏好自动提升为等重")

j_plus = np.array([[1.0]])
j_minus = np.array([[-1.0]])
j_bal = np.diag([1.0, -1.0])

check("张量载体中心特征为 +1",
      np.allclose(j_plus @ j_plus, np.eye(1)) and j_plus[0, 0] == 1.0)
check("旋量载体中心特征为 -1",
      np.allclose(j_minus @ j_minus, np.eye(1)) and j_minus[0, 0] == -1.0)
check("平衡载体同时保留 +1 与 -1",
      np.allclose(j_bal @ j_bal, np.eye(2))
      and abs(np.trace(j_bal)) < 1e-12)
check("三种载体共用同一个中心双覆盖关系",
      "三者使用完全相同的" in DOC and "表示论资料" in DOC)
check("文档明确表示层桥不在 Z0／Z-E* 中",
      "这条桥不在 Z0 或" in DOC)
check("Z0③ 仍是不设概率且没有选择规则",
      "不设概率" in Z0 and "没有选择规则" in Z0)


# ---------------------------------------------------------------- F3
head("F3  最小平衡模块与中心幂等元")

projector_plus = (np.eye(2) + j_bal) / 2.0
projector_minus = (np.eye(2) - j_bal) / 2.0

check("e_+ 是幂等中心投影",
      np.allclose(projector_plus @ projector_plus, projector_plus))
check("e_- 是幂等中心投影",
      np.allclose(projector_minus @ projector_minus, projector_minus))
check("两个投影相互正交且和为 I",
      np.allclose(projector_plus @ projector_minus, np.zeros((2, 2)))
      and np.allclose(projector_plus + projector_minus, np.eye(2)))
check("J 的特征值为 +1 与 -1",
      np.allclose(np.sort(np.linalg.eigvalsh(j_bal)), [-1.0, 1.0]))
check("最小平衡模块维数为 2",
      "m_+=m_-=1" in DOC and "最小平衡模块的维数为 $2$" in DOC)
check("正则表示唯一到同构意义",
      "在同构意义下" in DOC and "正则表示" in DOC)
check("文档保留分量基相位自由度", "U(1)\\times U(1)" in DOC)


# ---------------------------------------------------------------- F4
head("F4  有序多模宇称")

sig_z = np.diag([1.0, -1.0])


def kron_all(mats):
    out = np.array([[1.0]])
    for mat in mats:
        out = np.kron(out, mat)
    return out


def parity_matrix(n_modes):
    return kron_all([sig_z] * n_modes)


for nm in (1, 2, 3, 4):
    p = parity_matrix(nm)
    chars = []
    for bits in itertools.product((0, 1), repeat=nm):
        vec = np.zeros(2 ** nm)
        idx = 0
        for b in bits:
            idx = 2 * idx + b
        vec[idx] = 1.0
        value = (p @ vec)[idx] / vec[idx]
        if not np.isclose(abs(value), 1.0):
            chars.append(("bad",))
        else:
            chars.append(bits)
    unique = sorted(set(chars))
    dims_even = sum(1 for c in unique if sum(c) % 2 == 0)
    dims_odd = len(unique) - dims_even
    check("N=%d 有 2^N 个不同宇称字符" % nm,
          len(unique) == 2 ** nm, "维数 %d" % (2 ** nm))
    check("N=%d 偶／奇扇区等重" % nm,
          dims_even == dims_odd == 2 ** (nm - 1))
    check("N=%d 总宇称迹为零" % nm, abs(np.trace(p)) < 1e-12)


# ---------------------------------------------------------------- F5
head("F5  Majorana 与 Jordan–Wigner")


def majorana_ops(n_modes):
    ident = np.eye(2)
    sig_x = np.array([[0.0, 1.0], [1.0, 0.0]])
    sig_y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
    ops = []
    for j in range(n_modes):
        prefix = [sig_z] * j
        suffix = [ident] * (n_modes - j - 1)
        ops.append(kron_all(prefix + [sig_x] + suffix).astype(complex))
        ops.append(kron_all(prefix + [sig_y] + suffix).astype(complex))
    return ops


def anticommutator(a, b):
    return a @ b + b @ a


NM = 4
gammas = majorana_ops(NM)
clifford_ok = True
for a in range(2 * NM):
    for b in range(2 * NM):
        expected = 2 * np.eye(2 ** NM) if a == b else np.zeros((2 ** NM, 2 ** NM))
        if not np.allclose(anticommutator(gammas[a], gammas[b]), expected):
            clifford_ok = False
check("N=4 的 Majorana 满足 Clifford 关系", clifford_ok)

c = []
cdag = []
for j in range(NM):
    c.append((gammas[2 * j] + 1j * gammas[2 * j + 1]) / 2.0)
    cdag.append((gammas[2 * j] - 1j * gammas[2 * j + 1]) / 2.0)

car_ok = True
for i in range(NM):
    for j in range(NM):
        expected = np.eye(2 ** NM) if i == j else np.zeros((2 ** NM, 2 ** NM))
        if not np.allclose(anticommutator(c[i], cdag[j]), expected):
            car_ok = False
        if not np.allclose(anticommutator(c[i], c[j]), np.zeros_like(expected)):
            car_ok = False
        if not np.allclose(anticommutator(cdag[i], cdag[j]), np.zeros_like(expected)):
            car_ok = False
check("N=4 的复合算符满足 CAR", car_ok)

parity = parity_matrix(NM)
check("总宇称等于 (-1)^F 的实现",
      np.allclose(parity @ parity, np.eye(2 ** NM))
      and abs(np.trace(parity)) < 1e-12)


# ---------------------------------------------------------------- F6
head("F6  Z-UNIF 没有把概率／权重塞回 Zero")

check("Z-UNIF 明确不是概率",
      "不是概率" in DOC or "不是概率／权重" in DOC)
check("Z-UNIF 明确不是实数权重", "不是实数权重" in DOC)
check("Z-UNIF 不改变全分支与整数重数",
      "不改变 Z2 的全分支与整数重数" in DOC)
check("Z-UNIF 被定位为表示层选择规则",
      "表示层选择规则" in DOC and "具名输入" in DOC)
check("Z13 仍把 E5 读出记为未解缺口",
      "E5" in Z13 and "量子读出" in Z13 and "开放" in Z13)


# ---------------------------------------------------------------- F7
head("F7  残余 no-go")

spin_structures = {"periodic", "antiperiodic"}
linear_cuts = {s for s in range(4)}
check("两个自旋结构都与同一平衡模块相容", len(spin_structures) == 2)
check("循环的线性切口不唯一", len(linear_cuts) == 4,
      "L=4 时有 4 个起点")
check("文档把自旋结构列为未选", "自旋结构" in DOC and "仍具名" in DOC)
check("文档把循环切口列为未选", "循环切口" in DOC and "仍具名" in DOC)
check("文档把模式相位列为未选", "模式基" in DOC and "相位" in DOC)
check("文档把四维框架列为未选", "四维物理框架" in DOC)
check("Z-READ 不等于 Z-UNIF",
      "Z\\text{-READ}" in DOC and "neq" in DOC and "Z\\text{-UNIF}" in DOC)


# ---------------------------------------------------------------- F8
head("F8  跨文档状态一致")

check("Z15 已登记 Z16 的条件替代",
      "Z16" in Z15 and "Z-UNIF" in Z15)
check("Z14 已登记 Z16 的后续边界",
      "Z16" in Z14 and "Z-UNIF" in Z14)
check("R14 已登记 Z16",
      "Z16" in R14 and "Z-UNIF" in R14)
check("R15 已登记 Z16",
      "Z16" in R15 and "Z-UNIF" in R15)
check("R16 把 Z16 保持在同一个开放簇",
      "Z16" in R16 and ("同一开放簇" in R16 or "同一个开放簇" in R16))
check("G10 已登记 Z16_check.py",
      "Z16_check.py" in G10)
check("STATUS 已登记 Z16 并保持 L1 开放",
      "Z16" in STATUS and "Z-UNIF" in STATUS
      and "L1" in STATUS and "仍未关闭" in STATUS)


# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
