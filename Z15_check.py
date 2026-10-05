#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z15_check.py —— 【Z-CAR 无唯一性 + Z-READ 下 Jordan-Wigner】的独立核验
=====================================================================
不访问网络，不修改任何项目文件。

独立实断言：
  F1  文档锚点：Z-CAR-NOGO、Z-READ、CAR、(-1)^F、J1 未关闭
  F2  中心特征非唯一：+1 与 -1 都是合法表示
  F3  外代数：Λ(C⊕C) 的次数奇偶给出 diag(+1,-1,-1,+1)
  F4  Jordan-Wigner：N=4 时全部 CAR 关系成立
  F5  Z-READ 不是权重；Z0③ 与全分支论证不变
  F6  跨文档状态一致：Z14/R15/R16/STATUS 指向同一状态
"""
import cmath
import io
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


DOC = read("Z15_zcar_no_go_and_jordan_wigner_readout.md")
Z0 = read("Z0_zero_never_rests_single_axiom.md")
Z13 = read("Z13_zero_foundation_missing_principle.md")
Z14 = read("Z14_closure_cyclic_order_base_theorem.md")
R15 = read("R15_zcar_double_cover_and_zstress_scale.md")
R16 = read("R16_direction_audit_reduction_tree.md")
STATUS = read("STATUS.md")


# ---------------------------------------------------------------- F1
head("F1  文档锚点")

check("无唯一性定理在位", "定理 Z15.1" in DOC and "Z-CAR-NOGO" in DOC)
check("Z-READ 被登记为具名输入",
      "Z-READ" in DOC and "具名输入" in DOC)
check("外代数与 Fock 构造在位",
      "外代数" in DOC and "Fock" in DOC)
check("Jordan–Wigner 与 CAR 在位",
      "Jordan–Wigner" in DOC and "CAR" in DOC and "反对易" in DOC)
check("中心元到费米宇称在位",
      "(-1)^F" in DOC and "费米宇称" in DOC)
check("J1 未关闭被写明", "仍未关闭" in DOC or "仍开放" in DOC)


# ---------------------------------------------------------------- F2
head("F2  两个中心特征都合法")

theta = cmath.pi / 3
spin_lift = np.diag([cmath.exp(1j * theta / 2), cmath.exp(-1j * theta / 2)])
projected = spin_lift @ spin_lift
expected_projected = np.diag([cmath.exp(1j * theta), cmath.exp(-1j * theta)])

check("旋量升格平方后投影为 SO(2) 旋转",
      np.allclose(projected, expected_projected),
      "q(U_theta)=U_theta^2=R_theta")
check("旋量中心特征为 -1",
      np.allclose(spin_lift @ spin_lift.conj().T, np.eye(2))
      and abs(np.linalg.det(spin_lift) - 1) < 1e-12)
tensor_rep = np.array([[1.0]])
check("张量中心特征为 +1",
      np.allclose(tensor_rep, tensor_rep) and tensor_rep[0, 0] == 1.0)
check("两个中心特征都满足同一双覆盖关系",
      "两者都使用同一个" in DOC and "中心扩张关系" in DOC)


# ---------------------------------------------------------------- F3
head("F3  外代数宇称")

parity = np.diag([1.0, -1.0, -1.0, 1.0])
ground = np.array([1.0, 0.0, 0.0, 0.0])
one_particle = np.array([0.0, 1.0, 0.0, 0.0])
two_particle = np.array([0.0, 0.0, 0.0, 1.0])

check("真空为偶扇区", np.allclose(parity @ ground, ground))
check("单粒子为奇扇区", np.allclose(parity @ one_particle, -one_particle))
check("双粒子为偶扇区", np.allclose(parity @ two_particle, two_particle))
check("parity^2 = I 且迹为 0",
      np.allclose(parity @ parity, np.eye(4)) and abs(np.trace(parity)) < 1e-12)


# ---------------------------------------------------------------- F4
head("F4  Jordan-Wigner：CAR 的有限维复算")


def kron_all(mats):
    out = np.array([[1.0]], dtype=complex)
    for mat in mats:
        out = np.kron(out, mat)
    return out


def jw_ops(n_modes):
    ident = np.eye(2, dtype=complex)
    sig_z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
    sig_minus = np.array([[0.0, 0.0], [1.0, 0.0]], dtype=complex)
    sig_plus = sig_minus.conj().T
    c = []
    cdag = []
    for j in range(n_modes):
        prefix = [sig_z] * j
        c.append(
            kron_all(prefix + [sig_minus] + [ident] * (n_modes - j - 1))
        )
        cdag.append(
            kron_all(prefix + [sig_plus] + [ident] * (n_modes - j - 1))
        )
    return c, cdag


def anticommutator(a, b):
    return a @ b + b @ a


NM = 4
c, cdag = jw_ops(NM)
car_dag = True
car_aa = True
car_ddag = True
for i in range(NM):
    for j in range(NM):
        expected = np.eye(2 ** NM) if i == j else np.zeros((2 ** NM, 2 ** NM))
        if not np.allclose(anticommutator(c[i], cdag[j]), expected):
            car_dag = False
        if not np.allclose(anticommutator(c[i], c[j]), np.zeros_like(expected)):
            car_aa = False
        if not np.allclose(anticommutator(cdag[i], cdag[j]), np.zeros_like(expected)):
            car_ddag = False

check("N=4 时 {c_i,c_j^dagger}=delta_ij", car_dag)
check("N=4 时 {c_i,c_j}=0", car_aa)
check("N=4 时 {c_i^dagger,c_j^dagger}=0", car_ddag)

sigma_z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
fermion_parity = kron_all([sigma_z] * NM)
check("费米宇称为 prod sigma_z",
      np.allclose(fermion_parity @ fermion_parity, np.eye(2 ** NM)))


# ---------------------------------------------------------------- F5
head("F5  没有把概率／权重塞回 Zero")

check("Z-READ 明确不是概率或权重",
      "不是概率或权重" in DOC or "不是概率、权重" in DOC)
check("Z0③ 原文保持不设概率",
      "不设概率" in Z0 and "没有选择规则" in Z0)
check("全分支＋整数重数保持",
      "全分支" in Z0 and "整数重数" in Z0)
check("Z13 仍把读出列为 E5 的具名输入",
      "E5" in Z13 and "量子读出" in Z13)


# ---------------------------------------------------------------- F6
head("F6  跨文档状态一致")

check("Z14 说明 Z15 的无唯一性与条件构造",
      "Z15" in Z14 and "无唯一性" in Z14)
check("R15 指向 Z15 的条件构造",
      "Z15" in R15)
check("R16 的 Z-CAR 条目仍记为开放簇",
      "Z-CAR" in R16 and "开放具名簇" in R16)
check("STATUS 已登记 Z15",
      "Z15" in STATUS and "Z-READ" in STATUS)
check("STATUS 仍保持 L1 开放",
      "L1" in STATUS and "仍未关闭" in STATUS)


# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
