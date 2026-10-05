#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R34_check.py -- 核验 P2 探针结果与"有限维不可能"定理的文档表述。

独立复算：so(1,3) 关系、K 的非 Hermitian 性、无实 Hermitian 解、boost 指数不酉、
K_omega 的 Hermitian 性与跨度 log50、so(4)/so(1,3) 的 Killing 形式符号。
"""

from __future__ import annotations

import io
import json
import os
import sys
from math import log

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


def head(title: str) -> None:
    print("")
    print("=" * 72)
    print(title)
    print("=" * 72)


def check(name: str, cond: bool, detail: str = "") -> None:
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def read(name: str) -> str:
    with io.open(os.path.join(HERE, name), encoding="utf-8") as handle:
        return handle.read()


DOC = read("R34_finite_dimensional_boost_obstruction.md")
R33 = read("R33_action_phase_match_project.md")
R19 = read("R19_L1_upstream_probability_phase_and_missing_boost.md")
R12 = read("R12_zero_native_gap_filling.md")
R13 = read("R13_L1_strong_resolvent_attempt.md")
G62 = read("G62_quantum_sector_from_GNS_modular_flow_gleason.md")
G72 = read("G72_kappa1_from_the_ledger.md")
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")
RES = json.load(io.open(os.path.join(HERE, "R34_mixing_generator_results.json"), encoding="utf-8"))

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
P = [sx, sy, sz]
EPSD = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1,
        (1, 0, 2): -1, (2, 1, 0): -1, (0, 2, 1): -1}


def comm(A, B):
    return A @ B - B @ A


def e(i, j, k):
    return EPSD.get((i, j, k), 0)


J = [s / 2 for s in P]
K = [1j * s / 2 for s in P]

# ======================================================================
head("F1  文档结构：定理、问法、目标重述")

check("标题与性质为探针结果＋不可能定理",
      "P2 探针结果" in DOC
      and "有限维**不可能**承载 boost" in DOC
      and "探针结果＋有限维不可能定理＋目标重述" in DOC)
check("四个问 Q1-Q4 与定理 R34.1 在位",
      all(("**Q%d**" % i) in DOC for i in (1, 2, 3, 4))
      and "定理 R34.1" in DOC)
check("目标重述 P2′／指纹 P2″ 在位",
      "P2′（重述）" in DOC and "P2″（可测指纹，新增）" in DOC)
check("没有写成已证极限或已证几何模流",
      "没有证明细化极限含 type III 因子" in DOC
      and "没有证明 $2\\pi$ 归一化" in DOC
      and "没有由四维标签推出 Lorentz" in DOC)

# ======================================================================
head("F2  独立复算：so(1,3) 关系与 Hermitian 性")

check("Q1: [J_i,J_j] = i eps J_k（旋转关系）",
      all(np.allclose(comm(J[i], J[j]), 1j * sum(e(i, j, k) * J[k] for k in range(3)))
          for i in range(3) for j in range(3)))
check("Q1: [J_i,K_j] = i eps K_k（混合关系）",
      all(np.allclose(comm(J[i], K[j]), 1j * sum(e(i, j, k) * K[k] for k in range(3)))
          for i in range(3) for j in range(3)))
check("Q1: [K_i,K_j] = -i eps J_k（boost-boost）",
      all(np.allclose(comm(K[i], K[j]), -1j * sum(e(i, j, k) * J[k] for k in range(3)))
          for i in range(3) for j in range(3)))
check("Q2: J Hermitian，K 非 Hermitian",
      all(np.allclose(J[i].conj().T, J[i]) for i in range(3))
      and all(not np.allclose(K[i].conj().T, K[i]) for i in range(3)))
check("Q2b: 不存在实系数 Hermitian boost（2c^2=-1/2 无实解）",
      not any(abs(1j * c * 0 + 2 * c * c + 0.5) < 1e-9 for c in np.linspace(-3, 3, 6001) if abs(c) > 1e-12))
check("Q2c: boost 指数不酉、旋转指数酉",
      not np.allclose(np.array([[np.cosh(0.5), -np.sinh(0.5)], [-np.sinh(0.5), np.cosh(0.5)]]).conj().T
                      @ np.array([[np.cosh(0.5), -np.sinh(0.5)], [-np.sinh(0.5), np.cosh(0.5)]]), np.eye(2))
      and np.allclose((lambda R: R.conj().T @ R)(
          np.array([[np.cos(0.5), -1j * np.sin(0.5)], [-1j * np.sin(0.5), np.cos(0.5)]])), np.eye(2)))
check("Q2c: boost 指数范数随 theta 增长（无界）",
      RES["Q3_boost_norm_growth"][-1] > 10 * RES["Q3_boost_norm_growth"][0],
      "norms=%s" % [round(x, 3) for x in RES["Q3_boost_norm_growth"]])

# ======================================================================
head("F3  独立复算：原生 K_omega 与紧性")

counts = np.array([1.0, 2.0, 10.0, 50.0])
K_omega = np.diag(-np.log(counts / counts.sum()))
check("K_omega Hermitian 且谱有限（4 个点）",
      np.allclose(K_omega.conj().T, K_omega) and len(np.diag(K_omega)) == 4)
check("K_omega 跨度 = log 50（与 G72 §1 一致）",
      abs((np.diag(K_omega).max() - np.diag(K_omega).min()) - log(50)) < 1e-9,
      "span=%.6f" % (np.diag(K_omega).max() - np.diag(K_omega).min()))
check("有限维 ⇒ 模流闭包紧（准周期，pair-periods 有限）",
      RES["Q3_modular_flow_compact"] and RES["Q3_pair_periods_count"] == 6)
check("G72 §1 的 3.912 原文在位",
      "3.912" in G72)

# ======================================================================
head("F4  独立复算：Killing 形式的紧／非紧分界")


def so_basis(sig):
    eta = np.diag(sig).astype(float)
    basis = []
    for a in range(4):
        for b in range(a + 1, 4):
            M = np.zeros((4, 4))
            M[a, b], M[b, a] = 1.0, -1.0
            basis.append(M @ eta)
    return basis


def killing_eigs(basis):
    n = len(basis)
    B = np.zeros((n, n))
    ads = []
    for i in range(n):
        ad = np.zeros((n, n))
        for k in range(n):
            for target in range(n):
                C = basis[i] @ basis[target] - basis[target] @ basis[i]
                # 投影到 basis 上（basis 正交性由构造保证到常数因子）
                for l in range(n):
                    den = np.vdot(basis[l], basis[l]).real
                    if den > 1e-12:
                        c = np.vdot(basis[l], C) / den
                        if np.allclose(C, c * basis[l], atol=1e-8):
                            ad[l, target] = c.real
        ads.append(ad)
    for i in range(n):
        for j in range(n):
            B[i, j] = np.trace(ads[i] @ ads[j])
    return np.linalg.eigvalsh((B + B.T) / 2)


ev4 = killing_eigs(so_basis((1, 1, 1, 1)))
ev13 = killing_eigs(so_basis((1, -1, -1, -1)))
check("so(4) Killing 形式负定（紧）", bool(np.all(ev4 < 1e-9)), "eigs=%s" % np.round(ev4, 3).tolist())
check("so(1,3) Killing 形式不定（非紧）",
      bool(np.any(ev13 > 1e-9) and np.any(ev13 < -1e-9)), "eigs=%s" % np.round(ev13, 3).tolist())

# ======================================================================
head("F5  与既有 no-go 的合流")

check("R19 的 boost 缺失原文在位",
      "旋转双覆盖不能提供 boost" in R19 and "非紧、非阿贝尔的 boost" in R19)
check("R12／R13 的原文在位",
      "常数剖面" in R12 and "交换代数" in R12 and "谱半径" in R13 and "线性发散" in R13)
check("G62 的原生 Pauli 关系在位",
      "[\\sigma_x,\\sigma_y]=2i\\sigma_z" in G62 or "\\sigma_x,\\sigma_y" in G62)
check("文档写明 R34.1 是维数性质而非 Zero 缺失",
      "不是\"$Zero$ 缺 boost\"" in DOC or "任何有限维代数都缺 boost" in DOC)

# ======================================================================
head("F6  账本登记")

check("STATUS 已登记 R34",
      "### 2.34" in STATUS
      and "R34_finite_dimensional_boost_obstruction.md" in STATUS
      and "R34_check.py" in STATUS
      and "定理 R34.1" in STATUS)
check("INDEX 已收录 R34",
      "R34_finite_dimensional_boost_obstruction.md" in INDEX
      and "R34_check.py" in INDEX)

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
