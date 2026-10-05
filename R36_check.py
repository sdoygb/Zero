#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R36_check.py -- 核验 CHSH 探针：机制校验、G68 非双体、G82 直积必 S=2、违反门槛量化。
"""

from __future__ import annotations

import io
import json
import os
import sys

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


DOC = read("R36_chsh_bell_locality.md")
G68 = read("G68_interference_from_coarse_graining.md")
G82 = read("G82_B_from_two_independent_Z2.md")
G88 = read("G88_measurement_as_typicality.md")
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")
RES = json.load(io.open(os.path.join(HERE, "R36_chsh_results.json"), encoding="utf-8"))

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)
P = [sx, sy, sz]


def s_max(rho):
    T = np.array([[np.trace(rho @ np.kron(P[i], P[j])).real for j in range(3)] for i in range(3)])
    u = np.linalg.svd(T, compute_uv=False)
    return 2.0 * float(np.sqrt(u[0] ** 2 + u[1] ** 2)), u


def bell():
    psi = np.zeros(4, dtype=complex)
    psi[0] = psi[3] = 1 / np.sqrt(2)
    return np.outer(psi, psi.conj())


# ======================================================================
head("F1  文档结构：判据、非双体判定、门槛、缺口名称")

check("标题与性质为 CHSH 探针＋量化门槛",
      "CHSH 探针" in DOC and "Bell 局域" in DOC and "违反所需的门槛已量化" in DOC)
check("Horodecki 判据 (R36-4) 与直积结论 (R36-2)(R36-3) 在位",
      "(R36-2)" in DOC and "(R36-3)" in DOC and "(R36-4)" in DOC
      and "S_{\\max}=2\\sqrt{u_1^2+u_2^2}" in DOC)
check("G68 非双体 (R36-1) 与缺口名称 (R36-5) 在位",
      "(R36-1)" in DOC and "(R36-5)" in DOC and "单系统多路径" in DOC)
check("边界写明：没有否定量子性",
      "**没有否定量子性**" in DOC and "未能检验" in DOC)
check("与 D31 的三项登记对上",
      "复合系统与张量积" in DOC and "参照时钟的选择" in DOC)

# ======================================================================
head("F2  独立复算：机制校验")

s_bell, _ = s_max(bell())
check("Bell 态给 2√2", abs(s_bell - 2 * np.sqrt(2)) < 1e-9, "S=%.6f" % s_bell)
rho_prod = np.kron(np.diag([1, 0]).astype(complex), np.diag([1, 0]).astype(complex))
s_prod, u_prod = s_max(rho_prod)
check("直积态给 2（T 秩 1）", abs(s_prod - 2) < 1e-9 and abs(u_prod[1]) < 1e-12)
rho_cc = 0.5 * (np.kron(np.diag([1, 0]).astype(complex), np.diag([1, 0]).astype(complex))
                + np.kron(np.diag([0, 1]).astype(complex), np.diag([0, 1]).astype(complex)))
s_cc, _ = s_max(rho_cc)
check("经典关联态给 2（取等）", abs(s_cc - 2) < 1e-9)
ps = np.linspace(0.0, 1.0, 4001)
thr = next((float(p) for p in ps if s_max(p * bell() + (1 - p) * np.eye(4, dtype=complex) / 4)[0] > 2 + 1e-9), None)
check("Werner 门槛 ≈ 1/√2", thr is not None and abs(thr - 1 / np.sqrt(2)) < 5e-3,
      "p*=%.4f" % thr)

# ======================================================================
head("F3  独立复算：G68 设定不是双体；G82 直积必 S=2")

a68 = np.array([0.6 + 0.3j, 0.5 - 0.2j, 0.2 - 0.1j, -0.35 + 0.15j])
check("G68 振幅复算：Σ|a|²=0.935、合并=1.245、交叉=+0.310",
      abs(np.sum(np.abs(a68) ** 2) - 0.935) < 1e-9
      and abs(abs(a68[0] + a68[1]) ** 2 + abs(a68[2] + a68[3]) ** 2 - 1.245) < 1e-9
      and abs(float(sum(abs(a68[i] + a68[j]) ** 2 - abs(a68[i]) ** 2 - abs(a68[j]) ** 2
                       for i, j in [(0, 1), (2, 3)])) - 0.310) < 1e-9)
check("G68 文档确认干涉=合并+线性（非 Bell 结构）",
      "合并路径" in G68 and "振幅线性" in G68)
check("G82 文档给出 B=2×2（两独立 Z2）", "2\\times2=4" in G82 or "B=2\\times2=4" in G82)
s_z2, u_z2 = s_max(np.kron(np.diag([1, 0]).astype(complex), np.diag([1, 0]).astype(complex)))
check("G82 直积态 S=2 且 T 秩 1", abs(s_z2 - 2) < 1e-9 and abs(u_z2[1]) < 1e-12)
rho_zz = np.eye(4, dtype=complex) / 4 + np.kron(sz, sz) / 4
s_zz, u_zz = s_max(rho_zz)
check("σ_z⊗σ_z 关联拉满仍 S=2（单方向关联不违反）",
      abs(s_zz - 2) < 1e-9 and abs(u_zz[1]) < 1e-12)

# ======================================================================
head("F4  违反门槛：T 需秩≥2 且 u1²+u2²>1")

def two_dir(r):
    return np.eye(4, dtype=complex) / 4 + r * np.kron(sx, sx) / 4 + r * np.kron(sy, sy) / 4

s05, _ = s_max(two_dir(0.5))
s08, u08 = s_max(two_dir(0.8))
check("r=0.5 不违反（S=√2）", abs(s05 - np.sqrt(2)) < 1e-9)
check("r=0.8 违反（S≈2.2627）", s08 > 2 and abs(s08 - 2.2627416997969525) < 1e-9,
      "S=%.4f" % s08)
check("两方向态给出秩 2 的 T", abs(u08[2]) < 1e-12 and u08[0] > 0 and u08[1] > 0)

# ======================================================================
head("F5  与结果 JSON 及账本一致")

check("JSON 判定 G68 非双体、G82 不纠缠",
      RES["verdict"]["G68_setting_is_bipartite"] is False
      and RES["verdict"]["G82_two_Z2_state_entangled"] is False
      and RES["verdict"]["any_native_bipartition_violates"] is False)
check("JSON 给出量化目标", "u1^2+u2^2>1" in RES["verdict"]["quantitative_target"])
check("G88 引用 D31 的未恢复项在位",
      "复合系统与张量积" in G88)
check("STATUS 已登记 R36",
      "### 2.36" in STATUS
      and "R36_chsh_bell_locality.md" in STATUS
      and "R36_check.py" in STATUS
      and "Bell 局域" in STATUS)
check("INDEX 已收录 R36",
      "R36_chsh_bell_locality.md" in INDEX and "R36_check.py" in INDEX)

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
