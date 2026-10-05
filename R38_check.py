#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R38_check.py -- 核验"共同起因 ＋ 未记录自由度 ⇒ 纠缠"的探针与文档。
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


DOC = read("R38_entanglement_from_shared_closure_origin.md")
R36 = read("R36_chsh_bell_locality.md")
R37 = read("R37_kcbs_contextuality.md")
R28 = read("R28_phase_identity_cluster_reduction.md")
Z15 = read("Z15_zcar_no_go_and_jordan_wigner_readout.md")
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")
RES = json.load(io.open(os.path.join(HERE, "R38_jw_entanglement_results.json"), encoding="utf-8"))

I2 = np.eye(2, dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
sp = np.array([[0, 0], [1, 0]], dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
P = [sx, sy, sz]
N = 4


def op(single, site):
    out = np.array([[1.0 + 0j]])
    for k in range(N):
        out = np.kron(out, single if k == site else I2)
    return out


def c_dag(j):
    out = np.array([[1.0 + 0j]])
    for k in range(N):
        out = np.kron(out, sz if k < j else (sp if k == j else I2))
    return out


def s_max(rho):
    T = np.array([[np.trace(rho @ np.kron(P[i], P[j])).real for j in range(3)] for i in range(3)])
    u = np.linalg.svd(T, compute_uv=False)
    return 2.0 * float(np.sqrt(u[0] ** 2 + u[1] ** 2))


def reduce2(psi):
    T = psi.reshape([2] * N)
    rho = np.tensordot(T, T.conj(), axes=([list(range(2, N)), list(range(2, N))]))
    return rho.reshape(4, 4)


vac = np.zeros(2 ** N, dtype=complex); vac[0] = 1.0

# ======================================================================
head("F1  文档结构：H 的判据、两条路径、边界")

check("标题与性质为机制立项＋正面验证",
      "纠缠的涌现机制" in DOC and "共同起因" in DOC and "探针验证（正面）" in DOC)
check("主判据 (R38-1) 与两条路径 (R38-2) 在位",
      "(R38-1)" in DOC and "(R38-2)" in DOC
      and "EDGE-CONNECTION" in DOC and "未记录自由度" in DOC)
check("决定性问题 (R38-3) 在位",
      "(R38-3)" in DOC and "失明还是记录" in DOC)
check("边界写明 Z-READ 未导出、链是 1+1 维",
      "`Z-READ`" in DOC and "1+1" in DOC)

# ======================================================================
head("F2  独立复算：JW 单粒子扇区的跨位点纠缠")

check("Z15 的 JW 构造与 CAR 在位",
      "Jordan" in Z15 and "CAR" in Z15)
for j in range(N):
    psi = c_dag(j) @ vac
    psi = psi / np.linalg.norm(psi)
    S = s_max(reduce2(psi))
    check("单位点 c_%d†|0> 不纠缠（S=2）" % j, abs(S - 2) < 1e-9, "S=%.6f" % S)
sup = (c_dag(0) + c_dag(1)) @ vac / np.sqrt(2)
sup = sup / np.linalg.norm(sup)
S_sup = s_max(reduce2(sup))
check("跨位点相干叠加给 S=2√2（最大纠缠）",
      abs(S_sup - 2 * np.sqrt(2)) < 1e-9, "S=%.6f" % S_sup)
# 经典混合（位点被记录）
rho_mix = 0.5 * reduce2(c_dag(0) @ vac) + 0.5 * reduce2(c_dag(1) @ vac)
S_mix = s_max(rho_mix)
check("位点被记录（经典混合）退回 S=2", abs(S_mix - 2) < 1e-9, "S=%.6f" % S_mix)
check("真空不纠缠", abs(s_max(reduce2(vac)) - 2) < 1e-9)

# ======================================================================
head("F3  与结果 JSON 及既有账本一致")

check("JSON 判定叠加态纠缠、混合态不纠缠",
      RES["verdict"]["single_particle_superposition_across_two_sites_is_entangled"] is True
      and RES["verdict"]["which_site_record_decides"] is True
      and abs(RES["verdict"]["S_superposition"] - 2 * np.sqrt(2)) < 1e-9)
check("JSON 写明共同起因与未记录自由度",
      "费米子" in RES["verdict"]["shared_origin"]
      and "位点" in RES["verdict"]["unrecorded_dof"])
check("R36 的 CHSH 判据被正确引用（u1^2+u2^2>1）", "u_1^2+u_2^2" in R36 or "u1^2+u2^2" in R36)
check("R37 的单体结论未被改动", "语境" in R37)
check("R28 的 EDGE-CONNECTION 在位", "EDGE-CONNECTION" in R28)
check("STATUS 已登记 R38",
      "### 2.38" in STATUS
      and "R38_entanglement_from_shared_closure_origin.md" in STATUS
      and "R38_check.py" in STATUS)
check("INDEX 已收录 R38",
      "R38_entanglement_from_shared_closure_origin.md" in INDEX
      and "R38_check.py" in INDEX)

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
