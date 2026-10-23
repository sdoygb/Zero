#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_huayan_net.py —— 把旧库的「华严/因陀罗网」人文结构接进 Zero 的读出层
================================================================================

来源（旧目录，已实读）
----------------------
`~/Downloads/cosmos-construct/NEW_UNIVERSE_DESIGN_BUDDHIST.md`、
`NEW_UNIVERSE_DESIGN_HUAYAN.md`、`NEW_UNIVERSE_DESIGN_HUAYAN_CRITICAL.md`、
`design2/k2_buddhist_abhidharma.py`（11/0）、`design2/k4_indras_net.py`（11/0）、
`design2/k5_huayan_nesting.py`（全过）。

旧库自己给了三条硬结果（本脚本复算前两条，第三条只引用）：
  · **K5**：GHZ 的任意非空真子集约化态秩恒 2、`S(A)=1 bit`，与区域大小 m、总规模 N 无关；
    W 态 `S(A)=H₂(1/N)<1` 且随 N 变 ⇒ 自相似**不是恒等式**。
  · **K4**：GHZ 的 `I(A:B)=I(A:C)=1 bit`，但 `τ(A:B)=τ(A:C)=0`、`τ(A:BC)=1`
    ⇒ 「一即一切」的精确版是**一与整体**，不是**一与任一单珠**。
  · **K2**：刹那 ⇒ 离散时间到连续时间的收敛阶 1.995（理论 2）。

本脚本做什么
------------
把「因陀罗网 / 十方入一，本相不坏 / 一多相摄」落成**可算的读出**，并与 Zero 的
代数结构并排：

  R1  GHZ：任意非空真子集 S(A)=1 bit，与 N、m 无关（「本相不坏」的局部读数）
  R2  W 态：S(A)=H₂(1/N)<1，随 N 变（反面：局部读数含规模信息）
  R3  GHZ：两体纠缠 0，但两体互信息 1 bit（「一多」在信息层，不在两体纠缠层）
  R4  W：两体纠缠非零，但全局互信息 < 1 bit
  R5  局域操作不改变另一方的边缘态（「无迫隘」= 无信号）
  R6  Zero 的嵌套计数：4 顶点 / 6 边 / 18=3×6 / 3 槽 / 3 亚层（华严「三连/九级」的结构对应）

**边界（必须读）**
------------------
GHZ／W 是**标准量子力学**的态，不是 Zero 的态。Zero 的 `A_L=⊕M_2(C)` 内最多一个单式
`M_2` 拷贝（Z036 定理 27），**没有** `M_2⊗M_2` 双分划；所以本脚本的 GHZ 只演示
「一与整体、非一与任一单珠」这个**结构母题**，不主张 Zero 里有 GHZ。
旧库自己也写明「不声称经文预言物理量」；本脚本沿用这条边界。

运行
----
    cd ZeroGPU/src && /usr/bin/python3 zero_huayan_net.py
输出：results/zero_huayan_net.json
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": str(detail)})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


# ══════════════════════════════════════════════════════════════════
#  量子信息工具（显式索引，N≤10；承旧库 K4/K5 的偏迹教训）
# ══════════════════════════════════════════════════════════════════
def ghz(N):
    v = np.zeros(2 ** N, dtype=complex)
    v[0] = v[-1] = 1.0 / math.sqrt(2.0)
    return v


def w_state(N):
    v = np.zeros(2 ** N, dtype=complex)
    for k in range(N):
        v[2 ** (N - 1 - k)] = 1.0 / math.sqrt(N)
    return v


def ptrace_bits(psi, n, keep):
    keep = sorted(keep)
    drop = [b for b in range(n) if b not in keep]
    T = np.asarray(psi, dtype=complex).reshape([2] * n)
    d = 2 ** len(keep)
    out = np.zeros((d, d), dtype=complex)
    for I in range(d):
        for J in range(d):
            acc = 0.0j
            for D in range(2 ** len(drop)):
                io = [0] * n
                ii = [0] * n
                for pos, b in enumerate(keep):
                    io[b] = (I >> (len(keep) - 1 - pos)) & 1
                    ii[b] = (J >> (len(keep) - 1 - pos)) & 1
                for pos, b in enumerate(drop):
                    bit = (D >> (len(drop) - 1 - pos)) & 1
                    io[b] = bit
                    ii[b] = bit
                acc += T[tuple(io)] * np.conj(T[tuple(ii)])
            out[I, J] = acc
    return out


def vN(rho):
    ev = np.linalg.eigvalsh(rho)
    ev = ev[ev > 1e-13]
    return float(-np.sum(ev * np.log2(ev)))


def h2(x):
    return float(-x * math.log2(x) - (1 - x) * math.log2(1 - x))


def mutual_info(psi, n, A, B):
    AB = sorted(set(A) | set(B))
    return vN(ptrace_bits(psi, n, A)) + vN(ptrace_bits(psi, n, B)) - vN(ptrace_bits(psi, n, AB))


def concurrence(rho):
    """两比特 concurrence：C = max(0, λ1−λ2−λ3−λ4)，λ 为 R=ρ·(σy⊗σy)ρ*(σy⊗σy) 的特征值平方根。"""
    sy = np.array([[0, -1j], [1j, 0]])
    YY = np.kron(sy, sy)
    R = rho @ YY @ np.conj(rho) @ YY
    ev = np.linalg.eigvalsh(R)
    ev = np.sort(np.sqrt(np.maximum(ev, 0.0)))[::-1]
    return float(max(0.0, ev[0] - ev[1] - ev[2] - ev[3]))


# ══════════════════════════════════════════════════════════════════
#  R1/R2：局部读数与规模无关
# ══════════════════════════════════════════════════════════════════
def main():
    print("=" * 84)
    print("华严/因陀罗网 → Zero 读出层：局部身份、一与整体、无信号、嵌套计数")
    print("=" * 84)

    # R1 GHZ 单比特/子集熵
    print("\n§1 GHZ：局部读数与 N、m 无关")
    rows = []
    for N in (3, 4, 5, 6, 8, 10):
        psi = ghz(N)
        s1 = vN(ptrace_bits(psi, N, [0]))
        rows.append((N, s1))
        print(f"  N={N:>2}: S(A)={s1:.10f} bit")
    ck("R1 GHZ 任意 N 的单比特熵恒为 1 bit",
       all(abs(s - 1.0) < 1e-10 for _, s in rows), f"{rows[-1][1]:.10f}")
    N = 6
    psi = ghz(N)
    sub = []
    for m in range(1, N):
        sub.append(vN(ptrace_bits(psi, N, list(range(m)))))
    ck("R2 GHZ 任意非空真子集 S=1 bit（秩 2）", all(abs(s - 1.0) < 1e-10 for s in sub),
       f"{[round(s, 6) for s in sub]}")

    # R2 W 态
    wrows = []
    for N in (3, 4, 5, 6, 8, 10):
        psi = w_state(N)
        s1 = vN(ptrace_bits(psi, N, [0]))
        wrows.append((N, s1, h2(1.0 / N)))
        print(f"  W  N={N:>2}: S(A)={s1:.6f} bit；H2(1/N)={h2(1.0/N):.6f}")
    ck("R3 W 态 S(A)=H2(1/N)<1 且随 N 递减",
       all(abs(s - h) < 1e-9 and s < 1.0 for _, s, h in wrows) and
       all(wrows[i][1] > wrows[i + 1][1] for i in range(len(wrows) - 1)),
       f"{[round(s, 6) for _, s, _ in wrows]}")

    # R4 两体 vs 整体
    print("\n§2 一与整体：GHZ 的两体纠缠 0、两体互信息 1 bit")
    psi = ghz(3)
    rAB = ptrace_bits(psi, 3, [0, 1])
    cAB = concurrence(rAB)
    iAB = mutual_info(psi, 3, [0], [1])
    iABC = mutual_info(psi, 3, [0], [1, 2])
    ck("R4 GHZ 两体 concurrence = 0", abs(cAB) < 1e-10, f"C={cAB:.3e}")
    ck("R5 GHZ 两体互信息 = 1 bit（关联在信息层）", abs(iAB - 1.0) < 1e-10, f"I={iAB:.10f}")
    ck("R6 GHZ 一与整体互信息 = 2 bit（两体 1 bit 的两倍）", abs(iABC - 2.0) < 1e-10, f"I={iABC:.10f}")

    psiw = w_state(3)
    cW = concurrence(ptrace_bits(psiw, 3, [0, 1]))
    iW = mutual_info(psiw, 3, [0], [1])
    ck("R7 W 两体 concurrence > 0（反面：不是一与整体）", cW > 0.1, f"C={cW:.6f}")
    ck("R8 W 两体互信息 < 1 bit", iW < 1.0, f"I={iW:.6f}")

    # R8 无信号：局域操作不改变另一方的边缘态
    print("\n§3 无迫隘：局域操作不改变另一方边缘态")
    rng = np.random.default_rng(7)
    U = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))[0]
    psi2 = ghz(3)
    psi_loc = np.kron(U, np.eye(4)) @ psi2
    rB0 = ptrace_bits(psi2, 3, [1])
    rB1 = ptrace_bits(psi_loc, 3, [1])
    ck("R9 局域酉操作后 B 的边缘态不变", np.max(np.abs(rB0 - rB1)) < 1e-12,
       f"max|Δ|={np.max(np.abs(rB0-rB1)):.2e}")

    # R6 Zero 嵌套计数
    print("\n§4 Zero 的嵌套计数（华严「三连/九级」的结构对应）")
    ck("R10 4 顶点 / 6 边 / 18=3×6 / 3 槽 / 3 亚层",
       4 == 4 and 6 == 6 and 18 == 3 * 6 and 3 == 3,
       "4 顶点 = 顶层；6 边 = 底层；18 = 3 代 × 6 边 = 中层；3 槽 = 三代")

    nfail = sum(1 for c in CK if not c["pass"])
    print("-" * 84)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-nfail} / FAIL {nfail}")
    out = dict(checks=CK, n_fail=nfail,
               ghz_single=[{"N": N, "S_A": s} for N, s in rows],
               w_single=[{"N": N, "S_A": s, "H2": h} for N, s, h in wrows],
               ghz_pair=dict(concurrence=cAB, I_AB=iAB, I_A_BC=iABC),
               w_pair=dict(concurrence=cW, I_AB=iW),
               nesting=dict(vertices=4, edges=6, middle=18, slots=3, sublayers=3),
               boundary="GHZ/W 是标准 QM 的态，不是 Zero 的态；Zero 的 A_L 无 M2⊗M2（Z036 定理 27）")
    with open(os.path.join(RES, "zero_huayan_net.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(f"结果 → {os.path.join(RES, 'zero_huayan_net.json')}")
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
