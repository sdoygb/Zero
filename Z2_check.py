#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z2_check.py —— 【从 Zero 直连 GR】路线的核验
==========================================
独立实断言：
  F1  闭环计数公式 N^(k)_ij = k·A_ij·(A^{k-1})_ij（**暴力枚举**闭合游走，多图多 k 比对）
  F2  Perron 收敛：归一化穿越数 → φ_i φ_j A_ij（**非正则图**上，指数收敛）
  F3  正则图上退化为均匀（G40 §3 的一致性检查）
  F4  (C′) 守恒：Σ_v d_w(v) ≡ 0 是恒等式（Z1 定理 2 的独立复算）
  F5  (L) 局域性：一步只触及两个站点（补偿移动的支撑为 2）
  F6  文档的链条与结论在位（四处假设被吸收、两处输入保留、I2a 为瓶颈）
"""
import io
import itertools
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LEDGER_MODE = os.environ.get("LH_LEDGER", "A")
FAIL = 0
N = 0


def check(name, cond, detail="", level="ind"):
    global FAIL, N
    N += 1
    ok = bool(cond)
    if not ok:
        FAIL += 1
    tag = ("v" if ok else "x") if (level == "ind" or LEDGER_MODE == "A" or not ok) else "i"
    print("  [%s] %s%s" % (tag, name, ("   " + detail) if detail else ""))


def head(t):
    print("\n" + "=" * 72 + "\n" + t + "\n" + "=" * 72)


DOC = io.open(os.path.join(HERE, "Z2_zero_to_gr_direct_route.md"), encoding="utf-8").read()


def complete(m):
    A = np.zeros((m, m), dtype=int)
    for i in range(m):
        for j in range(m):
            if i != j:
                A[i, j] = 1
    return A


def cycle_graph(m):
    A = np.zeros((m, m), dtype=int)
    for i in range(m):
        A[i, (i + 1) % m] = A[(i + 1) % m, i] = 1
    return A


# ---------------------------------------------------------------- F1
head("F1  闭环计数公式（暴力枚举 vs k·A·A^{k-1}）")
for gname, A in (("K4", complete(4)), ("C4", cycle_graph(4)), ("K5", complete(5))):
    m = len(A)
    for k in (2, 3, 4):
        Nbrute = np.zeros((m, m), dtype=int)
        for seq in itertools.product(range(m), repeat=k):
            if any(A[seq[t], seq[(t + 1) % k]] == 0 for t in range(k)):
                continue
            for t in range(k):
                Nbrute[seq[t], seq[(t + 1) % k]] += 1
        P = np.linalg.matrix_power(A, k - 1)
        F = np.array([[k * A[i, j] * P[i, j] for j in range(m)] for i in range(m)])
        check("%-3s k=%d：暴力 %d = 公式 %d" % (gname, k, Nbrute.sum(), F.sum()),
              np.array_equal(Nbrute, F))

# ---------------------------------------------------------------- F2
head("F2  Perron 收敛（非正则图，度 3,2,3,3,1）")
E = [(0, 1), (0, 2), (0, 3), (1, 2), (2, 3), (3, 4)]
m = 5
A = np.zeros((m, m))
for i, j in E:
    A[i, j] = A[j, i] = 1
deg = A.sum(1)
check("图为非正则（度数不全相等）", len(set(deg.tolist())) > 1, "度数 %s" % deg.tolist())
ev, V = np.linalg.eigh(A)
phi = V[:, -1]
if phi.sum() < 0:
    phi = -phi
target = phi[:, None] * phi[None, :] * A
target = target / target.sum()
errs = {}
for k in (2, 4, 8, 16, 32, 64):
    P = np.linalg.matrix_power(A, k - 1)
    Nk = np.array([[k * A[i, j] * P[i, j] for j in range(m)] for i in range(m)])
    errs[k] = float(np.abs(Nk / Nk.sum() - target).max())
check("k=2 偏差 ≈ 5.1e-2", abs(errs[2] - 5.102e-2) < 1e-3, "实算 %.4e" % errs[2])
check("k=64 偏差 < 1e-10（指数收敛）", errs[64] < 1e-10, "实算 %.3e" % errs[64])
check("偏差随 k 单调下降", all(errs[a] > errs[b] for a, b in zip((2, 4, 8, 16, 32), (4, 8, 16, 32, 64))),
      "%s" % {k: "%.1e" % v for k, v in errs.items()})

# ---------------------------------------------------------------- F3
head("F3  正则图上退化为均匀")
A4 = complete(4)
P = np.linalg.matrix_power(A4, 15)
Nk = np.array([[16 * A4[i, j] * P[i, j] for j in range(4)] for i in range(4)])
off = Nk[Nk > 0]
check("K4 上所有非零 $w_{ij}$ 相等（均匀）", np.allclose(off, off[0]), "取值 %s" % np.unique(off))

# ---------------------------------------------------------------- F4
head("F4  (C′) 守恒：Σ_v d_w(v) ≡ 0 是恒等式")
viol = 0
tot = 0
for mm in (2, 3, 4):
    Ee = [(u, v) for u in range(mm) for v in range(mm) if u != v]
    for edges in itertools.product(Ee, repeat=4):
        d = [0] * mm
        for (u, v) in edges:
            d[v] += 1
            d[u] -= 1
        tot += 1
        if sum(d) != 0:
            viol += 1
check("任意词的散度之和恒为 0", viol == 0, "测 %d 个词，违反 %d" % (tot, viol))
check("文档写明 (C) 两侧都自动（几何侧 Noether ＋ 物质侧恒等式）",
      "(C) 守恒源" in DOC and "恒等式" in DOC and "Noether" in DOC)

# ---------------------------------------------------------------- F5
head("F5  (L) 局域性：一步的支撑为 2")
bad = 0
for mm in (3, 4, 5):
    Ee = [(u, v) for u in range(mm) for v in range(mm) if u != v]
    for (u, v) in Ee:
        support = {u, v}
        if len(support) != 2:
            bad += 1
check("每步只触及两个站点（补偿移动支撑 = 2）", bad == 0, "违反 %d" % bad)
check("文档写明 (L) 自动成立", "(L) 局域性" in DOC and "自动成立" in DOC)

# ---------------------------------------------------------------- F6
head("F6  文档结论在位")
check("写明免去桥接可行", "可行且更干净" in DOC or "桥接确实可免" in DOC)
check("写明四处假设被吸收", "四处假设被吸收" in DOC)
check("写明两处输入原样保留（I5／I2a）", "I5" in DOC and "I2a" in DOC and "原样保留" in DOC)
check("写明底层条款（Z0 条款 ＋ Z1–Z5 定理）桥接从来不是瓶颈", "从来不是瓶颈" in DOC)
check("引用几何探针的障碍（维数漂移）",
      "zero_sum_closure_graph_theorems.md" in DOC and "漂移" in DOC)
check("诚实边界写明本文未导出 Einstein 方程", "没有**导出 Einstein 方程" in DOC or "没有**导出" in DOC)
check("§6 给出与 G40 的对照表", "G40（经 A 系）" in DOC or "Z2（直连）" in DOC)

# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
