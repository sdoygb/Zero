#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z13_check.py —— 【Zero 基础缺什么】的独立核验
==============================================
不访问网络，不修改任何项目文件。

独立实断言：
  F1  文档锚点：选择／读出原理、五个“哪一个”、E5、Z13-OPEN、两条路径
  F2  与 G62／G75／G76／G77 的原文自认逐字对上
  F3  Z0③ 原文在位（不设概率／没有选择规则）
  F4  E1–E4 账本在位，且 STATUS 已登记 Z13 评估
  F5  数值：分支均匀（Z0③ 无偏好）但闭环计数在非正则图上给出**非均匀**边权，
      且归一化后收敛到 Perron 形式（“计数权重 ≠ 概率测度”的可复算证据）
  F6  R13 的 L1 判决与 Z-CRIT-DER 被正确引用
"""
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


def head(t):
    print("\n" + "=" * 72 + "\n" + t + "\n" + "=" * 72)


def read(name):
    return io.open(os.path.join(HERE, name), encoding="utf-8").read()


DOC = read("Z13_zero_foundation_missing_principle.md")
STATUS = read("STATUS.md")
Z0 = read("Z0_zero_never_rests_single_axiom.md")
G62 = read("G62_quantum_sector_from_GNS_modular_flow_gleason.md")
G75 = read("G75_quantum_geometry_modular_readout.md")
G76 = read("G76_area_law_in_2d.md")
G77 = read("G77_staggered_coupling_from_A5.md")
R13 = read("R13_L1_strong_resolvent_attempt.md")


# ---------------------------------------------------------------- F1
head("F1  文档锚点")

check("标题指出缺的是选择／读出原理",
      "选择／读出原理" in DOC or "选择/读出原理" in DOC)
check("五个“哪一个”都在位",
      all(x in DOC for x in ["E1", "E2", "E3", "E4", "E5"]))
check("E5 被明确写为账本漏记项",
      "账本漏记" in DOC and "量子读出" in DOC)
check("两条路径（政策 vs 基础改写）都在位",
      "政策路径" in DOC and "基础改写路径" in DOC)
check("内生候选 Z13-OPEN 在位",
      "Z13-OPEN" in DOC and "计数典型性" in DOC)
check("三角困境与‘计数权重 ≠ 概率测度’被引用",
      "三角困境" in DOC and "计数权重" in DOC and "概率测度" in DOC)


# ---------------------------------------------------------------- F2
head("F2  与底层文档的原文自认逐字对上")

check("G77 §4 自认自由费米形式未从 Z3 推出",
      "未" in G77 and "从 Z3 严格推出自由费米形式" in G77)
check("G75／G76 自认交错质量是识别",
      "交错质量是" in G75 and "识别" in G75
      and "交错质量仍是" in G76 and "识别" in G76)
check("G62 自认加性域输入未从 Z0 条款导出",
      "未从 Z0 条款导出" in G62)
check("G62 自认粗粒化 π 同时承重宏观与量子",
      "同时承重宏观动力学与量子动力学" in G62)
check("R13 的第一缺口即 E5（Z-CRIT-DER）",
      "Z-CRIT-DER" in R13 and "自由费米" in R13)


# ---------------------------------------------------------------- F3
head("F3  Z0③ 原文")

check("Z0③ 写明不设概率",
      "不设概率" in Z0)
check("Z0③ 写明没有选择规则",
      "没有选择规则" in Z0)
check("Z0 把全分支当推论（分水岭）",
      "分水岭" in Z0 and "全分支" in Z0)


# ---------------------------------------------------------------- F4
head("F4  账本与 STATUS 登记")

check("E1–E4 账本在 STATUS／Z13 中都被引",
      all(("**%s**" % e) in STATUS for e in ("E1", "E2", "E3", "E4"))
      and all(x in DOC for x in ["E1", "E2", "E3", "E4"]))
check("STATUS 已登记 Z13 底层评估",
      "Z13" in STATUS and "选择／读出原理" in STATUS)
check("STATUS 把 E5 登记为账本漏记／建议具名",
      "E5" in STATUS and "量子读出" in STATUS)


# ---------------------------------------------------------------- F5
head("F5  数值：无偏好 ⇒ 分支均匀；计数 ⇒ 边权非均匀（Perron）")


def nonregular_graph():
    """K4 去一边 (1,2)，再在节点 1 挂一悬挂点 4。度序列 (3,3,3,2,1)。"""
    A = np.zeros((5, 5), dtype=float)
    for i in range(4):
        for j in range(i + 1, 4):
            if (i, j) == (1, 2):
                continue
            A[i, j] = A[j, i] = 1.0
    A[1, 4] = A[4, 1] = 1.0
    return A


def perron(A):
    vals, vecs = np.linalg.eigh(A)
    phi = np.abs(vecs[:, -1])
    return phi / np.linalg.norm(phi)


def closed_walk_edge_weights(A, k):
    """A_ij * (A^{k-1})_ij：长度 k 闭合游走对无向边的穿越强度（差常数因子 k）。"""
    P = np.linalg.matrix_power(A, k - 1)
    return A * P


A = nonregular_graph()
deg = np.diag(A).copy()
deg = A.sum(axis=1)

# (i) 分支层：Z0③ 无偏好 ⇒ 每条有向移动重数相同（都是 1）。等权 ⇒ 度差异只能来自图结构。
branch_weights = A[A > 0]
check("分支层无偏好：每条允许移动的重数相同",
      np.allclose(branch_weights, 1.0),
      "边数=%d，全部重数=1" % int(A.sum() / 2))
check("图确实非正则（否则计数退化为均匀）",
      float(deg.max() - deg.min()) > 0.5,
      "度序列=%s" % deg.astype(int).tolist())

# (ii) 计数层：闭环计数权重在非正则图上**非均匀**。
k = 40
W = closed_walk_edge_weights(A, k)
Wn = W / W.sum()
edge_mask = A > 0
spread = float(Wn[edge_mask].max() / Wn[edge_mask].min())
check("闭环计数给出非均匀边权（计数权重 ≠ 均匀测度）",
      spread > 1.05,
      "边权 max/min = %.4f" % spread)

# (iii) 归一化后收敛到 Perron 形式 φ_i φ_j A_ij。
phi = perron(A)
Per = np.outer(phi, phi) * A
Per = Per / Per.sum()
rel = float(np.abs(Wn - Per)[edge_mask].max() / Per[edge_mask].max())
check("归一化闭环计数收敛到 Perron 形式 φ_iφ_jA_ij",
      rel < 1e-6,
      "k=%d 时最大相对偏差 = %.2e" % (k, rel))

# (iv) 稳定性：k 加倍，偏差继续下降（指数收敛）。
W2 = closed_walk_edge_weights(A, 2 * k)
W2 = W2 / W2.sum()
rel2 = float(np.abs(W2 - Per)[edge_mask].max() / Per[edge_mask].max())
check("Perron 收敛随 k 指数稳定",
      rel2 < rel,
      "k=%d: %.2e → k=%d: %.2e" % (k, rel, 2 * k, rel2))


# ---------------------------------------------------------------- F6
head("F6  L1 判决的引用一致性")

check("Z13 指向 R13 的 L1 判决",
      "R13" in DOC and "排除了原样强预解收敛" in DOC)
check("STATUS 仍保持 L1 开放／条件恢复",
      "L1 几何模极限" in STATUS and "条件恢复" in STATUS)


# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
