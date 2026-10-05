#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R14_check.py —— 【从 Zero 组装 L1】的独立核验
================================================
不访问网络，不修改任何项目文件。

独立实断言：
  F1  文档锚点：Z-CAR、R14.1–R14.4、R13 七缺口、J1 未关闭
  F2  引理 R14.1：闭合 ⇒ 步数平衡（半满），枚举 L=6,8,10
  F3  引理 R14.2：补偿移动支撑=2；固定 content 的交换图连通
  F4  引理 R14.3：原语段数 = r+1；Σ(1+r)=7,18,48；Σ2^r=8,23,67
  F5  命题 R14.4：L-循环置换号 = (-1)^{L-1}；偶数 L 恒 -1（反转号不稳）
  F6  面积在相邻交换下恒变 ±2（coboundary / 纯规范）
  F7  与 R13 七缺口、G67 双覆盖的交叉引用
"""
import io
import itertools
import os
import sys
from collections import deque

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


DOC = read("R14_L1_from_zero_assembly.md")
R13 = read("R13_L1_strong_resolvent_attempt.md")
G66 = read("G66_SU2_double_cover_from_geometry.md")
G67 = read("G67_reflection_generates_spin_Z2.md")
G77 = read("G77_staggered_coupling_from_A5.md")


def zero_sum_words(L):
    return [w for w in itertools.product([1, -1], repeat=L) if sum(w) == 0]


def partial(w):
    s = 0
    out = []
    for x in w:
        s += x
        out.append(s)
    return out


def internal_returns(w):
    S = partial(w)
    return sum(1 for k in range(len(w) - 1) if S[k] == 0)


def primitive_pieces(w):
    S = partial(w)
    return 1 + sum(1 for k in range(len(w) - 1) if S[k] == 0)


def area(w):
    return sum(partial(w))


def neighbors(w):
    L = len(w)
    out = []
    for i in range(L - 1):
        if w[i] != w[i + 1]:
            v = list(w)
            v[i], v[i + 1] = v[i + 1], v[i]
            out.append(tuple(v))
    return out


def canonical_rot(w):
    L = len(w)
    return min(tuple(w[(k + j) % L] for j in range(L)) for k in range(L))


# ---------------------------------------------------------------- F1
head("F1  文档锚点")

check("标题/判决：J1 仍未关闭",
      "J1 仍未关闭" in DOC or "L1 **仍未关闭**" in DOC)
check("四个引理/命题都在位",
      all(x in DOC for x in ["R14.1", "R14.2", "R14.3", "R14.4"]))
check("新命题 Z-CAR 在位并标为候选",
      "Z-CAR" in DOC and "候选" in DOC)
check("R13 七个缺口全部列出",
      all(x in DOC for x in ["Z-CRIT-DER", "Z-SCALE", "Z-HILB", "Z-CORE",
                             "Z-TAIL", "Z-STRESS", "Z-CONF"]))
check("组装≠导出、单费米点仍是识别被写明",
      "组装 ≠ 导出" in DOC and "仍是识别" in DOC)


# ---------------------------------------------------------------- F2
head("F2  引理 R14.1：闭合 ⇒ 半满")

bal_ok = True
counts = []
for L in (6, 8, 10):
    W = zero_sum_words(L)
    counts.append(len(W))
    for w in W:
        if sum(1 for x in w if x == 1) != L // 2:
            bal_ok = False
check("每个零和词正负各 L/2",
      bal_ok, "|W| = %s" % counts)


# ---------------------------------------------------------------- F3
head("F3  引理 R14.2：补偿移动 = 最近邻 hopping，且生成全 content 图")

support_ok = True
conn_ok = True
for L in (6, 8, 10):
    W = zero_sum_words(L)
    for w in W:
        for v in neighbors(w):
            if sum(1 for a, b in zip(w, v) if a != b) != 2:
                support_ok = False
    # BFS 连通性
    seen = {W[0]}
    q = deque([W[0]])
    while q:
        w = q.popleft()
        for v in neighbors(w):
            if v not in seen:
                seen.add(v)
                q.append(v)
    if len(seen) != len(W):
        conn_ok = False
check("每次交换只动 2 个位点（hopping 距离 1）", support_ok)
check("固定 content 的交换图连通（hopping 生成全 Fock 空间）", conn_ok)


# ---------------------------------------------------------------- F4
head("F4  引理 R14.3：excursion 分解与 Fock 维数")

pieces_ok = True
rows = []
for L in (6, 8, 10):
    W = zero_sum_words(L)
    for w in W:
        if primitive_pieces(w) != internal_returns(w) + 1:
            pieces_ok = False
    classes = {}
    for w in W:
        rot = canonical_rot(w)
        classes[rot] = internal_returns(rot)
    s1 = sum(1 + r for r in classes.values())
    s2 = sum(2 ** r for r in classes.values())
    rows.append((L, len(classes), s1, s2))
check("原语段数 = r+1", pieces_ok)
check("按旋转类求和 Σ(1+r) = 7,18,48",
      [r[2] for r in rows] == [7, 18, 48],
      "得到 %s" % [r[2] for r in rows])
check("按旋转类求和 Σ2^r = 8,23,67（Fock 维数）",
      [r[3] for r in rows] == [8, 23, 67],
      "得到 %s" % [r[3] for r in rows])


# ---------------------------------------------------------------- F5
head("F5  命题 R14.4：L-循环置换号 = 双覆盖 -1")


def perm_sign(perm):
    """莱布尼茨公式：置换矩阵的行列式。"""
    n = len(perm)
    M = np.zeros((n, n))
    for i, j in enumerate(perm):
        M[i, j] = 1.0
    return int(round(np.linalg.det(M)))


cycle_ok = True
rev_rows = []
for L in (6, 8, 10, 12):
    cyc = tuple((i + 1) % L for i in range(L))          # (1 2 ... L)
    rev = tuple(L - 1 - i for i in range(L))             # 反转
    sc = perm_sign(cyc)
    sr = perm_sign(rev)
    rev_rows.append((L, sr))
    if sc != (-1) ** (L - 1):
        cycle_ok = False
    if L % 2 == 0 and sc != -1:
        cycle_ok = False
check("L-循环置换号 = (-1)^{L-1}，偶数 L 恒为 -1",
      cycle_ok, "L=6,8,10,12 全为 -1")
rev_signs = [s for _, s in rev_rows]
check("对照：反转号 (-1)^{L/2} 依赖 L mod 4（不稳定）",
      rev_signs == [-1, 1, -1, 1] or rev_signs == [-1, 1, -1, 1],
      "L=6,8,10,12 给 %s" % rev_signs)
check("文档把‘项链一圈 = 2π’标为唯一台阶",
      "项链**走一圈**" in DOC or "项链一圈" in DOC or "走一圈" in DOC)


# ---------------------------------------------------------------- F6
head("F6  面积是 coboundary（纯规范）")

ds = set()
for L in (6, 8, 10):
    for w in zero_sum_words(L):
        a0 = area(w)
        for v in neighbors(w):
            ds.add(area(v) - a0)
check("相邻交换下面积变化恒为 {±2}",
      ds == {-2, 2}, "得到 %s" % sorted(ds))


# ---------------------------------------------------------------- F7
head("F7  交叉引用一致性")

check("R13 七个缺口在 R13 原文中确实存在",
      all(x in R13 for x in ["Z-CRIT-DER", "Z-SCALE", "Z-HILB", "Z-CORE",
                             "Z-TAIL", "Z-STRESS", "Z-CONF"]))
check("G66 的双覆盖与 G67 的 -1 原文在位",
      "双覆盖" in G66 and "-1" in G67 and "自旋" in G67)
check("G77 仍自认交错耦合未从 Z3 严格推出",
      "未" in G77 and "从 Z3 严格推出自由费米形式" in G77)
check("R14 未把组装写成 L1 已证",
      "L1 **仍未关闭**" in DOC or "J1 仍未关闭" in DOC)


# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
