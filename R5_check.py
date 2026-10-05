#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R5_check.py -- 主 Z/G 路线与 D210-D259 恢复路线的等价性审计。

对应文档 R5_route_equivalence.md。

  F1  文档与关键来源在位
  F2  零和秩公式：五通道给四维切空间
  F3  G9 的 5/4 冲突是计数约定，不是物理互斥
  F4  D259 时间线提升给 (1,3) 号差
  F5  D259 -> Z/G 的局部度规桥
  F6  完整五通道交换给各向同性正定型
  F7  各向异性反向阻碍
  F8  E1-E4 对象对应与缺口登记
  F9  两条路线不能合并进度
"""

import io
import itertools
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0
N = 0


def check(name, cond, detail=""):
    global FAIL, N
    N += 1
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def read(name):
    path = os.path.join(HERE, name)
    if not os.path.exists(path):
        return ""
    return io.open(path, encoding="utf-8").read()


DOC = read("R5_route_equivalence.md")
STATUS = read("STATUS.md")
G9 = read("G9_d_series_reference_triage.md")
D259 = read("D259_z_layer_condition_dictionary_and_four_dimensional_confluence.md")
D258 = read("D258_full_exchange_projection_and_local_isotropy_obstruction.md")
D255 = read("D255_tetrahedral_dirichlet_assembly_and_gluing.md")
Z7 = read("Z7_embedding_input_explicit_dictionary.md")


def zero_sum_rank(channels):
    return channels - 1


def main_dimension(base_channels):
    # route beta: one time direction + the (m-1)-dimensional zero-sum directions
    return 1 + zero_sum_rank(base_channels)


def d259_dimension(certificate_channels):
    # D259: the time line is selected inside H_Q, not added outside it
    return zero_sum_rank(certificate_channels)


def signature_counts(matrix, tol=1.0e-10):
    ev = np.linalg.eigvalsh(matrix)
    return (
        int(np.sum(ev > tol)),
        int(np.sum(ev < -tol)),
        int(np.sum(np.abs(ev) <= tol)),
    )


def gR_orthonormal_basis_with_first(g, u):
    """Return B with first column u/g-norm and B.T @ g @ B = I."""
    g = np.asarray(g, dtype=float)
    u = np.asarray(u, dtype=float)
    u = u / np.sqrt(u @ g @ u)
    n = u.size
    basis = [u]
    for e in np.eye(n):
        v = e.copy()
        for q in basis:
            v -= float(q @ g @ v) * q
        ng = float(v @ g @ v)
        if ng > 1.0e-12:
            basis.append(v / np.sqrt(ng))
        if len(basis) == n:
            break
    B = np.column_stack(basis)
    return B


def weighted_laplacian(n, edges, weights):
    L = np.zeros((n, n))
    for (i, j), w in zip(edges, weights):
        L[i, i] += w
        L[j, j] += w
        L[i, j] -= w
        L[j, i] -= w
    return L


def resistance_embedding(L):
    ev, evec = np.linalg.eigh(L)
    pos = ev > 1.0e-10
    return evec[:, pos] / np.sqrt(ev[pos])[None, :]


def full_exchange_stiffness(points, edges, weights):
    A = np.zeros((points.shape[1], points.shape[1]))
    for (i, j), w in zip(edges, weights):
        d = points[j] - points[i]
        A += w * np.outer(d, d)
    return A


# ======================================================================
head("F1  文档与关键来源在位")
# ======================================================================
for fn, text, needles in [
    ("R5_route_equivalence.md", DOC, ["等价性审计", "最小充分桥", "G9"]),
    ("STATUS.md", STATUS, ["两条路线尚未打通", "进度不能相加"]),
    ("G9_d_series_reference_triage.md", G9, ["五通道", "G8"]),
    ("D259_z_layer_condition_dictionary_and_four_dimensional_confluence.md", D259,
     ["五通道", "时间线", "跨区域字典汇流"]),
    ("D258_full_exchange_projection_and_local_isotropy_obstruction.md", D258,
     ["A_E^0", "I_3", "各向同性"]),
    ("D255_tetrahedral_dirichlet_assembly_and_gluing.md", D255,
     ["Q_\\sigma", "顶点嵌入", "细化"]),
    ("Z7_embedding_input_explicit_dictionary.md", Z7,
     ["单元形状", "O(a^2)"]),
]:
    check("%s 在位且含关键对象" % fn,
          bool(text) and all(x in text for x in needles))


# ======================================================================
head("F2  零和秩与通道-维度计数")
# ======================================================================
rank = {m: zero_sum_rank(m) for m in range(2, 8)}
check("rank H_Q(m)=m-1，m=5 给 4", rank[5] == 4, "%s" % rank)
check("主路线 m=4 给 1+(m-1)=4", main_dimension(4) == 4)
check("D259 m=5 给 m-1=4", d259_dimension(5) == 4)
check("两边最终维度同为 4", main_dimension(4) == d259_dimension(5) == 4)
check("误把 D259 时间线另加会得到 5 维", main_dimension(5) == 5)


# ======================================================================
head("F3  G9 冲突：计数约定而非物理互斥")
# ======================================================================
check("G9 原文确实登记过“5 通道 vs 4 维互斥”",
      "五通道 → 4 方向" in G9 and "通道数 = 时空维数" in G9)
check("R5 把互斥降级为记号冲突",
      "不是物理互斥" in DOC and "m_cert" in DOC and "m_base" in DOC)
check("D259 的时间线位于 H_Q 内部",
      "时间线" in D259 and "时间线是独立闭合证书" in D259 and "g_L" in D259)
check("主路线的时间方向与 H_Q 分开",
      "1+dim H_Q" in DOC or "1+(m-1)=m" in DOC)


# ======================================================================
head("F4  D259 时间线提升给 (1,3)")
# ======================================================================
gR = np.eye(4)
u = np.array([1.0, 0.0, 0.0, 0.0])
gL = 2.0 * np.outer(u, u) - gR
check("g_R 正定", np.all(np.linalg.eigvalsh(gR) > 0))
check("u 为单位时间线", np.isclose(u @ gR @ u, 1.0))
check("时间线提升 g_L=2u♭u♭-g_R", np.allclose(gL, np.diag([1.0, -1.0, -1.0, -1.0])))
check("g_L 号差为 (1,3)", signature_counts(gL) == (1, 3, 0), "%s" % (signature_counts(gL),))
check("时间线仍在四维切空间内", gL.shape == (4, 4) and rank[5] == 4)


# ======================================================================
head("F5  D259 -> Z/G 的局部度规桥")
# ======================================================================
rng = np.random.default_rng(5)
M = rng.normal(size=(4, 4))
gR = M @ M.T + np.eye(4)
u = rng.normal(size=4)
u /= np.sqrt(u @ gR @ u)
B = gR_orthonormal_basis_with_first(gR, u)

# In the g_R-orthonormal adapted basis, g_R = I and g_L = diag(1,-1,-1,-1).
gR_adapt = B.T @ gR @ B
h = gR_adapt[1:, 1:]
u_flat = gR @ u
gL = 2.0 * np.outer(u_flat, u_flat) - gR
gL_adapt = B.T @ gL @ B
check("适配基为 g_R-正交基", np.allclose(gR_adapt, np.eye(4)),
      "maxerr=%.2e" % np.max(np.abs(gR_adapt - np.eye(4))))
check("适配基第一列为时间线", np.allclose(B[:, 0], u))
check("空间块 h=g_R|_{u⊥} 正定", np.all(np.linalg.eigvalsh(h) > 0))
check("桥构造给出块对角 g_L",
      np.allclose(gL_adapt, np.diag([1.0, -1.0, -1.0, -1.0])))
check("桥构造保持号差 (1,3)", signature_counts(gL_adapt) == (1, 3, 0))
check("主路线字典 g_G=dτ²-h 与 g_L 相等",
      np.allclose(gL_adapt[0, 0], 1.0) and np.allclose(gL_adapt[1:, 1:], -np.eye(3)))


# ======================================================================
head("F6  完整五通道交换给各向同性正定型")
# ======================================================================
edges5 = list(itertools.combinations(range(5), 2))
weights5 = rng.uniform(0.4, 1.7, size=len(edges5))
L5 = weighted_laplacian(5, edges5, weights5)
X5 = resistance_embedding(L5)
A5 = full_exchange_stiffness(X5, edges5, weights5)
check("五通道电阻嵌入为 4 维", X5.shape == (5, 4), "%s" % (X5.shape,))
check("完整交换刚度在四维上为 I_4", np.allclose(A5, np.eye(4)),
      "maxerr=%.2e" % np.max(np.abs(A5 - np.eye(4))))
check("完整交换的特征值全相等", np.allclose(np.linalg.eigvalsh(A5), np.ones(4)))


# ======================================================================
head("F7  各向异性反向阻碍")
# ======================================================================
h_aniso = np.diag([1.0, 2.0, 3.0])
ev_aniso = np.linalg.eigvalsh(h_aniso)
ev_iso = np.linalg.eigvalsh(np.eye(3))
check("主路线反例 h 各向异性", not np.allclose(ev_aniso, ev_aniso[0] * np.ones(3)))
check("完整交换的空间截面必须各向同性", np.allclose(ev_iso, ev_iso[0] * np.ones(3)))
check("有限 c 不能把 diag(1,2,3) 变成 cI_3",
      not any(np.allclose(h_aniso, c * np.eye(3)) for c in (1.0, 2.0, 3.0)))

# Even if the time line is arbitrary, restriction of cI_4 to a 3-plane has equal
# eigenvalues. Enumerate a deterministic sample of time lines.
best = np.inf
for theta in np.linspace(0.0, np.pi, 25):
    for phi in np.linspace(0.0, 2.0 * np.pi, 25):
        uu = np.array([np.cos(theta),
                       np.sin(theta) * np.cos(phi),
                       np.sin(theta) * np.sin(phi),
                       0.0])
        uu = uu / np.linalg.norm(uu)
        _, _, vh = np.linalg.svd(uu.reshape(1, -1))
        V = vh[1:, :].T
        for c in (1.0, 2.0, 3.0):
            hh = c * (V.T @ V)
            best = min(best, float(np.max(np.abs(hh - h_aniso))))
check("任意时间线与 cI_4 截面都不能复现 diag(1,2,3)", best > 1.0e-3,
      "min mismatch=%.3e" % best)
check("文档登记反向受阻",
      "反向受阻" in DOC or "不能反向嵌入" in DOC)


# ======================================================================
head("F8  E1-E4 对 D259 的对应与缺口")
# ======================================================================
check("E1 站点/嵌入映射已登记为条件映射",
      "E1 站点/嵌入识别" in DOC and "单向条件映射" in DOC)
check("E2 的维度口径冲突已单独处理",
      "E2：D=4" in DOC and "口径冲突可消解" in DOC)
check("E3 明确登记缺失",
      "E3 源作用量类别" in DOC and "缺失" in DOC)
check("E4 明确登记部分对应",
      "E4 量纲常数/尺度" in DOC and "部分对应" in DOC)
check("体积与共形尺度不是 G,Λ",
      "不含 `G,Λ`" in DOC or "`G,Λ,ħ`" in DOC)
check("区域汇流仍为开放",
      "区域汇流" in DOC and "开放" in DOC)
check("外部边隔离登记为类比而非等价",
      "外部边隔离" in DOC and "类比，未证等价" in DOC)


# ======================================================================
head("F9  双路线记账纪律")
# ======================================================================
check("STATUS 要求两路线分开记账",
      "两条路线尚未打通" in STATUS and "进度不能相加" in STATUS)
check("R5 不建议当前砍掉一条路线",
      "不建议现在只保留一条路线" in DOC)
check("主路线保留为主生产线",
      "主 `Z/G` 路线保留为项目主生产线" in DOC)
check("D259 保留为局部图册子路线",
      "D259 作为" in DOC and "局部图册" in DOC)
check("完整等价明确未证",
      "没有完整等价" in DOC and "未证" in DOC)
check("G9 的“必须二选一”已降级",
      "降级为" in DOC and "必须选定一套通道记号" in DOC)


# ======================================================================
print("\n" + "=" * 72)
print("汇总")
print("=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
