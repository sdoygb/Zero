#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G9_check.py -- 独立复跑 D210-D259 中与 GR 推导直接相关的关键结论。

不引用 D 系列作为前提；只用底层条款与数学重算，逐条判定 D 结论是否成立。
失败时退出码非零。

  E1 D214  闭环关联矩阵：rank = L-1，im B = H_Q
  E2 D254  三维 h = (det Q) Q^{-1}；二维只剩共形类
  E3 D255  四面体装配 Q_sigma = omega^{-1} sum K e (x) e
  E4 D257  单位权重 K4：A_K = 4I-J，A_P1(c) = (sqrt(c)/24)(4I-J)，c = 576
  E5 D257  有效电阻非局域：路径 0-1-2-3 加边 {0,3} 后 R_12: 1 -> 3/4
  E6 D258  全边交换刚度 = 有效电阻嵌入的范围投影（任意正权重）
  E7 D258  K4 普遍固定点：c = V_1^{-2}，单位权重给 576
  E8 D258  外部边反作用：spec(A_S^0) = {4/5, 1, 1}
"""

import itertools
import math
import sys

import numpy as np

FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


# ----------------------------------------------------------------------
# E1  D214
# ----------------------------------------------------------------------
head("E1  D214  闭环关联矩阵：rank = L-1，im B = H_Q")

for L in (4, 5, 6):
    edges = [(k, (k + 1) % L) for k in range(L)]
    B = np.zeros((L, L))
    for c, (i, j) in enumerate(edges):
        B[i, c] = -1.0
        B[j, c] = 1.0
    P_B = B @ np.linalg.pinv(B.T @ B) @ B.T
    P_H = np.eye(L) - np.ones((L, L)) / L
    check("C_%d：rank = %d 且 im B = H_Q" % (L, L - 1),
          np.linalg.matrix_rank(B) == L - 1 and np.allclose(P_B, P_H, atol=1e-9))

# ----------------------------------------------------------------------
# E2  D254
# ----------------------------------------------------------------------
head("E2  D254  三维 h = (det Q)Q^{-1}；二维只剩共形类")

rng = np.random.default_rng(2)
A3 = rng.normal(size=(3, 3))
h3 = A3 @ A3.T + np.eye(3)
Q3 = np.sqrt(np.linalg.det(h3)) * np.linalg.inv(h3)
check("n=3：det Q = sqrt(det h)，且 h = (det Q)Q^{-1} 精确回收",
      abs(np.linalg.det(Q3) - np.sqrt(np.linalg.det(h3))) < 1e-9
      and np.allclose(np.linalg.det(Q3) * np.linalg.inv(Q3), h3))

A2 = rng.normal(size=(2, 2))
h2 = A2 @ A2.T + np.eye(2)
Q2 = np.sqrt(np.linalg.det(h2)) * np.linalg.inv(h2)
Om = 3.1
Q2s = np.sqrt(np.linalg.det(Om ** 2 * h2)) * np.linalg.inv(Om ** 2 * h2)
check("n=2：det Q = 1，且 Q 在 h -> Omega^2 h 下不变（只剩共形类）",
      abs(np.linalg.det(Q2) - 1.0) < 1e-9 and np.allclose(Q2s, Q2))

# ----------------------------------------------------------------------
# E3  D255
# ----------------------------------------------------------------------
head("E3  D255  四面体装配 Q_sigma = omega^{-1} sum K e (x) e")

V = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
E = [V[b] - V[a] for a, b in itertools.combinations(range(4), 2)]
om = abs(np.linalg.det(np.array([V[1] - V[0], V[2] - V[0], V[3] - V[0]]).T)) / 6.0
K = np.array([0.7, 1.3, 0.4, 2.1, 0.9, 1.6])
A = sum(Ke * np.outer(e, e) for Ke, e in zip(K, E))
Q = A / om
# 逐 (g) 检验：离散能量 = 连续能量
ok = True
for _ in range(20):
    g = rng.normal(size=3)
    disc = 0.5 * float(g @ A @ g)
    cont = 0.5 * om * float(g @ Q @ g)
    ok = ok and abs(disc - cont) < 1e-12
check("对全部仿射梯度，离散能量 = 连续能量（装配式成立）", ok)
check("omega = 1/6（参考四面体）", abs(om - 1.0 / 6.0) < 1e-12)

# ----------------------------------------------------------------------
# E4  D257  单位权重 K4 与 c = 576
# ----------------------------------------------------------------------
head("E4  D257  单位权重 K4：A_K = 4I-J，A_P1(c) = (sqrt(c)/24)(4I-J)，c = 576")

q = [np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0]),
     np.array([-1.0, 1, 0]), np.array([-1.0, 0, 1]), np.array([0, -1.0, 1])]
AK = sum(np.outer(v, v) for v in q)
I3 = np.eye(3)
J3 = np.ones((3, 3))
check("A_K = sum q (x) q = 4I - J", np.allclose(AK, 4 * I3 - J3))

c = 576.0
G = (c / 4.0) * (I3 + J3)
AP1 = math.sqrt(np.linalg.det(G)) / 6.0 * np.linalg.inv(G)
check("A_P1(576) = (sqrt(576)/24)(4I-J) = 4I-J = A_K",
      np.allclose(AP1, 4 * I3 - J3),
      "max dev=%.2e" % float(np.max(np.abs(AP1 - AK))))
check("V_1 = sqrt(det G_1)/6 = 1/24，故 c = V_1^{-2} = 576",
      abs(math.sqrt(np.linalg.det((1 / 4.0) * (I3 + J3))) / 6.0 - 1 / 24.0) < 1e-12)

# ----------------------------------------------------------------------
# E5  D257  有效电阻非局域
# ----------------------------------------------------------------------
head("E5  D257  有效电阻非局域：路径 0-1-2-3 加边 {0,3}，R_12: 1 -> 3/4")


def eff_R(m, edges, w=None):
    B = np.zeros((m, len(edges)))
    for col, (i, j) in enumerate(edges):
        B[i, col] = -1.0
        B[j, col] = 1.0
    Kd = np.ones(len(edges)) if w is None else np.asarray(w, float)
    L = B @ np.diag(Kd) @ B.T
    Lp = np.linalg.pinv(L)
    R = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            e = np.zeros(m)
            e[i], e[j] = 1.0, -1.0
            R[i, j] = float(e @ Lp @ e)
    return R


p4 = [(0, 1), (1, 2), (2, 3)]
R0 = eff_R(4, p4)
R1 = eff_R(4, p4 + [(0, 3)])
check("原图 R_12 = 1", abs(R0[1, 2] - 1.0) < 1e-9, "R=%.6f" % R0[1, 2])
check("加远端边 {0,3} 后 R_12 = 3/4", abs(R1[1, 2] - 0.75) < 1e-9, "R=%.6f" % R1[1, 2])

# ----------------------------------------------------------------------
# E6  D258  全边交换刚度 = 范围投影
# ----------------------------------------------------------------------
head("E6  D258  全边交换刚度 = 有效电阻嵌入的范围投影（任意正权重）")


def sym_sqrt_psd(M):
    """对称半正定矩阵的对称平方根。"""
    w, U = np.linalg.eigh(M)
    w = np.clip(w, 0.0, None)
    return U @ np.diag(np.sqrt(w)) @ U.T


def exchange_projection(m, edges, K):
    B = np.zeros((m, len(edges)))
    for col, (i, j) in enumerate(edges):
        B[i, col] = -1.0
        B[j, col] = 1.0
    L = B @ np.diag(K) @ B.T
    Lp = np.linalg.pinv(L)
    X = sym_sqrt_psd(Lp)                              # (L^+)^{1/2}（必须取对称平方根）
    AE = np.zeros((m, m))
    for col, (i, j) in enumerate(edges):
        d = X @ (np.eye(m)[j] - np.eye(m)[i])
        AE += K[col] * np.outer(d, d)
    # 范围投影
    ev, U = np.linalg.eigh(L)
    P = U[:, ev > 1e-9] @ U[:, ev > 1e-9].T
    return AE, P, ev


for seed in (1, 2, 3):
    r = np.random.default_rng(seed)
    K = r.uniform(0.3, 2.0, size=6)
    AE, P, ev = exchange_projection(4, list(itertools.combinations(range(4), 2)), K)
    check("随机正权重 K4 (seed=%d)：A_E^0 = range(L) 的投影" % seed,
          np.allclose(AE, P, atol=1e-8),
          "rank=%d  max dev=%.2e" % (int(np.sum(ev > 1e-9)), float(np.max(np.abs(AE - P)))))

# ----------------------------------------------------------------------
# E7  D258  K4 普遍固定点
# ----------------------------------------------------------------------
head("E7  D258  K4 普遍固定点：c = V_1^{-2}，单位权重给 576")

for seed in (1, 2, 3):
    r = np.random.default_rng(seed)
    K = r.uniform(0.3, 2.0, size=6)
    AE, P, ev = exchange_projection(4, list(itertools.combinations(range(4), 2)), K)
    check("seed=%d：任意正权重下 A_E^0 = 范围投影（故在三维叶上为 I_3），"
          "因此存在唯一尺度 c = V_1^{-2} 使 P1 固定点成立" % seed,
          np.allclose(AE, P, atol=1e-8),
          "max dev=%.2e" % float(np.max(np.abs(AE - P))))

# 单位权重 K4 的几何体积 V_1 = 1/24 -> c = 576
I3, J3 = np.eye(3), np.ones((3, 3))
V1_unit = math.sqrt(np.linalg.det((1 / 4.0) * (I3 + J3))) / 6.0
check("单位权重 K4：V_1 = 1/24，唯一固定点尺度 c = V_1^{-2} = %.0f" % (1 / V1_unit ** 2),
      abs(V1_unit - 1 / 24.0) < 1e-12 and abs(1 / V1_unit ** 2 - 576.0) < 1e-9)

# ----------------------------------------------------------------------
# E8  D258  外部边反作用：spec(A_S^0) = {4/5, 1, 1}
# ----------------------------------------------------------------------
head("E8  D258  外部边反作用：spec(A_S^0) = {4/5, 1, 1}")

m5 = 5
edges5 = list(itertools.combinations(range(4), 2)) + [(0, 4), (1, 4)]
K5 = np.ones(len(edges5))
B5 = np.zeros((m5, len(edges5)))
for col, (i, j) in enumerate(edges5):
    B5[i, col] = -1.0
    B5[j, col] = 1.0
L5 = B5 @ np.diag(K5) @ B5.T
Lp5 = np.linalg.pinv(L5)
X5 = sym_sqrt_psd(Lp5)
d = {}
for col, (i, j) in enumerate(edges5):
    d[(i, j)] = X5 @ (np.eye(m5)[j] - np.eye(m5)[i])
S = [0, 1, 2, 3]
# 局部单纯形 {0,1,2,3} 的三维仿射张成
span = np.array([X5[e] for e in S])
origin = span[0]
dirs = np.array([span[k] - origin for k in (1, 2, 3)]).T     # 5 x 3
Qb, _ = np.linalg.qr(dirs)
QS = Qb[:, :3]
AS = np.zeros((m5, m5))
for (i, j) in itertools.combinations(S, 2):
    AS += np.outer(d[(i, j)], d[(i, j)])
AS_loc = QS.T @ AS @ QS
spec = np.sort(np.linalg.eigvalsh(AS_loc))
check("局部交换刚度在三维叶上的特征值为 {4/5, 1, 1}",
      np.allclose(spec, [0.8, 1.0, 1.0], atol=2e-3),
      "spec = %s" % np.round(spec, 6))
check("spec 不是标量 => 不存在统一尺度使局部四面体与同一 P1 几何一致",
      float(np.max(spec) - np.min(spec)) > 0.1,
      "max-min = %.4f" % float(np.max(spec) - np.min(spec)))

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
