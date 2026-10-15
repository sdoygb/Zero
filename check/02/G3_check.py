#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G3_check.py -- 核验「导纳自洽不动点」(G2 引理 11) 的齐次度、显式解与稳定性。

对应文档 G3_admittance_fixed_point.md。
只使用底层条款（Z0 条款 ＋ Z1–Z5 定理）与数学，不引用任何 D* 结论。
失败时退出码非零。

单元模型：d 维单纯形，参考坐标顶点 {x_a}，边向量 e = x_b - x_a，
参考体积 omega = |det|/d!。给定权重 K > 0，
  A = sum_e K_e e_e (x) e_e,   Q = A/omega,
  h = (det Q)^{1/(d-2)} Q^{-1}          (d != 2；d=2 只有共形类)
  l_e = sqrt(e_e^T h e_e),     Phi(K)_e = l_e^{d-2}
不动点方程 K_e = kappa * l_e^{d-2}，等价于 Phi 的射线不动点。
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
# 单元模型
# ----------------------------------------------------------------------

def make_cell(d, seed, regular=False):
    if regular:
        verts = [np.zeros(d)] + [np.eye(d)[i] for i in range(d)]
    else:
        rng = np.random.default_rng(seed)
        verts = [np.zeros(d)] + [rng.normal(size=d) for _ in range(d)]
        M = np.array([verts[i] - verts[0] for i in range(1, d + 1)]).T
        if abs(np.linalg.det(M)) < 1e-3:
            return make_cell(d, seed + 1000)
    verts = np.array(verts)
    edges = [verts[b] - verts[a] for a, b in itertools.combinations(range(d + 1), 2)]
    M = np.array([verts[i] - verts[0] for i in range(1, d + 1)]).T
    vol = abs(np.linalg.det(M)) / math.factorial(d)
    return edges, vol


def Q_of_K(d, edges, vol, K):
    A = np.zeros((d, d))
    for Ke, e in zip(K, edges):
        A = A + Ke * np.outer(e, e)
    return A / vol


def h_of_Q(d, Q):
    if d == 2:
        raise ValueError("d=2 退化：只有共形类")
    detQ = float(np.linalg.det(Q))
    return (detQ ** (1.0 / (d - 2))) * np.linalg.inv(Q)


def Phi(d, edges, vol, K):
    Q = Q_of_K(d, edges, vol, K)
    h = h_of_Q(d, Q)
    return np.array([float(e @ h @ e) ** ((d - 2) / 2.0) for e in edges])


def normalize(v):
    return v / float(np.linalg.norm(v))


def solve_ray(d, edges, vol, K0, steps=2000, damp=0.5):
    K = normalize(np.abs(K0))
    for _ in range(steps):
        K = normalize((1.0 - damp) * K + damp * normalize(Phi(d, edges, vol, K)))
    return K


# ----------------------------------------------------------------------
# E1  齐次度：Phi(tK) = t Phi(K)，对每个 d 都是 1 次
# ----------------------------------------------------------------------
head("E1  引理 12  Phi 在 K 上齐次度为 1（对所有 d）")

degs = []
for d in (3, 4, 5, 6):
    edges, vol = make_cell(d, 7)
    rng = np.random.default_rng(3)
    K = rng.uniform(0.5, 2.0, size=len(edges))
    base = Phi(d, edges, vol, K)
    local = []
    for t in (0.25, 2.3, 7.0):
        deg = math.log(
            float(np.linalg.norm(Phi(d, edges, vol, t * K)) / np.linalg.norm(base))
        ) / math.log(t)
        local.append(deg)
        check("d=%d  t=%.2f：齐次度 = 1" % (d, t), abs(deg - 1.0) < 1e-8,
              "measured=%.10f" % deg)
    degs.append((d, local))

check("负结果：齐次度在 d=3,4,5,6 上全部为 1 => 齐次度不选择维数",
      all(abs(x - 1.0) < 1e-8 for _, loc in degs for x in loc))

head("E1b  h 在 K 上的尺度：h(tK) = t^{2/(d-2)} h(K)")

for d in (3, 4, 5):
    edges, vol = make_cell(d, 11)
    rng = np.random.default_rng(5)
    K = rng.uniform(0.5, 2.0, size=len(edges))
    h1 = h_of_Q(d, Q_of_K(d, edges, vol, K))
    t = 3.0
    h2 = h_of_Q(d, Q_of_K(d, edges, vol, t * K))
    pred = t ** (2.0 / (d - 2))
    ratio = float(np.linalg.norm(h2) / np.linalg.norm(h1))
    check("d=%d  h(tK)/h(K) = t^{2/(d-2)} = %.4f" % (d, pred),
          abs(ratio / pred - 1.0) < 1e-6, "measured=%.6f" % ratio)

# ----------------------------------------------------------------------
# E2  参考坐标缩放下 Phi 不变（l 是坐标不变量）
# ----------------------------------------------------------------------
head("E2  参考坐标缩放下 Phi 不变（l 是坐标不变量）")

for d in (3, 4):
    edges, vol = make_cell(d, 13)
    rng = np.random.default_rng(9)
    K = rng.uniform(0.5, 2.0, size=len(edges))
    base = Phi(d, edges, vol, K)
    for s in (0.3, 4.0):
        got = Phi(d, [s * e for e in edges], vol * s ** d, K)
        check("d=%d  s=%.1f：Phi 不变" % (d, s),
              float(np.linalg.norm(got - base) / np.linalg.norm(base)) < 1e-9)

# ----------------------------------------------------------------------
# E3  单纯形协方差恒等式（引理 13 的证明核心）
# ----------------------------------------------------------------------
head("E3  引理 13  单纯形协方差恒等式：白化后成对距离恒为 2")

for d in (2, 3, 4, 5):
    for seed in (0, 1, 2):
        rng = np.random.default_rng(100 * d + seed)
        V = np.array([np.zeros(d)] + [rng.normal(size=d) for _ in range(d)])
        y = V - V.mean(axis=0)
        G = y.T @ y
        Gi = np.linalg.inv(G)
        M = np.array([[float(y[a] @ Gi @ y[b]) for b in range(d + 1)] for a in range(d + 1)])
        dists = [float((V[b] - V[a]) @ Gi @ (V[b] - V[a]))
                 for a, b in itertools.combinations(range(d + 1), 2)]
        diag = np.diag(M)
        off = M[~np.eye(d + 1, dtype=bool)]
        ok = (max(dists) - min(dists) < 1e-10
              and diag.max() - diag.min() < 1e-10
              and off.max() - off.min() < 1e-10
              and abs(float(diag[0]) - d / (d + 1)) < 1e-10
              and abs(float(off[0]) + 1.0 / (d + 1)) < 1e-10)
        check("d=%d seed=%d：对角 = d/(d+1)，非对角 = -1/(d+1)，成对距离 = 2" % (d, seed),
              ok, "dist spread=%.1e" % (max(dists) - min(dists)))

# ----------------------------------------------------------------------
# E4  均匀 K 是任意单纯形的射线不动点
# ----------------------------------------------------------------------
head("E4  引理 13 续  均匀 K 对任意单元形状都是射线不动点")

for d in (3, 4, 5):
    worst = 0.0
    for seed in (0, 1, 2, 3, 4, 5):
        edges, vol = make_cell(d, seed)
        phi = Phi(d, edges, vol, np.ones(len(edges)))
        ratio = phi / phi[0]
        worst = max(worst, float(np.max(np.abs(ratio - 1.0))))
    check("d=%d  6 个随机单元上 Phi(uniform) ∝ uniform" % d, worst < 1e-10,
          "max rel spread=%.2e" % worst)

edges, vol = make_cell(3, 0, regular=True)
phi1 = Phi(3, edges, vol, np.ones(6))
check("d=3  正四面体：l_e = sqrt(288) = %.4f（解析值）" % math.sqrt(288.0),
      abs(float(phi1[0]) - math.sqrt(288.0)) < 1e-8, "measured=%.8f" % float(phi1[0]))

# ----------------------------------------------------------------------
# E5  迭代收敛：随机初值 -> 均匀射线（唯一吸引子）
# ----------------------------------------------------------------------
head("E5  引理 13 续  迭代收敛到均匀射线，与单元形状无关")

for d, seed in ((3, 21), (3, 22), (4, 31)):
    edges, vol = make_cell(d, seed)
    E = len(edges)
    target = normalize(np.ones(E))
    rng = np.random.default_rng(seed)
    starts = [np.abs(rng.normal(size=E)) + 0.1 for _ in range(24)]
    limits = [solve_ray(d, edges, vol, s) for s in starts]
    residuals = [float(np.linalg.norm(normalize(Phi(d, edges, vol, L)) - L)) for L in limits]
    dist_to_uniform = max(float(np.linalg.norm(L - target)) for L in limits)
    check("d=%d seed=%d：射线残差 < 1e-8（不动点存在）" % (d, seed),
          max(residuals) < 1e-8, "max residual=%.3e" % max(residuals))
    check("d=%d seed=%d：24 个初值全部收敛到均匀射线" % (d, seed),
          dist_to_uniform < 1e-8, "max |L - uniform|=%.3e" % dist_to_uniform)
    Q = Q_of_K(d, edges, vol, limits[0])
    check("d=%d seed=%d：不动点处 Q 正定" % (d, seed),
          float(np.min(np.linalg.eigvalsh(Q))) > 0)

# ----------------------------------------------------------------------
# E6  局部稳定性：均匀不动点附近的扰动衰减
# ----------------------------------------------------------------------
head("E6  均匀不动点的局部稳定性（扰动衰减）")

for d, seed in ((3, 21), (3, 23), (4, 31)):
    edges, vol = make_cell(d, seed)
    E = len(edges)
    s0 = normalize(np.ones(E))
    rng = np.random.default_rng(seed + 500)
    factors = []
    for _ in range(8):
        pert = rng.normal(size=E)
        pert = pert - (pert @ s0) * s0
        pert = pert / float(np.linalg.norm(pert))
        K = normalize(s0 + 1e-2 * pert)
        for _ in range(150):
            K = normalize(0.5 * K + 0.5 * normalize(Phi(d, edges, vol, K)))
        factors.append(float(np.linalg.norm(K - s0)))
    check("d=%d seed=%d：8 个扰动方向全部衰减（局部稳定）" % (d, seed),
          max(factors) < 1e-3, "max final deviation=%.3e" % max(factors))

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
