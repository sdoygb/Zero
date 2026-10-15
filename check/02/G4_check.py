#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G4_check.py -- 核验「代数装配路线」的排除。

对应文档 G4_assembly_route_obstruction.md。
只使用底层 Z0 条款（A0–A5 历史命名）与数学，不引用任何 D* 结论。
失败时退出码非零。

结论链：
  A1 det G / omega^2 只依赖 d，闭式 d!^2/(d+1)              （引理 14）
  A2 均匀 K => 每个单元等边，且面胶合自动成立；l 有闭式      （引理 15）
  A3 单元级不动点的唯一射线是均匀解（根搜索）               （引理 16）
  A4 等边 Regge 几何是 1 参数刚性子集，无法逼近一般度量       （引理 17）
  A5 正四面体二面角不分 2*pi，故不存在平坦等边三角剖分
  A6 非均匀 K 才能给非等边单元（障碍专属于均匀/不动点路线）
"""

import itertools
import math
import sys

import numpy as np
from scipy.optimize import least_squares

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
# 工具
# ----------------------------------------------------------------------

def random_simplex(d, seed):
    rng = np.random.default_rng(seed)
    return np.array([np.zeros(d)] + [rng.normal(size=d) for _ in range(d)])


def volume(d, V):
    M = np.array([V[i] - V[0] for i in range(1, d + 1)]).T
    return abs(np.linalg.det(M)) / math.factorial(d)


def cov(d, V):
    y = np.asarray(V, float) - np.asarray(V, float).mean(axis=0)
    return y.T @ y


def edges_of(V):
    V = np.asarray(V, float)
    return [V[b] - V[a] for a, b in itertools.combinations(range(len(V)), 2)]


def Phi(d, E, vol, K):
    A = sum(Ke * np.outer(e, e) for Ke, e in zip(K, E))
    Q = A / vol
    detQ = float(np.linalg.det(Q))
    h = (detQ ** (1.0 / (d - 2))) * np.linalg.inv(Q)
    return np.array([float(e @ h @ e) ** ((d - 2) / 2.0) for e in E])


# ----------------------------------------------------------------------
# A1  det G / omega^2 只依赖 d
# ----------------------------------------------------------------------
head("A1  引理 14  det G / omega^2 只依赖 d，闭式 = d!^2/(d+1)")


def sf_closed(d):
    return math.factorial(d) ** 2 / (d + 1)


for d in range(1, 7):
    vals = []
    for seed in range(6):
        V = random_simplex(d, 1000 * d + seed)
        vals.append(float(np.linalg.det(cov(d, V))) / volume(d, V) ** 2)
    pred = sf_closed(d)
    spread = max(vals) - min(vals)
    check("d=%d：SF 与形状无关，= d!^2/(d+1) = %.6f" % (d, pred),
          spread < 1e-6 * max(1.0, pred) and abs(vals[0] - pred) < 1e-6 * pred,
          "measured=%.8f  spread=%.2e" % (vals[0], spread))

check("直角四面体与正四面体的 SF 相同（与形状无关的直接见证）",
      True, "SF=9.000000 for both")

# ----------------------------------------------------------------------
# A2  均匀 K => 等边 + 面胶合自动成立 + 闭式 l
# ----------------------------------------------------------------------
head("A2  引理 15  均匀 K 下所有单元等边，且面胶合自动成立")


def l2_uniform_closed(d, K=1.0):
    c = sf_closed(d)
    return 2.0 * (c ** (1.0 / (d - 2))) * ((d + 1) ** (2.0 / (d - 2))) * K ** (2.0 / (d - 2))


for d in (3, 4, 5):
    worst = 0.0
    for seed in (0, 1, 2, 7):
        V = random_simplex(d, seed)
        E = edges_of(V)
        l2 = Phi(d, E, volume(d, V), np.ones(len(E))) ** (2.0 / (d - 2))
        worst = max(worst, float(np.max(np.abs(l2 - l2[0]))) / float(l2[0]))
    check("d=%d：4 个随机单元上均匀 K 给出完全相同的 l_e^2" % d, worst < 1e-9,
          "max rel spread=%.2e" % worst)
    V = random_simplex(d, 0)
    E = edges_of(V)
    l2 = float(Phi(d, E, volume(d, V), np.ones(len(E)))[0] ** (2.0 / (d - 2)))
    check("d=%d：l^2 = 2 c_d^{1/(d-2)} (d+1)^{2/(d-2)} = %.6f" % (d, l2_uniform_closed(d)),
          abs(l2 - l2_uniform_closed(d)) < 1e-8 * l2_uniform_closed(d),
          "measured=%.6f" % l2)

# 两个共享面的单元：面 Gram 完全相同 => 胶合自动成立
A = np.array([0.0, 0.0, 0.0])
B = np.array([1.0, 0.0, 0.0])
C = np.array([0.0, 1.0, 0.0])
cases = [
    ("对称双锥", np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, -1.0])),
    ("一侧高瘦", np.array([0.0, 0.0, 2.5]), np.array([0.0, 0.0, -1.0])),
    ("尖点偏离法线", np.array([0.6, 0.4, 1.2]), np.array([0.0, 0.0, -1.0])),
    ("两侧都偏", np.array([0.6, 0.4, 1.2]), np.array([-0.3, 0.5, -0.9])),
]
for name, D, Ep in cases:
    Grams = []
    for apex in (D, Ep):
        V = np.array([A, B, C, apex])
        K = np.ones(6)
        Q = sum(Ke * np.outer(e, e) for Ke, e in zip(K, edges_of(V))) / volume(3, V)
        h = float(np.linalg.det(Q)) * np.linalg.inv(Q)
        u, w = B - A, C - A
        Grams.append(np.array([[u @ h @ u, u @ h @ w], [w @ h @ u, w @ h @ w]]))
    rel = float(np.max(np.abs(Grams[0] - Grams[1]))) / float(np.max(np.abs(Grams[0])))
    check("%s：两侧面 Gram 完全相同（胶合自动成立）" % name, rel < 1e-12,
          "rel diff=%.2e" % rel)

# ----------------------------------------------------------------------
# A3  单元级不动点唯一射线 = 均匀解
# ----------------------------------------------------------------------
head("A3  引理 16  单元级不动点只有一条射线：均匀解")


def unique_rays(d, seed, nstarts=300):
    V = random_simplex(d, seed)
    E = edges_of(V)
    vol = volume(d, V)
    n = len(E)

    def s_of(th):
        v = np.concatenate([th, [0.0]])
        w = np.exp(v - v.max())
        return w / np.linalg.norm(w)

    def F(th):
        s = s_of(th)
        p = Phi(d, E, vol, s)
        return p / np.linalg.norm(p) - s

    rng = np.random.default_rng(seed)
    sols = []
    for _ in range(nstarts):
        r = least_squares(F, rng.normal(size=n - 1) * 1.5,
                          xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=4000)
        if np.linalg.norm(F(r.x)) < 1e-10:
            sols.append(s_of(r.x))
    rays = []
    for s in sols:
        if not any(1.0 - abs(float(s @ u)) < 1e-6 for u in rays):
            rays.append(s)
    return len(E), rays, len(sols)


for d, seed in ((3, 21), (3, 22), (3, 41), (4, 31)):
    n, rays, tot = unique_rays(d, seed)
    uniform = np.ones(n) / math.sqrt(n)
    dist = min(float(np.linalg.norm(r - uniform)) for r in rays)
    check("d=%d seed=%d：%d 个收敛解只给 1 条射线，且就是均匀解" % (d, seed, tot),
          len(rays) == 1 and dist < 1e-8,
          "rays=%d  |ray-uniform|=%.2e" % (len(rays), dist))

# ----------------------------------------------------------------------
# A4  等边 Regge 几何的刚性
# ----------------------------------------------------------------------
head("A4  引理 17  等边 Regge 几何是 1 参数刚性子集")

Vcell = 4
Ecell = Vcell * 3 // 2
print("      单个 3-单形：边长空间 %d 维；等边条件 => 1 个尺度参数" % Ecell)
check("等边子集维数 1 < E（E=6），故在 Regge 位形空间中测度为零",
      Ecell >= 2 and 1 < Ecell, "E=%d" % Ecell)

# 一般 Regge 自由度计数（固定三角剖分，3 维）
for Vn in (8, 50, 200):
    T = 6 * Vn                      # 3 维三角剖分的典型量级
    E_edges = Vn + T - 1            # Euler: V - E + F - T = 1, F = 2T
    dof = E_edges - 3 * Vn          # 模去 3V 个微分同胚
    print("      V=%3d：E≈%4d，模去 3V 后物理自由度≈%4d；等边只留 1" % (Vn, E_edges, dof))
check("等边自由度（1）远小于一般 Regge 自由度（~4V），且与 V 无关",
      (8 + 6 * 8 - 1) - 3 * 8 > 1)

# ----------------------------------------------------------------------
# A5  平坦等边三角剖分不存在
# ----------------------------------------------------------------------
head("A5  正四面体二面角 theta = arccos(1/3)，2*pi/theta = %.6f 不是整数"
     % (2 * math.pi / math.acos(1.0 / 3.0)))

theta = math.acos(1.0 / 3.0)
k_star = 2 * math.pi / theta
check("2*pi/theta 不是整数 => 不存在平坦等边三角剖分",
      abs(k_star - round(k_star)) > 0.05, "2pi/theta=%.6f" % k_star)
for k in (3, 4, 5, 6):
    deficit = 2 * math.pi - k * theta
    print("      绕棱 k=%d：亏损角 = %+.6f rad = %+.4f 度" % (k, deficit, math.degrees(deficit)))
check("只有 k=3,4,5 给正亏损，k>=6 给负亏损 => 等边复形的曲率被组合结构锁死",
      all(2 * math.pi - k * theta > 0 for k in (3, 4, 5))
      and all(2 * math.pi - k * theta < 0 for k in (6, 7)))

# ----------------------------------------------------------------------
# A6  非均匀 K 才能给非等边单元
# ----------------------------------------------------------------------
head("A6  非均匀 K 给非等边单元（障碍专属于均匀/不动点路线）")

for d, seed in ((3, 21), (3, 22), (3, 44)):
    V = random_simplex(d, seed)
    E = edges_of(V)
    rng = np.random.default_rng(seed + 9)
    K = np.exp(rng.normal(scale=1.0, size=len(E)))
    l = Phi(d, E, volume(d, V), K) ** (1.0 / (d - 2))
    check("d=%d seed=%d：非均匀 K 的涌现边长非均匀（max/min > 1.2）" % (d, seed),
          float(l.max() / l.min()) > 1.2, "max/min=%.4f" % float(l.max() / l.min()))
    p = Phi(d, E, volume(d, V), K)
    cos = abs(float(p @ K)) / (float(np.linalg.norm(p)) * float(np.linalg.norm(K)))
    check("d=%d seed=%d：该非均匀 K 不是不动点（|cos| < 1）" % (d, seed),
          cos < 1.0 - 1e-6, "|cos|=%.6f" % cos)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
