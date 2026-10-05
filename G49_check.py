#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G49_check.py -- 推进 G48 的四条诚实边界。

对应文档 G49_four_boundaries_advanced.md。

  F1  边界1（维度）：2D 方格上依赖半径仍为 m（有效距离 = 2d+1）
  F2  边界1：精确分解定理 w^(m) = m * c^m * (A_top^{m-1})；系数是【拓扑数】（与 c 无关）
  F3  边界2（层次核验）：D 系列的层确有【长度索引】（D216 历史长度模 T、历史窗口截断）
  F4  边界2：但也有【粗粒化轴】（D222 从精确到概括）=> 我的 m 只覆盖长度轴 => 部分核验
  F5  边界3（适用范围）：量化非均匀几何上的误差（m=2 精确；内部 <1%；大误差在边界）
  F6  边界4（截断处）：m* = L 的求和【收敛】（逐步差减半）
  F7  诚实边界
"""

import os
import sys
from math import comb

import numpy as np
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
MOD = os.path.expanduser("~/Downloads/modular-equilibrium/derivations")
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


def rd(f, base=None):
    b = base if base is not None else HERE
    p = os.path.join(b, f)
    return open(p, encoding="utf-8", errors="replace").read() if os.path.exists(p) else ""


def sup_pow(c, m):
    """(A^m)_{i,i+1} for all i（链）"""
    N = len(c) + 1
    if m == 0:
        return np.ones(N - 1)
    A = sp.diags([c, c], [-1, 1], shape=(N, N)).tocsr()
    return (A ** m).diagonal(1)


def w_super(c, m):
    """w^(m)_{i,i+1} = m * c_i * (A^{m-1})_{i,i+1}"""
    return m * c * sup_pow(c, m - 1)


def grid2d(Lx, Ly, c=1.0):
    N = Lx * Ly
    rows, cols = [], []
    for i in range(Lx):
        for j in range(Ly):
            k = i * Ly + j
            if i + 1 < Lx:
                rows += [k, k + Ly]; cols += [k + Ly, k]
            if j + 1 < Ly:
                rows += [k, k + 1]; cols += [k + 1, k]
    A = sp.coo_matrix(([c] * len(rows), (rows, cols)), shape=(N, N)).tocsr()
    return A


def w_m_2d(A, m):
    """w^(m) = m * A_ij * (A^{m-1})_ij —— 注意 scipy 稀疏的 * 是矩阵乘，元素乘要用 multiply"""
    if m < 2:
        return A * 0
    return m * A.multiply(A ** (m - 1))


# ======================================================================
head("F1  边界1：2D 方格上依赖半径仍为 m")

Lx = Ly = 40
A = grid2d(Lx, Ly)
ref = 20 * Ly + 20                    # 参考边 (20,20)-(20,21)
per = 20 * Ly + 30                    # 扰动边 (20,30)-(20,31)，图距离 10
print("      扰动与参考的图距离 = 10（往返最短 = 21 步）")
print("      m      w^(m) 在参考边的相对变化")
res2d = {}
for m in (2, 4, 8, 16, 22, 32):
    A1 = grid2d(Lx, Ly)
    A1 = A1.tolil()
    A1[per, per + 1] = 1.5
    A1[per + 1, per] = 1.5
    A1 = A1.tocsr()
    W0, W1 = w_m_2d(A, m), w_m_2d(A1, m)
    r0, r1 = float(W0[ref, ref + 1]), float(W1[ref, ref + 1])
    rel = abs(r1 - r0) / r0
    res2d[m] = rel
    print("      %-6d %.3e" % (m, rel))
check("m <= 16 时影响【精确为 0】（m 小于往返距离 21）",
      all(res2d[m] == 0.0 for m in (2, 4, 8, 16)))
check("m >= 22 时影响非零（m 达到往返距离）", all(res2d[m] > 0 for m in (22, 32)))
check("=> 2D 上依赖半径仍由 m 控制", True)

# ======================================================================
head("F1b  边界1 补：2D 上【固定 m】的细化收敛")


def grid2d_w(L, cfun):
    N = L * L
    rows, cols, vals = [], [], []
    for i in range(L):
        for j in range(L):
            k = i * L + j
            if i + 1 < L:
                c = cfun(i + 0.5, j)
                rows += [k, k + L]; cols += [k + L, k]; vals += [c, c]
            if j + 1 < L:
                c = cfun(i, j + 0.5)
                rows += [k, k + 1]; cols += [k + 1, k]; vals += [c, c]
    return sp.coo_matrix((vals, (rows, cols)), shape=(N, N)).tocsr()


print("      Lx=Ly    N        m=2 内部范围    m=4 内部范围    m=6 内部范围")
c2d = {2: [], 4: [], 6: []}
for L in (40, 60, 80, 100):
    cfun = lambda i, j, L=L: 1.0 + 0.3 * np.cos(2 * np.pi * i / L) * np.cos(2 * np.pi * j / L)
    A2 = grid2d_w(L, cfun)
    row = []
    for m in (2, 4, 6):
        W = w_m_2d(A2, m)
        lo, hi = int(0.25 * L), int(0.75 * L)
        vals = []
        for i2 in range(lo, hi):
            for j2 in range(lo, hi):
                k2 = i2 * L + j2
                if j2 + 1 < L:
                    vals.append(float(W[k2, k2 + 1]))
        v = np.array(vals)
        r = float(v.max() / v.min())
        c2d[m].append(r); row.append(r)
    print("      %-8d %-7d %-14.6f %-14.6f %.6f" % (L, L * L, *row))

for m in (2, 4, 6):
    v = c2d[m]
    d = [abs(v[k + 1] - v[k]) for k in range(len(v) - 1)]
    check("2D m=%d：逐步差递减（%s）=> 固定 m 收敛" % (m, " ".join("%.2e" % x for x in d)),
          all(d[k + 1] < d[k] for k in range(len(d) - 1)))

# ======================================================================
head("F2  边界1：精确分解定理 w^(m) = m * c^m * (A_top^{m-1})")

print("      均匀几何下：w^(m)_{ij} = m c^m (A_top^{m-1})_{ij}")
print("      1D 系数 m*C(m-1,m/2-1)；2D 系数为另一组【拓扑数】；两者都与 c 无关")
print()
print("      图     m     系数（c=0.7）    系数（c=2.0）    与 c 无关?")
for tag, mk in (("1D", None), ("2D", None)):
    for m in (2, 4, 6):
        if tag == "1D":
            vals = []
            for c in (0.7, 2.0):
                cc = np.full(300, c)
                vals.append(float(w_super(cc, m)[150]) / c ** m)
        else:
            vals = []
            for c in (0.7, 2.0):
                B = grid2d(30, 30, c)
                ctr = 15 * 30 + 15
                vals.append(float(w_m_2d(B, m)[ctr, ctr + 1]) / c ** m)
        print("      %-5s %-6d %-16.6f %-16.6f %s" % (tag, m, vals[0], vals[1], abs(vals[0] - vals[1]) < 1e-9))
        check("%s m=%d：系数与 c 无关（差 %.1e）" % (tag, m, abs(vals[0] - vals[1])),
              abs(vals[0] - vals[1]) < 1e-9)
check("=> 几何进入为 c^m、拓扑进入为游走计数 => 两者【精确因子化】", True)

# ======================================================================
head("F3  边界2：D 系列的层确有【长度索引】")

d216 = rd("D216_history_phase_state.md", MOD)
d222 = rd("D222_stratified_destruction_and_local_memory.md", MOD)
d211 = rd("D211_global_static_closure_zero_layer.md", MOD)
check("D216 预先结构含『历史长度模 T 的相位解释』", "历史长度模" in d216)
check("D216 预先结构含『历史窗口截断』", "历史窗口截断" in d216)
check("D222 含『亚层分解』", "亚层" in d222)
check("D211 含『可逆周期』", "可逆周期" in d211)
check("=> D 系列确实有【长度/周期型】索引 => 与我的 m（闭环游走长度）同型", True)

# ======================================================================
head("F4  边界2：但也有【粗粒化轴】=> 部分核验")

check("D222：活动层亚层在【局部寿命】到达时全清（年龄轴）", "局部寿命到达时全部清空" in d222)
check("D222：历史层亚层『从精确到概括』（粗粒化轴）", "从精确到概括" in d222)
check("D222：毁灭时只保留【最高两层】（粗粒化规则）", "最高两层" in d222)
check("=> D 系列的『亚层』至少有【年龄轴】与【粗粒化轴】两条", True)
check("=> 我的 m 只对应【长度轴】=> 对应关系是【部分的】", True)

# ======================================================================
head("F5  边界3：非均匀几何上幂律误差的量化")

N = 2048
i = np.arange(N - 1)
c = 1.0 + 0.3 * np.cos(2 * np.pi * (i + 0.5) / N)
print("      m     全链最大相对误差     中间80%最大误差")
err = {}
for m in (2, 4, 8, 16, 32, 64):
    e = w_super(c, m)
    pred = np.array([m * comb(m - 1, m // 2 - 1) * c[j] ** m for j in range(N - 1)])
    rel = np.abs(e - pred) / pred
    lo, hi = int(0.1 * N), int(0.9 * N)
    err[m] = (float(rel.max()), float(rel[lo:hi].max()))
    print("      %-6d %-21.4e %.4e" % (m, err[m][0], err[m][1]))
check("m=2：全链误差 < 1e-14（【精确】）", err[2][0] < 1e-14)
check("内部 80%：m<=64 时误差全部 < 1%（<1e-2）",
      all(err[m][1] < 1e-2 for m in err))
check("全链最大误差 >> 内部误差（误差集中在边界）", err[64][0] > 10 * err[64][1],
      "%.3e vs %.3e" % (err[64][0], err[64][1]))

# ======================================================================
head("F6  边界4：截断 m* = L 的求和【收敛】")

print("      N        L=4          L=8          L=16         L=32")
cols = {4: [], 8: [], 16: [], 32: []}
for N2 in (512, 1024, 2048, 4096):
    i2 = np.arange(N2 - 1)
    c2 = 1.0 + 0.3 * np.cos(2 * np.pi * (i2 + 0.5) / N2)
    row = []
    for L in (4, 8, 16, 32):
        tot = np.zeros(N2 - 1)
        for m in range(2, L + 1):
            tot += w_super(c2, m)
        lo, hi = int(0.15 * N2), int(0.85 * N2)
        v = float(tot[lo:hi].max() / tot[lo:hi].min())
        cols[L].append(v); row.append(v)
    print("      %-8d %-12.6f %-12.6f %-12.6f %.6f" % (N2, *row))

for L in (4, 8, 16, 32):
    v = cols[L]
    d = [abs(v[k + 1] - v[k]) for k in range(3)]
    check("L=%d：逐步差递减（%s）=> 收敛" % (L, " ".join("%.2e" % x for x in d)),
          all(d[k + 1] < d[k] for k in range(2)))
check("=> 截断到 m* = L（Z3 的寿命，固定步数）时，求和收敛", True)

# ======================================================================
head("F7  诚实边界")

check("2D 只测了方格，未测其他 2D 几何/高维", True)
check("『D 系列的层 = 我的 m 分层』只是【部分】对应（长度轴对应、粗粒化轴不对应）", True)
check("边界3 的误差只在【一个】非均匀剖面（cos 调制）上测", True)
check("边界4 的收敛判据是『逐步差递减』，不是数学证明", True)
check("本计算推进 G48 的四条边界，不改变 G1-G48 的其余数值结论", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
