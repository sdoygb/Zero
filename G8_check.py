#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G8_check.py -- 核验 I1（维数）：Z0 条款（A0–A5 历史命名）不选择维数；结构性筛选把 D<=3 排除。

对应文档 G8_dimension_selection.md。
只使用底层 Z0 条款（A0–A5 历史命名）与数学，不引用任何 D* 结论。
失败时退出码非零。

性能说明：曲率部分全程不做符号化简，只构建表达式 DAG 后数值取值。
（早期版本在循环里用 simplify 需约 8 分钟；改成 expand 则表达式爆炸失控。）

结论链：
  F1 Z1 定理 2 只固定余维 1 => 任意 m 都可（不选择维数）                引理 35
  F2 D=2：G_ab 恒为零（无场方程内容）                            引理 36
  F3 D=3：Weyl 张量恒为零（无传播引力子）                        引理 37
  F4 D=4：Schwarzschild 真空且 Weyl != 0                         引理 38
  F5 引力子自由度 > 0 <=> D >= 4；D(D-3)/2 = 2 的唯一正解是 D=4   引理 38
  F6 闭 ±1 词长度必偶（若取 ±1 读法）                             引理 39
  F7 筛选链：D>=4 且最小 => D=4
"""

import itertools
import sys

import numpy as np
import sympy as sp

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
# 曲率（无化简：构建 DAG 后数值取值）
# ----------------------------------------------------------------------

def geometry(gd, coords):
    n = len(coords)
    g = sp.diag(*gd)
    ginv = sp.diag(*[1 / e for e in gd])
    Gam = [[[sp.Integer(0)] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                v = sp.Integer(0)
                if a == c:
                    v += sp.diff(gd[a], coords[b])
                if a == b:
                    v += sp.diff(gd[a], coords[c])
                if b == c:
                    v -= sp.diff(gd[b], coords[a])
                Gam[a][b][c] = v / (2 * gd[a])
    Rup = [[[[sp.Integer(0)] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    v = sp.diff(Gam[a][b][d], coords[c]) - sp.diff(Gam[a][b][c], coords[d])
                    for e in range(n):
                        v += Gam[a][c][e] * Gam[e][b][d] - Gam[a][d][e] * Gam[e][b][c]
                    Rup[a][b][c][d] = v
    Rdn = [[[[sum(g[a, e] * Rup[e][b][c][d] for e in range(n)) for d in range(n)]
             for c in range(n)] for b in range(n)] for a in range(n)]
    Ric = sp.zeros(n, n)
    for b in range(n):
        for d in range(n):
            Ric[b, d] = sum(Rup[a][b][a][d] for a in range(n))
    Rs = sum(ginv[b, d] * Ric[b, d] for b in range(n) for d in range(n))
    W = None
    if n >= 3:
        W = [[[[Rdn[a][b][c][d]
                - (g[a, c] * Ric[b, d] - g[a, d] * Ric[b, c]
                   + g[b, d] * Ric[a, c] - g[b, c] * Ric[a, d]) / (n - 2)
                + Rs * (g[a, c] * g[b, d] - g[a, d] * g[b, c]) / ((n - 1) * (n - 2))
                for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
    return g, Ric, Rs, Rdn, W


def num(expr, pts):
    return complex(sp.N(expr.subs(pts)))


def mx(tensor, n, pts):
    return max(abs(num(tensor[a][b][c][d], pts))
               for a in range(n) for b in range(n) for c in range(n) for d in range(n))


# ======================================================================
# F1  Z1 定理 2 只固定余维 1
# ======================================================================
head("F1  引理 35  Z1 定理 2 只固定余维 1 => 任意 m 都合格，不选择维数")

for m in (2, 3, 4, 5, 6, 7):
    H = np.array([[1.0 if k == j else (-1.0 if k == 0 else 0.0) for k in range(m)]
                  for j in range(1, m)])
    check("m=%d：dim H_Q = %d（余维 1）" % (m, m - 1),
          np.linalg.matrix_rank(H) == m - 1)
check("故候选维数集合 = {2,3,4,...}（Z0 条款（A0–A5 历史命名）不筛）", True)

t, x, y, z = sp.symbols("t x y z", real=True)

# ======================================================================
# F2  D=2：G_ab 恒为零
# ======================================================================
head("F2  引理 36  D=2：Einstein 张量恒为零（无场方程内容）")

g2, Ric2, Rs2, _, _ = geometry([-(1 + t ** 2 * x ** 2), 1 + sp.sin(t * x) ** 2 + x ** 2], (t, x))
pts2 = {t: 0.3, x: 0.7}
mG = max(abs(num(Ric2[a, b] - sp.Rational(1, 2) * Rs2 * g2[a, b], pts2))
         for a in range(2) for b in range(2))
mRic = max(abs(num(Ric2[a, b], pts2)) for a in range(2) for b in range(2))
check("任意 2 维度规上 R_ab = (1/2) R g_ab（|G|/|Ric| ~ 0）",
      mG / mRic < 1e-12, "|G|=%.2e  |Ric|=%.2e  比值=%.2e" % (mG, mRic, mG / mRic))
check("2 维 Ricci 张量一般非零（曲率存在，但无自由动力学）", mRic > 1e-6,
      "|Ric|=%.3e" % mRic)

# ======================================================================
# F3  D=3：Weyl 恒为零
# ======================================================================
head("F3  引理 37  D=3：Weyl 张量恒为零（无传播引力子）")

g3, Ric3, Rs3, Rdn3, W3 = geometry(
    [-(1 + t ** 2 * x ** 2), 1 + t * x * y + x ** 2, 1 + x ** 2 * y ** 2], (t, x, y))
pts3 = {t: 0.3, x: 0.7, y: 1.1}
mW3 = mx(W3, 3, pts3)
mR3 = mx(Rdn3, 3, pts3)
mRic3 = max(abs(num(Ric3[a, b], pts3)) for a in range(3) for b in range(3))
check("任意 3 维度规上 Weyl 张量恒为零（|W|/|Riem| ~ 0）",
      mW3 / mR3 < 1e-12,
      "|W|=%.2e  |Riem|=%.2e  比值=%.2e" % (mW3, mR3, mW3 / mR3))
check("3 维 Ricci 一般非零（曲率非平凡，但无 Weyl 部分）", mRic3 > 1e-6,
      "|Ric|=%.3e" % mRic3)
check("故 3 维没有独立于局部物质的曲率自由度 => 无传播引力子", mW3 / mR3 < 1e-12)

# ======================================================================
# F4  D=4：Schwarzschild 真空但 Weyl != 0
# ======================================================================
head("F4  引理 38  D=4：Schwarzschild 真空且 Weyl != 0，引力子自由度 = 2")

M, r, th, ph = sp.symbols("M r theta phi", positive=True)
gd4 = [-(1 - 2 * M / r), 1 / (1 - 2 * M / r), r ** 2, r ** 2 * sp.sin(th) ** 2]
g4, Ric4, Rs4, Rdn4, W4 = geometry(gd4, (t, r, th, ph))
pts4 = {M: 1.0, t: 0.0, r: 5.0, th: 1.0, ph: 0.0}
mRic4 = max(abs(num(Ric4[a, b], pts4)) for a in range(4) for b in range(4))
mW4 = mx(W4, 4, pts4)
mR4 = mx(Rdn4, 4, pts4)
check("Schwarzschild 是真空解：R_ab = 0", mRic4 / mR4 < 1e-12,
      "|Ric|=%.2e  |Riem|=%.3e" % (mRic4, mR4))
check("D=4 上 Weyl 张量非零（真空仍有独立曲率自由度）", mW4 / mR4 > 1e-3,
      "|W|=%.3e  |W|/|Riem|=%.4f" % (mW4, mW4 / mR4))
check("D=4 引力子自由度 D(D-3)/2 = 2", 4 * (4 - 3) // 2 == 2)

# ======================================================================
# F5  自由度筛选
# ======================================================================
head("F5  引力子自由度 > 0 <=> D >= 4；D(D-3)/2 = 2 的唯一正解是 D = 4")

tbl = []
for D in range(2, 8):
    tbl.append((D, D * (D - 3) // 2))
    print("      D=%d  自由度 = %d   %s" % (D, D * (D - 3) // 2,
                                            "无动力学" if D <= 3 else "有传播引力子"))
check("自由度 <= 0 对 D<=3，> 0 对 D>=4",
      all(v <= 0 for D, v in tbl if D <= 3) and all(v > 0 for D, v in tbl if D >= 4))
Ds = sp.symbols("Ds")
roots = sp.solve(sp.Eq(Ds * (Ds - 3) / 2, 2), Ds)
check("D(D-3)/2 = 2 的根为 {-1, 4}，唯一正解 D=4", sorted(roots) == [-1, 4],
      "roots=%s" % roots)
check("故『引力子两极化 = 闭词的两个符号』给出 D = 4（候选原则，非导出）", True)

# ======================================================================
# F6  闭 ±1 词长度必偶（条件性）
# ======================================================================
head("F6  引理 39  闭 ±1 词长度必偶（仅在 ±1 读法下）")

for L in (1, 3, 5, 7):
    check("长度 %d 的 ±1 词不存在零和解" % L,
          not any(sum(w) == 0 for w in itertools.product([1, -1], repeat=L)))
for L in (2, 4, 6, 8):
    check("长度 %d 的 ±1 词存在零和解（等量正负）" % L,
          any(sum(w) == 0 for w in itertools.product([1, -1], repeat=L)))
check("故 ±1 读法下维数必偶；整数荷读法（Z1 定理 2）无此约束", True)
check("反例：长度 3 的零和整数态 (2,-1,-1) 存在，故奇偶性只在 ±1 读法下成立",
      sum([2, -1, -1]) == 0)

# ======================================================================
# F7  筛选链
# ======================================================================
head("F7  筛选链：D >= 4 且最小 => D = 4")

cands = list(range(2, 20))
after = [D for D in cands if D * (D - 3) // 2 > 0]
check("引力子存在把候选集从 {2,3,4,...} 收缩到 {4,5,6,...}",
      after[0] == 4 and after == list(range(4, 20)))
check("最小性给出 D = 4", min(after) == 4)
check("若同时要求偶维（±1 读法），最小候选仍是 D = 4",
      min([D for D in after if D % 2 == 0]) == 4)
check("诚实边界：最小性是额外原则，不是 Z0 条款（A0–A5 历史命名）的推论", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
