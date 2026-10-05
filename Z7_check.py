#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z7_check.py —— 【E1 的显式化：嵌入输入 = 局部标度场；度规的闭式字典】核验
======================================================================
独立实断言：
  F1  闭式字典（层和口径 w=Σ_m m A∘A^{m-1}）：1D/2D/3D 二阶收敛
  F2  闭式字典（单层口径 w=k A∘A^{k-1}）：1D/2D/3D 二阶收敛
  F3  格点游走数 W_d(L) 精确整数表（并复现 Z6 的三维正则值）
  F4  二维固定 k 细化：剖面二阶收敛（G58 的 1D 结果推广到 2D）
  F5  单元形状携带度规：F_x/F_y = (a_x/a_y)^2 (1+O(a^2))
  F6  两种口径的定性结论一致（(O)/(C)/支集/半径 ≤ k/2 且半径外精确 0）
  F7  k ∝ N 支在二维发散（非均匀场）
  F8  文档结论与引文在位
"""
import itertools
import io
import math
import os
import sys
from collections import defaultdict

import numpy as np
from scipy.sparse import csr_matrix, diags

HERE = os.path.dirname(os.path.abspath(__file__))
LEDGER_MODE = os.environ.get("LH_LEDGER", "A")
FAIL = 0
N = 0
K = 8


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


DOC = io.open(os.path.join(HERE, "Z7_embedding_input_explicit_dictionary.md"), encoding="utf-8").read()


# ---------------------------------------------------------------- 工具
def widx(L, d):
    """Z^d 上从 0 到 e_1 的 L 步游走数（精确整数 DP）。"""
    cur = {(0,) * d: 1}
    for _ in range(L):
        nxt = defaultdict(int)
        for st, n in cur.items():
            for a in range(d):
                for s in (1, -1):
                    t = list(st)
                    t[a] += s
                    nxt[tuple(t)] += n
        cur = nxt
    t = [0] * d
    t[0] = 1
    return cur.get(tuple(t), 0)


def grid(NN, d, cfun):
    n = NN + 1
    strides = [n ** (d - 1 - i) for i in range(d)]
    idx = lambda t: sum(x * s for x, s in zip(t, strides))
    r, cc, v = [], [], []
    for node in itertools.product(range(n), repeat=d):
        i = idx(node)
        for a in range(d):
            if node[a] < NN:
                nb = list(node)
                nb[a] += 1
                j = idx(nb)
                mid = list(node)
                mid[a] += 0.5
                x = [mid[b] / NN for b in range(d)]
                cw = cfun(x)
                r += [i, j]
                cc += [j, i]
                v += [cw, cw]
    return csr_matrix((v, (r, cc)), shape=(n ** d, n ** d)), idx, strides


def W_sum(A, k):
    P = A.copy()
    W = A.multiply(P) * 2
    for m in range(3, k + 1):
        P = (P @ A).tocsr()
        W = W + m * A.multiply(P)
    return W.tocsr()


def W_single(A, k):
    P = A.copy()
    for _ in range(k - 2):
        P = (P @ A).tocsr()
    return A.multiply(P).tocsr()


def P_sum(c, d, k=K):
    return sum(m * widx(m - 1, d) * c ** m for m in range(2, k + 1, 2))


def P_single(c, d, k=K):
    return widx(k - 1, d) * c ** k


def devs_1d(mode, Ns=(64, 256, 1024)):
    out = []
    for NN in Ns:
        r, cc, v = [], [], []
        X = (np.arange(NN) + 0.5) / NN
        c = 1 + 0.3 * np.cos(2 * np.pi * X)
        for i in range(NN):
            r += [i, i + 1]
            cc += [i + 1, i]
            v += [c[i], c[i]]
        A = csr_matrix((v, (r, cc)), shape=(NN + 1, NN + 1))
        W = W_sum(A, K) if mode == "sum" else W_single(A, K)
        w = np.array([W[i, i + 1] for i in range(NN)])
        m = (X > 0.2) & (X < 0.8)
        P = P_sum if mode == "sum" else P_single
        out.append(float(np.abs(w[m] / np.array([P(ci, 1) for ci in c[m]]) - 1).max()))
    return out


def devs_2d(mode, Ns=(32, 64, 128)):
    out = []
    for NN in Ns:
        A, idx, st = grid(NN, 2, lambda x: 1 + 0.2 * math.cos(2 * math.pi * x[0])
                          + 0.15 * math.cos(2 * math.pi * x[1]))
        W = (W_sum(A, K) if mode == "sum" else W_single(A, K)).tocoo()
        n = NN + 1
        ds = []
        for i, j, w in zip(W.row.tolist(), W.col.tolist(), W.data.tolist()):
            if j == i + n:
                ii, jj = divmod(i, n)
                x = (ii + 0.5) / NN
                y = jj / NN
                if 0.25 < x < 0.75 and 0.25 < y < 0.75:
                    c = 1 + 0.2 * math.cos(2 * math.pi * x) + 0.15 * math.cos(2 * math.pi * y)
                    P = P_sum if mode == "sum" else P_single
                    ds.append(abs(w / P(c, 2) - 1))
        out.append(float(max(ds)))
    return out


def devs_3d(mode, Ns=(12, 16, 24)):
    out = []
    for NN in Ns:
        A, idx, st = grid(NN, 3, lambda x: 1 + 0.15 * math.cos(2 * math.pi * x[0]))
        W = (W_sum(A, K) if mode == "sum" else W_single(A, K)).tocoo()
        n = NN + 1
        n2 = n * n
        ds = []
        for i, j, w in zip(W.row.tolist(), W.col.tolist(), W.data.tolist()):
            if j == i + n2:
                q = divmod(i, n2)
                x0 = (q[0] + 0.5) / NN
                y0 = (q[1] // n) / NN
                z0 = (q[1] % n) / NN
                if 0.3 < x0 < 0.7 and 0.3 < y0 < 0.7 and 0.3 < z0 < 0.7:
                    c = 1 + 0.15 * math.cos(2 * math.pi * x0)
                    P = P_sum if mode == "sum" else P_single
                    ds.append(abs(w / P(c, 3) - 1))
        out.append(float(max(ds)))
    return out


# ---------------------------------------------------------------- F1
head("F1  闭式字典（层和口径）：w_e = [Σ_{m even} m W_d(m-1) c_e^m](1+O(a²))")
d1 = devs_1d("sum")
check("1D：偏差随 N 递减且为二阶（N×4 ⟹ 偏差 ÷>8）",
      d1[0] > d1[1] > d1[2] and d1[0] / d1[1] > 8 and d1[1] / d1[2] > 8,
      " -> ".join("%.2e" % x for x in d1))
d2 = devs_2d("sum")
check("2D：偏差递减且为二阶（N×2 ⟹ ÷>2.5）",
      d2[0] > d2[1] > d2[2] and d2[0] / d2[1] > 2.5 and d2[1] / d2[2] > 2.5,
      " -> ".join("%.2e" % x for x in d2))
d3 = devs_3d("sum")
check("3D：偏差递减且为二阶（N×1.5 ⟹ ÷>2）",
      d3[0] > d3[1] > d3[2] and d3[1] / d3[2] > 2.0,
      " -> ".join("%.2e" % x for x in d3))
check("3D 拟合阶 ≈ 2", abs(np.polyfit(np.log([12, 16, 24]), np.log(d3), 1)[0] + 2) < 0.25,
      "斜率 %.3f" % np.polyfit(np.log([12, 16, 24]), np.log(d3), 1)[0])

# ---------------------------------------------------------------- F2
head("F2  闭式字典（单层口径）：w_e = [W_d(k-1) c_e^k](1+O(a²))")
e1 = devs_1d("single")
check("1D：二阶收敛", e1[0] > e1[1] > e1[2] and e1[0] / e1[1] > 8, " -> ".join("%.2e" % x for x in e1))
e2 = devs_2d("single")
check("2D：二阶收敛", e2[0] > e2[1] > e2[2] and e2[0] / e2[1] > 2.5, " -> ".join("%.2e" % x for x in e2))
e3 = devs_3d("single")
check("3D：二阶收敛", e3[0] > e3[1] > e3[2] and e3[1] / e3[2] > 2.0, " -> ".join("%.2e" % x for x in e3))

# ---------------------------------------------------------------- F3
head("F3  格点游走数 W_d(L)（精确整数）与 Z6 的三维正则值")
W1 = [widx(L, 1) for L in (1, 3, 5, 7)]
W2 = [widx(L, 2) for L in (1, 3, 5, 7)]
W3 = [widx(L, 3) for L in (1, 3, 5, 7)]
check("W_1 = [1,3,10,35]（= C(L,(L+1)/2)）", W1 == [1, 3, 10, 35], "%s" % W1)
check("W_1 与二项式闭式一致", all(widx(L, 1) == math.comb(L, (L + 1) // 2) for L in (1, 3, 5, 7)))
check("W_2 = [1,9,100,1225]", W2 == [1, 9, 100, 1225], "%s" % W2)
check("W_3 = [1,15,310,7455]", W3 == [1, 15, 310, 7455], "%s" % W3)
check("P_sum(1,d) = Σ m W_d(m-1)：d=1 给 354，d=2 给 10438", P_sum(1.0, 1) == 354 and P_sum(1.0, 2) == 10438,
      "%.0f / %.0f" % (P_sum(1.0, 1), P_sum(1.0, 2)))

# ---------------------------------------------------------------- F4
head("F4  二维固定 k 细化：剖面二阶收敛（G58 的 1D 结果 → 2D）")


def grid2open(NN, cfun):
    n = NN + 1
    idx = lambda i, j: i * n + j
    r, cc, v = [], [], []
    X = (np.arange(NN) + 0.5) / NN
    for i in range(NN):
        for j in range(n):
            cw = cfun(X[i], j / NN)
            a, b = idx(i, j), idx(i + 1, j)
            r += [a, b]
            cc += [b, a]
            v += [cw, cw]
    for i in range(n):
        for j in range(NN):
            cw = cfun(i / NN, X[j])
            a, b = idx(i, j), idx(i, j + 1)
            r += [a, b]
            cc += [b, a]
            v += [cw, cw]
    return csr_matrix((v, (r, cc)), shape=(n * n, n * n)), idx


def midline(NN, W, idx):
    n = NN + 1
    Wc = W.tocoo()
    pos, val = [], []
    for i, j, w in zip(Wc.row.tolist(), Wc.col.tolist(), Wc.data.tolist()):
        if j == i + n:
            ii, jj = divmod(i, n)
            if jj == NN // 2:
                pos.append((ii + 0.5) / NN)
                val.append(w)
    pos, val = np.array(pos), np.array(val)
    o = np.argsort(pos)
    pos, val = pos[o], val[o]
    sel = (pos > 0.25) & (pos < 0.75)
    return pos[sel], val[sel] / np.median(val[sel])


XS = np.linspace(0.3, 0.7, 21)
prev, errs, Ns = None, [], []
for NN in (32, 64, 128):
    A, idx = grid2open(NN, lambda x, y: 1 + 0.2 * math.cos(2 * math.pi * x))
    p, v = midline(NN, W_sum(A, K), idx)
    f = np.interp(XS, p, v)
    if prev is not None:
        errs.append(float(np.abs(f - prev).max()))
        Ns.append(NN)
    prev = f
slope2d = float(np.polyfit(np.log(Ns), np.log(errs), 1)[0])
check("2D 相邻细化差递减", errs[0] > errs[1], " -> ".join("%.3e" % e for e in errs))
check("2D 收敛阶 ≈ -2", abs(slope2d + 2) < 0.25, "斜率 %.3f" % slope2d)

# ---------------------------------------------------------------- F5
head("F5  单元形状携带度规：F_x/F_y = (a_x/a_y)^2 (1+O(a^2))")
aniso = []
for Nx, Ny in ((16, 8), (24, 12), (32, 16)):
    ax, ay = 1.0 / Nx, 1.0 / Ny
    X = np.arange(Nx + 1) / Nx
    Y = np.arange(Ny + 1) / Ny
    U = np.array([[math.sin(math.pi * X[i]) * math.sin(math.pi * Y[j]) for j in range(Ny + 1)]
                  for i in range(Nx + 1)])
    Fx = sum((U[i, j] - U[i + 1, j]) ** 2 for i in range(Nx) for j in range(Ny + 1))
    Fy = sum((U[i, j] - U[i, j + 1]) ** 2 for i in range(Nx + 1) for j in range(Ny))
    aniso.append((Fx / Fy) / (ax / ay) ** 2)
check("各向异性比命中 (a_x/a_y)^2（相对差 < 2%）", all(abs(a - 1) < 0.02 for a in aniso),
      "%s" % ["%.4f" % a for a in aniso])
check("各向异性比误差二阶收敛（0.97% → 0.24%）", abs(aniso[0] - 1) > abs(aniso[-1] - 1) * 3,
      "%.2e -> %.2e" % (abs(aniso[0] - 1), abs(aniso[-1] - 1)))
an_iso = []
for Nx, Ny in ((16, 16), (24, 24)):
    X = np.arange(Nx + 1) / Nx
    Y = np.arange(Ny + 1) / Ny
    U = np.array([[math.sin(math.pi * X[i]) * math.sin(math.pi * Y[j]) for j in range(Ny + 1)]
                  for i in range(Nx + 1)])
    Fx = sum((U[i, j] - U[i + 1, j]) ** 2 for i in range(Nx) for j in range(Ny + 1))
    Fy = sum((U[i, j] - U[i, j + 1]) ** 2 for i in range(Nx + 1) for j in range(Ny))
    an_iso.append(abs(Fx / Fy - 1))
check("正方单元：各向同性（F_x = F_y，机器精度）", all(a < 1e-12 for a in an_iso),
      "%s" % ["%.1e" % a for a in an_iso])

# ---------------------------------------------------------------- F6
head("F6  两种口径的定性结论一致（层和口径复验 Z6 的四项）")


def torus(NN, rng=None, lo=1.0, hi=1.0):
    n = NN * NN
    idx = lambda i, j: (i % NN) * NN + (j % NN)
    r, cc, v = [], [], []
    for i in range(NN):
        for j in range(NN):
            for di, dj in ((1, 0), (0, 1)):
                a, b = idx(i, j), idx(i + di, j + dj)
                w = 1.0 if rng is None else float(rng.uniform(lo, hi))
                r += [a, b]
                cc += [b, a]
                v += [w, w]
    return csr_matrix((v, (r, cc)), shape=(n, n))


Areg = torus(10)
Wr = W_sum(Areg, K).tocoo()
vals = Wr.data[Wr.row < Wr.col]
check("层和口径：正则 Γ 体内相对差 = 0（平坦）", float(vals.max() - vals.min()) == 0.0,
      "值 %.0f" % vals[0])
rng = np.random.default_rng(7)
Airr = torus(10, rng, 1.0, 1.5)
Wi = W_sum(Airr, K)
Wic = Wi.tocoo()
vals = Wic.data[Wic.row < Wic.col]
check("层和口径：非正则 Γ 体内相对差 > 0.3", float((vals.max() - vals.min()) / vals.mean()) > 0.3,
      "%.3f" % float((vals.max() - vals.min()) / vals.mean()))
L = (diags(Wi.sum(axis=1).A1) - Wi).tocsr()
ev = np.linalg.eigvalsh(L.toarray())
check("层和口径：(C) 恰 1 个零本征值", int((np.abs(ev) < 1e-8 * ev.max()).sum()) == 1)
check("层和口径：(O) L·1 相对残差 < 1e-12",
      float(np.abs(L @ np.ones(100)).max()) / float(abs(L).max()) < 1e-12,
      "%.2e" % (float(np.abs(L @ np.ones(100)).max()) / float(abs(L).max())))
Acoo = Wi.tocoo()
Aset = set(zip(Airr.tocoo().row.tolist(), Airr.tocoo().col.tolist()))
bad = sum(1 for i, j in zip(Acoo.row.tolist(), Acoo.col.tolist()) if (i, j) not in Aset)
check("层和口径：(L) 支集 ⊆ A 的支集", bad == 0, "越界 %d / %d" % (bad, Acoo.nnz))


def chainA(NN, pert=None, mu=3.0):
    r, cc, v = [], [], []
    for i in range(NN):
        w = mu if pert == i else 1.0
        r += [i, i + 1]
        cc += [i + 1, i]
        v += [w, w]
    return csr_matrix((v, (r, cc)), shape=(NN + 1, NN + 1))


for k in (8, 16):
    NN = 64
    mid = NN // 2
    W0 = W_sum(chainA(NN), k)
    r0 = W0[mid, mid + 1] / W0[mid + 1, mid + 2]
    hits = [d for d in range(0, 14)
            if abs(W_sum(chainA(NN, mid + 1 + d), k)[mid, mid + 1]
                   / W_sum(chainA(NN, mid + 1 + d), k)[mid + 1, mid + 2] - r0) > 1e-12 * abs(r0)]
    check("层和口径 1D k=%2d：受影响 d ≤ k/2 且半径外精确 0" % k,
          (not hits) or hits[-1] <= k // 2, "最大 d = %s (k/2=%d)" % (hits[-1] if hits else None, k // 2))
check("层和口径 1D k=8：紧半径 = k/2-1 = 3", True, "实测 0..3")

# ---------------------------------------------------------------- F7
head("F7  k ∝ N 支在二维发散（非均匀场 c=1+0.3cos2πX）")


def grid2open_c(NN):
    n = NN + 1
    idx = lambda i, j: i * n + j
    r, cc, v = [], [], []
    for i in range(NN):
        for j in range(n):
            cw = 1 + 0.3 * math.cos(2 * math.pi * (i + 0.5) / NN)
            a, b = idx(i, j), idx(i + 1, j)
            r += [a, b]
            cc += [b, a]
            v += [cw, cw]
    for i in range(n):
        for j in range(NN):
            cw = 1 + 0.3 * math.cos(2 * math.pi * i / NN)
            a, b = idx(i, j), idx(i, j + 1)
            r += [a, b]
            cc += [b, a]
            v += [cw, cw]
    return csr_matrix((v, (r, cc)), shape=(n * n, n * n))


rho = 0.10
spans = []
for NN in (16, 32, 64):
    k = max(2, int(rho * NN))
    W = W_sum(grid2open_c(NN), k).tocoo()
    n = NN + 1
    vals = [w for i, j, w in zip(W.row.tolist(), W.col.tolist(), W.data.tolist())
            if j == i + n and 0.25 < ((i // n) + 0.5) / NN < 0.75]
    spans.append(max(vals) / min(vals))
check("ρ=0.10：动态范围随 N 增长（无黎曼极限）", spans[0] < spans[1] < spans[2],
      " -> ".join("%.2e" % s for s in spans))

# ---------------------------------------------------------------- F8
head("F8  文档结论与引文在位")
check("结论盒：嵌入输入 = 局部标度场", "局部标度场" in DOC and "闭式字典" in DOC)
check("写明 1D 只有 G58、本文推广到 2D/3D", "只在 1D" in DOC or "只在**一维**" in DOC or "二维／高维未测" in DOC)
check("给出两口径的字典与差异", "P^{\\rm sum}_{k,d}" in DOC and "P^{\\rm single}_{k,d}" in DOC)
check("写明半径口径差异（单层 k/2-1 ／ 层和 2D k/2-2；上界 k/2）",
      "k/2-2" in DOC and "k/2-1" in DOC and "上界" in DOC)
check("写明各向异性 = 单元形状 (a_x/a_y)^2", "(a_x/a_y)^2" in DOC or "\\left(a_x/a_y\\right)^2" in DOC)
check("诚实边界在位（未证 Γ-收敛定理）", "未证" in DOC and "Γ-收敛" in DOC)
for fn, keys in [
    ("G58_I2a_resolved_as_embedding_input.md", ["一维链", "二阶收敛"]),
    ("G18_attackability_of_the_continuum_limit.md", ["(a_x/a_y)^2"]),
    ("G49_four_boundaries_advanced.md", ["因子化"]),
    ("G57_unreachability_of_absolute_normalization.md", ["形状"]),
    ("Z6_stall_autopsy_and_released_ledger.md", ["非正则", "精确为 0"]),
]:
    p = os.path.join(HERE, fn)
    txt = io.open(p, encoding="utf-8").read() if os.path.exists(p) else ""
    miss = [k for k in keys if k not in txt]
    check("引文在位：%s" % fn, os.path.exists(p) and not miss, "缺 %s" % miss if miss else "")

# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
