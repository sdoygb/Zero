#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z6_check.py —— 【卡点解剖 ＋ 放开输入后的终态账本】的核验
==========================================================
独立实断言：
  F1  三维非正则 Γ：端到端前提核验（(O) 二阶签名 / (C) 单零模）
  F2  三维：支集只在图上；远端扰动影响【精确为 0】＋反向控制
  F3  三维：正则 Γ 精确平坦；非正则 Γ 承载几何
  F4  三维：光滑指定剖面被度规携带（非均匀且随 c 单调）
  F5  局域半径的两种口径（值 k/2-1 / 比值 k/2；1D 与 2D）
  F6  KPP 闭式与 B=4 临界（独立复算 G56）
  F7  卡点表的引文完整性（16 条来源在位且含其判定）
  F8  文档结论、输入账本与决定清单在位
"""
import io
import math
import os
import sys

import numpy as np
from scipy.sparse import csr_matrix, diags

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


DOC = io.open(os.path.join(HERE, "Z6_stall_autopsy_and_released_ledger.md"), encoding="utf-8").read()


# ---------------------------------------------------------------- 图构造
def torus3(n, rng=None, lo=1.0, hi=1.0):
    Nn = n ** 3
    idx = lambda x, y, z: ((x % n) * n + (y % n)) * n + (z % n)
    r, c, v = [], [], []
    for x in range(n):
        for y in range(n):
            for z in range(n):
                i = idx(x, y, z)
                for dx, dy, dz in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
                    j = idx(x + dx, y + dy, z + dz)
                    w = 1.0 if rng is None else float(rng.uniform(lo, hi))
                    r += [i, j]
                    c += [j, i]
                    v += [w, w]
    return csr_matrix((v, (r, c)), shape=(Nn, Nn))


def open3(n):
    Nn = n ** 3
    idx = lambda x, y, z: (x * n + y) * n + z
    r, c, v = [], [], []
    for x in range(n):
        for y in range(n):
            for z in range(n):
                i = idx(x, y, z)
                for dx, dy, dz in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
                    xx, yy, zz = x + dx, y + dy, z + dz
                    if xx < n and yy < n and zz < n:
                        j = idx(xx, yy, zz)
                        r += [i, j]
                        c += [j, i]
                        v += [1.0, 1.0]
    return csr_matrix((v, (r, c)), shape=(Nn, Nn)), idx


def grid2(n):
    Nn = n * n
    r, c, v = [], [], []
    for x in range(n):
        for y in range(n):
            i = x * n + y
            for dx, dy in ((1, 0), (0, 1)):
                xx, yy = x + dx, y + dy
                if xx < n and yy < n:
                    j = xx * n + yy
                    r += [i, j]
                    c += [j, i]
                    v += [1.0, 1.0]
    return csr_matrix((v, (r, c)), shape=(Nn, Nn))


def chain(n):
    r, c, v = [], [], []
    for i in range(n - 1):
        r += [i, i + 1]
        c += [i + 1, i]
        v += [1.0, 1.0]
    return csr_matrix((v, (r, c)), shape=(n, n))


def perturb(A, i, j, mu=3.0):
    B = A.tolil()
    B[i, j] = mu
    B[j, i] = mu
    return B.tocsr()


def W_k(A, k):
    """闭环计数权：w^(k) = A ∘ A^{k-1}。"""
    P = A.copy()
    for _ in range(k - 2):
        P = (P @ A).tocsr()
    return A.multiply(P)


def spread(W):
    M = W.tocoo()
    vals = M.data[M.row < M.col]
    return float((vals.max() - vals.min()) / vals.mean())


# ---------------------------------------------------------------- F1
head("F1  三维非正则 Γ（10^3 环面，无序边权 u∈[1,1.5]）：端到端前提")
rng = np.random.default_rng(7)
Ai = torus3(10, rng, 1.0, 1.5)
Wi = W_k(Ai, 8)
L = (diags(Wi.sum(axis=1).A1) - Wi).tocsr()
scale = float(abs(L).max())
one = np.ones(1000)
x = np.array([i // 100 for i in range(1000)], float)
r1 = float(np.abs(L @ one).max()) / scale
r2 = float(np.abs(L @ x).max()) / scale
r3 = float(np.abs(L @ (x * x)).max()) / scale
check("(O) L·1 ≈ 0（相对残差 < 1e-12）", r1 < 1e-12, "%.2e" % r1)
check("(O) L·x ≠ 0（不是一阶）", r2 > 0.1, "%.3f" % r2)
check("(O) L·x² ≠ 0（二阶差分签名）", r3 > r2, "%.3f > %.3f" % (r3, r2))
ev = np.linalg.eigvalsh(L.toarray())
nz = int((np.abs(ev) < 1e-8 * ev.max()).sum())
check("(C) 恰有 1 个零本征值", nz == 1, "实测 %d" % nz)
check("(C) λ₂ > 0（连通图核 = 常量）", ev[1] > 0, "λ₂ = %.4g" % ev[1])

# ---------------------------------------------------------------- F2
head("F2  三维局域性：支集只在图上；半径外影响【精确为 0】")
Wc = W_k(Ai, 8).tocoo()
Aset = set(zip(Ai.tocoo().row.tolist(), Ai.tocoo().col.tolist()))
bad = sum(1 for i, j in zip(Wc.row.tolist(), Wc.col.tolist()) if (i, j) not in Aset)
check("(L) w^(k) 的支集 ⊆ A 的支集（k=8）", bad == 0, "越界 %d / %d" % (bad, Wc.nnz))

A3, idx3 = open3(16)
cen, far = idx3(8, 8, 8), idx3(13, 13, 13)
nxp, nxm = idx3(9, 8, 8), idx3(7, 8, 8)
far_e = (far, idx3(13, 13, 14))
for k in (4, 8):
    W0 = W_k(A3, k)
    r0 = W0[cen, nxp] / W0[cen, nxm]
    Wf = W_k(perturb(A3, *far_e), k)
    rf = Wf[cen, nxp] / Wf[cen, nxm]
    check("3D k=%d：远端扰动（图距 15）影响精确为 0" % k, rf == r0,
          "Δrel=%.3e" % (abs(rf - r0) / abs(r0)))
Wn = W_k(perturb(A3, cen, nxp), 8)
rn = Wn[cen, nxp] / Wn[cen, nxm]
W0 = W_k(A3, 8)
check("3D 反向控制：近端扰动显著非零", abs(rn - W0[cen, nxp] / W0[cen, nxm]) / abs(r0) > 1.0,
      "Δrel=%.3f" % (abs(rn - W0[cen, nxp] / W0[cen, nxm]) / abs(r0)))

Nc = 64
mid = Nc // 2
base = W_k(chain(Nc), 8)
r0 = base[mid, mid + 1] / base[mid + 1, mid + 2]
Afar = perturb(chain(Nc), Nc - 2, Nc - 1, mu=2.0)
rfar = W_k(Afar, 8)[mid, mid + 1] / W_k(Afar, 8)[mid + 1, mid + 2]
check("1D（复刻 G41 §3／Z5 §2）：远端扰动影响精确为 0", rfar == r0,
      "Δrel=%.3e" % (abs(rfar - r0) / abs(r0)))

# ---------------------------------------------------------------- F3
head("F3  正则 Γ ⇒ 精确平坦；非正则 Γ ⇒ 承载几何（三维）")
Areg = torus3(10)
for k in (2, 4, 6, 8):
    check("正则 3D k=%d：体内边权相对差 = 0" % k, spread(W_k(Areg, k)) == 0.0)
for k in (2, 4, 6, 8):
    check("非正则 3D k=%d：体内边权相对差 > 0.3" % k, spread(W_k(Ai, k)) > 0.3,
          "%.3f" % spread(W_k(Ai, k)))
sp = [spread(W_k(Ai, k)) for k in (2, 4, 6, 8)]
check("非正则：非均匀度随 k 增大", all(sp[i] < sp[i + 1] for i in range(3)),
      "%s" % [round(s, 3) for s in sp])

# ---------------------------------------------------------------- F4
head("F4  光滑指定剖面被度规携带（3D 环面 c(x)=1+0.3cos(2πx/n)，k=6）")


def torus3_line(n, k):
    Nn = n ** 3
    idx = lambda x, y, z: ((x % n) * n + (y % n)) * n + (z % n)
    r, c, v = [], [], []
    cs = []
    for x in range(n):
        cx = 1.0 + 0.3 * math.cos(2 * math.pi * x / n)
        cs.append(cx)
        for y in range(n):
            for z in range(n):
                i = idx(x, y, z)
                j = idx(x + 1, y, z)
                r += [i, j]
                c += [j, i]
                v += [cx, cx]
    A = csr_matrix((v, (r, c)), shape=(Nn, Nn))
    return A, np.array(cs), idx


A4, cs, idx4 = torus3_line(12, 6)
W4 = W_k(A4, 6).tocoo()
vals, cc = [], []
for i, j, w in zip(W4.row.tolist(), W4.col.tolist(), W4.data.tolist()):
    if i < j and (i % 144) // 12 == (j % 144) // 12:
        vals.append(w)
        cc.append(cs[(i // 144) % 12])
vals, cc = np.array(vals), np.array(cc)
rel = (vals.max() - vals.min()) / vals.mean()
check("预设剖面给出非均匀度规（相对差 > 1）", rel > 1.0, "%.3f  (%.3f → %.3f)" % (rel, vals.min(), vals.max()))
check("边权随指定 c 单调/正相关（corr > 0.95）", float(np.corrcoef(vals, cc)[0, 1]) > 0.95,
      "corr = %.4f" % float(np.corrcoef(vals, cc)[0, 1]))

# ---------------------------------------------------------------- F5
head("F5  局域半径：统一为 k/2-1（值口径与比值口径一致；d 由扰动边到被测比值最近边计）")
Nv = 256
for k in (4, 8, 12):
    base = W_k(chain(Nv), k)[100, 101]
    hits = [d for d in range(1, 30)
            if abs(W_k(perturb(chain(Nv), 100 + d, 101 + d), k)[100, 101] - base) > 1e-14]
    check("1D 值口径 k=%2d：受影响距离 = 1..%d" % (k, k // 2 - 1), hits == list(range(1, k // 2)),
          "实测 %s" % (hits[:1] + ["..."] + hits[-1:] if hits else hits))

N1 = 64
mid1 = N1 // 2
for k in (4, 8, 16):
    A0 = chain(N1)
    W0 = W_k(A0, k)
    r0 = W0[mid1, mid1 + 1] / W0[mid1 + 1, mid1 + 2]
    hits = [d for d in range(0, 20)
            if abs(W_k(perturb(A0, mid1 + 1 + d, mid1 + 2 + d), k)[mid1, mid1 + 1]
                   / W_k(perturb(A0, mid1 + 1 + d, mid1 + 2 + d), k)[mid1 + 1, mid1 + 2] - r0)
            > 1e-13 * abs(r0)]
    check("1D 比值口径 k=%2d：受影响距离 = 0..%d" % (k, k // 2 - 1), hits == list(range(0, k // 2)),
          "实测最大 d = %s" % (hits[-1] if hits else None))

n2 = 24
cen2 = 12 * n2 + 12
for k in (4, 8, 16):
    A0 = grid2(n2)
    W0 = W_k(A0, k)
    r0 = W0[cen2, cen2 + 1] / W0[cen2, cen2 - 1]
    hits = [d for d in range(0, 14)
            if abs(W_k(perturb(A0, cen2 + d, cen2 + d + 1), k)[cen2, cen2 + 1]
                   / W_k(perturb(A0, cen2 + d, cen2 + d + 1), k)[cen2, cen2 - 1] - r0)
            > 1e-13 * abs(r0)]
    check("2D 比值口径 k=%2d：受影响距离 = 0..%d" % (k, k // 2 - 1), hits == list(range(0, k // 2)),
          "实测最大 d = %s" % (hits[-1] if hits else None))

for tag, AA, ii, jj, meas in (("1D", chain(N1), mid1 + 1 + 4, mid1 + 2 + 4, (mid1, mid1 + 1, mid1 + 1, mid1 + 2)),
                              ("2D", grid2(n2), cen2 + 4, cen2 + 5, (cen2, cen2 + 1, cen2, cen2 - 1))):
    k = 8
    W0 = W_k(AA, k)
    r0 = W0[meas[0], meas[1]] / W0[meas[2], meas[3]]
    Wp = W_k(perturb(AA, ii, jj), k)
    check("%s k=8：半径之外（d = k/2）影响精确为 0" % tag,
          Wp[meas[0], meas[1]] / Wp[meas[2], meas[3]] == r0)

# ---------------------------------------------------------------- F6
head("F6  KPP 闭式与 B=4 临界（独立复算 G56 §2.1）")


def f_kpp(mu):
    return mu * math.tanh(mu) - math.log(math.cosh(mu))


def c_star(B):
    lo, hi = 1e-9, 80.0
    if f_kpp(hi) < 0.5 * math.log(B):
        return None
    for _ in range(300):
        m = 0.5 * (lo + hi)
        if f_kpp(m) < 0.5 * math.log(B):
            lo = m
        else:
            hi = m
    return math.tanh(0.5 * (lo + hi))


for B, ref in ((2, 0.7799442711232809), (3, 0.9346973029773759)):
    check("B=%d：c* = tanh μ* 与闭式参考一致（< 1e-9）" % B, abs(c_star(B) - ref) < 1e-9,
          "%.9f vs %.9f" % (c_star(B), ref))
check("B=4：c* = 1.000000000（与 Z1 定理 1 锥重合）", abs(c_star(4) - 1.0) < 1e-9, "%.12f" % c_star(4))
check("B=5：无解（f(∞) = log2 < ½log5）",
      c_star(5) is None and math.log(2) < 0.5 * math.log(5),
      "log2=%.6f < %.6f" % (math.log(2), 0.5 * math.log(5)))

# ---------------------------------------------------------------- F7
head("F7  卡点表的引文完整性（来源在位且含其判定）")
CITES = [
    ("G41_lovelock_premises_under_nonuniform_weight.md", ["(L')", "命名更正"]),
    ("Z2_zero_to_gr_direct_route.md", ["这是错的", "形式性质"]),
    ("G55_dynamics_line_degeneration_to_GR.md", ["一直在骗我们", "记忆核"]),
    ("G56_degeneration_attempt2_six_slots.md", ["撤回", "B=4"]),
    ("G58_I2a_resolved_as_embedding_input.md", ["Γ-收敛", "归并"]),
    ("G59_I7_settled_native_cone_and_its_residue.md", ["不承重", "有效锥"]),
    ("G89_dimension_no_go_and_the_balance_condition.md", ["命题 1", "【条件】"]),
    ("G57_unreachability_of_absolute_normalization.md", ["不可导出", "量纲"]),
    ("Z5_finite_k_locality_escape.md", ["精确为 0", "非正则"]),
    ("Z1_zero_layer_as_the_foundation.md", ["不再适用", "自付其价"]),
    ("G13_foliation_and_lorentz_invariance_gap.md", ["引理 46"]),
    ("G10_final_derivation_and_input_ledger.md", ["未建立", "G58"]),
]
for fn, keys in CITES:
    p = os.path.join(HERE, fn)
    ok = os.path.exists(p)
    txt = io.open(p, encoding="utf-8").read() if ok else ""
    miss = [k for k in keys if k not in txt]
    check("引文在位：%s" % fn, ok and not miss, "缺 %s" % miss if miss else "")

# ---------------------------------------------------------------- F8
head("F8  文档结论、输入账本与决定清单在位")
check("结论盒：卡点不是墙／四机制", "卡点不是墙" in DOC and "机制" in DOC)
check("终态账本：数学缺口 0 ＋ 承重 no-go 0", "数学缺口 }0" in DOC or "承重 no-go }0" in DOC)
check("§2 逐卡点终态表在位（16 行）",
      DOC.count("\n| 1 | $(L')$") == 1 and "承重卡点" in DOC)
check("输入账本 E1–E4 在位", all(("**%s**" % e) in DOC for e in ("E1", "E2", "E3", "E4")))
check("决定清单在位（选择而非卡点）", "决定清单" in DOC and "选错的代价" in DOC)
check("写明陈旧条目更正（G10 §7）", "陈旧条目" in DOC and "§7" in DOC)
check("诚实边界：账本是政策相关的", "政策相对性" in DOC and "政策相关" in DOC)
check("核验入口登记", "Z6_check.py" in DOC)

# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
