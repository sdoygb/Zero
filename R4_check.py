#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R4_check.py -- 条件恢复之外的候选可检验预言

对应文档 R4_falsifiable_prediction.md。只读文件，不写回任何既有文件。

  F1  文档入口、状态边界与四候选字段
  F2  KPP 闭式、尾速率与当前误差
  F3  游走数、O(a^2) 序列与 L=4 二值标度
  F4  四维 Lovelock 动态项计数
  F5  量子谱、维数与 Casimir
  F6  2D/3D 面积律与离散修正重算
  F7  主预言、失败阈值、改动清单与未解项在位
"""

import io
import os
import sys
from itertools import product
from math import cosh, exp, log, tanh

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


DOC = io.open(
    os.path.join(HERE, "R4_falsifiable_prediction.md"), encoding="utf-8"
).read()


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


def anchor(*tokens):
    return all(tok in DOC for tok in tokens)


def relerr(a, b):
    return abs(a - b) / abs(b)


def cstar(B):
    """mu* tanh(mu*) - log cosh(mu*) = log(B)/2; return (c*, mu*)."""
    target = 0.5 * log(B)
    f = lambda mu: mu * tanh(mu) - log(cosh(mu))
    lo, hi = 1e-12, 200.0
    if f(hi) < target:
        return None, None
    for _ in range(220):
        mid = 0.5 * (lo + hi)
        if f(mid) < target:
            lo = mid
        else:
            hi = mid
    mu = 0.5 * (lo + hi)
    return tanh(mu), mu


def walk_count_to_neighbor(dim, steps):
    """Number of exact lattice walks from 0 to e_0 in Z^dim after `steps` steps."""
    target = (1,) + (0,) * (dim - 1)
    cur = {(0,) * dim: 1}
    for _ in range(steps):
        nxt = {}
        for pos, count in cur.items():
            for axis in range(dim):
                for delta in (-1, 1):
                    q = list(pos)
                    q[axis] += delta
                    q = tuple(q)
                    nxt[q] = nxt.get(q, 0) + count
        cur = nxt
    return cur.get(target, 0)


def canonical_rotation(word):
    return min(word[i:] + word[:i] for i in range(len(word)))


def build_2d(N, mass=0.0):
    idx = lambda i, j: (i % N) * N + (j % N)
    H = np.zeros((N * N, N * N))
    for i in range(N):
        for j in range(N):
            a = idx(i, j)
            for di, dj in ((1, 0), (0, 1)):
                b = idx(i + di, j + dj)
                H[a, b] = -1.0
                H[b, a] = -1.0
            if mass:
                H[a, a] += mass * ((-1) ** (i + j))
    return H


def build_3d(N, mass=0.0):
    idx = lambda i, j, k: ((i % N) * N + (j % N)) * N + (k % N)
    D = N ** 3
    H = np.zeros((D, D))
    for i in range(N):
        for j in range(N):
            for k in range(N):
                a = idx(i, j, k)
                for di, dj, dk in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
                    b = idx(i + di, j + dj, k + dk)
                    H[a, b] = -1.0
                    H[b, a] = -1.0
                if mass:
                    H[a, a] += mass * ((-1) ** (i + j + k))
    return H


def correlation_matrix(H):
    ev, U = np.linalg.eigh(H)
    occupied = ev < 0
    return U[:, occupied] @ U[:, occupied].conj().T, ev


def entropy_2d(C, L, N):
    sel = [(i * N + j) for i in range(L) for j in range(L)]
    CA = C[np.ix_(sel, sel)]
    nu = np.clip(np.linalg.eigvalsh((CA + CA.conj().T) / 2), 1e-13, 1 - 1e-13)
    return float(-np.sum(nu * np.log(nu) + (1 - nu) * np.log(1 - nu)))


def entropy_3d(C, L, N):
    sel = [
        ((i * N) + j) * N + k
        for i in range(L)
        for j in range(L)
        for k in range(L)
    ]
    CA = C[np.ix_(sel, sel)]
    nu = np.clip(np.linalg.eigvalsh((CA + CA.conj().T) / 2), 1e-13, 1 - 1e-13)
    return float(-np.sum(nu * np.log(nu) + (1 - nu) * np.log(1 - nu)))


def area_fit_3d(x, y):
    A = np.column_stack([x ** 2, x, np.ones_like(x)])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    rms = float(np.sqrt(np.mean((A @ coef - y) ** 2)))
    return float(coef[0]), float(coef[1]), float(coef[2]), rms


# ======================================================================
head("F1  文档入口、状态边界与四候选字段")

check("R4 文档存在并链接核验脚本", os.path.exists(os.path.join(HERE, "R4_falsifiable_prediction.md"))
      and "[`R4_check.py`](R4_check.py)" in DOC)
check("明确条件恢复，不冒充无条件 GR",
      anchor("当前状态边界", "条件恢复", "不从 `Z0` 无条件导出 GR", "E1/E2/E3/E4"))
check("四个候选均在位",
      anchor("C1 有限前沿锥与指数尾", "C2 离散与有限尺寸偏差",
             "C3 Lovelock 高曲率项尺度", "C4 量子谱与面积律"))
check("每类候选的必填字段均在位",
      all(anchor(x) for x in [
          "公式", "参数依赖", "当前误差", "实验或观测窗口", "可证伪条件",
          "函数形式", "数值输入",
      ]))
check("等级边界单列物理单位不可导出", anchor("物理单位", "E4 约定", "不可导出"))

# ======================================================================
head("F2  KPP 闭式、尾速率与当前误差")

expected = {
    2: (0.779944, 1.0452, 1.0167, 0.0273),
    3: (0.934697, 1.6943, 1.5509, 0.0846),
}
roots = {}
for B, (c_exp, mu_exp, slope_abs, err_exp) in expected.items():
    c, mu = cstar(B)
    roots[B] = (c, mu)
    err = abs(slope_abs - mu) / mu
    print("      B=%d  c*=%.6f  mu*=%.4f  实测斜率=-%.4f  相对差=%.2f%%"
          % (B, c, mu, slope_abs, 100 * err))
    check("B=%d：c* 与语料一致" % B, relerr(c, c_exp) < 2e-4)
    check("B=%d：mu* 与语料一致" % B, relerr(mu, mu_exp) < 2e-3)
    check("B=%d：尾速率误差与 R4 表一致" % B, abs(err - err_exp) < 8e-4)

c4, mu4 = cstar(4)
print("      B=4  c*=%.9f  mu*=%.4f  每格衰减=%.3e"
      % (c4, mu4, exp(-mu4)))
check("B=4：c*=1", abs(c4 - 1.0) < 1e-9)
check("B=4：mu*>30", mu4 > 30.0, "mu*=%.4f" % mu4)
check("B=4：每格衰减在 1e-15 量级", 5e-15 < exp(-mu4) < 8e-15)
check("文档登记 B=3 浅窗 13.6% 与 B 未定", anchor("13.6%", "B∈{2,3,4}"))

# ======================================================================
head("F3  游走数、O(a^2) 序列与 L=4 二值标度")

w1 = [walk_count_to_neighbor(1, L) for L in (1, 3, 5, 7)]
w2 = [walk_count_to_neighbor(2, L) for L in (1, 3, 5, 7)]
w3 = [walk_count_to_neighbor(3, L) for L in (1, 3, 5, 7)]
print("      W1=%s  W2=%s  W3=%s" % (w1, w2, w3))
check("W1 精确游走数 = [1,3,10,35]", w1 == [1, 3, 10, 35])
check("W2 精确游走数 = [1,9,100,1225]", w2 == [1, 9, 100, 1225])
check("W3 精确游走数 = [1,15,310,7455]", w3 == [1, 15, 310, 7455])

err_1d = [1.47e-2, 9.17e-4, 5.73e-5]
err_2d = [3.94e-2, 9.75e-3, 2.43e-3]
err_3d = [7.41e-2, 4.17e-2, 1.85e-2]
r1 = [err_1d[i] / err_1d[i + 1] for i in range(len(err_1d) - 1)]
r2 = [err_2d[i] / err_2d[i + 1] for i in range(len(err_2d) - 1)]
r3 = [err_3d[i] / err_3d[i + 1] for i in range(len(err_3d) - 1)]
print("      1D 比值=%s（预期 16）  2D 比值=%s（预期 4）  3D 比值=%s（预期 1.78,2.25）"
      % (["%.2f" % x for x in r1], ["%.2f" % x for x in r2], ["%.2f" % x for x in r3]))
check("1D 字典偏差按 O(a^2) 缩小", all(abs(x / 16.0 - 1.0) < 0.02 for x in r1))
check("2D 字典偏差按 O(a^2) 缩小", all(abs(x / 4.0 - 1.0) < 0.03 for x in r2))
check("3D 字典偏差比例命中 (16/12)^2 与 (24/16)^2",
      abs(r3[0] / (16.0 / 12.0) ** 2 - 1.0) < 0.02
      and abs(r3[1] / (24.0 / 16.0) ** 2 - 1.0) < 0.02)

classes = {}
for T_age in (2, 4):
    for w in product((-1, 1), repeat=T_age):
        if sum(w) == 0:
            classes.setdefault(canonical_rotation(w), []).append(w)
sizes = sorted(len(v) for v in classes.values())
mean_size = float(np.mean(sizes))
norm = sorted(set(round(size / mean_size, 10) for size in sizes))
print("      T_age=2,4 合并旋转类规模=%s  归一化标度=%s" % (sizes, norm))
check("T_age=2,4 合并旋转类规模 = {2,2,4}", sizes == [2, 2, 4])
check("归一化标度只有 {0.75,1.5}", norm == [0.75, 1.5])
check("二值标度对比度为 2，度规对比度上界为 16",
      abs(max(norm) / min(norm) - 2.0) < 1e-12 and 2 ** 4 == 16)

# ======================================================================
head("F4  四维 Lovelock 动态项计数")

D4 = 4
m_max_4 = (D4 - 1) // 2
D5 = 5
m_max_5 = (D5 - 1) // 2
print("      D=4: 动态 m=0,...,%d；D=5: 动态 m=0,...,%d"
      % (m_max_4, m_max_5))
check("D=4 动态 Lovelock 项只有 m=0,1", m_max_4 == 1)
check("D=5 才允许独立高曲率 m=2", m_max_5 == 2)
check("文档明确 Gauss-Bonnet 在 D=4 是拓扑项且不进入场方程",
      anchor("Gauss-Bonnet", "拓扑项", "不进入四维场方程"))
check("文档区分 D=4 零假设、D>4 的 alpha2 输入与无法给出的 SI 界限",
      anchor("零假设", "D>4", "α₂", "不能给出"))

# ======================================================================
head("F5  量子谱、维数与 Casimir")

sx = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
sy = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
sz = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
J = [0.5 * sx, 0.5 * sy, 0.5 * sz]
J2 = sum(mat @ mat for mat in J)
j3_spectrum = np.linalg.eigvalsh(J[2])
T_age = 4
dim_A = 4 * (T_age + 1)
dim_H = 2 * (T_age + 1)
print("      T_age=4: dim A=%d  dim H=%d  J3=%s  J2=%s"
      % (dim_A, dim_H, np.array2string(j3_spectrum, precision=3), np.array2string(J2.diagonal(offset=0), precision=3)))
check("SU(2) 对易关系成立", np.allclose(J[0] @ J[1] - J[1] @ J[0], 1j * J[2]))
check("Casimir J^2 = 3/4 I", np.allclose(J2, 0.75 * np.eye(2)))
check("J3 谱 = {-1/2,+1/2}", np.allclose(sorted(j3_spectrum), [-0.5, 0.5]))
check("T_age=L=4 给出 dim A=20, dim H=10", dim_A == 20 and dim_H == 10)
check("文档登记 KMS 残差、加性域与测量诠释输入",
      anchor("3.6e-16", "加性域", "测量诠释"))

# ======================================================================
head("F6  2D/3D 面积律与离散修正重算")

N2 = 24
LS2 = np.arange(2, 10, dtype=float)
b2 = {}
for m in (0.0, 0.25, 0.5, 1.0):
    C, _ = correlation_matrix(build_2d(N2, mass=m))
    ss = np.array([entropy_2d(C, int(L), N2) for L in LS2])
    slopes = np.diff(ss)
    xm = LS2[:-1]
    X = np.column_stack([np.ones_like(xm), np.log(xm)])
    coef, *_ = np.linalg.lstsq(X, slopes, rcond=None)
    b2[m] = float(coef[1])
    print("      2D m=%.2f  b=%.3f" % (m, b2[m]))
check("2D 无 gap 的 b≈0.549", abs(b2[0.0] - 0.549) < 2e-3)
check("2D m=0.25 的 b≈0.150", abs(b2[0.25] - 0.150) < 2e-3)
check("2D m=0.5 的 b≈0.051", abs(b2[0.5] - 0.051) < 2e-3)
check("2D m=1 的 b≈0.007", abs(b2[1.0] - 0.007) < 2e-3)
check("2D b(m) 单调递减且跨度 > 50 倍",
      b2[0.0] > b2[0.25] > b2[0.5] > b2[1.0] and b2[0.0] / b2[1.0] > 50.0)

N3 = 12
LS3 = np.arange(2, 7, dtype=float)
a3 = {}
rms3 = {}
for m in (2.0, 3.0, 4.0):
    C, _ = correlation_matrix(build_3d(N3, mass=m))
    ss = np.array([entropy_3d(C, int(L), N3) for L in LS3])
    a, b, c, rms = area_fit_3d(LS3, ss)
    a3[m] = a
    rms3[m] = rms
    print("      3D m=%.1f  a=%.4f  b=%.4f  c=%.4f  rms=%.2e"
          % (m, a, b, c, rms))
    check("3D m=%.1f 面系数与语料一致" % m, abs(a - {2.0: 0.6929, 3.0: 0.4770, 4.0: 0.3475}[m]) < 2e-4)
    check("3D m=%.1f rms<1e-4" % m, rms < 1e-4)

prods = np.array([(a3[m] / 3.0) * m for m in (2.0, 3.0, 4.0)])
print("      3D a_face*m = %.4f +- %.4f（相对 %.2f%%）"
      % (prods.mean(), prods.std(), 100 * prods.std() / prods.mean()))
check("3D a_face*m = 0.4674±0.0068", abs(prods.mean() - 0.4674) < 2e-4
      and abs(prods.std() - 0.0068) < 1e-3)
check("强 gap rms 序列随 m 递减",
      rms3[2.0] > rms3[3.0] > rms3[4.0])

L_a5 = np.array([8.0, 16.0, 32.0, 64.0, 128.0, 256.0])
gap_a5 = np.array([0.23266113, 0.11601535, 0.05796885, 0.02897959, 0.01448919, 0.00724452])
xi_over_L = (2.0 / gap_a5) / L_a5
print("      Z3 点汇 xi/L = %s" % np.array2string(xi_over_L, precision=4))
check("Z3 点汇给出 xi≈1.078L（相对涨落 <0.5%）",
      abs(xi_over_L.mean() - 1.078) < 2e-3 and xi_over_L.std() / xi_over_L.mean() < 0.005)
check("文档登记完整面积律前提永不满足",
      anchor("ξ≈1.078L", "前提", "不满足"))

# ======================================================================
head("F7  主预言、失败阈值、改动清单与未解项在位")

check("主预言 R4-P1 在位",
      anchor("主预言 R4-P1", "S_2(L,m)", "b_2(0.25)", "a_{\\rm face}(m)m"))
check("主预言写明 2D/3D 的可执行拟合方案",
      anchor("扩展二维拟合", "扩展三维拟合", "AIC/BIC", "自助法"))
check("成功与失败阈值在位",
      anchor("成功与失败阈值", "95% 区间排除", "rms<1e-4", "出现第三标度"))
check("实验窗口明确当前只能先数值/量子模拟",
      anchor("量子模拟", "不能给出米制", "E1", "E4"))
check("改动文件清单仅列两个新文件",
      anchor("R4_falsifiable_prediction.md", "R4_check.py", "没有修改任何既有文件"))
check("未解项列全至少 E1、Γ-收敛、E3、E2、B/L、精确锥",
      anchor("E1 类到站点", "Γ-收敛", "E3 源作用量", "E2 `D=4`", "`B` 和 `L`", "精确锥"))

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
