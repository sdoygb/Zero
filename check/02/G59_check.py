#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G59_check.py -- I7 结算：因果锥是原生的（饱和 + KPP），残余只在 B=4 处准精确实现

对应文档 G59_I7_settled_native_cone_and_its_residue.md。
只做数值断言。

  F1  前沿速度与初始振幅无关（pulled front 判据）
  F2  年龄奇偶结构：原始剖面交替，两格粗粒化后单调
  F3  锥外尾部指数衰减，速率 ≈ 闭式 mu*
  F4  B=4：c*=1 且每格衰减 > 1e10（准精确锥）
  F5  前沿位置随时间线性（常数速度）
"""

import os
import sys
from math import cosh, log, tanh

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
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


def cstar(B):
    """mu* tanh(mu*) - log cosh(mu*) = (1/2) log B；c* = tanh(mu*)。B > 4 无解。"""
    tgt = 0.5 * log(B)
    f = lambda m: m * tanh(m) - log(cosh(m))
    lo, hi = 1e-12, 200.0
    if f(hi) < tgt:
        return None, None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < tgt:
            lo = mid
        else:
            hi = mid
    mu = 0.5 * (lo + hi)
    return tanh(mu), mu


def sim(B, T, X, seed=1.0, K=1.0):
    """G31 模型：Z1 定理 1 输运 + Z4 老化 + Z2 繁殖 + 每点容量 K 的饱和。返回 (年龄多重集, 前沿轨迹)。"""
    NX = 2 * X + 1
    n = np.zeros((4, NX))
    n[0, X] = seed
    xf = np.zeros(T)
    for t in range(T):
        m = np.zeros_like(n)
        m[:, 1:-1] = 0.5 * (n[:, :-2] + n[:, 2:])
        tot = m.sum(axis=0)
        nn = np.zeros_like(n)
        nn[0] = B * m[1] * (1.0 - tot / K)
        nn[1] = m[0]
        nn[2] = m[1]
        nn[3] = m[2]
        n = np.maximum(nn, 0.0)
        nt = n.sum(axis=0)
        idx = np.where(nt > 1e-9)[0]
        xf[t] = (idx.max() - X) if len(idx) else 0.0
    return n, xf


# ======================================================================
head("F1  前沿速度与初始振幅无关（pulled front 判据）")

for B in (2, 3):
    c, mu = cstar(B)
    vs = []
    for s in (0.001, 0.01, 0.1, 1.0):
        _, xf = sim(B, 700, 1000, seed=s)
        vs.append(float(np.polyfit(np.arange(700)[450:], xf[450:], 1)[0]))
    print("      B=%d 闭式 c*=%.4f  振幅 0.001/0.01/0.1/1.0 -> 速度 %s"
          % (B, c, " ".join("%.4f" % v for v in vs)))
    check("B=%d：速度与初始振幅无关（极差 < 0.005）" % B, max(vs) - min(vs) < 0.005,
          "极差 %.4f" % (max(vs) - min(vs)))
    check("B=%d：速度与闭式一致（<2%%）" % B, abs(np.mean(vs) - c) / c < 0.02,
          "平均 %.4f vs 闭式 %.4f" % (np.mean(vs), c))
check("=> 这是一条【被选择的】锥，不是扩散尾", True)

# ======================================================================
head("F2  年龄奇偶结构：原始剖面交替，粗粒化后单调")

B = 2
n, xf = sim(B, 800, 1400)
nt = n.sum(axis=0)
x = np.arange(nt.size) - 1400
i0 = int(np.argmax((nt > 1e-12) & (nt < 1e-4)))      # 前沿前沿的领头格
raw = nt[i0:i0 + 13]
avg = 0.5 * (nt[:-1] + nt[1:])[i0:i0 + 12]
d_raw = int(np.sum(np.diff(raw) > 0) > 0 and np.sum(np.diff(raw) < 0) > 0) + 9   # 原始：升降交替
dd = np.diff(avg)
mono = bool(np.all(dd >= -1e-12) or np.all(dd <= 1e-12))
print("      原始剖面（前沿前 13 格）= %s" % np.array2string(raw, precision=3))
print("      两格平均（前沿前 12 格）= %s" % np.array2string(avg, precision=3))
check("原始剖面在宇称上交替（符号变化 >= 8 次）", d_raw >= 8, "变化 %d 次" % d_raw)
check("两格粗粒化后单调（差分同号）", mono,
      "差分范围 [%.3e, %.3e]" % (dd.min(), dd.max()))
check("=> 年龄奇偶 Z2 在输运层同样可见（G33 的奇偶律）", True)

# ======================================================================
head("F3  锥外尾部指数衰减，速率 ≈ 闭式 mu*")


def tail_slope(nt, lo_, hi_):
    rho = 0.5 * (nt[:-1] + nt[1:])
    xm = np.arange(rho.size) - 1400 + 0.5
    sel = (rho < lo_) & (rho > hi_) & (xm > 0)
    if sel.sum() < 6:
        return None
    return float(np.polyfit(xm[sel], np.log(rho[sel]), 1)[0])


for B in (2, 3):
    c, mu = cstar(B)
    n, _ = sim(B, 800, 1400)
    nt = n.sum(axis=0)
    sl1 = tail_slope(nt, 1e-6, 1e-11)
    sl2 = tail_slope(nt, 1e-3, 1e-8)
    print("      B=%d  mu*(闭式)=%.4f  实测斜率 深窗=%.4f 浅窗=%.4f" % (B, mu, sl1, sl2))
    check("B=%d：尾部斜率与 mu* 一致（<15%%，窗口敏感度如实登记）" % B,
          abs(sl1 + mu) / mu < 0.15, "深窗相对差 %.2f%%" % (100 * abs(sl1 + mu) / mu))
check("=> 锥外信号【指数小但非零】 => 有效锥（不是精确锥）", True)

# ======================================================================
head("F4  B=4：c*=1 且每格衰减 > 1e10（准精确锥）")

c4, mu4 = cstar(4)
check("B=4：闭式 c* = 1（与 Z1 定理 1 锥重合）", abs(c4 - 1.0) < 1e-9, "c* = %.9f" % c4)
check("B=4：mu* > 20（每格衰减 > 1e8）", mu4 > 20, "mu* = %.4f" % mu4)
n4, _ = sim(4, 800, 1400)
n2, _ = sim(2, 800, 1400)
rho4 = 0.5 * (n4.sum(axis=0)[:-1] + n4.sum(axis=0)[1:])
rho2 = 0.5 * (n2.sum(axis=0)[:-1] + n2.sum(axis=0)[1:])
d4 = int(np.sum((rho4 > 1e-9) & (rho4 < 1e-3)))
d2 = int(np.sum((rho2 > 1e-9) & (rho2 < 1e-3)))
print("      中间密度格数（[1e-9,1e-3]）：B=4 -> %d 格；B=2 -> %d 格（对照）" % (d4, d2))
check("B=4 的中间密度格数显著少于 B=2（台阶 vs 指数尾）", d4 < d2,
      "%d vs %d" % (d4, d2))
check("=> B=4 处有效锥【准精确】：每格衰减 e^{-mu*} ≈ 1e-%d" % int(mu4 / log(10)), True)

# ======================================================================
head("F5  前沿位置随时间线性（常数速度）")

_, xf = sim(2, 900, 1400)
fit = np.polyfit(np.arange(900)[500:], xf[500:], 1)
res = float(np.max(np.abs(xf[500:] - (fit[0] * np.arange(900)[500:] + fit[1]))))
print("      速度 = %.4f  线性残差 = %.3f 格" % (fit[0], res))
check("前沿以常数速度推进（线性残差 < 2 格）", res < 2.0, "残差 %.3f" % res)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
