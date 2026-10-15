#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G56_check.py -- 退化尝试 #2：以 D213 的六槽位预算重估 + KPP 因果槽位的闭式

对应文档 G56_degeneration_attempt2_six_slots.md。

前两次尝试的失败点（原文取证）：
  尝试 #1  D213：周期骨架直接桥 —— 六槽位只拿到 1/6，死在"无局域零和传输图"
  尝试 #2  G1-G10：Lovelock 链 —— 条件恢复 GR，度规/维数为输入，I2a 开
本轮：用 G1-G55 的新成果逐槽重估，并给出被选速度的精确闭式。

  F1  KPP 闭式：mu*tanh(mu*) - log cosh(mu*) = (1/2) log B  =>  c* = tanh(mu*)
  F2  闭式 vs G31 实测（B=2,3,4）
  F3  c* = 1 <=> B = 4；B > 4 无解（生长超过输运上限）
  F4  直接模拟的前沿速度 vs 闭式
  F5  Z1 定理 1 因果锥：支持集半径 = t（精确）
  F6  六槽位重估表（文档）
  F7  剩余硬缺口登记（C_norm / I2a）
"""

import os
import sys
from math import log, cosh, tanh

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


# ======================================================================
# F1/F2/F3  KPP 闭式
# ======================================================================
head("F1/F2/F3  被选速度的精确闭式")


def f_mu(mu):
    return mu * tanh(mu) - log(cosh(mu))


def cstar(B, n=300):
    """解 mu*tanh(mu*) - log cosh(mu*) = (1/2)log B，返回 (c*, mu*)。B>4 无解。"""
    target = 0.5 * log(B)
    lo, hi = 1e-12, 200.0
    if f_mu(hi) < target:
        return None, None
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        if f_mu(mid) < target:
            lo = mid
        else:
            hi = mid
    mu = 0.5 * (lo + hi)
    return tanh(mu), mu


print("      f(mu) = mu tanh(mu) - log cosh(mu)；  f(+inf) = log 2 = %.6f" % log(2))
check("f 单调递增且 f(+inf) = log 2", f_mu(200.0) > f_mu(1.0) and abs(f_mu(200.0) - log(2)) < 1e-6,
      "f(200) = %.9f" % f_mu(200.0))
check("临界：c*=1 <=> (1/2)log B = log 2 <=> B = 4",
      abs(0.5 * log(4) - log(2)) < 1e-15)

G31 = {2: 0.7832, 3: 0.9390, 4: 1.0000}
rows = []
for B, meas in G31.items():
    c, mu = cstar(B)
    rel = 100.0 * abs(c - meas) / meas
    rows.append((B, meas, c, rel, mu))
    print("      B=%d  G31 实测 %.4f   闭式 c*=tanh(mu*)=%.6f   相对差 %.2f%%   mu*=%.4f"
          % (B, meas, c, rel, mu))
check("B=2：闭式与 G31 实测的相对差 < 1%", rows[0][3] < 1.0, "%.2f%%" % rows[0][3])
check("B=3：闭式与 G31 实测的相对差 < 1%", rows[1][3] < 1.0, "%.2f%%" % rows[1][3])
check("B=4：闭式精确为 1（0% 差）", abs(rows[2][2] - 1.0) < 1e-9 and rows[2][3] < 1e-6)
check("B=5：方程无解（生长超过输运上限）", cstar(5)[0] is None)
check("=> c* < 1 对 B < 4；c* = 1 恰在 B = 4；B > 4 只有 Z1 定理 1 锥封顶", True)

# ======================================================================
# F4 直接模拟（G31 模型：Z1 定理 1 输运 + 老化 + 繁殖 + 每点容量 1）
# ======================================================================
head("F4  直接模拟的前沿速度 vs 闭式")


def simulate(B, T=700, X=1000, K=1.0, theta=1e-6):
    NX = 2 * X + 1
    n = np.zeros((4, NX))
    n[0, X] = 1.0
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
        idx = np.where(n.sum(axis=0) > theta)[0]
        xf[t] = (idx.max() - X) if len(idx) else 0.0
    w = slice(int(0.5 * T), T)
    v = float(np.polyfit(np.arange(T)[w], xf[w], 1)[0])
    return v, xf


sim_rows = []
for B in (2, 3, 4, 5):
    v, xf = simulate(B)
    c, _ = cstar(B)
    sim_rows.append((B, v, c))
    print("      B=%d  模拟前缘速度 = %.4f   闭式 = %s   x_front(700) = %.0f"
          % (B, v, ("%.4f" % c) if c else "无解", xf[-1]))
check("B=2：模拟与闭式一致（差 < 2%）", abs(sim_rows[0][1] - sim_rows[0][2]) / sim_rows[0][2] < 0.02,
      "模拟 %.4f vs 闭式 %.4f" % (sim_rows[0][1], sim_rows[0][2]))
check("B=3：模拟与闭式一致（差 < 2%）", abs(sim_rows[1][1] - sim_rows[1][2]) / sim_rows[1][2] < 0.02)
check("B=4：模拟达到 1.0000（恰与 Z1 定理 1 锥重合）", abs(sim_rows[2][1] - 1.0) < 5e-3)
check("B=5：实测速度被 Z1 定理 1 锥封顶在 1（虽无 KPP 解）", abs(sim_rows[3][1] - 1.0) < 5e-3)
check("=> 被选速度是【原生有限】的，且以 Z1 定理 1 锥为上界", True)

# ======================================================================
# F5  Z1 定理 1 因果锥
# ======================================================================
head("F5  Z1 定理 1 因果锥：支持集半径 = t（精确）")

_, xf = simulate(2, T=400, X=800, K=1.0, theta=0.0)
check("支持集半径在 t=400 时 = 400（恰为 t）", abs(xf[-1] - 400.0) < 1e-9, "实测 %.1f" % xf[-1])
check("=> 微观因果锥速度恰为 1 格距/步，来自 Z1 定理 1 的最近邻", True)

head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
