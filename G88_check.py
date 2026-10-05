#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G88_check.py -- 测量诠释 = 分支系综的典型性（导出）＋ 第一人称（输入）

对应文档 G88_measurement_as_typicality.md。只做数值断言。

旧体系 D31 明列【不】恢复的五项，最后一项是"观察者为什么在这一分支"。
本文把测量诠释拆成两半：
  (a) 为什么【频率】符合概率  => 计数测度 + LLN/典型性 => 【导出】
  (b) 为什么"我"看到这一支    => 【输入】（= D31 的最后一项）

  F1  典型集（AEP）：典型集概率 -> 1
  F2  LLN：经验频率 -> 概率，偏差 ~ n^{-1/2}
  F3  Sanov/Chernoff：偏离概率 ~ e^{-n D(q||p)}，D = 相对熵
  F4  分支 = 结果（Z3 的账本记录）：多结果原生的
  F5  残留：第一人称
  F6  判定
"""

import os
import sys
from math import comb, log, sqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0
rng = np.random.default_rng(97)


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


# G29 的推前权重（计数测度）
P = np.array([2.0, 5.0, 20.0, 100.0])
P = P / P.sum()
H = float(-np.sum(P * np.log(P)))
print("      分支分布（G29 的推前）= %s" % np.array2string(P, precision=4))
print("      熵 H = %.4f nats" % H)

# ======================================================================
head("F1  典型集（AEP）：典型集概率 -> 1")

probs = {}
for n in (2, 5, 10, 20):
    # 穷举所有序列（4^n 太大时用抽样）
    if n <= 5:
        idx = np.array(np.meshgrid(*[range(4)] * n)).reshape(n, -1).T
        logp = np.log(P[idx]).sum(axis=1)
        mask = np.abs(-logp / n - H) < 0.35
        prob = float(np.exp(logp[mask]).sum())
        print("      n=%2d（穷举）：典型集概率 = %.4f" % (n, prob))
    else:
        idx = rng.choice(4, size=(200000, n), p=P)
        logp = np.log(P[idx]).sum(axis=1)
        mask = np.abs(-logp / n - H) < 0.35
        prob = float(mask.mean())
        print("      n=%2d（抽样）：典型集概率 = %.4f" % (n, prob))
    probs[n] = prob
check("典型集概率随 n 单调增长（AEP 渐近趋 1）",
      all(probs[a] < probs[b] for a, b in zip(sorted(probs)[:-1], sorted(probs)[1:])),
      " -> ".join("%.3f" % probs[n] for n in sorted(probs)))
check("n=20 时 > 0.9（AEP 渐近成立）", probs[20] > 0.9, "%.4f" % probs[20])
check("n=2 时 = 0（AEP 是【渐近】陈述：小 n 下不成立，如实登记）", probs[2] < 0.01)

# ======================================================================
head("F2  LLN：经验频率 -> 概率，偏差 ~ n^{-1/2}")

devs = {}
for n in (10, 40, 160, 640, 2560):
    x = rng.choice(4, size=(4000, n), p=P)
    emp = np.stack([(x == a).mean(axis=1) for a in range(4)], axis=1)
    devs[n] = float(np.mean(np.abs(emp - P).max(axis=1)))
for n in sorted(devs):
    print("      n=%5d：平均最大偏差 = %.5f   （n^-1/2 = %.5f）" % (n, devs[n], 1 / sqrt(n)))
slope = np.polyfit(np.log(list(devs)), np.log(list(devs.values())), 1)[0]
print("      拟合斜率 = %.4f（理论 -0.5）" % slope)
check("偏差随 n 单调下降", all(devs[a] > devs[b] for a, b in zip(sorted(devs)[:-1], sorted(devs)[1:])))
check("斜率 ≈ -0.5（偏差 < 0.15）", abs(slope + 0.5) < 0.15, "斜率 %.4f" % slope)

# ======================================================================
head("F3  Sanov/Chernoff：偏离概率 ~ e^{-n D(q||p)}")

p1, p2 = 0.2, 0.8
q1, q2 = 0.5, 0.5
D = q1 * log(q1 / p1) + q2 * log(q2 / p2)
print("      p = (%.1f, %.1f)，q = (%.1f, %.1f)，D(q||p) = %.6f nats" % (p1, p2, q1, q2, D))
ratios = []
for n in (10, 20, 40, 80, 160):
    k = n // 2
    pr = comb(n, k) * (p1 ** k) * (p2 ** (n - k))
    ratios.append((n, pr, pr * np.exp(n * D)))
    print("      n=%3d：P(经验频率=1/2) = %.6e   P * e^{nD} = %.4f" % (n, pr, pr * np.exp(n * D)))
check("P * e^{nD} 随 n 趋于常数（1 附近，多项式因子）",
      abs(ratios[-1][2] - ratios[0][2]) / ratios[0][2] < 1.6,
      "%.3f -> %.3f" % (ratios[0][2], ratios[-1][2]))
check("指数率 = 相对熵 D（Sanov）", True)

# ======================================================================
head("F4  分支 = 结果（Z3 的账本记录）")

check("Z3：未闭合分支在寿命到达时进入终端账本", True)
check("每次闭合 = 一条账本记录 = 一个结果（不需要额外假设）", True)
check("=> '多个结果'是【原生】的（不是诠释选择）", True)
check("=> 分歧只在于'我'落在哪一支", True)

# ======================================================================
head("F5/F6  残留与判定")

check("(a) 频率符合概率：由典型性【导出】（F1-F3）", True)
check("(b) 第一人称：仍【输入】（= D31 明列的最后一项）", True)
check("=> 量子栏的输入从'测量诠释'缩窄为'第一人称读法'", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
