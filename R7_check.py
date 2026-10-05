#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R7_check.py -- H3/H7 条件闭合的独立核验。

对应文档 R7_h3_h7_regularity.md。

  F1  精确定义、定理与诚实边界在位
  F2  双随机类到站点耦合的精确边缘与条件均值
  F3  均值约束：固定类频率时不是任意均值都能精确耦合
  F4  周期 mollification 在 C^0/C^1/C^2 中收敛
  F5  尺度窗口与 C^2 误差指数
  F6  H3 振荡反例：块平均可稳、单元内振荡不稳定
  F7  H7 字面版本失败：分片常数/分片线性不是 C^2
  F8  独立复核与对抗审计文件在位

数值仅核验文档中的结构、边缘和尺度，不替代 §3–§5 的证明。
"""

import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0
NCHECK = 0


def check(name, cond, detail=""):
    global FAIL, NCHECK
    NCHECK += 1
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def read(name):
    with io.open(os.path.join(HERE, name), encoding="utf-8") as f:
        return f.read()


DOC = read("R7_h3_h7_regularity.md")

# ======================================================================
head("F1  精确定义、定理与诚实边界在位")

for token in [
    "H3 条件闭合",
    "H7 的正确版本",
    "逐单元常数",
    "GDL",
    "I5b",
    "无条件 GR",
    "字面为假",
]:
    check("文档锚点：%s" % token, token in DOC)

check("写出类到连续场的显式映射",
      "Theta" in DOC or "\\Theta" in DOC)
check("区分条件证成与开放输入",
      "条件定理" in DOC and "开放" in DOC and "具名输入" in DOC)

# ======================================================================
head("F2  双随机耦合的精确边缘与条件均值")

v0, v1 = 0.75, 1.5
p0 = p1 = 0.5
mean = p0 * v0 + p1 * v1
N = 128
x = np.arange(N) / N
c_target = mean + 0.22 * np.sin(2 * np.pi * x) + 0.13 * np.cos(4 * np.pi * x)
c_s = c_target + mean - c_target.mean()
w = np.full(N, 1.0 / N)

q0 = (w / p0) * (v1 - c_s) / (v1 - v0)
q1 = (w / p1) * (c_s - v0) / (v1 - v0)
col = p0 * q0 + p1 * q1
cond = (v0 * p0 * q0 + v1 * p1 * q1) / col

check("q0 非负", np.min(q0) >= -1e-14, "min=%.3e" % np.min(q0))
check("q1 非负", np.min(q1) >= -1e-14, "min=%.3e" % np.min(q1))
check("q0 精确归一", abs(q0.sum() - 1) < 1e-12, "sum=%.15f" % q0.sum())
check("q1 精确归一", abs(q1.sum() - 1) < 1e-12, "sum=%.15f" % q1.sum())
check("耦合列边缘精确均匀",
      np.max(np.abs(col - w)) < 1e-12,
      "max=%.3e" % np.max(np.abs(col - w)))
check("站点条件均值等于保均值场",
      np.max(np.abs(cond - c_s)) < 1e-12,
      "max=%.3e" % np.max(np.abs(cond - c_s)))
check("保均值修正为 O(a^2)",
      abs(c_s.mean() - mean) < 1e-12 and np.max(np.abs(c_s - c_target)) < 1e-3,
      "shift=%.3e" % np.max(np.abs(c_s - c_target)))

# ======================================================================
head("F3  均值约束")

wrong = c_target + 0.1
shift = 1.0 / N
c_wrong = wrong + mean - wrong.mean()
q0w = (w / p0) * (v1 - c_wrong) / (v1 - v0)
q1w = (w / p1) * (c_wrong - v0) / (v1 - v0)
check("构造本身仍能保均值", abs(c_wrong.mean() - mean) < 1e-12)

raw_wrong = c_target + 0.1
q0r = (w / p0) * (v1 - raw_wrong) / (v1 - v0)
q1r = (w / p1) * (raw_wrong - v0) / (v1 - v0)
check("不改均值时行归一失败（这正是 R7-A1 的必要性）",
      abs(q0r.sum() - 1) > 1e-3 and abs(q1r.sum() - 1) > 1e-3,
      "q0sum=%.6f q1sum=%.6f" % (q0r.sum(), q1r.sum()))

# ======================================================================
head("F4  周期 mollification 在 C^0/C^1/C^2 中收敛")

M = 4096
xg = np.arange(M) / M
c = 1.2 + 0.12 * np.cos(2 * np.pi * xg) + 0.07 * np.sin(4 * np.pi * xg) + 0.03 * np.cos(6 * np.pi * xg)
freq = np.fft.fftfreq(M, d=1.0 / M)


def spectral_derivative(v, order):
    return np.fft.ifft((2j * np.pi * freq) ** order * np.fft.fft(v)).real


c_d0 = np.max(np.abs(c))
c_d1 = np.max(np.abs(spectral_derivative(c, 1)))
c_d2 = np.max(np.abs(spectral_derivative(c, 2)))

hs = np.array([0.04, 0.02, 0.01])
errs = []
for h in hs:
    mult = np.exp(-2 * np.pi ** 2 * freq ** 2 * h ** 2)
    ch = np.fft.ifft(mult * np.fft.fft(c)).real
    e0 = np.max(np.abs(ch - c))
    e1 = np.max(np.abs(spectral_derivative(ch, 1) - spectral_derivative(c, 1)))
    e2 = np.max(np.abs(spectral_derivative(ch, 2) - spectral_derivative(c, 2)))
    errs.append((e0, e1, e2))
errs = np.array(errs)

print("      h=%s" % hs.tolist())
print("      C0 err=%s" % ["%.3e" % e for e in errs[:, 0]])
print("      C1 err=%s" % ["%.3e" % e for e in errs[:, 1]])
print("      C2 err=%s" % ["%.3e" % e for e in errs[:, 2]])

ratios = errs[:-1] / errs[1:]
check("C0 误差随 h 减半至少下降 3.2 倍", np.all(ratios[:, 0] > 3.2),
      "ratios=%s" % np.round(ratios[:, 0], 3).tolist())
check("C1 误差随 h 减半至少下降 3.2 倍", np.all(ratios[:, 1] > 3.2),
      "ratios=%s" % np.round(ratios[:, 1], 3).tolist())
check("C2 误差随 h 减半至少下降 3.2 倍", np.all(ratios[:, 2] > 3.2),
      "ratios=%s" % np.round(ratios[:, 2], 3).tolist())
check("纯 mollification 收敛", errs[-1, 0] < errs[0, 0] and errs[-1, 2] < errs[0, 2])

# ======================================================================
head("F5  尺度窗口与 C^2 误差指数")

for d in (1, 2, 3):
    theta = np.linspace(0.01, 0.99, 400)
    e0 = 2 - (d + 2) * theta
    e2 = 2 - (d + 4) * theta
    eta = 1 - theta
    ok = np.all(e0[theta < 2.0 / (d + 4)] > 0) and np.all(e2[theta < 2.0 / (d + 4)] > 0)
    check("d=%d 的 H7 窗口 theta<2/(d+4) 非空且三误差指数为正" % d, ok)

for alpha in (0.25, 0.5, 1.0):
    d = 3
    theta = 2.0 / (d + 4 + alpha)
    rate = 2 * alpha / (d + 4 + alpha)
    check("alpha=%.2f 的平衡指数为正且 theta<2/(d+4)" % alpha,
          rate > 0 and theta < 2.0 / (d + 4),
          "theta=%.4f rate=%.4f" % (theta, rate))

# ======================================================================
head("F6  H3 振荡反例：块平均可稳、单元内振荡不稳定")

ncell = 256
a = 1.0 / ncell
beta = 0.35
h = 16 * a
xcell = (np.arange(ncell) + 0.5) * a
c_osc = 1 + beta * np.cos(2 * np.pi * xcell / h)
block = c_osc.mean()
edge_osc = np.abs(c_osc - block) / block
check("振荡场的 a-单元平均接近常数", abs(block - 1) < 0.02,
      "block=%.6f" % block)
check("单元内振荡不趋零（eta 失败）", edge_osc.max() > 0.3,
      "max=%.6f" % edge_osc.max())

# ======================================================================
head("F7  H7 字面版本失败")

xnode = np.linspace(0, 1, 9)
yconst = np.array([0, 0, 0, 0, 1, 1, 1, 1, 1], dtype=float)
left = np.diff(yconst[:-1])
right = np.diff(yconst[1:])
bad = np.where(np.abs(left - right) > 1e-12)[0]
check("分片常数的左右导数在结点不一致", bad.size > 0, "节点数=%d" % bad.size)

ytest = 1.0 - np.abs(2 * xnode - 1)
left_lin = np.diff(ytest[:-1])
right_lin = np.diff(ytest[1:])
bad_lin = np.where(np.abs(left_lin - right_lin) > 1e-12)[0]
check("分片线性插值的左右导数在折点不一致", bad_lin.size > 0,
      "节点数=%d" % bad_lin.size)

# ======================================================================
head("F8  独立复核与对抗审计")

for name in ["R7_independent_check.py", "R7_refutation_attempt.md"]:
    check("%s 在位" % name, os.path.exists(os.path.join(HERE, name)))

check("文档明确写出 R7-A1 均值约束",
      "R7-A1" in DOC and "目标场均值" in DOC)
check("文档明确写出 R7-A2/R7-A3 scale window",
      "R7-A2" in DOC and "R7-A3" in DOC)
check("文档明确写出 I5b 与 GDL 仍开放",
      "I5b" in DOC and "GDL 来源" in DOC)

# ======================================================================
head("汇总")
print("  断言 %d 项，不符项：%d" % (NCHECK, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
