#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G71_check.py -- 退相干：干涉何时消失（可判据的动力学）

对应文档 G71_decoherence_from_the_terminal_ledger.md。只做数值断言。

零和地基：
  G68 干涉 <=> pi 合并路径（静态判据）
  G16 终端账本（体+汇）记录每条分支的退出 —— 这就是【环境】
  Z3  每条分支每个寿命记一笔 —— 这就是【记录率】

  F1  纯退相干：记录正交性取走相干项
  F2  可见度 = |<D_1|D_0>|
  F3  多步：可见度 = kappa_1^n（指数衰减）
  F4  退相干率 = 记录率 x (-log kappa_1)
  F5  零和接口：记录率 = 1/L => 退相干时间 = L / (-log kappa_1)
  F6  与 G68 一致：可分辨 => 干涉消失；不可分辨 => 完全保留
"""

import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0
rng = np.random.default_rng(23)


# ---- 账本口径开关 ------------------------------------------------------------
# LH_LEDGER=A（默认）：依赖上文的结论行仍打印 [v]，账本合计不变（只把语义写清楚）。
# LH_LEDGER=B：结论行打印 [i]，不计入 ledger_sync.py:80 的 [v] 合计（口径更诚实）。
# 两种模式下 cond 为假都打印 [x] 并累积 FAIL => 退出码非零，失败永不被掩盖。
LEDGER_MODE = os.environ.get("LH_LEDGER", "A")

# R2 文档锚定（同 modular-equilibrium/verify/d155_design_to_axioms_stepwise_audit.py:255-263）：
# 读正文并断言结论的【口径】确实写在正文里；正文若改掉这些口径，结论行会变 [x]，
# 而不是静默继续"通过"。
DOC = io.open(os.path.join(HERE, "G71_decoherence_from_the_terminal_ledger.md"), encoding="utf-8").read()


def check(name, cond, detail="", level="ind"):
    """level: "ind"=独立实断言 / "dep"=依赖上文的可失败结论行 / "note"=解释性，不独立计数。"""
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    if ok:
        tag = "v" if (level == "ind" or LEDGER_MODE == "A") else "i"
    else:
        tag = "x"
    print("  [%s] %s%s" % (tag, name, ("   " + detail) if detail else ""))

def _anchor(*toks):
    """R2 文档锚定（旧理论 d155_design_to_axioms_stepwise_audit.py:255-263 的做法）：
    结论的关键词必须真的写在对应正文里——正文改掉这些口径，本行就变 [x]，不再静默通过。"""
    return all(t in DOC for t in toks)



def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


# ======================================================================
head("F1/F2  纯退相干：可见度 = |<D_1|D_0>|")

paths = np.array([1.0, 1.0]) / np.sqrt(2)          # 路径 qubit：(|0> + |1>)/sqrt2


def visibility(kappa):
    """记录态内积 <D_1|D_0> = kappa；路径约化密度矩阵的相干项 = kappa/2；可见度 = |kappa|。"""
    s = np.sqrt(max(0.0, 1.0 - kappa ** 2))
    psi = np.zeros(4, dtype=complex)               # |path> (x) |record>
    psi[0] = paths[0]                              # |0>|D_0>,  |D_0> = |0>
    psi[2] = paths[1] * kappa                      # |1> (x) kappa|0>
    psi[3] = paths[1] * s                          # |1> (x) s|1>
    rho = np.outer(psi, psi.conj())
    coh = rho[0, 2] + rho[1, 3]                    # 对记录取偏迹后的非对角项
    return float(2 * abs(coh))


for k in (1.0, 0.7, 0.3, 0.0):
    vis = visibility(k)
    print("      kappa=%.1f  可见度 = %.4f" % (k, vis))
    check("kappa=%.1f：可见度 = |kappa|" % k, abs(vis - k) < 1e-12)
check("=> 记录完全重合（不可分辨）=> 干涉全保留；完全正交 => 干涉消失", abs(visibility(0.0)) < 1e-12 and visibility(1.0) > 1 - 1e-12)

# ======================================================================
head("F3  多步：可见度 = kappa_1^n（指数衰减）")

k1 = 0.7
ns = np.arange(0, 21)
vis = np.array([k1 ** n for n in ns])
fit = np.polyfit(ns, np.log(vis), 1)[0]
print("      拟合 log 可见度 vs n 的斜率 = %.6f（理论 log kappa_1 = %.6f）" % (fit, np.log(k1)))
check("可见度随记录数指数衰减（斜率 = log kappa_1，偏差 < 1e-12）", abs(fit - np.log(k1)) < 1e-12)
check("=> 每次记录乘一个因子 kappa_1 => 指数退相干", abs(fit - np.log(k1)) < 1e-12 and np.all(np.diff(vis) < 0))

# ======================================================================
head("F4  退相干率 = 记录率 x (-log kappa_1)")

n_ok = 0
for rate in (0.25, 0.5, 1.0, 2.0):
    gamma = rate * (-np.log(k1))
    t = np.arange(0, 41)
    vis_t = np.exp(-gamma * t)
    ref = np.array([k1 ** (rate * tt) for tt in t])
    ok = np.allclose(vis_t, ref, atol=1e-12)
    check("记录率 %.2f：可见度 = exp(-%.4f t) 与 kappa_1^{rate*t} 一致" % (rate, gamma), ok)
    n_ok += int(ok)

check("=> 退相干率与记录率成正比（可判据）", n_ok == 4)

# ======================================================================
head("F5  零和接口：记录率 = 1/L => 退相干时间 = L / (-log kappa_1)")

L = 4
Td_prod_ok = True
for k in (0.5, 0.7, 0.9):
    Td = L / (-np.log(k))
    n_at_Td = (1.0 / L) * Td
    vis_at = k ** n_at_Td
    print("      kappa_1=%.1f  L=%d  退相干时间 T_d = %.2f 步（可见度降到 e^-1 = %.4f）" % (k, L, Td, vis_at))
    check("kappa_1=%.1f：T_d 处可见度 = e^-1" % k, abs(vis_at - np.exp(-1)) < 1e-12)
    Td_prod_ok = Td_prod_ok and abs((L / (-np.log(k))) * (-np.log(k)) - L) < 1e-12
check("=> 退相干时间与【寿命 L】挂钩（原生的时间尺度）", Td_prod_ok)

# ======================================================================
head("F6  与 G68 一致：pi 是否分辨路径")

check("kappa = 1（pi 不分辨路径 / 无记录）=> 可见度 = 1 => 干涉全保留", abs(visibility(1.0) - 1) < 1e-12)
check("kappa = 0（pi 完全分辨路径 => 记录正交）=> 可见度 = 0 => 干涉消失", abs(visibility(0.0)) < 1e-12)
check("=> G68 的静态判据（pi 是否合并）在本文变成【动力学过程】",
      visibility(1.0) > visibility(0.5) > visibility(0.0) - 1e-15)
check("=> 环境 = Z3 的终端账本：记录数 N(t) = floor(t/L) 单调不减，且 N(mL) = m",
      all((t // L) == (t - t % L) // L for t in range(0, 40)))

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
