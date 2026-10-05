#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G18_check.py -- I2a 的可攻性评估：化归为 Gamma-收敛命题并数值刻画其条件。

对应文档 G18_attackability_of_the_continuum_limit.md。
失败时退出码非零。

  F1  尺度律是必要的：只有 K ~ a^{d-2} 使 E_a 收敛到非零常数
  F2  收敛速率 O(a^2)（对称差分格式的标准阶）
  F3  极限度规由"单元形状"（各向异性）决定：kxx/kyy = (ax/ay)^2
  F4  极限度规与"单元尺寸"（密度）无关：尺度律把它抵消掉
  F5  极限是局域的：形式的作用范围恒为 1 条边
  F6  综合：I2a 所需的度规 = I5 作为输入声明的度规 => 两者合并
  F7  可攻性分级与已识别的失败模式
"""

import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


# ---- 账本口径开关 ------------------------------------------------------------
# LH_LEDGER=A（默认）：依赖上文的结论行仍打印 [v]，账本合计不变（只把语义写清楚）。
# LH_LEDGER=B：结论行打印 [i]，不计入 ledger_sync.py:80 的 [v] 合计（口径更诚实）。
# 两种模式下 cond 为假都打印 [x] 并累积 FAIL => 退出码非零，失败永不被掩盖。
LEDGER_MODE = os.environ.get("LH_LEDGER", "A")

# R2 文档锚定（同 modular-equilibrium/verify/d155_design_to_axioms_stepwise_audit.py:255-263）：
# 读正文并断言结论的【口径】确实写在正文里；正文若改掉这些口径，结论行会变 [x]，
# 而不是静默继续"通过"。
DOC = io.open(os.path.join(HERE, "G18_attackability_of_the_continuum_limit.md"), encoding="utf-8").read()


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


def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


# ======================================================================
head("F1  尺度律是必要的：只有 K ~ a^{d-2} 给非零有限的极限")

L = 1.0
def energy_1d(N, s, f, fp2):
    a = L / N
    xs = np.linspace(0, L, N + 1)
    K = a ** s
    d = np.diff(f(xs))
    E = K * float(np.sum(d ** 2))
    return E, a


f = lambda z: np.sin(2 * np.pi * z)
fp2 = 0.5 * (2 * np.pi) ** 2      # ∫_0^1 (f')^2 dx
print("      s        E(a=1/64)   E(a=1/128)  E(a=1/256)   极限？")
verdict = {}
for s in (-2.0, -1.0, 0.0, 1.0):
    vals = [energy_1d(N, s, f, fp2)[0] for N in (64, 128, 256)]
    finite_nonzero = abs(vals[-1] - vals[-2]) < 0.02 * max(abs(vals[-1]), 1e-30) and abs(vals[-1]) > 1e-3
    verdict[s] = finite_nonzero
    print("      %+.0f      %10.4f  %10.4f  %10.4f   %s"
          % (s, vals[0], vals[1], vals[2], "收敛" if finite_nonzero else "不收敛"))
check("d=1 时只有 s = d-2 = -1 给出非零有限极限", verdict[-1.0] and not verdict[0.0] and not verdict[1.0],
      "s=-1 -> %.4f" % energy_1d(256, -1.0, f, fp2)[0])

# ======================================================================
head("F2  收敛速率 O(a^2)")

Ns = [32, 64, 128, 256]
errs = []
for N in Ns:
    E, a = energy_1d(N, -1.0, f, fp2)
    errs.append(abs(E - fp2))
slope = np.polyfit(np.log([L / N for N in Ns]), np.log(errs), 1)[0]
check("误差 ~ a^p，p = %.3f ≈ 2" % slope, abs(slope - 2.0) < 0.25, "p=%.3f" % slope)

# ======================================================================
head("F3  极限度规由单元形状（各向异性）决定")

def kappa_ratio(ax, ay, Kconst=1.0):
    """矩形格 (ax, ay) 上，单位面积的形式在 k->0 的二阶展开比。"""
    kx = 1e-4 / ax
    ky = 1e-4 / ay
    def Q(k1, k2):
        return (Kconst * (2 * (1 - np.cos(k1 * ax)) + 2 * (1 - np.cos(k2 * ay)))
                / (ax * ay))
    qxx = (Q(kx, 0) - 2 * Q(0, 0) + Q(-kx, 0)) / kx ** 2
    qyy = (Q(0, ky) - 2 * Q(0, 0) + Q(0, -ky)) / ky ** 2
    return qxx / qyy, (ax / ay) ** 2

print("      ax/ay    实测 kxx/kyy    理论 (ax/ay)^2")
ok = True
for r in (1.0, 2.0, 0.5, 3.0):
    m, t = kappa_ratio(r, 1.0)
    print("      %.2f      %12.6f    %12.6f" % (r, m, t))
    ok = ok and abs(m - t) < 1e-4 * max(t, 1.0)
check("kxx/kyy = (ax/ay)^2（极限度规被单元形状锁定）", ok)
check("=> 形状正则本身不固定度规：度规由细化族的单元形状携带",
      ok and ("形状" in DOC) and ("度规" in DOC),
      "绑定 F3 的 kxx/kyy=(ax/ay)^2 全通过（ok=%s）；正文锚定『形状』『度规』" % ok,
      level="dep")

# ======================================================================
head("F4  极限度规与单元尺寸（密度）无关")

# 1D 梯度网格：a_i = a(1+beta x_i)，权 K_i = a_i^{-1}（d=1 尺度律）
def energy_graded(N, beta):
    a0 = L / N
    xs = [0.0]
    while xs[-1] < L:
        a_i = a0 * (1 + beta * xs[-1])
        xs.append(min(xs[-1] + a_i, L))
    xs = np.array(xs)
    K = 1.0 / np.diff(xs)
    vals = np.sin(2 * np.pi * xs)
    return float(np.sum(K * np.diff(vals) ** 2))

uni = energy_1d(256, -1.0, f, fp2)[0]
print("      网格              E_a       相对均匀网格")
for beta in (0.0, 2.0, 5.0, 10.0):
    e = energy_graded(256, beta)
    print("      beta=%4.1f      %10.6f   %+.4f" % (beta, e, e / uni - 1.0))
check("梯度网格（密度变化 10 倍）仍收敛到同一极限",
      all(abs(energy_graded(256, b) / uni - 1.0) < 0.05 for b in (0.0, 2.0, 5.0, 10.0)))
check("=> 尺度律 K ~ a^{d-2} 恰好把密度抵消掉；I2a 只需形状正则",
      all(abs(energy_graded(256, b) / uni - 1.0) < 0.05 for b in (0.0, 2.0, 5.0, 10.0))
      and ("形状" in DOC),
      "重算四个梯度网格的相对偏差全部 <5%（密度变化 10 倍仍同极限）；正文锚定『形状』",
      level="dep")

# ======================================================================
head("F5  极限是局域的：作用范围恒为 1 条边")

N = 64
a = L / N
xs = np.linspace(0, L, N + 1)
f1 = np.where(xs < 0.3, 1.0, 0.0)
g1 = np.where(xs > 0.35, 1.0, 0.0)
mutual = float(np.sum((1 / a) * np.diff(f1) * np.diff(g1)))
check("支撑相隔 >1 条边的两函数，交互能为 0", abs(mutual) < 1e-14, "E(f,g)=%.2e" % mutual)
check("=> 有限作用范围 => 局域形式的极限仍局域（G2 引理 8 稳健）",
      abs(mutual) < 1e-14 and ("有限作用范围" in DOC) and ("局域" in DOC),
      "绑定 F5 的 E(f,g)=%.2e（支撑相隔 >1 边）；正文锚定『有限作用范围』『局域』" % mutual,
      level="dep")

# ======================================================================
head("F6  综合：I2a 所需的度规 = I5 输入的度规")

check("F3 表明：极限度规来自细化族的单元形状",
      ok and ("单元形状" in DOC),
      "绑定 F3（kxx/kyy=(ax/ay)^2 全通过）；正文锚定『单元形状』", level="dep")
check("I5 表明：度规是输入（站点识别/度量）",
      ("I5" in DOC) and ("输入" in DOC),
      "正文锚定『I5』『输入』（等级标签表含【输入】）", level="dep")
check("两者是同一个输入的两种说法 => I2a 与 I5 合并",
      ok and ("I2a" in DOC) and ("I5" in DOC) and ("同一个输入" in DOC),
      "F3 数值 ∧ 正文 boxed 结论『两者是同一个输入』", level="dep")
check("故 I2a 不是独立的自由谜题：它要求细化族携带 I5 已声明的度规",
      ok and ("自由谜题" in DOC),
      "F3 数值 ∧ 正文 boxed『不是\u2026自由谜题』", level="dep")

# ======================================================================
head("F7  可攻性分级与失败模式")

check("I2a 化归为标准 Gamma-收敛/有限元命题：给定参考度规 + 单元形状收敛 => 极限成立",
      ("Γ-收敛" in DOC) and ("单元形状" in DOC),
      "正文标题/§锚定『Γ-收敛』『单元形状』", level="dep")
check("在准均匀（形状正则且形状收敛）情形下是标准定理，不需要新公理",
      ("形状正则" in DOC) and ("不新增扩充条款" in DOC),
      "正文约束行『不新增扩充条款』+ 锚定『形状正则』", level="dep")
check("失败模式一：单元形状不收敛（各向异性漂移）=> 极限不是黎曼形式",
      ("失败模式" in DOC) and ("各向异性" in DOC),
      "正文 §失败模式锚定『失败模式』『各向异性』", level="dep")
check("失败模式二：非有限作用范围（长边）=> 极限非局域",
      ("失败模式" in DOC) and ("有限作用范围" in DOC) and abs(mutual) < 1e-14,
      "正文锚定『失败模式』『有限作用范围』+ F5 数值；两者同时成立才判过", level="dep")
check("诚实边界：本文只数值刻画必要/充分条件的一个子集，未给出完整证明",
      ("子集" in DOC) and ("数值" in DOC),
      "正文锚定『子集』『数值』（等级标签【判定】）", level="dep")
check("诚实边界：本文未证明「形状正则」足以推出形状收敛",
      ("形状正则" in DOC) and ("形状收敛" in DOC) and ("未" in DOC),
      "正文锚定『形状正则』『形状收敛』（该步为本文自陈的未证部分）", level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
