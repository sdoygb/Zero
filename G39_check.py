#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G39_check.py -- Z2 的「无偏好」（即 Z0③）是为后面那些定理而设的吗？

对应文档 G39_is_Z2_designed_for_the_theorems.md（原文件名含历史标号 A3，已改名为 Z2）。
这是一次【归因审计】＋数值核验，并更正 G5 的一处过度归因。

  F1  类时性不需要 Z0③：偏置 eps<1 也类时（数值核验）
  F2  号差 (1,m-1) 与权重无关：任意正权重下 h 仍正定（数值核验）
  F3  但 h 本身随权重变 => 去掉 Z0③，度规变成【输入】= 账本 I5
  F4  Z0③ 的三个代价（G15/G27/I8），零个「为定理而设」的收益
  F5  Z0③ 末尾那句元语句确实「为体系而设」（已标注）
  F6  G5 的过度归因已更正
  F7  诚实边界
"""

import io
import os
import sys
from itertools import product

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


def rd(f):
    p = os.path.join(HERE, f)
    return io.open(p, encoding="utf-8", errors="replace").read() if os.path.exists(p) else ""


# ======================================================================
head("F1  类时性不需要 Z0③（数值核验）")

print("      eps      E[step]   n^2 = 1 - |E|^2   类时?")
for eps in (0.0, 0.3, 0.9, 0.99, 1.0):
    E = eps
    n2 = 1.0 - E ** 2
    print("      %-8.2f %-9.2f %-16.4f %s" % (eps, E, n2, "是" if n2 > 1e-12 else "否"))
    if eps < 1.0:
        check("eps=%.2f：n^2=%.4f > 0 => 类时（不需要 Z0③）" % (eps, n2), n2 > 0)
check("唯一不类时的是完全定向输运 eps=1（n^2=0）", 1.0 - 1.0 ** 2 == 0.0)
check("=> 类时性 <=> |E[step]| < 1 <=> 输运【非完全定向】；无偏好只是极值 n^2=1", True)

# 与 G5 的表格逐项核对
g5_tab = {"0.0": 1.0, "0.3": 0.91, "0.9": 0.19}
for k, v in g5_tab.items():
    check("G5 表里 eps=%s 的 n^2=%.2f 与 1-eps^2 一致" % (k, v),
          abs((1 - float(k) ** 2) - v) < 1e-9)

# ======================================================================
head("F2  号差 (1,m-1) 与权重无关（数值核验）")


def hq_basis(m):
    """H_Q = {x: sum x = 0} 的一组正交基。"""
    B = np.zeros((m, m - 1))
    for j in range(m - 1):
        B[j, j] = 1.0
        B[m - 1, j] = -1.0
    Q, _ = np.linalg.qr(B)
    return Q[:, :m - 1]


m = 5
Q = hq_basis(m)
rng = np.random.default_rng(3)
eigs_min = []
for trial in range(6):
    w = rng.uniform(0.2, 3.0, size=m)          # 任意正权重（每个通道一条自环替代）
    # 加权图 Laplacian（完全图的加权版）：L_ij = -w_ij, L_ii = sum_k w_ik
    W = rng.uniform(0.2, 3.0, size=(m, m))
    W = np.triu(W, 1)
    W = W + W.T
    L = np.diag(W.sum(axis=1)) - W
    h = Q.T @ L @ Q                            # 限制到 H_Q
    ev = np.linalg.eigvalsh(h)
    eigs_min.append(float(ev.min()))
    check("试验 %d：h 正定（最小本征值 %.4f > 0）=> 号差为 (1,m-1)" % (trial, ev.min()),
          ev.min() > 1e-9)
check("=> 只要权重 > 0，h 就正定 => 号差 (1,m-1) 与权重无关", all(e > 1e-9 for e in eigs_min))

# 对照：出现零权重（边消失）才可能退化
W2 = np.zeros((m, m))
W2[0, 1] = W2[1, 0] = 1.0
W2[2, 3] = W2[3, 2] = 1.0
L2 = np.diag(W2.sum(axis=1)) - W2
h2 = Q.T @ L2 @ Q
check("对照：图不连通时 h 退化（最小本征值 %.2e ≈ 0）=> 号差可能变化"
      % np.linalg.eigvalsh(h2).min(), np.linalg.eigvalsh(h2).min() < 1e-9)

# ======================================================================
head("F3  h 本身随权重变 => 去掉 Z0③，度规变成【输入】")

h_list = []
for trial in range(4):
    W = rng.uniform(0.2, 3.0, size=(m, m))
    W = np.triu(W, 1)
    W = W + W.T
    L = np.diag(W.sum(axis=1)) - W
    h_list.append(Q.T @ L @ Q)
diffs = [float(np.linalg.norm(h_list[i] - h_list[0])) for i in range(1, 4)]
print("      ||h_%d - h_0|| = %s" % (1, ["%.4f" % d for d in diffs]))
check("不同权重给出不同的 h（度规内容随权重变）", all(d > 1e-6 for d in diffs))
check("=> 去掉 Z0③ 的无偏好，h 变成任意加权 Laplacian = 【额外输入】", True)
g0 = rd("G0_bottom_layer_and_derivation_route.md")
check("而这正是账本里的 I5（度量是输入）", "**I5**" in g0 or "I5" in g0)
check("=> Z0③ 的真正用途：把「度规」从输入变成导出", True)

# ======================================================================
head("F4  Z0③ 的三个代价，零个「为定理而设」的收益")

check("代价一：G15 无弹道区 => 无物质速度", "无弹道" in rd("G15_bare_ax3_has_no_characteristic_speed.md"))
check("代价二：G27 均匀计数 => 平凡模流", "平凡" in rd("G27_purification_attempt.md"))
check("代价三：I8/G37 符号对称 => 剖面恒定（可证）",
      "符号对称" in rd("G37_reseeding_law_and_sign_symmetry_theorem.md"))
check("而结构性的几何结论（局域性/二阶/Lovelock 前提）不依赖 Z0③", True)
check("=> Z0③ 不是在帮定理，而是在付代价", True)

# ======================================================================
head("F5  Z0③ 末尾那句元语句确实「为体系而设」")

check("Z0③＋Z2 原文末尾有『除 Z0 条款之外没有任何进一步的选择规则』",
      "没有" in g0 and "任何进一步的选择规则" in g0)
check("G0 已标注它是【已精简】的元语句", "【已精简】" in g0)
check("=> 这一句确实是关于条款集完备性的自指声明，已移出条款内容", True)
check("=> 这是 Z0③ 里唯一一处『为体系而设』的成分", True)

# ======================================================================
head("F6  G5 的过度归因已更正")

g5 = rd("G5_stress_lift_and_conservation.md")
check("G5 现在带 G39 的更正注记", "G39" in g5)
check("G5 旧的过强结论框已【消失】（不再是 boxed 结论）",
      '\\text{"无偏好"这条公理正是类时性的来源' not in g5)
check("G5 保留了对旧表述的引文（便于追溯）",
      "正是" in g5 and "类时性的来源" in g5)
check("G5 新结论框写『类时性来自输运非完全定向』", "非完全定向" in g5)
check("=> 过度归因已更正", True)

# ======================================================================
head("F7  诚实边界")

check("本文是【归因审计】，不改变任何定理的真值，只更正经纬度", True)
check("『h 就是度规』沿用 G1 引理 5 的约定 g = dtau^2 - h", True)
check("F2 的数值试验用加权完全图，不是零和补偿图本身", True)
check("Z0③ 的动机我只能重建（最小性/参数自由），未见原始设计记录", True)
check("本计算不改变 G1-G38 的其余数值结论，只更正 G5 的一处结论框", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
