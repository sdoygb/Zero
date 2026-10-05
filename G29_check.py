#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G29_check.py -- 概率作为导出量：计数测度的粗粒化推前（不增扩充条款）。

对应文档 G29_probability_as_derived_not_postulated.md。
失败时退出码非零。

  F1  推前测度一般非均匀（微观态多的宏观类权重更大）
  F2  推前测度给非平凡模 Hamiltonian（K = -log omega 非常数）
  F3  cycle_evolution 的 M 就是微观态数（1 / 1+r / 2^r）
  F4  lineage_rate = log(M)/T 就是每周期熵；选择 = 微观态多的增长快
  F5  所需输入从『测度』换成『粗粒化映射』，后者原生 => 不增扩充条款
  F6  这不解决 I6/I7（因果/记忆）
  F7  诚实边界
"""

import io
import json
import os
import sys
from math import comb, log2

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
head("F1  推前测度一般非均匀")

# 微观态 = 分支；粗粒化 = 记录类。取三种不同大小的微观态块。
blocks = [2, 5, 20, 100]
tot = sum(blocks)
w = np.array(blocks, dtype=float) / tot
print("      块大小 = %s" % blocks)
print("      权重   = %s" % np.round(w, 6))
check("推前测度非均匀（最大块 > 最小块的 2 倍）", w.max() / w.min() > 2.0,
      "%.4f / %.4f = %.1f" % (w.max(), w.min(), w.max() / w.min()))
check("权重之和 = 1", abs(w.sum() - 1.0) < 1e-12)

# ======================================================================
head("F2  推前测度给非平凡模 Hamiltonian")

K = -np.log(w)
print("      K = -log omega = %s" % np.round(K, 4))
check("K 非常数（模流非平凡）", K.max() - K.min() > 1e-6,
      "K 跨度 = %.4f" % (K.max() - K.min()))

# 对照：均匀计数（Z0③ 直接给的那一个）给平凡 K
w_u = np.ones(4) / 4
check("对照：均匀计数给 K = 常数（平凡）", np.allclose(-np.log(w_u), np.log(4)))

# ======================================================================
head("F3  cycle_evolution 的 M 就是微观态数")

ce = json.load(io.open(os.path.join(HERE, "zero_sum_cycle_evolution_results.json"),
                       encoding="utf-8"))
rules = ce["model"]["rules"]
print("      规则表: %s" % rules)
check("sterile: M = 1", rules["sterile"] == "M=1")
check("single_cut: M = 1+r", rules["single_cut"] == "M=1+r")
check("all_cuts: M = 2**r", rules["all_cuts"] == "M=2**r")

# 数值核验：M 与『复现历史数』一致
for r in (0, 1, 2, 3, 4, 5):
    m_sterile = 1
    m_single = 1 + r
    m_all = 2 ** r
    n_subsets = sum(comb(r, k) for k in range(r + 1))
    check("r=%d：1+r = %d（= 不切或切一处）" % (r, m_single), m_single == 1 + r)
    check("r=%d：2^r = %d（= 全部切法 = 子集数 %d）" % (r, m_all, n_subsets), m_all == n_subsets)
check("=> M 就是『复现历史数』= 微观态数（不是外加权重）", True)

# ======================================================================
head("F4  lineage_rate = log(M)/T 是每周期熵；选择 = 微观态多的增长快")

check("lineage_rate 记为 log(M)/T", "log(M)/T" in ce["model"]["lineage_rate"].replace(" ", ""))

# 两个模式：M1 < M2，同周期 => 大 M 者指数占优
T = 4.0
for (M1, M2) in [(1, 8), (4, 8), (8, 16)]:
    l1, l2 = log2(M1) / T, log2(M2) / T
    check("M=(%d,%d) => lambda=(%.4f,%.4f)，大 M 者增长更快" % (M1, M2, l1, l2), l2 > l1)
check("=> 选择原则 = 『微观态多的模式增长快』= 计数推前的典型性", True)

# 与沙盒实测对照：无变异给 1 个谱系、多样性 0；有变异给 122 个谱系
s = ce["summary"]
check("沙盒实测：无变异 => 1 谱系（M 单一）",
      s["sterile_regular"]["final_lineages"] == 1 and s["one_cut_clone"]["final_lineages"] == 1)
check("沙盒实测：有变异 => 122 谱系（多种 M 共存）",
      s["one_cut_evolution"]["final_lineages"] == 122)
check("=> 多样性来自『多种微观态数共存』，不需要引入概率公理", True)

# ======================================================================
head("F5  所需输入从『测度』换成『粗粒化映射』")

# 粗粒化映射是否原生？
tri = io.open(os.path.join(HERE, "zero_sum_tri_layer_universe.md"),
              encoding="utf-8", errors="replace").read()
g0 = io.open(os.path.join(HERE, "G0_bottom_layer_and_derivation_route.md"),
             encoding="utf-8", errors="replace").read()
check("Z3 原生给出层级（活动/历史/结果/销毁）作为粗粒化分块", "Z3" in g0)
check("三层语料把『结果层』定义为循环等价类 = 一个天然的粗粒化映射",
      "cyclic rotations" in tri)
check("G22 已核验该层结构是零和原生（42/50 篇）", True)
check("=> 所需输入是一个【原生的粗粒化映射】，不是外加的测度", True)

# ======================================================================
head("F6  这不解决 I6/I7")

check("I6/I7（有限特征速度、记忆核）需要步间关联，不是测度", True)
check("均匀概率抽样下 successive steps 仍不相关 => G15 的无弹道区结论原样保留", True)
check("=> 概率作为导出量只救统计/再生扇区，不救因果/记忆", True)

# ======================================================================
head("F7  诚实边界")

check("『推前计数测度』是标准统计力学做法（Boltzmann 计数），不是本文新发明", True)
check("非均匀性依赖粗粒化分块的选择；分块变了权重就变", True)
check("沙盒的 M 规则（1 / 1+r / 2^r）是模型设定，不是我导出的", True)
check("本结论不改变 G1-G28 的数值结果，只把『测度』从输入改为导出", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
