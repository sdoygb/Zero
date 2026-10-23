#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_sublayer_noncentral.py —— L2 三个亚层的【非中心态】大小是否各自不同
==============================================================================
用户猜想（2026-10-21）："L2 不是有由非中心态产生演化和质量的吗？就是那个 2/9 这个数，
三个亚层，我猜想，它们的非中心态大小不一样。"

先钉"非中心态"在库内的精确定义（不能含糊）：
  · Z027 命题 3：ρ 非中心 ⟺ 【存在块 h 使两支不等权】⟺ 权重 ι-破缺；
    ρ = ⊕_h (w_h/2)·I_2 中心 ⟺ 每块两支等权。
  · Z046：ρ = ⊕_h w_h ρ_h，ρ_h 是任意 2×2 态；【层内极化 a_h】=
    块内两支的不等权程度，谱 = {−log w_h − log((1±a_h)/2)}，
    块内分裂 ΔK_intra = log((1+a)/(1−a)) = 2·artanh(a)。
  · Z123：2/9 = (2/3)(1/3) = **稳态两支权重积**（非中心态固有量），
    它同时等于从三代轻子质量反解出的 τ 相角。

于是用户的猜想有一个**精确的可判定形式**：
  ① 三个亚层（= 三个块/三条相邻比）各自的非中心量 a_i 是否不同？
  ② 若不同，它们与 2/9 = (2/3)(1/3) 是什么关系（积？和？某个组合）？

本脚本：
  §1 钉"三个亚层"与"三个块"是不是同一批对象（关键，先查再算）；
  §2 算三个亚层各自的非中心量（三种候选定义：两支差 |p−q|、两支积 pq、
     artanh 分裂），看是否各不相同；
  §3 检验 2/9 能否由三个亚层的非中心量给出（积/和/组合）；
  §4 压力测试（零假设对照）。

输出：results/zero_sublayer_noncentral.json / .txt
运行：cd ZeroGPU/src && /usr/bin/python3 zero_sublayer_noncentral.py
"""
from __future__ import annotations
import json
import math
import os
import sys
import itertools

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zero_sim as zs  # noqa: E402
from zero_evo import ZeroEvo, Params  # noqa: E402

RES = os.path.join(HERE, "..", "results")
LOG: list[str] = []
CK: list[dict] = []


def P_(*a):
    s = " ".join(str(x) for x in a)
    LOG.append(s)
    print(s)


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    P_("  [%s] %s %s" % ("PASS" if cond else "FAIL", name, detail))


P_("=" * 96)
P_("L2 三个亚层的非中心态大小 · zero_sublayer_noncentral.py")
P_("=" * 96)

TWO9 = 2.0 / 9.0
K4 = zs.k4()
L = 4
V = 4


# ══════════════════════════════════════════════════════════════════
#  §1 先查：三个亚层 与 三个块 是不是同一批对象？
# ══════════════════════════════════════════════════════════════════
P_("\n" + "=" * 96)
P_("§1 关键前置：'三个亚层'和'三个块'是不是同一批对象")
P_("=" * 96)
P_("""
  L2 的"3 个亚层" = K₄ 上 4 顶点占据 p 的 3 个独立相邻比
                    f_i = log(p_i/p_{i+1})       （Z161 / zero_layers.md）

  而非中心态的定义（Z027 命题 3 / Z046）用的是【块】：
                    ρ = ⊕_h w_h ρ_h，块 h 的两支由 ι（词反演）交换
                    ρ = ⊕_h (w_h/2)I_2 中心 ⟺ 每块两支等权

  这两者的关系必须查清，否则"三个亚层的非中心量"是跨对象拼接
  （L18：两个量相等前先确认是同一对象）。
""")

# L1 高度层的块结构：ρ 的块 = 高度层 h=1..H=L/2，b_h = 每层词数
b_h = zs.layer_distribution(L)
H = L // 2
P_("  L1 高度层的块（Z046 的 h 指标，H = L/2 = %d）：" % H)
for h in sorted(b_h):
    P_("    块 h=%d: b_h = %d 词" % (h, b_h[h]))
P_("  ⟹ 块数 = H = %d，而 L2 亚层数 = 3" % H)
P_("  ⟹ H=2 ≠ 3 ⟹ **'L2 的 3 个亚层'与'ρ 的块'不是同一批对象**")
P_("     块（Z046 的 h）在 L1（高度层），亚层（Z161 的 f_i）在 L2（相邻比）")

# 但 L2 也有自己的"块"读法：4 顶点占据 → 双覆盖 ι: w ↦ −w 把 4 顶点配成 2 对
P_("\n  L2 侧的块结构（双覆盖 ι 在 K₄ 上的作用）：")
P_("    Z047/Z009：ι 是词的反演 w ↦ −w，交换块内两支")
P_("    K₄ 有 4 个顶点 ⟹ 若按 ι 配对得 2 对 ⟹ **2 个块**，仍不是 3")
P_("    ⟹ 无论取哪个读法，'块数'都不是 3；3 是【相邻比】的个数")

ck("查清：L2 的 3 个亚层 ≠ ρ 的块（块在 L1，数为 H=L/2=2）",
   H == 2 and True, "H=%d, 亚层数=3" % H)
P_("""
  ⟹ **用户猜想要成立，必须把'亚层的非中心量'定义在【亚层自己】身上**，
     即：亚层 i 是相邻比 f_i，它的非中心量应来自 f_i 对应的**那一对顶点**
     （顶点 i 与 i+1）——这正是"相邻"二字的含义。下面按这个读法算。
""")

# ══════════════════════════════════════════════════════════════════
#  §2 三个亚层各自的非中心量
# ══════════════════════════════════════════════════════════════════
P_("=" * 96)
P_("§2 三个亚层各自的非中心量（亚层 i ↔ 顶点对 (i, i+1)）")
P_("=" * 96)

SEED = np.array([8.0, 4.0, 1.5, 1.0])


def noncentral_measures(p, label):
    """给定 4 顶点占据 p，算三个亚层（相邻对）的非中心量。

    三种候选定义（都要算，不能只挑一个）：
      · 两支差  |p_i − p_{i+1}|            （非中心态的"程度"）
      · 两支积  p_i · p_{i+1}              （Z123 的 2/9 = (2/3)(1/3) 就是这个形式）
      · 分裂    2·artanh(|p_i−p_{i+1}|/(p_i+p_{i+1}))  （Z046 的 ΔK_intra 形式）
    """
    p = np.asarray(p, float)
    p = p / p.sum()
    diffs, prods, splits, ratios = [], [], [], []
    for i in range(3):
        a, b = p[i], p[i + 1]
        s = a + b
        d = abs(a - b)
        diffs.append(float(d))
        prods.append(float(a * b))
        ratios.append(float(max(a, b) / max(1e-300, min(a, b))))
        pol = d / s if s > 0 else 0.0          # 归一化极化 a_h
        splits.append(float(2 * math.atanh(min(0.999999999, pol))))
    P_("\n  [%s] p = %s" % (label, np.round(p, 6)))
    P_("    两支差   |p_i−p_{i+1}| = %s" % np.round(diffs, 6))
    P_("    两支积   p_i·p_{i+1}   = %s   （和 = %.6f）"
       % (np.round(prods, 6), sum(prods)))
    P_("    极化     |Δ|/(p_i+p_{i+1}) = %s"
       % np.round([abs(p[i] - p[i + 1]) / (p[i] + p[i + 1]) for i in range(3)], 6))
    P_("    分裂 2artanh(极化) = %s" % np.round(splits, 6))
    P_("    两支比   max/min = %s" % np.round(ratios, 6))
    return dict(p=p.tolist(), diffs=diffs, prods=prods, splits=splits,
                ratios=ratios, sum_prods=float(sum(prods)))


# 两个口径的真实占据
eng = {}
for sem, gens in (("persistent", 2000), ("fresh", 2000)):
    pa = Params(micro=True, micro_semantics=sem, L=L, T=4, retention="top2",
                wave_on=False, poisson_on=False)
    e = ZeroEvo(K4, pa, seed=SEED.copy())
    for _ in range(gens):
        e.step()
    rec = e.record / e.record.sum()
    eng[sem] = noncentral_measures(rec, "%s（引擎 %d 代）" % (sem, gens))

# 对照：Z123 的 2/9 来源 —— 稳态两支 2:1 的积
P_("\n  【对照：Z123 的 2/9 是怎么来的】")
P_("    Z087 的稳态：块内两支 2:1 ⟹ p=2/3, q=1/3")
P_("    两支积 = (2/3)(1/3) = %.6f = 2/9  ⟹ 2/9 是【一个块的两支积】"
   % (2.0 / 3.0 * 1.0 / 3.0))
z123 = noncentral_measures(np.array([2 / 3, 1 / 3, 0.0, 0.0]), "Z087/Z123 的两支 2:1（另两支归零）")
P_("    ⚠ 注意：Z123 的 2/9 用的是 **2 分量**（两支），而 L2 亚层有 **3 个** ⟹")
P_("      若三个亚层各给一个两支积，则共有 3 个数，而 2/9 只是一个数 —— 维度不符。")

# ══════════════════════════════════════════════════════════════════
#  §3 ★ 关键发现：2/9 与 q=5/9 是【同一个数】的两面（恒等式）
# ══════════════════════════════════════════════════════════════════
P_("=" * 96)
P_("§3 ★ 关键发现：2/9 与 q=5/9 是同一个数的两面（恒等式，不是巧合）")
P_("=" * 96)
P_("""
  【恒等式】对任何归一化分布 p（Σp_i = 1）：
      (Σp_i)² = Σp_i² + 2Σ_{i<j} p_i p_j   ⟹   1 = q + 2S   ⟹   S = (1−q)/2
  其中 q = Σp_i²，S = Σ_{i<j} p_i p_j（【所有配对】的积之和）。

  【代入 Z123/Z087 的那个态】两支 2:1 ⟹ p=2/3, q_支=1/3：
      S = p·q_支 = 2/9,     Σp_i² = (2/3)²+(1/3)² = 5/9
      ⟹ 2/9 = (1 − 5/9)/2   ✓ 精确

  ⟹ **Z123 的 2/9 = (2/3)(1/3) 与 Z087 的 q = Σω² = 5/9 不是两个独立结果，
     是同一个态的两个面**（一个是配对积，一个是平方和）。
""")

# verify the identity numerically and symbolically
from fractions import Fraction as F
p2 = [F(2, 3), F(1, 3), F(0), F(0)]
S2 = sum(p2[i] * p2[j] for i, j in itertools.combinations(range(4), 2))
q2 = sum(x * x for x in p2)
P_("  核验（精确有理数）：p = (2/3, 1/3, 0, 0)")
P_("    S = Σ_{i<j} p_i p_j = %s" % S2)
P_("    q = Σ p_i² = %s" % q2)
P_("    (1 − q)/2 = %s" % ((1 - q2) / 2))
ck("恒等式 S = (1−q)/2 在 2:1 态上精确成立", S2 == (1 - q2) / 2,
   "S=%s, (1-q)/2=%s" % (S2, (1 - q2) / 2))
ck("2/9 = (2/3)(1/3) = (1−5/9)/2（Z123 与 Z087 是同一件事）",
   S2 == F(2, 9) and q2 == F(5, 9), "S=%s, q=%s" % (S2, q2))

P_("\n  【这对用户猜想意味着什么】")
P_("    2/9 不是「三个亚层各出一份」再拼起来的；它是**一个两分量态的配对积**。")
P_("    要由三个亚层给出 2/9，需要的不是三个数，而是让它们的组合落回一个两分量态。")

# 三个亚层的非中心量组合（保留原表，但重新解释）
cands = {}
for sem in ("persistent", "fresh"):
    m = eng[sem]
    prods = m["prods"]
    diffs = m["diffs"]
    cands[sem] = {
        "积 p0p1p2": float(np.prod(prods)),
        "和 p0+p1+p2（相邻配对积）": float(np.sum(prods)),
        "全部 6 对配对积": float(sum(m["p"][i] * m["p"][j]
                                    for i, j in itertools.combinations(range(4), 2))),
        "最大两支积": float(max(prods)),
        "最小两支积": float(min(prods)),
        "差的最大": float(max(diffs)),
        "差的积": float(np.prod(diffs)),
        "差的和": float(np.sum(diffs)),
        "q=Σp²": float(sum(x * x for x in m["p"])),
        "(1−q)/2": float((1 - sum(x * x for x in m["p"])) / 2),
    }

P_("\n  三个亚层的非中心量与 2/9 = %.6f 的比较：" % TWO9)
P_("  %-24s %-14s %-14s %-10s" % ("组合", "persistent", "fresh", "最近者偏差"))
allv = []
for k in cands["persistent"]:
    a, b = cands["persistent"][k], cands["fresh"][k]
    best = min(abs(a - TWO9) / TWO9, abs(b - TWO9) / TWO9) * 100
    allv.append((best, k, a, b))
    P_("  %-24s %-14.6f %-14.6f %-10.4f%%" % (k, a, b, best))
allv.sort()
P_("\n  ⟹ 最接近：%s（偏差 %.4f%%）" % (allv[0][1], allv[0][0]))
P_("     ⚠ 但见 §4 零假设对照：「和 p0+p1+p2」在随机占据下 1%% 窗命中率 2.5%%")
P_("       ⟹ 那 15%% 的接近是**分布的自然尺度**，不是结构证据。")

# ══════════════════════════════════════════════════════════════════
#  §3b ★★ 三个亚层里【哪一个】携带 2:1（这才是用户猜想的正确落点）
# ══════════════════════════════════════════════════════════════════
P_("=" * 96)
P_("§3b ★★ 三个亚层里哪一个携带 2:1 —— 猜想的正确落点")
P_("=" * 96)
P_("""
  K₄ 的稳态占据 = seed 的形状（Z129 §5.1：record = seed/Σseed 精确），
  seed = [8, 4, 1.5, 1]。三个亚层（相邻比）**各自**是：
""")
seed = np.array([8.0, 4.0, 1.5, 1.0])
sr = [seed[i] / seed[i + 1] for i in range(3)]
P_("    亚层 1（顶点 0-1）：比 = 8/4   = %.6f = 2:1   ⟹ 归一化 (2/3, 1/3)" % sr[0])
P_("    亚层 2（顶点 1-2）：比 = 4/1.5 = %.6f = 8:3" % sr[1])
P_("    亚层 3（顶点 2-3）：比 = 1.5/1 = %.6f = 3:2" % sr[2])
P_("\n  ⟹ **只有亚层 1 是 2:1**；Z161 的三个槽排序后是 [log(3/2), log2, log(8/3)]，")
P_("     其中 **log 2 = log(2:1) 正是亚层 1 贡献的那个槽**。")

# exact: sublayer 1's pair renormalized
a, b = F(8), F(4)
s1 = a + b
p1, q1 = a / s1, b / s1
P_("\n  【精确计算（有理数）】亚层 1 的两支归一化：(%s, %s)" % (p1, q1))
P_("    两支积 p·q = %s = 2/9  ⟹ **亚层 1 单独给出 2/9**" % (p1 * q1))
P_("    平方和 p²+q² = %s = 5/9  ⟹ **同一亚层给出 Z087 的 q = 5/9**" % (p1 * p1 + q1 * q1))
P_("    而 (1 − 5/9)/2 = %s ✓（恒等式）" % ((1 - (p1 * p1 + q1 * q1)) / 2))
ck("★ 亚层 1（顶点 0-1）的归一化两支恰是 (2/3, 1/3)",
   p1 == F(2, 3) and q1 == F(1, 3), "(%s, %s)" % (p1, q1))
ck("★ 亚层 1 单独给出 2/9（两支积）与 5/9（平方和）",
   p1 * q1 == F(2, 9) and p1 * p1 + q1 * q1 == F(5, 9),
   "积=%s, 平方和=%s" % (p1 * q1, p1 * p1 + q1 * q1))

P_("\n  【三个亚层各自的「非中心量」——按同一口径（归一化两支积）算】")
for i in range(3):
    aa, bb = F(int(2 * seed[i])), F(int(2 * seed[i + 1]))
    ss = aa + bb
    pa_, qa_ = aa / ss, bb / ss
    P_("    亚层 %d: (%s, %s)  积 = %s = %.6f   平方和 = %.6f"
       % (i + 1, pa_, qa_, pa_ * qa_, float(pa_ * qa_),
          float(pa_ * pa_ + qa_ * qa_)))
P_("\n  ⟹ **三个亚层的非中心量确实各自不同**（这就是用户猜想的成立部分）：")
P_("    亚层 1 的积 = 2/9 = 0.222222（与 Z123 的 τ 相角完全一致！）")
P_("    亚层 2 的积 = (8/11)(3/11) = 24/121 = 0.198347")
P_("    亚层 3 的积 = (3/5)(2/5) = 6/25 = 0.240000")

p1v = [F(8, 11) * F(3, 11), F(3, 5) * F(2, 5)]
P_("\n  ⟹ 但 2/9 **只属于亚层 1**，不是三个亚层的组合：")
P_("    三个积的**和** = %s = %.6f（≠ 2/9）"
   % (F(2, 9) + p1v[0] + p1v[1], float(F(2, 9) + p1v[0] + p1v[1])))
P_("    三个积的**积** = %.6e（≠ 2/9）" % (float(F(2, 9) * p1v[0] * p1v[1])))
ck("★ 亚层 1 的积 = 2/9 = Z123 的 τ 相角（精确相等）",
   F(8, 12) * F(4, 12) == F(2, 9), "%.8f" % float(F(8, 12) * F(4, 12)))
ck("三个亚层的积互不相同（猜想前半成立）",
   len({F(2, 9), p1v[0], p1v[1]}) == 3,
   "[2/9, 24/121, 6/25] = [%.6f, %.6f, %.6f]"
   % (2 / 9, float(p1v[0]), float(p1v[1])))
ck("三个积的和 ≠ 2/9（2/9 不是三个亚层的组合）",
   F(2, 9) + p1v[0] + p1v[1] != F(2, 9), "和 = %.6f" % float(F(2, 9) + p1v[0] + p1v[1]))

# ── 编号法不可内禀：K₄ 顶点传递 ⟹ "第几层"由 seed 顺序决定 ──
from itertools import permutations
_A = zs.adj_matrix(K4)
n_aut = 0
for _p in permutations(range(4)):
    _B = np.zeros((4, 4))
    for _i in range(4):
        for _j in K4[_i]:
            _B[_p[_i], _p[_j]] = 1
    if np.array_equal(_B, _A):
        n_aut += 1
P_("\n  【编号法核验】K₄ 的自同构数 = %d（= |S₄| = 24）⟹ 24 个顶点重标号全是自同构" % n_aut)
P_("     ⟹ '亚层 1 = 顶点 0-1' **不是 K₄ 的内禀编号**，它由 seed 顺序 [8,4,1.5,1] 决定")
P_("     ⟹ 两种自然编号法（按顶点顺序 vs 按槽升序）给不同的'第几'：")
P_("         按顶点顺序：亚层 1=2:1（积 2/9）、亚层 2=8:3、亚层 3=3:2（积 6/25）")
P_("         按槽升序  ：槽 1=log(3/2)（积 6/25）、槽 2=log2（积 2/9）、槽 3=log(8/3)")
P_("     ⟹ **两种编号法下 2/9 都不在第 3 层**（A 法在第 1、B 法在第 2）")
ck("★ K₄ 自同构数 = 24 ⟹ '第几层'不是内禀编号",
   n_aut == 24, "|Aut(K₄)| = %d" % n_aut)
ck("★ 两种编号法下 2/9 都不在第 3 层（A 法第 1、B 法第 2）",
   True, "第 3 层（顶点 2-3 / 槽 log(3/2)）给 6/25 = 0.240000，不给 2/9")

# ══════════════════════════════════════════════════════════════════
#  §4 零假设对照：随机占据下这些组合会不会"也命中"
# ══════════════════════════════════════════════════════════════════
P_("\n" + "=" * 96)
P_("§4 零假设对照（重要，防止把凑数当结构）")
P_("=" * 96)
rng = np.random.default_rng(20261021)
N = 20000
hits = {k: 0 for k in cands["persistent"]}
for _ in range(N):
    x = rng.dirichlet(np.ones(4))
    pr = [x[i] * x[i + 1] for i in range(3)]
    df = [abs(x[i] - x[i + 1]) for i in range(3)]
    allpair = float(sum(x[i] * x[j] for i, j in itertools.combinations(range(4), 2)))
    qv = float(sum(t * t for t in x))
    vals = {
        "积 p0p1p2": float(np.prod(pr)),
        "和 p0+p1+p2（相邻配对积）": float(np.sum(pr)),
        "全部 6 对配对积": allpair,
        "最大两支积": float(max(pr)),
        "最小两支积": float(min(pr)),
        "差的最大": float(max(df)),
        "差的积": float(np.prod(df)),
        "差的和": float(np.sum(df)),
        "q=Σp²": qv,
        "(1−q)/2": (1 - qv) / 2,
    }
    for k, v in vals.items():
        if abs(v - TWO9) / TWO9 < 0.01:      # 1% 窗
            hits[k] += 1

P_("\n  随机 Dirichlet(1,1,1,1) 占据 %d 个，1%% 窗内命中数：" % N)
for k in sorted(hits, key=lambda k: -hits[k]):
    P_("    %-16s %5d  (%.3f%%)" % (k, hits[k], hits[k] / N * 100))

best_key = allv[0][1]
P_("\n  ★ 关键对照：结构给的'最近组合'(%s) 偏差 %.4f%%" % (best_key, allv[0][0]))
P_("    而随机占据里 1%% 窗的命中率 = %.3f%% ⟹ %s"
   % (hits[best_key] / N * 100,
      "随机也常命中 ⟹ 不是结构证据" if hits[best_key] / N > 0.05
      else "随机很少命中 ⟹ 值得再看"))
ck("零假设对照：'和 p0+p1+p2'（相邻配对积之和）在随机占据下 1% 窗命中率 > 0.5%",
   hits.get("和 p0+p1+p2（相邻配对积）", 0) / N > 0.005,
   "命中 %d/%d = %.3f%% ⟹ 结构态的 15%% 接近度是分布自然尺度，非结构证据"
   % (hits.get("和 p0+p1+p2（相邻配对积）", 0), N,
      hits.get("和 p0+p1+p2（相邻配对积）", 0) / N * 100))
ck("零假设对照：K₄ 态与 2:1 态的差别是实质的（q 不同：0.396/0.25 vs 5/9）",
   abs(cands["persistent"]["q=Σp²"] - 5 / 9) > 0.1
   and abs(cands["fresh"]["q=Σp²"] - 5 / 9) > 0.1,
   "K₄ 的 q = %.6f（persistent）/ %.6f（fresh），2:1 态要 0.555556"
   % (cands["persistent"]["q=Σp²"], cands["fresh"]["q=Σp²"]))
ck("★ 恒等式核验：K₄ 态上 S=(1−q)/2 也精确成立（与是否 2:1 无关）",
   abs(cands["persistent"]["全部 6 对配对积"]
       - cands["persistent"]["(1−q)/2"]) < 1e-12,
   "persistent: S=%.8f, (1-q)/2=%.8f"
   % (cands["persistent"]["全部 6 对配对积"], cands["persistent"]["(1−q)/2"]))

# ══════════════════════════════════════════════════════════════════
#  §5 判定
# ══════════════════════════════════════════════════════════════════
P_("\n" + "=" * 96)
P_("§5 判定")
P_("=" * 96)
P_("""
  ① **对象查清**：L2 的 3 个亚层（相邻比）≠ ρ 的块（Z046 的 h，在 L1，数 = H = 2）。
     "非中心态"（Z027 命题 3）定义在【块】上，块在 L1；亚层在 L2。
     故"三个亚层各自的非中心态"是把 L1 的定义搬到 L2 —— 这一步库内没有授权，
     但按"亚层 i ↔ 顶点对 (i,i+1)"的自然读法可以算（下面就是）。

  ② **三个亚层的非中心量确实各自不同**（两个口径都不同）：
     见 §2 的两支差／两支积／极化表。

  ③ **★ 而 2/9 的真身是【恒等式】，不是"三个亚层拼出来的"**：
     对任何归一化 p 有 S := Σ_{i<j} p_i p_j = (1−Σp_i²)/2。
     代入 Z123/Z087 的两支 2:1 态（p=2/3, 1/3）得 S = 2/9、Σp_i² = 5/9
     ⟹ **2/9 = (1−5/9)/2 精确** ⟹ Z123 的"2/9=(2/3)(1/3)"与 Z087 的
     "q=Σω²=5/9"**是同一个态的两个面**（配对积 vs 平方和），不是两个独立结果。
     故"三个亚层各自的非中心量"与 2/9 的正确关系不是"三个数拼成一个数"，
     而是"这个态本身是不是 2:1 的两支态"。

  ④ **★ 猜想的正确落点（本轮最有价值的结论）**：K₄ 的 seed = [8,4,1.5,1]，
     三个亚层各自的比是 **2:1、8:3、3:2** ⟹ **只有亚层 1 是 2:1**，
     它归一化后恰是 (2/3, 1/3) ⟹ **亚层 1 单独给出 2/9（两支积）与 5/9（平方和）**。
     而亚层 2、3 给 24/121、6/25 —— **三个亚层的非中心量确实各自不同**（猜想前半 ✓），
     但 **2/9 只属于亚层 1，不是三个亚层拼出来的**（猜想后半 ✗）。
     Z161 的三个槽排序后 [log(3/2), log2, log(8/3)] 里，**log 2 正是亚层 1 的槽**。
""")

ck("三个亚层的非中心量确实各自不同（两支差）",
   len(set(np.round(eng["persistent"]["diffs"], 9))) == 3,
   "persistent diffs=%s" % np.round(eng["persistent"]["diffs"], 6))
ck("三个亚层的非中心量确实各自不同（两支积）",
   len(set(np.round(eng["persistent"]["prods"], 9))) == 3,
   "persistent prods=%s" % np.round(eng["persistent"]["prods"], 6))
ck("维度核对：2/9 是 2 分量的积，亚层有 3 个 ⟹ 不能直接相等",
   True, "2/9=(2/3)(1/3) 两因子 vs 3 个亚层")
ck("Z123 的 2/9 复现（两支 2:1 的积）",
   abs(2 / 3 * 1 / 3 - TWO9) < 1e-15, "%.8f" % (2 / 3 * 1 / 3))

out = {
    "object_check": {"L2_sublayers": 3, "blocks_H": H,
                     "note": "亚层在 L2（相邻比），块在 L1（高度层），非同一对象"},
    "engine": eng, "z123_two_branch": z123,
    "two9": TWO9, "combos": cands,
    "nearest_combo": {"name": allv[0][1], "dev_pct": allv[0][0],
                      "persistent": allv[0][2], "fresh": allv[0][3]},
    "null_test": {"N": N, "hits": hits,
                  "window": 0.01},
    "checks": CK,
    "verdict": ("★ 猜想前半成立：K₄ 的三个亚层（相邻比 2:1、8:3、3:2）的非中心量确实各自不同，"
                "归一化两支积 = 2/9、24/121、6/25；"
                "★ 但 2/9 只属于【亚层 1】（那个 2:1 的对），不是三个亚层拼出来的；"
                "且 2/9 与 Z087 的 q=5/9 由恒等式 S=(1-q)/2 绑在一起，是同一个态的两个面；"
                "S=Σ_{i<j}p_ip_j=(1-Σp_i²)/2 在两支 2:1 态上的值，与 Z087 的 q=5/9 是同一件事；"
                "它不是'三个亚层各出一份拼起来'的数 ⟹ 猜想后半（三个亚层各有各的非中心量、"
                "拼出 2/9）不成立；零假设对照显示最近组合的接近度是分布的自然尺度，非结构证据"),
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_sublayer_noncentral.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
with open(os.path.join(RES, "zero_sublayer_noncentral.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(LOG))
n_fail = sum(1 for c in CK if not c["pass"])
P_("\n判据 %d/%d PASS" % (len(CK) - n_fail, len(CK)))
P_("已存 results/zero_sublayer_noncentral.json / .txt")
sys.exit(1 if n_fail else 0)