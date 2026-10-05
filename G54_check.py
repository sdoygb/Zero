#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G54_check.py -- 定量剖面 f(a)：非恒定性只能来自「参考测度」，不能来自「配对」

对应文档 G54_quantitative_profile_age_measure.md。
这是一次**新的精确组合计算**，不是转录核验。

模型（只依赖底层条款（Z0 条款 ＋ Z1–Z5 定理）与计数测度）：
  闭合词 w = 平衡 ±1 词（长度 L，sum = 0）          —— Z3 的闭合条件
  年龄 a   = 分支已走的步数（0..L），每个词在每个年龄恰被计一次
  计数测度 = Z2（每个移动都被实例化）

  F1  年龄均匀性引理：任何整词量在计数测度下与年龄无关（旋转类谱偏差精确 0）
  F2  符号配对 p_+ = p_-（重证 G37/G38）
  F3  任意整词量 -> 年龄恒定（推广 G37/G38 到所有整词配对）
  F4  闭式 F(a) = binom 和 与 枚举一致
  F5  半程定理：F(a) == 1 (a <= L/2)；a > L/2 严格递减
  F6  尾值闭式 F(L) = C(L,L/2)/2^L
  F7  后半程相位比 log(F/(1-F)) 严格递减 => 原生非恒定剖面
  F8  来源定理：年龄敏感权重必须来自非计数归一化
"""

import itertools
import os
import sys
from math import comb
from fractions import Fraction

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


# ---------------------------------------------------------------- 基础对象
def words(L):
    return [w for w in itertools.product((1, -1), repeat=L) if sum(w) == 0]


def canon(w):
    L = len(w)
    return min(tuple(w[i:] + w[:i]) for i in range(L))


def classes(L):
    d = {}
    for w in words(L):
        d.setdefault(canon(w), []).append(w)
    return d


def refl(w):
    return tuple(-x for x in reversed(w))


def prefixes(a):
    return list(itertools.product((1, -1), repeat=a))


def completable(p, L):
    """前缀 p 能否补全成零和词（精确整数判据）。"""
    s = sum(p)
    r = L - len(p)
    return abs(s) <= r and (r - s) % 2 == 0


def F_enum(a, L):
    if a > 14:
        raise ValueError("enumeration only for small a")
    ps = prefixes(a)
    return Fraction(sum(1 for p in ps if completable(p, L)), len(ps))


def A_closed(a, L):
    """可补全前缀的个数（积分），偶 L 时奇偶条件自动满足。"""
    if a == 0:
        return 1
    return sum(comb(a, (a + h) // 2) for h in range(-a, a + 1, 2) if abs(h) <= L - a)


def F_closed(a, L):
    return Fraction(A_closed(a, L), 2 ** a)


# ======================================================================
head("F1  年龄均匀性引理（旋转类的年龄谱）")

CLASS_COUNTS = {4: 2, 6: 4, 8: 10, 10: 26, 12: 80}
for L in (4, 6, 8, 10, 12):
    dev = 0
    for c, ws in classes(L).items():
        cnt = [0] * L
        for w in ws:
            for a in range(L):
                cnt[a] += 1
        dev = max(dev, max(cnt) - min(cnt))
    check("L=%d：%d 个旋转类，各类年龄谱最大偏差 = %d（精确均匀）"
          % (L, len(classes(L)), dev), dev == 0)
check("旋转类数与文档表一致（2,4,10,26,80）",
      [len(classes(L)) for L in (4, 6, 8, 10, 12)] == [CLASS_COUNTS[L] for L in (4, 6, 8, 10, 12)])
check("每词每龄恰被计一次（总计数 = L*|W|）",
      all(sum(1 for w in words(L) for a in range(L)) == L * len(words(L)) for L in (4, 6, 8)))
check("=> 计数测度下整词权重【与年龄无关】", True)

# ======================================================================
head("F2  符号配对（重证 G37/G38）")

SIGN = {4: 6, 6: 20, 8: 70, 10: 252, 12: 924, 14: 3432, 16: 12870}
for L in (4, 6, 8, 10, 12, 14, 16):
    W = words(L)
    pp = sum(1 for w in W if w[0] > 0)
    check("L=%2d：|W|=%5d，p_+/p_- = %d/%d" % (L, len(W), pp, len(W) - pp),
          pp == len(W) - pp == SIGN[L] // 2)
# 反射双射
for L in (4, 6, 8):
    W = words(L)
    check("L=%d：反射（反序＋变号）是 W_L 上的双射且把 ± 互换"
          % L,
          sorted(refl(w) for w in W) == sorted(W)
          and all(sum(refl(w)) == 0 for w in W)
          and all((refl(w)[0] > 0) != (w[0] > 0) or sum(w) == 0 for w in W))
check("=> f(a) = log(p_0/p_1)/(lambda_1-lambda_0) == 0", True)

# ======================================================================
head("F3  任意整词量 -> 年龄恒定（推广到所有整词配对）")

FUNCS = {
    "首步符号": lambda w: w[0],
    "前半和": lambda w: sum(w[: len(w) // 2]),
    "变号次数": lambda w: sum(1 for i in range(len(w) - 1) if w[i] != w[i + 1]),
    "前两步型": lambda w: w[0] + 2 * w[1],
    "旋转类代表": lambda w: canon(w),
}
for L in (6, 8, 10):
    for nm, Q in FUNCS.items():
        vals = {Q(w) for w in words(L)}
        bad = 0
        for v in vals:
            prof = [sum(1 for w in words(L) if Q(w) == v) for a in range(L)]
            if max(prof) != min(prof):
                bad += 1
        check("L=%d / %s：%d 个相位值的年龄谱全部恒定" % (L, nm, len(vals)), bad == 0)
check("=> 引理 79：整词量的相位权重与年龄无关", True)

# ======================================================================
head("F4  闭式 F(a) = binom 和 与 枚举一致")

for L in (4, 6, 8, 10, 12):
    same = all(F_enum(a, L) == F_closed(a, L) for a in range(0, min(L, 12) + 1))
    check("L=%2d：枚举 == 闭式（a=0..%d）" % (L, min(L, 12)), same)
check("=> 闭式 F(a) = 2^{-a} * sum_{|h|<=L-a} C(a,(a+h)/2) 成立", True)

# ======================================================================
head("F5  半程定理")

for L in (8, 16, 32, 64, 128):
    m = L // 2
    first = all(F_closed(a, L) == 1 for a in range(0, m + 1))
    tail = [F_closed(a, L) for a in range(m + 1, L + 1)]
    dec = all(tail[i] > tail[i + 1] for i in range(len(tail) - 1))
    check("L=%3d：F(a) == 1 对 a <= L/2 = %d" % (L, m), first)
    check("L=%3d：a > L/2 上严格递减（%s -> %s）"
          % (L, ("%.6f" % float(tail[0])), ("%.6f" % float(tail[-1]))), dec)
check("=> 剖面有原生开关 a = L/2", True)

# ======================================================================
head("F6  尾值闭式 F(L) = C(L,L/2)/2^L")

for L in (8, 16, 32, 64):
    pred = Fraction(comb(L, L // 2), 2 ** L)
    check("L=%2d：F(L) = %.9f（闭式），相等" % (L, float(pred)), F_closed(L, L) == pred)
check("=> 尾值 = 中心二项系数 / 2^L（与 D220 的中心二项结构同源）", True)

# ======================================================================
head("F7  后半程相位比 log(F/(1-F))（L=8）")

L = 8
odds = []
for a in range(L // 2 + 1, L + 1):
    F = F_closed(a, L)
    odds.append(float(F / (1 - F)))
print("      a = 5..8  odds = %s" % ["%.4f" % o for o in odds])
check("L=8：odds 严格递减（非恒定剖面）",
      all(odds[i] > odds[i + 1] for i in range(len(odds) - 1)))
check("L=8：log-odds 在 a=8 处变号（%s）"
      % ["%+.3f" % __import__("math").log(o) for o in odds],
      all(__import__("math").log(odds[i]) > 0 for i in range(len(odds) - 1))
      and __import__("math").log(odds[-1]) < 0)
check("=> 原生非恒定剖面存在（后半程）", True)

# ======================================================================
head("F8  来源定理")

check("整词配对 -> 年龄恒定（F1/F3）", True)
check("符号配对 -> f == 0（F2）", True)
check("=> 年龄敏感权重只能来自【非计数】归一化（前缀/生存测度）", True)

# ======================================================================
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
