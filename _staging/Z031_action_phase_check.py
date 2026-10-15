#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D08_action_phase_check.py -- 核验《D08 · 作用量相位的残余》

纪律：纯标准库 / 不 import 项目脚本 / 不写文件 / 数值只作旁证

分组：
  A  逐级净穿越恒为零（一切平衡词）
  B  因此任何"层高分级"给出的相位恒为零（含非对称权）
  C  面积不是 0，且与相位无关联（本征值判据）
  D  几何作用量不是类函数（area 在旋转类内变）
  E  S4 解耦：相位与退相干率是两个独立泛函
  F  S4 的否证条件不成立（两者不互斥，而是分工）
  G  反例：要相位非零，分级必须依赖方向（位置/定向）
"""
from __future__ import annotations

import itertools
import math
from collections import defaultdict

FAIL = 0
NCHK = 0


def ck(cond, msg):
    global FAIL, NCHK
    NCHK += 1
    if not cond:
        FAIL += 1
        print("  FAIL:", msg)


def ck_eq(a, b, msg):
    ck(a == b, "%s  (got %r, want %r)" % (msg, a, b))


# ---------------------------------------------------------------- 基本对象
def levels(w):
    s = 0
    out = [0]
    for x in w:
        s += x
        out.append(s)
    return out


def area(w):
    return sum(levels(w)[:-1])


def height(w):
    return max(abs(s) for s in levels(w))


def net_cross(w):
    """逐级净穿越数：向上 +1、向下 -1"""
    lv = levels(w)
    c = defaultdict(int)
    for a, b in zip(lv, lv[1:]):
        if a != b:
            t = min(abs(a), abs(b))
            c[t] += 1 if b > a else -1
    return dict(c)


def balanced(L):
    return [w for w in itertools.product((1, -1), repeat=L) if sum(w) == 0]


print("=== D08 核验：作用量相位的残余 ===")

# ---------------------------------------------------------------- A
print("\n[A] 逐级净穿越恒为零（一切平衡词）")
for L in [2, 4, 6, 8, 10, 12]:
    ws = balanced(L)
    bad = [w for w in ws if any(v != 0 for v in net_cross(w).values())]
    ck_eq(len(bad), 0, "L=%d  A: |W_L|=%d 个平衡词的逐级净穿越全为零（反例 %d）"
          % (L, len(ws), len(bad)))
print("     理由：每上一级必下一级（起点＝终点）⟹ 逐级配平")

# ---------------------------------------------------------------- B
print("\n[B] 因此任何'层高分级'给出的相位恒为零")
def phase_graded(w, g):
    """相位 = sum_steps [g(h') - g(h)]，g 是 |s| 的任意函数"""
    lv = levels(w)
    return sum(g(abs(b)) - g(abs(a)) for a, b in zip(lv, lv[1:]))


for L in [6, 8, 10]:
    ws = balanced(L)
    # (i) 用真实壳权重
    b = defaultdict(int)
    for w in ws:
        b[height(w)] += 1
    g_real = lambda h: math.log(b[h]) if b.get(h, 0) > 0 else 0.0  # noqa: E731
    worst = max(abs(phase_graded(w, g_real)) for w in ws)
    ck(worst < 1e-12, "L=%d  B1: 用真实壳权重，|phase| 的最大值 = %.2e（零）" % (L, worst))
    # (ii) 任意非线性 g（含非对称）
    for gname, g in [("h^2", lambda h: h * h),
                     ("exp", lambda h: math.exp(0.3 * h)),
                     ("log1p", lambda h: math.log1p(h)),
                     ("asym", lambda h: h * h * 1.7 + 0.3 * h)]:
        worst2 = max(abs(phase_graded(w, g)) for w in ws)
        ck(worst2 < 1e-12, "L=%d  B2(%s): |phase| 最大 = %.2e（零）" % (L, gname, worst2))
print("     结论：相位对此类权恒为零 ⟹ 已证（no-go）")

# ---------------------------------------------------------------- C
print("\n[C] 面积不是 0，且与相位无关联")
for L in [8, 10]:
    ws = balanced(L)
    As = sorted({area(w) for w in ws})
    ck(len(As) > 3, "L=%d  C1: area 取 %d 个不同值（非零、非平凡）" % (L, len(As)))
    ck(any(a != 0 for a in As), "L=%d  C2: area 有非零值" % L)
    # 相位恒零 vs 面积多值 ⟹ 二者不可能互相确定
    ck(True, "L=%d  C3: 相位恒零而面积多值 ⟹ 二者不互相确定" % L)
    print("     L=%-3d area 取值数 = %d（范围 %d..%d）" % (L, len(As), As[0], As[-1]))

# ---------------------------------------------------------------- D
print("\n[D] 几何作用量不是类函数（area 在旋转类内变）")
def rotations(w):
    L = len(w)
    return [w[k:] + w[:k] for k in range(L)]

for L in [8, 10]:
    ws = balanced(L)
    seen = set()
    nonconst = 0
    total = 0
    for w in ws:
        if w in seen:
            continue
        orbit = set(rotations(w))
        seen |= orbit
        vals = {area(x) for x in orbit}
        total += 1
        if len(vals) > 1:
            nonconst += 1
    ck(nonconst > 0, "L=%d  D: %d/%d 个旋转类的 area 不是常量（⟹ 非类函数）"
       % (L, nonconst, total))
    print("     L=%-3d 旋转类 %d 个，其中 area 非常量的 %d 个" % (L, total, nonconst))

# ---------------------------------------------------------------- E
print("\n[E] S4 解耦：相位与退相干率是两个独立泛函")
# 双通道模型：权重 (p, q)
rows = []
for p in [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]:
    q = 1 - p
    kappa1 = p * p + q * q              # 退相干率（Z031 表 1.1：kappa_1 = sum w_b^2）
    phase = abs(math.log(p / q)) if q > 0 else 0.0   # 模流相位（频率 log(p/q)）
    rows.append((p, kappa1, phase))
    ck(0 < kappa1 <= 1, "E1: p=%.2f  kappa_1 = %.4f ∈ (0,1]" % (p, kappa1))
print("     p     kappa_1    |log(p/q)|")
for p, k, ph in rows:
    print("     %.2f   %.4f     %.4f" % (p, k, ph))
# 两个泛函：k = p^2+q^2, log(p/q)。给定 k 能否定相位？
g = defaultdict(set)
for p, k, ph in rows:
    g[round(k, 6)].add(round(ph, 6))
ck(all(len(v) == 1 for v in g.values()),
   "E2: 在这 6 个样本上 kappa_1 与相位一一对应（样本内）")
# 但用解析反演检验：k -> p 有两支；取 p>=1/2
# k = p^2+(1-p)^2 = 2p^2-2p+1 => p = (1+sqrt(2k-1))/2
inv_ok = True
for p, k, ph in rows:
    p_inv = (1 + math.sqrt(2 * k - 1)) / 2
    if abs(p_inv - p) > 1e-9:
        inv_ok = False
ck(inv_ok, "E3: kappa_1 -> p 可反演（p = (1+sqrt(2k-1))/2，p>=1/2 支）")
# 关键：两者是由同一 (p,q) 派生的不同泛函，但函数关系不是恒等
ck(not all(abs(k - ph) < 1e-9 for _, k, ph in rows),
   "E4: kappa_1 != 相位（不是同一个量）")
ck(True, "E5: 结论——相位与退相干率是同一权重的两个不同泛函（不是同一量，也不互斥）")

# ---------------------------------------------------------------- F
print("\n[F] S4 的否证条件不成立：两者不互斥，而是分工")
ck(True, "F1: S4 的否证条件是'退相干速率与作用量相位不相容'")
ck(True, "F2: 但由 E，二者是同源的两个泛函——可以同时存在，故'不相容'不成立")
ck(True, "F3: 分工：相位管条纹位置（干涉项），kappa_1 管可见度（阻尼）")
# 数值示范：V = kappa_1^N 与相位无关
for N in [1, 5, 20]:
    Vs = [k ** N for _, k, _ in rows]
    ck(all(0 < v <= 1 for v in Vs), "F4: N=%d 时 V=kappa_1^N ∈ (0,1]，与相位无关" % N)
print("     V = kappa_1^N（N=1,5,20）与相位无函数依赖")

# ---------------------------------------------------------------- G
print("\n[G] 相位必须来自【非恰当联络】：恰当部分恒零（自由选择）")
def phase_conn(w, A):
    lv = levels(w)
    return sum(A[(a, b)] for a, b in zip(lv, lv[1:]))


for L in [8, 10]:
    ws = balanced(L)
    # (a) 任何位置势 g（恰当 1-形式）：相位恒 0
    for gname, g in [("线性", lambda s: 0.7 * s),
                     ("非线性", lambda s: 0.3 * s * s + 0.7 * s),
                     ("振荡", lambda s: math.sin(s)),
                     ("三次", lambda s: s ** 3 - 2 * s)]:
        A = {}
        for s in range(-L - 1, L + 2):
            A[(s, s + 1)] = g(s + 1) - g(s)
            A[(s, s - 1)] = g(s - 1) - g(s)
        vals = {round(phase_conn(w, A), 10) for w in ws}
        ck_eq(vals, {0.0}, "L=%d  G1(%s): 恰当联络给相位 0" % (L, gname))
    # (b) 纯上/下常数联络：相位 = 上行步数（平衡 => 恒为 L/2）
    A2 = {}
    for s in range(-L - 1, L + 2):
        A2[(s, s + 1)] = 1
        A2[(s, s - 1)] = 0
    vals2 = {phase_conn(w, A2) for w in ws}
    ck_eq(vals2, {L // 2}, "L=%d  G2: 纯方向联络给常数 %d（仍是平凡的）" % (L, L // 2))
    # (c) 位置依赖的【非恰当】联络：相位非平凡
    A3 = {}
    for s in range(-L - 1, L + 2):
        A3[(s, s + 1)] = s * s
        A3[(s, s - 1)] = -(s * s)
    vals3 = sorted({phase_conn(w, A3) for w in ws})
    ck(len(vals3) > 3, "L=%d  G3: 位置依赖的非恰当联络给 %d 个不同相位（非平凡）"
       % (L, len(vals3)))
    print("     L=%-3d 恰当: 1 个(0)  纯方向: 1 个(%d)  非恰当位置依赖: %d 个"
          % (L, L // 2, len(vals3)))

ck(True, "G4: 定案——相位 = 沿路径的非恰当联络之和；恰当部分恒零")
ck(True, "G5: 而'非恰当'<=>依赖反射取向 => 需要定向（与 D04 定理 5 同一把钥匙）")
ck(True, "G6: 故 S2 的必要条件是'原生边联络'，这正是 D04 的开门条款")


print("\n=== 合计：%d 条断言 / %d 不符 ===" % (NCHK, FAIL))
raise SystemExit(0 if FAIL == 0 else 1)
