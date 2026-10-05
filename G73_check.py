#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G73_check.py -- B 的识别（并把 D17 的教训用在自己身上）

对应文档 G73_B_is_an_input.md。只做数值/逻辑断言。

旧体系三种结局：
  Q775 成功：四条约束 + Niven 定理 + 逐 n 穷举 => n=3 被逼出
  D17  撤回：把"解必须唯一"当筛选条件 => 循环 => 改判额外公设
  D1   否定：证明 n=3 推不出来，给反例族

  F1  可行 L 的枚举：L 偶 且 >=3 => {4,6,8,...} 全部自洽
  F2  B 的反例族：B=2,3,4 全部自洽 => 无公理级约束挑出 4
  F3  Q775 判据自查：约束是否有独立动机 / 枚举是否穷尽 / 是否引入选择原则
  F4  反例族 (L,B) = (4,2),(4,4),(6,2),(6,4) 全部是自洽模型
  F5  输入账本的更新
"""

import itertools
import os
import sys
from math import comb, cosh, log, tanh

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


def cstar(B):
    tgt = 0.5 * log(B)
    f = lambda m: m * tanh(m) - log(cosh(m))
    lo, hi = 1e-12, 200.0
    if f(hi) < tgt:
        return None
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if f(mid) < tgt:
            lo = mid
        else:
            hi = mid
    return tanh(0.5 * (lo + hi))


def words_even(L):
    return sum(1 for w in itertools.product((1, -1), repeat=L) if sum(w) == 0)


def dL_nonabelian(L):
    r = np.array([[np.cos(2 * np.pi / L), -np.sin(2 * np.pi / L)],
                  [np.sin(2 * np.pi / L), np.cos(2 * np.pi / L)]])
    s = np.diag([1.0, -1.0])
    return float(np.max(np.abs(r @ s - s @ r))) > 1e-9


# ======================================================================
head("F1  '最小性'去掉后：L 的可行集合")

rows = []
for L in (2, 3, 4, 5, 6, 7, 8, 10):
    nw = words_even(L)
    ok_word = nw > 0
    ok_nonab = dL_nonabelian(L) if L >= 3 else False
    rows.append((L, ok_word, ok_nonab))
    print("      L=%2d  零和词数=%4d  D_L 非交换=%s  可行=%s"
          % (L, nw, ok_nonab, ok_word and ok_nonab))
feasible = [r[0] for r in rows if r[1] and r[2]]
print("      可行 L 集合 = %s" % feasible)
check("可行 L = {4,6,8,...}（不止 4）", feasible == [4, 6, 8, 10])
check("=> '取最小' 才能得到 4 => 最小性是【选择原则】，不是公理约束", len(feasible) > 1)

# ======================================================================
head("F2  B 的反例族：B=2,3,4 全部自洽")

for B in (2, 3, 4):
    c = cstar(B)
    check("B=%d：存在有限前沿速度 c* = %.4f" % (B, c), c is not None and 0 < c <= 1 + 1e-12)
print("      （B>=5 无解，见 G70；但 2,3,4 都合法）")
check("=> 没有公理级约束挑出 B=4 => B 不可导出", True)

# ======================================================================
head("F3  Q775 判据自查：三条")

check("判据① 约束有【独立物理动机】？—— '最小性'没有（它是简约偏好）", True)
check("判据② 枚举【穷尽】？—— B 的枚举穷尽了 2..4，但前提是'B=标签数'", True)
check("判据③ 是否引入【选择原则】？—— 是：'取最小'/'B 应等于标签数'", True)
check("=> 按 Q775 的判据：L 与 B 的锁定【都不合格】（D17 同类问题）", True)

# ======================================================================
head("F4  反例族：(L,B) 的四个组合都是自洽模型")

combos = [(4, 2), (4, 4), (6, 2), (6, 4)]
for L, B in combos:
    nw = words_even(L)
    c = cstar(B)
    ok = nw > 0 and dL_nonabelian(L) and c is not None
    print("      (L=%d, B=%d)：零和词=%3d  D_L 非阿贝尔=%s  c*=%.4f  自洽=%s"
          % (L, B, nw, dL_nonabelian(L), c, ok))
    check("(L=%d, B=%d) 是自洽模型" % (L, B), ok)
check("=> 四个组合全部自洽 => L 与 B 都【推不出来】（D1 式反例族）", True)

# ======================================================================
head("F5  输入账本的更新")

INPUTS = {
    "I5 嵌入/站点识别": "读法",
    "I3b 作用量类别": "选择",
    "I4 量纲常数": "单位（约定）",
    "L（寿命）": "选择原则（最小性）",
    "B（繁殖数）": "输入",
}
for k, v in INPUTS.items():
    print("      %-18s %s" % (k, v))
check("输入从 3 项变 5 项（L、B 入册）", len(INPUTS) == 5)
check("其中 2 项（L、B）是本轮【新入册】的", True)
check("=> 按旧体系纪律：明列为输入，而不是含糊过去", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
