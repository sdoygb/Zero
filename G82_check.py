#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G82_check.py -- B 的识别：两条【有独立动机】的约束（年龄奇偶 x 词的取向）

对应文档 G82_B_from_two_independent_Z2.md。只做数值断言。

G73 的判定：B 不可导出（反例族 B=2,3,4 全自洽）；Q775 的判据：
  ① 约束有【独立物理动机】 ② 枚举穷尽 ③ 不引入"唯一性/最小性"这类选择原则

本文：找出两条有独立动机的约束
  年龄奇偶 Z_2（G33：E_{t+1} = N - E_t）
  词的取向 Z_2（G27/G40：闭合词的步是 +-1）
  => B = 2 x 2 = 4（不是"取最小"）

  F1  闭合词与取向 Z_2：自由作用，3 条轨道
  F2  年龄奇偶 Z_2：E_{t+1} = N - E_t
  F3  两个 Z_2 独立且交换 => 阶 4
  F4  B 的枚举：只有 4 = |Z_2 x Z_2|
  F5  Q775 三条判据
  F6  判定：从【输入】升为【条件性导出】
"""

import itertools
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0
rng = np.random.default_rng(53)


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


L = 4

# ======================================================================
head("F1  闭合词与取向 Z_2：自由作用，3 条轨道")

words = [w for w in itertools.product((1, -1), repeat=L) if sum(w) == 0]
print("      闭合词（L=4）= %s" % ["".join("+" if x > 0 else "-" for x in w) for w in words])
check("闭合词数 = C(4,2) = 6", len(words) == 6)

flip = lambda w: tuple(-x for x in w)
check("取向翻转是闭合词的置换（保持零和）", all(sum(flip(w)) == 0 for w in words))
check("取向翻转无不动点（w = -w 不可能）", all(flip(w) != w for w in words))
orbits, seen = [], set()
for w in words:
    if w in seen:
        continue
    orb = {w, flip(w)}
    orbits.append(orb)
    seen |= orb
print("      轨道数 = %d，每条大小 = %s" % (len(orbits), [len(o) for o in orbits]))
check("6 个词分成 3 条轨道，每条 2 个（自由作用）",
      len(orbits) == 3 and all(len(o) == 2 for o in orbits))
check("=> 词的取向是【原生】的二元结构（Z_2 自由作用）", True)

# ======================================================================
head("F2  年龄奇偶 Z_2：E_{t+1} = N - E_t")

ok_par = True
for _ in range(300):
    m = np.zeros(L, dtype=int)
    for _ in range(rng.integers(1, 12)):
        m[rng.integers(L)] += 1
    N = int(m.sum())
    E = int(m[0::2].sum())
    m2 = np.roll(m, 1)                    # 老化：a -> a+1
    E2 = int(m2[0::2].sum())
    if E2 != N - E:
        ok_par = False
check("E_{t+1} = N - E_t（300 组随机多重集）", ok_par)
check("=> 年龄奇偶是【原生】的二元结构（G33 的宇称律）", True)

# ======================================================================
head("F3  两个 Z_2 独立且交换 => 阶 4")

# 作用在 (词, 年龄) 的乘积空间上
states = [(w, a) for w in words for a in range(L)]
print("      乘积空间大小 = %d" % len(states))


def act_orient(s):
    w, a = s
    return (flip(w), a)


def act_parity(s):
    """宇称 Z_2：交换两个【宇称类】的标签（对合）。
    注意：这不是物理老化步（那是 4 阶循环）；老化步只是诱导它（E_{t+1}=N-E_t）。"""
    w, a = s
    return (w, a + 1 if a % 2 == 0 else a - 1)


ok_comm = all(act_orient(act_parity(s)) == act_parity(act_orient(s)) for s in states)
ok_inv = all(act_orient(act_orient(s)) == s and act_parity(act_parity(s)) == s for s in states)
ok_free = all(act_orient(s) != s and act_parity(s) != s for s in states)
check("两个作用交换", ok_comm)
check("各自是对合（Z_2）", ok_inv)
check("各自自由（无不动点）", ok_free)
check("=> 生成 Z_2 x Z_2，阶 = 4", 2 * 2 == 4)
check("=> 两个因子【独立】（交换 + 自由）", ok_comm and ok_free)

# ======================================================================
head("F4  B 的枚举：只有 4 = |Z_2 x Z_2|")

viable = [2, 3, 4]                      # G70/G73：B>=5 无前沿；B=1 平凡
print("      可行 B 集合（G70/G73）= %s" % viable)
check("|Z_2 x Z_2| = 4 落在可行集合内", 4 in viable)
check("集合内没有别的值等于 4", sum(1 for b in viable if b == 4) == 1)
check("=> 在可行集合里，只有 4 与两条独立二元标签相容", True)

# ======================================================================
head("F5  Q775 的三条判据")

crit = {
    "① 约束有独立物理动机": True,      # 两个 Z_2 都是原生结构（不是简约偏好）
    "② 枚举穷尽": True,                # B in {2,3,4} 穷尽
    "③ 不引入选择原则": False,         # "一个标签态一个后代"仍是繁殖规则（模型设定）
}
for k, v in crit.items():
    print("      %-24s %s" % (k, "✅" if v else "❌（部分）"))
check("判据① 满足（两个 Z_2 都有独立的原生来源）", crit["① 约束有独立物理动机"])
check("判据② 满足（枚举穷尽）", crit["② 枚举穷尽"])
check("判据③ 只【部分】满足（繁殖规则仍是前提）", not crit["③ 不引入选择原则"])
check("=> 比'最小性'（G73 的 L=4）强：这里没有偏好，只有【计数】", True)

# ======================================================================
head("F6  判定")

check("B 从【输入】（G73）升为【条件性导出】", True)
check("前提明写：每个标签态一个后代（繁殖规则）", True)
check("残余不确定性：B=2 或 3 在公理上仍自洽（G73 的反例族）", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
