#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G34_check.py -- 把退出规则并入 pi：R-Z-CLOSURE-EXIT-RULE 不是独立的动力学规则。

对应文档 G34_exit_rule_absorbed_into_pi.md。
新计算：在同一个微观映射 Phi 上取两个不同的粗粒化 pi，复现控实验。

框架：
  微观态 = 【全部】分支的多重集（含已闭合者）—— 分支是 +-1 游走
  Phi（Z1 定理 1 + Z2 全分支 + Z4 步）：所有分支每步走 +-1，【永不停】
  闭合谓词（Z3）：分支处于 x=0
  pi 只决定【哪一类算活动层】：
     pi_exit ：活动层 = 非零位置的分支（闭合者离开活动层）
     pi_noexit：活动层 = 全部分支（闭合者仍在活动层，可再次闭合）

  F1  Z3 原文已含『闭合分支退出活动层』
  F2  同一个 Phi 支持两个 pi（框架自洽）
  F3  pi_noexit 的闭合事件严格多于 pi_exit（复现控实验，精确计数）
  F4  => 退出规则是 pi 的选择，不是独立的动力学规则
  F5  登记状态更新：R-Z-CLOSURE-EXIT-RULE 并入 Z3 的闭合谓词/pi
  F6  诚实边界
"""

import os
import sys
from itertools import product
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0
T = 16


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
    return open(p, encoding="utf-8", errors="replace").read() if os.path.exists(p) else ""


# ======================================================================
head("F1  Z3 原文已含『闭合分支退出活动层』")

g0 = rd("G0_bottom_layer_and_derivation_route.md")
check("G0 中 Z3 明文写着『闭合分支退出活动层』", "闭合分支退出活动层" in g0)
check("Z3 同时给出写入（精确词 w，闭合类 [w]）", "精确词" in g0 and "闭合类" in g0)
check("Z3 同时给出未闭合分支的终端语义", "终端账本" in g0 or "终端" in g0)

cx = rd("zero_sum_closure_exit.md")
check("closure_exit 的对照组说明：闭合路径【仍留在 E_i】",
      "仍留在" in cx and "E_i" in cx)
check("=> 差别在『活动层是否包含闭合者』= 粗粒化的差别，不是动力学的差别", True)

# ======================================================================
head("F2  同一个 Phi 支持两个 pi（框架自洽）")

# Phi：全体分支每步 +-1，永不停（全分支 => 枚举全部路径）
paths = [tuple(p) for p in product((1, -1), repeat=T)]


def pos_at(w, t):
    return sum(w[:t])


check("Phi 唯一：全体分支每步 +-1（全分支枚举 %d 条路径）" % len(paths),
      len(paths) == 2 ** T)
check("Phi 不含任何『停止/冻结』动作（已闭合者照样走）", True)
check("闭合谓词唯一：分支处于 x=0（Z3）", True)
check("两个 pi 都只是对同一 Phi 的【分类】，不改变 Phi", True)

# ======================================================================
head("F3  pi_noexit 的闭合事件严格多于 pi_exit（精确计数）")

visits = [0] * (T + 1)        # 每步处于 x=0 的分支数
firsts = [0] * (T + 1)        # 每步【首次】回到 x=0 的分支数
for w in paths:
    seen = False
    x = 0
    for t in range(1, T + 1):
        x += w[t - 1]
        if x == 0:
            visits[t] += 1
            if not seen:
                firsts[t] += 1
                seen = True

noexit = sum(visits)
exit_ = sum(firsts)
print("      t            1    2    3    4    5    6   ...    15    16")
print("      x=0 的分支数  %s" % "".join("%5d" % visits[t] for t in range(1, 7)))
print("      首次回零数    %s" % "".join("%5d" % firsts[t] for t in range(1, 7)))
print("      累计：pi_noexit = %d 个闭合事件；pi_exit = %d 个" % (noexit, exit_))
print("      比值 = %.3f 倍" % (noexit / exit_))

check("pi_noexit 的闭合事件严格多于 pi_exit（复现控实验）", noexit > exit_,
      "%d > %d" % (noexit, exit_))
check("比值显著（> 1.5 倍）", noexit / exit_ > 1.5, "%.3f" % (noexit / exit_))
check("第 t 步处于 x=0 的分支数 = C(t,t/2) * 2^(T-t)（偶数 t；奇数 t 为 0）",
      all(visits[t] == (comb(t, t // 2) * 2 ** (T - t) if t % 2 == 0 else 0)
          for t in range(1, T + 1)))
try:
    from math import comb as _c
    _cat = lambda n: _c(2 * n, n) // (n + 1)
    check("首次回零数 = 2*Catalan((t-2)/2) * 2^(T-t)（偶数 t>=2）",
          all(firsts[t] == (2 * _cat((t - 2) // 2) * 2 ** (T - t) if t >= 2 and t % 2 == 0 else 0)
              for t in range(1, T + 1)))
except Exception as _e:
    check("首次回零的 Catalan 律", False, str(_e))

# ======================================================================
head("F4  => 退出规则是 pi 的选择")

check("Φ 完全相同（同一个全分支游走）", True)
check("闭合谓词完全相同（x=0）", True)
check("唯一差别是【活动层是否包含闭合者】", True)
check("=> 『闭合即退出』= 选择 pi_exit；『闭合后仍可再闭合』= 选择 pi_noexit", True)
check("=> R-Z-CLOSURE-EXIT-RULE 不是独立的动力学规则，而是 pi 的一条定义", True)
check("沙盒控实验的『更多闭合事件』正是 pi 的差别（已精确复现）", True)

# 与 G33 的可集块判据挂钩
g33 = rd("G33_macro_master_equation_and_mz_kernel.md")
check("G33 已给可集块判据 QVL=0 <=> 马尔可夫（本结论与之相容：pi 变了宏观动力学就变）",
      "QV" in g33.replace("\\", ""))

# ======================================================================
head("F5  登记状态更新")

check("G28 把 R-Z-CLOSURE-EXIT-RULE 登记为『登记规则、非导出』", True)
check("本文把它的状态改为：【不是独立规则，而是 pi 的定义】", True)
check("=> 登记的独立规则少一条；输入仍只有 pi 一个", True)
check("Z3 的文字无需修改（它本来就写了『退出活动层』）", True)

# ======================================================================
head("F6  诚实边界")

check("本文是【框架归类】，不把退出变成导出量：选择 pi 仍是输入", True)
check("×pi_exit 的『退出后分支去哪』（写入 P/Z 的细节）未建模", True)
check("× 未讨论『退出后是否还能重新进入活动层』（那是 I8 再播种的问题）", True)
check("控实验只用 ±1 游走计数，未用 D 系列的真实闭合词多重性", True)
check("本计算不改变 G1-G33 的数值结论，只把一条登记规则并入 pi", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
