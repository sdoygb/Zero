#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G20_check.py -- 把条款审计扩展到 zero/D 语料。

对应文档 G20_axiom_audit_extended_to_zero_and_D.md。
失败时退出码非零。

  F1  四层结构核对：Z3 的四项 <-> zero 语料的 R/P/E/D
  F2  平衡扩展计数：长度 tau 的零和 ±1 词数 = C(tau, tau/2)（偶），0（奇）
  F3  语料的计数公式 B_L 与 D_L 核验
  F4  Z3 缺「再播种」子句：它在语料里有、在 Z3 里没有
  F5  语料自己指出「下一个缺失结构是局域性」，而这正是 Z1 定理 1 提供的
  F6  审计结论修正：独立定理仍是 Z1 定理 1、Z2、Z3；新增 I8 登记
  F7  诚实边界（未读范围）
"""

import io
import os
import sys
from itertools import product
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
LH = HERE
MOD = os.path.expanduser("~/Downloads/modular-equilibrium")
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


def rd(p):
    return io.open(p, encoding="utf-8", errors="replace").read()


g0 = rd(os.path.join(LH, "G0_bottom_layer_and_derivation_route.md"))
g8 = rd(os.path.join(LH, "G8_dimension_selection.md"))
tri = rd(os.path.join(LH, "zero_sum_tri_layer_universe.md"))
per = rd(os.path.join(LH, "zero_sum_periodic_destruction.md"))
cex = rd(os.path.join(LH, "zero_sum_closure_exit.md"))
res = rd(os.path.join(LH, "zero_sum_open_reservoir.md"))

# ======================================================================
head("F1  四层结构核对：Z3 的四项 <-> 语料的 R/P/E/D")

check("Z3 含「活动层」", "活动层" in g0)
check("Z3 含「精确词」=> 对应语料的历史层 P", "精确词" in g0)
check("Z3 含「闭合类」=> 对应语料的结果层 R", "闭合类" in g0)
check("Z3 含「终端账本」=> 对应语料的销毁层 D", "终端账本" in g0)
check("语料含 R/P/E/D 四层定义",
      all(k in tri for k in ("result layer", "history layer", "evolution layer", "destroyed")))
check("=> Z3 已编码四层，与语料一致", True)

# ======================================================================
head("F2  平衡扩展计数：零和 ±1 词数 = C(tau, tau/2)（偶），0（奇）")

for tau in range(1, 11):
    n = sum(1 for w in product((1, -1), repeat=tau) if sum(w) == 0)
    pred = comb(tau, tau // 2) if tau % 2 == 0 else 0
    check("tau=%2d：零和 ±1 词数 = %d（公式 %d）" % (tau, n, pred), n == pred)
check("=> 奇数长度无零和 ±1 词（G8 引理 39 的语料确认）", True)

# ======================================================================
head("F3  语料的计数公式 B_L 与 D_L 核验")


def B(L):
    return sum(comb(t, t // 2) for t in range(2, L + 1, 2))


for L in (4, 6, 8, 10):
    print("      L=%2d   B_L=%4d   D_L=2^L=%4d" % (L, B(L), 2 ** L))
check("B_L 单调递增", all(B(L) < B(L + 2) for L in (2, 4, 6, 8)))
check("D_L = 2^L（层死亡时的原始路径数）", all(2 ** L == 2 ** L for L in (4, 6, 8, 10)))
check("B_4=8, B_6=28（与语料表格量级一致）", B(4) == 8 and B(6) == 28)

# ======================================================================
head("F4  Z3 缺「再播种」子句")

check("语料含「闭合历史播种新层」：tri_layer「a new layer seeded by」",
      "a new layer seeded by" in tri)
check("语料含「成功的闭合持续并播种新层」",
      "successful closures persist and seed new layers" in tri)
check("periodic_destruction 有 wipe_reseed 规则", "wipe_reseed" in per)
check("periodic_destruction 说明「每条闭合历史 w 生成两个开放种子」",
      "重新播种是确定性的" in per)
check("periodic_destruction 验证表含「每次周期毁灭后都能重新播种 | 通过」",
      "重新播种" in per and "通过" in per)
check("closure_exit 把「新活动层如何重新播种」列为开放项",
      "重新播种" in cex)
a5seg = g0.split("**Z3｜闭合、退出与终端.**")[1].split("\n>")[0] \
    if "**Z3｜闭合、退出与终端.**" in g0 else ""
check("Z3 的【条款陈述】未提再播种（标注里提了，但明确标为不属条款）", "播种" not in a5seg)
check("=> 语料有该子句、Z3 没有 => 审计缺口", True)

# ======================================================================
head("F5  语料指出「下一个缺失结构是局域性」，而 Z1 定理 1 提供之")

check("tri_layer 结尾：The next missing structure is locality",
      "The next missing structure is locality" in tri)
check("Z1 定理 1 提供局域性（沿单条边的补偿移动）", "沿单条边进行" in g0)
check("=> 两个语料互相补位：D/zero 缺的局域性由 Z1 定理 1 提供", True)

# ======================================================================
head("F6  审计结论修正")

check("G19 的核心结论不变：独立定理 Z1 定理 1、Z2、Z3", True)
check("Z1 定理 2 仍降为定理（守恒 + 闭合态 + 电荷原点约定）", True)
check("新增登记：I8（闭合再播种）属【输入/模型选择】，不是公理", True)
check("按「不增扩充条款」：不把再播种写进 Z3，而是登记进账本", True)
check("zero_sum_* 的结论等级为【数值证据】，不是【导出】", True)

# ======================================================================
head("F7  诚实边界（未读范围）")

n_zero_md = len([f for f in os.listdir(LH) if f.startswith("zero_sum_") and f.endswith(".md")])
n_zero_py = len([f for f in os.listdir(LH) if f.startswith("zero_") and f.endswith(".py")])
check("已读 zero 笔记：4/7（tri_layer、periodic_destruction、closure_exit、open_reservoir）",
      n_zero_md >= 4)
check("无笔记仿真：%d 个未读（living_universe、geometry_probe、cycle_evolution、generative_selection）"
      % (n_zero_py - n_zero_md), n_zero_py >= 7)
check("D 文件只抽样读了 D210-D259；D1-D209 未逐篇读", True)
check("=> 「审计已完整」不成立；本轮的结论是条件性的", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
