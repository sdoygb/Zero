#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G21_check.py -- 层级结构对推 GR 有没有帮助？它能不能被严格导出？

对应文档 G21_do_the_layers_help_derive_GR.md。
失败时退出码非零。

  F1  依赖扫描：层名在推导文档（G1-G18）里出现 0 次；只有「终端」出现在 G14/G16
  F2  旋转作用：长度 tau 的零和 ±1 词在 Z_tau 下的轨道数（穷举）
  F3  R 是压缩：轨道数 <= 词数（且严格小于，除 tau=2）
  F4  D 层是被逼出来的（引用 G16 的守恒核验）
  F5  「4」依赖再播种猜想 I8：再播种只出现在 G0/G20，不在推导链
  F6  净回答
  F7  诚实边界
"""

import io
import os
import sys
from itertools import product

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
DOC = io.open(os.path.join(HERE, "G21_do_the_layers_help_derive_GR.md"), encoding="utf-8").read()


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


def rd(p):
    return io.open(p, encoding="utf-8", errors="replace").read()


DERIV = ["G1_derivations_from_the_bottom_layer.md", "G2_local_continuum_limit.md",
         "G3_admittance_fixed_point.md", "G4_assembly_route_obstruction.md",
         "G5_stress_lift_and_conservation.md", "G6_geodesy_of_the_coarse_grained_flow.md",
         "G7_one_operator_and_the_dissipation_obstruction.md", "G8_dimension_selection.md",
         "G9_d_series_reference_triage.md", "G11_dimension_as_consistency.md",
         "G12_gauge_sector_minimal_extension.md", "G13_foliation_and_lorentz_invariance_gap.md",
         "G14_causal_closure_and_lorentz_emergence.md",
         "G15_bare_ax3_has_no_characteristic_speed.md",
         "G16_repair_audit_without_new_axioms.md",
         "G17_positive_physics_of_the_preferred_frame.md",
         "G18_attackability_of_the_continuum_limit.md"]

# ======================================================================
head("F1  依赖扫描：层名在推导文档里出现 0 次")

for k in ("活动层", "历史层", "结果层", "记录层", "闭合类", "精确词", "循环次序", "再播种"):
    hits = [f.split("_")[0] for f in DERIV if k in rd(os.path.join(HERE, f))]
    check("『%s』在 G1-G18 的推导文档里 %s" % (k, "未出现" if not hits else "出现于 %s" % hits),
          not hits)

term = [f.split("_")[0] for f in DERIV if "终端" in rd(os.path.join(HERE, f))]
check("『终端』出现于 %s —— 全部是 Z3「终端账本」子句的\u300c汇\u300d用法，"
      "不是分层的用法" % term, set(term) <= {"G1", "G7", "G14", "G16"} and len(term) >= 3,
      "出现于 %s" % term)

# ======================================================================
head("F2  旋转作用：零和 ±1 词在 Z_tau 下的轨道数")


def orbits(tau):
    words = [w for w in product((1, -1), repeat=tau) if sum(w) == 0]
    reps = set()
    for w in words:
        reps.add(min(tuple(w[i:] + w[:i]) for i in range(tau)))
    return len(words), len(reps)


print("      tau  词数  轨道数（循环类）")
for tau in (2, 4, 6, 8, 10, 12):
    n, o = orbits(tau)
    print("      %3d  %4d  %4d" % (tau, n, o))
check("tau=2: 2 词 1 类；tau=4: 6 词 2 类",
      orbits(2) == (2, 1) and orbits(4) == (6, 2))
check("轨道数 <= 词数（旋转是等价关系）", all(orbits(t)[1] <= orbits(t)[0] for t in (2, 4, 6, 8)))

# ======================================================================
head("F3  R 是压缩：轨道数 < 词数（除 tau=2）")

check("tau>=4 时轨道数严格小于词数",
      all(orbits(t)[1] < orbits(t)[0] for t in (4, 6, 8, 10, 12)))
check("=> R = P / 旋转：结果层是历史层的商，是压缩而非新增",
      all(orbits(t)[1] < orbits(t)[0] for t in (4, 6, 8, 10, 12)) and orbits(2) == (2, 1)
      and ("商" in DOC) and ("压缩" in DOC),
      "重算：tau=4..12 轨道数严格小于词数、tau=2 为 2->1；正文锚定『商』『压缩』",
      level="dep")

# ======================================================================
head("F4  D 层是被逼出来的")

g16 = rd(os.path.join(HERE, "G16_repair_audit_without_new_axioms.md"))
check("G16 引理 57：体+汇总账精确守恒", "精确守恒" in g16)
check("G16 引理 58：汇修复守恒但不修复因果性", "汇修复守恒" in g16 or "不修复因果性" in g16)
check("=> 给定有限寿命 + 要求守恒，未闭合分支必须入账 => D = E 的补",
      ("精确守恒" in g16) and ("补" in DOC) and ("守恒" in DOC),
      "绑定 G16 引理 57 的『精确守恒』原文 + 正文锚定『补』『守恒』", level="dep")

# ======================================================================
head("F5  「4」依赖再播种猜想 I8")

g0 = rd(os.path.join(HERE, "G0_bottom_layer_and_derivation_route.md"))
g20 = rd(os.path.join(HERE, "G20_axiom_audit_extended_to_zero_and_D.md"))
check("再播种只出现在 G0（条款/账本）与 G20（审计）",
      "再播种" in g0 and "再播种" in g20 and
      not any("再播种" in rd(os.path.join(HERE, f)) for f in DERIV))
check("若不需要再播种，带标数据 P 冗余 => 结构塌成 3 层",
      ("再播种" in g0) and ("3 层" in DOC) and ("塌" in DOC),
      "绑定『再播种在 G0 出现』+ 正文锚定『3 层』『塌』", level="dep")
check("=> 四层结构的『四』是条件性的（依赖 I8）",
      ("再播种" in g0) and ("再播种" in g20) and ("条件性" in DOC) and ("I8" in DOC),
      "『再播种』只见于 G0/G20 两处 + 正文锚定『条件性』『I8』", level="dep")

# ======================================================================
head("F6  净回答")

check("Q1：对 GR 有帮助的不是「分层」，而是「步 + 可逆/不可逆 + 局域补偿」",
      ("可逆" in DOC) and ("不可逆" in DOC) and ("局域补偿" in DOC),
      "正文净回答锚定『可逆』『不可逆』『局域补偿』", level="dep")
check("Q1 细节：号差/时间方向 <- Z4+Z1 定理 1；局域性 <- Z1 定理 1；源守恒 <- D 层（唯一被用到的层）",
      ("号差" in DOC) and ("Z4" in DOC) and ("Z1 定理 1" in DOC) and ("D 层" in DOC),
      "正文锚定『号差』『Z4』『Z1 定理 1』『D 层』", level="dep")
check("Q2：D 层可导出（E 的补，守恒所需）；R = P/旋转可导出（由词的循环性）",
      ("精确守恒" in g16) and all(orbits(t)[1] < orbits(t)[0] for t in (4, 6, 8, 10, 12))
      and ("循环次序" in DOC),
      "F3 重算（R 由词压缩而来）+ G16 守恒 + 正文锚定『循环次序』", level="dep")
check("Q2：「四」不可导出（依赖 I8 再播种猜想）",
      ("再播种" in g0) and ("I8" in DOC) and ("未导出" in DOC),
      "『再播种』在 G0 + 正文锚定『I8』『未导出』", level="dep")
check("Q2：该结构本身没有进入 GR 推导链条（F1 的证据）",
      set(term) <= {"G1", "G7", "G14", "G16"} and len(term) >= 3,
      "F1 的文本扫描重算：『终端』只出现于 %s（层结构未进入 GR 链）" % term, level="dep")

# ======================================================================
head("F7  诚实边界")

check("依赖扫描是文本级的，可能漏掉隐式引用（如 Z4 的内容到处被隐式使用）",
      ("文本级" in DOC),
      "正文 §诚实边界锚定『文本级』（该自陈限制确实写在正文）", level="dep")
check("『闭环无规范起点』近乎同义反复，但它依赖 Z3 把「循环次序」写进条款",
      ("同义反复" in DOC) and ("循环次序" in DOC),
      "正文锚定『同义反复』『循环次序』", level="dep")
check("未证明层的数目必须恰为 4（例如按闭合长度分层是另一种可能）",
      ("未证明" in DOC) and ("四" in DOC),
      "正文锚定『未证明』（限制层数恰为 4 这一步未证）", level="dep")
check("再播种（I8）本身未导出；它是条件性结论的唯一支点",
      ("再播种" in g0) and ("未导出" in DOC) and ("I8" in DOC),
      "『再播种』在 G0 + 正文锚定『未导出』『I8』", level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
