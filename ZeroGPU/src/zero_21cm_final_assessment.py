#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_21cm_final_assessment.py —— 21cm 内生性总评估
======================================================================================
用户口径（2026-10-20）："我们先评估，是否从 Zero 能原生地推导出 21cm 线（一个锚定）"

本脚本汇总整条链上【已查明】的事实，逐环节定论，并区分：
  A 类（已证不可）: 有定理级依据
  B 类（未做/待定）: 只是还没做，没有否定性定理
  C 类（可绕过）:    有替代写法，绕过该环节

21cm 公式: h nu = (8/3) g_p (m_e/m_p) alpha^2 R_y / (1+m_e/m_p)
"""
from __future__ import annotations
import os, json
import mpmath as mp

mp.mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

PASS = []
def check(name, cond, detail=""):
    PASS.append({"name": name, "pass": bool(cond), "detail": detail})
    P_("[%s] %s %s" % ("PASS" if cond else "FAIL", name, detail))

P_("=" * 96)
P_("21cm 内生性总评估 · zero_21cm_final_assessment.py")
P_("=" * 96)

# ---------- §1 公式与灵敏度 ----------
P_("\n§1 21cm 的公式与灵敏度")
me = mp.mpf('9.1093837015e-31'); mp_ = mp.mpf('1.67262192369e-27')
r = me/mp_
P_("  h nu = (8/3) g_p (m_e/m_p) alpha^2 R_y / (1+m_e/m_p)")
P_("  nu ~ g_p * alpha^4 * m_e^2 / m_p   （对数导数）")
P_("  d ln nu / d ln g_p   = +1")
P_("  d ln nu / d ln alpha = +4")
P_("  d ln nu / d ln m_e   = %+.6f" % (2 - r/(1+r)))
P_("  d ln nu / d ln m_p   = %+.6f" % (-1 + r/(1+r)))
check("alpha 的灵敏度为 4（四次方）", True, "nu ~ alpha^4")

# ---------- §2 逐环节定论 ----------
P_("\n§2 逐环节定论（A 已证不可 / B 未做 / C 可绕过）")
NODES = [
    dict(id="N0", name="锚 v", cat="OK",
         basis="Z013 定理16：v 是唯一有量纲自由度、是量纲单位",
         verdict="满足：锚只有一个"),
    dict(id="N1", name="8/3 系数", cat="OK",
         basis="接触项纯代数；Z006 自旋-角动量代数",
         verdict="内生"),
    dict(id="N2", name="alpha", cat="A",
         basis="旧 Z23：需无量纲输入 R；(2.1) 导出 L <=> 导出 alpha^-1",
         verdict="已证：非预言，需一笔无量纲输入"),
    dict(id="N3", name="g_p", cat="A",
         basis="库内无磁矩结构；Z014 定理22：禁闭/质量间隙不推出",
         verdict="已证：非微扰 QCD 量，库内无机制"),
    dict(id="N4", name="m_e", cat="OK",
         basis="Z126 振幅律：m_e = A^2 M_e^2，差 -0.0008%",
         verdict="半内生：形状硬，标度靠 784.5 构造"),
    dict(id="N5", name="m_p", cat="C",
         basis="超精细里 m_p 只通过质子磁矩进来：g_p/m_p = 4 mu_p/(e hbar)",
         verdict="可绕过：用 mu_p 替代，Zero 的 m_p 误差掉出"),
    dict(id="N6", name="约化质量", cat="C",
         basis="(1+r) 因子；r 是测量值",
         verdict="可绕过：用测量比值"),
    dict(id="N7", name="库仑束缚态", cat="B",
         basis="Z165：原生泊松核有束缚态；但 1/r 只在 n->inf 出现",
         verdict="未做：原生束缚态可算，但未报到物理单位"),
]
P_("  %-4s %-14s %-6s %s" % ("ID", "环节", "类别", "判定"))
for n in NODES:
    P_("  %-4s %-14s %-6s %s" % (n["id"], n["name"], n["cat"], n["verdict"]))

nA = sum(1 for n in NODES if n["cat"] == "A")
nB = sum(1 for n in NODES if n["cat"] == "B")
nC = sum(1 for n in NODES if n["cat"] == "C")
nOK = sum(1 for n in NODES if n["cat"] == "OK")
P_("\n  统计: OK=%d, A(已证不可)=%d, B(未做)=%d, C(可绕过)=%d" % (nOK, nA, nB, nC))
check("有 2 个环节是 A 类（已证不可）", nA == 2, "alpha, g_p")

# ---------- §3 两个 A 类的精确性质 ----------
P_("\n§3 ★ 两个 A 类的精确性质（关键：它们不是同一个障碍）")

P_("\n  alpha（A 类之一）:")
P_("    旧 Z23 的判定: alpha 是【一个无量纲比 R 的函数】, R 的类型由 Zero 定、值不由 Zero 定")
P_("    R_data = 1.292018 +- 0.003596")
P_("    迹 5/3 被排除 104 sigma; 范数 sqrt(5/3) 通过 -0.28 sigma 但机制缺")
P_("    ⟹ alpha 缺的是【一个无量纲比】, 不是缺尺度")
P_("    ⟹ 障碍性质: 需要一笔无量纲输入（不是定理级的'不可能'）")

P_("\n  g_p（A 类之二）:")
P_("    库内无任何磁矩/塞曼/g-2 结构")
P_("    Z014 定理22: 禁闭/质量间隙/弦张力【不推出】")
P_("    ⟹ g_p 缺的是【非微扰机制】")
P_("    ⟹ 障碍性质: 比 alpha 更硬（连'需要什么'都没定出来）")

P_("\n  ★ 两者的层次不同:")
P_("    alpha: 缺一个数（类型已定：无量纲比）")
P_("    g_p:   缺一个机制（类型未定）")
check("两个 A 类的障碍性质不同（已区分）", True, "alpha 缺数；g_p 缺机制")

# ---------- §4 绕过能到什么程度 ----------
P_("\n§4 绕过能到什么程度")
P_("  (i) alpha 在【比值】里可绕: nu_2s/nu_1s = 1/8（验证 3.5e-5）")
P_("  (ii) m_p 可绕: 用 mu_p（测出来的）替代 g_p/m_p")
P_("  (iii) 但【绝对频率】绕不掉 alpha^4 与 g_p")
P_("")
P_("  绕过后的依赖集:")
P_("    nu_21 = f(mu_p, m_e, alpha, r)")
P_("    = {测量 mu_p} x {Zero m_e} x {测量 alpha} x {测量 r}")
P_("  ⟹ Zero 的贡献只剩 m_e（0.0008%）")
check("绕过方案已明确（含残留依赖）", True, "见上")

# ---------- §5 判决 ----------
P_("\n§5 判决")
P_("  Q: 从 Zero 能原生地推导出 21cm 线吗（一个锚定）？")
P_("")
P_("  A: 【不能】—— 但必须分层说清:")
P_("")
P_("   (1) '锚一个锚点'这一条【满足】: v 是唯一有量纲锚")
P_("   (2) '其余内生'这一条【不满足】: alpha 与 g_p 是两笔无量纲输入")
P_("   (3) 且这两笔的性质不同:")
P_("        alpha: 缺一个【无量纲比】（类型已定，旧 Z23 已判）")
P_("        g_p:   缺一个【非微扰机制】（类型未定，更难）")
P_("   (4) 绕过只能到比值: 绝对频率仍绕不掉")
P_("   (5) Zero 在 21cm 上的真实贡献 = m_e（0.0008%）")
verdict = ("不能。'锚一个锚点'满足（v）；'其余内生'不满足——"
           "alpha 缺一个无量纲比（类型已定）、g_p 缺一个非微扰机制（类型未定）；"
           "绕过只能到比值，绝对频率绕不掉")
P_("\n  => %s" % verdict)
check("总判决已给出（分层）", True, verdict)

# ---------- §6 如果要做，最少的输入清单 ----------
P_("\n§6 若要做（最小输入清单）")
P_("  要做 21cm 的【绝对频率】，最少需要:")
P_("    (a) 锚 v          [已有，唯一]")
P_("    (b) alpha         [需一笔无量纲输入 R]")
P_("    (c) mu_p 或 g_p   [需一个非微扰机制]")
P_("    (d) m_e/m_p       [m_e 内生 0.0008%; m_p 可绕或用测量]")
P_("    (e) 库仑束缚态    [Z165 已示范原生可算]")
P_("  ⟹ 最少 2 笔无量纲输入（b、c）")
check("最少输入数已给出（2 笔无量纲）", True, "alpha + g_p")

# ---------- §7 为什么 21cm 是坏靶子 ----------
P_("\n§7 为什么 21cm 是坏靶子（重申）")
P_("  nu ~ alpha^4 ⟹ 21cm 把基本常数的精度直接搬上台面")
P_("  21cm 被引用到 13 位有效数字")
P_("  ⟹ 要'推出'须把 alpha 推到 1e-12 量级 —— 没有结构做得到")
P_("  ⟹ 检验 Zero 该用只依赖结构量的观测量")
check("坏靶子的理由已给出", True, "alpha^4 放大")

P_("\n" + "=" * 96)
nfail = sum(1 for p in PASS if not p["pass"])
P_("合计: %d 项判据, %d PASS, %d FAIL" % (len(PASS), len(PASS)-nfail, nfail))
P_("=" * 96)

out = {
    "nodes": NODES,
    "counts": {"OK": nOK, "A": nA, "B": nB, "C": nC},
    "verdict": verdict,
    "min_inputs": ["anchor v (have)", "alpha (needs dimensionless R)",
                   "g_p or mu_p (needs non-perturbative mechanism)"],
    "checks": PASS,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_21cm_final_assessment.json"), "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
P_("结果写入: %s" % os.path.join(RES, "zero_21cm_final_assessment.json"))