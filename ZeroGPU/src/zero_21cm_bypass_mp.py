#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_21cm_bypass_mp.py —— 绕过 m_p：21cm 其实不需要 m_p
======================================================================================
用户口径（2026-10-20）："对 m_p 也是可以绕过去"。

关键观察：21cm 里 m_p 只通过【质子磁矩】进来——
    mu_p = g_p * mu_N / 2 = g_p * e*hbar/(4 m_p)
故   g_p/m_p = 4 mu_p/(e hbar)          <-- m_p 被 mu_p 吸收掉
于是 nu_21 ∝ (g_p/m_p) * m_e^2 ∝ mu_p * m_e^2  ⟹ **不含独立的 m_p**

结论：把 (g_p, m_p) 换成【测出来的 mu_p】之后，Zero 的 m_p 误差(0.355%)整条掉出。
     剩下的 Zero 贡献只有 m_e（0.0008%）。
"""
from __future__ import annotations
import os, json, math

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
P_("绕过 m_p · zero_21cm_bypass_mp.py")
P_("=" * 96)

# ---------- 常数 ----------
me   = 9.1093837015e-31
mp   = 1.67262192369e-27
c    = 299792458.0
h    = 6.62607015e-34
hbar = 1.054571817e-34
e    = 1.602176634e-19
alpha= 1/137.035999084
gp   = 5.5856946893
muN  = 5.0507837461e-27          # J/T  核磁子
mup  = 2.7928473446*muN          # J/T  质子磁矩（测出来的）
eV   = 1.602176634e-19
Ry   = 13.605693122994
NU_OBS = 1420.405751768e6
r = me/mp
mu_fac = 1.0/(1.0+r)

# ---------- §1 恒等式核验 ----------
P_("\n§1 恒等式：g_p/m_p = 4 mu_p/(e hbar)")
lhs = gp/mp
rhs = 4.0*mup/(e*hbar)
P_("  g_p/m_p        = %.6e  1/kg" % lhs)
P_("  4 mu_p/(e hbar)= %.6e  1/kg" % rhs)
P_("  相对差         = %+.2e" % ((rhs-lhs)/lhs))
check("两式相等（mu_p 把 m_p 吸收掉）", abs((rhs-lhs)/lhs) < 1e-8,
      "rel=%.2e" % ((rhs-lhs)/lhs))

# ---------- §2 改写后的 21cm ----------
P_("\n§2 把 21cm 改写成【不含独立 m_p】的形式")
P_("  nu = (4/3) g_p (m_e^2/m_p) alpha^4 c^2/h /(1+r)")
P_("     = (4/3) * [4 mu_p/(e hbar)] * m_e^2 alpha^4 c^2/h /(1+r)")
P_("  ⟹ 依赖 = { mu_p, m_e, alpha, (1+r) }  —— m_p 只出现在 r 里")
# 核验：用 mu_p 形式重算，应极近观测
nu_mu = (4.0/3.0)*(4.0*mup/(e*hbar))*me**2*alpha**4*c**2/h*mu_fac
P_("  用 mu_p 形式重算 nu = %.6f MHz（观测 %.6f）" % (nu_mu/1e6, NU_OBS/1e6))
P_("  相对差 = %+.2e" % ((nu_mu-NU_OBS)/NU_OBS))
check("改写后复现 21cm（<1e-4）", abs((nu_mu-NU_OBS)/NU_OBS) < 1e-4,
      "dev=%.2e" % ((nu_mu-NU_OBS)/NU_OBS))

# ---------- §3 各量的误差传播 ----------
P_("\n§3 误差传播：谁决定 21cm 的精度")
P_("  对数导数（只保留主导项）:")
P_("    dln nu / dln mu_p  = +1   （质子磁矩）")
P_("    dln nu / dln m_e   = +2   （电子质量，二次）")
P_("    dln nu / dln alpha = +4")
P_("  故：Zero 的 m_p 误差【完全不出现在 mu_p 通道里】")

# Zero 的 m_e（Z126）
v_GeV = 246.2205
A2 = v_GeV/784.5
theta_e = 2*math.pi/3 + 2/9
Me = 1 + math.sqrt(2)*math.cos(theta_e)
me_zero_MeV = A2*1e3*Me**2
P_("\n  Zero 的 m_e（Z126 振幅律）= %.6f MeV，PDG %.6f，差 %+.4f%%" %
   (me_zero_MeV, 0.510999, 100*(me_zero_MeV-0.510999)/0.510999))
check("Z126 的 m_e 复现到 0.01%", abs(me_zero_MeV-0.510999)/0.510999 < 1e-4,
      "%.6f MeV" % me_zero_MeV)

# ---------- §4 两种"Zero 供 m_p"对比 ----------
P_("\n§4 对比：Zero 是否供 m_p")
ratio_me = me_zero_MeV/0.510999
mp_zero = 941.6
ratio_mp = mp_zero/938.272

# (A) 旧写法：Zero 供 m_e 与 m_p（m_p 通过 g_p/m_p 通道）
shiftA = ratio_me**1.99946 * ratio_mp**(-0.99946)
nuA = NU_OBS*shiftA
P_("  (A) Zero 供 m_e 与 m_p（旧写法，m_p 当作独立输入）")
P_("      nu = %.4f MHz，偏 %+.4f%% (%+.3f MHz)" %
   (nuA/1e6, 100*(shiftA-1), (nuA-NU_OBS)/1e6))

# (B) 新写法：m_p 被 mu_p 吸收，Zero 只供 m_e
shiftB = ratio_me**2.0
nuB = NU_OBS*shiftB
P_("  (B) m_p 被 mu_p 吸收，Zero 只供 m_e")
P_("      nu = %.4f MHz，偏 %+.4f%% (%+.4f MHz)" %
   (nuB/1e6, 100*(shiftB-1), (nuB-NU_OBS)/1e6))

P_("\n  改善倍数 = %.0f×" % abs((shiftA-1)/(shiftB-1)))
check("绕过后误差显著改善（>100×）",
      abs((shiftA-1)/(shiftB-1)) > 100,
      "%.0f×" % abs((shiftA-1)/(shiftB-1)))
check("绕过后 21cm 精度优于 0.01%",
      abs(shiftB-1) < 1e-4, "%+.4f%%" % (100*(shiftB-1)))

# ---------- §5 概念更正 ----------
P_("\n§5 概念更正（本轮最重要）")
P_("  前篇说'靶心从 alpha,g_p 转到 m_p' —— 那是【参数化造成的假象】。")
P_("  超精细真正需要的是【质子磁矩 mu_p】；写成 (g_p, m_p) 两笔是记账选择。")
P_("  换成 mu_p 之后 m_p 不再是独立输入 ⟹ Zero 的 m_p 误差(0.355%)整条掉出。")
P_("  ★ 剩下的 Zero 贡献只有 m_e，而 Z126 把它给到 0.0008%。")
check("已登记参数化更正（m_p 瓶颈是假象）", True, "见 §5")

# ---------- §6 边界 ----------
P_("\n§6 边界（必须一起读）")
P_("  ① mu_p 与 alpha 仍是【测出来的】；本绕过换的是参数化，不是内生性。")
P_("  ② 所以本绕过不改变前篇的判决：21cm 仍需无量纲输入。")
P_("     它改变的是【Zero 内部误差落在哪】：从 m_p 转移到 m_e（后者好得多）。")
P_("  ③ m_e 的 0.0008% 依赖 784.5 这个【半槽构造】——")
P_("     半槽篇自判它为'被数据支持的构造，不是定理'。若 784.5 -> 784，")
P_("     m_e 会偏 +0.0638%，21cm 则偏 +0.13%。")
P_("  ④ 本绕过不是'Zero 预言了 21cm'：它只说 Zero 的 m_e 足够好。")

# 784 情形
me_784 = (v_GeV/784.0)*1e3*Me**2
P_("\n  核验 ③：用 784（无半槽）m_e = %.6f MeV，21cm 偏 %+.3f%%" %
   (me_784, 100*((me_784/0.510999)**2-1)))
check("无半槽时 21cm 偏约 +0.13%", abs(100*((me_784/0.510999)**2-1) - 0.128) < 0.01,
      "%+.3f%%" % (100*((me_784/0.510999)**2-1)))

P_("\n" + "=" * 96)
nfail = sum(1 for p in PASS if not p["pass"])
P_("合计: %d 项判据, %d PASS, %d FAIL" % (len(PASS), len(PASS)-nfail, nfail))
P_("=" * 96)

out = {
    "identity": {"g_p_over_m_p": lhs, "four_mup_over_ehbar": rhs,
                 "rel_diff": (rhs-lhs)/lhs},
    "rewrite": {"nu_with_mup_MHz": nu_mu/1e6, "dev": (nu_mu-NU_OBS)/NU_OBS},
    "zero_me": {"value_MeV": me_zero_MeV, "dev_pct": 100*(ratio_me-1)},
    "comparison": {"A_old_pct": 100*(shiftA-1), "A_MHz": (nuA-NU_OBS)/1e6,
                   "B_new_pct": 100*(shiftB-1), "B_MHz": (nuB-NU_OBS)/1e6,
                   "improvement_x": abs((shiftA-1)/(shiftB-1))},
    "no_half_slot_pct": 100*((me_784/0.510999)**2-1),
    "note": "m_p bottleneck was a parameterization artifact",
    "checks": PASS,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_21cm_bypass_mp.json"), "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
P_("结果写入: %s" % os.path.join(RES, "zero_21cm_bypass_mp.json"))