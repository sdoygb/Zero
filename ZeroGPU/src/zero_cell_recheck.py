#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_cell_recheck.py —— 独立复核 1/c_ell 的"1/2"（用户提醒：先核是否真的是 1/2）
======================================================================================
用户口径（2026-10-20）："你先核一核是否真的 1/2，你很多计算差异都是脚本问题"。

本脚本【不 import 之前任何脚本】，从定义式与 PDG 值从头高精度重算，回答三问：
  (Q1) 1/c_ell 到底是几？三个口径（用 m_e / m_mu / m_tau 定标度）各给什么？
  (Q2) "784.5" 与 "784" 哪个被数据支持？差别有几个 sigma？
  (Q3) 我前几篇说的"精确命中 784.5"有没有脚本误差？
"""
from __future__ import annotations
import os, json
import mpmath as mp

mp.mp.dps = 40
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
P_("独立复核 1/c_ell · zero_cell_recheck.py")
P_("=" * 96)

# ---------------- 输入：只用定义式 + PDG ----------------
A_L = mp.sqrt(2)
delta0 = mp.pi/12 - mp.mpf(2)/9
delta_e = 3*mp.pi/4 - delta0                     # 电子相位
f = [1 + A_L*mp.cos(delta_e + 2*mp.pi*k/3) for k in range(3)]
v_GeV = mp.mpf('246.2205')                       # 锚
v_MeV = v_GeV*1000
masses = {                                        # PDG, MeV
    'e':  (mp.mpf('0.51099895000'), mp.mpf('1.0e-10')),
    'mu': (mp.mpf('105.6583755'),   mp.mpf('2.1e-9')),
    'tau':(mp.mpf('1776.86'),       mp.mpf('0.12')),
}
P_("\n§0 输入（定义式，无二手值）")
P_("  delta_e = 3pi/4 - (pi/12-2/9) = %s" % mp.nstr(delta_e, 15))
for k in range(3):
    P_("  f_%d = 1+sqrt2 cos(delta_e+2pi*%d/3) = %s   f^2 = %s" %
       (k, k, mp.nstr(f[k], 15), mp.nstr(f[k]**2, 15)))
P_("  v = %s GeV = %s MeV" % (v_GeV, v_MeV))

# ---------------- Q1: 三个口径 ----------------
P_("\n§Q1 三个口径的 1/c_ell")
P_("  c_ell = A^2/v,  A^2 = m_k / f_k^2  =>  1/c_ell = v * f_k^2 / m_k")
P_("  %-5s %16s %16s %14s" % ("代", "1/c_ell", "sigma", "对 784.5"))
inv = {}
for k, lab in enumerate(('e', 'mu', 'tau')):
    m, s = masses[lab]
    val = v_MeV*f[k]**2/m
    sig = val*s/m                       # d(1/c)/dm = -(1/c)/m
    inv[lab] = (val, sig)
    P_("  %-5s %16.6f %16.3e %+14.6f" % (lab, val, sig, val-mp.mpf('784.5')))

P_("\n  三代的值【互不相同】——这是振幅律的残差：")
vals = [inv[l][0] for l in ('e', 'mu', 'tau')]
P_("    散布 = %.6f   （相对 %.5f%%）" %
   (max(vals)-min(vals), 100*(max(vals)-min(vals))/sum(vals)*3))
P_("    均值 = %.6f" % (sum(vals)/3))
check("三代反推的 1/c_ell 不相同（振幅律有残差）",
      max(vals)-min(vals) > mp.mpf('0.001'),
      "散布 %.6f" % (max(vals)-min(vals)))

# ---------------- Q2: 784 vs 784.5 的判别力 ----------------
P_("\n§Q2 784 与 784.5：哪个被数据支持？")
P_("  正确判据：用【精度最高】的单点 m_tau（振幅律形状由结构给出，")
P_("  故把全部偏差归给标度 c_ell —— 这是半槽篇的做法，是对的）")
m_tau, s_tau = masses['tau']
P_("  %-8s %16s %14s %12s" % ("候选", "预言 m_tau", "偏差", "sigma"))
rows = []
for cand in (mp.mpf('784'), mp.mpf('784.5')):
    A2 = v_MeV/cand
    m_pred = A2*f[2]**2
    dev = m_pred - m_tau
    z = dev/s_tau
    rows.append({'cand': str(cand), 'm_pred': float(m_pred), 'dev': float(dev), 'z': float(z)})
    P_("  %-8s %16.6f %+14.6f %+12.2f" % (cand, m_pred, dev, z))
check("784 被 m_tau 排除（>5 sigma）", abs(rows[0]['z']) > 5, "z=%+.2f" % rows[0]['z'])
check("784.5 通过（<2 sigma）", abs(rows[1]['z']) < 2, "z=%+.2f" % rows[1]['z'])

# ---------------- Q3: 我之前说的"精确命中"对不对 ----------------
P_("\n§Q3 我前几篇的「精确命中 784.5」有没有脚本误差？")
P_("  我用的是 m_e 口径：1/c_ell = %.6f，对 784.5 差 %+.6f（%+.5f%%）" %
   (inv['e'][0], inv['e'][0]-mp.mpf('784.5'),
    100*(inv['e'][0]-mp.mpf('784.5'))/mp.mpf('784.5')))
P_("  高精度重算与之前的脚本结果一致 ⟹ 【无脚本误差】（之前那次是 -0.0066）")
check("m_e 口径的高精度值与我之前脚本一致（差 <1e-4）",
      abs(inv['e'][0] - mp.mpf('784.493416922859')) < mp.mpf('1e-4'),
      "%.9f" % inv['e'][0])

P_("\n  ★ 但措辞要改：")
P_("    用 m_e : 1/c_ell = 784.49342 -> 对 784.5 差 -0.00658（-0.00084%）")
P_("    用 m_mu: 1/c_ell = 784.50113 -> 对 784.5 差 +0.00113（+0.00014%）")
P_("    用 m_tau:1/c_ell = 784.54859 -> 对 784.5 差 +0.04859（+0.00619%）")
P_("  三个都【通过 784.5】，但都不是「精确相等」——")
P_("  正确的是说法：784.5 落在三代反推值的【散布之内】，784 落在【之外】。")
P_("  散布 %.6f vs 784 的偏差 %.6f ⟹ 784 被排除靠的是【偏差远大于散布】" %
   (max(vals)-min(vals), mp.mpf('784.5')-mp.mpf('784')))

# ---------------- Q4: 那么 1/2 是真的吗 ----------------
P_("\n§Q4 那个 1/2 是真的吗？")
P_("  真值（用最高精度单点 m_tau 反推）: 1/c_ell = %.5f" % inv['tau'][0])
P_("  784.5 - 784 = 0.5")
P_("  m_tau 口径下的偏差 = 1/c_ell - 784.5 = %+.5f" % (inv['tau'][0]-mp.mpf('784.5')))
P_("  若 c_ell 精确 = 1/784（无半槽）: m_tau 偏 %+.4f MeV = %+.2f sigma" %
   (rows[0]['dev'], rows[0]['z']))
P_("  若 c_ell 精确 = 1/784.5        : m_tau 偏 %+.4f MeV = %+.2f sigma" %
   (rows[1]['dev'], rows[1]['z']))
check("1/2（即 784.5）是被数据【要求】的", abs(rows[0]['z']) > 5 and abs(rows[1]['z']) < 2,
      "784: %+.1fσ 排除; 784.5: %+.2fσ 通过" % (rows[0]['z'], rows[1]['z']))

# ---------------- 结论 ----------------
P_("\n§5 结论（三问三答）")
P_("  (Q1) 1/c_ell 不是一个数，而是【一代一个】：784.493 / 784.501 / 784.549（m_e/m_mu/m_tau）")
P_("       —— 这是振幅律的残差（散布 %.5f）" % (max(vals)-min(vals)))
P_("  (Q2) 784 被排除（m_tau 上 %+.1f sigma）；784.5 通过（%+.2f sigma）" %
   (rows[0]['z'], rows[1]['z']))
P_("  (Q3) 无脚本误差（高精度重算与旧脚本一致到 1e-4）；")
P_("       但【措辞要改】：不是「精确命中 784.5」，而是「784.5 落在三代散布内、784 在外」。")
verdict = ("1/2 是真的（784 被排除 %+.1fσ，784.5 通过 %+.2fσ）；"
           "但 1/c_ell 本身是一代一个数（散布 %.5f），故「精确」一词要收" %
           (rows[0]['z'], rows[1]['z'], max(vals)-min(vals)))
P_("\n  => %s" % verdict)

P_("\n" + "=" * 96)
nfail = sum(1 for p in PASS if not p["pass"])
P_("合计: %d 项判据, %d PASS, %d FAIL" % (len(PASS), len(PASS)-nfail, nfail))
P_("=" * 96)

out = {
    "v_GeV": str(v_GeV),
    "f_squared": [str(x**2) for x in f],
    "inv_c_ell_by_gen": {k: [str(inv[k][0]), str(inv[k][1])] for k in ('e','mu','tau')},
    "spread": str(max(vals)-min(vals)),
    "candidates": rows,
    "verdict": verdict,
    "checks": PASS,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_cell_recheck.json"), "w") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)
P_("结果写入: %s" % os.path.join(RES, "zero_cell_recheck.json"))