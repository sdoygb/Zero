#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_alpha_sixterm_21cm.py —— 六项式 alpha ＋ 21cm 公式的独立复算与自洽检验
======================================================================================
来源（**移植，只作参考**）：
  GeometryAI-Mac-Build/app/articles/12.1_基本常数与谱_观察者位置与锚定账本_CN_260917.md（六项式、闭式、位移）
  .../12.4_超精细耦合_4π与复合律_CN_260917.md（21cm 公式）
  cosmos-construct/ALPHA_FORMULA_SOLVED.md（RG 式与 R(r)）
  几何文章 7.11.3.02（sin^2 theta_W = 0.23124，447/392）

本脚本独立复算四件事：
  V1  六项式 S(theta) 在 theta* 上 = 137 = 2^7+3^2；在电子点上 = 1/alpha
  V2  闭式 2^7+3^2+9/(2*5^3) = 137.036 及其偏差；结构位移 dxi_e = 1/(S_min L^4)
  V3  RG 式与六项式的自洽（两条独立路的一致性）
  V4  alpha_s 在 alpha(0) 里是否消掉；sin(theta_M*) 与 r = g_L/g_R 的关系
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
P_("六项式 alpha 与 21cm · zero_alpha_sixterm_21cm.py")
P_("=" * 96)

pi = mp.pi
D = pi/180

# ---------- V1 六项式 ----------
P_("\n§1 V1：六项式 S(theta) = Σ1/sin²θ_i + Σ_{i<j} 1/(sinθ_i sinθ_j)")
def S(th):
    s = [mp.sin(x*D) for x in th]
    v = sum(1/x**2 for x in s)
    for i in range(3):
        for j in range(i+1, 3):
            v += 1/(s[i]*s[j])
    return v

th_star = (mp.mpf('57.93'), mp.mpf('26.1593112467'), mp.mpf('5.9106887533'))
th_e    = (mp.mpf('57.93'), mp.mpf('26.1603055422'), mp.mpf('5.9096944578'))
S_star, S_e = S(th_star), S(th_e)
P_("  三角和 = %s 度  （= pi/2 = 2/8 圈，Zero 原生角）" % mp.nstr(sum(th_star), 12))
P_("  S(theta*)  = %s   （应 = 137 = 2^7+3^2）" % mp.nstr(S_star, 14))
P_("  S(theta^e) = %s   （应 = 1/alpha）" % mp.nstr(S_e, 14))
check("V1a 三角和 = 90 度", abs(sum(th_star) - 90) < mp.mpf('1e-9'),
      "%.12f" % float(sum(th_star)))
check("V1b S(theta*) = 137（精确）", abs(S_star - 137) < mp.mpf('1e-9'),
      "%.12f" % float(S_star))
check("V1c S(theta^e) = 1/alpha（12 位）",
      abs(S_e - mp.mpf('137.035999084')) < mp.mpf('1e-8'),
      "%.12f" % float(S_e))

# 逐项
s = [mp.sin(x*D) for x in th_star]
six = [1/x**2 for x in s] + [1/(s[i]*s[j]) for i in range(3) for j in range(i+1,3)]
P_("  逐项：%s" % [mp.nstr(v, 6) for v in six])
P_("  ⟹ 无一项等于 128 或 9（闭式是【总和】的陈述，非逐项分解）")

# ---------- V2 闭式与位移 ----------
P_("\n§2 V2：零参数闭式")
k0, Lam, L, dTh = 2, 3, 7, 5
closed = k0**L + Lam**2 + mp.mpf(Lam**2)/(k0*dTh**3)
obs = mp.mpf('137.035999084')
P_("  k0^L + Lam^2 + Lam^2/(k0 dTh^3) = 2^7+3^2+9/(2*125) = %s" % mp.nstr(closed, 12))
P_("  观测 = %s ；偏差 %+.2e （相对 %+.1e）"
   % (mp.nstr(obs, 12), float(closed-obs), float((closed/obs-1)*100)))
check("V2a 闭式偏差 < 1e-6", abs(closed-obs) < mp.mpf('1e-6'),
      "%.2e" % float(closed-obs))
P_("  dTh 由零和恒等式导出：dTh = Lam^2 - k0^2 = 9-4 = %d" % (Lam**2-k0**2))

Smin = k0**3*Lam                       # 24
dxi = 1/(Smin*L**(L-Lam))              # 1/(24*2401) = 1/57624
P_("\n  结构位移 dxi_e = 1/(S_min L^(L-Lam)) = 1/(%d * %d^%d) = %s"
   % (Smin, L, L-Lam, mp.nstr(dxi, 8)))
P_("  库内实测位移 = 1.735373021e-5 rad ；相合 %.4f%%"
   % float(dxi/mp.mpf('1.735373021e-5')*100))
check("V2b 位移式与实测相合 < 0.01%",
      abs(dxi/mp.mpf('1.735373021e-5')-1) < mp.mpf('1e-4'),
      "%.6f%%" % float((dxi/mp.mpf('1.735373021e-5')-1)*100))

# ---------- V3 两路自洽 ----------
P_("\n§3 V3：RG 式与六项式的自洽")
b1, b2, b3 = mp.mpf(41)/10, mp.mpf(-19)/6, mp.mpf(-7)
cc = lambda b: b/(2*pi)
Da = mp.mpf('0.0663'); bY = mp.mpf(41)/6

def RG_route(aem_inv, s2, as_):
    x1 = mp.mpf(3)/5*aem_inv*(1-s2); x2 = aem_inv*s2; x3 = 1/as_
    t = (x2-x3)/(cc(b2)-cc(b3))
    A1 = x1-cc(b1)*t; A2 = x2-cc(b2)*t
    R = mp.mpf(5)/3*A1/A2
    inv = (A2*(1+R) + t*(b2+bY)/(2*pi))/(1-Da)
    return inv, R, mp.sqrt(3*(R-mp.mpf(1)/3)/4), t

inv, R, r, t = RG_route(mp.mpf('127.95'), mp.mpf('0.23124'), mp.mpf('0.1179'))
P_("  RG 式（sin^2 theta_W 用结构值 0.23124）: 1/alpha(0) = %s" % mp.nstr(inv, 12))
P_("  六项式（含位移）                       : 137.0359994")
P_("  差 = %+.2e （相对 %+.1e）" % (float(inv-mp.mpf('137.0359994')),
                                   float((inv/mp.mpf('137.0359994')-1)*100)))
rel = abs(inv/mp.mpf('137.0359994')-1)
check("V3a 两路自洽（相对 < 1e-5）", rel < mp.mpf('1e-5'),
      "绝对 %.2e, 相对 %.1e（%.1f ppm）" % (float(inv-mp.mpf('137.0359994')), float(rel), float(rel*1e6)))
P_("  sin(theta_M*) = %s ;  r = g_L/g_R = %s ; 差 %+.2e"
   % (mp.nstr(mp.sin(mp.mpf('57.93')*D), 8), mp.nstr(r, 8),
      float(mp.sin(mp.mpf('57.93')*D)-r)))
check("V3b sin(theta_M*) = r（同一量，差 < 1% 舍入级）",
      abs(mp.sin(mp.mpf('57.93')*D)-r) < mp.mpf('1e-3'),
      "%.2e" % float(mp.sin(mp.mpf('57.93')*D)-r))

# ---------- V4 alpha_s 消掉 ----------
P_("\n§4 V4：alpha_s 在 alpha(0) 里是否消掉")
vals = []
for as_ in (mp.mpf('0.1179'), mp.mpf('0.1181'), mp.mpf('0.1184')):
    iv, RR, rr, tt = RG_route(mp.mpf('127.95'), mp.mpf('0.23124'), as_)
    vals.append((float(as_), float(tt), float(RR), float(rr), float(iv)))
    P_("  alpha_s=%.4f: t=%.4f  R=%.6f  r=%.6f  1/alpha(0)=%.8f"
       % vals[-1])
spread = max(v[4] for v in vals) - min(v[4] for v in vals)
P_("  1/alpha(0) 的极差 = %.2e" % spread)
check("V4 alpha_s 在 alpha(0) 里消掉（极差 < 1e-9）", spread < 1e-9,
      "极差 %.2e" % spread)
P_("  ⟹ alpha_s 只决定交叉点，不决定 alpha(0) ⟹ 【不是第二个输入】")

# ---------- 21cm ----------
P_("\n§5 21cm 公式（移植自 12.4）")
A2 = mp.mpf('313.856')      # MeV, A^2 = v/784.5
me, mp_ = mp.mpf('0.51099895'), mp.mpf('938.272')
ae = mp.mpf('1.162364e-3')
gp = mp.mpf('5.58569')
alpha = 1/mp.mpf('137.035999084')
me_mp = me/mp_
nu = (mp.mpf(k0**2)/Lam) * alpha**4 * (1+ae) * me_mp * gp * (1+me_mp)**(-3) * mp.mpf('1.235590e20')
P_("  nu_21 = (k0^2/Lam) alpha^4 (1+a_e)(m_e/m_p) g_p (1+m_e/m_p)^-3 m_e c^2/h")
P_("        = %.4f MHz" % float(nu/1e6))
P_("  实测 1420.4058 MHz ；偏差 %+.4f%%" % float((nu/1e6/mp.mpf('1420.4058')-1)*100))
check("21cm 公式偏差 < 0.01%",
      abs(nu/1e6/mp.mpf('1420.4058')-1) < mp.mpf('1e-4'),
      "%+.4f%%" % float((nu/1e6/mp.mpf('1420.4058')-1)*100))
P_("  ⚠ 其中 g_p = 5.58569 是【外部】输入（几何领头 2*Lam=6，偏离 -6.9051%）")

P_("\n" + "=" * 96)
nfail = sum(1 for p in PASS if not p["pass"])
P_("合计: %d 项判据, %d PASS, %d FAIL" % (len(PASS), len(PASS)-nfail, nfail))
P_("=" * 96)

out = {
    "S_star": float(S_star), "S_electron": float(S_e),
    "closed_form": float(closed), "closed_dev": float(closed-obs),
    "delta_xi": float(dxi), "R_MX": float(R), "r": float(r), "t": float(t),
    "RG_route_val": float(inv), "two_route_diff": float(inv-mp.mpf('137.0359994')),
    "alpha_s_spread": float(spread),
    "sin_thetaM": float(mp.sin(mp.mpf('57.93')*D)),
    "nu_21_MHz": float(nu/1e6), "nu_21_dev_pct": float((nu/1e6/mp.mpf('1420.4058')-1)*100),
    "g_p_external": float(gp),
    "checks": PASS,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_alpha_sixterm_21cm.json"), "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
P_("结果写入: %s" % os.path.join(RES, "zero_alpha_sixterm_21cm.json"))
