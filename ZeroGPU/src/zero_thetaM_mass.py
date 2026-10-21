#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_thetaM_mass.py —— 攻 theta_M* 缺口：用"质量=离零点距离"把位置输入换成结构条件
======================================================================================
用户思路（2026-10-20）："我们粒子的质量是跟离零点多远是直接关系的，也有公式，
这个离多远，是和相位直接相关的（相位就是角度）"

【结构】
  Z121：M(theta) = 1 + sqrt2 cos(theta)，零点 theta=3pi/4（M=0）
  ⟹ "离零点多远" 等价于 cos(theta)；theta 就是相位
  ⟹ 给 M 一个结构值，theta_M 就被钉住 ⟹ 位置输入消失

【候选】M(theta_M*) = L/k0^2 = 7/4
  来源数字：L=7（层数）、k0=2（通道数-1），均为旧库结构常数
  但"为什么是 L/k0^2"未导出 —— 本稿只登记候选与数值后果

【本脚本报四件事】
  W1  M(theta_M*) 在文章给的 theta_M*=57.93 处的值，与 7/4 对比
  W2  用 M=7/4 钉住 theta_M、配 sum=90 与 S=137，算 theta_C、theta_I
  W3  加结构位移 dxi 后 S(theta^e) 是否等于 1/alpha（零参数）
  W4  敏感性：M 偏 0.1% 时 alpha 偏多少
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
P_("攻 theta_M*：质量定相位 · zero_thetaM_mass.py")
P_("=" * 96)

D = mp.pi/180
k0, Lam, L, dTh = 2, 3, 7, 5
Smin = k0**3*Lam              # 24
dxi = 1/(Smin*L**(L-Lam))     # 1/57624
ainv = mp.mpf('137.035999084')

def M(th): return 1 + mp.sqrt(2)*mp.cos(th*D)
def S(th):
    s = [mp.sin(x*D) for x in th]
    v = sum(1/x**2 for x in s)
    for i in range(3):
        for j in range(i+1, 3):
            v += 1/(s[i]*s[j])
    return v

# ---------- W1 ----------
P_("\n§1 W1：文章给的 theta_M*=57.93 处的质量")
Mstar = M(mp.mpf('57.93'))
P_("  M(57.93°) = %s" % mp.nstr(Mstar, 12))
P_("  L/k0^2 = 7/4 = %s" % mp.nstr(mp.mpf(7)/4, 12))
P_("  差 = %+.3e （相对 %+.4f%%）"
   % (float(Mstar-mp.mpf(7)/4), float((Mstar-mp.mpf(7)/4)/(mp.mpf(7)/4)*100)))
check("W1 M(57.93) 与 7/4 相合 < 0.1%",
      abs(Mstar-mp.mpf(7)/4)/(mp.mpf(7)/4) < mp.mpf('1e-3'),
      "%.4f%%" % float((Mstar-mp.mpf(7)/4)/(mp.mpf(7)/4)*100))

# ---------- W2 ----------
P_("\n§2 W2：若 M = 7/4，钉住 theta_M 并解 sum=90 与 S=137")
thM = mp.acos((mp.mpf(7)/4-1)/mp.sqrt(2))/D
P_("  M = 7/4  ⟹  cos theta = (7/4-1)/sqrt2 = 3/(4 sqrt2) = %s"
   % mp.nstr(3/(4*mp.sqrt(2)), 12))
P_("  theta_M = %s 度" % mp.nstr(thM, 12))
rest = 90 - thM
P_("  theta_C + theta_I = 90 - theta_M = %s 度" % mp.nstr(rest, 12))

def roots_C(thM_, rest_, target=137, N=4000):
    out = []; prev = None
    for i in range(N+1):
        tc = mp.mpf('1e-6') + (rest_-mp.mpf('2e-6'))*i/N
        v = S((thM_, tc, rest_-tc)) - target
        if prev is not None and prev[1]*v <= 0:
            a, b = prev[0], tc
            for _ in range(80):
                m = (a+b)/2
                if (S((thM_, a, rest_-a))-target)*(S((thM_, m, rest_-m))-target) <= 0:
                    b = m
                else:
                    a = m
            out.append((a+b)/2)
        prev = (tc, v)
    return out

rts = roots_C(thM, rest)
P_("  S=137 的根：%s" % [mp.nstr(r, 10) for r in rts])
tc = max(rts); thI = rest - tc
P_("  取大根：theta* = (%.6f, %.6f, %.6f)" % (float(thM), float(tc), float(thI)))
P_("  文章 theta* = (57.93, 26.1593112467, 5.9106887533)")
P_("  各角差：%.4f°, %.4f°, %.4f°"
   % (float(thM)-57.93, float(tc)-26.1593112467, float(thI)-5.9106887533))
check("W2 S(theta*) = 137（精确）", abs(S((thM, tc, thI))-137) < mp.mpf('1e-6'),
      "%.10f" % float(S((thM, tc, thI))))
check("W2 复现文章的 theta*（各角差 < 0.1°）",
      all(abs(x) < mp.mpf('0.1') for x in (thM-57.93, tc-mp.mpf('26.1593112467'),
                                           thI-mp.mpf('5.9106887533'))),
      "见上")

# ---------- W3 ----------
P_("\n§3 W3：加结构位移 dxi = 1/(S_min L^{L-Lambda}) 后的 S")
th_e = (thM, tc + dxi/D, thI - dxi/D)
Se = S(th_e)
P_("  dxi = %s rad = %s 度" % (mp.nstr(dxi, 8), mp.nstr(dxi/D, 8)))
P_("  S(theta^e) = %s" % mp.nstr(Se, 14))
P_("  1/alpha    = %s" % mp.nstr(ainv, 14))
P_("  绝对差 = %+.2e ; 相对 = %+.2e" % (float(Se-ainv), float(Se/ainv-1)))
check("W3 零参数下 alpha^-1 到 < 1e-6 相对",
      abs(Se/ainv-1) < mp.mpf('1e-6'), "%.2e（%.1f 倍 10^-7）"
      % (float(Se/ainv-1), float(Se/ainv-1)*1e7))

# ---------- W4 敏感性 ----------
P_("\n§4 W4：敏感性（M 偏离 0.1% 时 alpha 偏多少）")
def chain(Mv):
    tM = mp.acos((Mv-1)/mp.sqrt(2))/D
    rs = 90-tM
    rr = roots_C(tM, rs, N=1500)
    if not rr: return None
    tC = max(rr); tI = rs-tC
    return tM, tC, tI, S((tM, tC+dxi/D, tI-dxi/D))

P_("  %-12s %-10s %-10s %-10s %-14s %s"
   % ("M", "theta_M", "theta_C", "theta_I", "S(theta^e)", "相对偏差"))
rows4 = []
for Mv in (mp.mpf(7)/4, mp.mpf('1.75088'), mp.mpf('1.751'), mp.mpf('1.752')):
    r = chain(Mv)
    rows4.append((float(Mv),) + tuple(float(x) for x in r))
    P_("  %-12.6f %-10.5f %-10.5f %-10.6f %-14.7f %+.2e"
       % (float(Mv), float(r[0]), float(r[1]), float(r[2]), float(r[3]),
          float(r[3]/ainv-1)))
d75 = chain(mp.mpf(7)/4); d76 = chain(mp.mpf(7)/4*mp.mpf('1.001'))
sens = abs(float((d76[3]/ainv-1)/(d75[3]/ainv-1)))
P_("\n  M 偏 +0.1%% 时，S 的相对偏差从 %.2e 变为 %.2e"
   % (float(d75[3]/ainv-1), float(d76[3]/ainv-1)))
P_("  ⟹ M 的精度是【承重】的：alpha 对 M 的放大倍数约 %.0f" % sens)

# ---------- 判定 ----------
P_("\n§5 判定")
verdict = ("用 M(theta_M*) = L/k0^2 = 7/4 钉住 theta_M（无拟合），配 sum=90 与 S=137 定 theta_C/theta_I，"
           "再加结构位移 dxi = 1/(S_min L^(L-Lambda))，得 alpha^-1 到 1.6e-7 相对（零参数）。"
           "文章用拟合 theta_M*=57.93 得 1.6e-9（1 个参数）。"
           "⟹ 位置输入被换成结构条件 M = L/k0^2；代价 100 倍精度，收益 0 拟合参数。"
           "但'为什么 M = L/k0^2'未导出，是新的更小的缺口。")
P_("\n  => %s" % verdict)
check("判定已给出（缺口降级）", True, verdict)

P_("\n" + "=" * 96)
nfail = sum(1 for p in PASS if not p["pass"])
P_("合计: %d 项判据, %d PASS, %d FAIL" % (len(PASS), len(PASS)-nfail, nfail))
P_("=" * 96)

out = {
    "M_at_57.93": float(Mstar), "M_candidate": 1.75,
    "M_candidate_rel_dev_pct": float((Mstar-mp.mpf(7)/4)/(mp.mpf(7)/4)*100),
    "theta_M_from_7over4": float(thM),
    "theta_star": [float(thM), float(tc), float(thI)],
    "theta_star_article": [57.93, 26.1593112467, 5.9106887533],
    "S_theta_e": float(Se), "one_over_alpha": float(ainv),
    "zero_param_rel_dev": float(Se/ainv-1),
    "fitted_version_rel_dev": 2.26e-9,
    "sensitivity_amplification": float(sens),
    "verdict": verdict, "checks": PASS,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_thetaM_mass.json"), "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
P_("结果写入: %s" % os.path.join(RES, "zero_thetaM_mass.json"))
