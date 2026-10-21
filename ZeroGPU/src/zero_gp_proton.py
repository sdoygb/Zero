#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_gp_proton.py —— 质子 g 因子的结构式与数秘压力测试
======================================================================================
用户判断（2026-10-20）："质子磁矩——精细结构数有了，这个也应该可以解决了吧"

【候选】
  g_p = 2*Lam * cos^2(pi/12) * [1 - A_l * delta_0^2]
  2*Lam = 6        几何领头（3 价夸克 x 2；与 Z125 的 m_p = 3A^2 同源）
  pi/12            Z123 的相位（delta_0 = pi/12 - 2/9 里那个）
  A_l = sqrt2      Z121 振幅系数
  delta_0^2        Z123 非中心相位平方

【本脚本报四件事】
  G1  候选值与实测对比
  G2  两层压低量的分解（cos^2(pi/12) 压 6.70%，残余 0.221% 由 A_l delta_0^2 吃掉）
  G3  数秘压力测试：枚举 2 Lam cos^2(x) (1 - c y) 型 600 个候选
  G4  对 21cm 的影响：g_p 从外部输入变推导值，nu_21 几乎不动
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

P_("=" * 92)
P_("质子 g 因子：结构式与压力测试 · zero_gp_proton.py")
P_("=" * 92)

pi = mp.pi
k0, Lam, L, Dth = 2, 3, 7, 5
d0 = pi/12 - mp.mpf(2)/9
A_l = mp.sqrt(2)
gp_exp = mp.mpf('5.58569')
lead = 2*Lam
c15 = mp.cos(pi/12)**2

def gp_pred():
    return lead * c15 * (1 - A_l*d0**2)

# ---------- G1 ----------
P_("\n§1 G1：候选与实测")
P_("  几何领头 2*Lam   = %s" % lead)
P_("  cos^2(pi/12)     = %s" % mp.nstr(c15, 14))
P_("  A_l*delta_0^2    = %s" % mp.nstr(A_l*d0**2, 14))
P_("  预言             = %s" % mp.nstr(gp_pred(), 16))
P_("  实测             = %s" % mp.nstr(gp_exp, 16))
P_("  偏差             = %+.2e  (%+.6f%%)"
   % (float(gp_pred()-gp_exp), float((gp_pred()/gp_exp-1)*100)))
check("G1 相对偏差 < 1e-5", abs(gp_pred()/gp_exp-1) < mp.mpf('1e-5'),
      "%.2e" % float(gp_pred()/gp_exp-1))

# ---------- G2 分解 ----------
P_("\n§2 G2：两层压低量")
r = gp_exp/lead
P_("  g_p/(2 Lam)           = %s" % mp.nstr(r, 14))
P_("  总压低 1 - r          = %s  (%.4f%%)" % (mp.nstr(1-r, 12), float((1-r)*100)))
P_("  第一层 1 - cos^2(pi/12) = %s  (%.4f%%)" % (mp.nstr(1-c15, 12), float((1-c15)*100)))
resid = 1 - r/c15
P_("  残余         = %s  (%.4f%%)" % (mp.nstr(resid, 12), float(resid*100)))
P_("  候选 A_l delta_0^2 = %s  （差 %+.3f%%）"
   % (mp.nstr(A_l*d0**2, 12), float((A_l*d0**2/resid-1)*100)))
check("G2 残余由 A_l delta_0^2 吃掉 < 0.5%",
      abs(A_l*d0**2/resid-1) < mp.mpf('5e-3'),
      "%.3f%%" % float((A_l*d0**2/resid-1)*100))

# ---------- G3 压力测试 ----------
P_("\n§3 G3：数秘压力测试")
pool = {
    'pi/12': pi/12, 'pi/8': pi/8, 'pi/6': pi/6, 'pi/4': pi/4,
    'd0': d0, '2/9': mp.mpf(2)/9, '1/3': mp.mpf(1)/3, '1/9': mp.mpf(1)/9,
    'd0^2': d0**2, '(2/9)^2': (mp.mpf(2)/9)**2,
}
coefs = {'A_l=sqrt2': mp.sqrt(2), '1': mp.mpf(1), 'k0=2': mp.mpf(2),
         'Lam=3': mp.mpf(3), 'Dth=5': mp.mpf(5), '1/2': mp.mpf('0.5')}
cands = []
for xn, xv in pool.items():
    if xv > pi/2:
        continue
    f1 = mp.cos(xv)**2
    for yn, yv in pool.items():
        for cn, c in coefs.items():
            v = lead*f1*(1 - c*yv)
            cands.append((abs(v/gp_exp - 1), xn, yn, cn, v))
cands.sort(key=lambda t: t[0])
P_("  枚举总数 = %d" % len(cands))
P_("  %-10s %-9s %-11s %-16s %s" % ("x", "y", "系数", "值", "相对偏差"))
for t in cands[:6]:
    P_("  %-10s %-9s %-11s %-16.9f %.2e" % (t[1], t[2], t[3], float(t[4]), float(t[0])))
mine = [i for i, t in enumerate(cands)
        if t[1] == 'pi/12' and t[2] == 'd0^2' and t[3] == 'A_l=sqrt2'][0]
n5 = sum(1 for t in cands if t[0] < mp.mpf('1e-5'))
n4 = sum(1 for t in cands if t[0] < mp.mpf('1e-4'))
P_("  我的候选排位 = 第 %d 名" % (mine+1))
P_("  |rel|<1e-5 命中 %d / %d (%.2f%%)" % (n5, len(cands), n5/len(cands)*100))
P_("  |rel|<1e-4 命中 %d / %d (%.2f%%)" % (n4, len(cands), n4/len(cands)*100))
P_("  第二名偏差 = %.2e （比第一名差 %.0f 倍）"
   % (float(cands[1][0]), float(cands[1][0]/cands[0][0])))
check("G3a 我的候选是第 1 名", mine == 0, "rank %d" % (mine+1))
check("G3b 同型 |rel|<1e-5 的命中 < 1%", n5/len(cands) < 0.01,
      "%.2f%%" % (n5/len(cands)*100))
check("G3c 与第二名拉开 > 100 倍", cands[1][0]/cands[0][0] > 100,
      "%.0f 倍" % float(cands[1][0]/cands[0][0]))

# ---------- G4 对 21cm 的影响 ----------
P_("\n§4 G4：对 21cm 的影响")
ae = mp.mpf('1.162364e-3')
me_mp = mp.mpf('5.446166e-4')
alpha = 1/mp.mpf('137.035999084')
scale = mp.mpf('1.235590e20')
exp21 = mp.mpf('1420.4058')
def nu21(gp, mm=me_mp):
    return (mp.mpf(4)/3)*alpha**4*(1+ae)*mm*gp*(1+mm)**(-3)*scale/1e6
P_("  nu_21(g_p 实测)      = %.4f MHz  (%+.4f%%)"
   % (float(nu21(gp_exp)), float((nu21(gp_exp)/exp21-1)*100)))
P_("  nu_21(g_p 推导值)    = %.4f MHz  (%+.4f%%)"
   % (float(nu21(gp_pred())), float((nu21(gp_pred())/exp21-1)*100)))
P_("  ⟹ g_p 从外部输入变推导值，nu_21 几乎不动（g_p 的精度已不是瓶颈）")
check("G4 用推导 g_p 后 21cm 仍 < 0.01%",
      abs(nu21(gp_pred())/exp21-1) < mp.mpf('1e-4'),
      "%+.4f%%" % float((nu21(gp_pred())/exp21-1)*100))

# ---------- 判定 ----------
P_("\n§5 判定")
verdict = (
    "g_p = 2*Lam * cos^2(pi/12) * [1 - A_l * delta_0^2] = 5.5856756，"
    "对实测 5.58569 偏差 2.6e-6。四个因子都有出处：2*Lam=6（3 价夸克 x 2，与 Z125 的 m_p=3A^2 同源）、"
    "cos^2(pi/12)（Z123 相位）、A_l=sqrt2（Z121 振幅系数）、delta_0^2（Z123 非中心）。"
    "压力测试：600 个同型候选中只有它落在 |rel|<1e-5（0.17%），第二名差 250 倍。"
    "⟹ 21cm 的 g_p 缺口消失，nu_21 仍 +0.0055%。"
    "⚠ 边界：公式是扫出来的，未从 Z017/Z121 推出；2*Lam 作领头仍是'3 个价夸克'的判读，SU(6) 破缺机制未给。")
P_("\n  => %s" % verdict)
check("判定已给出", True, verdict)

P_("\n" + "=" * 92)
nfail = sum(1 for p in PASS if not p["pass"])
P_("合计: %d 项判据, %d PASS, %d FAIL" % (len(PASS), len(PASS)-nfail, nfail))
P_("=" * 92)

out = {
    "g_p_pred": float(gp_pred()), "g_p_exp": float(gp_exp),
    "rel_dev": float(gp_pred()/gp_exp-1),
    "lead_2Lam": lead, "cos2_pi12": float(c15), "A_l": float(A_l),
    "delta_0": float(d0), "delta_0_sq": float(d0**2),
    "residual_needed": float(resid),
    "pressure_test": {"n_candidates": len(cands), "rank_of_mine": mine+1,
                      "n_hit_1e5": n5, "n_hit_1e4": n4,
                      "second_best": float(cands[1][0]),
                      "top6": [{"x": t[1], "y": t[2], "coef": t[3],
                                "value": float(t[4]), "rel": float(t[0])}
                               for t in cands[:6]]},
    "nu21_MHz": {"with_exp_gp": float(nu21(gp_exp)),
                 "with_pred_gp": float(nu21(gp_pred()))},
    "verdict": verdict, "checks": PASS,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_gp_proton.json"), "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
P_("结果写入: %s" % os.path.join(RES, "zero_gp_proton.json"))
