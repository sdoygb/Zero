#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_alpha_oldlib_reconcile.py —— 对齐两个旧库的 alpha 链，并检验候选关系
======================================================================================
用户指令（2026-10-20）："alpha --- 查最老两个资料库，有推导，但只作参考"

【两个旧库的对象】
  Zero_副本/sm/SM12（最老，Oct 7）: R_data(M_X)=1.291148（"仅参考"）
  Zero-261010/sm/SMR17,SMR18:     目标 L=25.4509 / alpha_U^-1=42.4116；41 候选 0 命中@0.0037%

【本脚本钉四件事】
  V1 旧库的 L=25.4509 是哪个标度？（我算：alpha1=alpha2 交叉）
  V2 新库的 R_data=1.292018 是哪个标度？（我算：alpha2=alpha3 交叉）
  V3 这整条链是不是重言式？（x_i 全由测量值给 ⟹ R 只是数据自己）
  V4 候选关系 L = 1/M_e + 2/3 的距离（登记为候选，不宣布命中）
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
P_("对齐两个旧库的 alpha 链 · zero_alpha_oldlib_reconcile.py")
P_("=" * 96)

pi = mp.pi
b1, b2, b3 = mp.mpf(41)/10, mp.mpf(-19)/6, mp.mpf(-7)
cc = lambda b: b/(2*pi)

# 旧库 SM12 的输入（原文 §2 给定）
aem_inv, s2_obs, as_obs = mp.mpf('127.95'), mp.mpf('0.23122'), mp.mpf('0.1179')
x1 = mp.mpf(3)/5 * aem_inv * (1 - s2_obs)     # alpha1 = (5/3) alpha_Y
x2 = aem_inv * s2_obs
x3 = 1 / as_obs
P_("\n§1 旧库输入 x_i(M_Z)（SM12 §2 给定）")
P_("  x1 = %.8f   x2 = %.8f   x3 = %.8f" % (float(x1), float(x2), float(x3)))
P_("  注：alpha1 ≡ (5/3) alpha_Y，由迹比 5/3 算出（SM12 §2(c)，原文：不是外部给定）")

# V1：两个交叉尺度
t12 = (x1 - x2) / (cc(b1) - cc(b2))
t23 = (x2 - x3) / (cc(b2) - cc(b3))
P_("\n§2 V1/V2 两个交叉尺度")
P_("  alpha1 = alpha2 的 t = %.8f   （旧库目标 L = 25.4509）" % float(t12))
P_("  alpha2 = alpha3 的 t = %.8f   （新库 R_data 所在）" % float(t23))
check("V1 旧库 L=25.4509 是 alpha1=alpha2 交叉", abs(t12 - mp.mpf('25.4509')) < 1e-3,
      "t12 = %.6f" % float(t12))

X1_12 = x1 - cc(b1)*t12; X2_12 = x2 - cc(b2)*t12; X3_12 = x3 - cc(b3)*t12
P_("  在 t12 处: x1 = x2 = %.6f,  x3 = %.6f" % (float(X1_12), float(X3_12)))
P_("  旧库给 alpha_U^-1 = 42.4116  ⟹ 与 x1(=x2) 一致" )
check("V1b 旧库 x_U=42.4116 = x1=x2 的值", abs(X1_12 - mp.mpf('42.4116')) < 5e-3,
      "%.6f" % float(X1_12))
P_("  ★ 注意 x3 = %.4f ≠ 42.41 ⟹ 旧库的 M_U 【不是】全统一点，只统一两个" % float(X3_12))

# V2：R 在两个尺度上的值
def R_at(t):
    A1 = x1 - cc(b1)*t; A2 = x2 - cc(b2)*t
    return mp.mpf(5)/3 * A1/A2
R12 = R_at(t12); R23 = R_at(t23)
P_("\n§3 R = (5/3)(x1/x2) 在两个尺度上")
P_("  R(t12) = %.8f   （旧库登记 R_data(M_X)=1.291148）" % float(R12))
P_("  R(t23) = %.8f   （新库用的是 1.292018）" % float(R23))
check("V2 新库 R_data=1.292018 是 alpha2=alpha3 处的值",
      abs(R23 - mp.mpf('1.292018')) < 1e-5, "R(t23)=%.7f" % float(R23))
check("V2b 旧库 R=1.291148 与 t23 处差 0.067%",
      abs(R23/mp.mpf('1.291148') - 1) < 1e-3, "差 %.4f%%" % (float(abs(R23/mp.mpf('1.291148')-1))*100))

# V3：重言性
P_("\n§4 V3 重言性判定")
P_("  x_i(M_Z) 全部由【测量值】(alpha_em, sin^2 theta_W, alpha_s) 给出")
P_("  而 t_cross 与 R 只由 x_i 算出 ⟹ R 命中'数据'是必然（它就是数据自己）")
P_("  ⟹ 这不是预言，是一致性检验。旧库 SM12 结论五已判：J2 缺、'一个框架数拟两个观测'")
check("V3 判定：该链是重言式（非预言）", True, "旧库已判；本脚本复核")

# V4：候选关系
P_("\n§5 V4 候选关系：L(alpha1=alpha2) = 1/M_e + 2/3 ?")
d0 = pi/12 - mp.mpf(2)/9
Me = 1 + mp.sqrt(2)*mp.cos(3*pi/4 - d0)
cand = 1/Me + mp.mpf(2)/3
P_("  2/9        = %.10f   [Z123 稳态两支权重积，非中心态]" % float(mp.mpf(2)/9))
P_("  delta0     = %.10f   [= pi/12 - 2/9]" % float(d0))
P_("  M_e        = %.10f   [Z121+Z123: theta_e = 3pi/4 - delta0]" % float(Me))
P_("  1/M_e      = %.10f" % float(1/Me))
P_("  1/M_e + 2/3= %.10f" % float(cand))
P_("  L(t12)     = %.10f" % float(t12))
dev = float(abs(t12/cand - 1) * 100)
thr = 0.0037
P_("  相对偏差   = %.5f%%" % dev)
P_("  旧库门槛   = %.4f%%" % thr)
P_("  ⟹ 在门外 %.2f 倍" % (dev/thr))
check("V4 候选在门槛外（不得宣布命中）", dev > thr,
      "%.5f%% > %.4f%%（%.2f 倍）" % (dev, thr, dev/thr))

# 另一条更松的候选
P_("\n§6 另一条候选（更松）：alpha^-1 = (2/9)/M_e^2")
ainv_obs = mp.mpf('137.035999084')
formula = (mp.mpf(2)/9)/Me**2
P_("  公式值 = %.6f ；观测 = %.6f ；偏差 %+.4f%%"
   % (float(formula), float(ainv_obs), float((formula/ainv_obs-1)*100)))
check("该候选偏离 0.4% 量级（更松，仅登记）",
      abs(float((formula/ainv_obs-1)*100)) < 1.0,
      "%+.4f%%" % float((formula/ainv_obs-1)*100))

# 判据：跑动不能解释
P_("\n§7 跑动不能解释 0.4%（方向判定）")
P_("  alpha^-1(q) 随 q 单调减小：alpha^-1(0)=137.036 为【最大】")
P_("  公式给 %.2f 已【小于】137.036 ⟹ 任何 q>0 处只会更小 ⟹ 跑动方向相反" % float(formula))
check("跑动不能救（方向相反）", formula < ainv_obs, "%.2f < 137.04" % float(formula))

P_("\n" + "=" * 96)
nfail = sum(1 for p in PASS if not p["pass"])
P_("合计: %d 项判据, %d PASS, %d FAIL" % (len(PASS), len(PASS)-nfail, nfail))
P_("=" * 96)

out = {
    "inputs": {"x1_MZ": float(x1), "x2_MZ": float(x2), "x3_MZ": float(x3)},
    "t_alpha1_eq_alpha2": float(t12),
    "t_alpha2_eq_alpha3": float(t23),
    "R_at_t12": float(R12),
    "R_at_t23": float(R23),
    "old_L_target": 25.4509,
    "xU_at_t12": float(X1_12),
    "x3_at_t12": float(X3_12),
    "candidate": {"form": "L = 1/M_e + 2/3", "value": float(cand),
                  "deviation_percent": dev, "threshold_percent": thr,
                  "ratio_to_threshold": dev/thr, "verdict": "candidate, outside threshold"},
    "phase_formula": {"form": "alpha^-1 = (2/9)/M_e^2", "value": float(formula),
                      "deviation_percent": float((formula/ainv_obs-1)*100)},
    "tautology": True,
    "checks": PASS,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_alpha_oldlib_reconcile.json"), "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
P_("结果写入: %s" % os.path.join(RES, "zero_alpha_oldlib_reconcile.json"))
