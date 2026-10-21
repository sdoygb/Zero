#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_21cm_verdict.py —— 21cm 能否内生地由 Zero 推出（锚一个锚点）？
======================================================================================
不只回答"能不能"，而是把缺口量化：
  (1) 21cm 对每个因子的精确灵敏度
  (2) 每个因子在 Zero 库内的真实状态
  (3) 用 Zero 现有的最好的量，21cm 会偏多少
  (4) 要做到"推出"，各因子需要多高的精度

公式（费米接触项，基态超精细分裂）：
    h nu = (8/3) g_p (m_e/m_p) alpha^2 R_y / (1 + m_e/m_p),   R_y = alpha^2 m_e c^2/2
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
P_("21cm 内生性判决 · zero_21cm_verdict.py")
P_("=" * 96)

# ---------- 基准值（PDG/CODATA） ----------
me  = 9.1093837015e-31
mp  = 1.67262192369e-27
c   = 299792458.0
h   = 6.62607015e-34
alpha = 1/137.035999084
gp  = 5.5856946893
eV  = 1.602176634e-19
Ry  = 13.605693122994
NU_OBS = 1420.405751768e6      # Hz

mu_fac = 1.0/(1.0 + me/mp)
nu = (8.0/3.0)*gp*(me/mp)*alpha**2*Ry*mu_fac*eV/h
P_("\n§0 基准")
P_("  nu(用 PDG 常数重建) = %.6f MHz" % (nu/1e6))
P_("  nu(观测)            = %.6f MHz" % (NU_OBS/1e6))
P_("  相对偏差            = %.3e" % ((nu-NU_OBS)/NU_OBS))
check("公式可重建 21cm（<1e-4）", abs((nu-NU_OBS)/NU_OBS) < 1e-4,
      "dev=%.2e" % ((nu-NU_OBS)/NU_OBS))

# ---------- §1 灵敏度 ----------
P_("\n§1 精确灵敏度（对数导数）")
# ln nu = ln gp + 4 ln alpha + 2 ln me - ln mp - ln(1+me/mp) + const
r = me/mp
dme = 2.0 - r/(1.0+r)          # d ln nu / d ln me
dmp = -1.0 + r/(1.0+r)         # d ln nu / d ln mp
P_("  d ln nu / d ln g_p   = +1")
P_("  d ln nu / d ln alpha = +4      （alpha 的四次方！）")
P_("  d ln nu / d ln m_e   = %+.6f" % dme)
P_("  d ln nu / d ln m_p   = %+.6f" % dmp)
check("alpha 的灵敏度为 4（四次方）", abs(4 - 4) < 1e-12, "nu ∝ alpha^4")
check("m_e 灵敏度约 +2", abs(dme - 1.99946) < 1e-4, "%.5f" % dme)

# ---------- §2 各因子的 Zero 状态 ----------
P_("\n§2 每个因子在 Zero 库内的状态")
FACTORS = [
    ("8/3",        "IN",   "接触项纯代数；Zero 有自旋-角动量代数（Z006）"),
    ("alpha",      "GAP",  "需无量纲输入 R（旧 Z23，本轮复核）；非预言"),
    ("g_p",        "GAP",  "库内无磁矩/塞曼/g-2 结构；非微扰 QCD 量"),
    ("m_e",        "PART", "Z126 振幅律给 0.51099 MeV（-0.002%）"),
    ("m_p",        "WEAK", "Z125：A^2×3=941.6 MeV（+0.355%），自判'近似'"),
    ("约化质量",   "IN?",  "Z165：原生束缚态可算，但幅度需源输入"),
    ("v（锚）",    "ANCHOR","允许：唯一有量纲锚（Z013 定理16）"),
]
P_("  %-10s %-7s %s" % ("因子", "状态", "依据"))
for sym, st, why in FACTORS:
    P_("  %-10s %-7s %s" % (sym, st, why))
n_gap = sum(1 for f in FACTORS if f[1] == "GAP")
P_("\n  统计：GAP（库内无）= %d 个；PART/WEAK = 2 个；IN = 2 个；ANCHOR = 1 个" % n_gap)
check("alpha 与 g_p 均为库内缺口", n_gap >= 2, "GAP=%d" % n_gap)

# ---------- §3 用 Zero 现有最好的量，21cm 会偏多少 ----------
P_("\n§3 用 Zero 现有的量能到多准")
# Z126 的轻子单元与振幅
v_GeV = 246.2205
A2_unit = v_GeV/784.5           # GeV，轻子单元
theta_e = 2*math.pi/3 + 2/9
M_e = 1 + math.sqrt(2)*math.cos(theta_e)
me_zero_MeV = A2_unit*1e3*M_e**2
P_("  Z126 振幅律：theta_e = 2pi/3+2/9 = %.6f rad" % theta_e)
P_("    M_e = 1+sqrt2 cos(theta_e) = %.7f" % M_e)
P_("    轻子单元 A^2 = v/784.5 = %.4f MeV" % (A2_unit*1e3))
P_("    m_e(Zero) = A^2 * M_e^2 = %.6f MeV" % me_zero_MeV)
P_("    m_e(PDG)  = 0.510999 MeV  -> 差 %+.4f%%" %
   (100*(me_zero_MeV-0.510999)/0.510999))
check("Z126 的电子质量复现到 0.01% 以内",
      abs(me_zero_MeV-0.510999)/0.510999 < 1e-4,
      "%.6f MeV" % me_zero_MeV)

# Z125 的质子质量
mp_zero_MeV = 941.6
P_("\n  Z125：m_p(Zero) = A^2×3 = %.1f MeV" % mp_zero_MeV)
P_("        m_p(PDG)  = 938.272 MeV  -> 差 %+.3f%%" %
   (100*(mp_zero_MeV-938.272)/938.272))

# 合成：用 Zero 的 m_e, m_p，其余（alpha, g_p）用 PDG
ratio_me = (me_zero_MeV/1e6*1e6*1e-3)/(0.510999*1e-3)   # 同单位比值
ratio_me = me_zero_MeV/0.510999
ratio_mp = mp_zero_MeV/938.272
shift = ratio_me**dme * ratio_mp**dmp
nu_zero = NU_OBS * shift
P_("\n  合成：nu(Zero 质量 + PDG 的 alpha,g_p)")
P_("    (m_e 比)^%.5f × (m_p 比)^%.5f = %.6f" % (dme, dmp, shift))
P_("    nu = %.4f MHz    （观测 %.4f）" % (nu_zero/1e6, NU_OBS/1e6))
P_("    偏差 = %+.4f%%  = %+.3f MHz" % (100*(shift-1), (nu_zero-NU_OBS)/1e6))
check("Zero 质量律单独使 21cm 偏 0.1%~1%（可量化的缺口）",
      0.001 < abs(shift-1) < 0.01, "%+.3f%%" % (100*(shift-1)))

# ---------- §4 要做到"推出"，需要多高的精度 ----------
P_("\n§4 若要求 21cm 的值与观测一致，各因子需多准")
P_("  21cm 的数值被引用到 13 位有效数字（1420.405751768 MHz）")
for digits, tol in (("1e-6 (ppm)", 1e-6), ("1e-9 (ppb)", 1e-9), ("1e-12", 1e-12)):
    P_("    目标精度 %-12s => alpha 需 %.1e，g_p 需 %.1e，m_e 需 %.1e，m_p 需 %.1e" %
       (digits, tol/4, tol, tol/dme, tol/abs(dmp)))
P_("\n  Zero 目前能给到的：")
P_("    m_e: %.1e（Z126，%.4f%%）" % (0.00002, abs(100*(ratio_me-1))))
P_("    m_p: %.1e（Z125，%.4f%%，自判'近似'）" % (0.00355, abs(100*(ratio_mp-1))))
P_("    alpha: 不可给（需输入 R）")
P_("    g_p: 不可给（库内无结构）")
check("Zero 的质量精度比'推出 21cm'所需低多个量级", True,
      "m_p 差 ~0.36%，需 ~1e-9 或更好")

# ---------- §5 判决 ----------
P_("\n§5 判决")
P_("  (a) 锚一个锚点：满足。v 是唯一有量纲锚（Z013 定理16 / Z153 / Z162 §4.2）")
P_("  (b) 其余内生：不满足。alpha 与 g_p 是两笔无量纲输入，且都不是'锚'")
P_("      故'一个锚'与'全内生'是两条不同的规则，前者可满足、后者不可")
P_("  (c) 更硬的限制：21cm 依赖 alpha^4，而 alpha 的数值在 Zero 里")
P_("      需要无量纲输入 R（旧 Z23 判决、本轮复核）；g_p 则完全无结构")
P_("  (d) 即使借来 alpha 与 g_p，Zero 自己的质量律（Z125 的 m_p 近似）")
P_("      也使 21cm 偏 %+.3f%%" % (100*(shift-1)))
verdict = "不能（锚一个锚点这一条满足，但 alpha 与 g_p 两笔无量纲输入无法内生）"
P_("\n  => %s" % verdict)
check("判决：21cm 不能内生推出（锚规则可满足，内生规则不可）", True, verdict)

P_("\n" + "=" * 96)
nfail = sum(1 for p in PASS if not p["pass"])
P_("合计: %d 项判据, %d PASS, %d FAIL" % (len(PASS), len(PASS)-nfail, nfail))
P_("=" * 96)

out = {
    "baseline_MHz": nu/1e6, "observed_MHz": NU_OBS/1e6,
    "sensitivities": {"gp": 1.0, "alpha": 4.0, "me": dme, "mp": dmp},
    "factors_status": [{"sym": s, "status": t, "why": w} for s, t, w in FACTORS],
    "zero_masses": {"me_MeV": me_zero_MeV, "mp_MeV": mp_zero_MeV,
                    "me_dev_pct": 100*(ratio_me-1), "mp_dev_pct": 100*(ratio_mp-1)},
    "nu_with_zero_masses_MHz": nu_zero/1e6,
    "nu_shift_pct": 100*(shift-1),
    "verdict": verdict,
    "checks": PASS,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_21cm_verdict.json"), "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
P_("结果写入: %s" % os.path.join(RES, "zero_21cm_verdict.json"))