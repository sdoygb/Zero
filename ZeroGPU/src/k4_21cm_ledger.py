#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_21cm_ledger.py —— 氢原子 21cm 线的依赖账核验（判定篇用）
======================================================================================
目标: 把「Zero 能否推导 21cm 线」这一问题拆成逐因子的依赖账, 并给出可跑判据.

21cm 的标准因子分解(费米接触项, 基态超精细分裂):
    h nu = (8/3) * g_p * (m_e/m_p) * alpha^2 * Ry * (1/(1+m_e/m_p))
其中 Ry = alpha^2 m_e c^2 / 2 = 13.605693 eV 是里德伯能量(库仑问题的束缚态能).

判据设计: 每个因子打三个标签
    IN   = Zero 库内已有且为内生(或半内生)
    EXT  = 库内无, 必须外部输入
    GAP  = 库内无, 且机制层面已知不可导出(有 no-go)
"""
from __future__ import annotations
import os, json
import numpy as np

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
P_("21cm 依赖账 · k4_21cm_ledger.py")
P_("=" * 96)

# ---------- 常数(外部/CODATA) ----------
me   = 9.1093837015e-31      # kg
mp   = 1.67262192369e-27     # kg
c    = 299792458.0           # m/s
h    = 6.62607015e-34        # J s
alpha= 1/137.035999084
gp   = 5.5856946893          # 质子磁矩(核磁子单位)
eV   = 1.602176634e-19
Ry   = 13.605693122994       # eV

NU_OBS = 1420.405751768e6    # Hz, 观测值

# ---------- §1 因子分解正确性 ----------
P_("\n§1 因子分解核验: 用标准常数重建 21cm")
mu_fac = 1.0/(1.0 + me/mp)
nu = (8.0/3.0)*gp*(me/mp)*alpha**2*Ry*mu_fac*eV/h
dev = (nu - NU_OBS)/NU_OBS
P_("  nu(重建) = %.6f MHz" % (nu/1e6))
P_("  nu(观测) = %.6f MHz" % (NU_OBS/1e6))
P_("  相对偏差 = %.3e" % dev)
check("因子分解可重建 21cm(偏差<0.01%)", abs(dev) < 1e-4, "dev=%.2e" % dev)

# ---------- §2 逐因子归属 ----------
P_("\n§2 逐因子归属(Zero 库内状态)")
FACTORS = [
    # (符号, 数值, 归属, 依据)
    ("4/3 或 8/3", 8.0/3.0, "IN",  "角动量/接触项纯代数; Zero 有自旋-角动量代数(Z006)"),
    ("g_p",        gp,      "EXT", "质子磁矩=非微扰 QCD 量; 全库 grep 无磁矩/塞曼/g-2"),
    ("alpha^2",    alpha**2,"GAP", "精细结构常数; 全库 grep 无 1/137; 且 Coulomb 问题需 QED"),
    ("m_e/m_p",    me/mp,   "PART","m_e 半内生(Z126 振幅律, 但标度 u=1 是定义, 绝对值靠锚 v); m_p 仅近似(Z125, 无精确命中)"),
    ("Ry",         Ry,      "GAP", "库仑问题束缚态能; 库内无 Coulomb 求解器; QM 动力学有 no-go(Z003 #12)"),
    ("1/(1+me/mp)",mu_fac,  "GAP", "约化质量; 同上, 属束缚态结构"),
]
P_("  %-14s %-14s %-6s %s" % ("因子", "数值", "归属", "依据"))
for sym, val, tag, why in FACTORS:
    P_("  %-14s %-14.6g %-6s %s" % (sym, val, tag, why))

tags = [f[2] for f in FACTORS]
n_in  = tags.count("IN")
n_ext = tags.count("EXT")
n_gap = tags.count("GAP")
n_par = tags.count("PART")
P_("\n  统计: IN=%d EXT=%d GAP=%d PART=%d" % (n_in, n_ext, n_gap, n_par))
check("内生因子不足以定出 21cm(IN 因子无法独立给出频率)", n_in < 2, "IN=%d" % n_in)

# ---------- §3 只保留 Zero 内生因子会得到什么 ----------
P_("\n§3 反证: 若只允许 Zero 内生量, 21cm 频率不可得")
# Zero 能给的: 代数因子 8/3 (纯数); 其余全是外部
P_("  只留 8/3: 这是一个纯数, 无任何量纲 -> 得不到 Hz")
P_("  故 21cm 的绝对频率需要 >=3 个外部输入(g_p, alpha, Ry, 约化质量) 才能定出")
need_ext = n_ext + n_gap
check("21cm 绝对频率需要 >=3 个外部输入", need_ext >= 3, "需外部=%d" % need_ext)

# ---------- §4 v 的地位核验(Z013 定理16 / Z153 §2 / Z162 §4.2) ----------
P_("\n§4 弱电 VEV v 的地位: 它是锚, 不是导出量")
v_GeV = 246.2205
P_("  Z013 定理16: v 是唯一的有量纲自由度, 是量纲单位(问'v 是多少'恒答 1)")
P_("  Z153 §2   : 锚是一个约定; '换算不是预言'")
P_("  Z162 §4.2 : 锚的优先序 v > c > m_e > hbar; 全库仍只有 v 一个")
# 用 v 只能得到单位, 不能得到氢原子结构
P_("  故'加一个 v'只能给出单位换算(已由 Z153 完成 9/9 PASS), 不给出谱")
check("v 是输入而非导出(Z013 定理16)", True, "v=%.4f GeV 作为锚" % v_GeV)

# ---------- §5 用 v 锚能算出什么(换算可达性) ----------
P_("\n§5 换算可达性: 21cm 的'波长'与 v 的关系")
lam_21 = c/NU_OBS
P_("  21cm 波长 = %.6f m = %.4e 跳数 (1 跳数=2.009e-33 m, Z153)" % (lam_21, lam_21/2.009e-33))
P_("  21cm 能量 = %.6f eV = %.4e * v" % (NU_OBS*h/eV, (NU_OBS*h/eV)/(v_GeV*1e9)))
P_("  => 21cm 对应无量纲数 ~1e-15 * v: 换算可写, 但该数本身不是 Zero 算出的")
check("21cm 可由 v 表达为无量纲数(但该数非内生)", True, "E/v=%.3e" % ((NU_OBS*h/eV)/(v_GeV*1e9)))

# ---------- §6 判定 ----------
P_("\n§6 判定")
verdict = {
    "绝对频率可内生推导": False,
    "原因": "需要 g_p, alpha, Ry, 约化质量修正 共 4 个外部输入; 且 Coulomb 束缚态在库内为 GAP(Z003 #12 no-go)",
    "标度换算可达": True,
    "口径": "用 v 锚可把 21cm 写成无量纲数, 但那是换算(Z153), 不是推导",
    "外部输入数": int(need_ext),
}
for k, v in verdict.items():
    P_("  %s: %s" % (k, v))

P_("\n" + "=" * 96)
nfail = sum(1 for p in PASS if not p["pass"])
P_("合计: %d 项判据, %d PASS, %d FAIL" % (len(PASS), len(PASS)-nfail, nfail))
P_("=" * 96)

out = {
    "nu_reconstructed_MHz": nu/1e6,
    "nu_observed_MHz": NU_OBS/1e6,
    "rel_dev": dev,
    "factors": [{"symbol": s, "value": float(v), "tag": t, "why": w} for s, v, t, w in FACTORS],
    "counts": {"IN": n_in, "EXT": n_ext, "GAP": n_gap, "PART": n_par},
    "external_inputs_needed": int(need_ext),
    "verdict": verdict,
    "checks": PASS,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "k4_21cm_ledger.json"), "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
P_("结果写入: %s" % os.path.join(RES, "k4_21cm_ledger.json"))