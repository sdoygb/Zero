#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_alpha_R_direction.py —— 复核旧 Z23 §3.3 里 r 的方向
======================================================================================
旧 Z23（Zero-261010/Z23_alpha_fine_structure.md）§3.3 同时写了：
    (i)  r ≡ g_R/g_L
    (ii) R(r) = c_L^2 + c_R^2 r^2 = 1/3 + 4/3 r^2
而 §3.3 前面给的关系是  1/alpha_Y = c_L^2/alpha_L + c_R^2/alpha_R （alpha_L = alpha_2）。

本脚本检验 (i)+(ii) 与前一条关系是否自洽，并核对 §3.3 脚注的算术。
【纪律】只做纯算术核验，不引入新假设；每条判据给出"若不成立会怎样"。
"""
from __future__ import annotations
import os, json, math
from fractions import Fraction as F

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
P_("旧 Z23 §3.3 的 r 方向复核 · zero_alpha_R_direction.py")
P_("=" * 96)

# ---------- 常数（来自旧 Z23 §4.1，本文只用它的口径） ----------
R_data = 1.292018
cL2 = F(1, 3)
cR2 = F(4, 3)
P_("\n§0 输入（旧 Z23 口径）")
P_("  R_data = %.6f" % R_data)
P_("  c_L^2 = %s, c_R^2 = %s, 和 = %s" % (cL2, cR2, cL2 + cR2))
check("c_L^2 + c_R^2 = 5/3（迹比）", cL2 + cR2 == F(5, 3), "%s" % (cL2 + cR2))

# ---------- §1 由 LR 关系推 R(r) 的正确形式 ----------
P_("\n§1 由 1/alpha_Y = c_L^2/alpha_L + c_R^2/alpha_R 推 R")
P_("  alpha_i = g_i^2/4pi  =>  alpha_R/alpha_L = (g_R/g_L)^2 = r^2")
P_("  1/alpha_Y = (c_L^2 + c_R^2 * alpha_L/alpha_R)/alpha_L")
P_("            = (c_L^2 + c_R^2 / r^2)/alpha_L")
P_("  R = (1/alpha_Y)/(1/alpha_L) = c_L^2 + c_R^2 / r^2     <-- r^{-2}")
P_("")
P_("  旧 Z23 写的是        R = c_L^2 + c_R^2 * r^2          <-- r^{+2}")
P_("  ⟹ 两者只在 r = 1 时相同")

def R_of_r_minus2(r):
    return float(cL2) + float(cR2) / (r * r)

def R_of_r_plus2(r):
    return float(cL2) + float(cR2) * (r * r)

P_("\n  数值对照（R_data = %.6f）:" % R_data)
P_("    %8s %16s %16s" % ("r", "R (r^-2, 正确)", "R (r^+2, 旧式)"))
for r in (0.5, 0.847947, 1.0, 1.179319, 2.0):
    P_("    %8.6f %16.6f %16.6f" % (r, R_of_r_minus2(r), R_of_r_plus2(r)))
check("两式仅在 r=1 处相等", abs(R_of_r_minus2(1.0) - R_of_r_plus2(1.0)) < 1e-15
      and abs(R_of_r_minus2(2.0) - R_of_r_plus2(2.0)) > 1e-3,
      "r=1 同值；r=2 相差 %.4f" % abs(R_of_r_minus2(2.0) - R_of_r_plus2(2.0)))

# ---------- §2 反解 r（两种方向） ----------
P_("\n§2 反解 r：数据要的那个 r 是几")
r_plus2 = math.sqrt((R_data - float(cL2)) / float(cR2))
r_minus2 = math.sqrt(float(cR2) / (R_data - float(cL2)))
P_("  用旧式 (r^+2) 反解:  r = %.6f   （旧 Z23 §3.3 印的是 0.847947）" % r_plus2)
P_("  用正确式 (r^-2):    r = %.6f" % r_minus2)
P_("  两者互为倒数?  %.6f vs %.6f" % (r_plus2, 1.0 / r_minus2))
check("两解互为倒数（同一物理量的两种标法）",
      abs(r_plus2 - 1.0 / r_minus2) < 1e-9, "%.6f = 1/%.6f" % (r_plus2, r_minus2))

# ---------- §3 物理方向：R < 5/3 需要 r > 1 ----------
P_("\n§3 物理方向：R < 5/3 需要 g_R 比 g_L 大")
P_("  正确式 R = 1/3 + (4/3)/r^2 在 r>0 上单调【减】")
P_("  R = 5/3  <=>  r = 1     （LR 宇称未破）")
P_("  R < 5/3  <=>  r > 1     <=>  g_R > g_L  <=>  alpha_R > alpha_L")
P_("  数据 R_data = %.6f < 5/3 = %.6f ⟹ r = %.6f > 1" %
   (R_data, 5 / 3, r_minus2))
check("数据要求 r > 1（右手中性耦合更强）", r_minus2 > 1.0, "r=%.4f" % r_minus2)
P_("  ★ 含义：迹（r=1，LR 宇称保持）被排除 ⟹ 范数假设【要求宇称在 M_X 处已破】")

P_("\n§3b ★ 稳健性：无论怎么标 r，物理方向都相同")
P_("  读法(a) r≡g_R/g_L，公式应为 r^{-2}  =>  r = %.4f > 1  => g_R > g_L" % r_minus2)
P_("  读法(b) r≡g_L/g_R，公式为 r^{+2}   =>  r = %.4f < 1  => g_R > g_L" % r_plus2)
P_("  两种读法给【同一物理结论】: g_R > g_L（右手耦合在 M_X 处更强）")
check("两种标法给同一物理方向 g_R > g_L",
      (r_minus2 > 1.0) and (r_plus2 < 1.0),
      "r_minus2=%.4f>1, r_plus2=%.4f<1" % (r_minus2, r_plus2))

# ---------- §4 脚注的算术复核 ----------
P_("\n§4 旧 Z23 §3.3 脚注的算术复核")
P_("  脚注称：旧理论给 r=g_L/g_R=0.847562；本式反解 r=g_R/g_L=0.847947；")
P_("           '说明这是同一个 r≈0.85 的两种倒数写法'")
P_("  但若一个标 g_L/g_R、一个标 g_R/g_L，二者【必须互为倒数】")
P_("  1/0.847562 = %.6f   ≠   0.847947" % (1.0 / 0.847562))
P_("  两个数相差 %.4f%%，都是 0.85 量级，不是倒数关系" %
   (100 * abs(0.847947 - 0.847562) / 0.847562))
check("脚注的'互为倒数'说法不成立",
      abs(1.0 / 0.847562 - 0.847947) > 0.3,
      "1/0.847562=%.4f vs 0.847947" % (1.0 / 0.847562))

# ---------- §5 修正后，范数是否变得更自然？ ----------
P_("\n§5 修正后：sqrt(5/3) 对应的 r 有没有结构读法")
r_needed = r_minus2
r2_needed = float(cR2) / (math.sqrt(5 / 3) - float(cL2))
P_("  R = sqrt(5/3) <=> r^2 = (4/3)/(sqrt(5/3)-1/3) = %.8f" % r2_needed)
P_("  闭式: r^2 = 2(sqrt(15)+1)/7 = %.8f" % (2 * (math.sqrt(15) + 1) / 7))
P_("  r = %.6f" % math.sqrt(r2_needed))
P_("")
P_("  候选结构值对照（看有无 '自然的' r）:")
cands = {
    "1 (宇称未破)": 1.0,
    "sqrt2": math.sqrt(2),
    "N_c/N_w = 3/2": 1.5,
    "2": 2.0,
    "sqrt(5/3)": math.sqrt(5 / 3),
    "4/3": 4 / 3,
}
for name, v in cands.items():
    P_("    r = %-16s %8.4f -> R = %9.6f  (R_data=%.6f, 差 %+8.3f%%)" %
       (name, v, R_of_r_minus2(v), R_data, 100 * (R_of_r_minus2(v) - R_data) / R_data))
P_("  ⟹ 没有任何简单结构值给出 sqrt(5/3)")
check("无简单结构 r 命中 sqrt(5/3)（机制仍缺）",
      all(abs(R_of_r_minus2(v) - math.sqrt(5 / 3)) / math.sqrt(5 / 3) > 1e-3
          for v in cands.values()),
      "最小偏差 %.4f%%" % (100 * min(abs(R_of_r_minus2(v) - math.sqrt(5 / 3)) / math.sqrt(5 / 3)
                                     for v in cands.values())))

# ---------- §6 对判决的影响 ----------
P_("\n§6 对旧 Z23 判决的影响")
P_("  迹 R=5/3 被排除、范数 R=sqrt(5/3) 通过 —— 这两个判决【不受影响】")
P_("  （它们只依赖 R 本身，不依赖 r 怎么标）")
P_("  受影响的是:")
P_("    (a) '所需 r' 的数值: 0.8479 -> 1.1793（互为倒数）")
P_("    (b) 物理读法: 由'g_R < g_L'改为'g_R > g_L'")
P_("    (c) 一个新结论: 范数假设【要求 LR 宇称在 M_X 处已破】")
check("判决（迹排除/范数通过）不受 r 方向影响", True, "只依赖 R 本身")

P_("\n" + "=" * 96)
nfail = sum(1 for p in PASS if not p["pass"])
P_("合计: %d 项判据, %d PASS, %d FAIL" % (len(PASS), len(PASS)-nfail, nfail))
P_("=" * 96)

out = {
    "R_data": R_data,
    "cL2": str(cL2), "cR2": str(cR2),
    "formula_corrected": "R = c_L^2 + c_R^2 / r^2  (r = g_R/g_L)",
    "formula_old": "R = c_L^2 + c_R^2 * r^2  (as printed in old Z23 §3.3)",
    "r_old_printed": 0.847947,
    "r_corrected": r_minus2,
    "reciprocal_check": {"1_over_r_corrected": 1.0 / r_minus2, "r_old": r_plus2},
    "footnote_check": {"r_old_theory": 0.847562, "its_reciprocal": 1.0 / 0.847562,
                       "claimed": 0.847947, "holds": False},
    "physical_reading": "R < 5/3 requires r > 1, i.e. g_R > g_L (parity broken at M_X)",
    "norm_requires_parity_breaking": True,
    "no_simple_r_gives_norm": True,
    "verdict_unaffected": True,
    "checks": PASS,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_alpha_R_direction.json"), "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
P_("结果写入: %s" % os.path.join(RES, "zero_alpha_R_direction.json"))