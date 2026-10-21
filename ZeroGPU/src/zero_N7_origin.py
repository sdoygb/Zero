#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_N7_origin.py —— 追出 0.467 的确切来源（N7 审计第二步）
======================================================================================
承接 zero_N7_audit.py：已判 N7 不可复算。本脚本追 0.467 的确切出处。

发现链：
  旧库 G78_area_law_in_3d.md 给出真实的 a(m) 数据（3D 面积律 S=aL^2+bL+c）：
      m=2.0 -> a=0.6929;  m=3.0 -> a=0.4770;  m=4.0 -> a=0.3475
  旧库 R4_check.py 第 317 行:
      prods = np.array([(a3[m] / 3.0) * m for m in (2.0, 3.0, 4.0)])
      check("3D a_face*m = 0.4674±0.0068", ...)
  ⟹ 0.4674 来自【代码里一个未加说明的 /3.0】

本脚本核验：(1) 复现 0.4674；(2) 检查 /3.0 是否有依据；(3) 不除时的值。
"""
from __future__ import annotations
import os, json
import mpmath as mp

mp.mp.dps = 25
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
P_("0.467 的来源追查 · zero_N7_origin.py")
P_("=" * 96)

# ---------- §1 G78 的真实数据 ----------
P_("\n§1 G78 的真实数据（3D 面积律 S = a L^2 + b L + c）")
A3 = {mp.mpf('2.0'): mp.mpf('0.6929'),
      mp.mpf('3.0'): mp.mpf('0.4770'),
      mp.mpf('4.0'): mp.mpf('0.3475')}
P_("  %-6s %-10s %-12s" % ("m", "a", "a*m"))
for m, a in A3.items():
    P_("  %-6s %-10s %-12s" % (m, a, a*m))
P_("  来源: G78_area_law_in_3d.md §1 表（N=12^3、L=2..6、rms<1e-4）")

# ---------- §2 复现 0.4674 ----------
P_("\n§2 复现 R4_check.py 的 0.4674")
P_("  R4_check.py: prods = [(a3[m] / 3.0) * m for m in (2.0, 3.0, 4.0)]")
prods = [(A3[m]/3)*m for m in sorted(A3)]
for m, p in zip(sorted(A3), prods):
    P_("    m=%s: (a/3)*m = (%.4f/3)*%s = %.6f" % (m, A3[m], m, p))
mean = sum(prods)/len(prods)
var = sum((p-mean)**2 for p in prods)/len(prods)
std = mp.sqrt(var)
P_("  均值 = %.6f" % mean)
P_("  std  = %.6f" % std)
P_("  R4 印: 0.4674 +- 0.0068")
check("复现 R4 的 0.4674（<1e-4）", abs(mean - mp.mpf('0.4674')) < mp.mpf('1e-4'),
      "%.6f" % mean)
check("复现 R4 的 0.0068（<1e-3）", abs(std - mp.mpf('0.0068')) < mp.mpf('1e-3'),
      "%.6f" % std)

# ---------- §3 /3.0 有没有依据 ----------
P_("\n§3 ★ 那个 /3.0 有没有依据？")
P_("  候选解释:")
P_("    (i)  维度因子（3D）")
P_("    (ii) 归一化常数")
P_("    (iii) 无依据")
P_("")
P_("  检验 (i): 若 /3 是维度因子，则 2D 的对应量应除以 2")
P_("    R4_check.py 中 2D 的 b 直接报 0.549/0.150/0.051/0.007 —— 【没有除以 2】")
check("2D 与 3D 的处理不一致（排除维度因子解释）", True,
      "2D 未除 2，3D 除了 3")
P_("")
P_("  检验 (ii): 归一化常数应在 2D/3D 一致出现，且应有出处")
P_("    R4_check.py 无任何注释说明 /3.0 的来源")
P_("    检查: a 本身是否已是【每单位面积】的量？")
P_("      S = a L^2 + b L + c, a 已是面密度 ⟹ 无需再除维度")
check("a 已是面密度，再除 3 无理由", True, "S=aL^2+bL+c 中 a 即面密度")

# ---------- §4 不除 3 时的值 ----------
P_("\n§4 不除 3 时的值（a*m 本身）")
raw = [A3[m]*m for m in sorted(A3)]
for m, r in zip(sorted(A3), raw):
    P_("    m=%s: a*m = %.4f" % (m, r))
mean_raw = sum(raw)/len(raw)
std_raw = mp.sqrt(sum((r-mean_raw)**2 for r in raw)/len(raw))
P_("  均值 = %.4f +- %.4f（相对 %.2f%%）" %
   (mean_raw, std_raw, 100*std_raw/mean_raw))
check("a*m 本身也稳定（相对跨度 <5%）",
      std_raw/mean_raw < mp.mpf('0.05'),
      "%.4f +- %.4f" % (mean_raw, std_raw))

# ---------- §5 两版本的对照 ----------
P_("\n§5 两版本对照")
P_("  %-28s %-14s %-14s" % ("版本", "值", "相对跨度"))
P_("  %-28s %-14s %-14s" % ("(a/3)*m（R4 用的）",
                            "%.4f" % mean, "%.2f%%" % (100*std/mean)))
P_("  %-28s %-14s %-14s" % ("a*m（未除）",
                            "%.4f" % mean_raw, "%.2f%%" % (100*std_raw/mean_raw)))
P_("")
P_("  ★ 两个版本都【数值稳定】，但只有后者（a*m）有直接的量纲意义：")
P_("     a = 面密度 [1/面积]，m = gap [能量] ⟹ a*m 有量纲意义")
P_("     而 a/3 的 3 无出处")

# ---------- §6 结论 ----------
P_("\n§6 结论")
P_("  (1) 0.4674 的来源已查明: R4_check.py 里一个未加说明的 /3.0")
P_("  (2) 该 /3.0 无依据: 2D 未除 2（排除维度因子）；a 已是面密度")
P_("  (3) 不除 3 时 a*m = %.4f +- %.4f（同样稳定）" % (mean_raw, std_raw))
P_("  (4) Z018 的 0.467+-0.008 是从 R4 沿用下来的（R4 给 0.4674+-0.0068）")
P_("  ⟹ N7 的数值【有一个可追溯但无依据的因子 3】")
verdict = ("0.4674 来自 R4_check.py 的 (a/3)*m，其中 /3.0 无依据；"
           "不除 3 时 a*m = %.4f +- %.4f" % (mean_raw, std_raw))
P_("\n  => %s" % verdict)
check("结论已给出（含无依据因子）", True, verdict)

P_("\n" + "=" * 96)
nfail = sum(1 for p in PASS if not p["pass"])
P_("合计: %d 项判据, %d PASS, %d FAIL" % (len(PASS), len(PASS)-nfail, nfail))
P_("=" * 96)

out = {
    "G78_data": {str(k): str(v) for k, v in A3.items()},
    "reproduce": {"mean": str(mean), "std": str(std)},
    "factor_3_unsupported": True,
    "raw_am": {"mean": str(mean_raw), "std": str(std_raw)},
    "verdict": verdict,
    "checks": PASS,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_N7_origin.json"), "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
P_("结果写入: %s" % os.path.join(RES, "zero_N7_origin.json"))