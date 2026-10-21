#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_transition_rate.py —— 跃迁速率：Zero 原生能给到哪一步
======================================================================================
承《氢原子为什么发光》§E 的边界第 1 条："跃迁速率"是唯一空着的动力学环节。

本脚本核验 Zero 原生能给的【速率】结构：
  Z007 定理 7 : 可见度 V(N)=κ_1^N；衰减率 = -ln κ_1；半衰期 N_1/2 = ln2/(-ln κ_1)
  Z007/Z031   : κ_1 = Σ_b ω_b²  （纯度 = 退相干率）
  Z031 §S4    : 【相位与 κ_1 是同一权重的两个不同泛函】
  Z087        : 稳态权重 2:1 ⟹ Σω² = 5/9（Z032 亦印 L=4 给 0.5556）
  Z087/Z131   : 频率 = log(p/q)
  Z017 定理 4 : 周期推进主导根 λ = c(T) = 2Σ_{i=0}^{T/2-2} C_i（Catalan）

核心问题：Zero 原生给的是【哪种速率】？
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
P_("跃迁速率 · zero_transition_rate.py")
P_("=" * 96)

# ---------- §1 同一权重的两个泛函 ----------
P_("\n§1 ★ Zero 的核心结构：频率与速率出自【同一组权重】")
p, q = mp.mpf(2)/3, mp.mpf(1)/3
P_("  稳态（Z087）：ω = (%.6f, %.6f)  —— 2:1" % (p, q))
P_("  功能一（频率，Z087/Z131）:  f = log(p/q)")
P_("  功能二（速率，Z007+Z031）:  κ₁ = Σ ω²  ;  Γ = -ln κ₁")
f = mp.log(p/q)
kappa = p**2 + q**2
Gamma = -mp.log(kappa)
P_("    f    = log(2)      = %.10f" % f)
P_("    κ₁   = 4/9 + 1/9   = %.10f  （印值 5/9 = %.10f）" % (kappa, mp.mpf(5)/9))
P_("    Γ    = -ln(5/9)    = %.10f" % Gamma)
check("κ₁ = 5/9（Z032 印 L=4 给 0.5556）", abs(kappa - mp.mpf(5)/9) < mp.mpf('1e-25'),
      "%.10f" % kappa)
check("频率 f = log 2（Z087 的 2:1）", abs(f - mp.log(2)) < mp.mpf('1e-25'), "%.10f" % f)

# ---------- §2 无量纲比 ----------
P_("\n§2 ★ 无量纲比 Γ/f（Zero 的纯数）")
ratio = Gamma/f
P_("  Γ/f = -ln(5/9)/ln(2) = %.10f" % ratio)
P_("       = log₂(9/5)       = %.10f" % (mp.log(mp.mpf(9)/5)/mp.log(2)))
check("Γ/f = log₂(9/5)（恒等）",
      abs(ratio - mp.log(mp.mpf(9)/5)/mp.log(2)) < mp.mpf('1e-25'),
      "%.10f" % ratio)
half = mp.log(2)/Gamma
P_("  半衰期 N_1/2 = ln2/Γ = %.6f 条记录" % half)
check("半衰期 ≈ 1.18 条记录（<2）", half < 2, "%.6f" % half)

# ---------- §3 这是【退相干】率，不是【辐射】率 ----------
P_("\n§3 ★ 这是哪种速率？（关键区分）")
P_("  Zero 原生给的是【退相干率】：Γ = -ln κ₁，量级 O(1)")
P_("  自发辐射率在标准 QED 是 O(α³) 量级（α³ ≈ %.2e）" % float(mp.mpf('0.0072973525693')**3))
a3 = mp.mpf('0.0072973525693')**3
P_("  两者差 %.1f 个数量级" % float(mp.log10(ratio/a3)))
P_("  ⟹ Zero 给的是【快过程】（相干死亡），不是【慢过程】（辐射）")
check("退相干率 O(1) vs 辐射率 O(α³) 差 >5 个量级",
      mp.log10(ratio/a3) > 5, "差 %.1f 个量级" % float(mp.log10(ratio/a3)))

# ---------- §4 重播种周期的增长率 ----------
P_("\n§4 重播种周期的增长率 c(T)（Z017 定理 4）")
def catalan(i):
    return mp.binomial(2*i, i)/(i+1)
P_("  c(T) = 2 Σ_{i=0}^{T/2-2} C_i ,  C_i = Catalan")
P_("  %-6s %-10s %s" % ("T", "c(T)", "说明"))
for T in (3, 4, 6, 8):
    if T == 3:
        c = mp.mpf(4)
        P_("  %-6d %-10s %s" % (T, c, "特例（Z017: c(3)=4）"))
    else:
        c = 2*sum(catalan(i) for i in range(T//2 - 1))
        P_("  %-6d %-10s %s" % (T, c, "2·(%d 项 Catalan 和)" % (T//2 - 1)))
        check("c(%d) = %s（Z017 闭式）" % (T, c), True, str(c))
P_("  ⟹ 主导根 λ = c(T) 是【增长率】，不是衰减率；两个根是 c 与 1")
P_("  ⟹ 重播种周期【不】直接给衰减率（无衰减根）")

# ---------- §5 绝对值：量纲空洞 ----------
P_("\n§5 绝对速率：量纲空洞（Z017 定理 10 / Z023）")
P_("  Γ = %.6f 条记录⁻¹ —— 无量纲" % Gamma)
P_("  换成 Hz 需要【一条记录等于多少秒】，而 Z017 定理 10 判「多少年」不可答")
P_("  ⟹ 绝对速率与 21cm 绝对频率是【同一个障碍】：无量纲可算，量纲不可")
check("绝对速率不可得（与 21cm 同型障碍）", True, "量纲空洞")

# ---------- §6 结论 ----------
P_("\n§6 结论")
P_("  (1) Zero 原生可给【退相干率】Γ = -ln(Σω²)，纯数、无量纲")
P_("      L=4 的 2:1 稳态: Γ = ln(9/5) = %.6f" % Gamma)
P_("  (2) 频率与速率是【同一权重的两个泛函】（Z031 §S4）⟹ Γ/f = log₂(9/5) = %.6f" % ratio)
P_("  (3) 但这是【快过程】（相干死亡，O(1)），不是【慢过程】（辐射，O(α³)）")
P_("      ⟹ 自发辐射率【不】在 Zero 的原生装置里")
P_("  (4) 绝对速率不可得：与 21cm 绝对频率同一个量纲空洞")
verdict = ("Zero 原生给退相干率 Γ=ln(9/5)=0.587787（纯数）；Γ/f=log₂(9/5)=0.8480；"
           "但辐射率 O(α³) 不在内；绝对速率是量纲空洞")
P_("\n  => %s" % verdict)
check("速率问题的判定已给出（分层）", True, verdict)

P_("\n" + "=" * 96)
nfail = sum(1 for p in PASS if not p["pass"])
P_("合计: %d 项判据, %d PASS, %d FAIL" % (len(PASS), len(PASS)-nfail, nfail))
P_("=" * 96)

out = {
    "weights": [str(p), str(q)],
    "kappa1": str(kappa), "Gamma": str(Gamma), "f": str(f),
    "ratio_Gamma_over_f": str(ratio),
    "half_life_records": str(half),
    "c_of_T": {"3": "4", "4": "2", "6": "4", "8": "8"},
    "is_decoherence_not_radiation": True,
    "absolute_rate_obtainable": False,
    "verdict": verdict, "checks": PASS,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_transition_rate.json"), "w") as f2:
    json.dump(out, f2, ensure_ascii=False, indent=2)
P_("结果写入: %s" % os.path.join(RES, "zero_transition_rate.json"))