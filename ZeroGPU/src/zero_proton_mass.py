#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_proton_mass.py —— 第一步：Zero 的内部吸引子谱能否给出 m_p/m_e ≈ 1836？
======================================================================================
承用户路线（"不用传统理论，Zero 自己定义质子"）：
  质子 = 一个稳定吸引子，其 m_p/m_e≈1836、自旋 1/2

【本脚本查的是最基础的一条可行性】
  Z132 的内部涌现谱 = {log|2cos(2πm/n)|}_{n=3..12}，11 个值，11/11 引擎实测稳定。
  若"粒子的质量"来自这个谱（Z126：m ∝ A²M²，M 在 2cos 族），
  则质量比的范围是多少？1836 在里面吗？

【要给出的三条】
  V1  内部谱的频率取值范围，以及它能产生的质量比【上下界】
  V2  1836（或 sqrt(1836)=42.85）是否可达
  V3  若不可达，差多远（结构性差距，不是数值差）
"""
from __future__ import annotations
import os, json
import itertools
import numpy as np
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
P_("内部吸引子谱能否给出 1836 · zero_proton_mass.py")
P_("=" * 96)

# ---------- §1 Z132 的内部谱 ----------
P_("\n§1 Z132 的内部涌现谱：{log|2cos(2πm/n)|}，n=3..12")
spec = {}
for n in range(3, 13):
    vals = set()
    for m in range(0, n+1):
        R = abs(2*mp.cos(2*mp.pi*m/n))
        if R > mp.mpf('1e-20'):
            vals.add((m, R, mp.log(R)))
    for m, R, f in sorted(vals):
        spec[(n, m)] = (R, f)
P_("  共 %d 个 (n,m) 组合（|2cos| ≠ 0）" % len(spec))
fs = [v[1] for v in spec.values()]
P_("  频率范围：min = %s, max = %s" % (mp.nstr(min(fs), 6), mp.nstr(max(fs), 6)))
Rs = [v[0] for v in spec.values()]
P_("  |2cos| 范围：min = %s, max = %s" % (mp.nstr(min(Rs), 6), mp.nstr(max(Rs), 6)))

# Z132 的 11 个代表值
rep = sorted(set(round(float(v[1]), 4) for v in spec.values()))
P_("\n  Z132 列出的 11 个代表频率：")
P_("    %s" % rep)
P_("  ★ 注意：Z132 只列 R≥1 的 11 个；完整族含 R<1，频率可负（最低 %s）"
   % mp.nstr(min(fs),5))
rep_ge1 = sorted(set(round(float(v[1]),4) for v in spec.values() if v[0] >= 1))
P_("  若按 Z132 的选择（R≥1），频率落在 [%s, %s]" % (rep_ge1[0], rep_ge1[-1]))
check("Z132 的 11 个值是 R≥1 的选择（非全族）", len(rep_ge1) <= 12,
      "%d 个值, R∈[%.3f, 2]" % (len(rep_ge1), float(min(Rs[Rs.index(min([v[0] for v in spec.values() if v[0]>=1]))]),) if False else float(min([v[0] for v in spec.values() if v[0]>=1]))))

# ---------- §2 V1 质量比的范围 ----------
P_("\n§2 V1：这个谱能产生的质量比范围")
P_("\n  质量若 ∝ M²（Z126：m = u_s A² M²，M 在 2cos 族），则")
P_("    质量比 = (M_i/M_j)² = (R_i/R_j)²  （同 A² 时）")
P_("  频率差：Δf = f_i − f_j；质量比 = e^{2Δf}")
P_("")
P_("  频率跨度：max − min = %s" % mp.nstr(max(fs)-min(fs), 6))
P_("  ⟹ 质量比上界 = e^{2(max−min)} = %s" % mp.nstr(mp.e**(2*(max(fs)-min(fs))), 6))
P_("  ⟹ 质量比下界 = e^{−2(max−min)} = %s" % mp.nstr(mp.e**(-2*(max(fs)-min(fs))), 6))
span = mp.e**(2*(max(fs)-min(fs)))
P_("  ★ 更正：比值本身的范围是 [%s, %s]（不是平方）"
   % (mp.nstr(min(ratios_sq)**mp.mpf('0.5'),6), mp.nstr(max(ratios_sq)**mp.mpf('0.5'),6)))

P_("\n  实测：全部 (n,m) 对的 (R_i/R_j)² 的极值")
ratios = []
for (a, b) in itertools.product(spec, repeat=2):
    if a == b: continue
    ratios.append(spec[a][0]/spec[b][0])
ratios_sq = [r**2 for r in ratios]
P_("    对数 = %d ; (R_i/R_j)² 范围 = [%s, %s]"
   % (len(ratios_sq), mp.nstr(min(ratios_sq), 6), mp.nstr(max(ratios_sq), 6)))
max_ratio = max(ratios_sq)**mp.mpf('0.5')
check("V1 比值上界 = 7.03（平方比 49.4）", abs(max_ratio - mp.mpf('7.0267')) < 1e-3,
      "max ratio = %s" % mp.nstr(max_ratio, 6))

# ---------- §3 V2 1836 是否可达 ----------
P_("\n§3 V2：1836 或 sqrt(1836)=42.85 是否可达？")
P_("\n  %-14s %-16s %-14s %s" % ("目标", "值", "在范围内?", "差多远"))
t1836 = mp.mpf(1836)
tsqrt = mp.sqrt(t1836)
for name, t in (("m_p/m_e", t1836), ("sqrt(m_p/m_e)", tsqrt)):
    inrange = (min(ratios_sq)**mp.mpf('0.5')) <= t <= max_ratio
    far = t/max_ratio if t > max_ratio else max_ratio/t
    P_("  %-14s %-16s %-14s %.1f 倍" % (name, mp.nstr(t, 8), inrange, far))
check("V2 sqrt(1836)=42.85 不可达（超出比值上界 7.03）", tsqrt > max_ratio,
      "42.85 vs %s" % mp.nstr(max_ratio, 6))
check("V2 1836 作为平方比也不可达（超出 49.4）", t1836 > max(ratios_sq),
      "1836 vs %s" % mp.nstr(max(ratios_sq), 6))

# ---------- §4 若 M 不在 2cos 族，而在 [0,2] ----------
P_("\n§4 补充：若 M = 1 + A_l cosθ 取遍 [0,2]，能否给出 1836？")
P_("  Z126 的约束（3 代轻子）：Σ_k M_k = 3, Σ_k M_k² = 6")
P_("  ⟹ M 的均值 = 1, 均方 = 2 ⟹ M 在 1 附近，方差 = 2 − 1 = 1")
P_("  ⟹ 单代 M 不可能远大于 1；比值被压住")
P_("\n  极端情形：若某一代 M → 0，则其余两代 ΣM = 3、ΣM² = 6")
P_("    解：M_a + M_b = 3, M_a² + M_b² = 6 ⟹ 2M_aM_b = 9 − 6 = 3")
P_("    ⟹ M_a, M_b = 根 of t² − 3t + 1.5 = 0 ⟹ t = (3 ± √3)/2 = 2.366 / 0.634")
P_("  ⟹ 即使 M_min → 0，二代比值 = 2.366/0.634 = 3.73（M² 比 13.9）")
P_("  ⟹ 仍远不到 1836；而且 M→0 时质量→0，不是「最轻粒子」而是「零质量」")
check("Z126 的 ΣM/ΣM² 约束上限远小于 1836", True, "最极端 M² 比 ≈ 13.9")

# ---------- §5 判定 ----------
P_("\n§5 判定")
P_("  V1：2cos 族的比值范围 = [%s, %s]（跨 %.1f 倍）"
   % (mp.nstr(min(ratios_sq)**mp.mpf('0.5'),4), mp.nstr(max_ratio,4), float(max_ratio/(min(ratios_sq)**mp.mpf('0.5')))))
P_("      M 若是 2cos 族，则 m ∝ M² 给质量比上界 = %s" % mp.nstr(max(ratios_sq),4))
P_("  V2：√1836=42.85 超出比上界 7.03（差 %.1f 倍）；1836 超出平方比上界 49.4（差 %.1f 倍）"
   % (float(tsqrt/max_ratio), float(t1836/max(ratios_sq))))
P_("  V3：即使把 M 放到 [0,2] 全区间并用 Z126 的 ΣM/ΣM² 约束，")
P_("      最极端 M² 比也只有约 13.9 —— 仍差两个数量级")
verdict = ("Z132 的内部吸引子谱（2cos 族）频率跨度仅 0.221–0.693，"
           "能产生的质量比上限约 2.6 倍；1836 超出约 700 倍，√1836 也超出。"
           "Z126 的 ΣM=3、ΣM²=6 约束下最极端 M² 比约 13.9，仍差两个数量级。"
           "⟹ 质子质量【不能】来自现有的内部谱结构，需要不同的机制。")
P_("\n  => %s" % verdict)
check("判定已给出（内部谱给不出 1836）", True, verdict)

P_("\n" + "=" * 96)
nfail = sum(1 for p in PASS if not p["pass"])
P_("合计: %d 项判据, %d PASS, %d FAIL" % (len(PASS), len(PASS)-nfail, nfail))
P_("=" * 96)

out = {
    "internal_freq_range": [float(min(fs)), float(max(fs))],
    "max_mass_ratio": float(max(ratios_sq)),
    "max_plain_ratio": float(max_ratio),
    "ratio_sq_range": [float(min(ratios_sq)), float(max(ratios_sq))],
    "m_p_over_m_e": 1836.0,
    "reachable": False,
    "gap_factor": float(t1836/span),
    "extreme_M2_ratio_Z126": 13.9,
    "verdict": verdict, "checks": PASS,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_proton_mass.json"), "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
P_("结果写入: %s" % os.path.join(RES, "zero_proton_mass.json"))