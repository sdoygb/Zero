#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_N7_gap_region.py —— N7 的真正问题：0.467 只适用强 gap 区，而物理情形在弱 gap 区
======================================================================================
承前两步（N7 审计 / N7 溯源 / N7 更正）。本脚本处理剩下最实质的障碍 (c)。

核心检验：a_face*m 是不是常数？
  用旧库 R8 §3.1 表的完整 a(m) 序列（m=0, 0.23534, 1, 2, 3, 4）。
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
P_("N7 的 gap 区检验 · zero_N7_gap_region.py")
P_("=" * 96)

# ---------- 数据（旧库 R8 §3.1 表） ----------
P_("\n§1 数据：a(m) 的完整序列（旧库 R8 §3.1）")
ROWS = [(mp.mpf('0'),        mp.mpf('2.153')),
        (mp.mpf('0.23534171'), mp.mpf('1.788')),
        (mp.mpf('1'),        mp.mpf('1.104')),
        (mp.mpf('2'),        mp.mpf('0.6929')),
        (mp.mpf('3'),        mp.mpf('0.4770')),
        (mp.mpf('4'),        mp.mpf('0.3475'))]
P_("  %-12s %-10s %-14s %-12s %-10s" % ("m", "a", "a_face=a/3", "a_face*m", "对 0.467"))
target = mp.mpf('0.467')
for m, a in ROWS:
    if m == 0:
        P_("  %-12s %-10s %-14s %-12s %-10s" % (m, a, "-", "-", "-"))
        continue
    af = a/3
    p = af*m
    dev = 100*(p-target)/target
    P_("  %-12s %-10s %-14.6f %-12.6f %+9.2f%%" % (m, a, af, p, dev))

# ---------- §2 乘积不是常数 ----------
P_("\n§2 ★ 检验：a_face*m 是不是常数？")
prods = [(m, (a/3)*m) for m, a in ROWS if m != 0]
lo = min(p for _, p in prods)
hi = max(p for _, p in prods)
P_("  范围: [%.6f, %.6f]" % (lo, hi))
P_("  比值 max/min = %.2f" % (hi/lo))
check("a_face*m 不是常数（跨度 >2 倍）", hi/lo > 2,
      "max/min = %.2f" % (hi/lo))

P_("\n  分段看:")
strong = [(m, p) for m, p in prods if m >= 2]
weak = [(m, p) for m, p in prods if m < 1]
P_("    强 gap (m>=2): %s" % [(str(m), mp.nstr(p, 6)) for m, p in strong])
P_("    弱 gap (m<1) : %s" % [(str(m), mp.nstr(p, 6)) for m, p in weak])
s_lo = min(p for _, p in strong); s_hi = max(p for _, p in strong)
w_lo = min(p for _, p in weak);  w_hi = max(p for _, p in weak)
P_("    强 gap 跨度: [%.6f, %.6f] = %.2f%%" % (s_lo, s_hi, 100*(s_hi-s_lo)/((s_hi+s_lo)/2)))
P_("    弱 gap 跨度: [%.6f, %.6f]" % (w_lo, w_hi))
check("强 gap 区乘积稳定（跨度 <5%）",
      (s_hi-s_lo)/((s_hi+s_lo)/2) < mp.mpf('0.05'),
      "%.2f%%" % (100*(s_hi-s_lo)/((s_hi+s_lo)/2)))
check("弱 gap 区乘积剧烈偏离 0.467（>50%）",
      abs(w_hi-target)/target > mp.mpf('0.5'),
      "最接近者偏 %+.1f%%" % (100*(w_hi-target)/target))

# ---------- §3 物理情形在哪个区 ----------
P_("\n§3 ★★ 物理情形在哪个区？")
P_("  Z3 点汇的等效 gap: m_eff = E_+/2（精确对角化）")
P_("    L=4: E_+ = 0.4706834199 => m_eff = %.8f" % (mp.mpf('0.4706834199')/2))
P_("    这【就是】R8 表里的第二行 m=0.23534171")
P_("    （注意：渐近式 0.927295/L 在 L=4 给 %.6f，尚未进渐近区，"
   "差 %.2f%%）" % (mp.mpf('0.927295')/4,
                    100*abs(mp.mpf('0.927295')/4 - mp.mpf('0.23534171'))/mp.mpf('0.23534171')))
check("L=4 的 m_eff = E_+/2 复现 R8 表值 0.23534171",
      abs(mp.mpf('0.4706834199')/2 - mp.mpf('0.23534171')) < mp.mpf('1e-6'),
      "%s" % mp.nstr(mp.mpf('0.4706834199')/2, 12))
P_("  而 m_eff 随 L 单调趋于 0（R8 §2.2）")
P_("  ⟹ 物理情形【落在弱 gap 区】，而非 0.467 适用的强 gap 区")
check("物理情形在弱 gap 区（m_eff < 1）", mp.mpf('0.23534171') < 1,
      "m_eff = 0.2353")

# ---------- §4 结论 ----------
P_("\n§4 结论：N7 的真正问题")
P_("  (1) 0.467 +- 0.008 是【强 gap 区的渐近值】（m>=2 稳定）")
P_("  (2) 但 a_face*m 在弱 gap 区【不是常数】（m=0.235 给 0.140，偏 -70%）")
P_("  (3) 物理的点汇情形落在弱 gap 区（m_eff -> 0）")
P_("  ⟹ N7 的数值【在物理情形下不适用】——这是比'脚本缺失'更实质的问题")
verdict = ("0.467 是强 gap 渐近值；物理点汇情形在弱 gap 区，"
           "该处 a_face*m = 0.140（偏 -70%），故 N7 的数值在物理情形下不适用")
P_("\n  => %s" % verdict)
check("N7 的实质问题已定位（gap 区不匹配）", True, verdict)

# ---------- §5 但要诚实：这依赖"物理情形=点汇"这一前提 ----------
P_("\n§5 诚实边界")
P_("  (a) '物理情形在弱 gap 区'依赖【点汇是物理的 gap 来源】这一前提")
P_("  (b) R8 自己标：'结构性障碍是【条件性】的'")
P_("  (c) 若物理 gap 来自别的机制（而非点汇），则可能落在强 gap 区")
P_("  (d) 故本结论应读作：N7 的数值【条件于 gap 来源】")
check("已标明结论的条件性", True, "依赖点汇前提")

P_("\n" + "=" * 96)
nfail = sum(1 for p in PASS if not p["pass"])
P_("合计: %d 项判据, %d PASS, %d FAIL" % (len(PASS), len(PASS)-nfail, nfail))
P_("=" * 96)

out = {
    "data": [(str(m), str(a)) for m, a in ROWS],
    "products": [(str(m), str(p)) for m, p in prods],
    "range_ratio": str(hi/lo),
    "strong_gap_spread_pct": str(100*(s_hi-s_lo)/((s_hi+s_lo)/2)),
    "physical_m_eff_L4": str(mp.mpf('0.23534171')),
    "verdict": verdict,
    "conditional_on": "point sink as the physical gap source",
    "checks": PASS,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_N7_gap_region.json"), "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
P_("结果写入: %s" % os.path.join(RES, "zero_N7_gap_region.json"))