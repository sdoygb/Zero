#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
up_amplitude_mechanism.py —— 256/81 是机制还是巧合？
================================================================================
问题：Z126 记 A_up² = 256/81 = (4/3)^4（差 −0.006%）。这个"命中"有多特殊？

两条独立的检验：
  【一】命中的**唯一性**：用与下型同一套 528 形态去搜 A_up²，看有几个命中。
         若只有 1 个 ⟹ 它是特殊读法；若有一把 ⟹ 只是众多巧合之一。
  【二】读法的**可推广性**：$4/3$ 是亚层比梯级 r=(3/2, 2, 8/3) 的公比。
         检验三个已知振幅是否都是同一个泛函 f(r) 的值。

【检验】
  【M1】亚层比构成等比梯级，公比 4/3
  【M2】三个已知振幅在 r 梯级上的表达
  【M3】256/81 在 528 形态里的命中数与偶然期望
  【M4】"取梯级跨度"这条读法能否推广（轻子/中微子是不是跨度）
  【M5】槽数不决定振幅（轻子与中微子同为每代 2 槽，振幅却不同）

用法：/usr/local/bin/python3 up_amplitude_mechanism.py
输出：ZeroGPU/results/up_amplitude_mechanism.json
"""
from __future__ import annotations
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "up_amplitude_mechanism.json")

ok, bad = 0, []


def chk(cond, msg, detail=""):
    global ok
    if cond:
        ok += 1
        print(f"  [v] {msg}" + (f"   {detail}" if detail else ""))
    else:
        bad.append(msg)
        print(f"  [x] {msg}" + (f"   {detail}" if detail else ""))


# ── 亚层比梯级（Z129 的 L2 权重 (16,8,3,2)/29 的相邻比）──
W = (16, 8, 3, 2)
# Z129 的 r = (3/2, 2, 8/3) 是升序（三个相邻比按升序排列）
_r = sorted([W[0] / W[1], W[1] / W[2], W[2] / W[3]])
R1, R2, R3 = _r
RHO = R2 / R1

print("=" * 92)
print("【M1】亚层比梯级是等比的")
print("=" * 92)
print(f"  权重 (16,8,3,2)/29（Z129：L2 的 4 个亚层，对角锁死的稳态）")
print(f"  相邻比 r = ({R1:.6f}, {R2:.6f}, {R3:.6f}) = (3/2, 2, 8/3)")
print(f"  r2/r1 = {R2/R1:.10f}，r3/r2 = {R3/R2:.10f}")
chk(abs(R2 / R1 - 4 / 3) < 1e-12 and abs(R3 / R2 - 4 / 3) < 1e-12,
    "梯级等比，公比 = 4/3", f"{RHO:.10f}")

print()
print("=" * 92)
print("【M2】三个已知振幅在 r 梯级上的表达")
print("=" * 92)
A_NU2 = 3 / 2
A_L2 = 2.0
A_UP2 = 256 / 81
rows = [
    ("中微子", math.sqrt(A_NU2), A_NU2, "r1"),
    ("轻子", math.sqrt(A_L2), A_L2, "r2 = r1*rho"),
    ("上型", 16 / 9, A_UP2, "(r3/r1)^2 = rho^4"),
]
print(f"  {'层':>6} {'A':>10} {'A²':>12} {'读法':>18} {'对 rho 的指数':>14}")
for nm, A, A2, rd in rows:
    # A2 = 3/2 * rho^k  →  k = ln(A2/(3/2))/ln(rho)
    k = math.log(A2 / (3 / 2)) / math.log(RHO)
    print(f"  {nm:>6} {A:>10.6f} {A2:>12.8f} {rd:>18} {k:>14.6f}")
print()
print("  ★ 若三层共用同一个泛函 f(r)，指数应相同；实测：")
ks = [math.log(A2 / (3 / 2)) / math.log(RHO) for _, _, A2, _ in rows]
print(f"     k = {['%.4f' % v for v in ks]}")
chk(len(set(round(k, 3) for k in ks)) == len(ks),
    "三个指数互不相同 ⟹ 三层用了 r 梯级的**三种不同泛函**", f"{['%.4f' % v for v in ks]}")

print()
print("=" * 92)
print("【M3】256/81 在 528 形态里的命中数")
print("=" * 92)
ATOMS = {
    "1": 1.0, "2": 2.0, "3": 3.0, "4": 4.0, "6": 6.0, "8": 8.0, "16": 16.0, "29": 29.0,
    "3/2": 1.5, "4/3": 4 / 3, "8/3": 8 / 3, "16/9": 16 / 9, "5/2": 2.5,
    "sqrt2": math.sqrt(2), "sqrt3": math.sqrt(3), "sqrt{3/2}": math.sqrt(1.5),
}
forms = [(k, v) for k, v in ATOMS.items()]
for k1, v1 in ATOMS.items():
    for k2, v2 in ATOMS.items():
        forms.append((f"{k1}/{k2}", v1 / v2))
        forms.append((f"{k1}*{k2}", v1 * v2))
vals = [(n, v) for n, v in forms]

print(f"  形态数 = {len(vals)}")
print(f"  {'窗口':>8} {'命中数':>8} {'命中形态':>40} {'偶然期望':>10}")
hits_by_w = {}
for w in (0.001, 0.006, 0.01, 0.1, 1.0):
    h = [(n, v) for n, v in vals if abs(v - A_UP2) / A_UP2 < w / 100]
    hits_by_w[w] = h
    exp = len(vals) * 2 * w / 100
    show = "、".join(f"{n}={v:.6f}" for n, v in h[:4]) + ("…" if len(h) > 4 else "")
    print(f"  ±{w:>6}% {len(h):>8} {show:>40} {exp:>10.2f}")
print()
h = hits_by_w[0.01]
chk(len(h) <= 1, "±0.01% 窗口内至多 1 个形态命中（故 256/81 是唯一读法）",
    f"{len(h)} 个：{[(n, round(v,6)) for n, v in h]}")

print()
print("=" * 92)
print("【M4】「取梯级跨度」能否推广")
print("=" * 92)
SPAN = (R3 / R1) ** 2          # 上型用的泛函
print(f"  上型的泛函：A² = (r3/r1)² = {SPAN:.10f}")
print(f"  {'层':>6} {'A² 实测':>12} {'若用同一泛函':>14} {'相对差':>10}")
for nm, _, A2, _ in rows:
    print(f"  {nm:>6} {A2:>12.8f} {SPAN:>14.8f} {100*(SPAN-A2)/A2:>+9.2f}%")
others = [(nm, A2) for nm, _, A2, _ in rows if nm != "上型"]
chk(all(abs(SPAN - A2) / A2 > 0.3 for _, A2 in others),
    "另两层若用同一泛函会差 30% 以上 ⟹ 该泛函不通用",
    "最小差 " + "%.1f%%" % min(100 * abs(SPAN - A2) / A2 for _, A2 in others))

print()
print("=" * 92)
print("【M5】槽数不决定振幅")
print("=" * 92)
# 每代槽数（Z136 槽表）：轻子 2、中微子 2、上型 6、下型 6
SLOTS = {"轻子": 2, "中微子": 2, "上型": 6, "下型": 6}
AMP2 = {"轻子": A_L2, "中微子": A_NU2, "上型": A_UP2, "下型": 7.197684}
print(f"  {'层':>6} {'每代槽数':>9} {'A²':>12}")
for nm in SLOTS:
    print(f"  {nm:>6} {SLOTS[nm]:>9} {AMP2[nm]:>12.8f}")
print()
chk(AMP2["轻子"] != AMP2["中微子"] and SLOTS["轻子"] == SLOTS["中微子"],
    "轻子与中微子同为每代 2 槽，A² 却不同 ⟹ 振幅不是槽数的函数",
    f"A²={AMP2['轻子']:.4f} vs {AMP2['中微子']:.4f}")
chk(SLOTS["上型"] == SLOTS["下型"] and AMP2["上型"] != AMP2["下型"],
    "上型与下型同为每代 6 槽，A² 也不同", f"{AMP2['上型']:.4f} vs {AMP2['下型']:.4f}")

out = {
    "r_ladder": [R1, R2, R3], "rho": RHO,
    "A2": {"neutrino": A_NU2, "lepton": A_L2, "up": A_UP2, "down": AMP2["下型"]},
    "rho_exponents": ks,
    "hits_by_window": {str(k): [{"form": n, "value": v} for n, v in v] for k, v in hits_by_w.items()},
    "form_count": len(vals),
    "assertions": {"passed": ok, "failed": len(bad), "bad": bad},
}
os.makedirs(RES, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print()
print(f"断言：{ok} 条通过 / {len(bad)} 条不符")
for x in bad:
    print("  ✗", x)
print(f"已存 {OUT}")
raise SystemExit(1 if bad else 0)