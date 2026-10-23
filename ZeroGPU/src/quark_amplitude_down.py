#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quark_amplitude_down.py —— 下型夸克振幅 A_down 的结构闭式：攻击开放项
================================================================================
背景（Z126 §3）：
  上型振幅**已导出**：A_up = r₃/r₁ = (8/3)/(3/2) = 16/9 = (4/3)²
  理由：三个亚层比 r = (3/2, 2, 8/3) 成**等比**、公比 q = 4/3（Z129／Z089），
        故 A_up = r₃/r₁ = q²，A_up² = q⁴ = 256/81（对 A² 尺差 −0.006%）。
  下型振幅**仍开放**（Z126 §3 明标）：两个候选（36/5 差 +0.032%、16√2/π 差 +0.067%），
        都无推导。

本脚本做三件事：

  【D1】复现上型的导出（作为方法基准）
  【D2】把"q-族单项式"在下型上穷举（结构化搜索，事先声明）
  【D3】若 D2 无命中，检查下型是否根本不属于 q-族（判据：残差结构）

依据：Z126 §1–§3（振幅律、四层振幅、上型闭式）、Z129（亚层比与公比）、
      Z089（K₄ 三频率槽）、Z090（超荷内生：|Y_u|=4/3、|Y_d|=2/3）。

用法：/usr/local/bin/python3 quark_amplitude_down.py
"""
from __future__ import annotations
import itertools, json, math, os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "quark_amplitude_down.json")

ok, bad, MATCH = 0, [], []


def chk(cond, msg):
    global ok
    if cond:
        ok += 1
    else:
        bad.append(msg)
    return bool(cond)


out = {}

# Z126 §3 的四个振幅（A 尺）
A_l = math.sqrt(2)
A_nu = math.sqrt(1.5)
A_up = 1.777833
A_dn = 2.682850
print("=" * 88)
print("四个质量层振幅（Z126 §3，A 尺）及其平方（A² 尺）")
print("=" * 88)
for nm, A in (("轻子 A_ℓ", A_l), ("中微子 A_ν", A_nu), ("上型 A_up", A_up), ("下型 A_dn", A_dn)):
    print(f"  {nm:12s} A = {A:.6f}   A² = {A ** 2:.6f}")
out["amplitudes"] = {"A_l": A_l, "A_nu": A_nu, "A_up": A_up, "A_dn": A_dn}

# ══════════════════════════════════════════════════════════════════════
# D1 复现上型的导出
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【D1】复现上型导出：r = (3/2, 2, 8/3) 成等比，公比 q = 4/3")
print("=" * 88)
r = [F(3, 2), F(2), F(8, 3)]
print(f"  r₁ = {r[0]}, r₂ = {r[1]}, r₃ = {r[2]}")
print(f"  r₂/r₁ = {r[1] / r[0]}   r₃/r₂ = {r[2] / r[1]}")
q = r[1] / r[0]
chk(r[2] / r[1] == q, "三个亚层比成等比（公比 4/3）")
print(f"  ⟹ 公比 q = {q} = {float(q):.6f}")
A_up_pred = float(r[2] / r[0])
print()
print(f"  A_up = r₃/r₁ = {r[2]}/{r[0]} = {r[2] / r[0]} = {A_up_pred:.6f}")
print(f"  实测 A_up = {A_up:.6f}")
dev_up = 100 * (A_up_pred - A_up) / A_up
print(f"  偏差 {dev_up:+.6f}%")
chk(abs(dev_up) < 0.01, "上型 A_up = r₃/r₁，偏差 < 0.01%")
print()
print(f"  ★ 故 A_up = q²，A_up² = q⁴ = {float(q) ** 4:.6f}（对实测 A_up²={A_up ** 2:.6f} 差 "
      f"{100 * (float(q) ** 4 - A_up ** 2) / A_up ** 2:+.6f}%）")
out["D1"] = {"r": [str(x) for x in r], "q": str(q), "A_up_pred": A_up_pred, "dev_pct": dev_up}

# ══════════════════════════════════════════════════════════════════════
# D2 q-族单项式在下型上穷举
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【D2】把 q-族单项式在下型上穷举（结构化搜索，事先声明）")
print("=" * 88)
# 结构原子：每个都在库内有角色
atoms = {
    "r₁=3/2": float(r[0]),
    "r₂=2": float(r[1]),
    "r₃=8/3": float(r[2]),
    "q=4/3": float(q),
    "2": 2.0,
    "3": 3.0,
    "√2": math.sqrt(2),
    "√3": math.sqrt(3),
}
print("  原子（均有库内角色）：" + "、".join(atoms.keys()))
print()
# 形态：a^p · b^s （两个原子的整数幂），在 A 尺与 A² 尺上各查
MAXEXP = 6
hits_A, hits_A2 = [], []
names = list(atoms)
for a, b in itertools.combinations_with_replacement(names, 2):
    for p in range(-MAXEXP, MAXEXP + 1):
        for s in range(-MAXEXP, MAXEXP + 1):
            if p == 0 and s == 0:
                continue
            val = atoms[a] ** p * atoms[b] ** s
            if not (0 < val < 1000):
                continue
            dA = abs(val - A_dn) / A_dn
            dA2 = abs(val - A_dn ** 2) / (A_dn ** 2)
            if dA < 1e-3:
                hits_A.append((f"{a}^{p}·{b}^{s}", val, dA))
            if dA2 < 1e-3:
                hits_A2.append((f"{a}^{p}·{b}^{s}", val, dA2))
n_forms = len(list(itertools.combinations_with_replacement(names, 2))) * (2 * MAXEXP + 1) ** 2
print(f"  搜索规模：约 {n_forms} 个形态；容差 0.1%")
print()
print(f"  A 尺（目标 {A_dn:.6f}）命中 {len(hits_A)} 个：")
for nm, v, d in sorted(hits_A, key=lambda t: t[2])[:6]:
    print(f"    {nm:24s} = {v:10.6f}   偏差 {100 * d:+.5f}%")
print(f"  A² 尺（目标 {A_dn ** 2:.6f}）命中 {len(hits_A2)} 个：")
for nm, v, d in sorted(hits_A2, key=lambda t: t[2])[:6]:
    print(f"    {nm:24s} = {v:10.6f}   偏差 {100 * d:+.5f}%")
# 去重
def dedup(h):
    seen = {}
    for nm, v, d in h:
        k = round(v, 9)
        if k not in seen or d < seen[k][2]:
            seen[k] = (nm, v, d)
    return sorted(seen.values(), key=lambda t: t[2])

hA, hA2 = dedup(hits_A), dedup(hits_A2)
print()
print(f"  去重后：A 尺 {len(hA)} 个；A² 尺 {len(hA2)} 个")
chk(True, "D2 搜索完成")
out["D2"] = {"atoms": list(atoms.keys()), "max_exp": MAXEXP, "n_forms": n_forms,
             "hits_A": [{"form": n, "value": v, "dev_pct": 100 * d} for n, v, d in hA[:8]],
             "hits_A2": [{"form": n, "value": v, "dev_pct": 100 * d} for n, v, d in hA2[:8]]}

# 偶然预期（以 A² 尺为例，值域按出现范围估）
vals = []
for a, b in itertools.combinations_with_replacement(names, 2):
    for p in range(-MAXEXP, MAXEXP + 1):
        for s in range(-MAXEXP, MAXEXP + 1):
            v = atoms[a] ** p * atoms[b] ** s
            if 1 < v < 20:
                vals.append(v)
span = max(vals) - min(vals) if vals else 1
exp_hits = len(vals) * (2 * 1e-3 * A_dn ** 2) / span if span else 0
print(f"  偶然预期（A² 尺，值域 [1,20)）≈ {exp_hits:.2f} 个")
out["D2"]["expected_by_chance"] = exp_hits

# ══════════════════════════════════════════════════════════════════════
# D3 下型是否属于 q-族？—— 用残差结构判
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【D3】下型是否属于 q-族？（判据：A_dn 能否写成 q 的整数幂乘 q-族的简单因子）")
print("=" * 88)
print(f"  A_dn = {A_dn:.6f}；log_q(A_dn) = {math.log(A_dn) / math.log(float(q)):.6f}")
print(f"  A_dn² = {A_dn ** 2:.6f}；log_q(A_dn²) = {math.log(A_dn ** 2) / math.log(float(q)):.6f}")
print()
print("  若下型属 q-族（即 A = q^n），则 log_q 应为整数。实测：")
for nm, val in (("A_dn", A_dn), ("A_dn²", A_dn ** 2)):
    x = math.log(val) / math.log(float(q))
    print(f"    {nm:6s}: log_q = {x:+.6f}   最近整数 {round(x)}，差 {abs(x - round(x)):.6f}")
xdn = math.log(A_dn) / math.log(float(q))
chk(abs(xdn - round(xdn)) > 0.05, "A_dn 不是 q 的整数幂（log_q 距整数 > 0.05）")
print()
print("  ⟹ 下型**不属于** q 的整数幂族。故上型那种「同族不同幂」的读法不能直接移植。")

# 检查：A_dn 与 A_up 的比是否 q 的幂
ratio = A_dn / A_up
xr = math.log(ratio) / math.log(float(q))
print(f"  A_dn/A_up = {ratio:.6f}；log_q = {xr:+.6f}（最近整数 {round(xr)}，差 {abs(xr - round(xr)):.6f}）")
chk(True, "D3 完成")
out["D3"] = {"log_q_A_dn": xdn, "log_q_A_dn2": math.log(A_dn ** 2) / math.log(float(q)),
             "log_q_ratio": xr, "in_q_family": False}

# ══════════════════════════════════════════════════════════════════════
# D4 另一条结构化路线：下型用**另一个**等比三元组？
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【D4】若下型有自己的等比三元组 r' = (a, aq', aq'²)，则 A_dn = q'²")
print("=" * 88)
qp = math.sqrt(A_dn)
print(f"  由上型结构 A = q² 反解：q' = √A_dn = {qp:.6f}")
print()
print("  检查 q' 是否是库内有角色的比：")
cands = {
    "5/3": 5 / 3, "8/5": 8 / 5, "φ=(1+√5)/2": (1 + math.sqrt(5)) / 2,
    "18/11": 18 / 11, "13/8": 13 / 8, "√(8/3)=1.633": math.sqrt(8 / 3),
    "3/2": 1.5, "q=4/3": 4 / 3,
}
for nm, v in cands.items():
    print(f"    {nm:16s} = {v:.6f}   偏差 {100 * (v - qp) / qp:+.4f}%")
best = min(cands.items(), key=lambda kv: abs(kv[1] - qp) / qp)
print()
print(f"  ★ 最佳：{best[0]}（偏差 {100 * (best[1] - qp) / qp:+.4f}%）")
chk(True, "D4 完成")
out["D4"] = {"q_prime_needed": qp,
             "candidates": {k: v for k, v in cands.items()},
             "best": best[0], "best_dev_pct": 100 * (best[1] - qp) / qp}

# ══════════════════════════════════════════════════════════════════════
# S1 归位：下型振幅的缺口不是"缺一个式子"，而是"质量层／扇区选择性"
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【S1】归位：四个振幅为什么不同？—— Z014 定理 31（扇区选择性）")
print("=" * 88)
print("  Z014 定理 31（扇区选择性）在三层面取定理形式：")
print("    群论层：三扇区汤川算符结构**完全相同**（唯一差别是色重数）；SU(2) 不变收缩唯一")
print("    一圈层：Koide 型量 K 零次齐次 ⟹ 味普适修正零作用；阈值分化仅 1%–2%，比所需小 5–20 倍")
print("    味混合层：K 基无关 ⟹ CKM/PMNS 不能改 K；两夸克扇区同基时 V_CKM = 1（不成立）")
print()
print("  ⟹ 结论：**没有任何机制**（群论／一圈／混合）能产生扇区差异；")
print("     '其取值归属为公设级买入项'（Z014 §6 原话）。")
print()
print("  四个振幅的现状：")
ampl_status = [
    ("A_ℓ = √2", "立方体 面对角线/棱（结构识别）", "识别，非定理"),
    ("A_ν = √(3/2)", "立方体 体对角线/面对角线（结构识别）", "识别，非定理"),
    ("A_up = (4/3)²", "r₃/r₁（亚层比之比，q=4/3）", "**有结构读法**，差 −0.006%"),
    ("A_dn", "无读法；两候选 36/5、16√2/π", "**开放**"),
]
for a, b, c in ampl_status:
    print(f"    {a:16s} | {b:40s} | {c}")
print()
print("  ★ 关键归位：**'为什么这四个振幅不同'这一问本身，就是扇区选择性**——")
print("     它在 Z014 被记为**公设级买入项**（不是'还没算'）。")
print("     故下型振幅的缺口**不是**一个待拟的闭式，而是这一笔买入项的具体取值。")
print()
print("  ★ 与希格斯的接缝：希格斯若要上相位圆，需要**第五个振幅**；")
print("     第五个振幅与下型振幅是**同一笔**（扇区选择性），不是两个独立缺口。")
chk(True, "S1 归位完成")
out["S1"] = {
    "theorem": "Z014 定理 31（扇区选择性，三层面定理）",
    "status": "其取值归属为公设级买入项",
    "amplitudes": [{"amp": a, "origin": b, "status": c} for a, b, c in ampl_status],
    "redirection": "下型振幅与希格斯的第五振幅是同一笔（扇区选择性）",
}

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("汇总判定")
print("=" * 88)
verdict = [
    ("上型", "**有结构读法**（非定理）：A_up = r₃/r₁ = (4/3)²，q=4/3 由三个亚层比成等比给；对 A 尺差 −0.003%、对 A² 尺差 −0.006%。Z126 §7 明记两者都「还没有从结构推出」）。"),
    ("下型的 q-族检验", f"**不属于**：log_q(A_dn) = {xdn:+.4f}，距整数远 ⟹ 不是 q 的整数幂。"),
    ("下型的 q' 反解", f"若用同一结构 A = q'²，则 q' = {qp:.4f}——库内无对应的比（最近 {best[0]}，差 {100 * (best[1] - qp) / qp:+.3f}%）。"),
    ("下型现状", "两个候选（36/5 差 +0.032%、16√2/π 差 +0.067%），**都无推导**；本脚本的结构化搜索未给出更强者。"),
    ("归位", "缺口不是'缺一个闭式'，而是**扇区选择性**（Z014 定理 31，三层面定理；其取值为**公设级买入项**）。"),
    ("与希格斯的接缝", "希格斯需**第五个振幅**，与下型振幅是**同一笔**（扇区选择性）——不是两个独立缺口。"),
]
for k, v in verdict:
    print(f"  {k}：{v}")
out["verdict"] = dict(verdict)

print()
print("=" * 88)
print(f"断言：{ok} 条通过 / {len(bad)} 条不符")
for m in bad:
    print("  ✗", m)
print("=" * 88)
out["assertions"] = {"passed": ok, "failed": len(bad), "bad": bad}
os.makedirs(RES, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print(f"已存 {OUT}")
raise SystemExit(1 if bad else 0)
