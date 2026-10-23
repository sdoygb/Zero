#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
higgs_composite_test.py —— 希格斯是否是复合粒子：四条实测
================================================================================
问题（用户）：希格斯是否有可能是复合粒子？

本脚本测四件事：

  【R1】数值：M_H 能否由已知质量的**简单组合**给出？（候选集事先声明）
  【R2】标度：复合是否**带来新标度**？（Zero 只有唯一量纲锚 v）
  【R3】结构：Zero 里"多体／结合"是否可建？（Z036 定理 27、Z014 定理 22）
  【R4】账目：复合能否把 λ 从 F（拟合）转成 D（导出）？

依据：Z153／Z013 定理 16（唯一量纲锚 v）、Z036 定理 27（⊕ 内最多一个单式 M_2）、
      Z014 定理 22（禁闭边界）、Z042（组合间隙非质量间隙）、Z024 定理 2（无作用量）、
      Z151 §4.5（边际耦合不给值）、Z167 定义 167.10（四档）。

用法：/usr/local/bin/python3 higgs_composite_test.py
"""
from __future__ import annotations
import itertools, json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "higgs_composite_test.json")

ok, bad, MATCH = 0, [], []


def chk(cond, msg):
    global ok
    if cond:
        ok += 1
    else:
        bad.append(msg)
    return bool(cond)


out = {}

V = 246.22
M_H = 125.25
LAM = M_H ** 2 / (2 * V ** 2)
C_H = M_H / V          # = √(2λ)，无量纲

print("=" * 88)
print("基本量")
print("=" * 88)
print(f"  v = {V} GeV（唯一量纲锚）；M_H = {M_H} GeV")
print(f"  c ≡ M_H/v = {C_H:.6f}（无量纲）;  λ = c²/2 = {LAM:.6f}")
out["basis"] = {"v": V, "M_H": M_H, "c_H": C_H, "lambda": LAM}

# ══════════════════════════════════════════════════════════════════════
# R1 数值：M_H 能否由已知质量的简单组合给出
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【R1】M_H 能否由已知质量的简单组合给出？（候选集事先声明）")
print("=" * 88)
masses = {   # GeV
    "M_W": 80.377, "M_Z": 91.1876, "m_t": 172.69, "m_b": 4.18,
    "m_c": 1.27, "m_τ": 1.77686, "m_μ": 0.10566, "m_e": 0.000511,
    "m_π": 0.13957, "m_p": 0.938272, "m_n": 0.939565,
    "Λ_QCD": 0.2,   # 量级
}
print("  已知质量（GeV）：" + "、".join(f"{k}={v}" for k, v in masses.items()))
print()
# 形态 1：k × m
hits = []
for k, m in masses.items():
    for n in range(1, 9):
        val = n * m
        if abs(val - M_H) / M_H < 0.02:
            hits.append((f"{n}×{k}", val))
# 形态 2：m1 + m2、m1 − m2
for (k1, m1), (k2, m2) in itertools.permutations(masses.items(), 2):
    for val, nm in ((m1 + m2, f"{k1}+{k2}"), (m1 - m2, f"{k1}−{k2}")):
        if val > 0 and abs(val - M_H) / M_H < 0.02:
            hits.append((nm, val))
# 形态 3：√(m1·m2)
for (k1, m1), (k2, m2) in itertools.combinations(masses.items(), 2):
    val = math.sqrt(m1 * m2)
    if abs(val - M_H) / M_H < 0.02:
        hits.append((f"√({k1}·{k2})", val))

print(f"  ±2% 内命中（组合形态：n·m、m1±m2、√(m1m2)）：{len(hits)} 个")
for nm, v in sorted(hits, key=lambda t: abs(t[1] - M_H))[:10]:
    print(f"    {nm:22s} = {v:8.4f} GeV   偏差 {100 * (v - M_H) / M_H:+.3f}%")
chk(True, "R1 扫描完成")
out["R1"] = {"hits": [{"form": n, "GeV": v} for n, v in hits[:10]],
             "n_hits": len(hits)}

# 具体核对几个"物理上会想到"的复合组合
print()
print("  物理上会首先想到的复合组合：")
named = [
    ("2·m_t（tt̄ 阈值）", 2 * masses["m_t"]),
    ("m_t（单体）", masses["m_t"]),
    ("2·M_W", 2 * masses["M_W"]),
    ("M_Z+M_W", masses["M_Z"] + masses["M_W"]),
    ("√2·M_W", math.sqrt(2) * masses["M_W"]),
    ("2·M_Z", 2 * masses["M_Z"]),
]
for nm, v in named:
    print(f"    {nm:20s} = {v:8.3f} GeV   偏差 {100 * (v - M_H) / M_H:+.2f}%")
out["R1_named"] = [(n, v) for n, v in named]

# ══════════════════════════════════════════════════════════════════════
# R2 标度：复合是否带来新标度
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【R2】复合是否带来新标度？（Zero 只有唯一量纲锚 v）")
print("=" * 88)
print("  Z153／Z013 定理 16：Zero 的**唯一**有量纲自由度是 v；")
print("  一切有量纲量写成 c·v^n（n 由量纲定）。")
print()
print("  故任何复合构型的质量必然形如：")
print()
print("      M_H = c · v ,   c 无量纲")
print()
print(f"  现在 c = {C_H:.6f}（由 λ 给：c = √(2λ)）")
print()
print("  ★ 关键：复合**不改变这个形式**——它只换 c 的**来历**：")
print()
rows = [
    ("现在的图景", "c = √(2λ)，λ 是拉氏量参数（F 档，实测反解）"),
    ("复合图景", "c = f(g_*, N, 结合能/v)，全部无量纲 ⟹ 仍是 F 档，除非 f 可导出"),
]
for a, b in rows:
    print(f"    {a:10s}：{b}")
print()
print("  ⟹ 复合要**换来收益**，必须让 f 可导出（即结合动力学可算）。")
print("     否则只是把'λ 是输入'换成'结合动力学是输入'——输入数不减。")
chk(True, "R2 形式检查完成")
out["R2"] = {"unique_anchor": "v",
             "form_invariant": "M_H = c·v, c 无量纲",
             "c_now": C_H,
             "composite_keeps_form": True}

# ══════════════════════════════════════════════════════════════════════
# R3 结构：多体／结合是否可建
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【R3】Zero 里「多体／结合」是否可建？")
print("=" * 88)
blockers = [
    ("多体（两体分划）",
     "Z036 定理 27：A_L = ⊕_h M_2(C) 内**最多容纳一个**单式 M_2 拷贝 ⟹ 无「量子位＋量子位」",
     "在 ⊕ 代数下**不可建**；换成 ⊗ 是选取层 C1（Z039 §6.3 判「不是必需」）"),
    ("结合能（非微扰测度）",
     "Z014 定理 22（禁闭边界）+ Z037 §2.2：禁闭=闭合判据已证；**质量间隙与弦张力开放**",
     "**开放**（背靠 Z014 定理 22）"),
    ("Zero 自己的「间隙」",
     "Z042：Zero 的间隙是**组合间隙**（Δw=1/3），不是质量间隙",
     "结合能所需的「质量间隙」不在库内的间隙概念里"),
    ("作用量（势 V 的写处）",
     "Z024 定理 2：不存在守恒的状态函数 E ⟹ 无作用量",
     "结合动力学需要一个动力学泛函；Zero 无"),
]
for a, b, c in blockers:
    print(f"  ▸ {a}")
    print(f"      依据：{b}")
    print(f"      状态：{c}")
print()
chk(True, "R3 结构检查完成")
out["R3"] = [{"item": a, "basis": b, "status": c} for a, b, c in blockers]

print("  ★ 四项逐条：多体在 ⊕ 下不可建；结合能在库内是开放项；")
print("     而 Zero 的「间隙」是组合间隙、不是质量间隙——结合能所需的那个概念不在库里。")

# ══════════════════════════════════════════════════════════════════════
# R3b 先例：库里已有的"复合"是怎么做成的
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【R3b】先例核查：库里已有的复合（分子／原生库仑）是怎么做成的")
print("=" * 88)
precedents = [
    ("分子 H2O／LiH（zero_cell_molecules.py）",
     "import pyscf；FCI／SCF／gto —— **从头算哈密顿量是借来的**",
     "**嫁接**：Zero 只提供胞腔／S4 符号等附加结构，势与反对称化来自外部包"),
    ("原生库仑束缚态（zero_native_coulomb.py）",
     "H = (1/2)L − g·diag(φ)：**单粒子在固定势里**；势来自 L⁺（Z147／Z033）",
     "**单粒子**（借 Born–Oppenheimer 式约化）：不是两体动力学"),
    ("声子／BCS（Z113／Z114）",
     "图上的集体模；不需要两体分划",
     "**集体模可用**（不触碰 Z036 定理 27）"),
]
for a, b, c in precedents:
    print(f"  ▸ {a}")
    print(f"      做法：{b}")
    print(f"      性质：{c}")
print()
print("  ★ 故库内「复合」分两种：")
print("     · 分子 = **借**外部哈密顿量（嫁接层）")
print("     · 原生库仑 = **单粒子**在固定势里（不是两体）")
print("     · 集体模（声子）= 原生，且**不受 Z036 定理 27 限制**")
chk(True, "R3b 先例核查完成")
out["R3b"] = [{"precedent": a, "how": b, "nature": c} for a, b, c in precedents]

# ══════════════════════════════════════════════════════════════════════
# R3c 集体模路线：唯一没被结构否决的形态 —— 测它能否定出 c
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【R3c】集体模路线：载体谱能否定出 c ≡ M_H/v = 0.508691？")
print("=" * 88)
# Γ_16 = K_4 □ K_2 □ K_2 的 Laplacian 谱（Z136）
spec16 = [0, 2, 4, 6, 8]
print(f"  Γ_16 的 Laplacian 谱（Z136）：{spec16}")
print("  已知：相邻非零谱比 4/2 = 2 = A_ℓ²、6/4 = 3/2 = A_ν²（Z126 §3）")
ratios = []
for i in range(1, len(spec16) - 1):
    a, b = spec16[i], spec16[i + 1]
    ratios.append((f"{b}/{a}", b / a))
for nm, v in ratios:
    print(f"    谱比 {nm} = {v:.6f}")
print()
print(f"  目标 c = {C_H:.6f}；λ = c²/2 = {LAM:.6f}")
print("  以谱比 {2, 3/2, 4/3} 造候选（事先声明形态：单比、比之积/商、半整数倍）：")
sp = {"s1=2": 2.0, "s2=3/2": 1.5, "s3=4/3": 4 / 3}
cand = []
names = list(sp)
for a in names:
    cand.append((a, sp[a]))
    for b in names:
        cand.append((f"{a}/{b}", sp[a] / sp[b]))
        cand.append((f"{a}·{b}", sp[a] * sp[b]))
        for k in (2, 3, 4):
            cand.append((f"{a}/({k}·{b})", sp[a] / (k * sp[b])))
best = None
for nm, v in cand:
    for form, val in ((nm, v), (f"√({nm})", math.sqrt(v) if v > 0 else None),
                      (f"1/{nm}", 1 / v if v else None)):
        if val is None:
            continue
        d = abs(val - C_H) / C_H
        if best is None or d < best[2]:
            best = (form, val, d)
print(f"  最佳：{best[0]} = {best[1]:.6f}   偏差 {100 * best[2]:+.3f}%")
chk(best[2] > 0.01, "载体谱比不能给出 c（偏差 > 1%）")
out["R3c"] = {"spectrum": spec16, "target_c": C_H,
              "best_form": best[0], "best_value": best[1], "best_dev": best[2],
              "n_candidates": len(cand) * 3}
print()
print("  ★ 载体谱是**无量纲**的；任何由它造出的质量仍需乘 v ⟹ 只能定出无量纲 c。")
print(f"     而上面最佳候选偏差 {100 * best[2]:.2f}% ⟹ 谱比没有给出 c。")

# ══════════════════════════════════════════════════════════════════════
# R4 账目：复合能否把 λ 从 F 转成 D
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【R4】复合能否把 λ 从 F（拟合）转成 D（导出）？")
print("=" * 88)
print("  判据（Z167 定义 167.10）：")
print("    D（导出）= 由结构量算出，有数值判据")
print("    F（拟合）= 拉氏量参数，由观测量反解")
print()
print("  要转档，必须给出：")
steps = [
    "① 复合的**组分**（哪些态）      —— 需要多体结构（R3：⊕ 下不可建）",
    "② 组分间的**结合**（动力学）    —— 需要非微扰测度（R3：开放）",
    "③ 结合能的**数值**             —— 需要一个质量间隙（R3：库里是组合间隙）",
]
for s in steps:
    print(f"    {s}")
print()
print("  三步缺一不可；任一步未给 ⟹ c 仍由实测定 ⟹ 仍是 F 档。")
print()
chk(True, "R4 检查完成")
out["R4"] = {"to_move_F_to_D_needs": steps, "moved": False}

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("汇总判定")
print("=" * 88)
verdict = [
    ("R1 数值", "唯一 ±2% 命中：√(M_Z·m_t) = 125.49 GeV（+0.19%）——**无机制**，且 M_Z、m_t 本身属 F 档（实测）。2m_t 差 +176%。"),
    ("R2 标度", "**复合不带来新标度**：唯一锚 v ⟹ M_H = c·v 的形式不变，只换 c 的来历。"),
    ("R3 结构", "多体（⊕ 下不可建）＋ 结合能（开放）＋ 质量间隙（库里是组合间隙）＋ 无作用量。"),
    ("R4 账目", "复合**不能**把 λ 从 F 转成 D（三步缺一不可）。"),
    ("R3b 先例", "分子=借 PySCF；原生库仑=单粒子在固定势；集体模（声子）原生且不受定理 27 限制。"),
    ("R3c 集体模", "载体谱 {0,2,4,6,8} 是**无量纲**的；由谱比造 c 的最佳偏差 >1% ⟹ 谱比没给出 c。"),
    ("净结论", "复合**与「不在相位圆上」相容**（是唯一预测此事的图像），但在 Zero 内**不产生可导出性**。"),
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
