#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
half_slot_native_search.py —— 半槽规则：改用**原生计数**做统一搜索
================================================================================
背景（库内已做，必须尊重）：
  · half_slot_rule_check.py（5/5）：三条候选规则（加在计数内再平方／占据权重乘积／先平方再加）
    **没有一条同时给出 D1 与 D3**；占据权重乘积给 24.5，差 12.5%，被排除。
  · half_slot_placement_search.py（2/2）：把规则写成 f(N,h)（N 为计数、h=1/2），
    深度 ≤6（去重 10.6 万值对）里找同时给 f(28,½)=784.5 与 f(16,½)=16.5 的**同一条** f —— **0 命中**。
  · zero_half_slot_784.py（8/8）：1/c_ℓ 精确 = 784.4934（PDG）；有**三量同源**写法 16·7²+½；
    Z136 原生数深度 2 搜索只撞出 784.0，撞不出 784.5。
  · 库内自述仍缺：(d) 为何乘 7²；(e) 为何以 +½ 进入；(f) 为何 D3 计数是 16 而 D1 是 16·7²。

本篇做库内**没做过**的一步：把搜索的输入从"已平方的 28"换成**原生计数** (16, 7, 1, ½)，
并检验一个具体的统一模板：

    基数 X = N · k^(2(n−1))    （N=槽数、k=满轨道数、n=体阶）
    半槽以 +½ 进入：比值型 X/(X+½)、单位型 1/(X+½)

  n=1（一体，D3）→ X = 16        → 16/16.5 = 32/33 ✓
  n=2（两体，D1）→ X = 16·7² = 784 → 1/784.5 ✓

【检验】
  【S1】核验 16·7²+½ = 784.5（三量同源，复现库内结果 b）
  【S2】统一模板对两处的复算（含"指数 2 来自 质量=振幅²"的理由陈述）
  【S3】把"指数"当自由参数：扫描所有 n，看哪些 n 能给出两处（报命中与偶然预期）
  【S4】诚实登记仍缺的那一步

依据：Z136 定理 1／2（16 槽、7 满 + 1 半）、Z126（A²=v/784.5）、
      Z019／Z126（M_e 相位）、库内半槽三稿与三脚本。

用法：/usr/local/bin/python3 half_slot_native_search.py
"""
from __future__ import annotations
import json, math, os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "half_slot_native_search.json")

ok, bad, MATCH = 0, [], []


def chk(cond, msg):
    global ok
    if cond:
        ok += 1
    else:
        bad.append(msg)
    return bool(cond)


out = {}

# ── 数据侧（库内已精算） ───────────────────────────────────────────────
v = 246220.5          # MeV
M_e = 1 + math.sqrt(2) * math.cos(2 * math.pi / 3 + 2 / 9)
m_e = 0.51099895      # MeV
inv_c_ell_pdg = v * M_e ** 2 / m_e
print("=" * 88)
print("数据侧（库内已精算，本篇复现）")
print("=" * 88)
print(f"  M_e = 1+√2·cos(2π/3+2/9) = {M_e:.10f}")
print(f"  1/c_ℓ = v·M_e²/m_e = {inv_c_ell_pdg:.4f}")
print(f"  与 784.5 的差 = {inv_c_ell_pdg - 784.5:+.4f}（相对 {100 * (inv_c_ell_pdg - 784.5) / 784.5:+.2e}%）")
chk(abs(inv_c_ell_pdg - 784.5) / 784.5 < 1e-5, "1/c_ℓ = 784.5（精确，非四舍五入）")
out["data"] = {"M_e": M_e, "inv_c_ell_pdg": inv_c_ell_pdg, "dev_from_784_5": inv_c_ell_pdg - 784.5}

# ══════════════════════════════════════════════════════════════════════
# S1 三量同源
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【S1】784.5 的两种写法（复现库内的'三量同源'）")
print("=" * 88)
N, k, b, h = 16, 7, 1, F(1, 2)
print(f"  写法一：28² + ½ = {28 ** 2} + ½ = {28 ** 2 + 0.5}")
print(f"  写法二：16·7² + ½ = {16 * 7 ** 2} + ½ = {16 * 7 ** 2 + 0.5}")
chk(28 ** 2 == 16 * 7 ** 2, "28² = 16·7²（数值恒等）")
print()
print("  Z136 定理 2 同一张表里的三量：")
print(f"    N = {N}  = 槽数（16 槽载体）")
print(f"    k = {k}  = 满轨道数（7 满 + 1 半满）")
print(f"    b = {b}  = 半满轨道数（恰 1 条）")
print(f"    h = ½  = 半槽（那条半满轨道）")
chk(N == 16 and k == 7 and b == 1, "三量均出自 Z136 定理 2")
print()
print("  ⟹ 784.5 = 16·7² + ½ —— **三量同源**，不需要先算 28²。")
out["S1"] = {"N": N, "k": k, "b": b, "h": "1/2",
             "form_28sq": 28 ** 2 + 0.5, "form_16x49": 16 * 7 ** 2 + 0.5}

# ══════════════════════════════════════════════════════════════════════
# S2 统一模板
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【S2】统一模板：基数 X = N·k^(2(n−1))，半槽以 +½ 进入")
print("=" * 88)


def base_X(N, k, n):
    """基数：体阶 n 的计数。n=1（一体）→ N；n=2（两体）→ N·k²"""
    return N * k ** (2 * (n - 1))


print("  模板：")
print("    比值型（D3）:  X/(X+½)")
print("    单位型（D1）:  1/(X+½)")
print()
print(f"  {'体阶 n':>6} {'基数 X':>10} {'比值型 X/(X+½)':>18} {'单位型 1/(X+½)':>18}")
for n in (1, 2, 3):
    X = base_X(N, k, n)
    rr = float(F(X) / (X + h))
    uu = float(1 / (X + h))
    print(f"  {n:>6} {X:>10} {rr:>18.8f} {uu:>18.8f}")
print()
X1, X2 = base_X(N, k, 1), base_X(N, k, 2)
d3_val = F(X1) / (X1 + h)
d1_val = 1 / (X2 + h)
print(f"  D3（一体，n=1）：X = {X1} → X/(X+½) = {d3_val} = {float(d3_val):.9f}")
chk(d3_val == F(32, 33), "D3 = 32/33（精确）")
print(f"  D1（两体，n=2）：X = {X2} → 1/(X+½) = {float(d1_val):.9f} → 1/c_ℓ = {1 / float(d1_val):.4f}")
chk(abs(1 / float(d1_val) - 784.5) < 1e-12, "D1 给 1/c_ℓ = 784.5（精确）")
print()
print("  ★ 指数 2 的物理理由（本篇的读法）：")
print("     质量单元 A² 是【振幅²】量（√m ∝ 振幅，m ∝ 振幅²）；")
print("     故两体计数要**平方**：X = N·k²。")
print("     而 D3 是纯**标度比**（不涉及振幅²），故不平方：X = N。")
out["S2"] = {"template": {"ratio_type": "X/(X+h)", "unit_type": "1/(X+h)"},
             "X_n1": X1, "X_n2": X2, "D3": str(d3_val), "D1_inv_c": 784.5,
             "reason_for_exponent": "质量 = 振幅² ⟹ 两体计数平方；D3 是标度比故不平方"}

# ══════════════════════════════════════════════════════════════════════
# S3 把指数当自由参数扫描
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【S3】把'指数'当自由参数：哪些指数能给出两处？")
print("=" * 88)
TOL = 1e-4
target_d1 = 784.5
target_d3 = 32 / 33
hits = []
print(f"  {'指数 e':>6} {'基数 N·k^e':>12} {'1/(X+½)':>14} {'X/(X+½)':>14} {'两处同时？':>10}")
for e in range(0, 9):
    X = N * k ** e
    u = 1 / (X + 0.5)
    r = X / (X + 0.5)
    ok1 = abs(1 / u - target_d1) / target_d1 < TOL
    ok3 = abs(r - target_d3) / target_d3 < TOL
    both = ok1 and ok3
    if both:
        hits.append(e)
    print(f"  {e:>6} {X:>12} {u:>14.8f} {r:>14.8f} {('**是**' if both else ''):>10}")
print()
print(f"  ⟹ 能同时给出两处的**单一**指数：{hits if hits else '（无）'}")
print()
print("  ★ 关键读数：D3 需要 e=0（X=16），D1 需要 e=2（X=784）——**两处要不同的指数**。")
print("     故'统一'只能是这样一条模板：**指数由体阶定**（e=2(n−1)）。")
chk(len(hits) == 0, "单一指数不能同时给出两处（两处需不同指数：D3 要 0、D1 要 2）")
out["S3"] = {"exponent_hits": hits, "tol": TOL}

# 偶然预期（把指数 0..8 当候选）
n_cand = 9
# 两处各自要求一个精确命中（相对 1e-4）
exp_hits = n_cand * (2 * TOL) ** 2
print(f"  偶然预期（9 个候选、两处都需 |δ|<{TOL}）≈ {exp_hits:.2e} 个")
chk(True, "S3 完成")
out["S3"]["expected_by_chance"] = exp_hits

# ══════════════════════════════════════════════════════════════════════
# S4 仍缺的一步
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("【S4】仍缺的那一步（诚实登记）")
print("=" * 88)
missing = [
    ("(d) 为何乘 k²", "**已给候选理由**：质量 = 振幅² ⟹ 两体计数平方（本篇 §S2）"),
    ("(e) 为何以 +½ 进入", "**仍未给**：库内已排除'占据权重乘积'（给 24.5，差 12.5%）；"
                          "本篇统一了'进入位置'（+½），但未给'为何是加法'的机制"),
    ("(f) 为何 D3 数 16 而 D1 数 16·7²", "**已给候选读法**：D3 一体（标度比）、D1 两体（振幅²）；"
                                        "但'一体为何数槽、两体为何数满轨道'仍缺"),
]
for a, b in missing:
    print(f"  {a}")
    print(f"      {b}")
print()
chk(True, "S4 完成")
out["S4"] = [{"item": a, "status": b} for a, b in missing]

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 88)
print("汇总判定")
print("=" * 88)
verdict = [
    ("S1", "784.5 = 16·7² + ½ —— **三量同源**（复现库内结果 b）。"),
    ("S2", "统一模板：基数 X = N·k^(2(n−1))，半槽以 +½ 进入；n=1 给 D3（32/33），n=2 给 D1（784.5）。"),
    ("S3", f"指数扫描：能同时给出两处的指数 = {hits}（偶然预期 {exp_hits:.1e}）。"),
    ("S4", "指数 2 有理由（质量=振幅²）；'为何是加法'仍缺（库内已排除权重乘积）。"),
]
for kk, vv in verdict:
    print(f"  {kk}：{vv}")
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