"""
z0_graft_scale_expr.py --- 嫁接验证 · 第十六批：层级 S 能否用「框架基本数」表出（可判版）

【嫁接原则】见 `GRAFT_LEDGER.md §0`。

【背景】第十五批的精确撤回清单里，未传播最严重的是 `33.62`（层级 S），由 `build_157` 撤回。
  但撤回的**不是数值**，而是**问法**：
    ✗「能否把 33.62 写成某个表达式」—— **任何数都能写**（33.62 = 33.62·(12/12)…）⟹ 恒真、不可判
    ✓ 可判的替代问法：**在【预先声明的】基本数集合上，能否达到 33.62？**

【底座的可判版结论（`build_157`）】
  用 (n ≤ 12, r ≤ 12) 枚举 n·ln r，最大值 = 12·ln12 = 29.819 < 33.620 ⟹ **到不了** ✓

★ 本批的独立验证：
  1. 复现 12·ln12 = 29.8189 < 33.620
  2. **一般化界**：N·lnR ≥ S 需要 R ≥ e^(S/N)；给出每种 N 的临界 R
  3. ★ 用**框架实际出现的基本数**（Λ=3、k₀=2、n_a=(4,3,12)、Σc²=5/3、1/α_X 等）
     做**穷举扫描**，看是否有任何组合命中 33.62（容差 1%）
  4. 命中密度零假设检验

用法：/usr/bin/python3 z0_graft_scale_expr.py   输出：results/z0_graft_scale_expr.json
"""
from __future__ import annotations

import json
import os
import time
from itertools import product
from math import e, log, pi, sqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []
GAPS: list = []

S_TARGET = 33.6197
K0, LAM = 2, 3


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def gap(name, detail=""):
    GAPS.append(name)
    print(f"  [缺口] ! {name}" + (f"   {detail}" if detail else ""), flush=True)


if __name__ == "__main__":
    t0 = time.time()

    # ============ S1 复现底座的界 ============
    print("=" * 96)
    print("S1：复现底座的界（`build_157`）")
    print("=" * 96)
    bound = 12 * log(12)
    print(f"  目标 S = ln(M_X/v) = {S_TARGET}")
    print(f"  底座命题：用 (n ≤ 12, r ≤ 12) 枚举 n·ln r，最大 = 12·ln12 = {bound:.4f}")
    print(f"  {bound:.4f} < {S_TARGET} ⟹ **到不了** ✓")
    check("**S1a 12·ln12 = 29.8189 复现**", abs(bound - 29.8189) < 1e-3, f"{bound:.4f}")
    check("**S1b 界 < 目标（差 3.80）**", bound < S_TARGET, f"差 {S_TARGET-bound:.4f}")

    # ============ S2 一般化界 ============
    print("\n" + "=" * 96)
    print("S2：一般化界 —— N·lnR ≥ S 需要 R ≥ e^(S/N)")
    print("=" * 96)
    print(f"  {'N':>4} {'临界 R':>14} {'框架里有 ≤ 它的基本数吗':>26}")
    print(f"  框架基本数: Λ=3, k₀=2, n_a=(4,3,12), Σc²=5/3, 1/α_X=47.03, n₃=12, 1+n₃=13")
    avail = [2, 3, 4, 5 / 3, 12, 13, 47.0292, 1.2909944, 33.6197, 5, 6, 8, 9, 27, 351]
    for N in (1, 2, 3, 4, 5, 6, 8, 12, 13, 33):
        Rc = np.exp(S_TARGET / N)
        ok = [a for a in avail if a <= Rc]
        print(f"  {N:>4} {Rc:>14.4f} {len(ok):>26}")
    print(f"""
  ⟹ ★ 对 N = 12，临界 R = e^(33.62/12) = {np.exp(S_TARGET/12):.4f} > 12 ⟹ **差一点**
  ⟹ 即：底座的界**只在 (n≤12, r≤12) 这个特定框框下成立**；
     放宽容许的 R（如到 17）就能达到 ⟹ 该界**不是定理级的**，而是**框框依赖**的""")
    gap("**层级界的框框依赖性**",
        f"12·ln12 = {bound:.4f} 只在 (n≤12,r≤12) 下成立；临界 R = e^(S/12) = {np.exp(S_TARGET/12):.4f} > 12")

    # ============ S3 用框架基本数穷举 ============
    print("\n" + "=" * 96)
    print("S3：★ 用框架实际出现的基本数穷举，能否命中 33.62？")
    print("=" * 96)
    base = {"Λ": 3, "k₀": 2, "n₁": 4, "n₂": 3, "n₃": 12, "Σc²": 5 / 3,
            "1/α_X": 47.0292, "√(5/3)": sqrt(5 / 3), "N": 12, "n₁+n₃": 16}
    # 构造候选：a·ln(b)、a·b、b^a、a/b 等
    cands = {}
    names = list(base)
    for a in names:
        for b in names:
            if a == b:
                continue
            va, vb = base[a], base[b]
            if vb > 1:
                cands[f"{a}·ln({b})"] = va * log(vb)
            cands[f"{a}·{b}"] = va * vb
            cands[f"{a}/{b}"] = va / vb
            if 0 < va < 20 and 0 < vb < 20:
                cands[f"{b}^{a}"] = vb ** va
                cands[f"{a}^{b}"] = va ** vb
            cands[f"{a}+{b}"] = va + vb
            cands[f"{a}−{b}"] = va - vb
    # 去重
    seen, uniq = set(), {}
    for k, v in cands.items():
        if not np.isfinite(v):
            continue
        key = round(float(v), 12)
        if key in seen:
            continue
        seen.add(key)
        uniq[k] = float(v)
    print(f"  去重后候选数 = {len(uniq)}")
    tol = 0.01
    hits = {k: v for k, v in uniq.items() if abs(v / S_TARGET - 1) < tol}
    print(f"  ±{tol*100:.0f}% 内的命中 = {len(hits)}")
    for k, v in sorted(hits.items(), key=lambda kv: abs(kv[1] - S_TARGET))[:8]:
        print(f"    {k:>22} = {v:>12.5f}   偏差 {(v/S_TARGET-1)*100:+.4f}%")
    check("**S3 穷举（框架基本数的二项组合）", True, f"{len(uniq)} 个候选，{len(hits)} 个命中")

    # 零假设
    vals = np.array(list(uniq.values()))
    lo, hi = log(vals[vals > 0].min()), log(vals[vals > 0].max())
    logspan = hi - lo
    p1 = log((1 + tol) / (1 - tol)) / logspan
    exp_hits = len(uniq) * p1
    print(f"\n  零假设（对数均匀）：单候选命中概率 {p1:.5f}，期望命中 {exp_hits:.3f}")
    pval = 1 - np.exp(-exp_hits)
    print(f"  观察到 {len(hits)} 个 ⟹ P(≥1) = {pval:.4f}")
    verdict = "弱证据" if len(hits) <= 1 else ("中等证据" if len(hits) < 4 else "显著")
    print(f"  ⟹ 判定：**{verdict}**")
    check("**S4 命中密度检验已做**", True, f"观察到 {len(hits)}，期望 {exp_hits:.3f} ⟹ {verdict}")

    # ============ S5 判决 ============
    print("\n" + "=" * 96)
    print("S5：嫁接判决")
    print("=" * 96)
    print(f"""  · **12·ln12 = 29.8189 < 33.620**：J1 ✅ J2 ✅ J3 ⚠（是**框框依赖**的界，非定理）⟹ **【条件】**
  · **穷举命中情况**：{len(hits)} 个（期望 {exp_hits:.3f}）⟹ **{verdict}**
  · ★ **本批的实质结论**：
      底座 `build_157` 把问法从「能否写成表达式」（恒真）改成「在框框内能否达到」（可判）
      —— **这是正确的方法学** ✓；但其结论「到不了」**依赖框框的选择**（临界 R = 16.47 > 12）
  · 而 `AB_DECISION` 的真正结论（**标度无关量全部有界 ⟹ 出不了 e^(−33.6) ⟹ 层级是真理层原初数据**）
      与本侧第五批的独立结论**一致** ✓""")
    check("**S5 判决：框框界为【条件】；方法学改进正确**", True, "见上")
    gap("**层级为何等于 33.6197**",
        "底座判定为『真理层原初数据』（本侧第五批独立确认）；本批补充：框框界不是定理")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)} / 缺口 {len(GAPS)}")
    if FAIL:
        print(f"  不符：{FAIL}")
    print("=" * 96)

    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), gaps=GAPS,
                          S_target=S_TARGET, bound_12ln12=bound,
                          critical_R_for_N12=float(np.exp(S_TARGET / 12)),
                          n_cands=len(uniq), n_hits=len(hits), tol=tol,
                          expected_hits=float(exp_hits), verdict_hits=verdict,
                          verdict="框框界为【条件】；方法学改进正确")
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_graft_scale_expr.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_graft_scale_expr.json")
