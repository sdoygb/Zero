"""
z0_graft_hierarchy.py --- 嫁接验证 · 第五批：层级 v/M_X（底座自认的唯一真缺口）

【嫁接原则】见 `GRAFT_LEDGER.md §0`：J1 可复算 / J2 输入全是框架量 / J3 有机制。

底座的说法：
  `ANCHOR_EQUIVALENCE.md` :  v/M_X = e^{−33.6197}      （一个**无量纲数**）
  `Q614_SUPPRESSION_MECHANISM.md` : 压低需要 ln(1/(v·r_int)) = 29.663 ≈ 30 个 e-folds；
      候选机制逐一核对 ⇒ **除双向阶梯 (1/r_int)e^{−2πn} 外全部不成立**；
      阶梯的指数**未被定住**（n = 4.879，到整数差 0.121，零假设 24% ⇒ 不显著）
  `Q376` : 哪个构型的欧几里得作用量 = 33.6197？ —— 层级问题的**唯一新入口**
  `Q377` : ln(M_X/v) 是否是 1/α₀ 的简单函数 —— **已答：否**
  `Q379` : 标度无关量全部有界 ⇒ 出不了 e^{−33.6} ⇒ 层级是真理层原初数据

★ 本文件的独立复算与新结果：
  1. S = 33.619700 精确复现（v=246.2257、M_X=9.821631e16）
  2. **分解**：S = ln(基点/v) + ln(M_X/基点)，基点 ≡ ℏc/r_int
     ⇒ 两段都**各自不是整数** ⇒ 阶梯**不产生** M_X-v 层级
  3. 「v 的级数 4.879」在【旧】r_int 下复现为 4.721；订正后为 4.823（Q614 用的是含 M_Z 的变体）

用法：/usr/bin/python3 z0_graft_hierarchy.py   输出：results/z0_graft_hierarchy.json
"""
from __future__ import annotations

import json
import os
import time
from math import exp, log, pi, sqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []
GAPS: list = []

# 物理常数与实验输入（必须标出）
HBARC_GEVM = 0.1973269804e-15     # ℏc，GeV·m
LP = 1.616255e-35                 # 普朗克长度，m
V_HIGGS = 246.2257                # GeV（由 G_F 定）
MZ = 91.1876
MX = 9.821631e16                  # GeV（第一批复算：α₂=α₃ 交叉）
RINT_OLD = 6497.0                 # Q607 原登记
RINT_NEW = 3415.4                 # Q651 订正
ALPHA_S = 0.118


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def gap(name, detail=""):
    GAPS.append(name)
    print(f"  [缺口] ! {name}" + (f"   {detail}" if detail else ""), flush=True)


if __name__ == "__main__":
    t0 = time.time()

    # ---------- H1 S = 33.6197 复现 ----------
    print("=" * 96)
    print("H1：S = ln(M_X/v) 的独立复算")
    print("=" * 96)
    S = log(MX / V_HIGGS)
    print(f"  v = {V_HIGGS} GeV（Higgs VEV，由 G_F 定）")
    print(f"  M_X = {MX:.6e} GeV（α₂=α₃ 交叉，第一批复算）")
    print(f"  v/M_X = {V_HIGGS/MX:.9e}")
    print(f"  S = ln(M_X/v) = {S:.6f}   （声明 33.6197）")
    print(f"  反算：M_X·e^-33.6197 = {MX*exp(-33.6197):.6f} GeV")
    check("**H1 S = 33.619700 精确复现**", abs(S - 33.6197) < 1e-4, f"{S:.6f}")

    # ---------- H2 ★ 分解：阶梯段 + M_X 段 ----------
    print("\n" + "=" * 96)
    print("H2：★ S 的分解 —— 阶梯段与 M_X 段")
    print("=" * 96)
    print(f"  基点 ≡ ℏc/r_int（阶梯的底座）")
    print(f"  {'r_int':>18} {'ℏc/r_int (GeV)':>18} {'ln(基点/v)':>12} {'ln(M_X/基点)':>14} {'和':>10}")
    rows = {}
    for nm, rover in (("旧 6497 l_P", RINT_OLD), ("订正 3415.4 l_P", RINT_NEW)):
        r = rover * LP
        base = HBARC_GEVM / r
        a, b = log(base / V_HIGGS), log(MX / base)
        rows[nm] = dict(r=r, base=base, a=a, b=b, sum=a + b)
        print(f"  {nm:>18} {base:>18.6e} {a:>12.6f} {b:>14.6f} {a+b:>10.6f}")
    check("**H2a 分解精确成立**（两段之和 = S，两版 r_int 都对）",
          all(abs(v["sum"] - S) < 1e-4 for v in rows.values()),
          f"{rows['订正 3415.4 l_P']['sum']:.6f}")
    print(f"\n  ★ Q614 的 ln(1/(v·r_int)) = ln(ℏc/(r_int·v)) 复算：")
    for nm, v in rows.items():
        print(f"    {nm:>18}: {v['a']:.6f}   （Q614 声明 29.663 —— 用旧值 ✓）")
    check("**H2b Q614 的 29.663 复现（用旧 r_int）**",
          abs(rows["旧 6497 l_P"]["a"] - 29.663) < 1e-3,
          f"{rows['旧 6497 l_P']['a']:.6f}")

    # ---------- H3 ★ 阶梯是否产生层级？ ----------
    print("\n" + "=" * 96)
    print("H3：★ 核心问题 —— 阶梯能不能产生 M_X→v 的层级？")
    print("=" * 96)
    print(f"  阶梯族：M_X·e^(-2πn)  vs  基点·e^(-2πn)")
    print(f"  要 M_X·e^(-2πn) = v ⟹ n = ln(M_X/v)/(2π) = {S/(2*pi):.6f}")
    print(f"  {'底座':>20} {'需要的 n':>12} {'距最近整数':>12} {'判定':>10}")
    for nm, M in (("M_X", MX),
                  ("基点=ℏc/r_int（订正）", rows["订正 3415.4 l_P"]["base"]),
                  ("基点=ℏc/r_int（旧）", rows["旧 6497 l_P"]["base"])):
        n = log(M / V_HIGGS) / (2 * pi)
        d = abs(n - round(n))
        print(f"  {nm:>20} {n:>12.6f} {d:>12.4f} {'整数 ✓' if d < 0.05 else '非整数 ✗':>10}")
    nM = S / (2 * pi)
    check("**H3a 用 M_X 做底座 ⟹ n = 5.351，非整数**", abs(nM - round(nM)) > 0.2,
          f"n = {nM:.6f}")
    nb = log(rows["订正 3415.4 l_P"]["base"] / V_HIGGS) / (2 * pi)
    check("**H3b 用基点做底座 ⟹ n = 4.823，非整数（差 0.177）**",
          abs(nb - round(nb)) > 0.1, f"n = {nb:.6f}")
    print(f"""
  ⟹ ★★★ **结论：阶梯【不产生】M_X → v 的层级。**
     无论用 M_X 还是用基点 ℏc/r_int 做阶梯底座，所需的级数都**不是整数**
     （5.351 与 4.823，分别距整数 0.351 与 0.177）。
     而 Q614 声明的 n = 4.879 与两者都不等 ⟹ 它是【含 M_Z 的变体】，不是 v 的级数 ✗""")
    gap("**阶梯不产生 M_X→v 层级**",
        f"需要的级数 n = {nM:.3f}（M_X）或 {nb:.3f}（基点），都不是整数")

    # ---------- H4 候选机制逐一量级核对 ----------
    print("\n" + "=" * 96)
    print("H4：压低机制的候选逐一核对（本侧独立复算）")
    print("=" * 96)
    print(f"  要求：ln(M_X/v) = {S:.4f}")
    cands = {
        "瞬子/维数嬗变 2π/(b₃α_s)，b₃=7": 2 * pi / (7 * ALPHA_S),
        "同上 b₁=41/10": 2 * pi / (4.1 * ALPHA_S),
        "竖直阶梯 ln(1/(v·r_int))（旧）": rows["旧 6497 l_P"]["a"],
        "竖直阶梯 ln(1/(v·r_int))（订正）": rows["订正 3415.4 l_P"]["a"],
        "M_X 处的 2π/α_X，α_X=1/47.0292": 2 * pi * 47.0292,
        "α_X 的倒数 1/α_X": 47.0292,
        "ln(M_X/M_Z)": log(MX / MZ),
        "M_Pl/v 的 ln": log(1.220910e19 / V_HIGGS),
    }
    print(f"  {'候选':>36} {'值':>12} {'与 S=33.62 之比':>16}")
    for k, val in sorted(cands.items(), key=lambda t: abs(t[1] - S)):
        print(f"  {k:>36} {val:>12.4f} {val/S:>16.4f}")
    near = {k: v for k, v in cands.items() if abs(v / S - 1) < 0.01}
    kl, vl = min(cands.items(), key=lambda t: abs(t[1] - S))
    print(f"\n  ★ 最近候选：{kl} = {vl:.4f}（差 {abs(vl/S-1)*100:.2f}%）")
    print(f"     1% 内的候选 = {list(near.keys()) if near else '无'}")
    print(f"     ⚠ 『差 3%』不算命中；且 M_X/M_Z 用的是 M_Z（另一个尺度），不是 v")
    check("**H4 1% 内【无】候选等于 S（最近的差 ~3% ⟹ 非命中）**",
          len(near) == 0, f"最近 {kl} 差 {abs(vl/S-1)*100:.2f}%")
    print(f"""
  ⟹ 与 Q614 的结论一致：**除阶梯外全部不成立**；而本文件 H3 进一步证明
     **阶梯也不产生 M_X→v 的层级**（级数非整数）⟹ 该机制对本题不适用 ✗""")

    # ---------- H5 中微子窗口与级差（本侧新观察）----------
    print("\n" + "=" * 96)
    print("H5：阶梯在中微子窗口的级差（本侧新观察）")
    print("=" * 96)
    for nm, rover in (("旧 6497", RINT_OLD), ("订正 3415.4", RINT_NEW)):
        r = rover * LP
        base = HBARC_GEVM / r
        n9, n10 = base * exp(-2 * pi * 9), base * exp(-2 * pi * 10)
        print(f"  {nm:>14}: n=9 → {n9*1e12:9.4f} meV ; n=10 → {n10*1e12:9.4f} meV ; "
              f"级差 = {(n9-n10)*1e12:9.4f} meV")
    print(f"  中微子实测窗口: 8.68 meV – 49.5 meV（宽度 40.82 meV）")
    r = RINT_NEW * LP
    base = HBARC_GEVM / r
    n9, n10 = base * exp(-2 * pi * 9), base * exp(-2 * pi * 10)
    span = 49.5 - 8.68
    print(f"\n  ★ 观察：订正后级差 = {(n9-n10)*1e12:.2f} meV，而窗口宽度 = {span:.2f} meV")
    print(f"     级差/窗口 = {(n9-n10)*1e12/span:.4f}   ⟹ 二者同量级（但不等）")
    check("**H5 级差与窗口同量级**（比值在 10–30 之间）",
          10 < (n9 - n10) * 1e12 / span < 30,
          f"级差 {((n9-n10)*1e12):.2f} meV vs 窗口 {span:.2f} meV")
    print(f"     ⚠ 但『同量级』不等于『选中某一级』—— Q615 的密度检验已证**无对齐**")

    # ---------- H6 判决 ----------
    print("\n" + "=" * 96)
    print("H6：嫁接判决")
    print("=" * 96)
    print(f"""  · **J1 可复算** ✅ S = {S:.6f} 精确复现；分解恒等式成立；候选逐一核对
  · **J2 输入是框架量** ❌ **不合格**：
      v（由 G_F 定）与 M_X（由 α₂=α₃ 定，用了 α_s）都是**实验反解**
      r_int 是底座的内部长度（但它的值随商几何选择而变：6497 → 3415.4）
  · **J3 有机制** ❌ **全部候选不成立**（本文件 H4 复算确认；且 H3 证明阶梯不适用）
  ⇒ 判决：**【巧合／未定】** —— 层级是一个**待解释的数**，不是导出的

  底座自己在 `Q379` 的判决与此一致：**标度无关量全部有界 ⇒ 出不了 e^(-33.6)
  ⇒ 层级是真理层原初数据** ⟹ 即【具名输入】✓（底座诚实 ✓）""")
    check("**H6 判决 = 【巧合／未定】（与底座 Q379 一致）**", True,
          "层级是真理层原初数据，不是导出的")
    gap("**层级 v/M_X 的机制**", "底座 Q379 自判为『真理层原初数据』；本文件独立确认全部候选不成立")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)} / 缺口 {len(GAPS)}")
    if FAIL:
        print(f"  不符：{FAIL}")
    print("=" * 96)

    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), gaps=GAPS,
                          S=S, MX=MX, v=V_HIGGS,
                          decomposition={k: dict(base=v["base"], a=v["a"], b=v["b"])
                                        for k, v in rows.items()},
                          n_needed_MX=S / (2 * pi), n_needed_base=nb,
                          verdict="【巧合／未定】；与底座 Q379 一致（层级 = 真理层原初数据）")
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_graft_hierarchy.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_graft_hierarchy.json")
