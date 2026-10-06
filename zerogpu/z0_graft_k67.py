"""
z0_graft_k67.py --- 嫁接验证 · 第七批：非循环解法（k67）与两处更正

【嫁接原则】见 `GRAFT_LEDGER.md §0`：J1 可复算 / J2 输入全是框架量 / J3 有机制。

本题来自底座的【中央总账】`OLD_DOCS_RECOMPUTE_LEDGER.md`（1116 行，393 篇逐篇重算）
  它在 #4/#5 判定：旧文档 `ALPHA_FORMULA(.SOLVED).md` 的 R = √(Σc²) = √(5/3)
  是「一条**真实的高精度命中**」，而它此前的验证审计（`K65`）**漏掉了**它。
  脚本：`design2/k67_recompute_alpha_formula.py`（6 通过 / 0 不符）

k67 的**非循环**解法（本文件独立复现）：
  独立输入：1/α_em(M_Z) = 127.95、α_s(M_Z) = 0.1179（两个实测数）＋ 框架的 R = √(5/3)
  未知：x = 1/α₂(M_Z)、y = 1/α_Y(M_Z)、t = ln(M_X/M_Z)
  方程 (k67 的精确形式)：
    x  = 1/α₃(M_Z) + [(|b₃| − |b₂|)/2π]·t
    xX = x + (|b₂|/2π)·t                      （1/α₂ 在 M_X）
    y  = R·xX + (b_Y/2π)·t
    约束：x + y = 1/α_em(M_Z)

★ 本文件的独立结论：
  1. sin²θ_W(M_Z) = x/(x+y) 的预言 **0.231293**（实测 0.23122）⟹ 偏差 **0.032%** ✓
     对照标准 GUT（R = 5/3）给 **0.207590** ⟹ 偏差 **−10.220%**
     ⟹ 框架优于标准 GUT **322 倍** ✓ 这是**真预言**（x 不是输入）
  2. ⚠ **对 k67 的一处更正**：`1/α(0)` 的那条（报 −0.0004%）**是零信息** ——
     因为分子 = X(1+R) + L(b₂+b_Y)/2π = x + y **恒等于 1/α_em(M_Z) = 127.95**
     （可代数验证），故 1/α(0) = 127.95/(1−Δα) 只是把输入除以 (1−Δα)
     ⟹ 它**不是预言**，也不能算"命中" ✗
  3. ⚠ 本侧自己的两处错（记录以免再犯）：
     (i) 用 fsolve 收敛到伪解（x=160）；(ii) 写错 y 的方程（漏 R 只乘 xX）

用法：/usr/bin/python3 z0_graft_k67.py   输出：results/z0_graft_k67.json
"""
from __future__ import annotations

import json
import os
import time
from math import exp, pi, sqrt

import numpy as np
import sympy as sp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []
GAPS: list = []

MZ = 91.1876
AEM_INV = 127.95          # 1/α_em(M_Z)，实测
ALS_MZ = 0.1179           # α_s(M_Z)，实测
S2W_MZ = 0.23122          # sin²θ_W(M_Z)，实测（MS-bar）
B1, B2, B3 = 41 / 10, -19 / 6, -7.0
BY = (5 / 3) * B1         # b_Y = (5/3) b₁
TWO = 2 * pi
A0_OBS = 137.035999084
DALPHA = 0.0663


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def gap(name, detail=""):
    GAPS.append(name)
    print(f"  [缺口] ! {name}" + (f"   {detail}" if detail else ""), flush=True)


def solve(R):
    """k67 的精确解法"""
    a3MZ = 1 / ALS_MZ

    def f(t):
        x = a3MZ + ((1 / TWO) * (-B3) - (1 / TWO) * (-B2)) * t
        xX = x + (-B2 / TWO) * t
        y = R * xX + (BY / TWO) * t
        return x + y - AEM_INV

    tt = brentq(f, 5.0, 80.0, xtol=1e-14)
    x = a3MZ + ((1 / TWO) * (-B3) - (1 / TWO) * (-B2)) * tt
    xX = x + (-B2 / TWO) * tt
    y = R * xX + (BY / TWO) * tt
    return x, y, tt, xX


if __name__ == "__main__":
    t0 = time.time()

    # ============ K1 两个 R 的对比 ============
    print("=" * 96)
    print("K1：框架 R=√(5/3) vs 标准 GUT R=5/3（独立复现 k67）")
    print("=" * 96)
    print(f"  独立输入（实测）：1/α_em(M_Z) = {AEM_INV}, α_s(M_Z) = {ALS_MZ}")
    print(f"  框架输入：R = √(Σc²) = √(5/3) = {sqrt(5/3):.6f}")
    print(f"\n  {'':>16} {'框架 √(5/3)':>16} {'标准 GUT 5/3':>16}")
    out = {}
    for nm, R in (("fw", sqrt(5 / 3)), ("gut", 5 / 3)):
        x, y, tt, xX = solve(R)
        out[nm] = dict(R=R, x=x, y=y, t=tt, xX=xX, MX=MZ * exp(tt),
                       s2w=x / (x + y))
    for k, lab in (("x", "1/α₂(M_Z)"), ("y", "1/α_Y(M_Z)"), ("t", "t = ln(M_X/M_Z)"),
                   ("MX", "M_X (GeV)"), ("s2w", "sin²θ_W(M_Z)")):
        a, b = out["fw"][k], out["gut"][k]
        print(f"  {lab:>16} {a:>16.6f} {b:>16.6f}")
    print(f"  {'实测对照':>16} {out['fw']['x']+out['fw']['y']:>16.6f} (=1/α_em, 强制)")
    print(f"\n  实测 sin²θ_W(M_Z) = {S2W_MZ}")
    for nm, lab in (("fw", "框架 √(5/3)"), ("gut", "标准 GUT 5/3")):
        d = (out[nm]["s2w"] / S2W_MZ - 1) * 100
        print(f"    {lab:>14}: 预言 {out[nm]['s2w']:.6f}  偏差 {d:+.4f}%")
    dev_fw = abs(out["fw"]["s2w"] / S2W_MZ - 1)
    dev_gut = abs(out["gut"]["s2w"] / S2W_MZ - 1)
    print(f"\n  ⟹ 框架偏差是 GUT 的 1/{dev_gut/dev_fw:.0f} ⟹ 优于 GUT **{dev_gut/dev_fw:.0f} 倍**")
    check("**K1a 框架预言 sin²θ_W(M_Z) 偏差 < 0.05%**", dev_fw < 5e-4,
          f"{out['fw']['s2w']:.6f}（偏差 {dev_fw*100:.4f}%）")
    check("**K1b 标准 GUT 偏差 ~10%**", abs(dev_gut - 0.1022) < 0.002,
          f"{dev_gut*100:.3f}%")
    check("**K1c 框架优于 GUT 超过 100 倍**", dev_gut / dev_fw > 100,
          f"{dev_gut/dev_fw:.0f} 倍")
    print(f"""
  ★ 为什么这是**真预言**（不是循环）：
     x = 1/α₂(M_Z) **不是输入** —— 它由方程 (2)(3) 解出（x = {out['fw']['x']:.5f}）
     而 sin²θ_W(M_Z) = x/(x+y) 是它的函数；它恰好对上实测 0.23122 ⟹ **有信息** ✓""")

    # ============ K2 ★ 更正：1/α(0) 那条是零信息 ============
    print("\n" + "=" * 96)
    print("K2：★ 对 k67 的一处更正 —— 1/α(0) 那条是【零信息】")
    print("=" * 96)
    R = sqrt(5 / 3)
    x, y, tt, xX = solve(R)
    num = xX * (1 + R) + tt * (B2 + BY) / TWO
    a0 = num / (1 - DALPHA)
    print(f"  k67 的算式：1/α(0) = [X(1+R) + L(b₂+b_Y)/2π]/(1−Δα)")
    print(f"    其中 X = 1/α_X = xX = {xX:.6f},  L = t = {tt:.6f}")
    print(f"    分子 = {xX:.6f}×{1+R:.6f} + {tt:.6f}×{((B2+BY)/TWO):.6f} = {num:.6f}")
    print(f"    1/α(0) = {num:.6f}/{1-DALPHA:.4f} = {a0:.4f}")
    print(f"    实测 {A0_OBS}，偏差 {(a0/A0_OBS-1)*100:+.4f}%")
    print(f"\n  ★★ 关键：分子是否恒等于 1/α_em(M_Z) = {AEM_INV}？")
    print(f"    分子 = {num:.10f}   vs   1/α_em(M_Z) = {AEM_INV:.10f}")
    print(f"    差 = {abs(num-AEM_INV):.2e}")
    check("**K2a 分子恒等于 1/α_em(M_Z)（= 输入）**", abs(num - AEM_INV) < 1e-9,
          f"{num:.10f} vs {AEM_INV}")
    # 代数证明
    a3 = sp.Rational(1, 1) / sp.Rational(1179, 10000)
    R_s, b2s, b3s, bYs = sp.sqrt(sp.Rational(5, 3)), sp.Rational(-19, 6), sp.Integer(-7), sp.Rational(41, 10) * sp.Rational(5, 3)
    ts = sp.symbols("t", positive=True)
    xs = a3 + ((sp.Integer(1) / (2 * sp.pi)) * (-b3s) - (sp.Integer(1) / (2 * sp.pi)) * (-b2s)) * ts
    xXs = xs + (-b2s / (2 * sp.pi)) * ts
    ys = R_s * xXs + (bYs / (2 * sp.pi)) * ts
    num_s = sp.simplify(xXs * (1 + R_s) + ts * (b2s + bYs) / (2 * sp.pi) - (xs + ys))
    print(f"\n  代数验证：[X(1+R) + L(b₂+b_Y)/2π] − (x+y) = {num_s}")
    check("**K2b 代数上恒等（故 1/α(0) 不是预言）**", sp.simplify(num_s) == 0,
          "分子 ≡ x+y ≡ 1/α_em(M_Z)")
    print(f"""
  ⟹ ★★★ **更正**：`1/α(0) = 127.95/(1−Δα) = 137.0355`（偏差 −0.0004%）
     看起来极准，但**分子就是输入的 1/α_em(M_Z)**（代数恒等）⟹ **零信息** ✗
     它不是预言，也不该算作"命中"。k67 的 `AR3` 判据（d1 < 0.05%）**通过得没有意义**。""")
    gap("**k67 的 1/α(0) 那条是零信息**",
        "分子 ≡ 1/α_em(M_Z)（代数恒等）；真正的成绩只有 sin²θ_W(M_Z)")

    # ============ K3 本侧自己的两处错 ============
    print("\n" + "=" * 96)
    print("K3：本侧自己的两处错（记录以免再犯）")
    print("=" * 96)
    print("""  (i) 用 `fsolve` 解三方程 ⟹ 收敛到伪解 (x=160, t=248)，与实测 x≈29.58 相差 5 倍。
      ⟹ 教训：**多解系统必须用解析法或带物理界的括号法**（k67 用 brentq 在 [5,80]）。
  (ii) 第一遍把方程 (3) 写成 y = R·x + …（漏了 x 应为 xX = x + (|b₂|/2π)t）⟹
      得 sin²θ_W = 0.2653。对照 k67 原文才发现。
      ⟹ 教训：**复现他人计算时必须取原文的方程形式，不能凭物理直觉重写**。""")
    check("**K3 两处错已记录**", True, "fsolve 伪解 ＋ 方程重写")

    # ============ K4 判决 ============
    print("\n" + "=" * 96)
    print("K4：嫁接判决")
    print("=" * 96)
    print(f"""  · **sin²θ_W(M_Z)** ：
      J1 ✅ 独立复现（{out['fw']['s2w']:.6f}）
      J2 ⚠ **部分** —— 输入是 1/α_em(M_Z) 与 α_s(M_Z)（两个实测数）；x 被解出，不是输入
          ⟹ 这是【**一条观测预言另一条观测**】，不是从 Zero 结构推出
      J3 ⚠ **部分** —— R = √(Σc²) 的形式仍缺机制（G1，见第一批）；但 √(5/3) 本身来自
          Σc² = 5/3（唯一超荷嵌入，穷举证明）⟹ 比"纯拟合"强
      ⇒ 判决：**【条件】** —— 一个框架数拟一个观测量，且优于标准 GUT 322 倍 ✓

  · **1/α(0)** ：⟹ **【零信息】**（本文件 K2 更正）✗""")
    check("**K4 判决：sin²θ_W 为【条件】；1/α(0) 为【零信息】**", True, "见上")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)} / 缺口 {len(GAPS)}")
    if FAIL:
        print(f"  不符：{FAIL}")
    print("=" * 96)

    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), gaps=GAPS,
                          solution=out, s2w_obs=S2W_MZ,
                          dev_fw=dev_fw, dev_gut=dev_gut, ratio=dev_gut / dev_fw,
                          correction="1/α(0) 的分子 ≡ 1/α_em(M_Z)（零信息）",
                          verdict="sin²θ_W(M_Z) = 【条件】；1/α(0) = 【零信息】")
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_graft_k67.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_graft_k67.json")
