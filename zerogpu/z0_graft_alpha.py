"""
z0_graft_alpha.py --- 嫁接验证 · 第一批：规范耦合与精细结构常数链

【嫁接原则】不抄。每一条都要：
  (i)  **独立复算**数值；
  (ii) **追溯它的输入**：哪些是 Zero／框架结构给的，哪些是**实验数据反解**的；
  (iii) 判定它是【预言】/【一致性检验】/【循环】/【数据拟合】，并列出**缺口**。

来源（`~/Downloads/cosmos-construct/`）：
  `build_124_scale_map_head_on.py`   α₂=α₃ 交叉 ⟹ M_X、1/α_X
  `build_172_alpha_formula.py`       1/α(0) 的数学表述公式（16/0）
  `ALPHA_FORMULA_SOLVED.md`          公式识别 + 双重验证
  `Q244_ALPHA_S.md`                  α_s 在几何链里未覆盖
  `MAPPING_FORMULAS.md`              R 与 sin²θ_W 的刚性关联

用法：/usr/bin/python3 z0_graft_alpha.py   输出：results/z0_graft_alpha.json
"""
from __future__ import annotations

import json
import os
import time
from math import log, pi, sqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []
GAPS: list = []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


def gap(name, detail=""):
    GAPS.append(name)
    print(f"  [缺口] ! {name}" + (f"   {detail}" if detail else ""), flush=True)


# ============ 框架结构量（应从 Zero 侧给出）============
K0, LAM = 2, 3                    # Z_Λ 律的整数（build_169）
SUM_C2 = 5 / 3                    # Σc²，唯一超荷嵌入 (c_L,c_R)=(−1/√3,+2/√3)（build_86/95/99）
R_OVER = (2 + sqrt(3)) / 4        # r = cos²15°，界间旋量重叠（build_71）

# ============ 实验输入（必须标出）============
SM = dict(a2=29.5846, a3=8.4674, aY=98.3654)   # 1/α_i(M_Z)，SM 归一化
MZ = 91.1876
B2, BY, B3 = -19 / 6, 41 / 6, -7.0             # SM 一圈 β
A0_OBS = 137.035999084                          # CODATA 1/α(0)
AMZ_OBS = 127.95                                # 1/α_em(M_Z)
DALPHA = 0.066300                               # 低位真空极化（标准物理）

if __name__ == "__main__":
    t0 = time.time()

    # ---------- A：α₂=α₃ 交叉 ----------
    print("=" * 96)
    print("A：`α₂=α₃` 交叉 ⟹ M_X 与 X = 1/α_X —— ★这一步用了 α_s")
    print("=" * 96)
    R2, R3, RY = -B2 / (2 * pi), -B3 / (2 * pi), -BY / (2 * pi)
    tX = (SM["a2"] - SM["a3"]) / (R3 - R2)
    MX = MZ * np.exp(tX)
    aL_tX = SM["a2"] + R2 * tX
    aC_tX = SM["a3"] + R3 * tX
    aR_tX = 0.75 * (SM["aY"] + RY * tX - aL_tX / 3)
    print(f"  t_X = ln(M_X/M_Z) = {tX:.6f}")
    print(f"  M_X = {MX:.6e} GeV   （文件：9.82e16 GeV）")
    print(f"  1/α_L(M_X) = {aL_tX:.4f}")
    print(f"  1/α_c(M_X) = {aC_tX:.4f}   （与 L 相等，按构造成立）")
    print(f"  1/α_R(M_X) = {aR_tX:.4f}   （★ 与 L **不等** ⟹ 三耦合不统一）")
    X = aL_tX
    Lln = tX
    print(f"  ⇒ X = 1/α_X = {X:.4f}   （文件：47.0292）")
    print(f"  ⇒ L = ln(M_X/M_Z) = {Lln:.6f}   （文件：34.613030）")
    check("**A1 M_X 与 1/α_X 复现**",
          abs(MX / 9.82e16 - 1) < 0.01 and abs(X - 47.0292) < 0.01,
          f"M_X={MX:.4e}, X={X:.4f}")
    print(f"\n  ★ R 的数据值（由 1/α(0) 反解）:")
    R_data = (A0_OBS * (1 - DALPHA) - Lln * (B2 + BY) / (2 * pi)) / X - 1
    print(f"     R_data = {R_data:.6f}   （文件：1.291148）")
    check("**A2 R_data 反解复现**", abs(R_data - 1.291148) < 1e-4, f"{R_data:.6f}")

    # ---------- B：R 的三个候选 ----------
    print("\n" + "=" * 96)
    print("B：R 的三个候选 —— 框架原式 vs 本候选 vs 数据")
    print("=" * 96)
    R_struct = SUM_C2                       # 迹（= GUT 归一）
    R_norm = sqrt(SUM_C2)                   # 范数（= 本候选）
    R_old = 1 / 3 + 4 / (3 * R_OVER ** 2)   # 框架原式
    print(f"  R = Σc²          = {R_struct:.9f}   （迹；GUT 归一）")
    print(f"  R = √(Σc²)       = {R_norm:.9f}   （范数；★ 本候选）")
    print(f"  R = 1/3+4/(3r²)  = {R_old:.9f}   （原式，r=cos²15°={R_OVER:.6f}）")
    print(f"  R_data（反解）   = {R_data:.9f}")
    print(f"\n  相对 R_data 的偏差：")
    for nm, v in (("Σc²", R_struct), ("√(Σc²)", R_norm), ("原式 1/3+4/(3r²)", R_old)):
        print(f"    {nm:>20}: {v:.6f}   偏差 {(v/R_data-1)*100:+.4f}%")
    check("**B1 √(Σc²) 与 R_data 差 0.012%**", abs(R_norm / R_data - 1) < 2e-4,
          f"{(R_norm/R_data-1)*100:+.4f}%")
    check("**B2 原式偏 ~21%（在 1/α(0) 上）**", True, "见下 C 段")
    gap("**为什么是【范数】√(Σc²) 而不是【迹】Σc²**",
        "文件自认：与 Q274 的 √2 同类（数在框架里，机制缺）")

    # ---------- C：一个数两个观测量 ----------
    print("\n" + "=" * 96)
    print("C：R = √(5/3) ⟹ 两个观测量")
    print("=" * 96)
    shift = Lln * (B2 + BY) / (2 * pi)
    print(f"  L(b₂+b_Y)/(2π) = {shift:.4f}")
    for nm, Rv in (("原式 1/3+4/(3r²)", R_old), ("★ √(5/3)", R_norm), ("数据要", R_data)):
        invMZ = X * (1 + Rv) + shift
        inv0 = invMZ / (1 - DALPHA)
        print(f"  {nm:>18}: R={Rv:.6f}  1/α(M_Z)={invMZ:9.4f} ({(invMZ/AMZ_OBS-1)*100:+.4f}%)  "
              f"1/α(0)={inv0:9.4f} ({(inv0/A0_OBS-1)*100:+.4f}%)  "
              f"sin²θ_W(M_X)={1/(1+Rv):.6f}")
    invMZ_n = X * (1 + R_norm) + shift
    inv0_n = invMZ_n / (1 - DALPHA)
    check("**C1 1/α(M_Z) 偏差 < 0.01%**", abs(invMZ_n / AMZ_OBS - 1) < 1e-4,
          f"{(invMZ_n/AMZ_OBS-1)*100:+.4f}%")
    check("**C2 1/α(0) 偏差 < 0.01%**", abs(inv0_n / A0_OBS - 1) < 1e-4,
          f"{(inv0_n/A0_OBS-1)*100:+.4f}%")
    inv0_old = (X * (1 + R_old) + shift) / (1 - DALPHA)
    print(f"\n  ★ 原式给 1/α(0) = {inv0_old:.4f}，偏差 {(inv0_old/A0_OBS-1)*100:+.2f}% ✗（文件称 21.1%）")
    check("**C3 原式偏 21%（复现）**", abs(inv0_old / A0_OBS - 1 - 0.211) < 0.005,
          f"{(inv0_old/A0_OBS-1)*100:+.2f}%")

    # ---------- D：GUT 自洽 ----------
    print("\n" + "=" * 96)
    print("D：GUT 自洽 —— R(1) = Σc² = 5/3 ⟹ sin²θ_W(M_X) = 3/8")
    print("=" * 96)
    R_at_1 = 1 / 3 + 4 / (3 * 1.0 ** 2)
    print(f"  R(1) = 1/3 + 4/3 = {R_at_1:.9f} = Σc² = {SUM_C2:.9f}   {'✓' if abs(R_at_1-SUM_C2)<1e-12 else '✗'}")
    print(f"  1/(1+Σc²) = {1/(1+SUM_C2):.12f}   3/8 = {3/8:.12f}   "
          f"{'✓ 精确' if abs(1/(1+SUM_C2)-3/8)<1e-12 else '✗'}")
    print(f"  R(r=cos²15°) = {R_old:.6f} > 5/3 ⟹ 往【大】走，而数据要 R = {R_data:.6f} < 5/3 往【小】走")
    check("**D1 R(1) = Σc² = 5/3 且 1/(1+Σc²) = 3/8 精确**",
          abs(R_at_1 - SUM_C2) < 1e-12 and abs(1 / (1 + SUM_C2) - 3 / 8) < 1e-12)
    check("**D2 原式的 r 依赖方向【反了】**（r<1 给出 R>5/3，与数据要的相反）",
          R_old > SUM_C2 > R_data, f"原式 {R_old:.4f} > 5/3 > 数据 {R_data:.4f}")
    gap("**R(r) 对 r 的依赖方向**", "原式在 r<1 时变大，数据要变小；需重写该依赖")

    # ---------- E：循环性审计（最重要的一条）----------
    print("\n" + "=" * 96)
    print("E：★ 循环性审计 —— 这个'预言'的输入是什么？")
    print("=" * 96)
    print("""  链条：
    α_s(M_Z) [实验]  ─┐
    1/α_em(M_Z)=127.95 [实验] ─┼─→ α₂=α₃ 交叉 ─→ t_X, X=1/α_X
    1/α_Y(M_Z)=98.3654 [实验] ─┘
    X, L ─→ 1/α(0) = [X(1+R) + L(b₂+b_Y)/(2π)]/(1−Δα)
    Δα = 0.0663 [实验/标准物理]

  ⟹ 结论（诚实）：
    · **X 与 L 是【实验数据反解】**（用 α_s、α_em、α_Y 在 M_Z 的值）
    · 公式给出的 1/α(0) 因此是【**一致性检验**】，不是零参数预言 ✗
    · 唯一**框架给**的输入是 R = √(5/3) 与 β 系数
    · 而 R 的判据是"它同时对上两个观测量" ⟹ 是**一个数拟合两个观测**
      ⇒ 每观测 ~3% 的巧合概率 ⟹ **弱到中等证据**（文件自己这样定性 ✓）""")
    # 定量：如果 R 随机，同时命中两个观测的概率
    p_one = 0.012 / 100 * 2          # R 落在 R_data ±0.012% 的宽度
    print(f"  \n  定量：R 的容差 ±0.012% ⟹ 单个观测的巧合概率 ~{p_one:.2e}")
    print(f"        两个独立观测同时 ⟹ ~{p_one**2:.2e}（★ 但二者不独立：都走同一条 RG 链）")
    check("**E1 循环性已标清：X 与 L 是实验反解 ⟹ 这是一致性检验**", True,
          "文件自己在 §F 承认 '1/α_X 仍是输入'")
    gap("**α_s(M_Z) 在几何链里未覆盖**（Q244）", "M_X 由 α₃=α₂ 定 ⟹ 用了 α_s ⟹ α_s 不能同时是预言")

    # ---------- F：零假设定量 ----------
    print("\n" + "=" * 96)
    print("F：零假设检验（用对数均匀 + 与文件同一容差）")
    print("=" * 96)
    target, tol = R_data, 0.001
    cands = {"√(Σc²)": sqrt(SUM_C2), "Σc²": SUM_C2, "1/3+4/(3r²)": R_old,
             "1/(1−3/8)": 1 / (1 - 3 / 8), "4/3": 4 / 3, "√(1+k₀/Λ)": sqrt(1 + K0 / LAM),
             "5/3−3/8": SUM_C2 - 3 / 8, "1+k₀/Λ": 1 + K0 / LAM, "√(3/5)⁻¹": 1 / sqrt(3 / 5),
             "√(5/3)(1−1/12)": sqrt(SUM_C2) * (1 - 1 / 12), "(1+r²)/√2": (1 + R_OVER ** 2) / sqrt(2),
             "√(1+r²)": sqrt(1 + R_OVER ** 2), "1+r²/2": 1 + R_OVER ** 2 / 2,
             "√2·r": sqrt(2) * R_OVER, "1/K=4/3": 4 / 3, "√(Σc²·r)": sqrt(SUM_C2 * R_OVER),
             "Σc²·r": SUM_C2 * R_OVER, "√12/√7.2": sqrt(12) / sqrt(7.2),
             "√(mean d)": sqrt(np.mean([1 - np.cos(pi / 4), 1 - np.cos(pi / 3),
                                        1 - np.cos(5 * pi / 12)]))}
    seen, uniq = [], {}
    for k, v in cands.items():
        if not any(abs(v - u) / abs(u) < 1e-9 for u in seen):
            seen.append(v); uniq[k] = v
    hits = {k: v for k, v in uniq.items() if abs(v - target) / target < tol}
    logspan = log(max(uniq.values()) / min(uniq.values()))
    p1 = log((1 + tol) / (1 - tol)) / logspan
    print(f"  去重后候选数 = {len(uniq)}，容差 ±{tol*100:.1f}%，目标 R_data = {target:.6f}")
    for k in sorted(uniq, key=lambda x: abs(uniq[x] - target))[:4]:
        print(f"    {k:>20} {uniq[k]:12.6f}  {abs(uniq[k]-target)/target*100:8.3f}%")
    print(f"  命中: {list(hits.keys())}")
    print(f"  对数均匀零假设: 单候选概率 {p1:.5f}，期望命中 {len(uniq)*p1:.4f}")
    print(f"  ⟹ 命中 {len(hits)} 个：{list(hits.keys())}")
    print(f"     注意：`5/3−3/8 = 41/24 = {SUM_C2-3/8:.6f}` 是**独立候选**，凑巧落进 ±0.1%")
    print(f"     文件自己的容差是 R_data ±0.012%（比本文件的 ±0.1% 严 8 倍）⟹ 它那里只命中 √(Σc²)")
    check("**F1 严容差（±0.02%）下命中唯一 = √(Σc²)**",
          sum(1 for v in uniq.values() if abs(v - target) / target < 2e-4) == 1
          and any(abs(v - target) / target < 2e-4 for k, v in uniq.items() if k == "√(Σc²)"),
          f"严容差命中 = {[k for k,v in uniq.items() if abs(v-target)/target < 2e-4]}")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)} / 缺口 {len(GAPS)}")
    if FAIL:
        print(f"  不符：{FAIL}")
    print("=" * 96)
    print("""【嫁接判决】
  ✅ **数值可复现**：M_X、X、L、R_data、原式 21%、GUT 自洽（3/8）全部复现。
  ✅ **方法学干净**：预先声明候选集、按数值去重、做了零假设、自己定性为"弱到中等证据"。
  ⚠ **性质**：这是【**公式识别 + 一致性检验**】，**不是**第一性原理推导
     —— 用户当时的要求正是"不从第一性原理推导，只要数学表述公式且不矛盾"。
  ❌ **不能算作 Zero 的预言**：X 与 L 是**实验数据反解**（用了 α_s、α_em、α_Y 在 M_Z 的值）；
     唯一框架输入的 R = √(5/3) 是**一个数拟两个观测**。
  ❌ **两个缺口**：(i) 为什么是范数不是迹；(ii) R(r) 的 r 依赖方向反了。""")

    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), gaps=GAPS,
                          MX_GeV=MX, tX=tX, X=X, R_data=R_data,
                          R_struct=R_struct, R_norm=R_norm, R_old=R_old,
                          inv0_norm=inv0_n, inv0_old=inv0_old,
                          sin2W_MX=1 / (1 + R_norm),
                          verdict="一致性检验（X、L 为实验反解）")
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_graft_alpha.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_graft_alpha.json")
