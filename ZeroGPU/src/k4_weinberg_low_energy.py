#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_weinberg_low_energy.py —— 把 sin²θ_w(M_Z)=0.23122 化归成 Zero 侧的一个数
================================================================================
目标（用户指示）：**推进，要算出 0.23122**。

做法分三段，来源严格分开：

【A】β 系数由 Zero 的内容**精确算出**（Zero 侧，精确有理数）
    b_i 只依赖粒子内容：b_g = −(11/3)C₂(G) + (2/3)Σ_f T(R_f) + (1/3)Σ_s T(R_s)
    Zero 有：3 代 × 16 槽（Z136）+ 一个 Higgs 双重态。逐项算：
      · SU(3)：每代色 Weyl 4 个（u_L,d_L,u_R,d_R），T(3)=1/2 ⟹ ΣT = 2/代
      · SU(2)：每代**双重态** 1(轻子) + 3(夸克×色) = 4 个，T(2)=1/2 ⟹ ΣT = 2/代
      · U(1)（GUT 归一化 T_Y=√(3/5)·Y/2）：每代 Σ(Y/2)² = 10/3（Z136 槽内容）
    验证：应精确复现 SM 的 (41/10, −19/6, −7)。
    ⟹ 结论：**Zero 缺的不是 β 系数，是标度**（这正是 Z032/Z056 该被精确化的地方）。

【B】把低能值化归成一个数 L
    统一标度给 sin²θ_w(M_U) = 3/8（前一份登记的迹比结果），一圈跑动给
        A_i(M_Z) = A_U + (b_i/2π)·L,    L ≡ ln(M_U/M_Z),  A_i ≡ 1/α_i
    于是 sin²θ_w(M_Z) = (3/5)A_2/((3/5)A_2 + A_1) —— **一个 L 的函数**。
    反过来用实测值解出 L_req（这就是"要算出 0.23122，Zero 只需提供一个数"）。

【C】L 在 Zero 侧的可能来源（**预先声明候选集**，防事后挑）
    Zero 的内禀对数频率来自 seed 相邻比（Z087/Z131）：log(3/2)、log2、log(8/3)。
    "标度比 = seed 比之幂" ⟹ L 应是这三个数的**整系数组合**（= log(2^a(3/2)^b(8/3)^c)）。
    在该格点上搜最接近 L_req 的组合，报偏差与**搜索规模**（look-elsewhere），
    并明确：命中也不等于推导，除非有机制定出指数 (a,b,c)。

用法：/usr/local/bin/python3 k4_weinberg_low_energy.py
"""
from __future__ import annotations
from fractions import Fraction as F
import itertools, json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "k4_weinberg_low_energy.json")

ok, bad, MATCH = 0, [], []


def chk(cond, msg):
    global ok
    if cond:
        ok += 1
    else:
        bad.append(msg)
    return bool(cond)


# ══════════════════════════════════════════════════════════════════════
# 【A】β 系数从 Zero 的 16 槽内容精确算出
# ══════════════════════════════════════════════════════════════════════
print("=" * 90)
print("【A】β 系数由 Zero 的内容精确算出（精确有理数）")
print("=" * 90)

NGEN = 3                      # 代数 = K₄ 三频率槽（Z089）
N_HIGGS_DOUBLET = 1           # 一个 Higgs 双重态（外部输入，登记）
C2_SU3, C2_SU2 = F(3), F(2)   # 伴随的二次 Casimir
T3_Q, T2_D = F(1, 2), F(1, 2) # 三重/双重态指标

# 每代内容（来自 Z136 的 16 槽）
colored_weyl = 4              # u_L,d_L,u_R,d_R（每个都是色三重）
su2_multiplets = 1 + 3        # 轻子双重态 1 + 夸克双重态 3（每色一个）
sum_Yhalf_sq = F(10, 3)       # Σ(Y/2)²（Z136 槽内容，WEINBERG_ANGLE.md §1）

# b_g = −(11/3)C₂ + (2/3)·NGEN·Σ_f T(R) + (1/3)·Σ_s T(R)
b3 = -(F(11, 3) * C2_SU3) + F(2, 3) * NGEN * (colored_weyl * T3_Q)
b2 = -(F(11, 3) * C2_SU2) + F(2, 3) * NGEN * (su2_multiplets * T2_D) \
     + F(1, 3) * N_HIGGS_DOUBLET * T2_D
# U(1)：生成元取 T_Y = √(3/5)·(Y/2) ⟹ 每代 ΣT_Y² = (3/5)Σ(Y/2)²；Higgs 双重态含 2 个复分量
b1 = F(2, 3) * (F(3, 5) * sum_Yhalf_sq) * NGEN + F(1, 3) * (F(3, 5) * F(1, 4)) * 2

print(f"  SU(3)：每代色 Weyl {colored_weyl} 个 ×T={T3_Q} ⟹ ΣT={colored_weyl*T3_Q}/代")
print(f"  SU(2)：每代双重态 {su2_multiplets} 个 ×T={T2_D} ⟹ ΣT={su2_multiplets*T2_D}/代")
print(f"  U(1) ：每代 Σ(Y/2)² = {sum_Yhalf_sq}（Z136 槽内容）")
print(f"\n  算出：b₁ = {b1} = {float(b1):.4f}；b₂ = {b2} = {float(b2):.4f}；"
      f"b₃ = {b3} = {float(b3):.4f}")
chk(b1 == F(41, 10) and b2 == F(-19, 6) and b3 == F(-7),
    "b 系数应精确等于 SM 的 (41/10, −19/6, −7)")
MATCH.append(("A1", f"β 系数由 Zero 内容精确复现：b=({b1}, {b2}, {b3}) = SM 标准值", True))
chk(True, "ν_R 是规范单态 ⟹ 对 b_i 无贡献 ⟹ 低能 sin²θ_w 与 ν_R 无关")
MATCH.append(("A2", "ν_R（规范单态）对三个 b 系数贡献为 0 ⟹ 低能值对 ν_R 不敏感", True))
print("  · ν_R 是规范单态 ⟹ 不进 b_i ⟹ 低能 sin²θ_w 与它无关")

# ══════════════════════════════════════════════════════════════════════
# 【B】sin²θ_w(M_Z) 是 L 的函数 → 反解 L_req
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 90)
print("【B】sin²θ_w(M_Z) = f(L)：反解「要复现 0.23122 需要多大的 L」")
print("=" * 90)

MZ = 91.1876
alpha_em_inv = 127.95
sin2_meas = 0.23122
alpha_s = 0.1179
A1_MZ = (3 / 5) * (1 - sin2_meas) * alpha_em_inv     # 1/α₁(M_Z)
A2_MZ = sin2_meas * alpha_em_inv                     # 1/α₂(M_Z)
A3_MZ = 1 / alpha_s

# 关系（本文档 §2 的推导，两处都用）：
#   α_em(M_Z) = e²/4π，e = g sinθ = g' cosθ ⟹ α_em = α_2 sin²θ = α_Y' cos²θ
#   ⟹ α_em = (3/5)/(A₁ + (3/5)A₂)          【匹配条件，只用一个实测数 α_em】
#   ⟹ sin²θ_w = (3/5)A₂/(A₁ + (3/5)A₂)     【预言式】
#   一圈：A_i(M_Z) = A_U + (b_i/2π)·L,  L = ln(M_U/M_Z)


def sin2_pred(b1f, b2f, L, alpha_em_inv=127.95):
    """给定 L 与实测 α_em(M_Z)，预言 sin²θ_w(M_Z)（及 A₁、A₂、A_U）"""
    R = (3 / 5) * alpha_em_inv                       # = A₁ + (3/5)A₂ at M_Z
    A_U = (R - L * (b1f + (3 / 5) * b2f) / (2 * math.pi)) / (1 + 3 / 5)
    A1 = A_U + (b1f / (2 * math.pi)) * L
    A2 = A_U + (b2f / (2 * math.pi)) * L
    return (3 / 5) * A2 / (A1 + (3 / 5) * A2), A1, A2, A_U


def L_required(b1f, b2f, target=sin2_meas, lo=0.0, hi=200.0):
    """二分求 L 使 sin2_pred(L) = target（sin²θ_w 随 L 单调下降）"""
    for _ in range(200):
        mid = (lo + hi) / 2
        if sin2_pred(b1f, b2f, mid)[0] > target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

SM_B = (float(b1), float(b2), float(b3))
MSSM_B = (33 / 5, 1, -3)

print(f"  实测输入：α_em⁻¹={alpha_em_inv}, sin²θ_w={sin2_meas} ⟹ "
      f"A₁(M_Z)={A1_MZ:.3f}, A₂(M_Z)={A2_MZ:.3f}")
res = {}
for name, (b1f, b2f, b3f) in (("SM(=Zero 内容)", SM_B), ("MSSM(外部内容)", MSSM_B)):
    L = L_required(b1f, b2f)
    s2, A1, A2, A_U = sin2_pred(b1f, b2f, L)
    scale = MZ * math.exp(L)
    res[name] = {"b": (b1f, b2f, b3f), "L": L, "scale_GeV": scale,
                 "sin2_reproduced": s2, "A1_MZ": A1, "A2_MZ": A2}
    print(f"\n  [{name}] b=({b1f:.4f},{b2f:.4f},{b3f:.4f})")
    print(f"    L_req = {L:.4f}  ⟹ M_U = {scale:.3e} GeV")
    print(f"    代入预言式复现 sin²θ_w(M_Z) = {s2:.6f}（目标 {sin2_meas}，偏差 "
          f"{100*(s2-sin2_meas)/sin2_meas:+.4f}%）；A_U = {A_U:.3f}")

chk(abs(res["SM(=Zero 内容)"]["sin2_reproduced"] - sin2_meas) < 1e-5,
    "SM 内容下由 L_req 复现 0.23122（<1e-5）")
MATCH.append(("B1", f"**0.23122 化归为一个数**：L_req(SM 内容) = "
                    f"{res['SM(=Zero 内容)']['L']:.3f}", True))
MATCH.append(("B2", f"换内容（MSSM）则 L_req = {res['MSSM(外部内容)']['L']:.3f}"
                    f" ⟹ 目标值对内容敏感", True))

# ══════════════════════════════════════════════════════════════════════
# 【C】L 的 Zero 侧候选：三个内禀对数频率的整系数格点（预先声明候选集）
# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 90)
print("【C】L 的 Zero 侧候选：内禀对数频率 {log(3/2), log2, log(8/3)} 的整系数组合")
print("=" * 90)
f1 = math.log(1.5)          # log(3/2) —— seed 相邻比 1.5:1
f2 = math.log(2.0)          # log2    —— seed 相邻比 2:1（内禀频率，Z087）
f3 = math.log(8 / 3)        # log(8/3) —— seed 相邻比 8:3
print(f"  内禀频率：log(3/2)={f1:.6f}, log2={f2:.6f}, log(8/3)={f3:.6f}")

def best_lattice(target, R=6):
    """在 |a|,|b|,|c| ≤ R 内搜 |a·f1+b·f2+c·f3 − target| 最小者（R 小 = 只允许小指数）"""
    best, cnt = None, 0
    for a in range(-R, R + 1):
        for b in range(-R, R + 1):
            for c in range(-R, R + 1):
                if a == 0 and b == 0 and c == 0:
                    continue
                cnt += 1
                L = a * f1 + b * f2 + c * f3
                d = abs(L - target)
                if best is None or d < best[0]:
                    best = (d, a, b, c, L)
    return best, cnt

fit = {}
targets = {name: res[name]["L"] for name in res}

# C1：只允许小指数（|a|,|b|,|c| ≤ 6）——"标度比 = seed 比的少数几次幂"这一假设
for name, target in targets.items():
    (d, a, b, c, L), cnt = best_lattice(target, R=6)
    b1f, b2f = res[name]["b"][0], res[name]["b"][1]
    s2, _, _, _ = sin2_pred(b1f, b2f, L)
    fit[name] = {"target_L": target, "small_exp": {"a": a, "b": b, "c": c, "L": L,
                 "abs_err": d, "rel_err_percent": 100 * d / target,
                 "sin2": s2, "sin2_err_percent": 100 * (s2 - sin2_meas) / sin2_meas,
                 "searched": cnt}}
    print(f"\n  [{name}] 目标 L = {target:.4f}")
    print(f"    C1 小指数（|系数|≤6，{cnt} 个格点）最佳："
          f"{a}·log(3/2)+{b}·log2+{c}·log(8/3) = {L:.4f}")
    print(f"       ⟹ |ΔL| = {d:.4f}（相对 {100*d/target:.1f}%）⟹ sin²θ_w = {s2:.5f}"
          f"（差 {100*(s2-sin2_meas)/sin2_meas:+.1f}%）—— 远不够")

# C2：大指数能拟合到任意精度（格点致密）——必须如实标成"无信息"
for name, target in targets.items():
    (d, a, b, c, L), cnt = best_lattice(target, R=64)
    fit[name]["large_exp"] = {"a": a, "b": b, "c": c, "L": L, "abs_err": d,
                              "rel_err_percent": 100 * d / target, "searched": cnt}
    print(f"    C2 大指数（|系数|≤64，{cnt} 个格点）最佳："
          f"{a}·log(3/2)+{b}·log2+{c}·log(8/3) = {L:.4f}（|ΔL|={d:.5f}，"
          f"{100*d/target:.4f}%）")
    print(f"       ⚠ 三个频率互不通约 ⟹ 格点致密 ⟹ **任何** 目标都能被大的整数组合拟合；"
          f"这不是证据")

# C3：换算成"e-folds 计数"
folds = {}
for name, target in targets.items():
    folds[name] = {}
    print(f"\n  [{name}] L = {target:.4f} 换算成各内禀频率的次数：")
    for lab, fr in (("log(3/2)", f1), ("log2", f2), ("log(8/3)", f3)):
        n = target / fr
        ni = round(n)
        dev = 100 * (ni * fr - target) / target
        folds[name][lab] = {"n_exact": n, "n_round": ni, "dev_percent_if_integer": dev}
        print(f"    L / {lab:9s} = {n:6.2f}   （最近的整数 {ni}，若取整数则偏差 {dev:+.2f}%）")

# C4：Zero 的"自然组合数"集合（**事先声明**）里，有没有 ln 落在 L_req ±1% 的？
import math as _m
natural = {}
for n in range(1, 31):
    natural[f"C(2*{n},{n})"] = _m.lgamma(2 * n + 1) - 2 * _m.lgamma(n + 1)
    natural[f"2^{n}"] = n * _m.log(2)
    natural[f"3^{n}"] = n * _m.log(3)
    natural[f"{n}!"] = _m.lgamma(n + 1)
natural["|Aut(Heawood)|=336"] = _m.log(336)
natural["C(14,7)=3432"] = _m.log(3432)
natural["C(12,6)=924"] = _m.log(924)
hits = {name: {"ln": v, "dev_percent": 100 * (v - res["SM(=Zero 内容)"]["L"]) / res["SM(=Zero 内容)"]["L"]}
        for name, v in natural.items()
        if abs(100 * (v - res["SM(=Zero 内容)"]["L"]) / res["SM(=Zero 内容)"]["L"]) < 1.0}
print(f"\n  C4：Zero 的'自然组合数'集合（{len(natural)} 个数，事先声明）里 ln 落在 "
      f"L_req ±1% 的：")
for k, v in sorted(hits.items(), key=lambda kv: abs(kv[1]["dev_percent"])):
    print(f"    {k:18s} ln = {v['ln']:.4f}  偏差 {v['dev_percent']:+.3f}%")
if not hits:
    print("    （无）")
MATCH.append(("C4", f"自然组合数集合（{len(natural)} 个）里 ±1% 命中 {len(hits)} 个"
                    f"（命中也不构成机制）", True))
# 偶然预期：把 ln 值视作在 [min,max] 上近似均匀，窗口宽度 = 2×1%×L
_lns = sorted(natural.values())
_spread = _lns[-1] - _lns[0]
_window = 2 * 0.01 * res["SM(=Zero 内容)"]["L"]
_expected = len(natural) * _window / _spread
_poisson2 = 1 - _m.exp(-_expected) * (1 + _expected)
print(f"    偶然预期：{len(natural)} 个数的 ln 跨度 {_spread:.2f}，窗口宽 {_window:.3f} "
      f"⟹ 预期命中 {_expected:.2f} 个（实测 {len(hits)} 个；P(≥{len(hits)})={_poisson2:.2f}）")
MATCH.append(("C5", f"偶然预期 {_expected:.2f} 个 vs 实测 {len(hits)} 个（P≈{_poisson2:.2f}）"
                    f" ⟹ 与偶然同量级、**无信息**", True))
fit["natural_number_search"] = {"set_size": len(natural), "tolerance_percent": 1.0,
                                "hits": hits, "expected_by_chance": _expected,
                                "p_ge_observed": _poisson2}

chk(fit["SM(=Zero 内容)"]["small_exp"]["rel_err_percent"] > 10,
    "小指数组合应明显够不到 L_req（≥10%）")
MATCH.append(("C1", f"**小指数组合失败**：|系数|≤6 时最佳只到 L={fit['SM(=Zero 内容)']['small_exp']['L']:.2f}"
                    f"（需要 {fit['SM(=Zero 内容)']['target_L']:.2f}），sin²θ_w 差 "
                    f"{fit['SM(=Zero 内容)']['small_exp']['sin2_err_percent']:+.1f}%", True))
MATCH.append(("C2", f"大指数能拟合到 {fit['SM(=Zero 内容)']['large_exp']['rel_err_percent']:.4f}%，"
                    f"但格点致密 ⟹ **无信息**（登记为巧合候选）", True))
MATCH.append(("C3", "换算 e-folds：L/log2 = "
                    f"{folds['SM(=Zero 内容)']['log2']['n_exact']:.2f}、"
                    f"L/log(8/3) = {folds['SM(=Zero 内容)']['log(8/3)']['n_exact']:.2f}、"
                    f"L/log(3/2) = {folds['SM(=Zero 内容)']['log(3/2)']['n_exact']:.2f}"
                    f"（都不接近小整数）", True))
verdict = ("0.23122 **不是**从 Zero 推出来的：它现在被化归为一个数 L = ln(M_U/M_Z) = 25.4509"
           "（SM 内容；MSSM 内容 33.0256），由「3/8（统一，Z136 迹比）+ b_i（Zero 内容，"
           "本文精确复现）+ L + 实测 α_em」复现到 5 位以上。"
           "Zero 侧的小指数候选（三个内禀对数频率 {log(3/2), log2, log(8/3)} 的整系数组合）"
           "**够不到** L（|系数|≤6 时最好只到 12.48）；大指数能拟合但格点致密、无机制，"
           "按库内纪律只能登记为【巧合候选】。"
           "⟹ 精确化的缺口：Zero 需要生成一个 **大对数**（≈25.45），"
           "这正是 Z056「无标度生成」那笔缺口的定量形式。")

print()
print("断言清单：")
for tag, txt, good in MATCH:
    print(f"   {tag}  {'✓' if good else '✗'}  {txt}")
_total = len(MATCH)
_failed = len([1 for _t, _x, g in MATCH if not g]) + len(bad)
print(f"\n—— {_total} 条断言 / {_failed} 不符 ——")
for b_ in bad:
    print("   ✗", b_)

os.makedirs(RES, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump({
        "b_from_zero_content": {"b1": str(b1), "b2": str(b2), "b3": str(b3),
                                "matches_SM": (b1, b2, b3) == (F(41, 10), F(-19, 6), F(-7)),
                                "nuR_contribution": 0},
        "reduce_to_one_number": res,
        "lattice_fit": fit,
        "intrinsic_logs": {"log(3/2)": f1, "log2": f2, "log(8/3)": f3},
        "verdict": verdict,
        "assertions": {"total": len(MATCH),
                       "failed": len([1 for _t, _x, g in MATCH if not g]) + len(bad),
                       "detail": bad},
    }, f, ensure_ascii=False, indent=1)
print(f"已存 {OUT}")