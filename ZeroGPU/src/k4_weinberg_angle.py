#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_weinberg_angle.py —— 弱混合角 sin²θ_w：载波迹范数之比 + 一圈跑动
================================================================================
两个部分，**来源分开标注**：

【A】Zero 侧（精确有理数，纯结构）：
    用 Z136 的 16 槽内容（每代：色单态顶点 0 + 色顶点 1,2,3；同位旋 ↑/↓；手征 L/R），
    取每槽的 (T₃, Y, Q)（Y 取 Q = T₃ + Y/2 约定，即 e_R 的 Y=−2），算三个迹：
        Σ T₃²,  Σ (Y/2)²,  Σ T₃·(Y/2),  Σ Q²
    判定①（完备性）：Σ T₃·(Y/2) = 0 —— 无非对角混合 ⟹ 迹公式干净；
    判定②（主结果）：sin²θ_w|_统一标度 = ΣT₃² / (ΣT₃² + Σ(Y/2)²) = **3/8**（精确）。

    这个 3/8 的**唯一非 Zero 输入**是：「T₃ 与 Y/2 两个生成元的**迹范数相等**」
    —— 即统一群里的同一归一化（= 库里 Z086／Z135 §判定2 第③项「连续载体 U(1)」那笔输入）。
    本脚本把该条件显式写成归一化因子，并给出**反事实**：换成别的归一化，「3/8」就变。

【B】外部物理侧（一圈 RGE，标准模型/MSSM 系数 + 实测 α_i(M_Z)）：
    把测得的三个耦合跑到高能，看 α₁ 与 α₂ 的交点、以及 α₃ 是否也在那儿 ——
    这一步**不是 Zero 的推导**（Zero 无普适 β 函数：Z032；无标度生成：Z056），
    它的用途只是说明「3/8」与「观测 0.23122」之间隔着什么。

用法：/usr/local/bin/python3 k4_weinberg_angle.py
"""
from __future__ import annotations
from fractions import Fraction as F
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
OUT = os.path.join(RES, "k4_weinberg_angle.json")

ok = 0
bad = []


def chk(cond, msg):
    global ok
    if cond:
        ok += 1
    else:
        bad.append(msg)
    return bool(cond)


MATCH = []

# ══════════════════════════════════════════════════════════════════════
# 【A】Zero 侧：16 槽上的三个迹
# ══════════════════════════════════════════════════════════════════════
# 每代 16 槽（Z136 §4）：色单态顶点 v=0（轻子 4 槽）+ 色顶点 v=1,2,3（夸克 4 槽 ×3 色）
# 字段： (名字, 色重数, T₃, Y)   约定 Q = T₃ + Y/2
SM15 = [
    # 色单态（轻子）——顶点 0
    ("nu_L", 1, F(1, 2), F(-1)),
    ("e_L",  1, F(-1, 2), F(-1)),
    ("e_R",  1, F(0),    F(-2)),
    # 色三重（夸克）——顶点 1,2,3
    ("u_L",  3, F(1, 2), F(1, 3)),
    ("d_L",  3, F(-1, 2), F(1, 3)),
    ("u_R",  3, F(0),    F(4, 3)),
    ("d_R",  3, F(0),    F(-2, 3)),
]
# Z136 的补位槽（补出来的那个 = ν_R）：T₃ = 0, Y = 0
NUR = ("nu_R", 1, F(0), F(0))

def traces(fields):
    """返回 (ΣT₃², ΣY'², ΣT₃Y', ΣQ²)，Y' = Y/2"""
    sT2 = sum(m * t3 * t3 for _, m, t3, _ in fields)
    sY2 = sum(m * (y / 2) * (y / 2) for _, m, _, y in fields)
    sTY = sum(m * t3 * (y / 2) for _, m, t3, y in fields)
    sQ2 = sum(m * (t3 + y / 2) ** 2 for _, m, t3, y in fields)
    return sT2, sY2, sTY, sQ2

def sin2(sT2, sY2):
    """sin²θ_w = ΣT₃²/(ΣT₃²+ΣY'²)（迹范数相等归一化下）"""
    return sT2 / (sT2 + sY2)

ALL16 = SM15 + [NUR]
sT2, sY2, sTY, sQ2 = traces(ALL16)

print("=" * 88)
print("【A】Zero 侧：16 槽上的迹（精确有理数）")
print("=" * 88)
print(f"  槽数（含色重数） = {sum(m for _, m, _, _ in ALL16)}")
print(f"  Σ T₃²            = {sT2}   （= {float(sT2)}）")
print(f"  Σ (Y/2)²         = {sY2}   （= {float(sY2)}）")
print(f"  Σ T₃·(Y/2)       = {sTY}   ← 完备性判据（应为 0）")
print(f"  Σ Q²             = {sQ2}   （= {float(sQ2)}）")

chk(sum(m for _, m, _, _ in ALL16) == 16, "16 槽（含色重数）")
chk(sTY == 0, "Σ T₃·(Y/2) 应为 0（完备多重态 ⟹ 无非对角混合）")
chk(sQ2 == sT2 + sY2 + 2 * sTY, "Q = T₃ + Y/2 ⟹ ΣQ² 恒等式")
MATCH.append(("A1", "Σ T₃·(Y/2) = 0：载体是完备多重态（迹公式干净）", sTY == 0))

s_w = sin2(sT2, sY2)
chk(s_w == F(3, 8), "sin²θ_w 应为 3/8")
chk(sT2 == 2 and sY2 == F(10, 3) and sQ2 == F(16, 3),
    "三迹应为 2、10/3、16/3")
print(f"\n  ⟹ sin²θ_w|_统一标度 = {sT2} / ({sT2} + {sY2}) = {s_w} = {float(s_w):.6f}")
MATCH.append(("A2", "sin²θ_w = ΣT₃²/(ΣT₃²+ΣY'²) = 3/8（精确）", s_w == F(3, 8)))
chk(s_w == sT2 / sQ2, "与另一种写法 ΣT₃²/ΣQ² 一致（因 ΣT₃Y'=0）")
MATCH.append(("A3", "3/8 = ΣT₃²/ΣQ²（两写法相等，靠正交性）", s_w == sT2 / sQ2))

# 反事实 1：去掉 ν_R（15 槽）
sT2_15, sY2_15, _, _ = traces(SM15)
chk(sin2(sT2_15, sY2_15) == F(3, 8),
    "去掉 ν_R 仍是 3/8（ν_R 的 T₃=Y=0 ⟹ 对迹无贡献）")
MATCH.append(("A4", "**3/8 对 ν_R 不敏感**（去掉它也一样）⟹ 它不能当 ν_R 的证据", True))

# 反事实 2：只用左手 8 槽
LEFT8 = [f for f in SM15 if not f[0].endswith("_R")] + [NUR]
sT2_L, sY2_L, _, _ = traces(LEFT8)
chk(sin2(sT2_L, sY2_L) == F(3, 4), "只取左手扇区应给 3/4（≠ 3/8）")
MATCH.append(("A5", "只取左手扇区（8 槽）给 3/4 ⟹ **3/8 需要两个手征扇区**（完备性）",
              sin2(sT2_L, sY2_L) == F(3, 4)))
print(f"  · 去掉 ν_R（15 槽）：sin²θ_w = {sin2(sT2_15, sY2_15)}（不变）")
print(f"  · 只手征 L（8 槽）：sin²θ_w = {sin2(sT2_L, sY2_L)}（≠ 3/8，说明需要 R 扇区）")

# 反事实 3：归一化换了（把 Y 当生成元而不是 Y/2；等价于迹范数不等）
def sin2_wrongnorm(fields, k):
    """生成元取 (T₃, k·Y/2)：k≠1 表示两生成元归一化不相等"""
    sT2 = sum(m * t3 * t3 for _, m, t3, _ in fields)
    sY2 = sum(m * (k * y / 2) ** 2 for _, m, _, y in fields)
    return sT2 / (sT2 + sY2)

alt = {str(k): str(sin2_wrongnorm(ALL16, k)) for k in (F(1), F(2), F(1, 2), F(3, 5))}
chk(sin2_wrongnorm(ALL16, F(1)) == F(3, 8) and sin2_wrongnorm(ALL16, F(2)) != F(3, 8),
    "换归一化因子 k≠1 则 3/8 不再成立")
MATCH.append(("A6", "反事实：归一化因子 k≠1 ⟹ 「3/8」变成别的数（k=2 给 3/23）",
              sin2_wrongnorm(ALL16, F(2)) == F(3, 23)))
print(f"  · 归一化因子 k=1 → 3/8；k=2 → {sin2_wrongnorm(ALL16, F(2))}；k=1/2 → "
      f"{sin2_wrongnorm(ALL16, F(1,2))}；k=3/5 → {sin2_wrongnorm(ALL16, F(3,5))}")

# ── A7：反常抵消（这才是 Z090「关系②」的真实依据，见 NUR_VERIFICATION §1 第 6 行）──
# ⚠ 约定坑（本脚本第一版就栽在这里）：反常必须按**左手 Weyl 等价表示**算 ——
#   右手场要取共轭，即 Y → −Y。若把右手场也当成 +Y，会得到
#   [SU(3)²U(1)] = 4、[U(1)³] = −32/9 的**假不抵消**。
def anomalies(fields_dict):
    """fields_dict: name -> (色重数, T₃, Y, 是否双重态成员, 是否右手)
       返回 ([SU(2)²U(1)], [SU(3)²U(1)], [U(1)], [U(1)³])，全部按左手等价表示"""
    def yeff(rec):
        m, _t3, y, _dbl, isR = rec
        return -y if isR else y
    A_SU2sq_U1 = sum(m * yeff(r) for m, r in ((v[0], v) for v in fields_dict.values()) if r[3])
    A_SU3sq_U1 = sum(m * yeff(r) for m, r in ((v[0], v) for v in fields_dict.values()) if r[0] == 3)
    A_U1 = sum(m * yeff(r) for m, r in ((v[0], v) for v in fields_dict.values()))
    A_U1cube = sum(m * yeff(r) ** 3 for m, r in ((v[0], v) for v in fields_dict.values()))
    return A_SU2sq_U1, A_SU3sq_U1, A_U1, A_U1cube

per_gen = {                      # (色重数, T₃, Y, 是否双重态成员, 是否右手)
    "nu_L": (1, F(1, 2), F(-1), True, False), "e_L": (1, F(-1, 2), F(-1), True, False),
    "e_R": (1, F(0), F(-2), False, True), "nu_R": (1, F(0), F(0), False, True),
    "u_L": (3, F(1, 2), F(1, 3), True, False), "d_L": (3, F(-1, 2), F(1, 3), True, False),
    "u_R": (3, F(0), F(4, 3), False, True), "d_R": (3, F(0), F(-2, 3), False, True),
}
a2u1, a3u1, a1, a1c = anomalies(per_gen)
print(f"\n  反常（每代，左手等价表示）：[SU(2)²U(1)]={a2u1}  [SU(3)²U(1)]={a3u1}  "
      f"[U(1)]={a1}  [U(1)³]={a1c}")
chk(a2u1 == 0 and a3u1 == 0 and a1 == 0 and a1c == 0, "四条反常都应抵消")
MATCH.append(("A7", "16 槽内容满足四条规范反常抵消（左手等价表示下全为 0）",
              (a2u1, a3u1, a1, a1c) == (0, 0, 0, 0)))
# 反事实：把 Y_L(夸克) 换成别的值 ⟹ [SU(2)²U(1)] 不再抵消
bad_gen = dict(per_gen)
bad_gen["u_L"] = (3, F(1, 2), F(1, 2), True, False)
bad_gen["d_L"] = (3, F(-1, 2), F(1, 2), True, False)
b2u1, *_ = anomalies(bad_gen)
chk(b2u1 != 0, "改掉 Y_L(夸克) ⟹ [SU(2)²U(1)] 不再抵消")
MATCH.append(("A8", f"反事实：Y_L(夸克)≠−Y_L(轻子)/N_c ⟹ [SU(2)²U(1)]={b2u1}≠0"
                    f"（故关系②由**反常抵消**强制，而非「色重数」独立给出）", b2u1 != 0))
print(f"  · 反事实：Y_L(夸克) 取 1/2 ⟹ [SU(2)²U(1)] = {b2u1}（不抵消）")
# 附带：ν_R 对反常的贡献为 0
chk(all(a == 0 for a in anomalies({k: v for k, v in per_gen.items() if k == "nu_R"}))
    , "ν_R 单独对各条反常贡献 0")
MATCH.append(("A9", "ν_R（Y=T₃=0）对四条反常贡献全 0（它不参与定超荷）", True))



# ══════════════════════════════════════════════════════════════════════
# 【B】外部物理：一圈 RGE（不是 Zero 推导，见文件头）
# ══════════════════════════════════════════════════════════════════════
MZ = 91.1876                 # GeV
alpha_em_inv = 127.95        # α_em^{-1}(M_Z)
sin2_meas = 0.23122          # sin²θ_w(M_Z), MS-bar
alpha_s = 0.1179             # α_s(M_Z)

aY_inv = (1 - sin2_meas) * alpha_em_inv          # α_Y^{-1} = cos²θ_W · α_em^{-1}
a1_inv = F(3, 5) * F(aY_inv).limit_denominator(10**6) if False else (3 / 5) * aY_inv
a2_inv = sin2_meas * alpha_em_inv
a3_inv = 1 / alpha_s

print()
print("=" * 88)
print("【B】外部物理：实测耦合跑到高能（一圈 RGE，系数为 SM / MSSM 标准值）")
print("=" * 88)
print(f"  α_em⁻¹(M_Z)={alpha_em_inv}, sin²θ_w(M_Z)={sin2_meas}, α_s(M_Z)={alpha_s}")
print(f"  ⟹ α₁⁻¹={a1_inv:.2f}, α₂⁻¹={a2_inv:.2f}, α₃⁻¹={a3_inv:.2f}"
      f"（α₁=(5/3)α_Y 归一化）")

def run(b, t):
    """一圈：1/α_i(μ) = 1/α_i(M_Z) - (b_i/2π) t, t = ln(μ/M_Z)"""
    return tuple(a - (bi / (2 * math.pi)) * t for a, bi in zip((a1_inv, a2_inv, a3_inv), b))

def crossing_12(b):
    """解 1/α₁ = 1/α₂ 得 t，并给出该标度下三耦合"""
    b1, b2, b3 = b
    # a1 - b1 t/2π = a2 - b2 t/2π
    t = (a1_inv - a2_inv) / ((b1 - b2) / (2 * math.pi))
    return t, run(b, t)

SM_B = (41 / 10, -19 / 6, -7)
MSSM_B = (33 / 5, 1, -3)

ext = {}
for name, b in (("SM", SM_B), ("MSSM", MSSM_B)):
    t, (i1, i2, i3) = crossing_12(b)
    mu = MZ * math.exp(t)
    a1, a2, a3 = 1 / i1, 1 / i2, 1 / i3
    mism = a3 / a1
    ext[name] = {"b": b, "t": t, "scale_GeV": mu, "a1_inv": i1, "a2_inv": i2,
                 "a3_inv": i3, "alpha1": a1, "alpha2": a2, "alpha3": a3,
                 "alpha3_over_alpha1": mism, "mismatch_percent": 100 * (mism - 1)}
    print(f"\n  [{name}] α₁=α₂ 交点：t={t:.2f} ⟹ μ={mu:.3e} GeV")
    print(f"        该标度：α₁⁻¹={i1:.2f}、α₂⁻¹={i2:.2f}（等）、α₃⁻¹={i3:.2f}")
    print(f"        α₃/α₁ = {mism:.4f} ⟹ 三方不重合度 {100*(mism-1):+.1f}%")

chk(abs(ext["SM"]["scale_GeV"] - 1.0e13) / 1.0e13 < 0.5,
    "SM：α₁=α₂ 交点在 ~10^13 GeV 量级")
MATCH.append(("B1", f"SM：α₁=α₂ 交点 ≈ {ext['SM']['scale_GeV']:.2e} GeV", True))
chk(ext["SM"]["alpha3_over_alpha1"] > 1.10,
    "SM：α₃ 在该标度明显偏离（>10%）⟹ 三耦合不重合")
MATCH.append(("B2", f"**SM 不统一**：交点上 α₃ 高出 {ext['SM']['mismatch_percent']:+.1f}%", True))
chk(abs(ext["MSSM"]["alpha3_over_alpha1"] - 1) < 0.02,
    "MSSM：三耦合在 1% 内重合")
MATCH.append(("B3", f"MSSM（外部内容）：交点上 α₃ 偏差 {ext['MSSM']['mismatch_percent']:+.2f}%"
                    f"（≈{ext['MSSM']['scale_GeV']:.1e} GeV）", True))

# 交点上 sin²θ_w 的值（按定义就是 3/8 —— 说明「跑到低能得到 0.23」是自指的）
sin2_at_cross = 1 / (1 + (a2_inv * 1.0) / a1_inv) if False else None
sw_cross = (3 / 5) * (1 / ext["SM"]["a1_inv"]) / ((3 / 5) * (1 / ext["SM"]["a1_inv"]) + 1 / ext["SM"]["a2_inv"])
print(f"\n  · 交点上按定义算 sin²θ_w = {sw_cross:.6f}（α₁=α₂ ⟹ 恰为 3/8：**恒等式，不是预言**）")
print(f"  · 所以「3/8 经 RGE 跑到低能得 0.23，与观测一致」这句在 SM 内容下是**自指**的：")
print(f"    交点标度是用实测值反过来定的；真正的检验是 α₃ 是否同时重合 —— SM 差 "
      f"{ext['SM']['mismatch_percent']:+.1f}%，**不重合**。")

print()
print("断言清单：")
for tag, txt, good in MATCH:
    print(f"   {tag}  {'✓' if good else '✗'}  {txt}")
print(f"\n—— {ok} 条断言 / {len(bad)} 不符 ——")
for b_ in bad:
    print("   ✗", b_)

os.makedirs(RES, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump({
        "zero_side": {
            "slots": 16, "sum_T3sq": str(sT2), "sum_Yhalf_sq": str(sY2),
            "sum_T3_Yhalf": str(sTY), "sum_Qsq": str(sQ2),
            "sin2_thetaW_GUT": str(s_w), "sin2_thetaW_GUT_float": float(s_w),
            "no_nuR": str(sin2(sT2_15, sY2_15)), "left_only": str(sin2(sT2_L, sY2_L)),
            "wrong_normalization_k": alt,
            "inputs_beyond_zero": ["T₃ 与 Y/2 的迹范数相等（统一归一化）= Z086/Z135③ 的 U(1) 输入"],
        },
        "external_one_loop": ext,
        "verdict": ("Zero 侧：sin²θ_w|统一标度 = 3/8 精确（载体上 ΣT₃²/(ΣT₃²+ΣY'²)，"
                    "完备性由 ΣT₃Y'=0 保证；3/8 对 ν_R 不敏感）。"
                    "外部侧：SM 内容下三耦合不重合（α₃ 偏高 ~15%），"
                    "「3/8 跑到 0.23」在 SM 下自指；MSSM 内容才重合。"),
        "assertions": {"total": ok, "failed": len(bad), "detail": bad},
    }, f, ensure_ascii=False, indent=1)
print(f"已存 {OUT}")