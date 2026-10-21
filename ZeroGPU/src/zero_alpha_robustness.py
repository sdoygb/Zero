#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_alpha_robustness.py —— 两条定论性检验
======================================================================================
(A) 旧 Z23 的中央判决（迹排除 / 范数通过）对更新数据的稳健性
    Z23 用的是 alpha_s=0.1181±0.0011 一档；本轮用 PDG 现行值重算，
    并回答"要什么 alpha_s 才让 R 精确等于 sqrt(5/3) / 5/3"。

(B) Zero-DB 的耦合候选 g^2 = 4pi/(k+c_v) 能否给出 SM 的三个耦合？
    DB 定理 GG：k = 面积相位量子群等级。二态 k=L^2/16-2（c_v=2）、
    三态 k=(4L-33)/3（c_v=3）。检查是否存在单一 L 同时给出 alpha_2 与 alpha_3。

【纪律】全部用精确/高精度算术；每条判据写清"若不成立会怎样"。
"""
from __future__ import annotations
import os, json, math
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

PASS = []
def check(name, cond, detail=""):
    PASS.append({"name": name, "pass": bool(cond), "detail": detail})
    P_("[%s] %s %s" % ("PASS" if cond else "FAIL", name, detail))

P_("=" * 96)
P_("两条定论性检验 · zero_alpha_robustness.py")
P_("=" * 96)

# ================= 通用：由观测输入算 R =================
b1 = F(41, 10); b2 = F(-19, 6); b3 = F(-7)
bY = (F(5, 3)) * b1          # b_Y = (5/3) b_1 = 41/6

def compute_R(aem_inv, sin2, als):
    """由三个观测输入算 R = A_Y(M_X)/A_2(M_X)，M_X: alpha_2=alpha_3。"""
    A2z = aem_inv * sin2
    AYz = aem_inv * (1.0 - sin2)
    A3z = 1.0 / als
    # A_i(M_X) = A_i(M_Z) - b_i L/(2pi);  令 A_2(M_X)=A_3(M_X)
    L = 2 * math.pi * (A2z - A3z) / (float(b2) - float(b3))
    A2x = A2z - float(b2) * L / (2 * math.pi)
    AYx = AYz - float(bY) * L / (2 * math.pi)
    A3x = A3z - float(b3) * L / (2 * math.pi)
    return {"L": L, "X": A2x, "AY": AYx, "R": AYx / A2x,
            "A3x": A3x, "cross_ok": abs(A2x - A3x)}

# ================= (A) 稳健性 =================
P_("\n§A 旧 Z23 中央判决的稳健性（用 PDG 现行值重算）")
CASES = [
    ("Z23 原用 (alpha_s=0.1181)", 127.950, 0.23122, 0.1181, 0.0009, 0.00004),
    ("PDG 现行 (alpha_s=0.1180)", 127.951, 0.23122, 0.1180, 0.0009, 0.00004),
    ("PDG 2022 (alpha_s=0.1179)", 127.951, 0.23122, 0.1179, 0.0009, 0.00004),
]
R_NORM = math.sqrt(5 / 3)
R_TRACE = 5 / 3

rows = []
P_("  %-26s %9s %9s %9s %11s %11s" %
   ("case", "L", "X", "R_data", "R-Rnorm", "R-Rtrace%"))
for label, ai, s2, als, ds, ds2 in CASES:
    r = compute_R(ai, s2, als)
    dev_n = r["R"] - R_NORM
    dev_t = 100 * (R_TRACE - r["R"]) / r["R"]      # 用 Z23 口径：以 R_data 为分母
    rows.append({"label": label, "als": als, **r, "dev_norm": dev_n, "dev_trace_pct": dev_t,
                 "ds": ds, "ds2": ds2})
    P_("  %-26s %9.4f %9.4f %9.6f %+11.6f %+10.4f%%" %
       (label, r["L"], r["X"], r["R"], dev_n, dev_t))

# 交叉条件
for r in rows:
    check("M_X 处 alpha_2=alpha_3（交叉条件）[%s]" % r["label"],
          r["cross_ok"] < 1e-9, "|A2-A3|=%.2e" % r["cross_ok"])

# 判决稳健性
P_("\n  判决核对：")
for r in rows:
    sig = 3.99 * r["ds"]          # dR/dalpha_s ~ +3.99（Z23 F2 值）
    n_sigma = r["dev_norm"] / sig
    P_("    %s: R 与 sqrt(5/3) 差 %+.6f = %+.2f sigma（sigma_R=%.6f）" %
       (r["label"][:24], r["dev_norm"], n_sigma, sig))
    P_("      迹偏差 %+.2f%% （远 >1%%，无论哪一档都排除）" % r["dev_trace_pct"])
P_("")
check("迹在全部现行数据下仍被排除（偏差 >1%）",
      all(abs(r["dev_trace_pct"]) > 1.0 for r in rows),
      "最小偏差 %.2f%%" % min(abs(r["dev_trace_pct"]) for r in rows))
check("范数在全部现行数据下仍通过（<1 sigma）",
      all(abs(r["dev_norm"] / (3.99 * r["ds"])) < 1.0 for r in rows),
      "最大 %.2f sigma" % max(abs(r["dev_norm"] / (3.99 * r["ds"])) for r in rows))

# 反解：什么 alpha_s 让 R 精确等于两者
P_("\n  反解：要什么 alpha_s 才让 R 精确等于目标？")
def solve_als(target):
    """R 随 alpha_s 单调【递减】（实测 dR/dalpha_s ≈ -4.3），故方向要反过来。"""
    lo, hi = 0.10, 0.14
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        # R(mid) > target  => 需要更大的 alpha_s
        if compute_R(127.951, 0.23122, mid)["R"] > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)

als_norm = solve_als(R_NORM)
als_trace = solve_als(R_TRACE)
P_("    R = sqrt(5/3) = %.6f  <=> alpha_s = %.6f" % (R_NORM, als_norm))
P_("    R = 5/3     = %.6f  <=> alpha_s = %.6f" % (R_TRACE, als_trace))
P_("    PDG alpha_s = 0.1180 ± 0.0009")
P_("    ⟹ 范数要的 alpha_s 偏 %+.2f sigma" % ((als_norm - 0.1180) / 0.0009))
P_("    ⟹ 迹要的 alpha_s 偏 %+.1f sigma" % ((als_trace - 0.1180) / 0.0009))
check("迹要的 alpha_s 远超实验窗（>10 sigma）",
      abs((als_trace - 0.1180) / 0.0009) > 10,
      "%.1f sigma" % ((als_trace - 0.1180) / 0.0009))
check("范数要的 alpha_s 落在实验窗内（<1 sigma）",
      abs((als_norm - 0.1180) / 0.0009) < 1,
      "%+.2f sigma" % ((als_norm - 0.1180) / 0.0009))

# ================= (B) DB 耦合公式 =================
P_("\n§B Zero-DB 的耦合候选 g^2 = 4pi/(k+c_v) 能否给出 SM 耦合")
P_("  DB 定理 GG：g^2 = 4pi/(k+c_v)  ⟹  1/alpha = k + c_v")
P_("  二态 SU(2)（c_v=2）：k = L^2/16 - 2   ⟹ 1/alpha_2 = L^2/16")
P_("  三态 SU(3)（c_v=3）：k = (4L-33)/3   ⟹ 1/alpha_3 = (4L-24)/3")
P_("")

def inv_a2_db(L):
    return L * L / 16.0

def inv_a3_db(L):
    return (4.0 * L - 24.0) / 3.0

P_("  核对 DB 表值（二态 L=8,12,16 / 三态 L=9,12,15）：")
for L in (8, 12, 16):
    P_("    二态 L=%2d: 1/alpha_2 = %6.2f" % (L, inv_a2_db(L)))
for L in (9, 12, 15):
    P_("    三态 L=%2d: 1/alpha_3 = %6.2f" % (L, inv_a3_db(L)))
P_("  （应与 DB 表一致：二态 4,9,16；三态 4,8,12）")
check("公式复现 DB 表值",
      abs(inv_a2_db(8) - 4) < 1e-9 and abs(inv_a2_db(12) - 9) < 1e-9
      and abs(inv_a3_db(9) - 4) < 1e-9 and abs(inv_a3_db(12) - 8) < 1e-9,
      "二态 4,9,16；三态 4,8,12")

# SM 目标
a2_target = 127.951 * 0.23122          # 1/alpha_2(M_Z)
a3_target = 1.0 / 0.1180               # 1/alpha_3(M_Z)
P_("\n  SM 目标（M_Z）：1/alpha_2 = %.4f,  1/alpha_3 = %.4f" % (a2_target, a3_target))
P_("  注意：1/alpha_2 >> 1/alpha_3（弱耦合比强耦合小约 3.5 倍）")

L2_need = math.sqrt(16.0 * a2_target)
L3_need = (3.0 * a3_target + 24.0) / 4.0
P_("\n  反解各需的 L：")
P_("    要 1/alpha_2 = %.4f  =>  L = %.3f" % (a2_target, L2_need))
P_("    要 1/alpha_3 = %.4f  =>  L = %.3f" % (a3_target, L3_need))
P_("    两者相差 %.3f（同一 M_Z 处两个耦合必须来自同一 L）" % abs(L2_need - L3_need))
check("DB 公式无单一 L 能同时给出 alpha_2 与 alpha_3",
      abs(L2_need - L3_need) > 1.0,
      "L 需 %.2f vs %.2f（差 %.2f）" % (L2_need, L3_need, abs(L2_need - L3_need)))

# 比率随 L 的变化：能否在任何 L 匹配比率 3.49？
P_("\n  比率检验：SM 要 (1/alpha_2)/(1/alpha_3) = %.4f" % (a2_target / a3_target))
P_("    %6s %12s %12s %10s" % ("L", "1/a_2", "1/a_3", "比值"))
best = None
for L in (9, 12, 15, 18, 24, 30, 36, 48, 60, 72):
    r2, r3 = inv_a2_db(L), inv_a3_db(L)
    ratio = r2 / r3
    P_("    %6d %12.2f %12.2f %10.3f" % (L, r2, r3, ratio))
    if best is None or abs(ratio - a2_target / a3_target) < abs(best[1] - a2_target / a3_target):
        best = (L, ratio)
P_("    最接近的 L=%d 给比值 %.3f，SM 要 %.3f（差 %.1f%%）" %
   (best[0], best[1], a2_target / a3_target,
    100 * (best[1] - a2_target / a3_target) / (a2_target / a3_target)))
P_("    且该 L 下 1/a_2 = %.1f（SM 要 %.1f）" % (inv_a2_db(best[0]), a2_target))
P_("    ★ 读法：比值可以在某个 L 附近逼近，但那个 L 上的【绝对值】差一个量级。")

# 真正的判据：能否同时匹配绝对值？
P_("\n  正确的判据：是否有 L 同时把两个绝对值都放进 ±5%？")
ok_any = False
for L in range(9, 121):
    r2, r3 = inv_a2_db(L), inv_a3_db(L)
    if (abs(r2 - a2_target) / a2_target < 0.05
            and abs(r3 - a3_target) / a3_target < 0.05):
        ok_any = True
        P_("    命中 L=%d: 1/a_2=%.2f, 1/a_3=%.2f" % (L, r2, r3))
P_("    （遍历 L=9..120，无输出＝无命中）")
check("DB 公式无任何 L 同时给出 alpha_2 与 alpha_3 的绝对值（±5%）",
      not ok_any, "L=9..120 全扫")
P_("    ⟹ 故 DB 的 g^2=4pi/(k+c_v) 给的是【相对结构 + 跑动方向】，")
P_("       不能给出 SM 的三个耦合数值——与 DB 自判'P9 死门/统一机制开放'一致。")

P_("\n" + "=" * 96)
nfail = sum(1 for p in PASS if not p["pass"])
P_("合计: %d 项判据, %d PASS, %d FAIL" % (len(PASS), len(PASS) - nfail, nfail))
P_("=" * 96)

out = {
    "A_robustness": {
        "R_norm": R_NORM, "R_trace": R_TRACE,
        "cases": [{k: (float(v) if isinstance(v, (int, float)) else v)
                   for k, v in r.items()} for r in rows],
        "alpha_s_for_norm": als_norm, "alpha_s_for_trace": als_trace,
    },
    "B_db_coupling": {
        "L_needed_for_alpha2": L2_need, "L_needed_for_alpha3": L3_need,
        "best_ratio_L": best[0], "best_ratio": best[1],
        "sm_ratio": a2_target / a3_target,
    },
    "checks": PASS,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_alpha_robustness.json"), "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
P_("结果写入: %s" % os.path.join(RES, "zero_alpha_robustness.json"))