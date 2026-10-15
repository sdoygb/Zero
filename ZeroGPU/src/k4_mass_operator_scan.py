#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_mass_operator_scan.py —— 试 Zero 内生质量算子 (底座必须是 Zero)
======================================================================================
目标: 质量在 Zero 里怎么来。不套 cosmos 的 θ_物质=π/2, 而是试 Zero 已知结构
      (K₄ 三槽频率、模流谱、log2 等) → 各种质量算子 F → 对照观测
      轻子 1:207:3477 + Koide 2/3.

忠实计算: 每种算子独立算质量比 + Koide, 与观测并排, 诚实报告成败.
"""
from __future__ import annotations
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

def koide(m):
    m = np.asarray(m, float)
    return float(m.sum() / np.sqrt(m).sum() ** 2)

def norm_ratio(m):
    m = np.asarray(m, float)
    return m / m.min()

# ============ Zero 已知结构 ============
f_slots = np.array([0.405, 0.693, 0.981])          # K₄ 三槽 (log seed 相邻比, 等差 0.288)
log2 = 0.6931471805599453                           # 模流频率 (Z119)
Hmod = np.array([0.4055, 1.0986])                  # 模流谱 = diag(log3/2, log3) (Z120); 差=log2
seed = np.array([8.0, 4.0, 1.5, 1.0])
seed_ratios = np.array([2.0, 8/3, 1.5])

# ============ 观测 ============
obs = np.array([0.511, 105.66, 1776.86])           # e:μ:τ MeV
OBS_R = norm_ratio(obs)                             # 1:207:3477
OBS_K = koide(obs)                                  # 0.6667

def main():
    P_("=" * 100)
    P_("试 Zero 内生质量算子 (底座=Zero) · k4_mass_operator_scan.py")
    P_("=" * 100)
    P_(f"Zero 已知: K₄ 三槽 f={f_slots} (等差 {f_slots[1]-f_slots[0]:.3f}); log2={log2:.4f}; 模流谱={Hmod}")
    P_(f"观测: 轻子质量比={OBS_R.round(1)}; Koide={OBS_K:.4f}")
    P_("")
    out = {"zero_inputs": {"f_slots": f_slots.round(4).tolist(), "log2": log2,
                           "Hmod": Hmod.round(4).tolist(), "seed_ratios": seed_ratios.round(4).tolist()},
           "observation": {"ratio": OBS_R.round(1).tolist(), "koide": round(OBS_K, 5)}}

    results = []
    def test(name, m, note=""):
        m = np.asarray(m, float)
        k = koide(m); r = norm_ratio(m)
        # 距离观测的质量比误差 (对数空间)
        err = float(np.log(r / OBS_R).std())
        koide_ok = abs(k - OBS_K) < 0.005
        ratio_ok = err < 0.5
        results.append({"name": name, "koide": round(k, 5), "ratio": r.round(3).tolist(),
                        "logerr": round(err, 3), "koide_ok": koide_ok, "ratio_ok": ratio_ok, "note": note})
        P_(f"  {'✓' if (koide_ok and ratio_ok) else '·'} {name:34s} Koide={k:.4f} 质量比={r.round(2)} "
           f"(logerr={err:.2f}) {note}")

    # ===== 族 A: 频率的直接/幂/指数映射 (全是温和, 预期失败) =====
    P_("\n【族 A: 频率线性/幂/指数 —— 三槽等差, 给不出极端比】")
    for n in [1, 2, 3, 4]:
        test(f"m ∝ f^{n}", f_slots ** n)
    test("m ∝ 1/f²", 1 / f_slots ** 2)
    for c in [3, 5, 8, 12]:
        test(f"m ∝ e^({c}·f)", np.exp(c * f_slots))
    # 真实比 1:207:3477 非等比 → e^{c·f}(等差f)给等比, 必失败

    # ===== 族 B: 模流谱 → 质量 (三槽 vs 模流谱) =====
    P_("\n【族 B: 模流谱 / 谱差】")
    test("m ∝ e^(Hmod 外推)", np.exp(Hmod[0] + (Hmod[1]-Hmod[0]) * np.arange(3)))
    test("m ∝ 1/(0.288·i + a)²", 1 / (0.288 * np.arange(3) + 0.405) ** 2)  # 槽差幂

    # ===== 族 C: 嫁接 cosmos "相角零点放大" (明示嫁接, 锚点未 Zero) =====
    P_("\n【族 C: 嫁接 '相角零点放大' M(θ)=1+√2cosθ (cosmos 结构, 底座非 Zero, 仅参考)】")
    # 若 θ_i=δ+2πi/3 且 δ 扫到零点 3π/4 附近 → e 轻, μτ 重
    # 反求: 要质量比 1:207:3477, 需 M_i² 比. 若 θ_i=δ+2πi/3, 看 δ 能否给出
    # 用最小二乘反解 δ (若三等分结构成立)
    def Mth(delta, i):  # √m ∝ 1+√2 cos(delta+2πi/3)
        return 1 + np.sqrt(2) * np.cos(delta + 2 * np.pi * i / 3)
    from scipy.optimize import least_squares
    # 反解: 找 δ 使 √m_i 比例 = (1+√2cos(delta+2πi/3))
    def resid(delta):
        v = np.array([Mth(delta[0], i) for i in range(3)])
        v = v / v.min()
        o = np.sqrt(obs) / np.sqrt(obs).min()
        return (v - o) / o
    sol = least_squares(resid, [2.3], bounds=([-np.pi], [np.pi]))
    P_(f"  三等分结构反解: δ={sol.x[0]:.4f} rad = {np.degrees(sol.x[0])%360:.2f}°; "
       f"最优残差={float(np.abs(resid(sol.x)).max()):.4f}  (若残差小则结构可重建)")
    test("三等分 M(θ)² (δ 反解)", np.array([Mth(sol.x[0], i) for i in range(3)]) ** 2,
         note=f"δ={np.degrees(sol.x[0])%360:.1f}° 残差={float(np.abs(resid(sol.x)).max()):.3f}")

    # ===== 族 D: 纯 Zero 尝试 —— 让三槽等差进"cos 零点放大" =====
    # 若把槽当相角 θ_i = c·f_i (等差), 质量 = (1+√2 cos θ_i)², 扫描 c 找 Koide 附近
    P_("\n【族 D: 槽→相角 θ_i=c·f_i, 质量=(1+√2cosθ_i)², 扫 c (纯 Zero 槽, 但 √2 是 cosmos 的)】")
    best = None
    for c in np.linspace(0.1, 8.0, 400):
        th = c * f_slots
        m = (1 + np.sqrt(2) * np.cos(th)) ** 2
        k = koide(m)
        if abs(k - OBS_K) < 0.005:
            if best is None or abs(k-OBS_K) < abs(best[0]-OBS_K):
                best = (k, c, norm_ratio(m))
    if best:
        test(f"θ=c·f, (1+√2cosθ)² (c={best[1]:.2f})", (1+np.sqrt(2)*np.cos(best[1]*f_slots))**2,
             note=f"Koide={best[0]:.4f}")
    else:
        test("θ=c·f 扫描: 无 c 使 Koide≈2/3", np.zeros(3), note="扫400点无命中")

    out["scans"] = results
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_mass_operator_scan.json"), "w", encoding="utf-8") as f_:
        json.dump(out, f_, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "k4_mass_operator_scan.txt"), "w", encoding="utf-8") as f_:
        f_.write("\n".join(LOG))
    P_("\n已存 k4_mass_operator_scan.json / .txt")

    # ===== 总结 =====
    P_("\n" + "=" * 100)
    ok1 = [r["name"] for r in results if r["koide_ok"] and r["ratio_ok"]]
    ok2 = [r["name"] for r in results if r["koide_ok"]]
    P_("完全命中(质量比+Koide 都符): " + (str(ok1) if ok1 else "无"))
    P_("仅 Koide 接近: " + (str(ok2) if ok2 else "无"))
    P_("诚实结论: Zero 已知结构(频率/谱)直接给不出轻子质量; 质量算子缺'角结构'.")

if __name__ == "__main__":
    main()
