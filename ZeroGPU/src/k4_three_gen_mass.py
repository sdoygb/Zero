#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_three_gen_mass.py —— K₄ 三频率槽 → 三代质量比 → 对照观测(CKM/Koide)
======================================================================================
Z089 边界 #8: 3 频率槽的质量比是否对应观测(CKM/Koide), 待算.
K₄ 三槽: [0.405, 0.693, 0.981] (trip seed [8,4,1.5,1] → log 相邻比: log2, log(8/3), log1.5)
        且恰好等差 (间隔 0.288).

本脚本: 各种质量假设(m∝f, e^f, 2^f, seed 相邻比) → 归一化质量比 → 对照观测
        轻子 e:μ:τ=0.511:105.66:1776.86 MeV → Koide 值 vs 2/3.
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
    return float((m.sum()) / (np.sqrt(m).sum()) ** 2)

def norm_ratio(m):
    m = np.asarray(m, float)
    return m / m.min()

def main():
    P_("=" * 100)
    P_("K₄ 三频率槽 → 三代质量比 → 对照观测(CKM/Koide) (k4_three_gen_mass.py)")
    P_("=" * 100)
    out = {}

    # K₄ 三槽 (Z089): trip seed [8,4,1.5,1] → log 相邻比
    f_slots = np.array([0.405, 0.693, 0.981])          # 频率槽 (等差 0.288)
    seed = np.array([8.0, 4.0, 1.5, 1.0])
    ratios = np.array([seed[0]/seed[1], seed[1]/seed[2], seed[2]/seed[3]])  # [2, 2.667, 1.5]
    P_(f"K₄ 三频率槽 = {f_slots}  (log seed 相邻比, 等差 {f_slots[1]-f_slots[0]:.3f})")
    P_(f"seed 相邻比 = {ratios}  (log = {np.log(ratios).round(3)})")
    P_("")
    out["k4"] = {"f_slots": f_slots.round(3).tolist(), "seed_ratios": ratios.round(3).tolist()}

    # 观测: 轻子质量 (MeV)
    m_e, m_mu, m_tau = 0.511, 105.66, 1776.86
    obs = np.array([m_e, m_mu, m_tau])
    P_("【观测】轻子 e:μ:τ = %.3f : %.3f : %.2f MeV" % (m_e, m_mu, m_tau))
    P_(f"  质量比(最小=1) = {norm_ratio(obs).round(1)}")
    P_(f"  Koide 观测值 = (Σm)/(Σ√m)² = {koide(obs):.4f} (应≈2/3=0.6667)")
    P_("")
    out["observation"] = {"m_MeV": [m_e, m_mu, m_tau],
                          "ratio": norm_ratio(obs).round(1).tolist(),
                          "koide": round(koide(obs), 5)}

    # 各质量假设
    hyp = {
        "m ∝ f": f_slots,
        "m ∝ e^f": np.exp(f_slots),
        "m ∝ 2^f": 2.0 ** f_slots,
        "m ∝ seed 相邻比": ratios,
    }
    P_("【质量假设 → Koide】(K₄ 三槽线性映射, 均不满足 2/3)")
    rows = {}
    for name, m in hyp.items():
        m = np.asarray(m, float)
        k = koide(m)
        r = norm_ratio(m)
        P_(f"  {name:16s}: Koide={k:.4f} (vs 0.667); 质量比={r.round(3)}")
        rows[name] = {"koide": round(k, 5), "ratio": r.round(3).tolist()}
    out["hypotheses"] = rows
    P_("")

    # 反向: 若要 Koide=2/3, 需要质量比满足什么
    P_("【反向】Koide=2/3 的约束: 三代质量比须满足 (Σm)/(Σ√m)²=2/3")
    P_("  真实轻子满足 (0.6667); K₄ 线性映射全不满足 (≈0.34) ⟹ 质量比极端性需新机制")
    out["reverse"] = "真实轻子满足 Koide; K₄ 线性映射全不满足(≈0.34), 需新机制"

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_three_gen_mass.json"), "w", encoding="utf-8") as f_:
        json.dump(out, f_, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "k4_three_gen_mass.txt"), "w", encoding="utf-8") as f_:
        f_.write("\n".join(LOG))
    P_("\n已存 k4_three_gen_mass.json / .txt")

if __name__ == "__main__":
    main()
