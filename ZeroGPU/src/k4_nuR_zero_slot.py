#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_nuR_zero_slot.py —— Zero 三代频率槽直接约束中微子扇区
====================================================================================
不套跷跷板公式, 用 Zero 内生三代槽(质量归一[0.41339,0.70670,1.0], 来自K₄谱)
直接给中微子三代质量结构, 检验是否适用(普适三代结构?):
  设 m_i = A·norm_i (Zero槽), 两种定标:
    定标1(大气): A²(norm3²−norm2²)=Δm²_atm → 预测 Δm²_solar
    定标2(太阳): A²(norm2²−norm1²)=Δm²_solar → 预测 Δm²_atm
  比值结构: Zero槽 Δm²比 = (norm2²−norm1²)/(norm3²−norm2²) vs 观测 Δm²_solar/Δm²_atm
  若 Zero 槽是三代普适结构 → 中微子质量比应≈Zero槽比 → 偏差≈0
  若偏差大 → Zero 槽不适用中微子(近简并), 需额外机制/不同槽
观测: Δm²_solar=7.39e-5, Δm²_atm=2.45e-3
"""
from __future__ import annotations
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))

dms2_solar = 7.39e-5
dms2_atm   = 2.45e-3
norm = np.array([0.41339, 0.70670, 1.0])

def main():
    out = {"config": {"zero_slot_norm": list(norm), "dm2_solar": dms2_solar, "dm2_atm": dms2_atm}}
    print("="*90)
    print("Zero 三代频率槽直接约束中微子扇区 (不套跷跷板)")
    print("="*90)

    # === 定标1: 大气定标 → 预测太阳 Δm² ===
    A1 = np.sqrt(dms2_atm / (norm[2]**2 - norm[1]**2))
    m1_1 = A1*norm
    pred_solar1 = m1_1[1]**2 - m1_1[0]**2
    print("\n[定标1] 用大气 Δm²_atm 定标 Zero 槽")
    print(f"  m1={m1_1[0]:.4e}, m2={m1_1[1]:.4e}, m3={m1_1[2]:.4e} eV")
    print(f"  预测 Δm²_solar = {pred_solar1:.3e} eV²  vs 观测 {dms2_solar:.3e}")
    print(f"  偏差 = {pred_solar1/dms2_solar:.1f}×")

    # === 定标2: 太阳定标 → 预测大气 Δm² ===
    A2 = np.sqrt(dms2_solar / (norm[1]**2 - norm[0]**2))
    m1_2 = A2*norm
    pred_atm2 = m1_2[2]**2 - m1_2[1]**2
    print("\n[定标2] 用太阳 Δm²_solar 定标 Zero 槽")
    print(f"  m1={m1_2[0]:.4e}, m2={m1_2[1]:.4e}, m3={m1_2[2]:.4e} eV")
    print(f"  预测 Δm²_atm = {pred_atm2:.3e} eV²  vs 观测 {dms2_atm:.3e}")
    print(f"  偏差 = {dms2_atm/pred_atm2:.1f}×")

    # === 比值结构 ===
    r_zero = (norm[1]**2 - norm[0]**2)/(norm[2]**2 - norm[1]**2)
    r_obs = dms2_solar/dms2_atm
    print("\n[比值结构] Δm²_solar/Δm²_atm")
    print(f"  Zero槽预测比值 = {r_zero:.4f}  观测比值 = {r_obs:.4f}")
    print(f"  偏差 = {r_zero/r_obs:.1f}×")

    # === 观测中微子质量比(近简并) vs Zero 槽 ===
    # 正常序观测: m1≈0, m2=√(solar), m3=√(atm+solar)
    m2o = np.sqrt(dms2_solar); m3o = np.sqrt(dms2_atm+m2o**2)
    obs_ratio = [0.0, m2o/m3o, 1.0]
    zero_ratio = norm/norm[2]
    print("\n[质量比] m1:m2:m3 (归一化到 m3)")
    print(f"  观测(正常序近简并): {[round(x,4) for x in obs_ratio]}")
    print(f"  Zero槽:             {[round(x,4) for x in zero_ratio]}")
    print(f"  差异: Zero槽展宽(1:1.71:2.42) vs 观测近简并(0:0.17:1)")

    out["scaling1_atm"] = {"m": [round(m,6) for m in m1_1],
                           "pred_solar": round(float(pred_solar1),6), "obs_solar": dms2_solar,
                           "bias_x": round(float(pred_solar1/dms2_solar),1)}
    out["scaling2_solar"] = {"m": [round(m,6) for m in m1_2],
                             "pred_atm": round(float(pred_atm2),6), "obs_atm": dms2_atm,
                             "bias_x": round(float(dms2_atm/pred_atm2),1)}
    out["ratio"] = {"zero_slot": round(float(r_zero),4), "obs": round(float(r_obs),4),
                    "bias_x": round(float(r_zero/r_obs),1)}
    out["mass_ratio"] = {"obs_normal_order": [round(x,4) for x in obs_ratio],
                         "zero_slot": [round(x,4) for x in zero_ratio]}

    # === 结论判断 ===
    bias = r_zero/r_obs
    if bias > 5:
        verdict = "Zero槽(带电轻子展宽1:1.7:2.4)不直接适用于中微子(观测近简并); 需额外机制或中微子用不同槽"
    else:
        verdict = "Zero槽与观测中微子质量结构相容"
    out["verdict"] = verdict
    print(f"\n[结论] {verdict}")

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_nuR_zero_slot.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已存 k4_nuR_zero_slot.json")

if __name__ == "__main__":
    main()
