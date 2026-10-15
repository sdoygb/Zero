#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_nuR_features.py —— ν_R 详细特征: 三代质量谱 + 跷跷板 M_R + sterile 角色
====================================================================================
基于 Z098(ν_R sterile Majorana 候选) + Zero 三代质量槽(质量归一[0.413,0.707,1.0])。
计算详细特征:
  1) 三代轻中微子质量谱: 两种定标
     A 观测定标(近简并): Δm²_solar=7.39e-5, Δm²_atm=2.45e-3 (分层序)
     B Zero三代槽: m_i∝[0.413,0.707,1.0], 预测 Δm² 对比观测(诚实偏差)
  2) 跷跷板反推右手中微子质量:  M_R = y_ν²v²/m_ν, v=246GeV
  3) sterile 角色: 无规范耦合; 两个质量尺度(keV暗物质 vs 10¹⁵GeV跷跷板)
观测值: Δm²_solar≈7.39e-5 eV²(太阳), Δm²_atm≈2.45e-3 eV²(大气), m_ν<0.8eV(KATRIN/宇宙学)
"""
from __future__ import annotations
import sys, os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))

# 物理常数与观测
V_EW = 246.0            # GeV 电弱标度
dms2_solar = 7.39e-5    # eV² 太阳中微子质量差平方
dms2_atm   = 2.45e-3    # eV² 大气中微子
KEV = 1e-6              # 1 keV = 1e-6 GeV

def main():
    out = {"config": {"v_EW": V_EW, "dm2_solar": dms2_solar, "dm2_atm": dms2_atm,
                      "zero_mass_norm": [0.41339, 0.70670, 1.0]}}
    print("="*90)
    print("ν_R 详细特征: 三代质量谱 + 跷跷板 M_R + sterile 角色")
    print("="*90)

    # === 1A: 观测定标三代质量(近简并, 正常序) ===
    m1a = 0.0
    m2a = np.sqrt(dms2_solar)
    m3a = np.sqrt(dms2_atm + m2a**2)
    masses_obs = [m1a, m2a, m3a]
    print("\n[1A] 观测定标三代轻中微子质量(正常序, 近简并)")
    print(f"  m1={m1a:.3e}, m2={m2a:.3e}, m3={m3a:.3e} eV")

    # === 1B: Zero 三代槽预测中微子质量 ===
    norm = np.array([0.41339, 0.70670, 1.0])
    # 定标 A: 使 m3²−m2²=dms2_atm → A²(norm3²−norm2²)=dms2_atm
    A2 = dms2_atm / (norm[2]**2 - norm[1]**2)
    A = np.sqrt(A2)
    m_zero = A * norm
    dms2_solar_pred = m_zero[1]**2 - m_zero[0]**2
    print("\n[1B] Zero 三代槽预测中微子质量 (定标 Δm²_atm)")
    print(f"  m1={m_zero[0]:.3e}, m2={m_zero[1]:.3e}, m3={m_zero[2]:.3e} eV")
    print(f"  预测 Δm²_solar={dms2_solar_pred:.2e} eV²  vs 观测 {dms2_solar:.2e}")
    print(f"  偏差 = {dms2_solar_pred/dms2_solar:.1f}×")

    # === 2: 跷跷板反推右手中微子质量 M_R ===
    print("\n[2] 跷跷板反推 M_R = y_ν²v²/m_ν (Type-I)")
    seesaw = {}
    for y in [1.0, 0.1]:
        row = []
        for m in masses_obs[1:]:           # 用 m2, m3 (m1=0 退化)
            # M_R = y²v²/m_ν: v²[GeV²]/m_ν[eV] 单位 GeV²/eV = GeV²/(1e-9 GeV)=1e9 GeV
            M_R = y**2 * V_EW**2 * 1e9 / m
            row.append(M_R)
        seesaw[f"y={y}"] = {"M_R_ν2": row[0], "M_R_ν3": row[1]}
        print(f"  y_ν={y}: M_R(ν2)={row[0]:.2e} GeV, M_R(ν3)={row[1]:.2e} GeV")
    print("  (单位: M_R=y²v²[GeV²]×1e9/m_ν[eV] = GeV; v=246GeV)")

    # === 3: sterile 角色 + 暗物质候选 ===
    print("\n[3] sterile 角色: 两个质量尺度两种角色")
    print("  规范耦合: Y=T₃=色单态 → U(1)/SU(2)/SU(3) 全零 = sterile")
    # 暗物质候选: 轻 sterile ν_s ~keV (Dodelson-Widrow 机制)
    nu_s_keV = 7.1*KEV      # 7.1 keV 暗物质候选(观测: 3.5keV X射线线提示 7.1keV)
    seesaw_y1_nu3 = seesaw["y=1.0"]["M_R_ν3"]
    out["masses"] = {"obs_normal_order": [round(m,6) for m in masses_obs],
                     "zero_slot": [round(m,6) for m in m_zero],
                     "pred_vs_obs_dms2_solar": {"pred": round(float(dms2_solar_pred),6),
                                                "obs": dms2_solar, "ratio": round(float(dms2_solar_pred/dms2_solar),1)}}
    out["seesaw_MR"] = {k: {"M_R_ν2_GeV": round(v["M_R_ν2"],3), "M_R_ν3_GeV": round(v["M_R_ν3"],3)} for k,v in seesaw.items()}
    out["sterile_roles"] = {"keV_sterile_darkmatter": nu_s_keV,
                            "heavy_seesaw": seesaw_y1_nu3,
                            "note": "keV尺度=暗物质候选(DW); 10^14-10^16 GeV尺度=跷跷板重伙伴"}
    print(f"  keV sterile 暗物质候选: m~{nu_s_keV:.2e} GeV ({nu_s_keV/KEV:.1f} keV, DW机制)")
    print(f"  跷跷板重伙伴: M_R~{seesaw_y1_nu3:.2e} GeV (对 m_ν3~{m3a:.3f} eV)")

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "k4_nuR_features.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已存 k4_nuR_features.json")

if __name__ == "__main__":
    main()
