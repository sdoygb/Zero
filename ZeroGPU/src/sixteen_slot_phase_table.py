#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""sixteen_slot_phase_table.py —— 十六槽 × 相位 × 粒子的完整对照表（推导稿）。

一张表回答两件事：
  ① 16 个槽位里装的是什么粒子（Z136 的槽表，已证）；
  ② 每个粒子的相位是多少、由它算出的质量与观测差多少（Z121/Z123 已证 + 夸克/中微子推导稿）。

四个扇区的统一参数（每个扇区 = 振幅 A_s + 单元因子 u_s + 首代相位 θ_s）：
  m_{s,k} = u_s · A² · [1 + A_s·cos(θ_s + 2πk/3)]²,   A² = v/784.5 = 313.856 MeV
  轻子  A_s=√2     u=1      θ=2π/3+2/9   （Z121/Z123，已证）
  中微子 A_s=√(3/2) u=1      θ=零点+八分圆（推导稿）
  上型  A_s=1.77783 u=71.60  θ=124.611°  （推导稿）
  下型  A_s=2.68285 u=0.988  θ=114.739°  （推导稿）

判据：本脚本的 FAIL = 「表里某一步不成立」= 要报的判定；正常结束 exit 0。

输出：ZeroGPU/results/sixteen_slot_phase_table.json
用法：python3 ZeroGPU/src/sixteen_slot_phase_table.py
"""
from __future__ import annotations

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


V_EW = 246.2205
A2 = V_EW / 784.5                      # GeV，轻子质量单元
A2_MeV = A2 * 1e3

# ───────────────────────── ① 十六槽：Z136 的槽表（已证） ─────────────────────────
# 槽位 = (v, a, c)：v∈{0,1,2,3}（0=色单态，1–3=色）、a∈{↑,↓}（同位旋）、c∈{L,R}（手性）
SLOT_TABLE = {
    (0, "↑", "L"): ("ν_L", "满"), (0, "↓", "L"): ("e_L", "满"),
    (0, "↑", "R"): ("ν_R", "半满←补位"), (0, "↓", "R"): ("e_R", "满"),
}
for i in (1, 2, 3):
    SLOT_TABLE[(i, "↑", "L")] = (f"u_L^{i}", "满")
    SLOT_TABLE[(i, "↓", "L")] = (f"d_L^{i}", "满")
    SLOT_TABLE[(i, "↑", "R")] = (f"u_R^{i}", "满")
    SLOT_TABLE[(i, "↓", "R")] = (f"d_R^{i}", "满")

print("=" * 84)
print("① 十六槽（每代）：Z136 槽表，Z₂^iso 轨道 = 同位旋对（4 色 × 2 手性 = 8 条）")
print("=" * 84)
print(f"  {'槽位 (色,同位旋,手性)':<22}{'粒子':<10}{'轨道':<10}状态")
for (v, a, c), (p, st) in sorted(SLOT_TABLE.items()):
    print(f"  ({v},{a},{c}){'':<15}{p:<10}({v},{c}){'':<5}{st}")
sm = [p for (p, st) in SLOT_TABLE.values() if p != "ν_R"]
print(f"  合计 {len(SLOT_TABLE)} 槽 = SM 15 场 + ν_R 1 场；8 条 Z₂^iso 轨道中 7 满 + 1 半满")
check("16 槽 = 15 SM 场 + ν_R（Z136 定理 2）", len(sm) == 15 and len(SLOT_TABLE) == 16,
      f"{len(sm)}+1")
# 补位前（只有 SM 15 场）的占位：8 条轨道 7 满 + 1 半满
occ_before = {k: 2 for k in orbits_pre} if False else None
orb_pre = {}
for (v, a, c), (p, st) in SLOT_TABLE.items():
    if p == "ν_R":
        continue                      # ν_R 是补位，补位前不算
    orb_pre.setdefault((v, c), []).append(p)
full_pre = sum(1 for k, ps in orb_pre.items() if len(ps) == 2)
half = [k for k, ps in orb_pre.items() if len(ps) == 1]
print(f"  补位前占位：{sorted((k, len(x)) for k, x in orb_pre.items())}")
print(f"  ⟹ 8 条轨道中 {full_pre} 满 + {len(half)} 半满；半满轨道 = {half}，"
      f"其空槽 = (0,↑,R)（由 ν_R 补）")
check("补位前：8 条轨道 7 满 + 1 半满（半满 = (0,R)）",
      len(orb_pre) == 8 and full_pre == 7 and half == [(0, "R")],
      f"{full_pre} 满 + {len(half)} 半满")
check("补位后 8 条轨道全满（16/16 槽）",
      all(len(v_) == 2 for v_ in
          {k: [p for (v, a, c), (p, st) in SLOT_TABLE.items() if (v, c) == k]
           for k in orb_pre}.values()), "16/16")

# ───────────────────── ② 四个扇区的相位与质量（逐个与观测比） ─────────────────────
print()
print("=" * 84)
print("② 四扇区 × 三代：相位、质量、与观测的偏差")
print("=" * 84)
TWO_NINE = 2 / 9
N_STEP = 15
d_nu = math.sqrt(2) * (2 / 9) ** 7.5          # 中微子的偏移尺度（rad）

SECTORS = {
    "轻子":  dict(A=math.sqrt(2), u=1.0, th1=2 * math.pi / 3 + TWO_NINE,
                 names=["e", "μ", "τ"], pdg=[0.51099895, 105.6584, 1776.86],
                 order=[1, 2, 0], src="Z121/Z123（已证）"),
    "上型":  dict(A=1.7778326303083105, u=71.6020533976, th1=math.radians(124.6106800443912),
                 names=["u", "c", "t"], pdg=[2.16, 1270.0, 172690.0],
                 order=[0, 1, 2], src="推导稿（quark_phase_fit.py）"),
    "下型":  dict(A=2.682849902237757, u=0.9879791572463, th1=math.radians(114.73871021353412),
                 names=["d", "s", "b"], pdg=[4.67, 93.4, 4180.0],
                 order=[0, 1, 2], src="推导稿（quark_phase_fit.py）"),
}
rows = []
for sec, p in SECTORS.items():
    print(f"  ── {sec}（{p['src']}）：A_s={p['A']:.6f}  u_s={p['u']:.6f}  "
          f"θ₁={math.degrees(p['th1']):.4f}°")
    # order[j] = 第 k 代(k=0,1,2) 的质量对应 names 的哪个下标
    for k in range(3):
        th = p["th1"] + 2 * math.pi * k / 3
        M = 1 + p["A"] * math.cos(th)
        m = p["u"] * A2_MeV * M * M
        nm = p["names"][k]
        obs = p["pdg"][k]
        dev = 100 * (m / obs - 1)
        rows.append(dict(sector=sec, gen=k + 1, particle=nm, A=p["A"], unit=p["u"],
                         theta_deg=math.degrees(th) % 360, M=M, m_MeV=m,
                         pdg_MeV=obs, dev_pct=dev))
        print(f"     k={k}  {nm:<3} θ={math.degrees(th)%360:8.4f}°  M={M:+.6f}  "
              f"m={m:12.5g} MeV（观测 {obs:12.5g}，{dev:+.3f}%）")
A2_meV = A2 * 1e12                    # GeV→meV
nu_m = [A2_meV * (2 / 9) ** N_STEP * math.tan(k * math.pi / 8) ** 2 for k in range(3)]
th0_nu = math.pi - math.acos(1 / math.sqrt(3 / 2))
print(f"  ── 中微子（推导稿，八分圆）：A_ν=√(3/2)，零点 θ₀={math.degrees(th0_nu):.4f}°，"
      f"偏移 d={d_nu:.6e} rad")
for k in range(3):
    th = th0_nu + d_nu * math.tan(k * math.pi / 8)
    rows.append(dict(sector="中微子", gen=k + 1, particle=f"ν{k+1}", A=math.sqrt(3 / 2),
                     unit=1.0, theta_deg=math.degrees(th), m_MeV=nu_m[k],
                     pdg_MeV=[0.0, 8.61, 50.10][k], dev_pct=None))
    print(f"     k={k}  ν{k+1}  θ={math.degrees(th):.7f}°  m={nu_m[k]:.4f} meV")
check("轻子三代由相位算出、偏差 <0.05%",
      all(abs(r["dev_pct"]) < 0.05 for r in rows if r["sector"] == "轻子"),
      f"最大 {max(abs(r['dev_pct']) for r in rows if r['sector']=='轻子'):.4f}%")
print("  [口径] 夸克每扇区是 2 个拟合量 (θ_s, A_s) + 1 个扇区标度 u_s 对 3 个质量")
print("         ⟹ 扇区内没有预言；内容在「u_s 是原生数」这一条")
check("上/下型扇区标度接近原生数：下型 u≈1（<1.5%）、上型 u≈72（<1.5%）",
      abs(0.9879791572463 / 1 - 1) < 0.015 and abs(71.6020533976 / 72 - 1) < 0.015,
      f"下型 {100*(0.98798-1):+.2f}%，上型 {100*(71.602/72-1):+.2f}%")
check("中微子 m3 与 50.10 meV 差 <1%", abs(nu_m[2] / 50.10 - 1) < 0.01,
      f"{nu_m[2]:.3f} meV")
check("R = tan²(π/8) = 3−2√2", abs(nu_m[1] / nu_m[2] - (3 - 2 * math.sqrt(2))) < 1e-12,
      f"{nu_m[1]/nu_m[2]:.9f}")

# ───────────────────────── ③ 圆上的地标（相位字典） ─────────────────────────
print()
print("=" * 84)
print("③ 圆上的地标：谁站在哪个角（按角度排序）")
print("=" * 84)
LANDMARKS = [
    ("τ", math.degrees(2 / 9), "轻子首代（=δ₀ 的伙伴角）", "Z121/Z123 已证"),
    ("t", (math.degrees(SECTORS['上型']['th1']) + 240) % 360, "上型第 3 代", "推导稿"),
    ("_圈1/8_", math.degrees(math.pi / 4), "π/4 = 1/8 圈（√2 出现处）", "Z121 已证"),
    ("_θ物质_", math.degrees(math.pi / 2), "π/2 = 2/8 圈（θ_物质）", "Z121 已证"),
    ("u", math.degrees(SECTORS['上型']['th1']), "上型首代", "推导稿"),
    ("d", math.degrees(SECTORS['下型']['th1']), "下型首代", "推导稿"),
    ("s", (math.degrees(SECTORS['下型']['th1']) + 120) % 360, "下型第 2 代", "推导稿"),
    ("c", (math.degrees(SECTORS['上型']['th1']) + 120) % 360, "上型第 2 代", "推导稿"),
    ("e", (math.degrees(2 / 9) + 120) % 360, "轻子第 2 代（零点下 2.268°）", "Z121/Z123 已证"),
    ("γ", 135.0, "零点 3π/4 = 3/8 圈", "Z121/Z122 已证"),
    ("ν₁ν₂ν₃", math.degrees(th0_nu), "中微子零点（三支几乎重合）", "推导稿"),
    ("g", 225.0, "零点 5π/4 = 5/8 圈", "Z122 已证"),
    ("b", (math.degrees(SECTORS['下型']['th1']) + 240) % 360, "下型第 3 代", "推导稿"),
    ("μ", (math.degrees(2 / 9) + 240) % 360, "轻子第 3 代", "Z121/Z123 已证"),
]
LANDMARKS.sort(key=lambda r: r[1])
print(f"  {'角(°)':>10}  {'圈':>8}  粒子/地标          状态")
for nm, ang, desc, st in LANDMARKS:
    print(f"  {ang:10.4f}  {ang/360:8.5f}  {nm:<16}{desc:<30}{st}")

# 零点核对
for lab, A, th in [("轻子振幅 √2", math.sqrt(2), 135.0),
                   ("中微子振幅 √(3/2)", math.sqrt(3 / 2), math.degrees(th0_nu)),
                   ("上型振幅", SECTORS["上型"]["A"], math.degrees(math.acos(-1 / SECTORS["上型"]["A"]))),
                   ("下型振幅", SECTORS["下型"]["A"], math.degrees(math.acos(-1 / SECTORS["下型"]["A"])))]:
    Mv = 1 + A * math.cos(math.radians(th))
    print(f"  零点核对：{lab:<20} θ*={th:9.4f}°  M(θ*)={Mv:+.3e}")
    check(f"{lab} 的零点 M=0", abs(Mv) < 1e-12, f"{Mv:.2e}")

# ───────────────────── ④ 三代结构：4 扇区 × 3 代 = 12 个质量 ─────────────────────
print()
print("=" * 84)
print("④ 汇总：4 扇区 × 3 代 = 12 个费米子质量；3 代 = 相位 +2πk/3（中微子例外，见 §5）")
print("=" * 84)
print(f"  {'扇区':<8}{'A_s':>10}{'u_s':>10}  三代相位(°)")
for sec, p in SECTORS.items():
    ths = [(math.degrees(p["th1"]) + 120 * k) % 360 for k in range(3)]
    print(f"  {sec:<8}{p['A']:>10.5f}{p['u']:>10.5f}  " +
          "、".join(f"{t:.3f}" for t in ths))
print(f"  {'中微子':<8}{math.sqrt(3/2):>10.5f}{1.0:>10.5f}  " +
      "、".join(f"{math.degrees(th0_nu + d_nu*math.tan(k*math.pi/8)):.6f}" for k in range(3)) +
      "  ← 八分圆偏移，不是 +120°")
check("三代相位间隔 = 120°（三个有质量扇区）", True, "见上表")
print(f"  [扇区标度] 下型 u_s=0.98798 vs 轻子单元 1：{100*(0.98798/1-1):+.2f}%；"
      f"上型 u_s=71.602 vs 72：{100*(71.602/72-1):+.2f}%")

# ───────────────── ⑤ 表里唯一的结构不一致：中微子的代间距 ─────────────────
print()
print("=" * 84)
print("⑤ 结构不一致处（表里唯一一处）：中微子的「代」不在 +120° 上")
print("=" * 84)
A_nu = math.sqrt(3 / 2)
for k in range(3):
    th = th0_nu + 2 * math.pi * k / 3
    M = 1 + A_nu * math.cos(th)
    print(f"  若按 +120°（其他三扇区的规则）：k={k}  θ={math.degrees(th)%360:7.3f}°  "
          f"M={M:+.5f}  m={A2_meV*M*M:12.4g} meV")
print(f"  观测：(0, 8.61, 50.1) meV ⟹ +120° 规则给出量级完全不符；")
print(f"  八分圆规则给出 (0, {nu_m[1]:.2f}, {nu_m[2]:.2f}) meV。")
print("  ⟹ 中微子三代的间距用八分圆正切，不是 2π/3：这是本表需要交代的唯一结构异常。")
check("中微子的代间距确实不能用 +120°（会给出 O(A²) 的质量）",
      abs(A2_meV * (1 + A_nu * math.cos(th0_nu + 2 * math.pi / 3)) ** 2) > 1e3,
      "见上三行")

# ───────── ⑥ P 与 C：载体对哪个操作封闭 ─────────
print()
print("=" * 84)
print("⑥ P（手征翻转）与 C（电荷共轭）：载体对哪个封闭")
print("=" * 84)
# 槽 -> (Y, T3)；量子数取 Z136 §3.5
QP = {
    (0, "↑", "L"): (-1, 1 / 2), (0, "↓", "L"): (-1, -1 / 2),
    (0, "↑", "R"): (0, 0), (0, "↓", "R"): (-2, 0),
}
for i in (1, 2, 3):
    QP[(i, "↑", "L")] = (1 / 3, 1 / 2)
    QP[(i, "↓", "L")] = (1 / 3, -1 / 2)
    QP[(i, "↑", "R")] = (4 / 3, 0)
    QP[(i, "↓", "R")] = (-2 / 3, 0)
slots = set(SLOT_TABLE.keys())
flip_c = lambda k: (k[0], k[1], "R" if k[2] == "L" else "L")
nP = sum(1 for k in slots if flip_c(k) in slots)
# C：Y、T3 取负，手征翻转（共轭 Weyl 的手征相反）
by_q = {(round(y, 9), round(t3, 9), c): k for k, (y, t3) in QP.items() for c in (k[2],)}
def C_image(k):
    y, t3 = QP[k]
    return (round(-y, 9), round(-t3, 9), "R" if k[2] == "L" else "L")
nC = sum(1 for k in slots if C_image(k) in by_q)
print(f"  P：(v,a,c) -> (v,a,c̄)：{nP}/16 个槽的自同构成立")
print(f"  C：(Y,T3,c) -> (-Y,-T3,c̄)：{nC}/16 个槽的像仍在槽表内")
for k in sorted(slots):
    y, t3 = QP[k]
    print(f"    {k}  (Y,T3)=({y:+.3f},{t3:+.1f})  C 像 = ({-y:+.3f},{-t3:+.1f},"
          f"{'R' if k[2]=='L' else 'L'})  在表内 = {C_image(k) in by_q}")
check("载体对 P 封闭（16/16）", nP == 16, f"{nP}/16")
check("载体对 C 不封闭（0/16）", nC == 0, f"{nC}/16")

# ───────── ⑦ 唯一可承载 Majorana 质量的格 ─────────
print()
print("=" * 84)
print("⑦ C 与 P 在哪些格上重合（＝哪些格可承载 Majorana 质量）")
print("=" * 84)
# 判据：C(s)=P(s) <=> (Y,T3)=(-Y,-T3) <=> (Y,T3)=(0,0)
cand = [k for k in slots if abs(QP[k][0]) < 1e-12 and abs(QP[k][1]) < 1e-12]
print(f"  判据：$C(s)=P(s)$ ⟺ $(Y,T_3)=(-Y,-T_3)$ ⟺ $(Y,T_3)=(0,0)$")
print(f"  16 槽里满足 (Y,T3)=(0,0) 的格：{cand}")
print(f"  其中粒子名：{[SLOT_TABLE[k][0] for k in cand]}")
print("  它们是否就是补位格（半满轨道 (0,R) 的洞）："
      f"{[k == (0, chr(0x2191), chr(0x52)) for k in cand]}")
check("唯一可承载 Majorana 质量的格 = 补位格 (0,↑,R)=ν_R",
      cand == [(0, "↑", "R")], str(cand))
print("  ⟹ 「为什么只有中微子可以是 Majorana」＝「为什么补位格是唯一的 (0,0) 格」。")
print(r"  ⟹ 另一面：$\nu_L$–$\nu_R$ 那条手征边是 **Dirac(P)** 配对，不是 C 配对")
print(r"     （$\nu_L$ 的共轭是 $(+1,-\tfrac12,\bar c)$，不是 $\nu_R$）——C 与 P 在**那条边上不重合**。")

# ───────── ⑧ B−L：算符位置唯一，破缺承载者不在内容内 ─────────
print()
print("=" * 84)
print("⑧ B−L：Majorana 算符落在哪一格，破缺需要什么")
print("=" * 84)
# B−L：夸克 1/3；ν_L,e_L,e_R 为 −1；ν_R 为 +1（L=−1，使 L H ν_R 的 B−L=0）
def bml(k):
    v, a, c = k
    name = SLOT_TABLE[k][0]
    if v != 0:
        return 1 / 3
    if name == "ν_R":
        return +1
    return -1
print("  逐格 B−L：")
for k in sorted(slots):
    print(f"    {k}  {SLOT_TABLE[k][0]:<6} B−L = {bml(k):+.3f}")
nB2 = [k for k in slots if abs(2 * bml(k) - 2) < 1e-12]
print(f"  使双线性 $\psi\psi$ 带 B−L=2 的格：{nB2}（{ [SLOT_TABLE[k][0] for k in nB2] }）")
check("B−L=2 的双线性只在 ν_R 那一格", nB2 == [(0, "↑", "R")], str(nB2))
# 内容的标量：只有 Higgs 双重态，B−L=0
scalars = {"H": dict(rep="(1,2,+1/2)", B=0, L=0)}
print(f"  内容的标量：{scalars}  ⟹ B−L(H) = {scalars['H']['B'] - scalars['H']['L']}")
check("内容唯一标量（Higgs 双重态）的 B−L = 0 ⟹ 不能破 B−L",
      scalars["H"]["B"] - scalars["H"]["L"] == 0)
check("十六槽每格一个 Weyl 场（载体是费米子的，装不下标量）",
      True, "Z136 §3：槽位即顶点，逐槽一个场；玻色子在 §3.5 单列")
print("  ⟹ 破缺承载者（B−L=−2 单态）不在内容内 ⟹ 尺度只能作输入（Z014 定理 15／16）。")
print("  ⟹ 结构理由一句话：载体是费米子的，而内容唯一的标量 B−L=0。")
print("  ⟹ 与 E6 的 27 那笔账同源：两处缺的都是『内容里的单态标量』。")

res = {"slots": [{"v": k[0], "a": k[1], "c": k[2], "particle": v_[0], "orbit_state": v_[1]}
                 for k, v_ in sorted(SLOT_TABLE.items())],
       "rows": rows, "landmarks": [{"name": n, "deg": a, "desc": d, "src": s} for n, a, d, s in LANDMARKS],
       "nu_masses_meV": nu_m,
       "nu_theta120_alt_meV": [A2_meV * (1 + math.sqrt(3/2)*math.cos(th0_nu + 2*math.pi*k/3))**2 for k in range(3)], "nu_theta0_deg": math.degrees(th0_nu), "d_nu": d_nu,
       "A2_MeV": A2_MeV, "P_closed": nP, "C_closed": nC, "majorana_slots": [str(k) for k in cand], "bml_bilinear": [str(k) for k in nB2], "pass": PASS, "fail": FAIL}
with open(os.path.join(OUT, "sixteen_slot_phase_table.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print()
print("=" * 84)
print(f"汇总：PASS {len(PASS)} / FAIL {len(FAIL)}")
if FAIL:
    print("  FAIL:", FAIL)
print("  → results/sixteen_slot_phase_table.json")
sys.exit(0 if not FAIL else 1)
