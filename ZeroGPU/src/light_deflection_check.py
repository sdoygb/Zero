#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
light_deflection_check.py —— 补 #3：光线偏折的**因子 2**
=========================================================
对象：D36 §3.8.1「光线偏折推导出来了」（机制：质量附近势深 ⟹ 波速慢 ⟹ 波被偏折）。

三件事：
  L1 台账核对：Z104 的波动方程里 c²=0.8 是**常数**、L 不依赖 φ
      ⟹ 台账内**没有**「波速随势」这一项，该机制是新增输入
  L2 装上新机制（c_eff² = 1+2Φ）后，偏折角是多少
  L3 与 GR 对照：GR 全偏折 = 2× 折射型（h_00 与 h_xx 各贡献一半）

双路径（否定性结论的纪律）：
  路径 A：二阶射线方程 d/ds(n dr/ds)=∇n，RK4 积分（n = 1/c_eff）
  路径 B：一阶 eikonal  α = −∫ ∂_⊥ Φ dx（数值求积）
  两路径须一致；"GR = 2× 折射型"是弱场解析事实

诚实记录：**直接波场实验**（波包质心 / 平面波相位）在本实现下精度不足
  （本底 ±0.1 rad 量级 > 信号 0.05），故不作为判据，仅记录。

运行：cd ZeroGPU/src && /usr/bin/python3 light_deflection_check.py
"""
import json
import os
import sys

import numpy as np

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
CK, OUT = [], []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def note(name, value, detail=""):
    OUT.append({"name": name, "value": value, "detail": detail})
    print(f"      · {name} = {value}" + (f"   {detail}" if detail else ""))


SIG, AMP = 60.0, 5e-3


def Phi(x, y, amp=AMP):
    return -amp * np.exp(-(x * x + y * y) / (2 * SIG ** 2))


def n_of(x, y, amp=AMP):
    return 1.0 / np.sqrt(1.0 + 2 * Phi(x, y, amp))


def gradn(x, y, amp=AMP):
    f = -amp * np.exp(-(x * x + y * y) / (2 * SIG ** 2))
    gx, gy = f * (-x / SIG ** 2), f * (-y / SIG ** 2)
    return -gx / (1.0 + 2 * Phi(x, y, amp)) ** 1.5, -gy / (1.0 + 2 * Phi(x, y, amp)) ** 1.5


def ray(b, amp=AMP, x0=-400.0, x1=400.0, h=0.05):
    """路径 A：二阶射线方程（x 为参数）y'' = (1+y'²)/n·(∂_y n − y'∂_x n)，RK4"""
    def f(xx, yy, pp):
        nx_, ny = gradn(xx, yy, amp)      # gradn 返回 (∂_x n, ∂_y n) —— 顺序不能颠倒
        return (1 + pp * pp) / n_of(xx, yy, amp) * (ny - pp * nx_)
    y, yp, x = b, 0.0, x0
    while x < x1:
        k1 = f(x, y, yp); m1 = yp
        k2 = f(x + h / 2, y + h / 2 * m1, yp + h / 2 * k1); m2 = yp + h / 2 * k1
        k3 = f(x + h / 2, y + h / 2 * m2, yp + h / 2 * k2); m3 = yp + h / 2 * k2
        k4 = f(x + h, y + h * m3, yp + h * k3); m4 = yp + h * k3
        y += h / 6 * (m1 + 2 * m2 + 2 * m3 + m4)
        yp += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        x += h
    return float(np.arctan(yp))


def eik(b, amp=AMP, N=400001):
    """路径 B：eikonal α = −∫ ∂_⊥ Φ dx（沿 y=b 的直线）"""
    x = np.linspace(-800, 800, N)
    dp = -amp * (b / SIG ** 2) * np.exp(-(x ** 2 + b ** 2) / (2 * SIG ** 2))
    return float(-np.trapz(dp, x))


# ============================================================
print("=" * 72)
print("L1  台账核对：Z104 的波动方程里 c 是否依赖势")
print("=" * 72)
print("  Z104 §3：∂ₜ²φ = −c²Lφ + x(t)，c²=0.8（常数；k4_L1_wave.py 注释：稳定 c²·λ_max<4）")
print("  L 是固定图 Laplacian（不依赖 φ）")
ck("L1a 『波速随势变慢』这个机制在台账里有依据",
   False,
   "台账的 c² 是常数、L 不依赖 φ ⟹ 无 c(φ) 项；要产生偏折必须新增 c_eff(Φ)（新输入，"
   "按 Z016 §6 规则 4/6 须写明『上到 L2 的额外一步』）")

# ============================================================
print()
print("=" * 72)
print("L2/L3  装上新机制后：偏折角（双路径）与 GR 对照")
print("=" * 72)
print(f"  势阱 Φ = −{AMP}·exp(−r²/2σ²)，σ = {SIG}")
print(f"  {'b':>5} | {'路径A 射线 RK4':>16} | {'路径B eikonal':>14} | {'A/B':>7} | {'GR（=2×）':>12}")
rows = []
for b in (30.0, 45.0, 60.0, 90.0, 120.0):
    aA, aB = ray(b), eik(b)
    rows.append((b, aA, aB, abs(aA) / abs(aB), 2 * abs(aB)))
    print(f"  {b:>5.0f} | {aA:>16.9f} | {aB:>14.9f} | {abs(aA)/abs(aB):>7.4f} | {2*abs(aB):>12.9f}")
ratio = float(np.mean([r[3] for r in rows]))
ck("L2a 两条独立路径（射线 RK4 / eikonal）给出同一个偏折角",
   abs(ratio - 1) < 0.05,
   f"五个碰撞参数的偏折**大小**平均比 = {ratio:.4f}（逐点 {[round(r[3],4) for r in rows]}）")
print("  （两路径的符号差来自 yp 与 eikonal 的方向约定；物理方向都是朝势阱）")
conv = [ray(60.0, h=h) for h in (0.2, 0.1, 0.05, 0.025)]
ck("L2b 射线积分已收敛（步长无关）",
   max(conv) - min(conv) < 1e-6,
   f"b=60 时 h=0.2/0.1/0.05/0.025 给 {[round(c,9) for c in conv]}")
ck("L3a 该机制给出 GR 的偏折（α = 2×折射型）",
   abs(ratio - 2) < 0.1,
   f"实测 = 折射型 = ½·α_GR（比值 {ratio:.4f}，不是 2）⟹ 只给 GR 的一半")

# 太阳尺度
G, c, Msun, Rsun = 6.67430e-11, 2.99792458e8, 1.98847e30, 6.957e8
arc_refr = 2 * G * Msun / (c ** 2 * Rsun) * 180 / np.pi * 3600
print()
print(f"  太阳：折射型 2GM/(c²R) = {arc_refr:.3f}″；GR 全值 = {2*arc_refr:.3f}″")
note("太阳偏折", f"折射型 {arc_refr:.3f}″ vs GR {2*arc_refr:.3f}″",
     "D36 §3.8.1 的机制给前者；这正是历史上标量引力（Nordström）被否的判据")

print()
print("  诚实记录：直接波场实验（波包质心 / 平面波相位）在本实现下精度不足——")
print("    波包质心有无势阱无关的 −28 格点系统偏差（自由衍射＋吸收层）；")
print("    平面波相位的本底达 ±0.1 rad 量级，大于信号 0.05 ⟹ 不作为判据。")
note("波场实验", "精度不足，未用作判据", "本底 ±0.1 > 信号 0.05；改用两条几何光学路径")

# ============================================================
print()
print("=" * 72)
n_fail = sum(1 for x in CK if not x["pass"])
print(f"合计 {len(CK)} 条：PASS {len(CK)-n_fail} / FAIL {n_fail}")
print("=" * 72)
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "light_deflection_check.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": CK, "notes": OUT,
               "rows": [{"b": r[0], "alpha_ray": r[1], "alpha_eik": r[2], "ratio_AB": r[3],
                         "alpha_GR": r[4]} for r in rows],
               "ratio_AB_mean": ratio, "convergence": conv,
               "sun_arcsec_refraction": float(arc_refr), "sun_arcsec_GR": float(2 * arc_refr),
               "config": {"sigma": SIG, "amp": AMP}},
              f, ensure_ascii=False, indent=1)
print("结果写入 results/light_deflection_check.json")
sys.exit(1 if n_fail else 0)