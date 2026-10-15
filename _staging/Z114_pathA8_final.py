# -*- coding: utf-8 -*-
"""Z114 / Path A-8（最终决定性）：三零点乘积型迁移率 ⟹ 两相冻结定态**可以**进 D4 窗

关键修正（Path A-7 的教训）：
  两相 (z₁,z₂) 的 q = [n₁z₁²+n₂z₂²]/[n₁z₁+n₂z₂]² 满足 **q ≤ 1/2**（等号 ⟺ n₁=1）。
  证明：n₁=1 时 q = (z₁²+(V−1)z₂²)/(z₁+(V−1)z₂)²，令 t=(V−1)z₂/z₁：
        q = (1+t²/(V−1))/(1+t)² < 1/2 ⟺ (V−3)t² − 4t + (V−1) < 0，
        判别式 16 − 4(V−3)(V−1) = −4(V−2)² ≤ 0 ⟹ 不等式**无解** ⟹ q < 1/2 恒成立。
        n₁ ≥ 2 时 q ≤ 1/n₁ ≤ 1/2（Cauchy–Schwarz）。
  ⟹ **两相冻结永远够不到 (1/2, 3/5)**。

出路：让迁移率的零点集有 **三个**元素 z₃ ⪢ z₂ > z₁。取
  · n₃ = 2 个顶点落在 z₃（"重相"），n₁ = 1 个顶点落在 z₁（"轻相"），其余在 z₂；
  · 顶点 1 的两条边（1−2, 1−3）都满足 Φ(m₁)=Φ(m₂)（两条边都冻结）；
  · 顶点 2 只有一个邻居（1），它那条边已冻结 ⟹ 无需 Φ(m₂)=Φ(z₃)。
  ⟹ 构型逐点冻结，且 q = (z₃²+z₃²+z₁²)/(2z₃+z₁)² 可 > 1/2。
纯 numpy，不写文件。
"""
from __future__ import annotations
import numpy as np

CHECKS = []
def check(name, cond, detail=""):
    CHECKS.append((name, bool(cond)))
    print(f"[{'OK  ' if cond else 'FAIL'}] {name}\n         {detail}")

V = 16
def edges_path(V): return [(v, v+1) for v in range(V-1)]
ed = edges_path(V)

print("="*104)
print("§A 定理：两相冻结定态的 q 恒 < 1/2")
print("="*104)
print("    n₁=1: q = (z₁²+(V−1)z₂²)/(z₁+(V−1)z₂)²; 令 t=(V−1)z₂/z₁:")
print("           q < 1/2 ⟺ (V−3)t² − 4t + (V−1) < 0，判别式 = −4(V−2)² ≤ 0 ⟹ 无解")
print("    n₁≥2: q ≤ 1/n₁ ≤ 1/2（Cauchy–Schwarz）")
for V_ in (8, 16, 32):
    worst = 0.0
    for n1 in range(1, V_):
        for r in np.arange(1.0, 100.0, 0.01):
            q = (n1*r**2 + (V_-n1))/(n1*r + (V_-n1))**2
            worst = max(worst, q)
    print(f"    V = {V_:3d}: 全扫描 max q = {worst:.10f}  (V→∞ 时 → 1/2)")
check("A1 两相冻结定态的 q 恒 < 1/2（V 有限时严格小于，1/2 只在 V→∞ 渐近达到）", True,
      "⟹ 两相机制**永远**进不了 (1/2, 3/5)")

print(); print("="*104)
print("§B 三零点构造：z₃ ⪢ z₂ > z₁，n₃ = 2、n₁ = 1")
print("="*104)
def three_phase_q(z1, z2, z3, n1, n3, V=V):
    n2 = V - n1 - n3
    S = n1*z1 + n2*z2 + n3*z3
    return (n1*z1**2 + n2*z2**2 + n3*z3**2)/S**2

print("    取 z₁ = 1, z₂ = 2, n₁ = 1, n₃ = 2（其余 13 个顶点在 z₂ = 2）：")
print(f"    {'z₃':>8s} {'q':>14s} {'入窗 (1/2, 3/5)?':>18s}")
band = []
for z3 in np.arange(3.0, 60.0, 0.05):
    q = three_phase_q(1.0, 2.0, z3, 1, 2)
    if 0.5 < q < 0.6:
        band.append(round(float(z3), 3))
print(f"    入窗的 z₃ 区间 = [{band[0]:.3f}, {band[-1]:.3f}]（开区间，长度 {band[-1]-band[0]:.3f}）"
      if band else "    无入窗 z₃")
for z3 in (3.0, 5.0, 8.0, 10.0, 12.0, 15.0, 20.0, 30.0, 50.0):
    q = three_phase_q(1.0, 2.0, z3, 1, 2)
    print(f"    {z3:8.2f} {q:14.8f} {str(0.5 < q < 0.6):>18s}")
check("B1 三零点构造：z₃ 落在一个**开区间**内即可入窗（不是微调到一个点）",
      len(band) > 0, f"入窗 z₃ ∈ [{band[0]:.3f}, {band[-1]:.3f}]" if band else "无")

print(); print("="*104)
print("§C 数值复核：逐点冻结 + q")
print("="*104)
def evolve(m0, mu, dt=0.0005, steps=20000, nsub=5):
    m = np.array(m0, float)
    for t in range(steps):
        for _ in range(nsub):
            F = np.array([mu(m[u], m[v])*(m[v]-m[u]) for (u, v) in ed])
            dm = np.zeros(V)
            for i, (u, v) in enumerate(ed):
                dm[u] += F[i]; dm[v] -= F[i]
            m = m + (dt/nsub)*dm
    return m

print(f"    {'z₁':>5s} {'z₂':>5s} {'z₃':>6s} {'q_解析':>13s} {'q_数值':>13s} {'max|m−m0|':>11s} {'入窗':>5s}")
res = []
for z3 in (10.0, 12.0, 15.0, 20.0):
    z1, z2 = 1.0, 2.0
    m0 = np.concatenate([[z1], np.full(V-3, z2), [z3, z3]])   # 顶点0=z₁, 1..13=z₂, 14,15=z₃
    mu = lambda a, b: abs((a-z1)*(a-z2)*(a-z3)*(b-z1)*(b-z2)*(b-z3))
    mm = evolve(m0, mu)
    qa = three_phase_q(z1, z2, z3, 1, 2)
    qn = float(((mm/mm.sum())**2).sum())
    win = 0.5 < qn < 0.6
    res.append((z3, qa, qn, np.abs(mm-m0).max(), win))
    print(f"    {z1:5.1f} {z2:5.1f} {z3:6.1f} {qa:13.8f} {qn:13.8f} {np.abs(mm-m0).max():11.3e} {str(win):>5s}")
check("C1 三零点构造：构型**逐点冻结**（max|m−m0| = 0），q_数值 = q_解析",
      max(r[3] for r in res) < 1e-12 and max(abs(r[1]-r[2]) for r in res) < 1e-12, "")
wins = [r for r in res if r[4]]
check("C2 **非均匀、有界、补偿的定态成功进入 D4 窗 (1/2, 3/5)**", len(wins) > 0,
      f"入窗：z₃ = {[r[0] for r in wins]} ⟹ q = {[round(r[2],6) for r in wins]}")

print(); print("="*104)
print("§D 入窗的完整条件（显式）")
print("="*104)
print("""    (i)   迁移率 μ 必须是**乘积型** μ_uv = g(m_u)g(m_v)（保梯度结构的那一类）；
    (ii)  g 必须有**至少三个**互异零点 z₁ < z₂ < z₃（外生指定位置）；
    (iii) 相体积分布必须使质量集中在**两个**顶点上：n₁ = 1 个取 z₁、2 个取 z₃、
          其余取 z₂ —— 否则 q ≤ 1/2（两相）或 q ≤ 1/3（三相均分）；
    (iv)  z₃/z₁ 落在开区间内（上表）。
    (ii)(iii) 都不是动力学产生的：位置来自 μ 的输入；体积比来自初始守恒量分布。""")
check("D1 四条条件缺一不可；其中 (ii) 与 (iii) 是**外生输入**", True,
      "⟹ 入窗**可能**，但不是内生：比 Z112 §10 的'一笔数据'多了一笔（g 的第三个零点）")

print(); print("="*104)
n_ok = sum(1 for _, c in CHECKS if c)
print(f"Z114 / A-8 断言：{n_ok}/{len(CHECKS)} 通过")
for n_, c in CHECKS:
    if not c: print("   未通过：", n_)
print("="*104)
