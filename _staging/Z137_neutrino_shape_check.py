#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Z137_neutrino_shape_check.py —— 中微子质量形状能否由现行质量机制（Z121/Z123）给出？

用户 2026-10-16 提问：Z019 是早期篇，此后理论有长足发展；中微子质量能否算准？
本文回答两件事：① 机制能否给出中微子的形状；② 若能，代价是什么。

【定义（先消歧——作者在本题上连续踩过同类坑，故此处把每个符号的"幂次"钉死）】
  模型量：M(θ) = 1 + √2·cos θ      —— 偏离振幅，**可负**；√m ∝ |M|
  形状量：s(θ) := M(θ)²            —— ∝ 质量 m（**≥0**，排序用它）
  质量：  m_i = A·s_i              —— A 为质量尺度（eV）
  观测量：X := Δm²21/Δm²31 = (s2²−s1²)/(s3²−s1²)   ← **含平方**（s1,s2,s3 升序）
          R := m2/m3 = s2/s3                        ← m1=0 时才有 X = R²
          K := Σm/(Σ√m)² = Σs/(Σ√s)²               ← Koide 比，尺度无关

【机制（Z121 立、Z123 补）】三代角间隔 2π/3（"三项"）；零点 3π/4、振幅 √2 来自闭合 2π 分 8 份。

纪律：纯标准库；不 import 项目脚本；不写文件。
"""
from __future__ import annotations
import math

FAIL = 0
NCHK = 0

def ck(c, m):
    global FAIL, NCHK
    NCHK += 1
    if not c:
        FAIL += 1
        print("  FAIL:", m)

S2 = math.sqrt(2)
PI = math.pi
def M(t): return 1 + S2*math.cos(t)
def s(t): return M(t)**2

D21, D31 = 7.42e-5, 2.510e-3          # eV²，正序全球拟合（与 Z019 同源）
X_OBS = D21/D31
R019 = (73 - 28*math.sqrt(6))/25
ZP = (3*PI/4, 5*PI/4)                 # 两个零点（M=0）
COSMO = 0.064                          # eV，Z019 引用的宇宙学上限

def shape(t1): return sorted(s(t1 + 2*PI*i/3) for i in range(3))
def X_of(t1):
    v = shape(t1)
    return (v[1]**2 - v[0]**2)/(v[2]**2 - v[0]**2) if v[2] > v[0] else float('inf')
def K_of(t1):
    v = shape(t1)
    return sum(v)/(sum(math.sqrt(x) for x in v))**2
def delta0(t1):
    return math.degrees(min(min(abs(t1-z), abs(t1-z+2*PI), abs(t1-z-2*PI)) for z in ZP))

# ── Z123：偏移 δ₀ 完全内生（这是本题的**权威来源**）──
#   δ₀ = π/12 − (2/3)(1/3) = π/12 − 2/9
#      π/12 = 3π/4 − 2π/3（零点 − 一代轮转）——seed 8 格点纯几何
#      2/9  = (2/3)(1/3)（稳态两支权重积）——非中心态固有量
#   两个分量**都与扇区无关** ⟹ 机制在带电轻子侧**没有拟合参数**。
D0_E = PI/12 - 2/9
TH_E = 2*PI/3 + 2/9          # Z123 §1.3：电子相角 δ_e = 2π/3 + 2/9
def lepton_angles(): return [TH_E + 2*PI*i/3 for i in range(3)]
def masses_of(t1):
    """返回 (m 升序 [eV], A)。A² = Δm²31/(s3²−s1²)。"""
    v = shape(t1)
    A = math.sqrt(D31/(v[2]**2 - v[0]**2))
    return [A*x for x in v], A

print("=== 0. 定义与量纲先对清楚 ===")
print("  X := Δm²21/Δm²31 = (s2²−s1²)/(s3²−s1²)  [含平方]")
print("  R := m2/m3 = s2/s3                        [m1=0 时 X = R²]")
print("  K := Σm/(Σ√m)² = Σs/(Σ√s)²")
ck(abs(R019**2 - 0.0311774904) < 1e-9, "0a: Z019 的 R² = %.9f" % R019**2)
ck(abs(X_OBS - 0.029562) < 5e-7, "0b: X_obs = %.6f" % X_OBS)
ck(abs(R019**2/X_OBS - 1 - 0.05466) < 1e-3,
   "0c: Z019 形状比实测高 %.2f%%（＝Z019 §3.3 的 [B] 张力 +5.47%%）" % ((R019**2/X_OBS-1)*100))

# ── A. 机制在带电轻子侧的自检 ──
print("\n[A] 机制在带电轻子侧（若这里不成立，后面无需谈）")
d0e = D0_E                       # Z123：完全内生，非拟合
te = TH_E                        # = 3π/4 − δ₀ = 2π/3 + 2/9
amps = [M(t) for t in lepton_angles()]   # 振幅（非角！）
r = [abs(x)/abs(amps[0]) for x in amps]          # √m ∝ |M|
obs = [1.0, 14.38, 58.968]
print("     √m 比 = [1, %.4f, %.4f]  vs 观测 [1, %.2f, %.3f]" % (r[1], r[2], obs[1], obs[2]))
ck(abs(r[1]/obs[1]-1) < 5e-5, "A1: μ 偏差 %.4f%%" % (abs(r[1]/obs[1]-1)*100))
ck(abs(r[2]/obs[2]-1) < 5e-5, "A2: τ 偏差 %.4f%%" % (abs(r[2]/obs[2]-1)*100))
ck(abs(d0e - 0.039577166) < 1e-8, "A3: δ₀^e = %.9f（带电轻子的输入）" % d0e)

# ── B. 机制的恒等签名：Koide ≡ 2/3（三个 M 同号时）──
# Σcos(θ+2πk/3)=0（1+ω+ω²=0）⟹ ΣM ≡ 3；Σcos²=3/2（1+ω²+ω⁴=0）⟹ ΣM² ≡ 6。
# 由 √m ∝ |M|：K = ΣM²/(Σ|M|)² = 6/(3+2|S₋|)²，S₋ = 负振幅之和。
print("\n[B] 机制的恒等签名（与偏移无关，但**只在三个振幅同号时**给出 Koide=2/3）")
for deg in (0.0, 90.0, 132.7324, 147.8055, 200.5):
    v = [M(math.radians(deg) + 2*PI*i/3) for i in range(3)]
    ck(abs(sum(v)-3) < 1e-12, "B: δ=%.1f° ΣM = %.12f ≡ 3" % (deg, sum(v)))
    ck(abs(sum(x*x for x in v)-6) < 1e-12, "B: δ=%.1f° ΣM² = %.12f ≡ 6" % (deg, sum(x*x for x in v)))
v0 = [M(0.0 + 2*PI*i/3) for i in range(3)]
ck(all(x > 0 for x in v0), "B1a: δ=0° 三振幅同正 ⟹ K = %.12f = 2/3"
   % (sum(x*x for x in v0)/(sum(abs(x) for x in v0))**2))
v9 = [M(math.radians(90) + 2*PI*i/3) for i in range(3)]
K9 = sum(x*x for x in v9)/(sum(abs(x) for x in v9))**2
Sneg = abs(sum(x for x in v9 if x < 0))
ck(abs(K9 - 6/(3+2*Sneg)**2) < 1e-12, "B1b: δ=90° 混合符号 ⟹ K = 6/(3+2|S₋|)² = %.6f ≠ 2/3" % K9)

# ── C. 机制能否给出中微子形状？反解 X ──
print("\n[C] 反解：X = X_obs 的 δ（正确口径）")
sols, prev, pt, N = [], None, 0.0, 400000
for k in range(N+1):
    t = 2*PI*k/N
    X = X_of(t)
    if prev is not None and math.isfinite(X) and math.isfinite(prev) and (prev-X_OBS < 0) != (X-X_OBS < 0):
        a, b = pt, t
        for _ in range(90):
            mid = (a+b)/2
            if (X_of(mid)-X_OBS < 0) == (X_of(a)-X_OBS < 0): a = mid
            else: b = mid
        if all(abs((a+b)/2 - x) > 1e-6 for x in sols): sols.append((a+b)/2)
    prev, pt = X, t
d0s = sorted({round(delta0(t), 4) for t in sols})
print("     解的 δ₀^ν 取值：%s" % "、".join("%.4f°" % d for d in d0s))
mm, A = masses_of(sols[0])
K_mech = K_of(sols[0])
print("     m = [%.4f, %.4f, %.4f] meV；Σm_ν = %.3f meV；m1/m3 = %.5f"
      % (mm[0]*1e3, mm[1]*1e3, mm[2]*1e3, sum(mm)*1e3, mm[0]/mm[2]))
print("     K = %.6f" % K_mech)
ck(len(d0s) == 3, "C1: 三个不同的 δ₀^ν：%s" % d0s)
ck(abs(K_mech - 0.524368) < 1e-5, "C2: 拟合处的 K = %.6f" % K_mech)
ck(all(abs(K_of(t) - K_mech) < 1e-9 for t in sols), "C3: 全部解给出同一个 K（振幅符号皆为混合）")

# ── D. 数据侧 K 的允许范围，看机制是否相容 ──
print("\n[D] 相容性：机制的 K 是否落在数据允许范围内")
def K_data(ml, order):
    if order == "NO":
        ms = sorted([ml, math.sqrt(ml*ml+D21), math.sqrt(ml*ml+D31)])
    else:
        m3 = ml; m1 = math.sqrt(m3*m3+D31)
        ms = sorted([m1, math.sqrt(m1*m1+D21), m3])
    return sum(ms)/(sum(math.sqrt(x) for x in ms))**2
K0 = K_data(0.0, "NO")
KMAX = max(max(K_data(10**(-9 + k*2.0/100000), o) for k in range(100001)) for o in ("NO", "IO"))
print("     数据侧：m1→0（正规序）K = %.6f；全域上界 = %.6f" % (K0, KMAX))
print("     机制：K = %.6f  ⟹ %s" % (K_mech, "落在允许范围内（**相容**）" if K_mech < KMAX else "超出（不相容）"))
ck(abs(K0 - 0.585605) < 1e-5, "D1: K_data(m1=0) = %.6f" % K0)
ck(K_mech < KMAX, "D2: **机制与数据相容**：%.6f < %.6f（机制未被排除）" % (K_mech, KMAX))
ck(abs(K_mech - 2/3) > 0.1, "D3: 拟合处的 K 远离 2/3（因振幅混合符号）")

# ── E. 条件预言 ──
print("\n[E] 条件预言（机制拟合 X_obs 后，其余量成为预言）")
Ue = [math.sqrt(1 - 0.5443**2 - 0.1492**2), 0.5443, 0.1492]
mbeta = math.sqrt(sum(u*u*x*x for u, x in zip(Ue, mm)))
print("     Σm_ν = %.3f meV（Z019[A] 定标给 57.40 meV；宇宙学上限 %.0f meV）"
      % (sum(mm)*1e3, COSMO*1e3))
print("     m_β = %.3f meV；m1/m3 = %.5f；K = %.6f" % (mbeta*1e3, mm[0]/mm[2], K_mech))
ck(abs(sum(mm) - 0.059080) < 1e-5, "E1: Σm_ν = %.4f eV" % sum(mm))
ck(sum(mm) < COSMO, "E2: Σm_ν 在宇宙学上限内")
ck(abs(mbeta*1e3 - 8.831) < 0.05, "E3: m_β = %.3f meV" % (mbeta*1e3))

# ── F. 代价：需要的偏移是否内生 ──
print("\n[F] 代价：δ₀^ν 是否内生角")
ENDO = {"π/16": 180/16, "π/14": 180/14, "π/12": 15.0, "3π/16": 33.75, "π/8": 22.5,
        "π/4": 45.0, "π/3": 60.0, "3π/8": 67.5, "π/5": 36.0, "2π/9": 40.0,
        "3π/5": 108.0, "π/6": 30.0}
for d in d0s:
    b = min(ENDO.items(), key=lambda kv: abs(d/kv[1]-1))
    print("     δ₀^ν = %8.4f°  → 最近内生角 %s = %.2f°（差 %.2f%%）"
          % (d, b[0], b[1], abs(d/b[1]-1)*100))
best = min((abs(d/v-1), nm, d, v) for d in d0s for nm, v in ENDO.items())
ck(best[0] < 0.01, "F1: 最近内生角 %s 差 %.2f%% ⟹ 落在 1%% 级（凑数区）" % (best[1], best[0]*100))
print("     ⚠ 该 %s 含因子 5，与本框架的 8 分／3 分结构无来源关系 ⟹ 不足以采信（见 G）" % best[1])
ck(5 in [5] and best[1] == "3π/5", "F2: 最近的恰是 3π/5（分母 5 无结构来源）")

# ── G. 选择性检验：把角限制在内生角网格 ──
print("\n[G] 选择性检验（纪律：1% 级接近不构成证据）")
cands = sorted({round(k*PI/n, 12) for n in (8, 12, 16, 3, 4, 6, 24) for k in range(1, n)})
xs = []
for d in cands:
    if not (0 < d < PI): continue
    for t in cands:
        a, b = M(t+d), M(t+2*d)
        if abs(b) < 1e-6: continue
        r2 = (a/b)**2
        if r2 > 0: xs.append(r2)
xs.sort()
def hits(tg, tol=0.01): return sum(1 for v in xs if abs(v/tg-1) < tol)
h = hits(X_OBS)
self_counts = sorted(hits(v) for v in xs)
med = self_counts[len(self_counts)//2]
q1 = self_counts[len(self_counts)//4]
q3 = self_counts[3*len(self_counts)//4]
pct = 100*sum(1 for c in self_counts if c >= h)/len(self_counts)
print("     网格 %d 组；X_obs 的 ±1%% 邻域命中 %d 组；零假设下中位数 %d（IQR %d–%d，最大 %d）"
      % (len(xs), h, med, q1, q3, self_counts[-1]))
print("     ⟹ X_obs 排在前 %.1f%%" % pct)
ck(pct > 10, "G1: X_obs 未超 10%% 分位 ⟹ 网格的接近属典型值，不采信")

# ── H. 决定性检验：把 Z123 的**内生**结构用于中微子 ──
print("\n[H] 把内生结构用于中微子（偏移不是自由参数 ⟹ 这是**否证**而非拟合）")
Xe, _ = X_of(TH_E), None
print("    带电轻子：√m 比 [1, %.2f, %.2f] vs 观测 [1, 14.38, 58.968]；Koide = %.6f"
      % (sorted(abs(M(t)) for t in lepton_angles())[1]/sorted(abs(M(t)) for t in lepton_angles())[0],
         sorted(abs(M(t)) for t in lepton_angles())[2]/sorted(abs(M(t)) for t in lepton_angles())[0],
         K_of(TH_E)))
ck(abs(K_of(TH_E) - 2/3) < 1e-12, "H1: 带电轻子 Koide = 2/3 精确（振幅同号）")
ck(abs(d0e - 0.039577166) < 1e-8, "H2: δ₀ = π/12 − 2/9 = %.9f（Z123，完全内生）" % d0e)
Xend = X_of(TH_E)
print("    中微子（原样）  X = %.6f  vs 实测 %.6f ⟹ 差 %.2f 倍" % (Xend, X_OBS, X_OBS/Xend))
ck(abs(Xend - 0.003535) < 1e-5, "H3: 内生结构给 X = %.6f" % Xend)
ck(X_OBS/Xend > 8, "H4: **否证**：实测比内生预言大 %.2f 倍" % (X_OBS/Xend))
for sh, tag in ((PI, "加 π（另一零点侧）"), (-PI, "减 π")):
    Xs = X_of(TH_E + sh)
    print("    中微子（%s）X = %.6f ⟹ 差 %.2f 倍" % (tag, Xs, X_OBS/Xs if Xs else float('inf')))
    ck(abs(Xs/X_OBS - 1) > 0.05, "H5: %s 亦不符（因子 %.2f）" % (tag, Xs/X_OBS))
# 拟合路线的代价：需要的 δ₀^ν
need = min((abs(delta0(t) - 0.0), delta0(t)) for t in sols)[1]
print("    若放弃内生偏移、改为拟合：需要 δ₀^ν ∈ {%s}（都不是框架量）"
      % "、".join("%.2f°" % d for d in d0s))
ck(all(abs(d - d0e) > 1e-6 for d in d0s), "H6: 拟合所需的三个偏移都不等于内生值 %.4f°" % math.degrees(d0e))

print("\n=== 合计：%d 条断言 / %d 不符 ===" % (NCHK, FAIL))
raise SystemExit(0 if FAIL == 0 else 1)
