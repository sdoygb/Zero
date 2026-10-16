#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Z019_neutrino_shape_check.py —— 中微子形状能否由后来的质量机制（Z121/Z123）给出？

问题（用户 2026-10-16）：Z019 是早期篇，此后理论有长足发展；中微子质量能否算准？

【术语（先消歧，否则极易混淆——本文作者第一版就在这里错了）】
  X := Δm²21/Δm²31 = (m2²−m1²)/(m3²−m1²)   ← **可观测量**，尺度无关，故是机制的**纯预言**
  R := m2/m3                                 ← **质量比**；仅当 m1=0 时 X = R²
  Z019 采用形状 R = (73−28√6)/25 = 0.176571488083（Z019 §8 边界一：**来源不在本系列**）

【机制（Z121 立、Z123 补）】
  m = A²·M(θ)²,  M(θ) = 1+√2·cos θ;  三代间隔 2π/3（来自"三项"）；
  零点 3π/4、振幅 √2 来自闭合 2π 分 8 份（seed 8）。

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
def M(t): return 1 + S2 * math.cos(t)

# ── 数据（正序全球拟合，量级与 Z019 同源）──
D21, D31 = 7.42e-5, 2.510e-3          # eV²
X_OBS = D21 / D31                     # 可观测量 0.029562
R019 = (73 - 28*math.sqrt(6)) / 25    # Z019 采用的质量比 0.176571

ZP = (3*PI/4, 5*PI/4)                 # 两个零点（M=0）

def masses(t1):
    """给定第一代角，按质量升序返回 m²（相对值）。"""
    return sorted(M(t1 + 2*PI*i/3)**2 for i in range(3))

def X_of(t1):
    ms = masses(t1)
    if ms[2] - ms[0] < 1e-15: return float('inf')
    return (ms[1] - ms[0]) / (ms[2] - ms[0])

def delta0(t1):
    """离最近零点的角距（度）。"""
    return math.degrees(min(min(abs(t1-z), abs(t1-z+2*PI), abs(t1-z-2*PI)) for z in ZP))

print("=== 0. 三个数先对清楚 ===")
print("  X_obs（可观测量）           = %.6f" % X_OBS)
print("  Z019 的形状 R               = %.9f  ⟹ X = %.6f" % (R019, R019**2))
try:
    from math import isclose
except ImportError:
    pass
ck(abs(R019**2 - 0.0311774904) < 1e-9, "0a: R019² = %.9f" % R019**2)
ck(abs(X_OBS - 0.029562) < 5e-7, "0b: X_obs = %.6f（实测）" % X_OBS)
ck(abs(R019**2/X_OBS - 1 - 0.05466) < 1e-3,
   "0c: Z019 形状比实测高 %.2f%%（＝Z019 §3.3 的 [B] 定标张力 +5.47%%）" % ((R019**2/X_OBS-1)*100))

# ── A. 机制在带电轻子侧的自检 ──
print("\n[A] 机制在带电轻子侧（若这里不成立，后面无需谈）")
d0e = PI/12 - 2/9                     # Z123：δ₀ 完全内生
te = 3*PI/4 - d0e
Mm = [M(te + 2*PI*i/3) for i in range(3)]
r = [x/Mm[0] for x in Mm]
obs = [1.0, 14.38, 58.968]
print("     √m 比 = [1, %.4f, %.4f]  vs 观测 [1, %.2f, %.3f]" % (r[1], r[2], obs[1], obs[2]))
ck(abs(r[1]/obs[1]-1) < 5e-5, "A1: μ 偏差 %.4f%%" % (abs(r[1]/obs[1]-1)*100))
ck(abs(r[2]/obs[2]-1) < 5e-5, "A2: τ 偏差 %.4f%%" % (abs(r[2]/obs[2]-1)*100))
ck(abs(d0e - 0.039577166) < 1e-8, "A3: δ₀^e = %.9f（带电轻子的输入，非内生角）" % d0e)

# ── B. 中微子侧：Z019 自己的设定 m1=0 ⟹ δ₀^ν=0 ──
print("\n[B] 保留 Z019 的设定 m1=0 ⟹ 最轻态落在零点 ⟹ δ₀^ν = 0")
mb = masses(3*PI/4)
Xm = X_of(3*PI/4)
print("     θ1=3π/4: M = [%.6f, %.6f, %.6f]" % (M(3*PI/4), M(3*PI/4+2*PI/3), M(3*PI/4+4*PI/3)))
print("     ⟹ X = %.6f，实测 %.6f ⟹ 机制比实测高 %.2f 倍" % (Xm, X_OBS, Xm/X_OBS))
ck(abs(Xm - 0.071797) < 1e-5, "B1: 机制(m1=0) 预言 X = %.6f" % Xm)
ck(Xm/X_OBS > 2, "B2: **被实测否证**：高 %.2f 倍（远超 8.14%% 那处张力）" % (Xm/X_OBS))
ck(abs(mb[0]) < 1e-12, "B3: 该设定下 m1 恰为 0（与 Z019 一致）")

# ── C. 一般情形：允许 δ₀^ν 自由，反解 ──
print("\n[C] 一般情形：δ₀^ν 自由，反解 X = X_obs")
sols, prev, N = [], 0.0, 400000
fp = X_of(1e-9) - X_OBS
for k in range(1, N+1):
    x = 2*PI*k/N
    fx = X_of(x) - X_OBS
    if math.isfinite(fx) and (fp < 0) != (fx < 0):
        a, b = prev, x
        for _ in range(90):
            mid = (a+b)/2
            if (X_of(a)-X_OBS < 0) != (X_of(mid)-X_OBS < 0): b = mid
            else: a = mid
        if all(abs((a+b)/2 - s) > 1e-6 for s in sols): sols.append((a+b)/2)
    prev, fp = x, fx
d0s = sorted({round(delta0(t), 6) for t in sols})
print("     解的 δ₀^ν 取值：%s" % "、".join("%.4f°" % d for d in d0s))
for t in sols:
    ms = masses(t)
    print("     θ1=%.4fπ  δ₀^ν=%7.4f°   m1/m3=%.4f" % (t/PI, delta0(t), math.sqrt(ms[0]/ms[2])))
ck(len(d0s) >= 1, "C1: 反解有解（δ₀^ν 取 %d 个不同值）" % len(d0s))
ENDO = {"π/16": 180/16, "π/14": 180/14, "π/12": 180/12, "π/8": 180/8,
        "3π/16": 3*180/16, "π/4": 180/4, "π/10": 180/10, "2π/27": 360/27}
best = min((abs(d/v - 1), nm, d, v) for d in d0s for nm, v in ENDO.items())
print("     最接近的内生角：%s（%s=%.4f° vs 需 %.4f°），差 %.2f%%" % (best[1], best[1], best[3], best[2], best[0]*100))
ck(best[0] > 0.005, "C2: 需要的 δ₀^ν 与任何内生角都差 >0.5%% ⟹ **不是内生角**（最好者差 %.2f%%）" % (best[0]*100))
ck(all(d > 0.05 for d in d0s), "C3: 所有解都要求 δ₀^ν ≠ 0 ⟹ 与 Z019 的 m1=0 不相容")

# ── D. 选择性检验：把角限制在内生角网格，看能否"凑上" ──
print("\n[D] 选择性检验（本项目纪律：1% 级接近不构成证据）")
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
bx = min(xs, key=lambda v: abs(v/X_OBS - 1))
print("     网格 %d 组；最接近 X_obs 者 X=%.6f（差 %+.2f%%）" % (len(xs), bx, (bx/X_OBS-1)*100))

def hits(tg, tol=0.01): return sum(1 for v in xs if abs(v/tg-1) < tol)
h = hits(X_OBS)
# 零假设：**拿网格值自身当目标**（不是从任意区间抽随机目标——那样基线会被
# 网格自身的集中度误导，本文作者第一版即栽在这里）
self_counts = sorted(hits(v) for v in xs)
med = self_counts[len(self_counts)//2]
q1 = self_counts[len(self_counts)//4]
q3 = self_counts[3*len(self_counts)//4]
pct = 100*sum(1 for c in self_counts if c >= h)/len(self_counts)
print("     X_obs 的 ±1%% 邻域命中 %d 组；零假设下 中位数 %d（IQR %d–%d，最大 %d）"
      % (h, med, q1, q3, self_counts[-1]))
print("     ⟹ X_obs 排在前 %.1f%%" % pct)
ck(pct > 10, "D1: X_obs 未超过 10%% 分位 ⟹ 网格里那 %.2f%% 的接近是**典型值，非证据（凑数）**"
   % abs(bx/X_OBS-1)*100)

# ── E. 若机制适用且偏移自由：给出可否证的预言 ──
print("\n[E] 条件预言（若机制适用于中微子扇区、δ₀^ν 取最小解 %.4f°）" % d0s[0])
t1 = min(sols, key=delta0)
vs = masses(t1)
A2 = D31 / (vs[2] - vs[0])
mm = [math.sqrt(A2*v) for v in vs]            # eV
Ue = [math.sqrt(1 - 0.5443**2 - 0.1492**2), 0.5443, 0.1492]
mbeta = math.sqrt(sum(u*u*x*x for u, x in zip(Ue, mm)))
summ = sum(mm)
print("     m = [%.4f, %.4f, %.4f] meV；Σm_ν = %.3f meV；m_β = %.4f meV"
      % (mm[0]*1e3, mm[1]*1e3, mm[2]*1e3, summ*1e3, mbeta*1e3))
print("     m1/m3 = %.4f ⟹ **最轻态不零质量**（与 Z019 的 m1=0 冲突）" % (mm[0]/mm[2]))
ck(abs(mm[0]/mm[2] - 0.0644) < 5e-4, "E1: m1/m3 = %.4f（非零）" % (mm[0]/mm[2]))
ck(abs(summ - 0.062636) < 1e-5, "E2: Σm_ν = %.4f eV ⟹ 紧贴宇宙学上限（Z019 引用 0.064 eV）" % summ)
ck(abs(mbeta - 0.0093969) < 1e-6, "E3: m_β = %.4f meV" % (mbeta*1e3))

# ── F. 机制的恒等签名：Koide ≡ 2/3 ──
# 因 θ_i 是间隔 2π/3 的等分点：Σcos(θ+2πk/3)=0（1+ω+ω²=0）⟹ ΣM ≡ 3；
# 又 Σcos²=3/2（1+ω²+ω⁴=0）⟹ ΣM² ≡ 6。故 m∝M² 时
#   K = Σm/(Σ√m)² = ΣM²/(Σ|M|)² ；当三个 M 同号时 Σ|M| = |ΣM| = 3 ⟹ K = 6/9 = 2/3。
print("\n[F] 机制的恒等签名：Koide ≡ 2/3（与偏移无关，只要三个 M 同号）")
for deg in (0.0, 20.0, 90.0, 132.7324, 147.8055, 200.5):
    v = [M(math.radians(deg) + 2*PI*i/3) for i in range(3)]
    ck(abs(sum(v) - 3) < 1e-12, "F: δ=%.1f° 时 ΣM = %.12f ≡ 3" % (deg, sum(v)))
    ck(abs(sum(x*x for x in v) - 6) < 1e-12, "F: δ=%.1f° 时 ΣM² = %.12f ≡ 6" % (deg, sum(x*x for x in v)))
v = [M(0.0 + 2*PI*i/3) for i in range(3)]          # δ=0：M=[2.414,0.293,0.293] 三正
ck(all(x > 0 for x in v), "F1a: δ=0° 时三个 M 同为正（M=[%.4f,%.4f,%.4f]）" % tuple(v))
ms = sorted(x*x for x in v)
Kd = sum(ms)/(sum(math.sqrt(x) for x in ms))**2
ck(abs(Kd - 2/3) < 1e-12, "F1b: 同号时 Koide = %.12f = 2/3" % Kd)
v = [M(math.radians(90) + 2*PI*i/3) for i in range(3)]   # δ=90°：混合符号
ms = sorted(x*x for x in v)
Km = sum(ms)/(sum(math.sqrt(x) for x in ms))**2
ck(abs(Km - 0.504244923464) < 1e-9,
   "F1c: δ=90° 混合符号（M=[%.4f,%.4f,%.4f]）⟹ Koide = %.6f ≠ 2/3" % (v[0], v[1], v[2], Km))
# 一般式 K = 6/(3+2|S₋|)²（S₋ = 负值之和）
for deg in (250.0, 300.0):
    v = [M(math.radians(deg) + 2*PI*i/3) for i in range(3)]
    Sneg = abs(sum(x for x in v if x < 0))
    ms = sorted(x*x for x in v)
    K = sum(ms)/(sum(math.sqrt(x) for x in ms))**2
    ck(abs(K - 6/(3+2*Sneg)**2) < 1e-12, "F2: δ=%.0f° 时 K = 6/(3+2|S₋|)² = %.6f" % (deg, K))

# ── G. 用签名检验中微子扇区（决定性）──
print("\n[G] 决定性检验：机制在 X 拟合处把 K 钉死在 2/3，而数据要求 K ≤ 0.5855")
def K_data(ml, order):
    if order == "NO":
        ms = sorted([ml, math.sqrt(ml*ml+D21), math.sqrt(ml*ml+D31)])
    else:
        m3 = ml; m1 = math.sqrt(m3*m3+D31)
        ms = sorted([m1, math.sqrt(m1*m1+D21), m3])
    return sum(ms)/(sum(math.sqrt(x) for x in ms))**2
KMAX = max(max(K_data(10**(-9 + k*2.0/100000), o) for k in range(100001)) for o in ("NO", "IO"))
print("     数据侧 K 上界（两序、最轻质量全域）= %.6f" % KMAX)
ck(abs(KMAX - 0.585489) < 1e-5, "G1: K_max(data) = %.6f" % KMAX)
ck(2/3 > KMAX, "G2: 2/3 = %.6f > %.6f，差 %+.1f%% ⟹ 签名被违反" % (2/3, KMAX, (2/3/KMAX-1)*100))
# 扫描：X 在 X_obs ±5% 内的全部 δ，机制的 K 是否恒定在 2/3
ks, t1s = [], []
NN = 400000
for k in range(NN+1):
    t = 2*PI*k/NN
    v = [M(t + 2*PI*i/3) for i in range(3)]
    mss = sorted(x*x for x in v)
    if mss[2] - mss[0] < 1e-15: continue
    Xv = (mss[1]-mss[0])/(mss[2]-mss[0])
    if abs(Xv/X_OBS - 1) < 0.05:
        ks.append(sum(mss)/(sum(math.sqrt(x) for x in mss))**2); t1s.append(t)
print("     X 落在 X_obs ±5%% 的 δ 共 %d 个采样：K ∈ [%.6f, %.6f]（恒定）" % (len(ks), min(ks), max(ks)))
ck(len(ks) > 1000, "G3: 采样点数 %d > 1000" % len(ks))
ck(max(ks) - min(ks) < 1e-12, "G4: 该区间的 K **恒定**（极差 %.2e）" % (max(ks)-min(ks)))
ck(min(ks) > KMAX, "G5: **排除**：机制能配 X 的每一处，其 K = %.6f 都超过数据上界 %.6f" % (min(ks), KMAX))
# 机制的 K 全域（说明它确实能给低 K，只是那时 X 不对）
allk = []
for k in range(0, 200000):
    t = 2*PI*k/200000
    v = [M(t + 2*PI*i/3) for i in range(3)]
    mss = sorted(x*x for x in v)
    if mss[2]-mss[0] < 1e-15: continue
    allk.append(sum(mss)/(sum(math.sqrt(x) for x in mss))**2)
print("     （机制 K 全域 = [%.6f, %.6f]：它也能给低 K，但那时 X 远离 X_obs）"
      % (min(allk), max(allk)))
ck(min(allk) < KMAX, "G6: 机制全域 K 低至 %.6f ⟹ 排除来自 X 与 K 的不相容，不是 K 值域不够" % min(allk))

print("\n=== 合计：%d 条断言 / %d 不符 ===" % (NCHK, FAIL))
raise SystemExit(0 if FAIL == 0 else 1)
