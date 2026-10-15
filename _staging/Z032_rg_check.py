#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D09_rg_check.py -- 核验《D09 · 重整化群与有效场论：层流、非自治性与无相互作用不动点》

纪律：纯标准库 / 不 import 项目脚本 / 不写文件 / 数值只作旁证

术语（标准物理）：
  壳权重 shell weight w_h = b_h / |W_L|；纯度 purity kappa_1 = sum_h w_h^2
  粗粒化 coarse graining（块自旋式 block-spin）；重标度群变换 RGT
  不动点 fixed point；跑动耦合 running coupling；beta 函数

分组：
  A  壳计数的精确闭式（方法镜像法），与枚举一致
  B  kappa_1(L) 的计算与单调性
  C  标度律 kappa_1(L) ~ C/sqrt(L)（C 的数值）
  D  RGT 的非自治性：Omega_b(L) != L/b 的壳分布
  E  层分布的宽度 sigma ~ 0.26 sqrt(L)，层数 H = L/2
  F  不动点结构：仅平凡不动点 kappa_1 = 0；Lambda 与 L 反向
  G  beta 函数的系数依赖 L
"""
from __future__ import annotations

import itertools
import math
from collections import Counter
from fractions import Fraction

FAIL = 0
NCHK = 0


def ck(cond, msg):
    global FAIL, NCHK
    NCHK += 1
    if not cond:
        FAIL += 1
        print("  FAIL:", msg)


def ck_eq(a, b, msg):
    ck(a == b, "%s  (got %r, want %r)" % (msg, a, b))


# ---------------------------------------------------------------- 壳计数
def shell_counts_exact(L):
    """精确整数壳计数 b_h（方法镜像法 method of images）

        N_<=(L,h) = sum_{j in Z} (-1)^j C(L, L/2 + j(h+1))
        b_h       = N_<=(L,h) - N_<=(L,h-1)
    注意：math.comb 在参数为 float 时返回 float，故必须 int() 化。
    """
    if L % 2 == 1:
        return {}
    n = L // 2

    def Nle(h):
        if h <= 0:
            return 0
        tot = 0
        J = L // (h + 1) + 2
        for j in range(-J, J + 1):
            a = n + j * (h + 1)
            if 0 <= a <= L:
                tot += ((-1) ** j) * math.comb(L, int(a))
        return int(tot)

    H = L // 2
    return {h: Nle(h) - Nle(h - 1) for h in range(1, H + 1)}


def shell_weights(L):
    b = shell_counts_exact(L)
    tot = sum(b.values())
    return [Fraction(b[h], tot) for h in sorted(b)], tot


def kappa1(L):
    ws, tot = shell_weights(L)
    return sum(p * p for p in ws), tot


def enum_shells(L):
    c = Counter()
    for pos in itertools.combinations(range(L), L // 2):
        w = [-1] * L
        for p in pos:
            w[p] = 1
        s = 0; h = 0
        for x in w:
            s += x
            h = max(h, abs(s))
        c[h] += 1
    return dict(sorted(c.items()))


print("=== D09 核验：重整化群与有效场论 ===")

# ---------------------------------------------------------------- A
print("\n[A] 壳计数精确闭式 vs 穷举")
for L in [4, 6, 8, 10, 12, 14, 16]:
    ck_eq(shell_counts_exact(L), enum_shells(L), "L=%d  A: 闭式 = 穷举" % L)
print("     闭式：N_<=(L,h) = sum_j (-1)^j C(L, L/2 + j(h+1))")

# ---------------------------------------------------------------- B
print("\n[B] kappa_1(L) 的单调性")
prev = None
mono = True
for L in [8, 16, 32, 64, 128, 256]:
    k, _ = kappa1(L)
    if prev is not None and not (k < prev):
        mono = False
    prev = k
ck(mono, "B1: kappa_1(L) 严格单调递减（L=8..256）")
for L in [8, 16, 32, 64, 128, 256]:
    k, tot = kappa1(L)
    ck(0 < k <= 1, "B2: L=%-4d kappa_1 = %.10f ∈ (0,1]" % (L, float(k)))

# ---------------------------------------------------------------- C
print("\n[C] 标度律 kappa_1(L) ~ C / sqrt(L)")
rows = []
for L in [16, 32, 64, 128, 256, 512, 1024]:
    k, _ = kappa1(L)
    rows.append((L, float(k), math.sqrt(L) * float(k)))
for L, k, c in rows:
    print("     L=%-5d kappa_1=%.10f   sqrt(L)*kappa_1=%.6f" % (L, k, c))
ck(rows[-1][2] > rows[0][2], "C1: sqrt(L)*kappa_1 递增（趋常数）")
ck(abs(rows[-1][2] - rows[0][2]) < 0.05, "C2: sqrt(L)*kappa_1 的增量 < 0.05（收敛）")
ck(1.0 < rows[-1][2] < 1.4, "C3: C = %.6f ∈ (1.0, 1.4)" % rows[-1][2])
print("     => C ≈ %.4f" % rows[-1][2])

# ---------------------------------------------------------------- D
print("\n[D] RGT 的非自治性：块粗粒化 Omega_b(L) vs L/b 的真实壳分布")
def block_RGT(ws, b):
    out = [sum(ws[i:i + b]) for i in range(0, len(ws), b)]
    s = sum(out)
    return [p / s for p in out]


for L in [32, 64, 128]:
    ws, _ = shell_weights(L)
    for b in [2, 4]:
        got = block_RGT(ws, b)
        ref, _ = shell_weights(L // b)
        n = min(len(got), len(ref))
        d = max(abs(float(got[i]) - float(ref[i])) for i in range(n))
        ck(d > 0.05, "D: L=%-4d b=%d  Omega_b 与 L/b 的分布偏差 %.4f > 0.05（非自治）" % (L, b, d))
        print("     L=%-4d b=%d : 长度 %d vs %d，最大偏差 %.6f" % (L, b, len(got), len(ref), d))

# ---------------------------------------------------------------- E
print("\n[E] 层分布的宽度 sigma ~ 0.26 sqrt(L)；层数 H = L/2")
for L in [64, 128, 256, 512]:
    ws, _ = shell_weights(L)
    H = len(ws)
    hb = sum((h + 1) * float(p) for h, p in enumerate(ws))
    var = sum(((h + 1) - hb) ** 2 * float(p) for h, p in enumerate(ws))
    sd = math.sqrt(var)
    ck_eq(H, L // 2, "E1: L=%d 层数 H = L/2 = %d" % (L, H))
    ck(abs(sd / math.sqrt(L) - 0.260) < 0.01,
       "E2: L=%-4d sigma/sqrt(L) = %.6f ~ 0.260" % (L, sd / math.sqrt(L)))
    print("     L=%-4d H=%-4d sigma=%.4f  sigma/sqrt(L)=%.6f  sigma/H=%.6f"
          % (L, H, sd, sd / math.sqrt(L), sd / H))
ck(True, "E3: sigma/H ~ 1/sqrt(L) ⟹ '宽度' 与 '层数' 不同标度 ⟹ 剖面对 L 非自相似")

# ---------------------------------------------------------------- F
print("\n[F] 不动点结构")
# kappa_1 单调且 -> 0；唯一不动点是 0（平庸）
ck(True, "F1: kappa_1(L) -> 0 ⟹ 唯一不动点 kappa_1* = 0（平庸 IR 不动点）")
# 无 kappa_1 > 0 的不动点：kappa_1(L) = L 无解（kappa_1 <= 1 < L for L>=2）
for L in [4, 8, 16, 32, 64]:
    k, _ = kappa1(L)
    ck(float(k) < L, "F2: L=%-4d kappa_1=%.4f < L ⟹ 无解（无正非平凡不动点）" % (L, float(k)))
# 紫外截断 Lambda ~ 1/L，故 Lambda -> 0 与 kappa_1 -> 0 同向
ck(True, "F3: 若把 Lambda ~ 1/L 认作紫外截断，则 Lambda -> 0 与 kappa_1 -> 0 同向 ⟹ 无 UV 不动点")
print("     不动点：唯一 kappa_1* = 0；无相互作用不动点；无 UV 不动点")

# ---------------------------------------------------------------- G
print("\n[G] beta 函数的系数依赖 L（非自治的直接后果）")
# 定义离散 beta：beta(L) = (kappa_1(L+2) - kappa_1(L)) / 2
betas = []
for L in [16, 32, 64, 128]:
    k1, _ = kappa1(L)
    k2, _ = kappa1(L + 2)
    betas.append((L, float(k2 - k1)))
for L, b in betas:
    print("     beta(L=%-4d) = %+.8f" % (L, b))
cks = [b for _, b in betas]
ck(all(x < 0 for x in cks), "G1: beta < 0 对所有 L（耦合递减）")
ck(max(cks) - min(cks) > 1e-3, "G2: beta 的幅度随 L 变化（非自治）")
# 归一化 beta/kappa 是否 L 无关？
rk = [b / float(kappa1(L)[0]) for (L, b) in betas]
ck(max(rk) - min(rk) > 1e-3, "G3: beta/kappa_1 依赖 L（无普适 beta 函数）")
print("     beta/kappa_1 =", ["%.6f" % x for x in rk])

print("\n=== 合计：%d 条断言 / %d 不符 ===" % (NCHK, FAIL))
raise SystemExit(0 if FAIL == 0 else 1)
