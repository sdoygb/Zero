#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D11_qinfo_check.py -- 核验《D11 · 量子信息与开放系统》

纪律：纯标准库 / 不 import 项目脚本 / 不写文件 / 数值只作旁证

术语（标准物理）：中心态 central state；层间相干 inter-layer coherence；
  完全正映射 completely positive map；Kraus 表示 Kraus representation；
  量子信道 quantum channel；Lindblad 主方程；马尔可夫近似 Markov approximation；
  熵 entropy；时间反演 time reversal；反幺正算符 antiunitary operator；
  模共轭 modular conjugation。

分组：
  A  中心态的结构：块内最大混合
  B  层保持信道：不产生层间相干
  C  随机 Kraus 映射破坏中心性（故"噪声"在本框架内不可自由定义）
  D  熵的可加分解 S = H(w) + ln2
  E  两种翻转变换对相干的作用（共轭 vs 恒等）
  F  时间反演：J0 对 sy 敏感；换基后 J 改变（不典范）
"""
from __future__ import annotations

import math
import random
import re
from fractions import Fraction

FAIL = 0
NCHK = 0


def ck(cond, msg):
    global FAIL, NCHK
    NCHK += 1
    if not cond:
        FAIL += 1
        print("  FAIL:", msg)


# ---------------------------------------------------------------- 线性代数（复数小矩阵）
def mm(X, Y):
    n, k, m = len(X), len(Y), len(Y[0])
    return [[sum(X[i][t] * Y[t][j] for t in range(k)) for j in range(m)] for i in range(n)]


def dag(X):
    n, m = len(X), len(X[0])
    return [[X[j][i].conjugate() for j in range(n)] for i in range(m)]


def rho_central(w, r=0):
    """4x4：rho = (w0/2) I2 (+) (w1/2) I2，块外相干 r（在 (0,2) 与 (1,3) 位）"""
    return [[w[0] / 2, 0, r, 0],
            [0, w[0] / 2, 0, r],
            [r.conjugate() if hasattr(r, "conjugate") else r, 0, w[1] / 2, 0],
            [0, r.conjugate() if hasattr(r, "conjugate") else r, 0, w[1] / 2]]


def off_block(M):
    """跨块元的最大模（块 1 = 指标 0,1；块 2 = 指标 2,3）"""
    return max(abs(M[i][j]) for i in (0, 1) for j in (2, 3))


def intra_imbalance(M):
    return (abs(M[0][0] - M[1][1]), abs(M[2][2] - M[3][3]))


print("=== D11 核验：量子信息与开放系统 ===")
print("说明：全部为 4x4（= I2 (+) I2，两层的最小模型）上的可复算核验。")

# ---------------------------------------------------------------- A
print("\n[A] 中心态的结构：块内最大混合")
for w in [(Fraction(1, 2), Fraction(1, 2)), (Fraction(3, 10), Fraction(7, 10)),
          (Fraction(9, 10), Fraction(1, 10))]:
    M = rho_central([float(w[0]), float(w[1])])
    ck(abs(M[0][0] - M[1][1]) < 1e-15, "A1: w=%s 第一块内布居相等（最大混合）" % (w,))
    ck(abs(M[2][2] - M[3][3]) < 1e-15, "A2: w=%s 第二块内布居相等（最大混合）" % (w,))
    ck(off_block(M) == 0, "A3: w=%s 无层间相干（中心态）" % (w,))
    tr = sum(M[i][i] for i in range(4))
    ck(abs(tr - 1) < 1e-15, "A4: w=%s 迹为 1" % (w,))

# ---------------------------------------------------------------- B
print("\n[B] 层保持信道不产生层间相干")
LAYER_KEEPING = {
    "块标号翻转 I2⊗X": [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]],
    "块内翻转 X⊗I2": [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]],
    "对角相位": [[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1j, 0], [0, 0, 0, -1j]],
}
for name, K in LAYER_KEEPING.items():
    for r in [0.0, 0.1]:
        M = rho_central([0.3, 0.7], r)
        out = mm(mm(K, M), dag(K))
        ck(off_block(out) <= abs(r) + 1e-12,
           "B: %s r=%.1f 不放大层间相干（出 %.4f <= 入 %.4f）" % (name, r, off_block(out), abs(r)))

# ---------------------------------------------------------------- C
print("\n[C] 随机 Kraus 映射破坏中心性")
random.seed(11)
broke = 0
for _ in range(20):
    K = [[complex(random.uniform(-1, 1), random.uniform(-1, 1)) for _ in range(4)] for _ in range(4)]
    M = rho_central([0.3, 0.7])
    out = mm(mm(K, M), dag(K))
    tr = sum(out[i][i] for i in range(4))
    out = [[out[i][j] / tr for j in range(4)] for i in range(4)]
    i1, i2 = intra_imbalance(out)
    if off_block(out) > 1e-9 or i1 > 1e-9 or i2 > 1e-9:
        broke += 1
ck(broke >= 18, "C1: 20 个随机 Kraus 映射中 %d 个破坏中心性（≥18）" % broke)
print("     20 个随机映射中 %d 个把态推出中心 ⟹ 一般噪声在本框架内不可自由定义" % broke)

# ---------------------------------------------------------------- D
print("\n[D] 熵的可加分解 S = H(w) + ln2")
def H_shannon(w):
    return -sum(x * math.log(x) for x in w if x > 0)


for w in [(0.5, 0.5), (0.3, 0.7), (0.9, 0.1)]:
    S = -sum((x / 2) * math.log(x / 2) * 2 for x in w)
    ck(abs(S - (H_shannon(w) + math.log(2))) < 1e-12,
       "D1: w=%s S = H(w) + ln2" % (w,))
    ck(abs((S - math.log(2)) - H_shannon(w)) < 1e-12,
       "D2: w=%s 块间贡献恰为 ln2" % (w,))
    print("     w=%s : S = %.6f ; H(w) = %.6f ; ln2 = %.6f" % (w, S, H_shannon(w), math.log(2)))

# ---------------------------------------------------------------- E
print("\n[E] 两种翻转变换对层间相干的作用")
r = 0.1 + 0.05j
for name, K in [("块内翻转 X⊗I2", LAYER_KEEPING["块内翻转 X⊗I2"]),
                ("块标号翻转 I2⊗X", LAYER_KEEPING["块标号翻转 I2⊗X"])]:
    out = mm(mm(K, rho_central([0.3, 0.7], r)), dag(K))
    coh = out[0][2]
    if name.startswith("块内"):
        ck(abs(coh - r.conjugate()) < 1e-12, "E1: %s 把 r 变为 conj(r)" % name)
    else:
        ck(abs(coh - r) < 1e-12, "E1: %s 保持 r 不变" % name)
    ck(abs(out[0][0] - out[1][1]) < 1e-12 and abs(out[2][2] - out[3][3]) < 1e-12,
       "E2: %s 后仍为中心（I 在两种翻转下不变）" % name)
    print("     %-18s r=%s -> %s" % (name, r, complex(round(coh.real, 4), round(coh.imag, 4))))

# ---------------------------------------------------------------- F
print("\n[F] 时间反演候选：不典范")
I2 = [[1, 0], [0, 1]]
sz = [[1, 0], [0, -1]]
sx = [[0, 1], [1, 0]]
sy = [[0, -1j], [1j, 0]]


def J0(A):
    """标准基下的复共轭"""
    return [[A[i][j].conjugate() for j in range(2)] for i in range(2)]


# J0 固定 sz、sx，但翻转 sy（=> 对 sy 敏感）
def same(A, B):
    return all(abs(A[i][j] - B[i][j]) < 1e-12 for i in range(2) for j in range(2))


ck(same(J0(sz), sz), "F1: J0 固定 sz")
ck(same(J0(sx), sx), "F2: J0 固定 sx")
ck(not same(J0(sy), sy), "F3: J0 翻转 sy（对 sy 敏感）")
# 换基 U = diag(1, i) 后 J 改变
U = [[1, 0], [0, 1j]]
Ud = dag(U)


def Jp(A):
    return mm(mm(U, J0(mm(mm(Ud, A), U))), Ud)


ck(same(Jp(sz), J0(sz)), "F4: 换基后 sz 仍固定（不区分不了）")
ck(not same(Jp(sx), J0(sx)), "F5: 换基后 J 改变 ⟹ 不典范（对 sx 可检出）")
print("     J0 固定 sz/sx、翻转 sy；换基 U=diag(1,i) 后 J 改变 ⟹ 依赖基的选取")

print("\n=== 合计：%d 条断言 / %d 不符 ===" % (NCHK, FAIL))
raise SystemExit(0 if FAIL == 0 else 1)
