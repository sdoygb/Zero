#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z5_check.py —— 【有限 k：逃出三角困境与三笔代价】的核验
=====================================================
独立实断言：
  F1  (L') 修好：半径 k 之外远端扰动的影响**精确为 0**（1D 链，复刻 G41 §3 判据）
  F2  同一判据在 2D 方格上同样精确为 0
  F3  (O) 保住：w^(k) = A·A^{k-1} 的支撑仍只在最近邻
  F4  正则 Γ：有限 k 的体内权重**精确均匀**（故度规退化为平坦）
  F5  非正则 Γ：有限 k 的体内权重非均匀（几何由非正则性承载）
  F6  Z4 的“域基态剖面”是 Perron 假象（有限 k 下消失；Perron 下非均匀）
  F7  文档的结论、三笔代价与红标在位
"""
import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LEDGER_MODE = os.environ.get("LH_LEDGER", "A")
FAIL = 0
N = 0


def check(name, cond, detail="", level="ind"):
    global FAIL, N
    N += 1
    ok = bool(cond)
    if not ok:
        FAIL += 1
    tag = ("v" if ok else "x") if (level == "ind" or LEDGER_MODE == "A" or not ok) else "i"
    print("  [%s] %s%s" % (tag, name, ("   " + detail) if detail else ""))


def head(t):
    print("\n" + "=" * 72 + "\n" + t + "\n" + "=" * 72)


DOC = io.open(os.path.join(HERE, "Z5_finite_k_locality_escape.md"), encoding="utf-8").read()


def chain(n, mu=1.0, d=None):
    A = np.zeros((n, n))
    for i in range(n - 1):
        w = mu if i == d else 1.0
        A[i, i + 1] = A[i + 1, i] = w
    return A


def grid(n):
    Nn = n * n
    A = np.zeros((Nn, Nn))
    for x in range(n):
        for y in range(n):
            i = x * n + y
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                xx, yy = x + dx, y + dy
                if 0 <= xx < n and 0 <= yy < n:
                    A[i, xx * n + yy] = 1
    return A


def W_k(A, k):
    P = np.linalg.matrix_power(A, k - 1) if k > 1 else np.eye(len(A))
    return A * P


# ---------------------------------------------------------------- F1
head("F1  (L') 修好：1D 链上半径 k 外影响精确为 0")
Nn = 64
mid = Nn // 2
for k in (2, 4, 8, 16):
    far = Nn - 2                                  # 最远端边
    r0 = W_k(chain(Nn), k)[mid, mid + 1] / W_k(chain(Nn), k)[mid + 1, mid + 2]
    r1 = W_k(chain(Nn, mu=2.0, d=far), k)[mid, mid + 1] / W_k(chain(Nn, mu=2.0, d=far), k)[mid + 1, mid + 2]
    check("k=%-3d 远端扰动对中部权比的影响 = 0" % k, r1 == r0, "Δ=%.3e" % abs(r1 - r0))
# 近端扰动（距中部 1）非零 —— 反向控制
d_near = mid - 1
r0 = W_k(chain(Nn), 8)[mid, mid + 1] / W_k(chain(Nn), 8)[mid + 1, mid + 2]
r1 = W_k(chain(Nn, mu=2.0, d=d_near), 8)[mid, mid + 1] / W_k(chain(Nn, mu=2.0, d=d_near), 8)[mid + 1, mid + 2]
check("反向控制：近端扰动确实有影响（非零）", abs(r1 - r0) > 1e-6, "Δ=%.3e" % abs(r1 - r0))

# ---------------------------------------------------------------- F2
head("F2  2D 方格上同样精确为 0")
n2 = 24
mid2 = (n2 // 2) * n2 + (n2 // 2)
far2 = (n2 - 2) * n2 + (n2 - 3)
for k in (2, 4, 8, 16):
    A0 = grid(n2)
    A1 = grid(n2)
    A1[far2, far2 + 1] = A1[far2 + 1, far2] = 3.0
    W0, W1 = W_k(A0, k), W_k(A1, k)
    r0 = W0[mid2, mid2 + 1] / W0[mid2 + 1, mid2 + 2]
    r1 = W1[mid2, mid2 + 1] / W1[mid2 + 1, mid2 + 2]
    check("2D k=%-3d 远端扰动影响 = 0" % k, abs(r1 - r0) / abs(r0) == 0.0,
          "相对 %.3e" % (abs(r1 - r0) / abs(r0)))

# ---------------------------------------------------------------- F3
head("F3  (O) 保住：支撑只在最近邻")
A32 = chain(32)
for k in (2, 4, 6, 8):
    W = W_k(A32, k)
    nz = np.argwhere(np.abs(W) > 1e-12)
    far = max(abs(i - j) for i, j in nz)
    check("k=%-3d 支撑最大 |i-j| = 1" % k, far == 1, "实测 %d" % far)

# ---------------------------------------------------------------- F4
head("F4  正则 Γ：有限 k 体内权重精确均匀（度规退化为平坦）")
g = grid(16)
for k in (2, 4, 8):
    W = W_k(g, k)
    vals = np.array([W[x * 16 + y, x * 16 + y + 1] for x in range(4, 12) for y in range(4, 12)])
    rel = (vals.max() - vals.min()) / vals.mean()
    check("2D 方格 k=%-3d 体内相对差 = 0" % k, rel == 0.0, "%.3e" % rel)
c = chain(32)
for k in (2, 4, 8, 16):
    W = W_k(c, k)
    vals = np.array([W[i, i + 1] for i in range(8, 24)])
    check("1D 链   k=%-3d 体内相对差 = 0" % k, (vals.max() - vals.min()) == 0.0)

# ---------------------------------------------------------------- F5
head("F5  非正则 Γ：有限 k 体内权重非均匀（几何由非正则性承载）")
r = np.random.default_rng(7)
Ad = np.zeros((64, 64))
for i in range(63):
    v = 1.0 + r.random() * 0.35
    Ad[i, i + 1] = Ad[i + 1, i] = v
prevs = []
for k in (2, 4, 8):
    W = W_k(Ad, k)
    vals = np.array([W[i, i + 1] for i in range(10, 54)])
    rel = (vals.max() - vals.min()) / vals.mean()
    prevs.append(rel)
    check("无序链 k=%-3d 体内相对差 > 0.5（非均匀）" % k, rel > 0.5, "%.3f" % rel)
check("非均匀度随 k 增大（几何信息增强）",
      all(prevs[i] < prevs[i + 1] for i in range(len(prevs) - 1)),
      "%s" % [round(x, 3) for x in prevs])

# ---------------------------------------------------------------- F6
head("F6  Z4 的域基态剖面是 Perron 假象")
A = chain(64)
ev, V = np.linalg.eigh(A)
phi = V[:, -1]
if phi.sum() < 0:
    phi = -phi
Wp = A * np.outer(phi, phi)
vals = np.array([Wp[i, i + 1] for i in range(12, 52)])
rel_p = (vals.max() - vals.min()) / vals.mean()
check("Perron：体内非均匀（相对差 > 0.3）", rel_p > 0.3, "%.3f" % rel_p)
W4 = W_k(A, 4)
vals4 = np.array([W4[i, i + 1] for i in range(12, 52)])
check("有限 k=4：体内精确均匀 ⇒ 剖面消失", (vals4.max() - vals4.min()) == 0.0)

# ---------------------------------------------------------------- F7
head("F7  文档结论与三笔代价在位")
check("写明 (L') 修好、半径 k 外精确为 0", "(L')" in DOC and "精确为 0" in DOC)
check("写明 (O) 保住（支撑仍在最近邻）", "(O)" in DOC and "支撑" in DOC)
check("三笔代价表在位（半径 k／长度标度／Γ 非正则）",
      "局部半径 $k$" in DOC and "长度标度" in DOC and "Γ 必须非正则" in DOC)
check("指明前两笔不是新账（L／I2b／G44）", "I2b" in DOC and "G44" in DOC and "不是新账" in DOC)
check("红标：G40 §3 的'一致性检查'在 d≠2 下是 G2 的障碍",
      "红标" in DOC and "G40" in DOC and "G2" in DOC)
check("写明 G41 的判决只在 Perron 取法下成立", "只在 Perron" in DOC or "Perron（$k\\to\\infty$）取法下成立" in DOC)
check("给出下一步三条（I5 升级／随机细化 GH 极限／标度律复核）",
      "I5" in DOC and "Benjamini" in DOC and "标度律" in DOC)
check("诚实边界写明 \(L'\)] 未作一般证明", "未作一般证明" in DOC)
check("写明 k~L 是建议而非定理", "建议" in DOC and "不把它升为公理" in DOC)

# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
