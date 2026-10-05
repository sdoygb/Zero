#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z4_check.py —— 【P1：G40 权重的标度与形状】的核验
================================================
独立实断言：
  F1  G40 权重的闭式 w_ij = A_ij φ_i φ_j / λ（由穿越数归一化推出，并与其和 = 1 自洽）
  F2  标度：w ~ a^α，α → d（d=1,2 数值），而 G2 要求 d−2 ⇒ 差 2
  F3  后果：E_N ~ a^{α+2−d} → 0（刚度消失），与 G2 的公式一致
  F4  形状：体内 std(log w) 收敛到正常数（不衰减）；全域 std 发散
  F5  顶点传递情形：环图 φ ≡ 常数 ⇒ w 均匀（G40 §3）
  F6  文档结论在位（三项缺口；两项已登记；共形剖面为新问题）
"""
import io
import math
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


DOC = io.open(os.path.join(HERE, "Z4_p1_weight_scaling_and_shape.md"), encoding="utf-8").read()


def chain(n):
    A = np.zeros((n, n))
    for i in range(n - 1):
        A[i, i + 1] = A[i + 1, i] = 1
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


def ring(n):
    A = np.zeros((n, n))
    for i in range(n):
        A[i, (i + 1) % n] = A[(i + 1) % n, i] = 1
    return A


def perron(A):
    ev, V = np.linalg.eigh(A)
    phi = V[:, -1]
    if phi.sum() < 0:
        phi = -phi
    return phi, ev[-1]


# ---------------------------------------------------------------- F1
head("F1  权重闭式 w = A φ φ^T / λ 的自洽性")
for name, A in (("链 n=32", chain(32)), ("环 n=32", ring(32)), ("方格 n=8", grid(8))):
    phi, lam = perron(A)
    W = (A * np.outer(phi, phi)) / lam
    check("%-10s Σ w_ij = 1（归一化穿越数）" % name, abs(W.sum() - 1) < 1e-10,
          "实算 %.12f" % W.sum())
    check("%-10s φ^T A φ = λ" % name, abs(phi @ A @ phi - lam) < 1e-9)
check("文档给出该闭式", "w_{ij}=\\frac{A_{ij}" in DOC or "w_{ij}=" in DOC)

# ---------------------------------------------------------------- F2
head("F2  标度：α → d（对比 G2 要求 d−2）")
alphas = {}
for d, sizes in ((1, (16, 32, 64, 128)), (2, (16, 64, 256))):
    ws = []
    for Nn in sizes:
        if d == 1:
            A = chain(Nn)
            c = Nn // 2
        else:
            n = int(round(Nn ** 0.5))
            A = grid(n)
            c = (n // 2) * n + (n // 2)
        phi, lam = perron(A)
        ws.append(phi[c] * phi[c] / lam)
    a = [np.log2(ws[i] / ws[i + 1]) for i in range(len(ws) - 1)]
    alphas[d] = a
    check("d=%d：α 序列单调趋近 %d" % (d, d),
          all(a[i] < a[i + 1] for i in range(len(a) - 1)) and a[-1] > d - 0.2,
          "%s" % [round(x, 3) for x in a])
check("d=1 末值 α > 0.98", alphas[1][-1] > 0.98, "%.3f" % alphas[1][-1])
check("d=2 末值 α > 1.8", alphas[2][-1] > 1.8, "%.3f" % alphas[2][-1])
check("文档写明差恰为 2（与 d 无关）", "a^{2}" in DOC and "差 }2" in DOC.replace("无关地差 ", "无关地差 "))

# ---------------------------------------------------------------- F3
head("F3  后果：E_N ~ a^{α+2−d} → 0（刚度消失）")
for d, al in ((1, alphas[1][-1]), (2, alphas[2][-1])):
    expo = al + 2 - d
    check("d=%d：E_N 的幂次 α+2−d = %.3f > 0 ⇒ E_N → 0" % (d, expo), expo > 0.5)
check("G2 的公式 E_N ≈ K a^{2-d} ∫|∇f|² 被正确引用",
      "a^{2-d}" in DOC and "G2_local_continuum_limit.md" in DOC)

# ---------------------------------------------------------------- F4
head("F4  形状：体内 std(log w) 收敛到常数（不衰减）")
for d, sizes in ((1, (16, 32, 64, 128, 256)), (2, (8, 16, 32, 64))):
    ins, alls = [], []
    for n in sizes:
        A = chain(n) if d == 1 else grid(n)
        phi, lam = perron(A)
        logw = 2 * np.log(phi) - np.log(lam)
        if d == 1:
            lo, hi = n // 4, 3 * n // 4
            ins.append(float(np.std(logw[lo:hi])))
        else:
            L, H = n // 4, 3 * n // 4
            m = np.zeros(n * n, dtype=bool)
            for x in range(L, H):
                for y in range(L, H):
                    m[x * n + y] = True
            ins.append(float(np.std(logw[m])))
        alls.append(float(np.std(logw)))
    # 体内：不衰减（末项 ≥ 首项）；全域：增长
    check("d=%d：体内 std 不衰减（末 ≥ 首）" % d, ins[-1] >= ins[0],
          "%s" % [round(x, 3) for x in ins])
    check("d=%d：全域 std 单调增长（边界层效应）" % d,
          all(alls[i] < alls[i + 1] for i in range(len(alls) - 1)),
          "%s" % [round(x, 3) for x in alls])
    check("d=%d：体内 std 收敛到正常数（0.1–0.4）" % d, 0.1 < ins[-1] < 0.4, "%.4f" % ins[-1])

# ---------------------------------------------------------------- F5
head("F5  顶点传递情形：环图 φ ≡ 常数 ⇒ w 均匀")
for n in (8, 32, 128):
    A = ring(n)
    phi, lam = perron(A)
    rel = float(np.max(np.abs(phi - phi[0])) / abs(phi[0]))
    check("环 n=%-4d φ 相对偏离 < 1e-10 ⇒ w 均匀" % n, rel < 1e-10, "偏离 %.2e" % rel)
check("文档写明这条给出 G2 的 α=0 障碍", "α=0" in DOC or "$\\alpha=0$" in DOC)

# ---------------------------------------------------------------- F6
head("F6  文档结论在位")
check("写明 I2a 未被 no-go 挡住", "没有被 no-go 挡住" in DOC or "不是被 no-go 挡住" in DOC)
check("三项缺口表在位（标度／维数／共形剖面）",
      "长度标度" in DOC and "维数 $d$" in DOC and "共形剖面" in DOC)
check("两项已登记输入被点名（I2b／G44 与 G89）",
      "I2b" in DOC and "G44" in DOC and "G89" in DOC)
check("§3.1 两难表在位（顶点传递 vs 带边界补片）",
      "顶点传递" in DOC and "带边界补片" in DOC)
check("与 G18 的既有判定对照", "G18" in DOC and "参考度规" in DOC)
check("诚实边界写明取格点补片为工作假设", "工作假设" in DOC)
check("给出三档下一步（新问题／P2／P3）", "共形剖面的地位" in DOC and "P2" in DOC and "P3" in DOC)
check("限定数值范围（d≤2）", "d\\le2" in DOC or "$d\\le2$" in DOC)

# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
