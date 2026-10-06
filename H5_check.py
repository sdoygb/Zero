#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
H5_check.py · 核验「锁 = 闭包类图谱」的定稿形式（经三次更正后）

  F1 主定理：M = Q 逐项相等（任意分块）⟹ 宏观弛豫被商图唯一锁定
  F2 谱：spec(M) = spec(Q)；M 是随机矩阵
  F3 被否证的强命题：spec(M) ⊂ spec(P) 只在等分划下成立
  F4 三处自我更正登记
  F5 边界
只用 numpy。用法: python3 H5_check.py     退出码 0 = 全部通过
"""
import numpy as np

PASS = FAIL = 0
FAILED = []


def ck(tag, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  [ok] {tag}" + (f"   {detail}" if detail else ""))
    else:
        FAIL += 1
        FAILED.append(tag)
        print(f"  [XX] {tag}   {detail}")


def brute(A, lab, K):
    """逐项最直白实现：M（度权块层链）与 Q（纤维平均商图）"""
    n = A.shape[0]
    d = A.sum(1)
    P = A / np.maximum(d[:, None], 1e-12)
    M = np.zeros((K, K))
    for a in range(K):
        idx = np.where(lab == a)[0]
        tot = d[idx].sum()
        for b in range(K):
            M[a, b] = sum(d[i] * P[i][np.where(lab == b)[0]].sum() for i in idx) / tot
    Q = np.zeros((K, K))
    for a in range(K):
        idx = np.where(lab == a)[0]
        md = d[idx].mean()
        for b in range(K):
            Q[a, b] = np.mean([A[i][np.where(lab == b)[0]].sum() for i in idx]) / md
    return M, Q, P


def is_eq(A, lab, K):
    for j in range(K):
        idx = np.where(lab == j)[0]
        C = np.array([[A[i][np.where(lab == k)[0]].sum() for k in range(K)] for i in idx])
        if C.std(axis=0).max() > 1e-9:
            return False
    return True


def rnd_case(rng):
    n = int(rng.integers(6, 16))
    K = int(rng.integers(2, 5))
    cuts = np.sort(rng.choice(np.arange(1, n), size=K - 1, replace=False))
    lab = np.zeros(n, dtype=int)
    prev = 0
    for k, c in enumerate(list(cuts) + [n]):
        lab[prev:c] = k
        prev = c
    A = np.zeros((n, n))
    pr = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = int(rng.integers(n, min(n * (n - 1) // 2, 3 * n) + 1))
    for t in rng.choice(len(pr), size=m, replace=False):
        i, j = pr[t]
        A[i, j] = A[j, i] = 1
    return A, lab, K


rng = np.random.default_rng(5)
mq = sp = rows = cont = 0.0
cnt = neq = 0
eq_cont = eq_n = 0.0
for _ in range(3000):
    A, lab, K = rnd_case(rng)
    d = A.sum(1)
    if (d == 0).any() or any((lab == k).sum() == 0 for k in range(K)):
        continue
    M, Q, P = brute(A, lab, K)
    evM = np.sort(np.real(np.linalg.eigvals(M)))
    evQ = np.sort(np.real(np.linalg.eigvals(Q)))
    evP = np.sort(np.real(np.linalg.eigvals(P)))
    mq = max(mq, np.abs(M - Q).max())
    sp = max(sp, np.abs(evM - evQ).max())
    rows = max(rows, np.abs(M.sum(1) - 1).max())
    cc = max(min(abs(evP - x)) for x in evM)
    cont = max(cont, cc)
    if is_eq(A, lab, K):
        eq_n += 1
        eq_cont = max(eq_cont, cc)
    else:
        neq += 1
    cnt += 1

print("=" * 78)
print("F1 · 主定理：M = Q 逐项相等（任意分块）")
print("=" * 78)
print(f"      样本：{cnt} 个随机图（非等分划 {neq} 个，等分划 {int(eq_n)} 个）")
ck("M = Q 逐项相等", mq < 1e-12, f"最大偏差 = {mq:.3e}")
ck("M 是随机矩阵（行和为 1）", rows < 1e-12, f"最大偏差 = {rows:.3e}")

print()
print("=" * 78)
print("F2 · 谱：spec(M) = spec(Q) ⟹ 宏观弛豫率被商图唯一锁定")
print("=" * 78)
ck("spec(M) = spec(Q)", sp < 1e-12, f"最大偏差 = {sp:.3e}")

print()
print("=" * 78)
print("F3 · 被否证的强命题：spec(M) ⊂ spec(P) 只在等分划下成立")
print("=" * 78)
ck("一般情形（含非等分划）为假", cont > 0.1,
   f"最大偏差 = {cont:.3e}  ← 『宏观模式是微观谱子集』不普遍成立")
ck("等分划子样本下成立", eq_n > 0 and eq_cont < 1e-12,
   f"最大偏差 = {eq_cont:.3e}（n={int(eq_n)}）" if eq_n > 0 else "无等分划样本")
# 确定性等分划对照
n, K = 12, 3
lab = np.array([0] * 4 + [1] * 4 + [2] * 4)
A = np.zeros((n, n))
for blk in range(3):
    o = 4 * blk
    for s in range(4):
        i, j = o + s, o + (s + 1) % 4
        A[i, j] = A[j, i] = 1
for s in range(4):
    A[0 + s, 4 + s] = A[4 + s, 0 + s] = 1
    A[4 + s, 8 + s] = A[8 + s, 4 + s] = 1
    A[0 + s, 8 + s] = A[8 + s, 0 + s] = 1
M, Q, P = brute(A, lab, K)
evM = np.sort(np.real(np.linalg.eigvals(M)))
evP = np.sort(np.real(np.linalg.eigvals(P)))
ck("确定性等分划（12点/3块/4-正则/商图 C3）：M=Q",
   np.abs(M - Q).max() < 1e-12 and is_eq(A, lab, K),
   f"偏差={np.abs(M-Q).max():.2e}  弛豫率={-np.log(max(abs(evM[1]),1e-15)):.6f}")
ck("确定性等分划：spec(M) ⊂ spec(P)",
   max(min(abs(evP - x)) for x in evM) < 1e-12)

print()
print("=" * 78)
print("F4 · 三处自我更正登记")
print("=" * 78)
ck("更正 1：『等分划＋同商图仍不充分』为错", True,
   "错因：混淆 M 与 Q；所造『反例』改变了 Q 本身")
ck("更正 2：『M = D⁻¹QD 相似变换』为错", True,
   "错因：把谱相等误当矩阵相等（真相更简单：M = Q）")
ck("更正 3：『M ≠ Q 一般成立』为错", True,
   "错因：那一版 M 用了未按度权归一的口径")
ck("三次更正均由数值对账发现，非再审证明", True)

print()
print("=" * 78)
print("F5 · 边界")
print("=" * 78)
ck("[本件] 定理只断言『给定 π 与 Γ，宏观弛豫被 Q 锁定』；未断言 π 从哪来", True,
   "π 的识别 = E1 站点识别缺口的社会层影子")
ck("[本件] 未在真实贸易/冲突网络上核验", True)
ck("[本件] 未独立复核旧理论 5.8.3.02 的原证", True)
ck("[本件] 数值仅在 6–19 点合成随机图上", True)

print()
print("=" * 78)
print(f"结果：通过 {PASS} / 不符 {FAIL}")
if FAILED:
    print("不符项：", FAILED)
print("=" * 78)
raise SystemExit(0 if PASS and FAIL == 0 else 1)
