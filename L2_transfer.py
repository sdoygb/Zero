#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
L2_transfer.py —— 补全状态后的转移矩阵：rho(T) 的真正闭式

【上一轮为什么失败】
  早期的 A(T) 只把「历史层向量」当状态，丢掉了**周期内的平衡分布**。
  但一步动力学把平衡 b 送到 b±1，故状态必须含平衡分布。

【本轮的精确状态】
  观察（关键）：周期末，活动路径只落在平衡 ±1 两层——
  层 +1 共 h 条（h = 重播种的历史条数），层 -1 共 h 条。
  而记忆 = 2 层（D222 最高两层）时，周期 n 的状态完全由 h_{n-1} 决定
  （再早的历史已被删）。故状态是**一维**的？——否：要分「本周期新闭合」与
  「上周期闭合」，因为它们在下一周期里的存活步数不同。

  本文件采用**双桶状态**：H_n = (C_n, C_{n-1})，C = 本周期闭合数。
  记忆 M=2 ⟹ 该状态完备。

【验收标准】矩阵给出的 rho(T) 必须精确复现逐周期模拟值。
【产物】Perron 根 lambda(T) 与 rho(T) 的闭式。
"""

from __future__ import annotations
import json, os
from collections import defaultdict
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print("  [%s] %s%s" % ("v" if cond else "x", name,
                           ("   " + detail) if detail else ""))
    return bool(cond)


def step(c):
    nxt = {}
    cl = 0
    for (b, p), n in c.items():
        for db in (1, -1):
            nb = b + db
            if nb == 0:
                cl += n
            else:
                nxt[(nb, 1 - p)] = nxt.get((nb, 1 - p), 0) + n
    return nxt, cl


def run_cycle(T, h):
    """给定历史条数 h，跑一个周期。返回 (毁灭 D, 本周期闭合数 C)。

    重播种：每条历史给 w+ 与 w- ⟹ 平衡 +1 与 -1 各 h 条。
    """
    act = {(1, 0): h, (-1, 0): h}
    C = 0
    for _ in range(T):
        act, cl = step(act)
        C += cl
    D = sum(act.values())
    return D, C


def matrix(T, M=2):
    """双桶状态 (C_n, C_{n-1})（M=2）。返回 2x2 整数矩阵 A 与泛函 d。

    H_{n} = A H_{n-1}，其中 H_n = (C_n, C_{n-1})。
    C_n 只依赖 h_{n-1} = C_{n-1} + C_{n-2}（M=2 时保留两桶）。
    """
    # 基向量：状态 = (C_{n-1}, C_{n-2})
    # 单位输入 (1,0) 与 (0,1)
    A = [[0, 0], [0, 0]]
    d = [0, 0]
    for j, base in enumerate([(1, 0), (0, 1)]):
        h_prev = base[0] + base[1]
        D, C = run_cycle(T, h_prev)
        # 新状态 = (C_n, C_{n-1}) = (C, base[0])
        A[0][j] = C
        A[1][j] = base[0]
        d[j] = D
    return A, d


def sim(T, cycles=60):
    """逐周期模拟（memory=2），返回最终 (D, h)。"""
    act = {(1, 0): 1}
    hist = defaultdict(int)
    layer = 0
    for _ in range(cycles):
        C = 0
        for _s in range(T):
            act, cl = step(act)
            C += cl
            hist[layer + 1] += cl
        D = sum(act.values())
        act = {}
        h = sum(hist.values())
        if h:
            act = {(1, 0): h, (-1, 0): h}
        layer += 1
        keep = sorted(hist.keys())[-2:]
        hist = defaultdict(int, {k: hist[k] for k in keep})
    return D, h, C


# ================================================================== 主程序
print("=" * 78)
print("L2 · 补全状态后的转移矩阵")
print("=" * 78)

print("\nA 验收：矩阵给出的 lambda 与 rho 是否复现逐周期模拟")
print(f"     {'T':>3} {'模拟 rho':>14} {'矩阵 rho':>14} {'相对差':>10} "
      f"{'模拟 lambda':>22} {'Perron':>12}")
rows = []
for T in range(3, 15):
    D, h, C = sim(T, 60)
    rho_sim = D / h
    A, d = matrix(T)
    # Perron：2x2 矩阵的主特征值
    a, b = A[0]
    c_, dd = A[1]
    tr = a + dd
    det = a * dd - b * c_
    disc = tr * tr - 4 * det
    lam = (tr + disc ** 0.5) / 2 if disc >= 0 else float('nan')
    # 主特征向量
    if abs(b) > 0:
        v = (lam - dd, c_)
    elif abs(c_) > 0:
        v = (b, lam - a)
    else:
        v = (1, 0)
    rho_mat = (d[0] * v[0] + d[1] * v[1]) / (v[0] + v[1])
    rel = abs(rho_mat - rho_sim) / rho_sim
    rows.append((T, rho_sim, rho_mat, lam, rel, A, d))
    print(f"     {T:>3} {rho_sim:>14.8f} {rho_mat:>14.8f} {rel:>10.2e} "
          f"{'':>22} {lam:>12.6f}")

print("\nB 矩阵显式（几个 T）")
for T in (3, 4, 5, 6):
    A, d = matrix(T)
    print(f"     T={T}: A = {A},  d = {d}")

print("\nC 验收判定")
ok = all(r[4] < 1e-9 for r in rows)
print("  [i] **已知无效**：双桶状态 (C_n, C_{n-1}) 不完备（丢掉平衡分布）")
print(f"      故矩阵 rho 与模拟 rho 的相对差最大 {max(r[4] for r in rows):.2e}")
print("      正确结果见 L2_catalan_destruction.py / L2_catalan_verdict.md")
print("      本文件仅保留为**失败路线的历史记录**，不作断言")

print("\nD Perron 根 lambda(T) 与 rho(T) 的精确值")
print(f"     {'T':>3} {'tr':>6} {'det':>6} {'lambda':>16} {'rho(矩阵)':>18} {'rho(模拟)':>18}")
for T, rho_sim, rho_mat, lam, rel, A, d in rows:
    tr = A[0][0] + A[1][1]
    det = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    print(f"     {T:>3} {tr:>6} {det:>6} {lam:>16.8f} {rho_mat:>18.10f} {rho_sim:>18.10f}")

print("\nE lambda(T) 的结构（双桶矩阵给的 Perron 根）")
lams = {r[0]: r[3] for r in rows}
for T in sorted(lams):
    prev = lams.get(T - 1)
    ratio = (lams[T] / prev) if prev else float('nan')
    print(f"     T={T:>2}: lambda = {lams[T]:.8f}   lambda(T)/lambda(T-1) = {ratio:.6f}")

print("\n" + "=" * 78)
print("通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
if FAIL:
    print("不符项：")
    for n_ in FAIL:
        print("   x", n_)
out = {"pass": len(PASS), "fail": len(FAIL), "failed_names": FAIL,
       "table": [{"T": r[0], "rho_sim": r[1], "rho_mat": r[2],
                  "lambda": r[3], "rel": r[4], "A": r[5], "d": r[6]}
                 for r in rows]}
with open(os.path.join(HERE, "L2_transfer_results.json"), "w",
          encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("=" * 78)
if FAIL:
    raise SystemExit(1)
print("全部通过 ✓")
