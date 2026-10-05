#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E1_NN_check.py —— Nielsen–Ninomiya 四前提核查（独立探针；不改动任何既有结论）

【背景】
  SYNTHESIS_zero_to_standard_model.md §5.5 的 E1：
    Nielsen–Ninomiya 定理：局域 ＋ 厄米 ＋ 平移不变 ＋ 双线性
    ⟹ 手征费米子必然加倍（不能只有一个手征费米子）。
  我们的费米扇区候选来自：
    R14.2  补偿移动 ⟹ 最近邻 hopping，且算符 = G40 的图 Laplacian；
    R14.4/R15.1  项链 L-循环 ⟹ 双覆盖 −1 ⟹ 反对称（费米）扇区。
  必须逐条核查这四个前提。

【本脚本只使用】
  Z1 定理 1 的补偿移动 ( +,− ) ↔ ( −,+ )、Z0③ 的全分支/整数重数、
  R14.2 的算子识别。不引用任何 U 系材料，不引入概率。

【方法学更正（第一版探针的两处错误，本版已改）】
  错误 A：把「环图平移 w ↦ w[-1:]+w[:-1]」当成了内容图 Γ_L 的对称。
          实际它 **不** 保度（边界位置少一个邻位），故不是对称。
  错误 B：只试了「所有位置同时移位」这一种候选，未穷举 {0,…,L−1} 个移位，
          因而把「无平移不变」误报成「平移不变被破坏 ⟹ 前提不成立」，
          却没排除「某个非零移位仍保图」。本版穷举全部 L 个移位。
  错误 C：色散对照用了固定 k 网格，长度不匹配时越界。本版改为按谱序对齐。
"""

import itertools
import numpy as np
from collections import Counter

TOL = 1e-9
results = []      # 断言（参与通过/失败计数）
findings = []     # 发现项（预期的负结果，不参与计数）


def check(name, ok, detail=""):
    """断言：应当成立。"""
    results.append((name, bool(ok), detail))
    print(f"  {'✓' if ok else '✗'} {name}" + (f"   [{detail}]" if detail else ""))
    return bool(ok)


def finding(name, detail=""):
    """发现项：本探针要报告的负结果，不是失败。"""
    findings.append((name, detail))
    print(f"  ◆ {name}" + (f"   [{detail}]" if detail else ""))


# ------------------------------------------------------------------ 图
def words(L):
    """零和词：长度 L，n_+ = n_- = L/2（Z1 定理 2）。"""
    if L % 2:
        raise ValueError("零和 ⟹ L 偶")
    return [tuple(+1 if i in pos else -1 for i in range(L))
            for pos in itertools.combinations(range(L), L // 2)]


def moves(w):
    """Z1 定理 1 的补偿移动：交换一对相邻异号步。"""
    return [tuple(w[:i] + (w[i + 1], w[i]) + w[i + 2:])
            for i in range(len(w) - 1) if w[i] != w[i + 1]]


def build(L):
    W = words(L)
    idx = {w: i for i, w in enumerate(W)}
    n = len(W)
    rows, cols = [], []
    for w in W:
        i = idx[w]
        for v in moves(w):
            rows.append(i)
            cols.append(idx[v])
    A = np.zeros((n, n))
    A[rows, cols] = 1.0
    A = np.maximum(A, A.T)
    return W, idx, A


def perm_matrix(n, sig):
    P = np.zeros((n, n))
    P[sig, np.arange(n)] = 1.0
    return P


# ================================================================== 主程序
print("=" * 78)
print("E1 · Nielsen–Ninomiya 四前提核查（Z1 定理 1 的补偿移动）")
print("=" * 78)

L_LIST = [4, 6, 8, 10, 12]
summary = {}

for L in L_LIST:
    print(f"\n--- L = {L} ---")
    W, idx, A = build(L)
    n = len(W)
    Lap = np.diag(A.sum(axis=1)) - A
    ev = np.linalg.eigvalsh(Lap)
    n_zero = int(np.sum(np.abs(ev) < 1e-8))

    # ---------------- 前提① 局域 ----------------
    supps = {sum(1 for a, b in zip(w, v) if a != b)
             for w in W for v in moves(w)}
    ok1 = supps == {2}
    check(f"① 局域：移动支撑集合 = {sorted(supps)}（应为 {{2}}）", ok1,
          "R14.2 第 1 条")

    # 连通性（R14.2 第 2 条）
    check(f"   图连通：Laplacian 零本征值重数 = {n_zero}（应为 1）",
          n_zero == 1, "R14.2 第 2 条")

    # ---------------- 前提③ 厄米 ----------------
    herm = np.max(np.abs(Lap - Lap.conj().T))
    check(f"③ 厄米：‖L − L†‖∞ = {herm:.3e}", herm < TOL)

    # ---------------- 前提② 双线性 ----------------
    asym = np.max(np.abs(Lap - Lap.T))
    check(f"② 双线性：L 对称实，‖L − Lᵀ‖∞ = {asym:.3e} ⟹ 二次型良定义",
          asym < TOL and herm < TOL)

    # ---------------- 前提④ 平移不变（穷举全部 L 个移位）----------------
    keeper = []
    for r in range(L):
        sig = [idx[w[r:] + w[:r]] for w in W]
        if np.array_equal(perm_matrix(n, sig) @ A @ perm_matrix(n, sig).T, A):
            keeper.append(r)
    ok4 = len(keeper) > 1          # 恒等之外还有别的移位
    finding(f"④ 平移不变：保持 Γ_L 的移位 = {keeper} ⟹ "
            f"{'成立' if ok4 else '【不成立】NN 前提被打破'}",
            "穷举全部 L 个移位")

    # 反射 w ↦ reverse(w)（这是真实的对称）
    rev = [idx[w[::-1]] for w in W]
    ref_ok = np.array_equal(perm_matrix(n, rev) @ A @ perm_matrix(n, rev).T, A)
    finding(f"（对照）反射 w ↦ reverse(w) 保持图：{ref_ok}",
            "我们的真实对称是反射，不是平移")

    # ---------------- 度分布 ----------------
    deg = A.sum(axis=1).astype(int)
    dd = dict(sorted(Counter(deg.tolist()).items()))
    print(f"     度分布 {dd}")
    check("   度为位置依赖（非齐次）⟹ 无平移不变性的直接证据",
          len(dd) > 1)

    # ---------------- 结论侧：加倍 ----------------
    finding(f"加倍检查：无能隙模数 = {n_zero}", "NN 结论侧：>1 ⟹ 加倍")
    spec = ev
    sym_err = np.max(np.abs(spec + spec[::-1]))
    finding(f"手征谱对称性 max|λ + λ_rev| = {sym_err:.3e}",
            "手征算子是否存在（我们尚无 γ）")

    summary[L] = dict(n=n, n_zero=n_zero, ok1=ok1, ok4=ok4,
                      keeper=keeper, ev=ev, deg=dd)

# ================================================================== 汇总
print("\n" + "=" * 78)
print("汇总")
print("=" * 78)
print(f"{'L':>4} {'态数':>7} {'零模数':>7} {'①局域':>7} {'②双线性':>8} "
      f"{'③厄米':>7} {'④平移不变':>10} {'保图移位':>18}")
for L in L_LIST:
    s = summary[L]
    print(f"{L:>4} {s['n']:>7} {s['n_zero']:>7} "
          f"{'✓' if s['ok1'] else '✗':>7} {'✓':>8} {'✓':>7} "
          f"{'✓' if s['ok4'] else '✗':>10} {str(s['keeper']):>18}")

# ================================================================== 色散
print("\n" + "=" * 78)
print("低能色散（按谱序对齐；L = 10）")
print("=" * 78)
L = 10
ev = np.sort(summary[L]['ev'])
k = np.arange(len(ev)) * 2 * np.pi / L
pred = 2 * (1 - np.cos(k))
print(f"{'序':>3} {'λ（本征值）':>14} {'2(1−cos k)':>14} {'差':>13}")
for i in range(len(ev)):
    print(f"{i:>3} {ev[i]:>14.6f} {pred[i]:>14.6f} {ev[i] - pred[i]:>13.3e}")
gap = ev[1] - ev[0]
print(f"\n谱隙 λ_1 − λ_0 = {gap:.6f}（>0 ⟹ 该单粒子谱有隙，非无能隙色散）")

# ================================================================== 判定
print("\n" + "=" * 78)
print("判定")
print("=" * 78)
all1 = all(summary[L]['ok1'] for L in L_LIST)
all4 = all(summary[L]['ok4'] for L in L_LIST)
allz = all(summary[L]['n_zero'] == 1 for L in L_LIST)
print(f"①局域（全部 L）      : {all1}")
print(f"④平移不变（全部 L）  : {all4}")
print(f"零模数恒为 1         : {allz}")
print()
if all1 and not all4:
    print("⟹ 结论：四前提中 **平移不变不成立**（Γ_L 的非平凡保图移位集为空），")
    print("   故 Nielsen–Ninomiya 定理**不适用**于该结构，")
    print("   它**不构成**「手征费米子必加倍」对本结构的否决。")
    print("   同时：该结论**也不**构成对费米扇区的正面支持——")
    print("   因为 NN 的适用条件是「有平移不变性可定义能带结构」，")
    print("   失去它意味着该算符**没有动量空间描述**，其低能极限需另行处理。")

n_pass = sum(1 for _, ok, _ in results if ok)
print(f"\n断言：{n_pass} / {len(results)} 通过")
print(f"发现项（预期的负结果）：{len(findings)} 条")
if n_pass != len(results):
    print("不符项：")
    for name, ok, _ in results:
        if not ok:
            print("   ✗", name)
    raise SystemExit(1)
print("断言全部通过 ✓")
print("（发现项不是失败：前提④为假正是本探针要报告的结论）")
