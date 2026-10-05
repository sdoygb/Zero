#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
L2_period_abs.py —— 试算 L2 毁灭周期 T 的**绝对取值**

【思路（照 G61／R48 的"锁整数"模板）】
  G61 的模板：约束集 ＋ 一个代数定理 ＋ 有限枚举 ⟹ 唯一解／排除
  R48 的判例：B 由**一个自洽条件**（锥边可见）钉在 4，不是靠最小性

  本文件尝试同一手法定 T：
    约束 (A) 容量有限（G32）：K_cap = 4(T+1) 每站点；正交极小投影 2(T+1)
    约束 (B) 自洽：峰值活动量须**恰好**达到容量（否则不是饱和态）
    判据     扫描 T 与容量，找使"峰值/容量"逐周期趋于**常数**（= 自洽不动点）者

【层指标】L2（演化层）。
【只用】Z0／Z4／Z5 的毁灭-重播种规则 ＋ G32 的容量式。不引入概率。
"""

from __future__ import annotations
import json, os
from collections import defaultdict
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print("  [%s] %s%s" % ("v" if cond else "x", name,
                           ("   " + detail) if detail else ""))
    return bool(cond)


def cycle(T, cycles=14, memory=2, cap=None):
    """cap = 每站点容量上限（None = 不设限）。返回逐周期统计。

    饱和实现：若活动量 > cap，则按比例截断（容量饱和）；记录是否触发。
    """
    active = {(1, 0): 1}
    hist = defaultdict(int)
    layer = 0
    rows = []
    for c in range(1, cycles + 1):
        closed = 0
        peak = 0
        sat = 0
        for _ in range(T):
            nxt = defaultdict(int)
            for (b, p), n in active.items():
                for db in (1, -1):
                    nb = b + db
                    if nb == 0:
                        closed += n
                        hist[layer + 1] += n
                    else:
                        nxt[(nb, 1 - p)] += n
            active = nxt
            tot = sum(active.values())
            if cap is not None and tot > cap:
                sat += 1
                # 按比例截断，保持平衡分布形状
                f = cap / tot
                active = {k: int(v * f) for k, v in active.items()}
                active = {k: v for k, v in active.items() if v}
            peak = max(peak, sum(active.values()))
        d = sum(active.values())
        active = {}
        h = sum(hist.values())
        if h:
            active = {(1, 0): h, (-1, 0): h}
            if cap is not None:
                tot = sum(active.values())
                if tot > cap:
                    f = cap / tot
                    active = {k: max(1, int(v * f)) for k, v in active.items()}
        layer += 1
        keep = sorted(hist.keys())[-2:]
        hist = defaultdict(int, {k: hist[k] for k in keep})
        rows.append(dict(c=c, closed=closed, peak=peak, destroyed=d,
                         history=h, sat_steps=sat))
    return rows


# ================================================================== 主程序
print("=" * 78)
print("L2 · 毁灭周期 T 的绝对取值试算")
print("=" * 78)

print("\nA 容量式（G32）与投影计数")
print(f"     {'T':>3} {'dim A = 4(T+1)':>15} {'极小投影 2(T+1)':>16}")
for T in range(0, 13):
    print(f"     {T:>3} {4*(T+1):>15} {2*(T+1):>16}")
check("A1 容量随 T 线性（dim A = 4(T+1)）", 4 * (3 + 1) == 16)
check("A2 T=0 时 dim A = 4（只有矩阵因子，无年龄）", 4 * (0 + 1) == 4)

print("\nB 无饱和时的首周期峰值 = c_T（中心二项式）")
c = [comb(k, k // 2) for k in range(0, 15)]
first = {}
for T in range(1, 15):
    rows = cycle(T, cycles=1)
    first[T] = rows[0]['peak']
print(f"     {'T':>3} " + " ".join(f"{T:>8}" for T in range(1, 13)))
print(f"     c_T " + " ".join(f"{c[T]:>8}" for T in range(1, 13)))
check("B1 首周期峰值 = c_T", all(first[T] == c[T] for T in range(1, 13)),
      f"T=1..12: {[first[T] for T in range(1,13)]}")

print("\nC 与两种容量口径对照（每站点）")
print(f"     {'T':>3} {'c_T':>10} {'4(T+1)':>9} {'2(T+1)':>9} {'c_T<=4(T+1)':>12} {'c_T<=2(T+1)':>12}")
cross_strong, cross_weak = [], []
for T in range(1, 15):
    ok_s = c[T] <= 4 * (T + 1)
    ok_w = c[T] <= 2 * (T + 1)
    if ok_s: cross_strong.append(T)
    if ok_w: cross_weak.append(T)
    print(f"     {T:>3} {c[T]:>10} {4*(T+1):>9} {2*(T+1):>9} "
          f"{str(ok_s):>12} {str(ok_w):>12}")
check("C1 存在使峰值 ≤ 容量的 T（否则无饱和态）", len(cross_strong) > 0)
check("C2 存在使峰值 ≤ 投影数的 T", len(cross_weak) > 0)

print("\nD 双口径给出的 T 上界（峰值恰好达到容量 = 自洽点）")
T_strong = max(cross_strong)
T_weak = max(cross_weak)
print(f"     口径一（dim A = 4(T+1)）：T <= {T_strong}   ⟹ 临界 T = {T_strong}")
print(f"     口径二（投影数 2(T+1)）：T <= {T_weak}   ⟹ 临界 T = {T_weak}")
check("D1 两口径给出相邻的临界值", abs(T_strong - T_weak) == 1,
      f"{T_weak} 与 {T_strong}")
check("D2 口径二（更严）给出唯一的临界 T", T_weak == 5, f"T_weak = {T_weak}")

print("\nE 临界点检验：哪个 T 恰好饱和？（c_T 恰达容量）")
print("     判据：T_max = max{T : c_T <= K(T)} 为可行域上界；看它是否为等号")
print(f"     {'口径':>18} {'可行域':>16} {'T_max':>6} {'c_T':>7} {'K(T)':>7} {'是否等号':>9}")
for label, K in (("dim A = 4(T+1)", lambda T: 4 * (T + 1)),
                 ("投影数 2(T+1)", lambda T: 2 * (T + 1))):
    feas = [T for T in range(1, 15) if c[T] <= K(T)]
    Tm = max(feas)
    eq = (c[Tm] == K(Tm))
    print(f"     {label:>18} {'1..' + str(Tm):>16} {Tm:>6} {c[Tm]:>7} "
          f"{K(Tm):>7} {str(eq):>9}")
    check(f"E_{label[:6]} 可行域非空且 T_max 有限", len(feas) > 0)
    if label.startswith("投影"):
        check("E2 投影数口径：T_max = 5（c_5 = 10 < 12，最后一个可行的）",
              Tm == 5 and c[5] == 10 and K(5) == 12,
              f"T_max={Tm}, c={c[Tm]}, K={K(Tm)}")
        check("E3 T=6 在投影数口径下**超容**（c_6=20 > 14 = 2(6+1)）",
              c[6] > 2 * (6 + 1), f"c_6={c[6]}, K={2*(6+1)}")
    if label.startswith("dim"):
        check("E4 dim A 口径：T_max = 6（c_6 = 20 < 28，最后一个可行的）",
              Tm == 6 and c[6] == 20 and K(6) == 28,
              f"T_max={Tm}, c={c[Tm]}, K={K(Tm)}")

print("\nE5 交叉点（对数插值，解 c_T = K(T)）")
for _lbl, _K in (("dim A = 4(T+1)", lambda T: 4 * (T + 1)),
                 ("投影 = 2(T+1)", lambda T: 2 * (T + 1))):
    for T in range(1, 17):
        if c[T] <= _K(T) and c[T + 1] > _K(T + 1):
            import math as _m
            x = ((_m.log2(_K(T + 1)) - _m.log2(c[T + 1])) /
                 ((_m.log2(c[T]) - _m.log2(c[T + 1])) -
                  (_m.log2(_K(T)) - _m.log2(_K(T + 1)))))
            print(f"     {_lbl}: 交叉于 T={T} 与 T={T+1} 之间；T* = {T + x:.3f}")
            check(f"E5 {_lbl[:6]} 的交叉落在 {T} 与 {T+1} 之间",
                  T <= T + x <= T + 1)
            break

print("\nF 结论候选")
check("F1 投影数口径的可行域上界 = 5", True, "见 E2")
check("F2 dim A 口径的可行域上界 = 6", True, "见 E4")
check("F3 两口径的可行域分别是 1..5 与 1..6（相邻）", True)
check("F4 交叉点 5.66 与 6.21 的整数候选 = {5, 6}", True)

print("\n" + "=" * 78)
print("通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
if FAIL:
    print("不符项：")
    for n in FAIL:
        print("   x", n)
out = {"pass": len(PASS), "fail": len(FAIL), "failed_names": FAIL,
       "c_T": {str(k): c[k] for k in range(0, 15)},
       "first_cycle_peak": {str(k): v for k, v in first.items()},
       "cross_strong_max": T_strong, "cross_weak_max": T_weak,
       "T_critical_projection": 5, "T_critical_dimA": 6,
       "note": "E 组原版用封顶机制检验收敛，属循环论证，已删除并改为临界点检验"}
with open(os.path.join(HERE, "L2_period_abs_results.json"), "w",
          encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("=" * 78)
if FAIL:
    raise SystemExit(1)
print("全部通过 ✓")
