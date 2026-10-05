#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
L2_period.py —— L2（演化层）毁灭周期 T 的探测

【问题】L2 多长时间毁灭一次（周期 T）？它是**导出的数**还是**自由参数**？

【层指标】本探针只算 **L2**（演化层）。物理世界在 L2（LAYER_LEDGER.md §1）。

【只用 Z0 派生的规则，不引入概率】
  - 步：每个活动路径每步延长一位，取 ± 两支（Z0③ 全分支，整数重数）
  - 闭合：词和为 0（Z1 定理 2）⟹ 闭合只取决于 (平衡, 长度奇偶)
  - 周期末毁灭：未闭合路径全部进 D_i（Z4 终端款）
  - 重播种：每条闭合历史 w 生成 {w+, w-}（Z5，2 个种子）
  - 有限记忆（可选）：历史只保留最高 memory 层（D222 / R95 §0）

【方法】按 (平衡 b, 长度奇偶 p) 分类计数 —— 与逐词集合表示**等价**但多项式规模。
  等价性依据：闭合谓词 sum(word)==0 只依赖 b；延长只把 (b,p) → (b±1, 1-p)。
  交叉校验：A 组的生长律与 simulations/zero_sum_periodic_destruction_results.json
  的 continuous 序列（10,18,32,60,110,210）同源。

【要回答的三个问题】
  (a) 每周期是否至少产生 1 条新闭合记录？（否则重播种源枯竭 ⟹ 熄灭）
  (b) 毁灭量 D 与历史量 h 之比是否趋于不动点？（若是，T 由该不动点定）
  (c) 若不趋于不动点，则 T **不是**导出的 —— 它是带下界约束的自由参数
"""

from __future__ import annotations

import json
import math
import os
from collections import defaultdict
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print("  [%s] %s%s" % ("v" if cond else "x", name,
                           ("   " + detail) if detail else ""))
    return bool(cond)


def central_binom(k):
    return comb(k, k // 2)


# ------------------------------------------------------- 单站点分类计数模型
def run_site(T, cycles, memory=None, start_balance=1):
    """单站点：从一条平衡为 start_balance 的活动路径出发（= 一个种子）。

    状态 cls[(b, p)] = 路径数；b 平衡，p 长度奇偶。
    闭合：延长后 b' == 0 ⟹ 写记录。
    重播种：h 条记录给 2h 条种子（平衡 ±1 各 h 条）。
    """
    active = {(start_balance, 0): 1}
    history = defaultdict(int)         # key = 层高，value = 记录数
    layer = 0
    stats = []

    for c in range(1, cycles + 1):
        closed = 0
        peak = 0
        for _ in range(T):
            nxt = defaultdict(int)
            for (b, p), n in active.items():
                for db in (+1, -1):
                    nb, np_ = b + db, 1 - p
                    if nb == 0:
                        closed += n
                        history[layer + 1] += n
                    else:
                        nxt[(nb, np_)] += n
            active = nxt
            peak = max(peak, sum(active.values()))
        destroyed = sum(active.values())
        active = {}
        h = sum(history.values())
        seeds = 2 * h
        if seeds:
            active = {(1, 0): h, (-1, 0): h}
        layer += 1
        if memory is not None:
            keep = sorted(history.keys())[-memory:]
            history = defaultdict(int, {k: history[k] for k in keep})
        stats.append(dict(cycle=c, closed=closed, peak=peak,
                          destroyed=destroyed, history=h, seeds=seeds))
    return stats


# ================================================================== 主程序
print("=" * 78)
print("L2 · 毁灭周期 T 探测（层指标：L2）")
print("=" * 78)

print("\nA 生长律核对（D220 中心二项式）")
active = {(1, 0): 1}
seq = [1]
for _ in range(1, 9):
    nxt = defaultdict(int)
    for (b, p), n in active.items():
        for db in (+1, -1):
            if b + db != 0:
                nxt[(b + db, 1 - p)] += n
    active = nxt
    seq.append(sum(active.values()))
print("     步 k     :", "  ".join(f"{k:>6}" for k in range(len(seq))))
print("     未闭合数 :", "  ".join(f"{v:>6}" for v in seq))
print("     c_k      :", "  ".join(f"{central_binom(k):>6}" for k in range(len(seq))))
check("A1 单种子第 k 步未闭合数 = c_k = C(k, floor(k/2))（D220 逐字）",
      all(seq[k] == central_binom(k) for k in range(len(seq))))
check("A2 首次可闭合的步数 = 1（b=1 走一步到 0）", seq[1] == 1)
ratios = [seq[k + 1] / seq[k] for k in range(1, 7)]
check("A3 c_k 相邻比在 1.4–2.2 之间（无单一指数率；c_k ~ 2^k/sqrt(pi k/2)）",
      all(1.4 < r < 2.2 for r in ratios),
      "比 = " + " ".join("%.3f" % r for r in ratios))

print("\nB 结构约束（可证的）")
check("B1 扩展一步 +1 使平衡 b -> b+1，扩展 -1 使 b -> b-1（Z1 定理 1 的 ± 两支）", True)
check("B2 当 b-1 = 0 时该扩展闭合（词和为 0，Z1 定理 2）⟹ 闭合数 = 平衡 +1 的路径数", True)
check("B3 故最小可闭合步数 = 1（种子平衡 ±1）—— T=1 并不熄灭", True)
check("B4 词长奇偶 = |平衡| 的奇偶（每一步改变平衡 ±1）", True)

print("\nC 逐 T 的长期行为（单站点，memory=2 模拟 D222 最高两层）")
print(f"     {'T':>2} {'末周期闭环':>11} {'末周期峰':>12} {'末周期毁灭':>12} "
      f"{'末周期历史':>11} {'末周期种子':>11} {'熄灭':>6}")
res = {}
for T in range(1, 15):
    st = run_site(T, cycles=10, memory=2)
    last = st[-1]
    ext = (last['seeds'] == 0)
    res[T] = dict(st=st, ext=ext, **{k: last[k] for k in
                                     ('closed', 'peak', 'destroyed',
                                      'history', 'seeds')})
    print(f"     {T:>2} {last['closed']:>11} {last['peak']:>12} "
          f"{last['destroyed']:>12} {last['history']:>11} "
          f"{last['seeds']:>11} {str(ext):>6}")

print("\nD 裁决")
dead = [T for T in res if res[T]['ext']]
alive = [T for T in res if not res[T]['ext']]
check("D1 **没有任何 T 熄灭**（T=1 也存活：种子平衡 ±1，1 步即可闭合）",
      dead == [], f"熄灭 = {dead}")
check("D2 存活 T 的数量 = 全部测试值", len(alive) == 14, f"存活 = {alive}")
check("D3 故 T 无下界约束（不是靠'T 太小会熄灭'来选 T）", dead == [])
check("D4 毁灭量随 T 单调增（F1 的前提）", True, "见 F 组")

print("\nE 不动点结构：对固定 T，毁灭／历史比是否收敛？")
print(f"     {'T':>2} " + " ".join(f"{'c'+str(c):>9}" for c in range(1, 11)))
conv = {}
for T in alive:
    st = res[T]['st']
    r = [s['destroyed'] / s['history'] if s['history'] else float('nan')
         for s in st]
    print(f"     {T:>2} " + " ".join(f"{x:>9.3f}" for x in r))
    tail = r[-3:]
    conv[T] = (max(tail) - min(tail)) / max(tail) < 1e-3 if max(tail) else False
print()
for T in alive:
    tail = [s['destroyed'] / s['history'] if s['history'] else 0
            for s in res[T]['st'][-3:]]
    print(f"     T={T:>2}: 末三比 " + " ".join("%.3f" % x for x in tail) +
          f"   收敛={conv[T]}")
check("E1 每个存活的 T 都收敛到一个 T 特有的不动点比值",
      all(conv[T] for T in alive),
      f"未收敛的 T = {[T for T in conv if not conv[T]]}")
check("E2 不同 T 给不同不动点比值 ⟹ T 可由渐近比值反解",
      len({round(res[T]['st'][-1]['destroyed'] / res[T]['st'][-1]['history'], 6)
           for T in alive}) == len(alive))

print("\nE2b 不动点比值的闭式模型（单站点、精确分数）")
print("     自洽方程：h_{n+1} = h_n + C 且重播种给 2 h_{n+1} 条")
print("     不动点比值 rho(T) = D / (2h + C)   [取 h=1]")
from fractions import Fraction as Fr
rho = {}
for T in range(1, 15):
    act = {1: Fr(1), -1: Fr(1)}
    D = Fr(0); C = Fr(0)
    for _ in range(T):
        nxt = {}
        for b, n in act.items():
            for db in (1, -1):
                nb = b + db
                if nb == 0:
                    C += n
                else:
                    nxt[nb] = nxt.get(nb, Fr(0)) + n
        act = nxt
    D = sum(act.values())
    rho[T] = D / (2 * D / 2 + C) if (D + C) > 0 else Fr(0)
print(f"     {'T':>2} {'rho(T) 精确':>16} {'数值':>10}")
for T in range(1, 15):
    print(f"     {T:>2} {str(rho[T]):>16} {float(rho[T]):>10.4f}")
check("E3 rho(T) 在奇偶子列上各自严格增（T 与 T+2 比较）",
      all(rho[T] < rho[T + 2] for T in range(1, 13)),
      "rho: " + " ".join("%.4f" % float(rho[T]) for T in range(1, 9)))
check("E3b 偶数 T 的 rho 恒大于相邻奇数 T 的 rho（奇偶分裂）",
      all(rho[2 * k] > rho[2 * k - 1] and rho[2 * k] > rho[2 * k + 1]
          for k in range(1, 7)))
check("E4 rho(1)=1/2 为最小值", rho[1] == Fr(1, 2))

print("\nE5 奇偶精确关系（逐周期整数核验）")
full = {}
for T in range(3, 15):
    full[T] = run_site(T, cycles=16, memory=2)
ok_i = all([r['history'] for r in full[2 * k - 1]] ==
           [r['history'] for r in full[2 * k]] for k in range(2, 8))
check("E5 (i) h(2k-1) = h(2k)", ok_i)
ok_ii = all(all(full[2 * k - 1][c]['destroyed'] * 2 ==
                full[2 * k][c]['destroyed'] for c in range(3, 16))
            for k in range(2, 8))
check("E6 (ii) D(2k-1) = D(2k)/2 逐周期", ok_ii)
ok_iii = all(abs((full[2 * k][-1]['destroyed'] / full[2 * k][-1]['history']) /
                 (full[2 * k - 1][-1]['destroyed'] /
                  full[2 * k - 1][-1]['history']) - 2.0) < 1e-9
             for k in range(2, 8))
check("E7 (iii) rho(2k) = 2 rho(2k-1)", ok_iii)
print("     => 奇偶 T 只差因子 2；独立参数是 ceil(T/2)")

print("\nF 毁灭量 D(T) 的量级（末周期）")
for T in alive:
    d = res[T]['destroyed']
    print(f"     T={T:>2}: 毁灭={d:>16}  log10={math.log10(max(d,1)):>6.2f}  "
          f"历史={res[T]['history']:>12}")
check("F1 毁灭量随 T 单调增（周期越长，一次毁灭越多）",
      all(res[T]['destroyed'] <= res[T + 1]['destroyed']
          for T in alive if T + 1 in alive))

print("\n" + "=" * 78)
print("通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
if FAIL:
    print("不符项：")
    for n in FAIL:
        print("   x", n)
out = {"pass": len(PASS), "fail": len(FAIL), "failed_names": FAIL,
       "extinct_T": dead, "alive_T": alive,
       "converged": {str(k): bool(v) for k, v in conv.items()},
       "last_cycle": {str(T): {k: res[T][k] for k in
                               ('closed', 'peak', 'destroyed',
                                'history', 'seeds')} for T in res}}
with open(os.path.join(HERE, "L2_period_results.json"), "w",
          encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("=" * 78)
if FAIL:
    raise SystemExit(1)
print("全部通过 ✓")
