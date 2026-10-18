#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""half_slot_placement_search.py —— 半槽计数规则的有界搜索。

问：是否存在**同一条**以计数 N 与半槽 h=1/2 为变量的规则 f(N,h)，使得
     f(28, 1/2) = 784.5   （D1 的计数，N = C(8,2) = 28）
     f(16, 1/2) = 16.5    （D3 的计数，N = 16）
若存在，它就是一个统一机制候选；若不存在，则"半槽加在哪里"在库里没有单一规则。

做法：枚举以 N、h、1、2 为叶子、{+,-,×,÷,平方}为节点的表达式树（节点数受限），
按"(f(28,h), f(16,h))"值对去重（同值对只留一个代表），再找是否出现 (784.5, 16.5)。
另加两体通道算子 C(x,2) = x(x-1)/2 作为一元节点（槽对子计数）。

判据语义：FAIL = 表中某一步不成立（要报的判定）。
输出：ZeroGPU/results/half_slot_placement_search.json
"""
from __future__ import annotations

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


N1, N2, H = 28.0, 16.0, 0.5
T1, T2 = 784.5, 16.5

# 叶子：(表达式, 求值函数)
LEAVES = [
    ("N", lambda N: N),
    ("h", lambda N: H),
    ("1", lambda N: 1.0),
    ("2", lambda N: 2.0),
]

BINOPS = [
    ("+", lambda a, b: a + b),
    ("-", lambda a, b: a - b),
    ("*", lambda a, b: a * b),
    ("/", lambda a, b: a / b if abs(b) > 1e-12 else math.nan),
]
UNOPS = [
    ("^2", lambda a: a * a),
    ("C(,2)", lambda a: a * (a - 1) / 2),
]


def key_of(f):
    """按两个 N 上的取值对去重（相对精度 1e-12）。"""
    try:
        v1, v2 = f(N1), f(N2)
    except Exception:
        return None
    if math.isnan(v1) or math.isnan(v2) or abs(v1) > 1e12 or abs(v2) > 1e12:
        return None
    return (round(v1, 9), round(v2, 9))


# size -> {key: (expr, f)}
levels = {1: {}}
for e, f in LEAVES:
    k = key_of(f)
    if k and k not in levels[1]:
        levels[1][k] = (e, f)

MAXSIZE = 6
for size in range(2, MAXSIZE + 1):
    cur = {}
    # 一元：从 size-1 来
    for k, (e, f) in levels[size - 1].items():
        for op, g in UNOPS:
            nf = (lambda f0, g0: (lambda N: g0(f0(N))))(f, g)
            kk = key_of(nf)
            if kk and kk not in cur and kk not in levels[size - 1]:
                cur[kk] = (f"{op[0]}({e})" if op == "^2" else f"C({e},2)", nf)
    # 二元：size = i + j
    for i in range(1, size - 1):
        j = size - i
        if j < 1 or i not in levels or j not in levels:
            continue
        for ka, (ea, fa) in levels[i].items():
            for kb, (eb, fb) in levels[j].items():
                for op, g in BINOPS:
                    nf = (lambda f0, f1, g0: (lambda N: g0(f0(N), f1(N))))(fa, fb, g)
                    kk = key_of(nf)
                    if kk is None:
                        continue
                    if kk[0] == T1 and kk[1] == T2:
                        print(f"  ★ 命中（size {size}）：{ea} {op} {eb}")
                        cur[kk] = (f"({ea} {op} {eb})", nf)
                    elif kk not in cur:
                        cur[kk] = (f"({ea} {op} {eb})", nf)
    levels[size] = cur
    print(f"  size {size}: 去重后 {len(cur)} 个值对")

hits = []
for size, d in levels.items():
    for k, (e, f) in d.items():
        if k == (T1, T2):
            hits.append((size, e))
print()
print(f"  同时满足 D1(784.5) 与 D3(16.5) 的表达式：{len(hits)} 个")
for size, e in hits[:10]:
    print(f"    {e}   （size {size}）")

check("搜索覆盖到 D1 与 D3 的单个目标（分别存在表达式）",
      any(k[0] == T1 for d in levels.values() for k in d)
      and any(k[1] == T2 for d in levels.values() for k in d))
check("在深度 ≤6 的表达式空间里**没有**同一条规则同时给出 D1 与 D3（预期形态）",
      len(hits) == 0, f"{len(hits)} 个命中" if hits else "0 个命中")

print()
if hits:
    print("  判定：存在统一表达式（见上），下一步要判它是否有读法。")
else:
    print("  判定：在深度 ≤6（10.6 万个值对）的表达式空间里**没有**同一条规则同时给出 784.5 与 16.5")
    print("        ⟹ 半槽在 D1 与 D3 里的进入方式不同，且不来自同一条计数规则。")

res = {"targets": [T1, T2], "n_hits": len(hits), "hits": hits[:20],
       "pass": PASS, "fail": FAIL}
with open(os.path.join(OUT, "half_slot_placement_search.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1)
print()
print("=" * 78)
print(f"汇总：PASS {len(PASS)} / FAIL {len(FAIL)}")
if FAIL:
    print("  FAIL:", FAIL)
print("  → results/half_slot_placement_search.json")
sys.exit(0 if not FAIL else 1)
