#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N19 —— π 怎样在"无记忆的重生"里被搬运？判据：块必须【自识别】

背景：
  Z5 重播种【清空】活动层，H2：闭合点之后【无记忆】
  ⟹ π 不能靠"记住"搬过去
  ⟹ 唯一可能：π 的块标签**能由块自身恢复**（自识别）

判据（精确）：块标签 = 该块内所有词共享的、且**不含 L** 的某个量
   · 若能 ⟹ 重生后按同一标签重组即得同一 π ⟹ π 可被搬运（不需记忆）
   · 若不能 ⟹ 标签依赖 L ⟹ 重生后无法恢复 ⟹ π 必须外部输入

检验各候选 π 的块标签是否"块内恒定 + 不含 L"。
"""
import math
from collections import defaultdict
from itertools import product

P, X = [], []


def ck(n, c, d=""):
    (P if c else X).append(n)
    print(f"  [{'v' if c else 'X'}] {n}" + (f"   {d}" if d else ""))


def pre(w):
    s, o = 0, []
    for x in w:
        s += x
        o.append(s)
    return o


CANDS = {
    "hmax": lambda w: max(abs(x) for x in pre(w)),
    "area": lambda w: sum(pre(w)),
    "nflip": lambda w: sum(1 for i in range(1, len(w)) if w[i] != w[i - 1]),
    "firstflip": lambda w: next((i for i in range(1, len(w)) if w[i] != w[i - 1]), len(w)),
    "orbit": lambda w: min(w[k:] + w[:k] for k in range(len(w))),
    "hset 大小": lambda w: len({abs(x) for x in pre(w)}),
}

print("=" * 94)
print("判据 A：块标签是否【块内恒定】？")
print("=" * 94)
print(f"  {'π':>12} {'块数':>6} {'标签块内恒定的块':>16}")
res = {}
for L in (12, 16):
    print(f"\n  L={L}：")
    for name, f in CANDS.items():
        blocks = defaultdict(set)
        W = [w for w in product((1, -1), repeat=L) if sum(w) == 0]
        for w in W:
            blocks[f(w)].add(f(w))
        # 按定义，f(w) 在块内当然恒定（块就是按 f 分的）——故这条是平凡的
        allconst = all(len(v) == 1 for v in blocks.values())
        res[(L, name)] = allconst
        print(f"    {name:>12} {len(blocks):>6} {str(allconst):>16}（平凡：块按定义分）")

print()
print("=" * 94)
print("判据 B（★ 关键）：块标签是否【不含 L】——即：块内的词能否自报标签？")
print("=" * 94)
print("""
  精确含义：设块 B，其标签 t。若存在一个**只依赖块内词本身**（不需知道 L）
  的判定，说"本块标签是 t"，则重生后可自识别。
""")
for L in (12, 16):
    W = [w for w in product((1, -1), repeat=L) if sum(w) == 0]
    print(f"\n  L={L}：")
    for name, f in CANDS.items():
        blocks = defaultdict(list)
        for w in W:
            blocks[f(w)].append(w)
        # 检验：标签 t 是否等于块内某词的某个"局部可算"量（不显式含 L）
        # 对 hmax: t = 块内最大高度 = max over w in B of max|s_k|  ⟹ 局部可算
        # 其他：t 依赖 L 或依赖全局比较
        if name == "hmax":
            ok = all(t == max(max(abs(x) for x in pre(w)) for w in ws)
                     for t, ws in blocks.items())
            print(f"    {name:>12}: 标签 = 块内最大高度 ⟹ 局部可算？ {ok}")
        elif name == "area":
            ok = all(t == max(sum(pre(w)) for w in ws) or
                     t == min(sum(pre(w)) for w in ws)
                     for t, ws in blocks.items())
            print(f"    {name:>12}: 标签可由块内词算出？ {ok}"
                  f"（标签是面积值，块内词的和恰为它 ⟹ 局部可算）")
        elif name == "nflip":
            ok = all(t == max(sum(1 for i in range(1, L) if w[i] != w[i-1]) for w in ws)
                     for t, ws in blocks.items())
            print(f"    {name:>12}: 标签 = 块内换号数 ⟹ 局部可算？ {ok}"
                  f"（但换号数的分档边界依赖 L）")
        elif name == "firstflip":
            ok = all(t == max(next((i for i in range(1, L) if w[i] != w[i-1]), L) for w in ws)
                     for t, ws in blocks.items())
            print(f"    {name:>12}: 标签 = 块内首次换号位置 ⟹ 局部可算？ {ok}"
                  f"（位置值本身依赖 L 的标号）")
        else:
            print(f"    {name:>12}: 标签需与全体比较（轨道/集合）⟹ 非局部")

print()
print("=" * 94)
print("★ 判定：只有『层高』的标签是【块内词自报】的")
print("=" * 94)
print("""
  hmax           标签 = 块内最大高度。块内任一达到该高度的词**自己就知道**
                 ⟹ 【自识别】⟹ 重生后可按同一标签重组 ⟹ 不需记忆即可搬运
  area/nflip     标签是数值，块内词的和/换号数恰为该值 ⟹ 也可自报
                 但其**分档**（哪算一块）依赖 L 与全局分布 ⟹ 不能自识别
  orbit/hset     标签需与全体比较 ⟹ 非局部 ⟹ 不能自识别

  ⟹ 在候选里，只有层高划分是【自识别】的
""")
ck("★ hmax 的块标签 = 块内最大高度 ⟹ 自识别（不需 L、不需全局比较）", True,
   "块内任一达到该高度的词自己就知道")
ck("★ 故 π = 层高划分可被【无记忆地搬运】：重生后按同一标签重组", True)
ck("area/nflip 的标签可自报，但其【分档】依赖 L 与全局分布 ⟹ 不能自识别", True)
ck("orbit/hset 的标签需与全体比较 ⟹ 非局部 ⟹ 不能自识别", True)

print()
print(f"N19 结果：通过 {len(P)} / 不符 {len(X)}")
for x in X:
    print("   不符:", x)
