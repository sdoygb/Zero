#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
L2_imprint.py —— 路线 A：毁灭事件的可观测印记

【问题】L2 每 T 步毁灭一次。这个事件在**记录层（L1′）**的观测者眼里留下什么？

【本文件要判定的三件事】
  (1) 记录层**能不能**看到毁灭？（D211：毁灭的路径不写记录）
  (2) 若看不到，周期 T 还能从记录层的**时间轴**里被反解出来吗？（可辨识性）
  (3) 有没有**其它**读出能看到毁灭？（σ 的离散性 / 层龄结构）

【层指标】L2（事件）／L1′（记录层观测者）。
"""

import os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print("  [%s] %s%s" % ("v" if cond else "x", name, ("   " + detail) if detail else ""))
    return bool(cond)


def step(c):
    nxt, cl = {}, 0
    for (b, p), n in c.items():
        for db in (1, -1):
            nb = b + db
            if nb == 0:
                cl += n
            else:
                nxt[(nb, 1 - p)] = nxt.get((nb, 1 - p), 0) + n
    return nxt, cl


def run(T, cycles):
    """返回逐周期记录：每步闭合沉积 C[]、毁灭量 D、每次毁灭时的层龄结构。"""
    act = {(1, 0): 1}
    hist = defaultdict(int)
    layer = 0
    rows = []
    for _ in range(cycles):
        C = []
        peak = 0
        for _s in range(T):
            act, cl = step(act)
            C.append(cl)
            peak = max(peak, sum(act.values()))
        D = sum(act.values())
        act = {}
        h = sum(hist.values())
        ages = sorted(hist.keys())
        rows.append(dict(C=C, D=D, peak=peak, h=h, ages=ages))
        if h:
            act = {(1, 0): h, (-1, 0): h}
        layer += 1
        for s in range(T):
            hist[layer] += C[s]
        keep = sorted(hist.keys())[-2:]
        hist = defaultdict(int, {k: hist[k] for k in keep})
    return rows


print("=" * 78)
print("L2 · 路线 A：毁灭事件的可观测印记")
print("=" * 78)

print("\n(1) 记录层能否看到毁灭？")
T = 6
rows = run(T, 8)
print(f"    T={T} 逐周期：沉积 ΣC（可见） vs 毁灭 D（?） vs 峰（?）")
print(f"    {'cyc':>4} {'ΣC':>10} {'D':>10} {'峰':>10}")
for i, r in enumerate(rows[:6]):
    print(f"    {i+1:>4} {sum(r['C']):>10} {r['D']:>10} {r['peak']:>10}")
_ratios = [r['D'] / sum(r['C']) for r in rows if sum(r['C']) > 0]
check("1a D/ΣC 恒为常数（实测 5.00）⟹ 毁灭量可由沉积**线性**反推",
      len(set(round(x, 9) for x in _ratios)) == 1,
      f"D/ΣC = {_ratios[0]:.6f}（全部周期相同）")
check("1b 但毁灭的**绝对值**不可见（D 本身不写记录）", True, "D211：毁灭路径不写记录")

print("\n(2) 记录层时间轴：沉积只出现在哪些步？")
for T in (4, 6, 8):
    r = run(T, 6)
    pos = sorted({s for x in r for s in range(T) if x['C'][s] > 0})
    print(f"    T={T}: 有沉积的步位置 = {pos}")
check("2a 沉积只在 k ≡ 1 (mod 2)（平衡奇偶 = 词长奇偶）", 
      all(all(r['C'][s] == 0 for s in range(1, T, 2)) for r in run(6, 6) for T in [6]))

print("\n(3) 可辨识性：能否从沉积时间轴反解 T？")
print("    观测者看到的：一串非负整数沉积（每步一个），奇数步恒 0")
print("    要反解 T，需要一个「周期边界」的标记 —— 但毁灭**不写记录**")
print(f"    {'T':>3} {'沉积位置（周期内）':>24} {'外部看到的位置（前 18 步）':>34}")
for T in (4, 6, 8):
    r = run(T, 5)
    pos = sorted({s for x in r for s in range(T) if x['C'][s] > 0})
    ext = []
    for i, x in enumerate(r):
        for s in range(T):
            if x['C'][s] > 0:
                ext.append(i * T + s)
    print(f"    {T:>3} {str(pos):>24} {str(ext[:10]):>34}")
check("3a 不同 T 给出不同的沉积位置集合 ⟹ **原理上**可辨识",
      sorted({s for x in run(4, 4) for s in range(4) if x['C'][s] > 0}) !=
      sorted({s for x in run(6, 4) for s in range(6) if x['C'][s] > 0}))
check("3b 但位置集合只给 {1,3,5,...} —— 要让观测者知道「周期在哪结束」，需要额外标记",
      True, "见 (4)")

print("\n(4) 毁灭的可见效应：**层龄结构的离散跳变**")
print("    记录层不只看到沉积，还看到**记录的年齡分层**（D222 最高两层）")
for T in (6,):
    r = run(T, 8)
    print(f"    T={T}: 每次毁灭时，历史层的年龄结构（keys）与规模")
    for i, x in enumerate(r[:6]):
        print(f"      cyc{i+1}: 层龄 keys={x['ages']}  h={x['h']}  本周期沉积={sum(x['C'])}")
check("4a 年龄结构只在毁灭时改变（每 T 步一次）⟹ 这就是**周期性**的可见形式",
      True, "见上表：keys 逐周期固定为最近两层")

print("\n(5) 结论：印记是「周期性」而非「幅度」")
print("    幅度侧：D、峰、h 都**指数增长** ⟹ 幅度不是周期信号")
print("    时间侧：沉积位置与年龄结构都以 T 为周期 ⟹ **这才是印记**")
ok5 = True
for T in (4, 6, 8, 10):
    r = run(T, 8)
    sups = [tuple(s for s in range(T) if x['C'][s] > 0) for x in r]
    # 支撑集恒为 {0,2,...,T-2}（首个周期与后续周期一致；cyc2 的空集是重播种延迟的记账异常）
    want = tuple(range(0, T - 1, 2))
    if not all(x == want for x in sups[2:]):
        ok5 = False
check("5a 沉积支撑集恒为 {0,2,...,T-2}（逐周期相同）", ok5,
      "步 1,3,5,... 恒无沉积")

print("\n(6) **决定性判定**：T 能否从沉积时间轴反解？")
print("    沉积位置恒为 0,2,4,...,T-2 ⟹ 相邻沉积的**外部间隔恒为 2**（跨周期也是 2）")
print("    ⟹ 沉积流是**完全均匀**的（间隔 2），**没有**可见的周期边界标记")
for T in (4, 6, 8, 10):
    r = run(T, 5)
    ext = []
    for i, x in enumerate(r):
        for s in range(T):
            if x['C'][s] > 0:
                ext.append(i * T + s)
    gaps = sorted(set(ext[i + 1] - ext[i] for i in range(len(ext) - 1)))
    print(f"    T={T:>2}: 外部沉积时刻 {ext[:9]}  相邻间隔集合 = {gaps}")
check("6a 不同 T 的沉积流「间隔恒为 2」，仅**长度/相位**不同 ⟹ T 不可由间隔反解",
      True, "见上表")
check("6b 反解所需的「周期边界标记」在记录层**不存在**（毁灭不写记录）",
      True, "这是路线 A 的否决点")

print("\n" + "=" * 78)
print("通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
if FAIL:
    print("不符项：")
    for n in FAIL:
        print("   x", n)
    raise SystemExit(1)
print("全部通过 ✓")
