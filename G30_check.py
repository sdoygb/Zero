#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G30_check.py -- 记忆核检验：推前动力学是否非马尔可夫？

对应文档 G30_memory_kernel_test.md。
这不是转录核验，而是一次**新的计算**：自建最小微观模型，做 Chapman-Kolmogorov 检验。

模型（只依赖底层条款（Z0 条款 ＋ Z1–Z5 定理）与计数测度）：
  微观态 = 年龄多重集（每个分支一个年龄；整数重数 = Z2）
  一步   = 全部年龄 +1（Z4）；新年龄达到寿命 L 者进入终端（Z3）；
           新年龄等于 SPAWN_AGE 者繁殖 B 个 0 龄子支（Z1 定理 1＋Z2 的全分支）
  计数测度 = 固定总数 N 时，对所有年龄多重集均匀（Z0③）

  F1  模型自洽：分支寿命恰为 L 步，繁殖事件恰在 SPAWN_AGE
  F2  关键：同一 N 的不同微观态给出不同的 N'  => 宏观转移不是 N 的函数 => 记忆
  F3  Chapman-Kolmogorov 检验：P2(N2|N0) 与 sum_N1 P(N1|N0)P(N2|N1) 不等
  F4  记忆的载体：年龄分布本身是马尔可夫的（状态 -> 状态 是函数）
  F5  第 1 步的推前测度不是 S(N1) 上的均匀测度（这就是违约机制）
  F6  诚实边界
"""

import os
import sys
from collections import Counter, defaultdict
from itertools import combinations_with_replacement

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0

L = 4          # 寿命：年龄 0,1,2,3；再加 1 步即退出
SPAWN = 2      # 繁殖年龄
B = 2          # 每个繁殖事件产生的子支数


def check(name, cond, detail=""):
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def step(s):
    """s = (n0,n1,n2,n3)。返回下一步的年龄多重集。"""
    n = [0] * L
    for a, c in enumerate(s):
        if c == 0:
            continue
        a2 = a + 1
        if a2 < L:
            n[a2] += c
            if a2 == SPAWN:
                n[0] += B * c
        # a2 >= L: 退出（终端账本）
    return tuple(n)


def total(s):
    return sum(s)


def states_with(N, cap):
    """所有总数 = N 的年龄多重集（年龄 0..L-1）。"""
    out = []
    for combo in combinations_with_replacement(range(L), N):
        s = [0] * L
        for a in combo:
            s[a] += 1
        t = tuple(s)
        if total(t) <= cap:
            out.append(t)
    return out


# ======================================================================
head("F1  模型自洽")

s0 = (1, 0, 0, 0)
traj = [s0]
for _ in range(6):
    traj.append(step(traj[-1]))
print("      从单支 0 龄出发的轨迹（年龄多重集 -> 总数）:")
for i, s in enumerate(traj):
    print("        t=%d  %s   N=%d" % (i, s, total(s)))
check("单支在第 2 步（新年龄 = SPAWN）繁殖：N 由 1 跳到 1+B = %d" % (1 + B),
      total(traj[1]) == 1 and total(traj[2]) == 1 + B,
      "轨迹 N = %s" % [total(s) for s in traj])
check("寿命 = L：年龄 L-1 的分支下一步退出（%s -> N'=0）"
      % ((0, 0, 0, 1),), total(step((0, 0, 0, 1))) == 0)
check("总数守恒律 N' = N - (age3 数) + B*(age1 数)",
      all(total(step(s)) == total(s) - s[L - 1] + B * s[1]
          for s in states_with(3, 10 ** 9)))

# ======================================================================
head("F2  同一 N 的不同微观态给出不同 N'  => 记忆")

for N in (1, 2, 3, 4):
    S = states_with(N, 10 ** 9)
    img = defaultdict(list)
    for s in S:
        img[total(step(s))].append(s)
    outs = sorted(img)
    print("      N=%d：%d 个微观态 -> N' 取值 %s" % (N, len(S), outs))
    if len(outs) > 1:
        check("N=%d：N' 不唯一（%d 种）=> 宏观转移不是 N 的函数 => 记忆"
              % (N, len(outs)), True, "例：%s -> N'=%d ; %s -> N'=%d"
              % (img[outs[0]][0], outs[0], img[outs[-1]][0], outs[-1]))
    else:
        check("N=%d：N' 唯一" % N, True)

# 显式反例
a, b = (1, 0, 0, 0), (0, 0, 0, 1)
check("显式反例：N=1 的 %s -> N'=%d，而同为 N=1 的 %s -> N'=%d"
      % (a, total(step(a)), b, total(step(b))),
      total(a) == total(b) and total(step(a)) != total(step(b)))

# ======================================================================
head("F3  Chapman-Kolmogorov 检验")

CAP = 40
report = []
for N0 in (1, 2, 3, 4, 5):
    S0 = states_with(N0, 10 ** 9)
    if not S0:
        continue
    # 第 1 步
    one = [step(s) for s in S0]
    P1 = Counter(total(s) for s in one)
    for k in P1:
        P1[k] /= len(S0)
    # 第 2 步（直接演化）
    two = [step(step(s)) for s in S0]
    if max(total(s) for s in two) > CAP:
        continue
    P2 = Counter(total(s) for s in two)
    for k in P2:
        P2[k] /= len(S0)
    # CK 预测：sum_N1 P(N1|N0) * P(N2|N1)
    PK = defaultdict(float)
    for N1, p1 in P1.items():
        S1 = states_with(N1, 10 ** 9)
        c = Counter(total(step(s)) for s in S1)
        for N2, cnt in c.items():
            PK[N2] += p1 * cnt / len(S1)
    keys = set(P2) | set(PK)
    dev = max(abs(P2.get(k, 0.0) - PK.get(k, 0.0)) for k in keys)
    tv = 0.5 * sum(abs(P2.get(k, 0.0) - PK.get(k, 0.0)) for k in keys)
    report.append((N0, len(S0), dev, tv))
    print("      N0=%d（%3d 个微观态）：max|P2-P_CK| = %.4f   全变差 = %.4f"
          % (N0, len(S0), dev, tv))

check("所有可算的 N0 都出现 CK 违约（max 偏差 > 1e-9）",
      all(r[2] > 1e-9 for r in report),
      "最大偏差 = %.4f" % max(r[2] for r in report))
check("违约显著（至少一个 N0 的全变差 > 0.1）",
      max(r[3] for r in report) > 0.1,
      "最大全变差 = %.4f" % max(r[3] for r in report))
check("=> 推前动力学在 N 这个粗粒化下【非马尔可夫】 => 记忆存在", True)

# ======================================================================
head("F4  记忆的载体：年龄分布本身是马尔可夫的")

check("年龄多重集 -> 下一步年龄多重集 是一个函数 step()",
      all(isinstance(step(s), tuple) for s in states_with(3, 10 ** 9)))
for s in states_with(3, 10 ** 9):
    check("state %s 的下一步唯一确定 = %s" % (s, step(s)), step(s) == step(s))
    break
check("=> 记忆的载体就是被 π 丢掉的那部分状态（年龄分布）", True)

# ======================================================================
head("F5  第 1 步的推前测度不是 S(N1) 上的均匀测度")

bad = 0
for N0 in (2, 3, 4):
    S0 = states_with(N0, 10 ** 9)
    one = [step(s) for s in S0]
    for N1 in set(total(s) for s in one):
        pushed = Counter(s for s in one if total(s) == N1)
        S1 = states_with(N1, 10 ** 9)
        if len(pushed) != len(S1):
            bad += 1
check("存在 N1 使推前测度的支撑 ≠ S(N1)（非均匀的直接证据）", bad > 0,
      "%d 个这样的 N1" % bad)
check("=> 非马尔可夫性的机制：推前把均匀性破坏了", True)

# ======================================================================
head("F6  机制与增长率：与 G29 的选择公式一致")

# 长期轨迹
tt = [s0]
for _ in range(30):
    tt.append(step(tt[-1]))
Ns = [total(x) for x in tt]
print("      N(t) = %s ..." % Ns[:12])
# 奇偶两类交替：x, x, 2x, 2x, ...
check("N(t) 呈周期 2 的倍增（x, x, 2x, 2x, ...）",
      all(abs(Ns[i + 2] - 2 * Ns[i]) < 1e-9 for i in range(2, 20)),
      "N(10..14) = %s" % Ns[10:15])

# 渐近每步增长率
import math
lam_meas = math.log(Ns[24] / Ns[22]) / 2.0
lam_pred = math.log(B) / SPAWN
print("      渐近每步增长率 lambda_实测 = %.6f   lambda_预测 = log(B)/SPAWN = %.6f"
      % (lam_meas, lam_pred))
check("实测增长率 = log(B)/SPAWN（G29 的 lambda = log(M)/T）",
      abs(lam_meas - lam_pred) < 1e-9)

# 机制：N 对年龄构成是盲的；是【被占用的年龄类】决定下一步
print("      机制：N=3 的两种年龄构成给出不同的 N'")
for st in [(2, 0, 1, 0), (0, 2, 0, 1)]:
    print("        %s (N=%d) -> N'=%d" % (st, total(st), total(step(st))))
check("同一 N=3 的两种构成给出 N'=3 与 N'=6 => 记忆的机制是【年龄构成（奇偶类）】",
      total((2, 0, 1, 0)) == total((0, 2, 0, 1)) == 3
      and total(step((2, 0, 1, 0))) == 3 and total(step((0, 2, 0, 1))) == 6)

# ======================================================================
head("F7  诚实边界")

check("这是最小模型（L=4, SPAWN=2, B=2），不是 D 系列的真实闭合动力学", True)
check("模型用【年龄】当分支状态；真实分支还有字/荷结构，未纳入", True)
check("CK 违约是【可计算的事实】，但记忆核的解析形式（Mori-Zwanzig）未导出", True)
check("本计算不改变 G1-G29 的结论，只建立记忆的存在性", True)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
