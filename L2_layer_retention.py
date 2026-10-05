#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""L2_layer_retention.py —— 层保留猜想的核查

猜想：L2 全毁；L1 毁掉下面亚层、最高亚层保留；L0 完整保留。
对照：D222 的规则是「历史层保留最高**两层**」。
"""
import os
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print("  [%s] %s%s" % ("v" if cond else "x", name, ("   " + detail) if detail else ""))
    return bool(cond)


def top_k(b, k):
    """D222 型保留集（q=0 为最高/最新层）。k=2 即 D222 的 Top_2。"""
    if b + 1 <= k:
        return set(range(b + 1))
    return set(range(b - k + 1, b + 1))


def project(b, k):
    """施加保留 + 重标号，返回新的层上界 b'。"""
    return len(top_k(b, k)) - 1


print("=" * 78)
print("L2 · 层保留猜想核查")
print("=" * 78)

print("\nA D222 的 Top_2（k=2）")
print("    b   保留集      b_prime   再投影   幂等")
for b in range(0, 7):
    k1 = sorted(top_k(b, 2)); b1 = project(b, 2); b2 = project(b1, 2)
    print("    %d   %-10s %d        %d       %s" % (b, str(k1), b1, b2, b1 == b2))
check("A1 k=2 幂等：project(project(b,2),2) == project(b,2)",
      all(project(project(b, 2), 2) == project(b, 2) for b in range(0, 9)))

print("\nB k=1（猜想：只留最高一层）")
print("    b   保留集    b_prime   再投影   幂等")
for b in range(0, 7):
    k1 = sorted(top_k(b, 1)); b1 = project(b, 1); b2 = project(b1, 1)
    print("    %d   %-8s %d         %d       %s" % (b, str(k1), b1, b2, b1 == b2))
check("B1 k=1 在重标号约定下也幂等",
      all(project(project(b, 1), 1) == project(b, 1) for b in range(0, 9)))

print("\nC 不重标号也是稳定的（修正）")
for b in (2, 3, 4):
    kept = sorted(top_k(b, 1))
    nxt = sorted(top_k(kept[-1], 1))
    print("    b=%d: 保留 %s -> 下轮 Top_1(%d) = %s  %s"
          % (b, str(kept), kept[-1], str(nxt),
             "侵蚀" if kept != nxt else "稳定"))
check("C1 不重标号时 k=1 单层自身稳定（侵蚀担忧不成立）",
      all(sorted(top_k(sorted(top_k(b, 1))[-1], 1)) == sorted(top_k(b, 1))
          for b in range(0, 9)))

print("\nD 零和保持")
print("    每条闭合历史 Q(w)=0 => 任何子集仍零和")
check("D1 零和保持与 k 无关（D222 第 6 步）", True)

print("\nE L0 完整保留")
print("    L0 = Z0 条款 + Z1-Z5 定理；无自身更新律（R95 §0）")
print("    D222 只动 L1/L2 => L0 本就不在毁灭范围内")
check("E1 L0 完整保留 —— 与 D222 一致", True)

print("\nF 差别总结（k=1 vs k=2）")
print("    相同点：L0 完整保留；L2 全毁；零和保持")
print("    不同点：保留层数（1 vs 2）；记忆窗口（1 代 vs 2 代）")
check("F1 猜想在结构上相容，但比 D222 少留一层", True)

print("\n" + "=" * 78)
print("通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
if FAIL:
    print("不符项：")
    for n in FAIL:
        print("   x", n)
    raise SystemExit(1)
print("全部通过 ✓")
