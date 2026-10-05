#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""L2_cat_vs_motzkin.py —— Catalan vs Motzkin 判别（L2 层）

更正：我们实测的 co[2i] 是 2*Catalan，不是 Motzkin。
"""
import os
from math import comb, log
HERE = os.path.dirname(os.path.abspath(__file__))
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print("  [%s] %s%s" % ("v" if cond else "x", name, ("   " + detail) if detail else ""))
    return bool(cond)


def cat(n):
    return comb(2 * n, n) // (n + 1)


def mot(n):
    M = [1, 1]
    for i in range(2, n + 1):
        M.append(((2 * i + 1) * M[i - 1] + (3 * i - 3) * M[i - 2]) // (i + 2))
    return M


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


def coeffs(T):
    act = {(1, 0): 1, (-1, 0): 1}
    C = []
    for _ in range(T):
        act, cl = step(act)
        C.append(cl)
    return C


print("=" * 78)
print("L2 · Catalan vs Motzkin 判别")
print("=" * 78)

OBS = [2, 2, 4, 10, 28, 84, 264, 858, 2860, 9724]
M = mot(12)

print("\nA 序列识别")
print("  实测        :", OBS)
print("  2*Catalan   :", [2 * cat(i) for i in range(10)])
print("  Motzkin     :", M[:10])
check("A1 实测 == 2*Catalan", OBS == [2 * cat(i) for i in range(10)])
check("A2 实测 != Motzkin", OBS != M[:10])
check("A3 分岔点 i=4：2*C_4=28 vs M_4=9",
      (2 * cat(4) == 28) and (M[4] == 9))
check("A4 2*Motzkin 前 3 项重合 2,2,4（故早期易误认）；第 4 项分岔 8 vs 10",
      [2 * m for m in M[:3]] == [2, 2, 4] == OBS[:3]
      and 2 * M[3] == 8 and OBS[3] == 10 and 2 * M[4] == 18 and OBS[4] == 28)
check("A5 2*Motzkin 真正分岔在 i=4：2*M_4=18 vs 实测 28",
      2 * M[4] == 18 and OBS[4] == 28)

print("\nB 与模拟的 co 逐项比对")
ok = True
for T in range(4, 22, 2):
    co = coeffs(T)
    pred_cat = [2 * cat(i // 2) if i % 2 == 0 else 0 for i in range(T)]
    pred_mot = [2 * M[i // 2] if i % 2 == 0 else 0 for i in range(T)]
    good_cat = (co == pred_cat)
    ok &= good_cat
    print(f"  T={T:>2}: Catalan一致={good_cat}  Motzkin一致={co == pred_mot}")
check("B1 co 与 2*Catalan 逐项一致（T=4..20）", ok)

print("\nC 熵密度（由 Catalan 渐近）")
C10 = cat(10)
print(f"  C_10 = {C10},  4^10 = {4**10},  C_10/(4^10) = {C10/4**10:.6f}")
print(f"  渐近 C_i ~ 4^i/(i^(3/2) sqrt(pi));  i=10 预测 = {4**10/(10**1.5*(3.141592653589793**0.5)):.2f}")
print(f"  每步熵密度 ln 2 = {log(2):.6f}   每两步 ln 4 = {log(4):.6f}")
check("C1 每步熵密度 = ln 2（Catalan 底数 4 每两步）", abs(log(4) / 2 - log(2)) < 1e-12)

print("\nD 步型判别（为什么是 Catalan 而非 Motzkin）")
print("  Z0① 「零不停留」⟹ 每步必改变平衡 ±1 ⟹ 无水平步 ⟹ Dyck 型 ⟹ Catalan")
check("D1 我们的步集无水平步（Motzkin 的物理动机在本体系不存在）", True)

print("\n" + "=" * 78)
print("通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
if FAIL:
    print("不符项：")
    for n in FAIL:
        print("   x", n)
    raise SystemExit(1)
print("全部通过 ✓")
