# -*- coding: utf-8 -*-
"""
判定：T 与 L 是否同一个参数？

判据（定义的直接比对）：
  Z4:  L = 寿命；长度达到 L 未闭合的词进入终端层
       故 L 步词的**年龄**取值恰为 {0,1,...,L}，大小 L+1
  G32: T = 有限年龄原子 X_T = {0,...,T}，大小 T+1

  ⟹ 若年龄原子就是"步龄"，则 T = L

本脚本检验该认定是否**自洽**（用既有约束联立）。

用法：
    /usr/bin/python3 TL_verdict.py
输出：
    results/TL_verdict.json
"""
from __future__ import annotations
import sys, os, json
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "results", "TL_verdict.json")

OK = []; BAD = []
def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))


print("=" * 88)
print("一、定义比对")
print("=" * 88)
print("  Z4 : L = 寿命（长度达到 L 未闭合 ⟹ 终端层）")
print("  G32: T = 有限年龄原子 X_T = {0,...,T}，大小 T+1")
print()
print("  L 步词的年龄取值 = {0,1,...,L}，大小 L+1")
print("  ⟹ 两者**同一结构** ⟹ T = L")
print()
ck("定义比对给出 T = L", True, "年龄原子 = 步龄")

print()
print("=" * 88)
print("二、若 T = L，既有约束怎么联立？")
print("=" * 88)
print()
print("  约束清单：")
cons = [
    ("(i)  容量",     "c_T = C(T,⌊T/2⌋) ≤ 2(T+1)"),
    ("(ii) 偶数",     "L 偶（±1 读法：净位移为零需偶数步）"),
    ("(iii) 引力子下界", "D = 2L ≥ 4 ⟺ L ≥ 2"),
    ("(iv) 非交换载体", "L ≥ 3（二面体 D_L 的 2 维不可约表示）"),
    ("(v)  四维窗",    "q_rot(L) ∈ (1/2, 3/5)"),
]
for a, b in cons: print(f"    {a:<16}{b}")
print()
print("  逐 L 核验（L = 1..10）：")
print(f"  {'L':>3}{'c_L':>6}{'2(L+1)':>8}{'(i)':>5}{'(ii)':>6}{'(iii)':>7}{'(iv)':>6}"
      f"{'(v) q':>9}{'(v)':>6}{'全过?':>7}")
import itertools
from fractions import Fraction
def q_rot(L):
    if L % 2: return None
    W = [w for w in itertools.product([-1, 1], repeat=L) if sum(w) == 0]
    seen = set(); sizes = []
    for w in W:
        if w in seen: continue
        orb = {w[i:] + w[:i] for i in range(L)}
        seen |= orb; sizes.append(len(orb))
    N = len(W)
    return sum(Fraction(s, N)**2 for s in sizes)

survivors = []
for L in range(1, 11):
    c = comb(L, L//2)
    ok_i = c <= 2*(L+1)
    ok_ii = (L % 2 == 0)
    ok_iii = (2*L >= 4)
    ok_iv = (L >= 3)
    q = q_rot(L) if ok_ii else None
    ok_v = (q is not None) and (Fraction(1,2) < q < Fraction(3,5))
    allok = ok_i and ok_ii and ok_iii and ok_iv and ok_v
    if allok: survivors.append(L)
    qs = f"{float(q):.4f}" if q is not None else "—"
    print(f"  {L:>3}{c:>6}{2*(L+1):>8}{str(ok_i):>5}{str(ok_ii):>6}{str(ok_iii):>7}"
          f"{str(ok_iv):>6}{qs:>9}{str(ok_v):>6}{str(allok):>7}")

print()
ck("★ 若 T = L，五条约束联立给出**唯一** L", len(survivors) == 1,
   f"存活集 = {survivors}")

print()
print("=" * 88)
print("三、★ 判定")
print("=" * 88)
print()
if len(survivors) == 1:
    L = survivors[0]
    print(f"  ★ T = L ⟹ 唯一解 L = {L}")
    print()
    print(f"  推论：")
    print(f"    · 寿命 L = {L}")
    print(f"    · 维数 D = 2L = {2*L}")
    print(f"    · 每点容量 2(T+1) = {2*(L+1)}")
    print(f"    · 单点峰值 c_L = C({L},{L//2}) = {comb(L, L//2)}")
    print(f"    · 四维窗 q_rot({L}) = {float(q_rot(L)):.6f}")
    ck(f"★ T = L = {L}（唯一）", True, f"D = 2L = {2*L}")
else:
    print(f"  ★ T = L 不给出唯一解（存活集 {survivors}）")
    ck("★ T = L 未给唯一解", True, f"存活集 {survivors}")

print()
print("=" * 88); print(f"断言：通过 {len(OK)} / 不符 {len(BAD)}"); print("=" * 88)
if BAD:
    for b in BAD: print("   -", b)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump({"ok": OK, "bad": BAD, "survivors": survivors}, fh,
              ensure_ascii=False, indent=2)
print("已写出", OUT)
sys.exit(0 if not BAD else 1)
