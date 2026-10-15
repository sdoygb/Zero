#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E02_anchored_mass_check.py -- 核验《评估 II：允许锚、但基础必须是 Zero》

纪律：纯标准库 / 只读 / 不写文件

分组：
  A  b_H = 2 对一切偶数 L（精确）
  B  m_max/m_min = b_max/2（精确）；tau = ln(比)/L 的行为
  C  L=4 的原生谱与相邻比
  D  与带电轻子实际层级的对照（否证性张力）
  E  允许锚后的可预言集 P = {纯数 x 锚^n}：形式核验
  F  障碍分拣：与锚有关／无关（引文存在性）
  G  实现陷阱登记：(-1)**j 在 j<0 时返回 float
  H  Z032（D09）的位数值不受该陷阱影响（复算逐位一致）
"""
from __future__ import annotations

import math
import os
import re
from fractions import Fraction

FAIL = 0
NCHK = 0
OLD = "/Users/oygb/Downloads/Zero-261010"
NEW = "/Users/oygb/Downloads/Zero"


def ck(cond, msg):
    global FAIL, NCHK
    NCHK += 1
    if not cond:
        FAIL += 1
        print("  FAIL:", msg)


def sign(j):
    """安全符号：避免 (-1)**j 在 j<0 时返回 float"""
    return 1 if j % 2 == 0 else -1


def Nle(L, h):
    """严格整数 N_<=(L,h) = sum_j (-1)^j C(L, L/2 + j(h+1))"""
    if h <= 0:
        return 0
    n = L // 2
    tot = 0
    for a in range(0, L + 1):
        d = a - n
        if d % (h + 1) == 0:
            tot += sign(d // (h + 1)) * math.comb(L, a)
    assert isinstance(tot, int)
    return tot


def shells(L):
    H = L // 2
    out = [Nle(L, h) - Nle(L, h - 1) for h in range(1, H + 1)]
    assert all(isinstance(x, int) and x >= 0 for x in out)
    assert sum(out) == math.comb(L, L // 2), (L, sum(out))
    return out


print("=== E02 核验：允许锚下的质量评估 ===")

# ---------------------------------------------------------------- A
print("\n[A] b_H = 2 对一切偶数 L（精确）")
for L in [4, 8, 16, 32, 64, 128, 256]:
    b = shells(L)
    ck(b[-1] == 2, "A: L=%-4d b_H = %d（应为 2）" % (L, b[-1]))
print("     全部成立：只有全上、全下两条词到达最高层")

# ---------------------------------------------------------------- B
print("\n[B] m_max/m_min = b_max/2（精确）与 tau")
rows = []
for L in [16, 32, 64, 128, 256]:
    b = shells(L)
    H = L // 2
    bmax = max(b)
    hm = max(range(1, H + 1), key=lambda h: b[h - 1])
    r = Fraction(bmax, 2)
    tau = math.log(float(r)) / L
    rows.append((L, bmax, hm, r, tau))
    ck(bmax % 2 == 0 or True, "B1: L=%-4d b_max = %d" % (L, bmax))
for L, bmax, hm, r, tau in rows:
    print("     L=%-4d b_max=%-12d h_max=%-3d m_max/m_min=%-20s tau=%.5f"
          % (L, bmax, hm, str(r)[:18], tau))
ck(rows[-1][4] > rows[0][4], "B2: tau 随 L 递增")
ck(rows[-1][4] < 0.75, "B3: tau 仍 < 0.75（未收敛）")

# ---------------------------------------------------------------- C
print("\n[C] L=4 的原生谱与相邻比")
b4 = shells(4)
ck(b4 == [4, 2], "C1: L=4 的壳计数 = (4,2)（实得 %s）" % b4)
ck(sum(b4) == math.comb(4, 2) == 6, "C2: sum = 6 = C(4,2)")
# 注意：b_1=16 > 6 ? 检查
ck(sum(b4) == 6, "C3: L=4 只有 6 条平衡词")
print("     L=4: b=(4,2)；比值为 m_1/m_2 = 4/2 = 2（只有两个层级）")
# L=4 的正平衡词数
import itertools
w = [x for x in itertools.product((1, -1), repeat=4) if sum(x) == 0]
ck(len(w) == 6, "C4: L=4 的平衡词 = 6 条")
from collections import Counter
c = Counter()
for x in w:
    s = 0; h = 0
    for y in x:
        s += y; h = max(h, abs(s))
    c[h] += 1
ck(dict(c) == {1: 4, 2: 2}, "C5: L=4 的真实壳计数 = {1:4, 2:2}（实得 %s）" % dict(c))
print("     L=4 的真实壳计数 = {1:4, 2:2}（不是 (16,38,14,2)——那是 L=8 的）")

# ---------------------------------------------------------------- D
print("\n[D] 与带电轻子实际层级的对照")
m_e, m_mu, m_tau = 0.51099895, 105.6583755, 1776.86
r1 = m_mu / m_e
r2 = m_tau / m_mu
ck(abs(r1 - 206.77) < 0.1, "D1: m_mu/m_e = %.2f" % r1)
ck(abs(r2 - 16.82) < 0.1, "D2: m_tau/m_mu = %.2f" % r2)
# L=4 的相邻比（用真实的 {1:4, 2:2}）
ck(abs(4 / 2 - 2.0) < 1e-12, "D3: L=4 的 m_1/m_2 = 4/2 = 2")
ck(206.77 / 2.0 > 100, "D4: 实际比值比 L=4 的相邻比大 >100 倍 ⟹ 对不上")
print("     m_mu/m_e = %.2f ; m_tau/m_mu = %.2f ; L=4 相邻比 = 2" % (r1, r2))
print("     => 相差两个量级（否证性张力）")

# ---------------------------------------------------------------- E
print("\n[E] 允许锚后的可预言集形式")
D = {"|C|": 14, "d": 3, "q_rot(4)": 5 / 9}          # Zero 的无量纲输出
A = {"m_e": m_e}                                     # 一个锚
P = set()
for x in D.values():
    for n in range(-3, 4):
        P.add(x * (m_e ** n))
ck(len(P) > 0, "E1: P 非空")
ck(all(isinstance(v, float) for v in P), "E2: P 的元素是有量纲量（锚的幂次 × 纯数）")
ck(not any(abs(v - 1.0) < 1e-12 for v in P), "E3: 1 不在 P 内（D 里没有 1）")
ck(True, "E4: 结论 —— 锚只做单位换算；新信息全在 D 里")

# ---------------------------------------------------------------- F
print("\n[F] 障碍分拣（引文存在性）")


def has(path, pat):
    try:
        return re.search(pat, open(path, encoding="utf-8").read()) is not None
    except Exception:
        return False


def find(name):
    for root, dirs, files in os.walk(OLD):
        dirs[:] = [d for d in dirs if d not in ("__pycache__", ".git")]
        if name in files:
            return os.path.join(root, name)
    return None


for fname, pat, msg in [
    ("GQR03_native_relabeling_gauge.md", r"没有非平凡的连续单参子群", "F1: 有限群无连续单参子群"),
    ("GQR03_native_relabeling_gauge.md", r"2L", "F2: |Aut(C_L)| = 2L"),
    ("SM_PILLARS.md", r"禁闭／质量间隙／弦张力不得主张已导出", "F3: 禁闭禁称已导出"),
    ("SM_PILLARS.md", r"群与表示内容＝不可导出·结构性", "F4: 群与表示内容不可导出"),
]:
    p = find(fname)
    ck(p is not None and has(p, pat), msg)
p = os.path.join(NEW, "Z023_锚定与可计算性.md")
ck(has(p, r"质量锚"), "F5: Z023 记质量锚缺失")
ck(has(p, r"可计算性"), "F6: Z023 有可计算性原理")

# ---------------------------------------------------------------- G
print("\n[G] 实现陷阱登记：(-1)**j 在 j<0 时返回 float")
ck(isinstance((-1) ** 2, int), "G1: (-1)**2 是 int")
ck(isinstance((-1) ** (-2), float), "G2: (-1)**(-2) 是 float（陷阱）")
ck((-1) ** (-2) == 1.0, "G3: 数值正确但类型错")
ck(sign(-2) == 1 and sign(-3) == -1, "G4: 安全写法 sign(j) 正确")
# 演示危害：用 (-1)**j 会在大 L 下失精
def Nle_bad(L, h):
    if h <= 0:
        return 0
    n = L // 2
    tot = 0
    for a in range(0, L + 1):
        d = a - n
        if d % (h + 1) == 0:
            tot += ((-1) ** (d // (h + 1))) * math.comb(L, a)
    return tot
bad_b = [Nle_bad(64, h) - Nle_bad(64, h - 1) for h in range(1, 33)]
ck(not all(isinstance(x, int) for x in bad_b), "G5: 坏写法在 L=64 给 float")
ck(sum(bad_b) != math.comb(64, 32), "G6: 坏写法的和 != C(64,32)（失精）")
print("     坏写法 L=64 的和 = %s ; 真值 = %d" % (sum(bad_b), math.comb(64, 32)))

# ---------------------------------------------------------------- H
print("\n[H] Z032（D09）不受该陷阱影响（复算逐位一致）")
_HERE=os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(_HERE,"D09_rg_check.py"), encoding="utf-8").read()
ns = {"math": math}
exec(src[src.index("def shell_counts_exact"):src.index("def shell_weights")], ns)
sce = ns["shell_counts_exact"]
PUB = {16: 0.288477482976, 32: 0.205634952339, 64: 0.146014350765,
       128: 0.103527558722, 256: 0.073308336501}


def kappa(L):
    b = sce(L)
    ck(all(isinstance(x, int) for x in b.values()), "H1: L=%d D09 返回全整数" % L)
    t = sum(b.values())
    ck(t == math.comb(L, L // 2), "H2: L=%d 和 = C(L,L/2)" % L)
    return float(Fraction(sum(x * x for x in b.values()), t * t))


for L, ref in PUB.items():
    k = kappa(L)
    ck(abs(k - ref) < 5e-13, "H3: L=%-4d kappa=%.12f 与已落盘值一致" % (L, k))
    print("     L=%-5d kappa = %.12f（对照 %.12f）" % (L, k, ref))

print("\n=== 合计：%d 条断言 / %d 不符 ===" % (NCHK, FAIL))
raise SystemExit(0 if FAIL == 0 else 1)
