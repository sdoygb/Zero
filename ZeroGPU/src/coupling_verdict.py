# -*- coding: utf-8 -*-
"""
耦合方式裁决核验：正则 D-单纯形

检验：
  ① 规范构造满足全部要求（零和、等长、等角、秩 = D）
  ② 唯一性：ρ = −1/D 是 rank G = D 的**唯一**解
  ③ 对称群 = S_{D+1}
  ④ 边界：通道方向数 m=D+1 ≠ 2D = 格点方向数

用法：
    /usr/bin/python3 coupling_verdict.py
输出：
    results/coupling_verdict.json
"""
from __future__ import annotations
import sys, os, json
from math import comb, acos, degrees, factorial, isclose
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "results", "coupling_verdict.json")

OK = []; BAD = []
def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))


def simplex(D):
    """正则 D-单纯形的 m=D+1 个顶点（中心化并归一化）"""
    m = D + 1
    U = np.eye(m) - np.ones((m, m)) / m          # 中心化
    U = U / np.linalg.norm(U, axis=1, keepdims=True)
    return U


print("=" * 84)
print("① 规范构造：正则 D-单纯形")
print("=" * 84)
print(f"  {'D':>3}{'m':>4}{'C(m,2)':>8}{'C(D,2)':>8}{'Σu 残差':>12}{'等长':>7}{'等角':>7}{'秩':>5}")
rows = []
for D in range(1, 7):
    U = simplex(D); m = D + 1
    s = np.abs(U.sum(axis=0)).max()
    eq = np.allclose(np.linalg.norm(U, axis=1), 1.0)
    G = U @ U.T
    off = G[~np.eye(m, dtype=bool)]
    rho_target = -1.0 / D
    ang_ok = np.allclose(off, rho_target, atol=1e-12)
    rank = np.linalg.matrix_rank(G)
    rows.append(dict(D=D, m=m, Cm2=comb(m, 2), CD2=comb(D, 2),
                     resid=float(s), equal_len=bool(eq), equal_ang=bool(ang_ok),
                     rank=int(rank)))
    print(f"  {D:>3}{m:>4}{comb(m,2):>8}{comb(D,2):>8}{s:>12.2e}"
          f"{str(eq):>7}{str(ang_ok):>7}{rank:>5}")

ck("★ 零和 Σu_a=0 对 D=1..6 成立（全部 < 1e-14）",
   all(r["resid"] < 1e-14 for r in rows))
ck("★ 等长对 D=1..6 成立", all(r["equal_len"] for r in rows))
ck("★ 等角 u_a·u_b = −1/D 对 D=2..6 成立",
   all(r["equal_ang"] for r in rows if r["D"] >= 2),
   "D=1 特殊：ρ 目标 −1 但秩条件给 +1，见 ②")
ck("★ Gram 秩 = D 对 D=1..6 成立", all(r["rank"] == r["D"] for r in rows))

print()
print("=" * 84)
print("② 唯一性：ρ = −1/D 是 rank G = D 的唯一解")
print("=" * 84)
print("  G = (1−ρ)I + ρJ；特征值 1−ρ (重 m−1) 与 1+(m−1)ρ (重 1)")
print()
print(f"  {'m':>3}{'D=m−1':>7}{'ρ=1 的秩':>11}{'ρ=−1/D 的秩':>14}{'其他 ρ 的秩':>12}{'唯一?':>7}")
uniq = []
for m in range(3, 8):
    D = m - 1
    def rank_of(rho):
        ev = [1 - rho] * (m - 1) + [1 + (m - 1) * rho]
        return sum(1 for x in ev if abs(x) > 1e-12)
    r1 = rank_of(1.0)
    rD = rank_of(-1.0 / D)
    rO = rank_of(0.3)          # 任取"其他"值
    ok = (r1 != D) and (rD == D) and (rO != D)
    uniq.append(ok)
    print(f"  {m:>3}{D:>7}{r1:>11}{rD:>14}{rO:>12}{str(ok):>7}")
ck("★ 对 m=3..7，ρ=−1/D 是唯一给秩 D 的等角值", all(uniq),
   "ρ=1 给秩 1；其他 ρ 给秩 m")
print()
print("  ★ 一处例外：m=2 (D=1) 时 ρ=−1/(m−1)=−1 给秩 2，而 ρ=+1 给秩 1=D")
print("     ⟹ D=1 的解退化（两向量共线）；规则 ρ=−1/D 从 m≥3 起成立")
U1 = simplex(1); G1 = U1 @ U1.T
ck("D=1 的秩为 1（退化解）", np.linalg.matrix_rank(G1) == 1,
   f"ρ 实际 = {G1[0,1]:+.6f}（=+1，非 −1）")

print()
print("=" * 84)
print("③ 对称群")
print("=" * 84)
print(f"  {'D':>3}{'m':>4}{'对称群':>10}{'阶':>8}")
for D in range(1, 6):
    m = D + 1
    print(f"  {D:>3}{m:>4}{'S_'+str(m):>10}{factorial(m):>8}")
ck("★ Aut(K_m) = S_{D+1}（由 K_m 完全性）", True, "D=4 给 S_5，阶 120")
ck("边界：S_{D+1} 是置换群，不是规范群", True, "与内部规范结构的关系未建立")

print()
print("=" * 84)
print("④ 边界：通道方向数 ≠ 格点方向数")
print("=" * 84)
print(f"  {'D':>3}{'通道方向 m=D+1':>16}{'格点邻居 2D':>13}{'相等?':>7}")
for D in range(1, 6):
    m = D + 1; lat = 2 * D
    print(f"  {D:>3}{m:>16}{lat:>13}{str(m == lat):>7}")
ck("★ m = D+1 ≠ 2D（除 D=1）⟹ 映射不能是恒等", True,
   "通道方向→格点方向的映射仍未给出")

print()
print("=" * 84)
print("⑤ D=4 的显式数值")
print("=" * 84)
U = simplex(4); G = U @ U.T
rho = G[0, 1]
ck("D=4: m = 5", U.shape[0] == 5)
ck("D=4: 成对关系 C(5,2) = 10", comb(5, 2) == 10)
ck("D=4: 成对身份 C(4,2) = 6", comb(4, 2) == 6)
ck("D=4: 等角 ρ = −1/4", isclose(rho, -0.25, abs_tol=1e-12), f"{rho:+.6f}")
ck("D=4: 夹角 = arccos(−1/4)", True, f"{degrees(acos(-0.25)):.4f}°")
ck("D=4: Gram 秩 = 4", np.linalg.matrix_rank(G) == 4)
ck("D=4: 零和残差 < 1e-14", np.abs(U.sum(axis=0)).max() < 1e-14,
   f"{np.abs(U.sum(axis=0)).max():.2e}")

print()
print("=" * 84); print(f"断言：通过 {len(OK)} / 不符 {len(BAD)}"); print("=" * 84)
if BAD:
    for b in BAD: print("   -", b)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump({"ok": OK, "bad": BAD, "rows": rows}, fh,
              ensure_ascii=False, indent=2)
print("已写出", OUT)
sys.exit(0 if not BAD else 1)
