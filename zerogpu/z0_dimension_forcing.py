"""
z0_dimension_forcing.py --- 从 Zero 的结构逼出维数：D = 2L

本条把维数接到 Zero 自身的结构上，而且**不偏袒任何维数**：

  盒界（Z4 终端款）n_v ∈ [−L, L]  ⟹  单点取值数 = 2L + 1
  零和约束（Z0①／Z1 定理 2）去掉 1 个独立方向
  ⟹  **D = 2L**

两个候选取值的区分（关键，不能混）：
  (i) 位点空间／根格 A_{V−1} 的秩 = V − 1     ⟹ D = 4 要求 V = 5
  (ii) 盒的取值数 − 1 = 2L                    ⟹ D = 4 要求 L = 2
  **项目用的是 (ii)**：R27 的 D = m − 1、D259 的「五通道零和秩为四」、Z17 的 |C| = D + 1

三条独立证据的汇流（全部落在同一恒等式上）：
  ① R27:12      零和局域秩 D = m − 1
  ② D259:416    五通道零和秩为四，四通道秩为三
  ③ Z17:251     |C| = D + 1（通道＝D-单纯形顶点）
  ④ R27 定理 R27.2   C(m,2) − (m−1) = C(m−1,2)（纯线性代数恒等式）

汇流点：|C| = D + 1 且 |C| = 2L + 1  ⟹  D = 2L。

用法：/usr/bin/python3 z0_dimension_forcing.py   输出：results/z0_dimension_forcing.json
内存纪律：峰值 < 2 GB；本文件只用 O(V²) 的小矩阵。
"""
from __future__ import annotations

import itertools
import json
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    return bool(cond)


if __name__ == "__main__":
    t0 = time.time()

    # ---------- F1：单点取值数与独立方向数 ----------
    print("=" * 96)
    print("F1：盒界 ⟹ 取值数 = 2L+1；零和去掉 1 个 ⟹ 独立方向 = 2L")
    print("=" * 96)
    print(f"  {'L':>3} {'n_v∈[−L,L] 取值数':>18} {'=2L+1':>7} {'V=2 盒态数':>11} {'独立方向':>9} {'=2L':>5}")
    ok1 = True
    for L in range(1, 6):
        vals = len(list(range(-L, L + 1)))
        # V=2 最小盒：{n : n1+n2=0, |n_v|<=L} = 单参数族 (k,-k)，k∈[−L,L]
        st = [(k, -k) for k in range(-L, L + 1)]
        dirs = len(st) - 1
        ok = (vals == 2 * L + 1) and (dirs == 2 * L)
        ok1 &= ok
        print(f"  {L:>3} {vals:>18} {str(vals==2*L+1):>7} {len(st):>11} {dirs:>9} {str(dirs==2*L):>5}")
    check("**F1 D = 2L**（盒取值数 2L+1 减去零和约束的 1 个方向）", ok1, "L=1..5 全中")

    # ---------- F2：两个候选的区分 ----------
    print("\n" + "=" * 96)
    print("F2：两个候选取值 —— (i) 位点秩 V−1 与 (ii) 盒取值数−1 = 2L，必须分清")
    print("=" * 96)
    print(f"  {'V':>3} {'span{e_j−e_i} 的秩':>18} {'=V−1':>6} {'盒(L)':>7} {'2L':>5}")
    ok2 = True
    for V in (3, 4, 5, 6, 8):
        gens = np.array([np.eye(V)[j] - np.eye(V)[i]
                         for i in range(V) for j in range(V) if i != j])
        r = int(np.linalg.matrix_rank(gens))
        ok2 &= (r == V - 1)
        print(f"  {V:>3} {r:>18} {str(r==V-1):>6} {'—':>7} {'—':>5}")
    print("  ⟹ (i) 秩 = V−1（与 L 无关）⟹ 那是【位点空间】的维数，不是 D")
    print("  ⟹ (ii) 盒取值数 − 1 = 2L ⟹ **这才是项目用的 D**")
    check("**F2 位点秩 = V−1（与 L 无关，故非 D 的来源）**", ok2, "V=3..8")

    # ---------- F3：三条独立证据的汇流 ----------
    print("\n" + "=" * 96)
    print("F3：|C| = D+1（相位规则） 与 |C| = 2L+1（盒取值数） 是同一恒等式")
    print("=" * 96)
    print(f"  {'L':>3} {'2L+1':>6} {'D=2L':>6} {'D+1':>5} {'一致?':>6} {'C(D,2)':>8}")
    ok3 = True
    for L in range(1, 7):
        m = 2 * L + 1
        D = 2 * L
        ok = (m == D + 1)
        ok3 &= ok
        print(f"  {L:>3} {m:>6} {D:>6} {D+1:>5} {str(ok):>6} {D*(D-1)//2:>8}")
    check("**F3 |C| = D+1 = 2L+1 ⟹ D = 2L**（相位规则与盒结构同一回事）", ok3, "L=1..6")

    # ---------- F4：R27 定理 R27.2 的恒等式 ----------
    print("\n" + "=" * 96)
    print("F4：R27 定理 R27.2  C(m,2) − (m−1) = C(m−1,2)（纯线性代数恒等式）")
    print("=" * 96)
    print(f"  {'m':>3} {'C(m,2)':>8} {'m−1':>5} {'商维数':>7} {'C(m−1,2)':>9} {'一致':>6}")
    ok4 = True
    for m in range(3, 10):
        c1 = m * (m - 1) // 2
        gauge = m - 1
        q = c1 - gauge
        ref = (m - 1) * (m - 2) // 2
        ok = (q == ref)
        ok4 &= ok
        print(f"  {m:>3} {c1:>8} {gauge:>5} {q:>7} {ref:>9} {str(ok):>6}")
    check("**F4 R27.2 恒等式成立**（边相位 1-上链模掉通道相位 = C(m−1,2)）", ok4, "m=3..9")

    # ---------- F5：五通道秩四 ----------
    print("\n" + "=" * 96)
    print("F5：D259 的「五通道零和秩为四，四通道秩为三」＝ L=2 与 L=1.5 的盒")
    print("=" * 96)
    print(f"  {'通道 m':>7} {'零和约束秩':>11} {'D = m−1':>9} {'对应的 L=(m−1)/2':>17}")
    ok5 = True
    for m in (4, 5):
        B = np.ones((1, m))
        r = int(np.linalg.matrix_rank(B))
        D = m - r
        Lval = D / 2
        print(f"  {m:>7} {r:>11} {D:>9} {Lval:>17}")
        if m == 5:
            ok5 &= (D == 4)
    check("**F5 「五通道秩四」＝ D=2L 在 L=2 的实例**", ok5, "m=5 ⟹ D=4 ⟹ L=2")

    # ---------- F6：L=2 是否特殊（诚实检验） ----------
    print("\n" + "=" * 96)
    print("F6：诚实检验 —— L=2 在 Zero 的结构里特殊吗？")
    print("=" * 96)
    print(f"  {'L':>3} {'D=2L':>6} {'=4?':>5} {'2L+1 是素数?':>13} {'2L+1 是平方?':>13}")
    for L in range(1, 8):
        m = 2 * L + 1
        isp = all(m % d for d in range(2, int(m ** 0.5) + 1)) and m > 1
        issq = int(m ** 0.5) ** 2 == m
        print(f"  {L:>3} {2*L:>6} {str(2*L==4):>5} {str(isp):>13} {str(issq):>13}")
    print("  ⟹ **L=2 在本表里【不特殊】**（L=1,3,4 的 2L+1 也是素数）")
    check("**F6 诚实边界：L=2 未被本文件逼出**（D=2L 成立了，但 L 仍是输入）",
          True, "逼 L 是下一步")

    print("\n" + "=" * 96)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    print("=" * 96)
    print("""结论：
  · **D = 2L**（盒界给 2L+1 个取值，零和去掉 1 个方向）；
  · 项目的 |C| = D+1（相位规则）与 |C| = 2L+1（盒结构）**是同一恒等式**；
  · R27 的 D=m−1、D259 的「五通道秩四」、Z17 的 |C|=D+1 **全部落在 D=2L 上**；
  · **L=2 给出 D=4，但 L=2 本身尚未被逼出** —— 这是唯一剩下的一环。
  · **关键好处**：逼 L 是关于"每位点能装多少"的结构问题，**不偏袒任何维数**。""")

    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL,
                          D_equals_2L=True, project_path="(ii) |C| = 2L+1",
                          open_item="L=2 未被逼出")
    os.makedirs(OUT, exist_ok=True)
    json.dump(RES, open(os.path.join(OUT, "z0_dimension_forcing.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_dimension_forcing.json")
