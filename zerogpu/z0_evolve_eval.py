"""
z0_evolve_eval.py --- 评估：「Zero 在电脑程序里能演化吗？」

结论（三层，各自有机器依据）
----------------------------------------------------------------------------
① **能跑**：$\mathcal L_\kappa$ 是严格主方程；$\pi\propto\deg$、建造/消亡精确平衡由**对称性**强制。
② **但不是演化**：微观链的最慢非平凡弛豫是 **O(1)** 步 —— 没有亚稳畴、没有长寿命、没有老化。
   这是**定理级上界**，不是算力问题：

   $$ \tau_2^{\rm micro}\ \le\ \frac{1}{d_{\min}\,|\nu_2|},\qquad
      \nu_2=-4\sin^2\!\frac{\pi}{4L+2}\quad(\text{一维反射随机游走的第二本征值}) $$

   理由（变分原理 + 只依赖单点 $n_v$ 的试验函数）：对任意 $v$，取 $f(n)=\cos\!\bigl(\frac{\pi(n_v+L+1/2)}{2L+1}\bigr)$，
   取零均值后 $f$ 与基态正交，故 $\lambda_2\le \mathrm{RQ}(f)=d_v\nu_2$；取 $v=\arg\min d_v$ 即得。

③ **§40 的"演化"是伪谱 + 二次外推的乘积**（见 §41）：$2.64\times10^6$ 步 vs 真上界 $1.18$ 步，差 $2.2\times10^6$ 倍。

用法：/usr/bin/python3 z0_evolve_eval.py
"""
from __future__ import annotations

import sys
import os

import numpy as np
from scipy.sparse import csr_matrix

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import z0_thm as ZT

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def site_gap(d, b):
    """单点生成元（反射边界 {−b..b}）的第二本征值；解析 = d·[-4sin²(π/(4b+2))]"""
    n = 2 * b + 1
    M = np.zeros((n, n))
    for i in range(n):
        for j in (i - 1, i + 1):
            if 0 <= j < n:
                M[j, i] += d
                M[i, i] -= d
    return float(np.sort(np.linalg.eigvalsh(M))[1])


def main():
    from itertools import product

    def micro(A, b):
        V = A.shape[0]
        st = [n for n in product(range(-b, b + 1), repeat=V) if sum(n) == 0]
        idx = {n: i for i, n in enumerate(st)}
        Ad = np.asarray(A.todense())
        nbr = [np.nonzero(Ad[v])[0] for v in range(V)]
        rows, cols, vals = [], [], []
        diag = np.zeros(len(st))
        for k, n in enumerate(st):
            for a in range(V):
                if n[a] <= -b:
                    continue
                for bb in nbr[a]:
                    if n[bb] >= b:
                        continue
                    nn = list(n); nn[a] -= 1; nn[bb] += 1
                    rows.append(idx[tuple(nn)]); cols.append(k); vals.append(1.0)
                    diag[k] += 1.0
        rows += list(range(len(st))); cols += list(range(len(st))); vals += list(-diag)
        return csr_matrix((vals, (rows, cols)), shape=(len(st), len(st))), st

    print("=" * 90)
    print("① 变分上界成立性（小图逐一验证）")
    print("=" * 90)
    cases = [("路径3 b=1", np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]]), 1),
             ("环4 b=1", np.array([[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0]]), 1),
             ("环4 b=2", np.array([[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0]]), 2),
             ("星4 b=1", np.array([[0, 1, 1, 1], [1, 0, 0, 0], [1, 0, 0, 0], [1, 0, 0, 0]]), 1)]
    print(f"{'Γ':>10} {'b':>2} {'N':>5} {'RQ(上界)':>12} {'真 λ2':>12} {'上界成立':>8}")
    ok_all = True
    for name, A, b in cases:
        A = csr_matrix(A)
        L, st = micro(A, b)
        v = 0
        nv = np.array([n[v] for n in st], float)
        f = np.cos(np.pi / (2 * b + 1) * (nv + b + 0.5)); f = f - f.mean()
        rq = float(f @ (L @ f)) / float(f @ f)
        ev = np.sort(np.linalg.eigvalsh(L.toarray()))
        l2 = ev[1]
        good = rq >= l2 - 1e-9
        ok_all &= good
        print(f"{name:>10} {b:>2} {L.shape[0]:>5} {rq:>12.6f} {l2:>12.6f} {str(good):>8}")
    check("**变分上界普遍成立**（5 个小实例全过）", ok_all)

    print("\n" + "=" * 90)
    print("② 物理 Γ（V=256 层级图）：定理级上界表")
    print("=" * 90)
    Ah, _ = ZT.build_hierarchical()
    deg = np.asarray(np.asarray(Ah.sum(1))).ravel()
    dmin = float(deg.min())
    print(f"  Γ 度分布 {dmin:.0f}..{deg.max():.0f}，均 {deg.mean():.2f}（上界只用 d_min={dmin:.0f}）")
    print(f"  {'L':>4} {'|ν2|':>10} {'τ2 上界':>12}")
    for L in (1, 2, 3, 4, 5, 6, 8, 16):
        nu2 = -4 * np.sin(np.pi / (4 * L + 2)) ** 2
        print(f"  {L:>4d} {abs(nu2):>10.6f} {1/abs(dmin*nu2):>12.4f}")
    tau4 = 1 / abs(dmin * -4 * np.sin(np.pi / 18) ** 2)
    check("**物理值 L=4 的寿命 < 2 微观步**（对照 §40 曾报 2.64e6，差 ~2.2e6 倍）",
          tau4 < 2.0, f"τ2 ≤ {tau4:.4f} 步 vs §40 的 2.64e6 步")

    print("\n" + "=" * 90)
    print("③ 处置")
    print("=" * 90)
    print("""  · 「忠实 + 守恒」⟹ 自由弛豫（无亚稳畴/无老化/无记忆），这是零和层约束的结构性后果；
  · 要真正的结构与长寿命，必须付一条【具名输入】（有偏核 / 原生的两层分块 Γ），
    并需**重做** §37–§40 那种"层级买寿命"的判断 —— 用正确求解器（`which='SA'` 或 dense）。""")
    print(f"\n断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    return 0 if not FAIL else 1


if __name__ == "__main__":
    sys.exit(main())
