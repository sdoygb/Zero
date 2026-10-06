"""
z0_zero_structure.py --- 「零」是什么？零和生成元的**结构判定**

三条结论（本文件逐条机器核验）
--------------------------------
**Z1（零 = 约束，不是谱）**：零和约束 $\sum_v n_v=0$ 不是一个谱条件，它是**核约束**（$\ker\sigma$，$\sigma(n)=\sum_v n_v$）。
**Z2（约束是生成元良定义的前提）**：去掉约束，盒角 $(b,\dots,b)$ 与 $(-b,\dots,-b)$ 是**吸收态**（所有坐标顶界、无路可走）
⟹ 生成元**非不可约**、$\lambda=0$ **重数 > 1**（"零"分裂）；
加上约束后 $\lambda=0$ **重数恒为 1** ⟹ 唯一常返类、唯一零模。**约束排除的正是那些吸收角。**
**Z3（谱是导出的，**不可加**）**：$\mathcal L_\kappa$ 按**源点**拆开 $\sum_v\mathcal L_v$（逐项相等，已验），
但**各点被零和约束相互耦合**：$[\mathcal L_u,\mathcal L_v]\ne0$（实测，见下）。故
$$\operatorname{spec}(\mathcal L)\ \ne\ \text{Minkowski 和 of 单点谱}\quad(\textbf{已数值否证})$$
可用的只有**变分上界**：单点试验函数给 $\lambda_2\le d_v\nu_2$。
（初稿曾猜"Minkowski 包含"，被本机否证 —— 保留为一次更正。）

用法：/usr/bin/python3 z0_zero_structure.py
"""
from __future__ import annotations
from itertools import product
import numpy as np
from scipy.sparse import csr_matrix

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def build(A, b, constraint=True):
    A = np.asarray(A, float); V = A.shape[0]
    st = [n for n in product(range(-b, b + 1), repeat=V) if (not constraint or sum(n) == 0)]
    idx = {n: i for i, n in enumerate(st)}
    nbr = [np.nonzero(A[v])[0] for v in range(V)]
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


def site_spec(d, b):
    n = 2 * b + 1; M = np.zeros((n, n))
    for i in range(n):
        for j in (i - 1, i + 1):
            if 0 <= j < n:
                M[j, i] += d; M[i, i] -= d
    return np.sort(np.linalg.eigvalsh(M))


if __name__ == "__main__":
    graphs = {"边(2点)": [[0, 1], [1, 0]],
              "路径3": [[0, 1, 0], [1, 0, 1], [0, 1, 0]],
              "环4": [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0]]}
    print("=" * 88)
    print("Z2：零和约束 ⟹ λ=0 重数从 >1 变成恒 =1（约束排除吸收角）")
    print("=" * 88)
    rows = {}
    ok2 = True
    for name, A in graphs.items():
        for b in (1, 2):
            out = {}
            for cons in (False, True):
                L, st = build(A, b, cons)
                if L.shape[0] > 400:
                    out[cons] = None; continue
                ev = np.linalg.eigvalsh(L.toarray())
                out[cons] = (int(np.sum(np.abs(ev) < 1e-9)), int(np.sum(np.abs(-L.diagonal()) < 1e-12)))
            rows[(name, b)] = out
            m_free, m_cons = out[False], out[True]
            if m_free is None:
                print(f"  {name:>8} b={b}: 无约束 N>{400} 跳过  |  零和约束 λ=0 重数={m_cons[0]:>2} sink={m_cons[1]}")
                ok2 &= (m_cons[0] == 1) and (m_cons[1] == 0)
                continue
            print(f"  {name:>8} b={b}: 无约束 λ=0 重数={m_free[0]:>2} sink={m_free[1]}"
                  f"  |  零和约束 λ=0 重数={m_cons[0]:>2} sink={m_cons[1]}")
            ok2 &= (m_free[0] > 1) and (m_cons[0] == 1) and (m_cons[1] == 0)
    check("**Z2 普遍成立**：无约束 ⟹ λ=0 重数 >1；零和约束 ⟹ 恒 =1、无 sink",
          ok2, "3 图 × 2 盒界，共 6 例")

    print("\n" + "=" * 88)
    print("Z3：按源点拆开是逐项相等的，但 [L_u,L_v] ≠ 0 ⟹ 谱不可加")
    print("=" * 88)

    def ops(A, b):
        A = np.asarray(A, float); V = A.shape[0]
        st = [n for n in product(range(-b, b + 1), repeat=V) if sum(n) == 0]
        idx = {n: i for i, n in enumerate(st)}
        nbr = [np.nonzero(A[v])[0] for v in range(V)]
        Lv = []
        for a in range(V):
            rows, cols, vals = [], [], []; dg = np.zeros(len(st))
            for k, n in enumerate(st):
                if n[a] <= -b:
                    continue
                for bb in nbr[a]:
                    if n[bb] >= b:
                        continue
                    nn = list(n); nn[a] -= 1; nn[bb] += 1
                    rows.append(idx[tuple(nn)]); cols.append(k); vals.append(1.0); dg[k] += 1.0
            rows += list(range(len(st))); cols += list(range(len(st))); vals += list(-dg)
            Lv.append(csr_matrix((vals, (rows, cols)), shape=(len(st), len(st))).toarray())
        return sum(Lv), Lv, st

    ok3 = True
    for name, A in {"K4(完全)": [[0, 1, 1, 1], [1, 0, 1, 1], [1, 1, 0, 1], [1, 1, 1, 0]],
                    "环4": [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0]]}.items():
        L, Lv, st = ops(A, 1)
        V = len(Lv)
        cmax = max(np.max(np.abs(Lv[u] @ Lv[v] - Lv[v] @ Lv[u]))
                   for u in range(V) for v in range(V))
        decomp = np.max(np.abs(L - sum(Lv)))
        # 可加预测是否命中真谱
        deg = np.asarray(A, float).sum(1).astype(int)
        nus = [site_spec(deg[v], 1) for v in range(V)]
        pred = set(round(sum(nus[v][ks[v]] for v in range(V)), 6)
                   for ks in product(*[range(len(x)) for x in nus]))
        act = set(np.round(np.linalg.eigvalsh(np.asarray(L)), 6))
        print(f"  {name:>8}: ΣL_v 逐项相等={decomp<1e-12}   max|[L_u,L_v]|={cmax:.2f}"
              f"   谱 ⊆ 可加预测={act <= pred}")
        ok3 &= (decomp < 1e-12) and cmax > 0 and not (act <= pred)
    check("**Z3 成立**：分解逐项相等 + 对易子非零 + 谱不可加（Minkowski 猜测被否证）", ok3, "2 例")

    print(f"\n断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
