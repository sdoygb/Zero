"""
z0_domains2.py --- Γ 的**模块度扫描**：尖峰 ⟶ 畴 的临界点在哪？

**上一轮的结论**：度成块不够（活动场与棋盘格不对齐，只有尖峰凝聚）。
本轮换成**模块度**这个真正的"成畴"控制量：

    环面 L×L 划成 blk×blk 的模块；**块内边全留**，**块间边以概率 ε 保留**。
    ε=0 ⟹ 16 个互不相连的模块（Q→1）；ε=1 ⟹ 完整环面（Q→0）。

**测**：活跃簇数 $N_{\rm clust}$、活动 Gini、有多少个模块"热"、图的模块度 $Q$、
以及 $N_{\rm clust}$ 随 ε 的过渡（找临界 ε_c）。

引擎与 `z0_local.py` 相同：全分支（归一化）＋ 闭合写记录 ＋ 局部异步毁灭（Kac 寿命 $\tau_v=2|E|/\deg v$）
＋ 从保留层重播种（重数封顶 2）。**全程零权重、零随机选择**。

用法：/usr/bin/python3 z0_domains2.py
输出：results/z0_domains2.json
"""
from __future__ import annotations

import json
import os
import time
from collections import deque

import numpy as np
from scipy.sparse import csr_matrix, identity, kron

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)


_EDGE_R = None


def modular_graph(L=24, blk=6, eps=1.0, seed=0):
    """
    ★ 嵌套族：每条**块间边**有一个固定的随机阈值 r_e（只依赖 seed），保留 ⟺ r_e < eps。
    eps 增大只**增加**边 ⟹ 各 eps 之间可比（上一版每个 eps 各砍各的，不可比）。
    """
    global _EDGE_R
    if _EDGE_R is None:
        rng = np.random.default_rng(seed)
        R = np.zeros((L, L, 2))
        for i in range(L):
            for j in range(L):
                for k in range(2):
                    R[i, j, k] = rng.random()
        _EDGE_R = R
    R = _EDGE_R
    rows, cols = [], []
    for i in range(L):
        for j in range(L):
            for k, (di, dj) in enumerate(((1, 0), (0, 1))):
                a, b = (i + di) % L, (j + dj) % L
                same = (i // blk == a // blk) and (j // blk == b // blk)
                if not same and R[i, j, k] >= eps:
                    continue
                rows.append(i * L + j); cols.append(a * L + b)
    rows = np.array(rows); cols = np.array(cols)
    A = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(L * L, L * L))
    A = A + A.T
    A.data[:] = 1.0
    return csr_matrix(A)


def clusters(occ, A):
    N = A.shape[0]
    indptr, indices = A.indptr, A.indices
    lab = np.full(N, -1, np.int64)
    c = 0
    for s in np.nonzero(occ)[0]:
        if lab[s] >= 0:
            continue
        lab[s] = c
        q = deque([int(s)])
        while q:
            u = q.popleft()
            for v in indices[indptr[u]:indptr[u + 1]]:
                if occ[v] and lab[v] < 0:
                    lab[v] = c
                    q.append(int(v))
        c += 1
    return lab, c


def gini(x):
    x = np.sort(np.asarray(x, float))
    return 0.0 if x.sum() == 0 else float((2 * np.arange(1, len(x) + 1) - len(x) - 1).dot(x)
                                          / (len(x) * x.sum()))


def graph_modularity(A, blk, L):
    """Newman 模块度（按块划分）。"""
    deg = np.asarray(A.sum(1)).ravel()
    m2 = deg.sum()
    comm = np.array([(i // L) // blk * (L // blk) + (i % L) // blk for i in range(A.shape[0])])
    coo = A.tocoo()
    same = comm[coo.row] == comm[coo.col]
    e = same.sum() / m2
    a = np.array([deg[comm == c].sum() for c in np.unique(comm)]) / m2
    return float(e - (a ** 2).sum())


def run(A, steps=3000, sample_every=100, seed=0):
    V = A.shape[0]
    deg = np.asarray(A.sum(1)).ravel()
    tau = np.maximum(1, np.round(0.6 * deg.sum() / np.maximum(deg, 1))).astype(np.int64)
    N = np.eye(V) * (1.0 / V)
    clock = np.zeros(V, np.int64)
    record = np.zeros(V)
    term = np.zeros(V)
    hist = []
    for t in range(steps):
        N = N @ A
        mx = N.max()
        if mx > 0:
            N /= mx
        cl = np.diag(N).copy()
        if cl.sum() > 0:
            record += cl
            N[np.arange(V), np.arange(V)] = 0.0
        clock += 1
        for v in np.nonzero(clock >= tau)[0]:
            col = N[:, v].copy()
            if col.sum() > 0:
                term[v] += col.sum()
                term -= col
                N[:, v] = 0.0
            N[v, v] += min(record[v], 2)
            clock[v] = 0
        if (t + 1) % sample_every == 0:
            act = N.sum(axis=1)
            hist.append((t + 1, gini(act), float(np.corrcoef(act, deg)[0, 1])
                         if np.std(deg) > 0 and np.std(act) > 0 else 0.0,
                         float((N.sum(0) - N.sum(1) + term).sum())))
    return N.sum(axis=1), hist, deg


def main():
    t0 = time.time()
    L, blk = 24, 6
    nmod = (L // blk) ** 2
    out = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
           "design": f"环面 {L}×{L}，{blk}×{blk} 模块 × {nmod}，块间边以概率 ε 保留；引擎同 z0_local",
           "runs": []}
    print(f"Γ: 环面 {L}×{L}，{nmod} 个 {blk}×{blk} 模块；块间边保留概率 ε")
    print(f"{'ε':>6} {'Q图':>7} {'度均值':>7} {'活动Gini':>9} {'簇数θ=1':>8} "
          f"{'热模块数':>8} {'corr(act,deg)':>13} {'Σd':>10}")
    for eps in (0.0, 0.05, 0.1, 0.2, 0.35, 0.6, 1.0):
        globals()['_EDGE_R'] = None
        A = modular_graph(L=L, blk=blk, eps=eps, seed=5)
        Q = graph_modularity(A, blk, L)
        act, hist, deg = run(A, steps=3000, sample_every=100)
        # 模块级活动
        modsum = np.zeros(nmod)
        modmax = np.zeros(nmod)
        for i in range(L):
            for j in range(L):
                c = (i // blk) * (L // blk) + (j // blk)
                modsum[c] += act[i * L + j]
                modmax[c] = max(modmax[c], act[i * L + j])
        hot = int(np.sum(modsum > modsum.mean()))
        lab, nc = clusters(act > act.mean(), A)
        # 簇与模块的对齐：每个簇落在几个模块里
        align = []
        for j in range(nc):
            mods = {(idx // L) // blk * (L // blk) + (idx % L) // blk
                    for idx in np.nonzero(lab == j)[0]}
            align.append(len(mods))
        row = {"eps": eps, "Q_graph": round(Q, 4), "deg_mean": round(float(deg.mean()), 3),
               "gini": round(float(hist[-1][1]), 4), "clusters_theta1": int(nc),
               "hot_modules": hot, "corr_act_deg": round(float(hist[-1][2]), 4),
               "sum_d": float(hist[-1][3]),
               "modules_per_cluster": align[:12],
               "hist": [{"t": h[0], "gini": round(h[1], 4)} for h in hist],
               "gini_last3": [round(h[1], 4) for h in hist[-3:]]}
        out["runs"].append(row)
        print(f"{eps:>6} {Q:>7.3f} {row['deg_mean']:>7.3f} {row['gini']:>9.4f} "
              f"{nc:>8} {hot:>8} {row['corr_act_deg']:>13.4f} {row['sum_d']:>10.2e}"
              f"   末3个Gini={row['gini_last3']}")
    print(f"\n用时 {round(time.time()-t0,1)}s → results/z0_domains2.json")
    print("读法：ε=0 时模块互不相连 ⟹ 活动被困在各模块内；ε 增大到某处，畴会併吞成一片。")
    print("      簇数从 ~nmod 掉到 1 的那个 ε_c 就是「成畴 ↦ 成坨」的临界点。")


if __name__ == "__main__":
    main()
