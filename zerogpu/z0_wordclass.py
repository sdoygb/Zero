"""
z0_wordclass.py --- 换状态表示：**词级**的闭合类，能不能区分顶点？

**为什么只测这一件事**：§8–§11 四次否定（度异质／度成块／模块度／类竞争保留）都失败，
根因是引擎只跟踪**计数与长度**。而长度谱是全局的 ⟹ 保留规则无从区分。
本轮把"闭合类"换成**真正的类**：闭合词（边序列）的**旋转类**。

**决定性的小问题**：每个顶点 $v$ 的**闭合类谱**（长度 $\le k$ 的不同旋转类的个数）是否因顶点而异？

$$
\begin{cases}
\text{谱因顶点而异} &\Longrightarrow\ \text{标签有区分度} \Longrightarrow\ D222\text{ 的保留机制有作用对象}\\
\text{谱逐点相同} &\Longrightarrow\ \text{连真词级标签也救不了（至少在这张 }\Gamma\text{ 上）}
\end{cases}
$$

**方法**：把每条走路的边序列编码成整数（底 $2d$），枚举长度 $\le k$ 的全部走路；
闭合的（$cur=start$）取**最小旋转**当类标签；按起点统计不同类的个数。

Γ：(a) $4\times4$ 环面（**顶点传递**，对照）；(b) 同图＋枢纽（**非传递**）；(c) 原生 $\mathcal G_{10}$。

用法：/usr/bin/python3 z0_wordclass.py  →  results/z0_wordclass.json
"""
from __future__ import annotations

import json
import os
import time
from collections import defaultdict

import numpy as np
from scipy.sparse import csr_matrix, identity, kron

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)


def torus_adj(D, L):
    r = np.arange(L)
    A1 = csr_matrix((np.ones(L), (r, (r + 1) % L)), shape=(L, L))
    A1 = A1 + A1.T
    A = A1
    for _ in range(D - 1):
        A = kron(A, identity(L, format="csr")) + kron(identity(A.shape[0], format="csr"), A1)
    return csr_matrix(A)


def add_hubs(A, frac=0.25, extra=2, seed=0):
    rng = np.random.default_rng(seed)
    V = A.shape[0]
    A = A.tolil()
    for u in rng.choice(V, size=max(1, int(frac * V)), replace=False):
        for _ in range(extra):
            v = int(rng.integers(0, V))
            if v != u:
                A[u, v] = 1.0
                A[v, u] = 1.0
    return csr_matrix(A)


def min_rotation(digits, k, base):
    """长度 k、底 base 的循环序列的最小旋转（整数）。"""
    best = cur = digits
    shift = base ** (k - 1)
    for _ in range(k - 1):
        cur = (cur % shift) * base + (cur // shift)
        if cur < best:
            best = cur
    return best


def class_spectra(A, kmax=8, cap=4_000_000):
    """
    枚举长度 <= kmax 的全部走路，返回每个顶点的闭合类集合。
    状态：(start, cur, seq)；seq 以 base=2d* 编码（用邻居槽位编号，保证每步一个数字）。
    """
    V = A.shape[0]
    indptr, indices = A.indptr, A.indices
    deg = np.diff(indptr)
    dmax = int(deg.max())
    base = dmax
    # 邻居槽位表
    slot = np.full((V, dmax), -1, np.int64)
    for v in range(V):
        slot[v, :deg[v]] = indices[indptr[v]:indptr[v + 1]]
    start = np.repeat(np.arange(V), dmax)
    sl = np.tile(np.arange(dmax), V)
    ok = slot.ravel() >= 0
    start, sl = start[ok], sl[ok]
    cur = slot.ravel()[ok]
    seq = sl.copy()
    start0 = np.repeat(np.arange(V), dmax)[ok]
    classes = defaultdict(set)
    counts = []
    for k in range(1, kmax + 1):
        cl = cur == start0
        if cl.any():
            for s, q in zip(start0[cl], seq[cl]):
                classes[int(s)].add(min_rotation(int(q), k, base))
        counts.append(int(cl.sum()))
        if len(cur) > cap:
            break
        nxt = slot[cur]                                   # (m, dmax)
        m = nxt.shape[0]
        valid = nxt >= 0
        digits = np.broadcast_to(np.arange(dmax)[None, :], (m, dmax))
        seq_rep = np.repeat(seq[:, None], dmax, axis=1)
        st_rep = np.repeat(start0[:, None], dmax, axis=1)
        cur2 = nxt[valid]
        seq2 = seq_rep[valid] * base + digits[valid]
        start2 = st_rep[valid]
        cur, seq, start0 = cur2, seq2, start2
    return classes, counts


def main():
    t0 = time.time()
    out = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
           "question": "真正的闭合类标签（旋转类）能否区分顶点？",
           "runs": []}

    L = 4
    base = torus_adj(2, L)
    graphs = [("torus4x4_transitive", base), ("torus4x4_hubs", add_hubs(base, 0.25, 2, seed=1))]
    try:
        import l0_closure as L0
        import observable_sweep as OS
        from zcl import Engine
        eng = Engine()
        reps, _ = L0.enumerate_necklaces(10, eng, chunk_bits=26, verbose=False)
        graphs.append(("G_10", OS.build_sparse(10, reps, eng)))
    except Exception as e:
        print("  [warn] 闭链图跳过:", e)

    for name, A in graphs:
        V = A.shape[0]
        cls, counts = class_spectra(A, kmax=8)
        ncls = np.array([len(cls.get(v, ())) for v in range(V)], float)
        print("=" * 88)
        print(f"### Γ={name}  V={V}  闭合词数按长度: {counts}")
        print(f"    闭合类个数：min={int(ncls.min())} max={int(ncls.max())} "
              f"均值={ncls.mean():.2f} std={ncls.std():.2f}  不同值的个数={len(set(ncls.tolist()))}")
        # 类标签本身的重叠：任取两点，类集合是否相同
        sets = [frozenset(cls.get(v, ())) for v in range(V)]
        uniq = len(set(sets))
        print(f"    不同的类集合个数 = {uniq} / {V}   "
              f"{'⟹ 顶点彼此可区分' if uniq > 1 else '⟹ 所有顶点类谱完全相同（无区分度）'}")
        # 共有类 vs 私有类
        allc = set().union(*sets) if sets else set()
        common = frozenset.intersection(*sets) if sets else frozenset()
        print(f"    总共 {len(allc)} 个类；所有顶点共有 {len(common)} 个；"
              f"私有/稀有类 {len(allc) - len(common)} 个")
        out["runs"].append({"graph": name, "V": V, "closed_counts": counts,
                            "n_classes_min": int(ncls.min()), "n_classes_max": int(ncls.max()),
                            "n_classes_std": round(float(ncls.std()), 3),
                            "distinct_class_sets": int(uniq),
                            "total_classes": len(allc), "common_classes": len(common),
                            "per_vertex_classes": [len(s) for s in sets][:40]})

    with open(os.path.join(OUT, "z0_wordclass.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("=" * 88)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_wordclass.json")
    print("读法：'不同的类集合个数 > 1' 才是 D222 保留机制能起作用的前提。")
    print("      顶点传递的图上必然 = 1（对照）；非传递的图上要看差多少。")


if __name__ == "__main__":
    main()
