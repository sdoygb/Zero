"""
z0_class_exit.py --- Z3/Z4 规则改造：**按闭合类分流的退出 ＋ 类竞争保留**

**动机**：§8–§10 连续三次否定证明 —— 忠实引擎里**没有成畴机制**，因为畴需要双稳/正反馈，
而正反馈不能用权重实现（Z0③）。所以只能到**规则**里找。`D222` 原文给了现成的：

> 「每个局部历史层 $P_i$ 有从精确到概括的亚层，毁灭时**只保留最高两层**」
> 「**被删去的低层不能再影响下一代**」

**本轮的机制**（全程零权重、零随机选择）：

1. 每个顶点 $v$ 记住自己**闭合类的直方图**（按闭合长度 $k$ 分类）；
2. 时钟到 $\tau_v$ 时：清空 $v$ 的活动（散度照旧送终端层），但**重播种只在 $v$ 的类获得保留名额时进行**；
3. **保留名额 = 全局人口最多的前两个类**（`D222` 的"最高两层"）；
4. 未获保留 ⟹ $v$ **熄灭**（不再被播种）。

$$
\Longrightarrow\ \textbf{这是一个}\textbf{存活阈值}\text{机制：顶点靠"自己的类是否热门"竞争存活}
\ \Longrightarrow\ \textbf{可能出现存活/熄灭的畴}
$$

**测**：存活顶点数、存活集的连通块数（= **畴数**）、块大小分布、活动 Gini、保留类的身份随时间的变化。

用法：/usr/bin/python3 z0_class_exit.py  →  results/z0_class_exit.json
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


def torus_adj(D, L):
    r = np.arange(L)
    A1 = csr_matrix((np.ones(L), (r, (r + 1) % L)), shape=(L, L))
    A1 = A1 + A1.T
    A = A1
    for _ in range(D - 1):
        A = kron(A, identity(L, format="csr")) + kron(identity(A.shape[0], format="csr"), A1)
    return csr_matrix(A)


def add_hubs(A, frac=0.1, extra=6, seed=0):
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


def clusters_of(mask, A):
    V = A.shape[0]
    indptr, indices = A.indptr, A.indices
    lab = np.full(V, -1, np.int64)
    c = 0
    for s in np.nonzero(mask)[0]:
        if lab[s] >= 0:
            continue
        lab[s] = c
        q = deque([int(s)])
        while q:
            u = q.popleft()
            for v in indices[indptr[u]:indptr[u + 1]]:
                if mask[v] and lab[v] < 0:
                    lab[v] = c
                    q.append(int(v))
        c += 1
    return np.bincount(lab[lab >= 0]).tolist() if c else []


def gini(x):
    x = np.sort(np.asarray(x, float))
    return 0.0 if x.sum() == 0 else float((2 * np.arange(1, len(x) + 1) - len(x) - 1).dot(x)
                                          / (len(x) * x.sum()))


def run(A, keep_top=2, steps=4000, sample_every=200, seed=0, maxL=40):
    V = A.shape[0]
    deg = np.asarray(A.sum(1)).ravel()
    tau = np.maximum(1, np.round(0.6 * deg.sum() / np.maximum(deg, 1))).astype(np.int64)
    N = np.eye(V) * (1.0 / V)
    clock = np.zeros(V, np.int64)
    record = np.zeros(V)
    term = np.zeros(V)
    cls_hist = np.zeros((V, maxL), bool)          # 顶点 × 闭合长度类
    lam = np.full(V, -1, np.int64)                # 顶点的类标签
    alive = np.ones(V, bool)                      # 是否获得保留
    hist = []
    for t in range(steps):
        N = N @ A
        mx = N.max()
        if mx > 0:
            N /= mx
        dg = np.diag(N).copy()
        if dg.sum() > 0:
            record += dg
            idx = np.nonzero(dg)[0]
            k = min(maxL - 1, (t % maxL))
            cls_hist[idx, k] = True               # 记录闭合类（按当前步长近似）
            N[idx, idx] = 0.0
        clock += 1
        fire = np.nonzero(clock >= tau)[0]
        for v in fire:
            col = N[:, v].copy()
            if col.sum() > 0:
                term[v] += col.sum()
                term -= col
                N[:, v] = 0.0
            if cls_hist[v].any():
                lam[v] = int(np.argmax(cls_hist[v]))      # 主导闭合类
                cls_hist[v] = False
            if alive[v]:                                   # ★ 只在获得保留时才重播种
                N[v, v] += min(record[v], 2)
            clock[v] = 0
        # 保留名额：全局最热的前 keep_top 个类
        assigned = lam >= 0
        if assigned.any():
            cnt = np.bincount(lam[assigned], minlength=maxL)
            top = set(np.argsort(-cnt)[:keep_top].tolist())
            alive = np.array([(lam[v] < 0) or (lam[v] in top) for v in range(V)])
        if (t + 1) % sample_every == 0:
            act = N.sum(axis=1)
            szs = clusters_of(alive, A)
            hist.append({"t": t + 1, "alive": int(alive.sum()),
                         "clusters": len(szs), "top_sizes": sorted(szs, reverse=True)[:5],
                         "gini_alive": round(gini(act[alive]), 4) if alive.any() else None,
                         "gini_all": round(gini(act), 4),
                         "sum_d": float((N.sum(0) - N.sum(1) + term).sum()),
                         "retained": sorted(int(x) for x in
                                            set(lam[alive][lam[alive] >= 0].tolist()))[:6]})
    return hist, alive


def main():
    t0 = time.time()
    out = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
           "mechanism": "退出按闭合类分流；重播种只在类获保留名额（全局前 keep_top）时进行；未获保留 ⟹ 熄灭",
           "runs": []}
    L = 16
    base = torus_adj(2, L)
    graphs = [("torus_regular", base), ("torus_hub10pct", add_hubs(base, 0.10, 6, seed=3))]
    for name, A in graphs:
        for keep in (1, 2, 3):
            hist, alive = run(A, keep_top=keep, steps=4000, sample_every=200)
            print("=" * 100)
            print(f"### Γ={name}  V={A.shape[0]}  保留名额 keep_top={keep}")
            print(f"{'t':>6} {'存活顶点':>8} {'畴数':>6} {'最大几块':>22} {'活动Gini(存活)':>14} {'Σd':>10}")
            for h in hist[::max(1, len(hist) // 9)]:
                print(f"{h['t']:>6} {h['alive']:>8} {h['clusters']:>6} "
                      f"{str(h['top_sizes']):>22} {str(h['gini_alive']):>14} {h['sum_d']:>10.2e}")
            out["runs"].append({"graph": name, "keep_top": keep, "hist": hist})
    with open(os.path.join(OUT, "z0_class_exit.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("=" * 100)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_class_exit.json")
    print("读法：畴数 = 存活集在 Γ 上的连通块数。keep_top 越小，竞争越狠 ⟹ 存活越少、畴应当越清晰。")


if __name__ == "__main__":
    main()
