"""
z0_local.py --- 局部**异步**毁灭 ＋ **Kac 寿命**（D222 的最后一件）

**为什么这样做**
    `Z0_CORE.md` §7 的结论：全局同步清空下，"稳定"与"有空间结构"互斥（要么爆炸，要么种子均匀）。
    而 `D222` 原文说：「活动亚层全清是**局部事件**；**全局同步不是前置条件**」，每个区域有自己的 $\\tau_i$。

**寿命从哪来（原生、无参数）**
    Kac 引理：从 $v$ 出发的随机游走，**平均首次返回时间** $=\\tau_v=1/\\pi(v)=2|E|/\\deg(v)$。
    于是
$$
\\tau_v=\\frac{2|E|}{\\deg(v)}
\\quad\\Longrightarrow\\quad
\\textbf{正则 }\\Gamma\\textbf{：所有 }\\tau_v\\textbf{ 相等（同步，退化成上一版）};\\qquad
\\textbf{非正则 }\\Gamma\\textbf{：这才真正异步}
$$
    **所以预言是尖锐的**：异步毁灭要产生新结构，$\\Gamma$ **必须非正则** —— 这直接接到 `L2-LATTICE-ORIGIN`。

**引擎（归一化避免溢出）**
    全分支 $N\\leftarrow N A$；每步除以最大值（**形状不变**，因为全分支归一化 = 简单随机游走）。
    闭合 $N[v,v]$ → 写记录（中性，I3）。
    每个顶点有自己的时钟：$c_v$ 到 $\\tau_v$ ⟹ **只清空该顶点的活动**（不是全图），
    把它们的散度送进终端层 $D$，再从该顶点**保留的记录层**重播种（Z5，重数封顶 $K$）。

**测什么**
    ① 种子/活动分布是否长出**空间结构**（Gini、与 $\\deg v$ 的相关、活跃簇数）
    ② 是否到**稳态**（统计量随时间的漂移）
    ③ 四条恒等式（I2 $\\sum_x d(x)=0$、I3 记录中性、零和源/汇补偿）

用法：/usr/bin/python3 z0_local.py
输出：results/z0_local.json
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


def add_hubs(A, frac=0.1, extra=4, seed=0):
    """把一部分顶点变成"枢纽"：加随机长边 ⟹ **非正则**（度分布有重尾）。"""
    rng = np.random.default_rng(seed)
    V = A.shape[0]
    A = A.tolil()
    hubs = rng.choice(V, size=max(1, int(frac * V)), replace=False)
    for u in hubs:
        for _ in range(extra):
            v = int(rng.integers(0, V))
            if v != u:
                A[u, v] = 1.0
                A[v, u] = 1.0
    return csr_matrix(A)


def gini(x):
    x = np.sort(np.asarray(x, float))
    n = len(x)
    if n == 0 or x.sum() == 0:
        return 0.0
    return float((2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


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
    return c


def run(A, tau_scale=1.0, K=2, steps=3000, sample_every=30, seed=0, name=""):
    """
    局部异步毁灭 ＋ Kac 寿命。tau_scale 只用来把 Kac 时间缩放到可跑范围（**比例保持**）。
    """
    V = A.shape[0]
    deg = np.asarray(A.sum(1)).ravel()
    E2 = float(deg.sum())                       # 2|E|
    tau = np.maximum(1, np.round(tau_scale * E2 / np.maximum(deg, 1))).astype(np.int64)
    rng = np.random.default_rng(seed)
    N = np.eye(V) * (1.0 / V)                   # seed：均匀（无偏好）
    clock = np.zeros(V, np.int64)
    record = np.zeros(V)
    term = np.zeros(V)
    rows = []
    for t in range(steps):
        N = N @ A
        mx = N.max()
        if mx > 0:
            N = N / mx                          # 归一化：形状不变（全分支归一 = 简单随机游走）
        cl = np.diag(N).copy()
        if cl.sum() > 0:
            record += cl
            N[np.arange(V), np.arange(V)] = 0.0
        clock += 1
        fire = np.nonzero(clock >= tau)[0]
        if len(fire):
            # 只清空这些顶点的活动；散度送终端层
            for v in fire:
                col = N[:, v].copy()
                if col.sum() > 0:
                    # ★ 正确写法：这些词 (u→v) 的散度 = e_v − e_u
                    #   端点边缘全在 v；起点边缘分布在各个 u 上
                    term[v] += col.sum()
                    term -= col
                    N[:, v] = 0.0
                # Z5：从该顶点保留的记录层重播种（重数封顶 K）
                N[v, v] += min(record[v], K)
                clock[v] = 0
        if (t + 1) % sample_every == 0:
            d = N.sum(axis=0) - N.sum(axis=1) + term
            act = N.sum(axis=1)                     # 每个起点的活动量 = 种子分布
            occ = act > act.mean() * 0.5
            rows.append({"t": t + 1, "sum_d": float(d.sum()),
                         "act_gini": round(gini(act), 5),
                         "act_maxfrac": round(float(act.max() / max(act.sum(), 1e-30)), 5),
                         "seed_record_gini": round(gini(record), 5),
                         "corr_act_deg": round(float(np.corrcoef(act, deg)[0, 1]), 4)
                         if (np.std(act) > 0 and np.std(deg) > 0) else None,
                         "clusters": int(clusters(occ, A)),
                         "fired_cum": int(np.sum(tau <= t + 1)),
                         "n_fire_events": int(len(fire))})
    return rows, {"name": name, "V": V, "deg_min": int(deg.min()), "deg_max": int(deg.max()),
                  "deg_std": round(float(deg.std()), 4),
                  "tau_min": int(tau.min()), "tau_max": int(tau.max()),
                  "tau_ratio": round(float(tau.max() / max(tau.min(), 1)), 3)}


def main():
    t0 = time.time()
    out = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
           "mechanism": "局部异步毁灭（每顶点自己的时钟）＋ Kac 寿命 τ_v = 2|E|/deg(v)",
           "prediction": "正则 Γ ⟹ 所有 τ_v 相等 ⟹ 退化为同步；非正则 Γ ⟹ 真异步 ⟹ 可能长出空间结构",
           "runs": []}

    L = 16
    base = torus_adj(2, L)
    graphs = [("torus_regular", base, 1.0)]
    for frac, extra in ((0.10, 4), (0.05, 12)):
        g = add_hubs(base, frac=frac, extra=extra, seed=3)
        graphs.append((f"torus_hub{int(frac*100)}pct_x{extra}", g, 0.25))

    for name, A, sc in graphs:
        rows, info = run(A, tau_scale=sc, K=2, steps=3000, sample_every=30, name=name)
        info["tau_scale"] = sc
        out["runs"].append({"info": info, "rows": rows})
        print("=" * 100)
        print(f"### Γ={name}  V={info['V']}  度 [{info['deg_min']},{info['deg_max']}] "
              f"std={info['deg_std']}   τ ∈ [{info['tau_min']},{info['tau_max']}] "
              f"比值={info['tau_ratio']}")
        print(f"{'t':>6} {'Σd':>10} {'活动Gini':>9} {'活动max':>9} {'记录Gini':>9} "
              f"{'corr(活动,度)':>13} {'活跃簇':>7}")
        for r in rows[::max(1, len(rows) // 10)]:
            cc = r["corr_act_deg"]
            print(f"{r['t']:>6} {r['sum_d']:>10.2e} {r['act_gini']:>9.4f} "
                  f"{r['act_maxfrac']:>9.4f} {r['seed_record_gini']:>9.4f} "
                  f"{(f'{cc:13.4f}' if cc is not None else '         —   ')} {r['clusters']:>7}")
        last = rows[-1]
        print(f"  → 末态：活动Gini={last['act_gini']}  corr(活动,度)={last['corr_act_deg']} "
              f"活跃簇={last['clusters']}  Σd={last['sum_d']:.2e}")

    with open(os.path.join(OUT, "z0_local.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("=" * 100)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_local.json")


if __name__ == "__main__":
    main()
