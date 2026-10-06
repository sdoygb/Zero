"""
z0_domains.py --- 平滑调制 还是 真的成畴？（下一步）

**上一轮留下的问题**：Kac 时钟 $\tau_v=2|E|/\deg(v)$ 让活动与度**反相关**（$-0.44$），
但"活跃簇"只有 1–3 个 ⟹ 看起来是**平滑调制**，不是**畴**。要判定，得做两件事：

1. **让度在空间上成块**：如果结构只是 $\deg$ 的平滑函数，那么度成块 ⟹ 活动也成块，
   **块界应当陡**；如果仍是平滑的，说明没有畴。
2. **判据要扫阈值**：簇数 $N_{\rm clust}(\theta)$ 对阈值 $\theta$ 扫一遍 ——
   有畴则存在一段 $\theta$ 区间里簇数**平台**（多块并存）；平滑场则恒为 1。

**Γ 的构造**（棋盘格成块加边）：环面 $L{=}32$ 划 $8{\times}8$ 的块，**棋盘格一半的块加密**
（块内随机加边，度 $4\to\sim12$）⟹ 度场**成块**、$\tau$ 比值 $\sim3$。

**测**：度场图、活动场图、$\mathrm{act}$–$\deg$ 散点、$N_{\rm clust}(\theta)$ 扫描、块内/块外对比、
以及**收敛性**（跑长看 Gini 是否停）。

用法：/usr/bin/python3 z0_domains.py
输出：results/z0_domains.json ＋ zerogpu/z0_domains.png
"""
from __future__ import annotations

import json
import os
import time
from collections import deque

import numpy as np
from scipy.sparse import csr_matrix, identity, kron

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams["font.sans-serif"] = ["Heiti TC", "Arial Unicode MS", "Hiragino Sans GB"]
matplotlib.rcParams["axes.unicode_minus"] = False

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


def blocky_graph(L=32, blk=8, extra=8, seed=0):
    """棋盘格一半的块**加密**：度场成块。返回 (A, degree, block_id, dense_mask)。"""
    rng = np.random.default_rng(seed)
    A = torus_adj(2, L).tolil()
    nb = L // blk
    dense = np.zeros((L, L), bool)
    for i in range(L):
        for j in range(L):
            if ((i // blk) + (j // blk)) % 2 == 0:
                dense[i, j] = True
    idx = lambda i, j: i * L + j
    for i in range(L):
        for j in range(L):
            if not dense[i, j]:
                continue
            for _ in range(extra):
                di, dj = rng.integers(0, blk, 2)
                a = (i // blk) * blk + di
                b = (j // blk) * blk + dj
                u, v = idx(i, j), idx(a, b)
                if u != v:
                    A[u, v] = 1.0
                    A[v, u] = 1.0
    A = csr_matrix(A)
    deg = np.asarray(A.sum(1)).ravel()
    return A, deg, dense


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


def gini(x):
    x = np.sort(np.asarray(x, float))
    n = len(x)
    return 0.0 if x.sum() == 0 else float((2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


def run(A, steps=6000, sample_every=100, seed=0):
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
        fire = np.nonzero(clock >= tau)[0]
        for v in fire:
            col = N[:, v].copy()
            if col.sum() > 0:
                term[v] += col.sum()
                term -= col
                N[:, v] = 0.0
            N[v, v] += min(record[v], 2)
            clock[v] = 0
        if (t + 1) % sample_every == 0:
            act = N.sum(axis=1)
            hist.append({"t": t + 1, "gini": round(gini(act), 5),
                         "corr": round(float(np.corrcoef(act, deg)[0, 1]), 4),
                         "sum_d": float((N.sum(0) - N.sum(1) + term).sum())})
    return N.sum(axis=1), hist, tau


def main():
    t0 = time.time()
    L, blk = 32, 8
    A, deg, dense = blocky_graph(L=L, blk=blk, extra=8, seed=0)
    V = A.shape[0]
    tau = np.maximum(1, np.round(0.6 * deg.sum() / np.maximum(deg, 1))).astype(np.int64)
    print(f"Γ: 环面 {L}×{L} 成块加密（块 {blk}×{blk}，棋盘格一半）")
    print(f"   度 [{deg.min()},{deg.max()}]  均值 {deg.mean():.2f}  "
          f"加密块内均值 {deg[dense.ravel()].mean():.2f} vs 块外 {deg[~dense.ravel()].mean():.2f}")
    print(f"   τ 比值 = {tau.max()/tau.min():.2f}")
    act, hist, _ = run(A, steps=6000, sample_every=100)
    print(f"\n{'t':>6} {'活动Gini':>9} {'corr(活动,度)':>13} {'Σd':>10}")
    for h in hist[::max(1, len(hist) // 8)]:
        print(f"{h['t']:>6} {h['gini']:>9.4f} {h['corr']:>13.4f} {h['sum_d']:>10.2e}")

    # 阈值扫描
    ths = [0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0, 1.2, 1.5, 2.0, 3.0]
    nc = [clusters(act > th * act.mean(), A) for th in ths]
    print("\n阈值扫描 N_clust(θ)：")
    for th, c in zip(ths, nc):
        print(f"   θ={th:<4} → {c} 块   （占比 {float((act>th*act.mean()).mean()):.3f}）")

    inside = act[dense.ravel()].mean()
    outside = act[~dense.ravel()].mean()
    print(f"\n块内平均活动 / 块外 = {inside/outside:.4f}   （度比 {deg[dense.ravel()].mean()/deg[~dense.ravel()].mean():.3f}）")
    print(f"末态 Gini={hist[-1]['gini']}  corr={hist[-1]['corr']}  Σd={hist[-1]['sum_d']:.2e}")

    res = {"graph": {"L": L, "blk": blk, "V": int(V), "deg_min": int(deg.min()),
                     "deg_max": int(deg.max()), "tau_ratio": round(float(tau.max() / tau.min()), 3)},
           "hist": hist, "theta_scan": {"thetas": ths, "clusters": nc},
           "inside_outside_activity_ratio": round(float(inside / outside), 4),
           "deg_ratio": round(float(deg[dense.ravel()].mean() / deg[~dense.ravel()].mean()), 4)}
    with open(os.path.join(OUT, "z0_domains.json"), "w") as f:
        json.dump(res, f, ensure_ascii=False, indent=1)

    # ---- 图 ----
    fig, ax = plt.subplots(1, 4, figsize=(19, 4.6))
    ax[0].imshow(deg.reshape(L, L), cmap="magma", interpolation="nearest")
    ax[0].set_title("Γ 的度场 $\\deg v$（成块加密）")
    ax[1].imshow(act.reshape(L, L), cmap="viridis", interpolation="nearest")
    ax[1].set_title(f"活动场 $\\mathrm{{act}}_v$（Gini={hist[-1]['gini']}）")
    ax[2].scatter(deg, act, s=4, alpha=0.4)
    ax[2].set_xlabel("$\\deg v$"); ax[2].set_ylabel("活动")
    ax[2].set_title(f"活动 vs 度（corr={hist[-1]['corr']}）")
    ax[3].semilogx(ths, nc, "o-")
    ax[3].set_xlabel("阈值 $\\theta$（相对均值）"); ax[3].set_ylabel("活跃簇数")
    ax[3].set_title("$N_{\\rm clust}(\\theta)$：有平台 ⟹ 成畴")
    for a in ax[:2]:
        a.set_xticks([]); a.set_yticks([])
    fig.suptitle("局部异步毁灭（Kac 寿命）：平滑调制 还是 真的成畴？", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    png = os.path.join(HERE, "z0_domains.png")
    fig.savefig(png, dpi=110)
    print(f"\n图：{png}   用时 {round(time.time()-t0,1)}s")


if __name__ == "__main__":
    main()
