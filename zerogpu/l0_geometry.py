"""
l0_geometry.py --- L0 闭链图的**真正谱维数**：文档只算到 T=18，本脚本推到 T=26。

**文档口径（simulations/zero_sum_geometry_probe.py:180-221，逐字）**

    graph_spectrum : N<=1000 用稠密全谱；N>1000 用 eigsh(k=min(300,N-2), which="SM")
                     然后丢掉 <=1e-10 的特征值（即丢掉零模）
    spectral_curve : P(t) = mean(exp(-eigenvalues * t))   <-- 注意分母是 len(eigenvalues)
                     d_s(t) = -2 * d log P / d log t
                     times = logspace(-1, 1.8, 180)

**本脚本指出的两处口径问题（并给出更正后的量）**

  (1) **热核被截断**：真正的回返概率是 P(t) = (1/N) * Σ_{全部 N 个特征值} e^{-λ_k t}。
      文档在 N>1000 时只取最低 300 个特征值，却除以 300。
      低特征值权重最大 => 该估计**系统性高估** P(t) => 斜率被压平 => **低估 d_s**。
      本脚本用**完整谱**（稠密或 Lanczos 求积）算真 P(t)。
      （文档 T<=16 时 N<=810 走稠密，是精确的；**只有 T=18(N=2704) 一行受此影响**。）

  (2) **窗口外的点没有意义**：对任何有限图，t 很大时 P(t) ~ (m1/N)e^{-λ1 t}
      => d_s ~ 2λ1 t **线性发散**。所以 d_s(t) 只在
      [~1/λmax, ~1/λ1] 这个窗口里有意义。文档在 times 到 63 全程取点，
      P=0.02 对应的点可能已落在窗口外。

**本脚本输出**：d_s(t) 曲线、三个文档概率点的 d_s、以及窗口内的平台判定 —— 全部到 T=26。
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh_tridiagonal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from zcl import Engine
import l0_closure as L0
from analyze_l0 import cyc_necklace_graph

SPECTRAL_POINTS = (0.10, 0.05, 0.02)
TIMES = np.logspace(-1.0, 1.8, 180)


# ------------------------------------------------------------------ 建图
def build(T: int, eng: Engine):
    reps, _ = L0.enumerate_necklaces(T, eng, verbose=False)
    n, r, c = cyc_necklace_graph(T, reps, eng)
    E = np.unique(np.sort(np.stack([r, c], 1), 1), axis=0)
    A = sp.coo_matrix((np.ones(len(E)), (E[:, 0], E[:, 1])), shape=(n, n)).tocsr()
    A = A.maximum(A.T)
    deg = np.asarray(A.sum(1)).ravel()
    Lap = (sp.diags(deg) - A).tocsr().astype(np.float64)
    return reps, n, len(E), A, Lap, deg


# -------------------------------------------------------------- 谱：稠密
def full_spectrum_dense(Lap, nmax_dense):
    n = Lap.shape[0]
    if n > nmax_dense:
        return None
    return np.linalg.eigvalsh(Lap.toarray())


# ------------------------------------------- 谱：Lanczos 求积（真 Tr e^{-Lt}/N）
def lanczos_quadrature(Lap, n, times, m=250, probes=48, seed=0):
    """
    P(t) = Tr(e^{-Lt})/N 的 Hutchinson 估计：
        P(t) ~ (1/R) Σ_r  e_1^T exp(-T_m^{(r)} t) e_1 * (||z_r||²/N)
    对 ±1 随机向量 ||z_r||² = N，故归一因子为 1。
    每个探测向量的 Gauss 求积节点/权重由 m 步 Lanczos 三对角化给出。

    **必须做完全重正交化**：不做的话丢正交性会让 beta 失真、三对角矩阵数值不定，
    `eigh_tridiagonal` 会直接不收敛（T=20 上实测 LAPACK info=22）。
    """
    rng = np.random.default_rng(seed)
    acc = np.zeros_like(times)
    used = 0
    for _ in range(probes):
        z = rng.integers(0, 2, size=n).astype(np.float64) * 2.0 - 1.0
        beta0 = np.linalg.norm(z)
        Q = np.zeros((m, n))
        Q[0] = z / beta0
        alpha = np.zeros(m)
        beta = np.zeros(max(m - 1, 1))
        steps = m
        for j in range(m):
            w = Lap @ Q[j]
            a = float(Q[j] @ w)
            alpha[j] = a
            w -= a * Q[j]
            if j > 0:
                w -= beta[j - 1] * Q[j - 1]
            Qj = Q[:j + 1]                       # 两遍完全重正交化
            for _ in range(2):
                w -= Qj.T @ (Qj @ w)
            if j < m - 1:
                b = float(np.linalg.norm(w))
                if b < 1e-10:                    # 真 breakdown：不变子空间已穷尽
                    steps = j + 1
                    break
                beta[j] = b
                Q[j + 1] = w / b
        ev, evec = eigh_tridiagonal(alpha[:steps], beta[:steps - 1])
        w0 = evec[0, :] ** 2
        acc += (w0[:, None] * np.exp(-np.outer(ev, times))).sum(axis=0) * (beta0 ** 2) / n
        used += 1
    acc /= max(used, 1)
    return acc, used


def spectral_dimension(times, P):
    return -2.0 * np.gradient(np.log(np.maximum(P, 1e-300)), np.log(times))


def points_at(times, P, ds, probs=SPECTRAL_POINTS):
    out = {}
    for p in probs:
        i = int(np.argmin(np.abs(P - p)))
        out[f"{p:.2f}"] = {"time": float(times[i]), "spectral_dimension": float(ds[i])}
    return out


# ------------------------------------------------- 邻域增长 / 局部增长维数
def _gather(A, frontier):
    """向量化 CSR 邻接收集：把 frontier 里每个点的邻居拼起来（不用 Python 循环）。"""
    indptr, indices = A.indptr, A.indices
    starts = indptr[frontier]
    counts = indptr[frontier + 1] - starts
    total = int(counts.sum())
    if total == 0:
        return np.zeros(0, dtype=indices.dtype)
    offs = np.repeat(starts - np.concatenate(([0], np.cumsum(counts)[:-1])), counts)
    return indices[offs + np.arange(total)]


def neighbourhood_growth_sampled(A, n, max_radius=20, sources=64, seed=0, cap=400_000):
    """采样式 BFS 球增长（文档对全部源点平均；大图上改为抽样平均）。"""
    rng = np.random.default_rng(seed)
    src = rng.choice(n, size=min(sources, n), replace=False)
    growth = np.zeros(max_radius + 1, dtype=float)
    for s in src:
        visited = np.zeros(n, dtype=bool)
        visited[s] = True
        nvis = 1
        frontier = np.array([s], dtype=np.int64)
        growth[0] += 1.0
        for r in range(1, max_radius + 1):
            if frontier.size == 0 or nvis > cap:
                break
            nbrs = _gather(A, frontier)
            nbrs = np.unique(nbrs[~visited[nbrs]])
            if nbrs.size == 0:
                break
            visited[nbrs] = True
            nvis += nbrs.size
            frontier = nbrs
            growth[r] += nvis
    growth /= len(src)
    return growth


def local_growth_dimension(growth, total_nodes):
    valid = [r for r in range(2, len(growth))
             if growth[r] > growth[r - 1] and growth[r] <= 0.8 * total_nodes]
    if len(valid) < 3:
        return None, None
    radii = np.asarray(valid, dtype=float)
    vals = growth[np.asarray(valid)]
    slope = np.polyfit(np.log(radii), np.log(vals), 1)[0]
    return float(slope), [int(r) for r in valid]


# ---------------------------------------------------------------------- 主
def main(Ts, nmax_dense=12000):
    eng = Engine()
    print("device:", eng.info()["name"], flush=True)
    rows = []
    curves = {}
    t_all = time.time()
    for T in Ts:
        t0 = time.time()
        reps, n, ne, A, Lap, deg = build(T, eng)
        exact = full_spectrum_dense(Lap, nmax_dense)
        if exact is not None:
            lam = exact
            method = "dense-full"
            P = np.array([np.mean(np.exp(-lam * t)) for t in TIMES])
            used = n
        else:
            P, used = lanczos_quadrature(Lap, n, TIMES, m=250, probes=48)
            method = f"lanczos-quad(m=250,R={used})"
            lam = None
        ds = spectral_dimension(TIMES, P)
        pts = points_at(TIMES, P, ds)
        # 窗口内平台：取 P 落在 [0.2, 0.02] 的区间
        win = (P <= 0.20) & (P >= 0.02)
        plateau = (float(np.median(ds[win])), int(win.sum())) if win.sum() >= 3 else (None, 0)
        gr = neighbourhood_growth_sampled(A, n, sources=48 if n > 5000 else min(n, 256))
        lgd, vradii = local_growth_dimension(gr, n)
        row = {
            "T": T, "nodes": int(n), "edges": int(ne),
            "mean_degree": round(float(deg.mean()), 4),
            "method": method,
            "lambda1": None if lam is None else float(lam[1]),
            "lambda_max": None if lam is None else float(lam[-1]),
            "spectral_points": pts,
            "plateau_median_ds_P020_002": None if plateau[0] is None else round(plateau[0], 4),
            "plateau_n_points": plateau[1],
            "local_growth_dimension": None if lgd is None else round(lgd, 4),
            "growth_radii_used": vradii,
            "seconds": round(time.time() - t0, 1),
        }
        rows.append(row)
        curves[T] = {"times": TIMES.tolist(), "P": P.tolist(), "ds": ds.tolist(),
                     "growth": gr.tolist()}
        print(f"T={T:3d} N={n:9,} [{method:22s}] "
              f"d_s@0.05={pts['0.05']['spectral_dimension']:7.4f} "
              f"window-median={row['plateau_median_ds_P020_002']} "
              f"growth-dim={row['local_growth_dimension']} [{row['seconds']}s]", flush=True)

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "l0_geometry.json"), "w") as f:
        json.dump({"device": eng.info(), "rows": rows,
                   "total_seconds": round(time.time() - t_all, 1)}, f, indent=2)

    # ---- 画图 ----
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    cmap = plt.get_cmap("viridis")
    for i, T in enumerate(Ts):
        c = cmap(i / max(len(Ts) - 1, 1))
        cv = curves[T]
        ax[0].loglog(cv["times"], cv["P"], color=c, lw=1.6, label=f"T={T}")
        ax[1].semilogx(cv["times"], cv["ds"], color=c, lw=1.6, label=f"T={T}")
        g = np.array(cv["growth"])
        r = np.arange(len(g))
        ax[2].loglog(r[1:], g[1:], color=c, lw=1.6, marker="o", ms=2.5, label=f"T={T}")
    for p in (0.10, 0.05, 0.02):
        ax[0].axhline(p, color="k", ls=":", lw=0.8)
    ax[0].set_xlabel("t"); ax[0].set_ylabel("P(t) = Tr e^{-Lt}/N")
    ax[0].set_title("真热核（完整谱）"); ax[0].legend(fontsize=7); ax[0].grid(alpha=.25)
    ax[1].axhline(4, color="r", ls="--", lw=1.2, label="d=4")
    ax[1].set_xlabel("t"); ax[1].set_ylabel("$d_s(t)=-2\\,d\\ln P/d\\ln t$")
    ax[1].set_title("谱维数"); ax[1].set_ylim(0, 14); ax[1].legend(fontsize=7); ax[1].grid(alpha=.25)
    ax[2].set_xlabel("半径 r"); ax[2].set_ylabel("球内节点数")
    ax[2].set_title("邻域增长"); ax[2].legend(fontsize=7); ax[2].grid(alpha=.25)
    fig.suptitle("L0 闭链图：几何长出来了吗？（T=8..%d）" % Ts[-1], fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    png = os.path.join(os.path.dirname(os.path.abspath(__file__)), "l0_geometry.png")
    fig.savefig(png, dpi=150)
    print("\nplot:", png)
    print("\n" + "=" * 96)
    print(f"{'T':>3} {'N':>10} {'d_s@0.10':>9} {'d_s@0.05':>9} {'d_s@0.02':>9} "
          f"{'窗口中位':>9} {'增长维数':>9} {'方法':>12}")
    for r in rows:
        sp_ = r["spectral_points"]
        print(f"{r['T']:>3} {r['nodes']:>10,} {sp_['0.10']['spectral_dimension']:>9.4f} "
              f"{sp_['0.05']['spectral_dimension']:>9.4f} "
              f"{sp_['0.02']['spectral_dimension']:>9.4f} "
              f"{str(r['plateau_median_ds_P020_002']):>9} "
              f"{str(r['local_growth_dimension']):>9} {r['method']:>12}")


if __name__ == "__main__":
    Ts = [int(a) for a in sys.argv[1:]] or [8, 10, 12, 14, 16, 18, 20, 22, 24]
    main(Ts)
