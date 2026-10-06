"""
evo_survey.py --- 把自持演化的产出**按层**摆出来看（L0 / L1' / L2）

用户："我们先看看各层产生的东西"。

对每个（基底 × 模式）重跑一遍并**保存场**，然后
  ① 按层打印汇总表（L0 / L1' / L2 / 跨层）
  ② 画图：每个组合两块场 —— **记录层 $R_v$（L1'）** 与 **持续占用场 $f_v$（L2）**

格子类基底直接画二维；图类基底用 **Laplacian 谱布局**（前两个非平凡特征向量当坐标）。

用法：/usr/bin/python3 evo_survey.py
输出：results/evo_survey.json ＋ zerogpu/evo_layers.png
"""
from __future__ import annotations

import json
import os
import time

import numpy as np
from scipy.sparse import csr_matrix, identity, kron
from scipy.sparse.linalg import eigsh

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams["font.sans-serif"] = ["Heiti TC", "Arial Unicode MS", "Hiragino Sans GB"]
matplotlib.rcParams["axes.unicode_minus"] = False

import z0_evo as EV

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)


def run_capture(A, mode, steps, win=500, seed=7):
    """重跑一次，返回 (场, 汇总)。场 = 末窗口的 R 与 f。"""
    N = A.shape[0]
    indptr, indices = A.indptr, A.indices
    deg = np.diff(indptr).astype(np.int64)
    dmax = int(deg.max())
    nb = np.zeros((N, dmax), np.int64)
    mask = np.zeros((N, dmax), bool)
    for v in range(N):
        k = deg[v]
        nb[v, :k] = indices[indptr[v]:indptr[v + 1]]
        mask[v, :k] = True
    rng = np.random.default_rng(seed)
    n_walk = N
    pos = rng.integers(0, N, n_walk).astype(np.int64)
    start = pos.copy()
    R = np.zeros(N, np.int64)
    T = np.zeros((N, dmax), np.int64)
    occ = np.zeros(N, np.int64)
    exc = []
    n_closed = 0
    for t in range(steps):
        m = len(pos)
        nbr = nb[pos]
        ok = mask[pos]
        if mode == "UNIF":
            w = ok.astype(float)
        else:
            w = ok.astype(float).copy()
            if mode in ("REC", "BOTH"):
                w = w * (1.0 + np.minimum(R[nbr], 2))
            if mode in ("TRAF", "BOTH"):
                slot = np.arange(dmax)[None, :]
                w = w * (1.0 + T[pos[:, None], slot])
            w = np.where(ok, w, 0.0)
            bad = w.sum(axis=1) <= 0
            if bad.any():
                w[bad] = ok[bad].astype(float)
        cw = np.cumsum(w, axis=1)
        u = rng.random(m)[:, None] * cw[:, -1:]
        pick = (u > cw).sum(axis=1)
        np.clip(pick, 0, dmax - 1, out=pick)
        np.putmask(pick, ~ok[np.arange(m), pick], 0)
        nxt = nbr[np.arange(m), pick]
        if mode in ("TRAF", "BOTH"):
            T[pos, pick] += 1
        pos = nxt
        cl = np.nonzero(pos == start)[0]
        if len(cl):
            np.add.at(R, pos[cl], 1)
            n_closed += len(cl)
            exc.extend([0])                     # 长度不细究，只数次数
            start[cl] = pos[cl]
        occ[pos] += 1
        if (t + 1) % win == 0:
            occ[:] = 0
    # 再跑一个窗口，拿持续占用场
    for t in range(win):
        m = len(pos)
        nbr = nb[pos]
        ok = mask[pos]
        w = ok.astype(float).copy()
        if mode in ("REC", "BOTH"):
            w = w * (1.0 + np.minimum(R[nbr], 2))
        if mode in ("TRAF", "BOTH"):
            slot = np.arange(dmax)[None, :]
            w = w * (1.0 + T[pos[:, None], slot])
        w = np.where(ok, w, 0.0)
        bad = w.sum(axis=1) <= 0
        if bad.any():
            w[bad] = ok[bad].astype(float)
        cw = np.cumsum(w, axis=1)
        u = rng.random(m)[:, None] * cw[:, -1:]
        pick = (u > cw).sum(axis=1)
        np.clip(pick, 0, dmax - 1, out=pick)
        np.putmask(pick, ~ok[np.arange(m), pick], 0)
        nxt = nbr[np.arange(m), pick]
        if mode in ("TRAF", "BOTH"):
            T[pos, pick] += 1
        pos = nxt
        cl = np.nonzero(pos == start)[0]
        if len(cl):
            np.add.at(R, pos[cl], 1)
            n_closed += len(cl)
            age_reset = True
            start[cl] = pos[cl]
        occ[pos] += 1
    f = occ / win
    keep = f >= 0.5
    lab, c = EV.components(keep, A)
    sz = np.bincount(lab[lab >= 0]) if c else np.array([0])
    summary = {
        "N": int(N), "mode": mode, "steps": steps,
        "L0": {"mean_degree": round(float(deg.mean()), 3),
               "clustering": None, "traffic_gini": round(EV.gini(T.ravel()), 4),
               "closures": int(n_closed)},
        "L1p": {"records_total": int(R.sum()), "record_rate_per_step": round(R.sum() / steps, 4),
                "gini": round(EV.gini(R), 4),
                "max_frac": round(float(R.max() / max(1, R.sum())), 5),
                "support_frac": round(float((R > 0).mean()), 4),
                "entropy": round(float(-((R[R > 0] / R.sum()) * np.log(R[R > 0] / R.sum())).sum()), 4)},
        "L2": {"regions": int(c), "occupied_frac": round(float(keep.mean()), 4),
               "mean_region_size": round(float(sz.mean()), 2), "max_region_size": int(sz.max()),
               "region_size_hist": np.bincount(sz).tolist()[:20] if len(sz) else []},
        "cross": {"corr_R_f": round(float(np.corrcoef(R, f)[0, 1]), 4)
                  if np.std(R) > 0 and np.std(f) > 0 else None},
    }
    return R, f, summary


def layout(A, seed=0):
    """图类基底布局：networkx spring_layout（谱布局在 G_18 上会塌成一条线）。"""
    import networkx as nx
    G = nx.from_scipy_sparse_array(A)
    pos = nx.spring_layout(G, seed=seed, iterations=60)
    x = np.array([pos[i][0] for i in range(A.shape[0])])
    y = np.array([pos[i][1] for i in range(A.shape[0])])
    return x, y


def main():
    t0 = time.time()
    subs = [("torus_L32", EV.torus_adj(2, 32), True),
            ("rand4reg_N1024", EV.random_regular(1024, 4, seed=1), False)]
    try:
        import l0_closure as L0
        import observable_sweep as OS
        from zcl import Engine
        eng = Engine()
        reps, _ = L0.enumerate_necklaces(18, eng, chunk_bits=26, verbose=False)
        subs.append(("G_18", OS.build_sparse(18, reps, eng), False))
    except Exception as e:
        print("  [warn] 闭链图跳过:", e)

    modes = ["REC", "TRAF", "BOTH"]
    res, fields = [], {}
    print("逐层汇总（每个组合重跑并保存场）")
    for sname, A, is_lat in subs:
        steps = 20000 if A.shape[0] <= 1200 else 10000
        for mode in modes:
            R, f, s = run_capture(A, mode, steps)
            s["substrate"] = sname
            res.append(s)
            fields[(sname, mode)] = (R, f, is_lat)
            print(f"\n### {sname} × {mode}   (N={s['N']}, {steps} 步)")
            print(f"  L0  : 平均度={s['L0']['mean_degree']}  边流量Gini={s['L0']['traffic_gini']}  "
                  f"闭合事件={s['L0']['closures']:,}")
            print(f"  L1' : 记录={s['L1p']['records_total']:,}  速率={s['L1p']['record_rate_per_step']}/步  "
                  f"Gini={s['L1p']['gini']}  最大占比={s['L1p']['max_frac']}  "
                  f"覆盖={s['L1p']['support_frac']}  熵={s['L1p']['entropy']}")
            print(f"  L2  : 区域数={s['L2']['regions']}  占用={s['L2']['occupied_frac']}  "
                  f"区域尺寸 均值={s['L2']['mean_region_size']} 最大={s['L2']['max_region_size']}")
            print(f"  跨层: corr(R, f) = {s['cross']['corr_R_f']}")

    np.savez_compressed(os.path.join(OUT, "evo_fields.npz"),
                        **{f"{k[0]}|{k[1]}|{w}": v for k, (R, f, _l) in fields.items()
                           for w, v in (("R", R), ("f", f))})
    with open(os.path.join(OUT, "evo_survey.json"), "w") as fp:
        json.dump({"generated": time.strftime("%Y-%m-%d %H:%M:%S"), "runs": res}, fp,
                  ensure_ascii=False, indent=1)

    # ---------------- 画图 ----------------
    nrow, ncol = len(subs), len(modes)
    fig, axes = plt.subplots(nrow, ncol * 2, figsize=(4.2 * ncol * 2, 4.0 * nrow))
    for i, (sname, A, is_lat) in enumerate(subs):
        for j, mode in enumerate(modes):
            R, f, _ = fields[(sname, mode)]
            ax1, ax2 = axes[i, j * 2], axes[i, j * 2 + 1]
            if is_lat:
                L = int(round(np.sqrt(len(R))))
                ax1.imshow(R.reshape(L, L), cmap="inferno", interpolation="nearest")
                ax2.imshow(f.reshape(L, L), cmap="viridis", interpolation="nearest")
                ax1.set_xlabel("x"); ax1.set_ylabel("y")
            else:
                x, y = layout(A)
                ax1.scatter(x, y, c=R, s=3, cmap="inferno")
                ax2.scatter(x, y, c=f, s=3, cmap="viridis")
            ax1.set_title(f"{sname} × {mode}\nL1' 记录场 $R_v$ (max={R.max()})", fontsize=9)
            ax2.set_title(f"L2 持续占用场 $f_v$  (占用 {np.mean(f>=0.5):.2f}, "
                          f"区域 {res[i*len(modes)+j]['L2']['regions']})", fontsize=9)
            for ax in (ax1, ax2):
                ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("自持演化的产出按层摆开：左=记录层 L1'，右=活动层 L2（无预设目标）", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    png = os.path.join(HERE, "evo_layers.png")
    fig.savefig(png, dpi=110)
    print(f"\n图：{png}   用时 {round(time.time()-t0,1)}s")


if __name__ == "__main__":
    main()
