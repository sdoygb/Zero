"""
l2_multi.py --- 让 L2 **自发演化成多区域**（用户的第一优先）

**诊断（来自本会话的实测）**：L2 混成一坨，是因为我用的分区判据——"两条远足曾经在同一个
顶点相撞 ⟹ 同一区域"——在均匀 $\Gamma$、$D\ge2$ 时**必然**把全部远足焊成一个连通块
（`L2_REGIONS.md` §2：$D{=}1$ 得 2 个区域，$D\ge2$ 得 1 个）。所以缺的不是"更多远足"，而是
**一个局部的正反馈**：让某处的活跃度能喂养它自己的后续活跃度。没有它，均匀态是唯一吸引子。

**语料里的出处（不是我造的）**：`D222` §第 8 步——
  「**新活动层只能由保留的历史顶层与零层播种**」，且「只保留全局零层会让区域差异消失，
    因此**局部保留层不可省**」。
把它写成最小形式：

$$
\text{闭合在 }v\ \Rightarrow\ R_v \mathrel{+}= 1\ (\text{写记录，L1'})\qquad
\text{重播种：随机挑一条走者，搬到 }\ v\sim R_v^{\gamma}\ (\text{读局部历史})
$$

$\gamma$ 的强度用一个旋钮表达：$\alpha=$ 每步做多少次"重播种搬运"（$\alpha{=}0$ ⟹ 退化成 `z0_loop` 基线）。

**区域判据（两条，互相独立）**
  A **活动集连通块**：某一时刻被"未闭合走者"占据的格点，在 $\Gamma$ 上的连通分量 = 区域
  B **时间窗碰撞图**：窗口内相撞过的远足聚类（对照"曾经"判据）

**测什么**
  ① 相图：平均区域数、最大区域占比 vs $\alpha$（均匀 ⟶ 碎裂 ⟶ 凝聚）
  ② 区域的**寿命**（跨时间片用 Jaccard 重叠配对）——区域是不是**持续**的，还是噪声
  ③ 记录层 $R_v$ 的**局域化**（Gini 系数）——L1′ 有没有跟着分区
  ④ 区域大小分布的结构（幂律？单一巨块？）

**等级**：【数值证据】。$\Gamma$、$\alpha$、$\gamma$ 是输入；判据 A/B 的差别要分开报。
用法：/usr/bin/python3 l2_multi.py
输出：results/l2_multi.json
"""
from __future__ import annotations

import json
import os
import sys
import time
from collections import deque

import numpy as np
from scipy.sparse import csr_matrix, identity, kron

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

RECORDS = []


def rec(sec, obj, readout, x, y, meta=None, note=None, kind="curve", maxpts=64):
    x = np.asarray(x, float).ravel()
    y = np.asarray(y, float).ravel()
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) > maxpts:
        idx = np.unique(np.r_[np.linspace(0, len(x) - 1, maxpts - 1).astype(int), len(x) - 1])
        x, y = x[idx], y[idx]
    RECORDS.append({"id": f"{sec}:{obj}:{readout}", "sector": sec, "object": obj,
                    "readout": readout, "kind": kind,
                    "x": [float(v) for v in x], "y": [float(v) for v in y],
                    "meta": meta or {}, **({"note": note} if note else {})})


def torus_adj(D, L):
    r = np.arange(L)
    A1 = csr_matrix((np.ones(L), (r, (r + 1) % L)), shape=(L, L))
    A1 = A1 + A1.T
    A = A1
    for _ in range(D - 1):
        A = kron(A, identity(L, format="csr")) + kron(identity(A.shape[0], format="csr"), A1)
    return csr_matrix(A)


def components(occ, A):
    """占用格点在 Γ 上的连通分量（只走占用点）。返回标签数组（-1 = 未占用）。"""
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
    n = len(x)
    if n == 0 or x.sum() == 0:
        return 0.0
    return float((2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


def run(A, n_walk, steps, alpha, win=200, dens_thresh=0.5, gamma=1.0,
        measure_every=100, seed=0):
    """
    保守走者 + 记录层 + **局部**趋性重播种（读邻域历史，D222 第 8 步的最小局部形式）：
        每步，以概率 alpha 走向**邻域中记录数最高**的邻居；否则走向均匀随机邻居。
        闭合 → R_v += 1（写记录）；起点重置（Z0① 继续走）。
    区域判据 A：**持续占用场** f_v（长度 win 的窗口内被占的时间比例）超过 dens_thresh 的
                格点集合，在 Γ 上的连通分量。基线（alpha=0）应当本来就是一整块。
    """
    N = A.shape[0]
    indptr, indices = A.indptr, A.indices
    deg = np.diff(indptr).astype(np.int64)
    rng = np.random.default_rng(seed)

    start = rng.integers(0, N, n_walk).astype(np.int64)
    pos = start.copy()
    R = np.zeros(N, np.int64)
    occ = np.zeros(N, np.int64)
    counts, biggest, empties = [], [], []
    R_hist, masks = [], []

    for t in range(steps):
        m = len(pos)
        # --- 局部趋性移动 ---
        if alpha > 0:
            biased = rng.random(m) < alpha
            nb = np.nonzero(biased)[0]
            if len(nb):
                cnt = deg[pos[nb]]
                st = np.repeat(indptr[pos[nb]], cnt)
                wid = np.repeat(np.arange(len(nb)), cnt)
                nbr = indices[st]
                rv = R[nbr]
                order = np.lexsort((-rv, wid))
                w_sorted = wid[order]
                first = np.ones(len(order), bool)
                first[1:] = w_sorted[1:] != w_sorted[:-1]
                pos[nb] = nbr[order][first]
            un = np.nonzero(~biased)[0]
            if len(un):
                off = (rng.random(len(un)) * deg[pos[un]]).astype(np.int64)
                pos[un] = indices[indptr[pos[un]] + off]
        else:
            off = (rng.random(m) * deg[pos]).astype(np.int64)
            pos = indices[indptr[pos] + off]
        # --- 闭合 → 写记录 ---
        cl = np.nonzero(pos == start)[0]
        if len(cl):
            np.add.at(R, pos[cl], 1)
            start[cl] = pos[cl]
        occ[pos] += 1
        # --- 每 win 步：算持续占用场的连通分量 ---
        if (t + 1) % win == 0:
            f = occ / win
            keep = f >= dens_thresh
            lab, c = components(keep, A)
            counts.append(int(c))
            sz = np.bincount(lab[lab >= 0]) if c else np.array([0])
            biggest.append(float(sz.max() / max(1, sz.sum())))
            empties.append(float(1.0 - keep.mean()))
            R_hist.append(R.copy())
            masks.append([np.nonzero(lab == j)[0] for j in range(c)])
            occ[:] = 0
    # ---- 区域寿命：相邻窗口按"重叠 > 50%"配对，数连续存活多少个窗口 ----
    lives, persist = [], []
    prev = {}
    for w in range(1, len(masks)):
        cur = {}
        matched = 0
        for j, mj in enumerate(masks[w]):
            best, bi = 0.0, -1
            for i, mi in enumerate(masks[w - 1]):
                inter = np.intersect1d(mj, mi, assume_unique=False).size
                ov = inter / max(1, len(mj))
                if ov > best:
                    best, bi = ov, i
            if best > 0.5:
                matched += 1
                cur[j] = prev.get(bi, 1) + 1
                lives.append(cur[j])
        persist.append(matched / max(1, len(masks[w])))
        prev = cur
    return {"alpha": alpha, "n_walk": int(n_walk), "steps": steps, "win": win,
            "region_persistence": round(float(np.mean(persist)), 4) if persist else None,
            "mean_region_lifetime_windows": round(float(np.mean(lives)), 3) if lives else None,
            "max_region_lifetime_windows": int(np.max(lives)) if lives else None,
            "dens_thresh": dens_thresh,
            "mean_regions": float(np.mean(counts)) if counts else 0.0,
            "regions_series": counts,
            "mean_largest_frac": round(float(np.mean(biggest)), 4) if biggest else None,
            "mean_empty_frac": round(float(np.mean(empties)), 4) if empties else None,
            "record_gini": round(gini(R_hist[-1]), 4) if R_hist else None,
            "record_max_frac": round(float(R_hist[-1].max() / max(1, R_hist[-1].sum())), 4)
                                if R_hist else None,
            "records_total": int(R.sum())}


def main(argv):
    t0 = time.time()
    out = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
           "mechanism": "D222 第 8 步：重播种读局部保留历史 ⟹ r_v ∝ R_v^gamma，强度 alpha",
           "criteria": {"A": "活动集在 Γ 上的连通分量", "B": "时间窗碰撞图（对照，见 L2_REGIONS.md）"},
           "runs": []}

    graphs = [("torus_D2_L16", torus_adj(2, 16)),
              ("torus_D2_L32", torus_adj(2, 32))]
    try:
        import l0_closure as L0
        import observable_sweep as OS
        from zcl import Engine
        eng = Engine()
        for T in (14, 16):
            reps, _ = L0.enumerate_necklaces(T, eng, chunk_bits=26, verbose=False)
            graphs.append((f"G_{T}", OS.build_sparse(T, reps, eng)))
    except Exception as e:                     # GPU 不可用时跳过
        print("  [warn] 闭链图跳过:", e)

    alphas = [0.0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.9]
    print("=" * 100)
    for gname, A in graphs:
        N = A.shape[0]
        n_walk = N          # **稠密**：基线必须本来就是一整块，碎裂才有意义
        print(f"\n### Γ = {gname}（N={N}，走者 {n_walk}）")
        print(f"{'α':>6} {'区域数(3种子)':>13} {'区域数序列':>24} "
              f"{'空置率':>8} {'记录集中度':>10} {'持续率':>8} {'寿命(窗)':>8}")
        rows = []
        for a in alphas:
            reps = [run(A, n_walk, 4000, alpha=a, seed=3 + 17 * k) for k in range(3)]
            mr = np.array([x["mean_regions"] for x in reps])
            r = reps[0]
            r["graph"] = gname
            r["mean_regions_mean3"] = round(float(mr.mean()), 3)
            r["mean_regions_std3"] = round(float(mr.std()), 3)
            r["replicates"] = [{k: x[k] for k in ("mean_regions", "region_persistence",
                                                  "mean_region_lifetime_windows",
                                                  "max_region_lifetime_windows",
                                                  "mean_empty_frac", "record_max_frac")}
                               for x in reps]
            rows.append(r)
            print(f"{a:>6} {r['mean_regions_mean3']:>7.2f}±{r['mean_regions_std3']:<4.2f} "
                  f"{str(r['regions_series'][:7]):>24} "
                  f"{r['mean_empty_frac']:>8.3f} {r['record_max_frac']:>10.3f} "
                  f"{str(r['region_persistence']):>8} {str(r['mean_region_lifetime_windows']):>8}")
            out["runs"].append(r)
        rec("MULTI", gname, "平均区域数 vs 反馈强度 α", [r["alpha"] for r in rows],
            [r["mean_regions"] for r in rows], {"N": N, "n_walk": n_walk})
        rec("MULTI", gname, "最大区域占比 vs α", [r["alpha"] for r in rows],
            [r["mean_largest_frac"] for r in rows], {"N": N})
        rec("MULTI", gname, "记录层 Gini vs α", [r["alpha"] for r in rows],
            [r["record_gini"] for r in rows], {"N": N})

    with open(os.path.join(OUT, "l2_multi.json"), "w") as f:
        json.dump({**out, "records": RECORDS}, f, ensure_ascii=False, indent=1)
    print("=" * 100)
    print(f"用时 {round(time.time()-t0,1)}s → results/l2_multi.json")
    print("读法：区域数 ↑ 且持续率高 ⟹ L2 自发分成多区域；最大区域占比 → 1 ⟹ 又凝成一坨。")


if __name__ == "__main__":
    main(sys.argv[1:])
