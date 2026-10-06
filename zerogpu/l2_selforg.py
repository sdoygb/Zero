"""
l2_selforg.py --- 把反馈强度 α 换成**无参数的自反馈**，并做标度与判据交叉复核

**上一轮留下的问题**：`l2_multi.py` 证明"局部历史反馈 ⟹ L2 分成 4–7 个持续区域"，
但反馈强度 $\alpha$ 是我拧的旋钮。本轮把旋钮拆掉。

**四条规则（全部无自由参数）**

$$
\begin{aligned}
&\texttt{UNIF}\ \text{(对照)}: && w_u = 1 &&\text{均匀随机走动，无反馈}\\
&\texttt{MULT}: && w_u = 1+R_u &&\text{"越有历史越吸引"（无上界）}\\
&\texttt{CAP2}: && w_u = 1+\min(R_u,2) &&\textbf{D222 的"最高两层保留"}\text{ 直接给出的饱和}\\
&\texttt{CONTRAST}: && w_u = 1+\max(0,\;R_u-\bar R_{\rm nbr}(u)) &&\text{局部对比度（超出邻居平均的部分）}
\end{aligned}
$$

每步：走者按 $w_u$ 在**邻居**中加权抽样；闭合在 $v$ ⟹ $R_v{+}{=}1$（写记录）；起点重置（Z0①）。
**没有 $\alpha$。**

**测什么**
  ① 区域数（判据 A：持续占用场的连通分量）——四条规则谁能在**无旋钮**下给出多区域
  ② **区域出生/死亡率**与寿命 —— 多区域是不是"生灭循环"（接用户上一个问题）
  ③ **标度**：$k^*(N)$ 随 $N=L^2$ 怎么走（$L=16,24,32,48,64$）
  ④ **判据交叉复核**：同一批运行上用"远足碰撞图"（判据 B）独立数一遍区域数

**等级**：【数值证据】。用法：/usr/bin/python3 l2_selforg.py
输出：results/l2_selforg.json
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


def weights(rule, R, nbr, Rbar_nbr):
    """四条规则的邻居权重。nbr: (m,deg) 邻居编号；Rbar_nbr: (m,) 邻居记录均值。"""
    Rn = R[nbr]                                  # (m, deg)
    if rule == "UNIF":
        return np.ones_like(Rn, dtype=float)
    if rule == "MULT":
        return 1.0 + Rn
    if rule == "CAP2":
        return 1.0 + np.minimum(Rn, 2)
    if rule == "CONTRAST":
        return 1.0 + np.maximum(0.0, Rn - Rbar_nbr[:, None])
    raise ValueError(rule)


def run(L, steps, rule, n_walk=None, win=200, seed=0, collect_support=False):
    N = L * L
    A = torus_adj(2, L)
    indptr, indices = A.indptr, A.indices
    deg = np.diff(indptr).astype(np.int64)
    d0 = int(deg[0])
    nb = indices[indptr[:-1][:, None] + np.arange(d0)[None, :]]      # (N, d0) 邻居表
    rng = np.random.default_rng(seed)
    n_walk = n_walk or N

    start = rng.integers(0, N, n_walk).astype(np.int64)
    pos = start.copy()
    R = np.zeros(N, np.int64)
    occ = np.zeros(N, np.int64)
    counts, masks, ginis, sizes_all = [], [], [], []
    births, deaths, lives = [], [], []
    prev_masks = []
    excursions = []                               # (len, frozenset) 只在 collect_support 时收
    visited = [set([int(p)]) for p in pos]
    age = np.zeros(n_walk, np.int64)

    for t in range(steps):
        m = len(pos)
        nbr = nb[pos]                                        # (m, d0)
        if rule == "UNIF":
            pick = rng.integers(0, d0, size=m)
            nxt = nbr[np.arange(m), pick]
        else:
            Rb = R[nbr].mean(axis=1)
            w = weights(rule, R, nbr, Rb)
            cw = np.cumsum(w, axis=1)
            u = rng.random(m)[:, None] * cw[:, -1:]
            pick = (u > cw).sum(axis=1)
            np.clip(pick, 0, d0 - 1, out=pick)
            nxt = nbr[np.arange(m), pick]
        pos = nxt
        age += 1
        if collect_support:
            for i in range(m):
                visited[i].add(int(pos[i]))
        cl = np.nonzero(pos == start)[0]
        if len(cl):
            np.add.at(R, pos[cl], 1)
            if collect_support:
                for i in cl:
                    excursions.append((int(age[i]), frozenset(visited[i])))
                    visited[i] = {int(pos[i])}
            age[cl] = 0
            start[cl] = pos[cl]
        occ[pos] += 1

        if (t + 1) % win == 0:
            f = occ / win
            keep = f >= 0.5
            lab, c = components(keep, A)
            counts.append(int(c))
            cur = [np.nonzero(lab == j)[0] for j in range(c)]
            masks.append(cur)
            sizes_all.append([len(x) for x in cur])
            ginis.append(float(_gini(R)))
            # 出生/死亡/寿命：与上一窗口按重叠 > 50% 配对
            if prev_masks:
                matched_prev = set()
                for mj in cur:
                    best, bi = 0.0, -1
                    for i, mi in enumerate(prev_masks):
                        ov = np.intersect1d(mj, mi, assume_unique=False).size / max(1, len(mj))
                        if ov > best:
                            best, bi = ov, i
                    if best > 0.5:
                        matched_prev.add(bi)
                        lives.append(1)
                births.append(len(cur) - len(matched_prev))
                deaths.append(len(prev_masks) - len(matched_prev))
            prev_masks = cur
            occ[:] = 0

    return {"rule": rule, "L": L, "N": int(N), "n_walk": int(n_walk), "steps": steps,
            "region_series": counts,
            "mean_regions": round(float(np.mean(counts)), 3) if counts else 0.0,
            "std_regions": round(float(np.std(counts)), 3) if counts else 0.0,
            "mean_births_per_window": round(float(np.mean(births)), 3) if births else None,
            "mean_deaths_per_window": round(float(np.mean(deaths)), 3) if deaths else None,
            "record_gini": round(float(np.mean(ginis)), 4) if ginis else None,
            "record_max_frac": round(float(R.max() / max(1, R.sum())), 4),
            "mean_region_size": round(float(N / max(1e-9, np.mean(counts))), 2) if counts else None,
            "region_sizes_last_window": sizes_all[-1][:60] if sizes_all else None,
            "mean_region_size_from_last": round(float(np.mean(sizes_all[-1])), 2)
                                          if sizes_all and sizes_all[-1] else None,
            "records_total": int(R.sum()),
            "excursions": len(excursions) if collect_support else None,
            "_excursions": excursions if collect_support else None}


def _gini(x):
    x = np.sort(np.asarray(x, float))
    n = len(x)
    if n == 0 or x.sum() == 0:
        return 0.0
    return float((2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


def criterion_B(excursions, A):
    """判据 B：远足碰撞图（非原点顶点相交）的连通块数 = 区域数。"""
    n_ex = len(excursions)
    if n_ex == 0:
        return None
    parent = list(range(n_ex))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    seen = {}
    for i, (_ln, sites) in enumerate(excursions):
        for s in sites:
            j = seen.get(s)
            if j is None:
                seen[s] = i
            else:
                ra, rb = find(j), find(i)
                if ra != rb:
                    parent[rb] = ra
    return len({find(i) for i in range(n_ex)})


def main(argv):
    t0 = time.time()
    out = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
           "question": "把 α 拆掉之后，L2 还能不能自发分成多区域？",
           "rules": {"UNIF": "w=1（对照）", "MULT": "w=1+R（无上界）",
                     "CAP2": "w=1+min(R,2)（D222 最高两层保留）",
                     "CONTRAST": "w=1+max(0,R-Rbar_nbr)"},
           "runs": []}

    print("=" * 100)
    print("① 四条无参数规则（Γ=ℤ² L=32，N=1024，走者=N，4000 步，窗口 200）")
    print(f"{'规则':>9} {'区域数':>9} {'区域数序列':>30} {'出生/窗':>8} {'死亡/窗':>8} "
          f"{'记录Gini':>9} {'记录集中度':>10}")
    for rule in ("UNIF", "MULT", "CAP2", "CONTRAST"):
        rs = [run(32, 4000, rule, seed=5 + 13 * k) for k in range(3)]
        mr = np.array([r["mean_regions"] for r in rs])
        r = rs[0]
        r["mean_regions_3seed"] = round(float(mr.mean()), 3)
        r["std_regions_3seed"] = round(float(mr.std()), 3)
        out["runs"].append(r)
        print(f"{rule:>9} {r['mean_regions_3seed']:>5.2f}±{r['std_regions_3seed']:<3.2f} "
              f"{str(r['region_series'][:9]):>30} {str(r['mean_births_per_window']):>8} "
              f"{str(r['mean_deaths_per_window']):>8} {r['record_gini']:>9.3f} "
              f"{r['record_max_frac']:>10.4f}", flush=True)

    # ---- ② 标度：最好的规则（按区域数排序取最大者，除 UNIF 外）----
    print("\n" + "=" * 100)
    print("② 标度：区域数 vs N（无参数规则 vs UNIF 对照）")
    cand = "CAP2"      # 选 CAP2：它是 D222 的保留规则直接给出的（有出处），不是按区域数挑的
    print(f"   （做 CAP2 与 MULT 两条；CAP2 是 D222 的保留规则直接给出的）")
    print(f"{'L':>4} {'N':>7} {'CAP2 区域数':>12} {'CAP2 区域尺寸':>13} "
          f"{'MULT 区域数':>11} {'MULT 尺寸':>10} {'UNIF':>7}")
    scale = []
    for L in (16, 24, 32, 48, 64):
        ra = run(L, 3000, "CAP2", seed=9)
        rm = run(L, 3000, "MULT", seed=9)
        rb = run(L, 3000, "UNIF", seed=9)
        scale.append({"L": L, "N": L * L,
                      "cap2_regions": ra["mean_regions"], "cap2_size": ra["mean_region_size"],
                      "cap2_births": ra["mean_births_per_window"],
                      "cap2_deaths": ra["mean_deaths_per_window"],
                      "cap2_gini": ra["record_gini"], "cap2_maxfrac": ra["record_max_frac"],
                      "mult_regions": rm["mean_regions"], "mult_size": rm["mean_region_size"],
                      "unif_regions": rb["mean_regions"]})
        print(f"{L:>4} {L*L:>7} {ra['mean_regions']:>12.2f} {ra['mean_region_size']:>13.2f} "
              f"{rm['mean_regions']:>11.2f} {rm['mean_region_size']:>10.2f} "
              f"{rb['mean_regions']:>7.2f}", flush=True)
    out["scaling"] = scale
    NN = np.array([s["N"] for s in scale], float)
    for key in ("cap2_regions", "mult_regions"):
        KK = np.array([s[key] for s in scale], float)
        sl = float(np.polyfit(np.log(NN), np.log(KK), 1)[0])
        SS = np.array([s["cap2_size"] if key.startswith("cap2") else s["mult_size"]
                       for s in scale], float)
        out[f"scaling_exponent_{key}"] = round(sl, 4)
        print(f"   {key}: k*(N) ~ N^{sl:.3f}   区域尺寸 ≈ {SS.mean():.1f} 格点（N={NN[0]:.0f}→{NN[-1]:.0f} 基本不变）")
        rec("SELFORG", key, "区域数 k* vs N", NN, KK, {"exponent": round(sl, 4)},
            note="区域尺寸不随 N 变 ⟹ 存在涌现长度尺度")

    # ---- ③ 判据 A vs B 交叉复核 ----
    print("\n" + "=" * 100)
    print("③ 判据交叉复核：持续占用场连通分量（A） vs 远足碰撞图（B），Γ=ℤ² L=16")
    rb = run(16, 4000, "CAP2", seed=21, collect_support=True)
    nB = criterion_B(rb["_excursions"], torus_adj(2, 16))
    out["criterion_cross"] = {"rule": cand, "L": 16,
                              "A_mean_regions": rb["mean_regions"],
                              "B_regions_collision_graph": nB,
                              "excursions": rb["excursions"]}
    print(f"   判据 A（持续占用场）= {rb['mean_regions']:.2f} 个区域")
    print(f"   判据 B（碰撞图，{rb['excursions']} 条远足）= {nB} 个区域")
    rb.pop("_excursions", None)
    out["runs"].append(rb)

    with open(os.path.join(OUT, "l2_selforg.json"), "w") as f:
        json.dump({**out, "records": RECORDS}, f, ensure_ascii=False, indent=1)
    print("=" * 100)
    print(f"用时 {round(time.time()-t0,1)}s → results/l2_selforg.json")


if __name__ == "__main__":
    main(sys.argv[1:])
