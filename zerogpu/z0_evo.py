"""
z0_evo.py --- **自持演化引擎 + 自动观测台**（去掉全部预设目标）

**用户的要求**：把"4 维/粒子/光锥"这类预选设计全部去掉。让它自己跑，产生垃圾也好黄金也好，
我们只要**支持演化的架构**；能生成就生成，不能生成就**分析数据**。

**架构（三个模式，全部无自由参数）**
    L0   : Γ（给定基底）＋ 走者 ＋ 整数重数
    L1'  : 闭合 $\sum w=0$ $\Rightarrow$ 写记录 $R_v{+}{=}1$
    回灌 : 走者按邻居权重加权抽样；起点重置（Z0① 闭合不停）
    自持 : 走者数**守恒**（每条走者每步走一条边）⟹ 精确自持（`LOOP.md` 已证 $r^*{=}1$）

| 模式 | 邻居权重 $w_{u\to v}$ | 读什么 |
|:--|:--|:--|
| `A: REC`  | $1+\min(R_v,\,2)$ | 记录层（`D222` 最高两层保留 ⟹ 饱和常数 2） |
| `B: TRAF` | $1+T_{uv}$ | **走过这条边的次数**（Hebbian：走过的路更常走）⟹ **图结构自己演化** |
| `C: BOTH` | $(1+\min(R_v,2))\,(1+T_{uv})$ | 两者相乘 |

**观测台（不挑、不预设）**：每个窗口记录 ~24 条时间序列，全部交给 `structures.classify_curve`
**自动做结构分类**（幂律/指数/对数/饱和/delta/周期/简并），并**自动检测突变**（相邻窗口的相对跳变、均值变点）。
输出是"**长出了什么**"的目录，不是"有没有我要的东西"。

**基底（不选优，三个一起报）**：ℤ² 环面、随机 4-正则图、原生闭链图 $\mathcal G_T$。

**等级**：【数值证据】。用法：/usr/bin/python3 z0_evo.py
输出：results/z0_evo.json
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


# ------------------------------------------------------------------ 基底
def torus_adj(D, L):
    r = np.arange(L)
    A1 = csr_matrix((np.ones(L), (r, (r + 1) % L)), shape=(L, L))
    A1 = A1 + A1.T
    A = A1
    for _ in range(D - 1):
        A = kron(A, identity(L, format="csr")) + kron(identity(A.shape[0], format="csr"), A1)
    return csr_matrix(A)


def random_regular(N, d, seed=0):
    """配置模型 + 去自环/重边（近似 d-正则）。"""
    rng = np.random.default_rng(seed)
    while True:
        stubs = np.repeat(np.arange(N), d)
        rng.shuffle(stubs)
        a, b = stubs[0::2], stubs[1::2]
        ok = a != b
        a, b = a[ok], b[ok]
        A = csr_matrix((np.ones(len(a)), (a, b)), shape=(N, N))
        A = A + A.T
        A.data[:] = 1.0
        if A.nnz >= N * d * 0.9:
            return A


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


def gini(x):
    x = np.sort(np.asarray(x, float))
    n = len(x)
    if n == 0 or x.sum() == 0:
        return 0.0
    return float((2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


# ------------------------------------------------------------------ 引擎
def evolve(A, mode, steps=40000, win=500, seed=0, n_walk=None):
    N = A.shape[0]
    indptr, indices = A.indptr, A.indices
    deg = np.diff(indptr).astype(np.int64)
    dmax = int(deg.max())
    # 邻居表（补齐到 dmax，重复项权重为 0 无效）
    nb = np.zeros((N, dmax), np.int64)
    mask = np.zeros((N, dmax), bool)
    for v in range(N):
        k = deg[v]
        nb[v, :k] = indices[indptr[v]:indptr[v + 1]]
        mask[v, :k] = True
    rng = np.random.default_rng(seed)
    n_walk = n_walk or N
    pos = rng.integers(0, N, n_walk).astype(np.int64)
    start = pos.copy()
    age = np.zeros(n_walk, np.int64)
    R = np.zeros(N, np.int64)                       # 记录层
    T = np.zeros((N, dmax), np.int64)               # 边流量（按邻居槽位存）
    occ = np.zeros(N, np.int64)
    exc_len_recent = deque(maxlen=4000)
    series = {k: [] for k in (
        "t", "occupied_frac", "empty_frac", "regions", "mean_region_size", "max_region_size",
        "births", "deaths", "records_total", "record_rate", "record_gini", "record_maxfrac",
        "record_entropy", "closed_total", "excursion_mean_len", "excursion_std_len",
        "ac1", "ac4", "traffic_gini", "hot_sites_frac", "mean_degree", "clustering",
        "rec_occ_corr", "largest_region_ever")}
    prev_masks = []
    n_closed = 0

    # 辅助：边流量索引（无向：两个方向分别记）
    for t in range(steps):
        m = len(pos)
        nbr = nb[pos]                                  # (m, dmax)
        ok = mask[pos]
        if mode == "UNIF":
            w = ok.astype(float)
        else:
            w = ok.astype(float).copy()
            if mode in ("REC", "BOTH"):
                w = w * (1.0 + np.minimum(R[nbr], 2))
            if mode in ("TRAF", "BOTH"):
                # 该走者到各邻居槽位的流量：用 (pos, slot) 索引
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
        age += 1
        cl = np.nonzero(pos == start)[0]
        if len(cl):
            np.add.at(R, pos[cl], 1)
            n_closed += len(cl)
            exc_len_recent.extend(age[cl].tolist())
            age[cl] = 0
            start[cl] = pos[cl]
        occ[pos] += 1

        if (t + 1) % win == 0:
            f = occ / win
            keep = f >= 0.5
            lab, c = components(keep, A)
            sz = np.bincount(lab[lab >= 0]) if c else np.array([0])
            cur = [np.nonzero(lab == j)[0] for j in range(c)]
            births = deaths = 0
            if prev_masks:
                matched = set()
                for mj in cur:
                    best = 0.0
                    bi = -1
                    for i, mi in enumerate(prev_masks):
                        ov = np.intersect1d(mj, mi, assume_unique=False).size / max(1, len(mj))
                        if ov > best:
                            best, bi = ov, i
                    if best > 0.5:
                        matched.add(bi)
                births = len(cur) - len(matched)
                deaths = len(prev_masks) - len(matched)
            prev_masks = cur
            occ[:] = 0
            # 空间自相关（只对格点类基底有意义；环面上用 (0,±1) 邻居近似）
            fr = f
            ac1 = float(np.corrcoef(fr, fr[nb[:, 0]])[0, 1]) if np.std(fr) > 0 else 0.0
            ac4 = float(np.corrcoef(fr, fr[nb[:, 3 % dmax]])[0, 1]) if np.std(fr) > 0 else 0.0
            el = np.array(exc_len_recent, float) if exc_len_recent else np.array([0.0])
            hot = np.sort(R)[::-1][: max(1, N // 100)]
            s = series
            s["t"].append((t + 1) / 1000.0)
            s["occupied_frac"].append(float(keep.mean()))
            s["empty_frac"].append(float(1 - keep.mean()))
            s["regions"].append(int(c))
            s["mean_region_size"].append(float(sz.mean()))
            s["max_region_size"].append(float(sz.max()))
            s["births"].append(int(births))
            s["deaths"].append(int(deaths))
            s["records_total"].append(int(R.sum()))
            s["record_rate"].append(int(R.sum()) / (t + 1))
            s["record_gini"].append(gini(R))
            s["record_maxfrac"].append(float(R.max() / max(1, R.sum())))
            p = R[R > 0] / max(1, R.sum())
            s["record_entropy"].append(float(-(p * np.log(p)).sum()) if len(p) else 0.0)
            s["closed_total"].append(int(n_closed))
            s["excursion_mean_len"].append(float(el.mean()))
            s["excursion_std_len"].append(float(el.std()))
            s["ac1"].append(ac1)
            s["ac4"].append(ac4)
            s["traffic_gini"].append(gini(T.ravel()))
            s["hot_sites_frac"].append(float(hot.sum() / max(1, R.sum())))
            deg_now = np.diff(A.indptr).astype(float)
            s["mean_degree"].append(float(deg_now.mean()))
            # 局部聚类系数（近似：抽 200 个点）
            samp = rng.choice(N, size=min(200, N), replace=False)
            cl_coef = []
            for v in samp:
                nbs = indices[indptr[v]:indptr[v + 1]]
                if len(nbs) < 2:
                    continue
                sub = A[nbs][:, nbs].toarray()
                cl_coef.append(sub.sum() / (len(nbs) * (len(nbs) - 1)))
            s["clustering"].append(float(np.mean(cl_coef)) if cl_coef else 0.0)
            s["rec_occ_corr"].append(float(np.corrcoef(R, f)[0, 1]) if np.std(R) > 0 and np.std(f) > 0 else 0.0)
            s["largest_region_ever"].append(float(max(s["max_region_size"])))
    return series, {"mode": mode, "N": int(N), "steps": steps, "win": win,
                    "n_walk": int(n_walk), "records": int(R.sum()), "closures": int(n_closed)}


def main(argv):
    t0 = time.time()
    out = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
           "architecture": "走者数守恒的自持环 + 无参数邻居权重（REC / TRAF / BOTH）",
           "target_removed": "没有维数扫描、没有缺陷势、没有粒子/光锥指标；只做中性普查",
           "runs": [], "catalogue": []}

    substrates = [("torus_L32", torus_adj(2, 32)),
                  ("rand4reg_N1024", random_regular(1024, 4, seed=1))]
    try:
        import l0_closure as L0
        import observable_sweep as OS
        from zcl import Engine
        eng = Engine()
        reps, _ = L0.enumerate_necklaces(18, eng, chunk_bits=26, verbose=False)
        substrates.append(("G_18", OS.build_sparse(18, reps, eng)))
    except Exception as e:
        print("  [warn] 闭链图跳过:", e)

    print("=" * 104)
    print("自持演化：三个模式 × 三个基底（无预设目标）")
    for sname, A in substrates:
        N = A.shape[0]
        for mode in ("A:REC", "B:TRAF", "C:BOTH"):
            m = mode.split(":")[1]
            steps = 40000 if N <= 1200 else 20000
            ser, info = evolve(A, m, steps=steps, win=500, seed=7)
            info["substrate"] = sname
            out["runs"].append({"info": info, "series": ser})
            print(f"\n### {sname}  N={N}  模式 {mode}   记录 {info['records']:,}  闭合 {info['closures']:,}",
                  flush=True)
            key = ["regions", "occupied_frac", "record_gini", "record_maxfrac",
                   "births", "deaths", "excursion_mean_len", "ac1", "clustering", "traffic_gini"]
            print("   " + "  ".join(f"{k}={ser[k][-1]:.3g}" for k in key))

    # ---------------- 自动普查：对每条时间序列做结构分类 ----------------
    print("\n" + "=" * 104)
    print("自动普查（结构由 structures.classify_curve 判定，无人工挑选）")
    import structures as ST
    from collections import Counter
    cnt = Counter()
    for run in out["runs"]:
        info, ser = run["info"], run["series"]
        for k, y in ser.items():
            if k == "t" or len(y) < 8:
                continue
            st = ST.classify_curve(np.array(ser["t"]), np.array(y, float))
            yv = np.array(y, float)
            jumps = np.abs(np.diff(yv)) / np.maximum(np.abs(yv[:-1]), 1e-12)
            run.setdefault("_struct", []).append({
                "obs": k, "label": st["label"], "r2": st.get("r2"),
                "params": st.get("params"),
                "max_rel_jump": round(float(jumps.max()), 4) if len(jumps) else None,
                "jump_at_window": int(np.argmax(jumps)) if len(jumps) else None,
                "first": round(float(yv[0]), 6), "last": round(float(yv[-1]), 6)})
            cnt[st["label"]] += 1
            out["catalogue"].append({"substrate": info["substrate"], "mode": info["mode"],
                                     "obs": k, "label": st["label"],
                                     "params": st.get("params"), "r2": st.get("r2")})
    print("\n结构聚类（全部 24 条 × 9 次运行）：")
    for k, v in cnt.most_common():
        print(f"   {k:<28} {v}")

    # ---------------- 突变检测 ----------------
    print("\n突变（相邻窗口相对跳变 > 50% 的观测）：")
    ev = 0
    for run in out["runs"]:
        for row in run.get("_struct", []):
            if row["max_rel_jump"] and row["max_rel_jump"] > 0.5:
                ev += 1
                if ev <= 14:
                    print(f"   {run['info']['substrate']:<16} {run['info']['mode']:<6} "
                          f"{row['obs']:<20} 跳变 {row['max_rel_jump']:.2f} "
                          f"@窗口{row['jump_at_window']}  {row['first']:.4g}→{row['last']:.4g}")
    print(f"   （共 {ev} 条）")
    out["events"] = ev

    with open(os.path.join(OUT, "z0_evo.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("=" * 104)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_evo.json")


if __name__ == "__main__":
    main(sys.argv[1:])
