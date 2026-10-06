"""
l2_regions.py --- 零乱动 → 自发分域 → 多个 L2（"我们的世界"是其中一个？）

**要回答的事**：不人工分区，让最小系统自己演化，看它**自发**分成几个区域层
（L2 是区域性的），再看有没有哪个区域长得像我们这个世界（尤其：**维数**）。

**"零乱动"是什么**
    Z0①（零从不停）＋ Z1 的图 Γ（有限连通、相邻两方向可通行）＋ Z2（每一步**所有可用边**
    都被实例化；同一词按整数计重数） ⟹ 词沿 Γ 走动：n_{k+1} = A n_k。
    分支数是 deg(v) 而不是 2（见 AUDIT_missing_graph.md）。

**分区凭什么不是人工的**
    **亚稳态／几乎不变集**：随机走动 P = D^{-1}A 的慢模。λ_2,λ_3,… 成群贴着 1，
    说明系统自己划出了几个"很久才漏出去"的区域。**区域数由谱隙决定，不是手选的。**
    区域 = 慢模谱嵌入上的 k-means 聚类（谱聚类）。

**每个区域自己的物理**
    · d_s  内部谱维数：内部 Laplacian 的热核 P(t)=⟨e^{-tμ}⟩ ⟹ d_s = -2 dlnP/dlnt
          窗口取 [3/μ_max, 0.3/μ_min]（μ_min/μ_max 用 Lanczos 的 Ritz 值估）
    · d_g  内部生长维数：vol(r) ~ r^{d_g}
    · 泄漏率、内部度、内部谱隙
    · **动力学核验**：真实随机轨迹在区域内的停留时间（应指数分布，mean≈std）
      ＋ 区域间转移矩阵 —— 区域是**动力学**对象，不只是谱的产物

**估计器标定**（诚实声明）：d_s 估计器在**已知答案**的 ℤ^D 环面上标定过，
误差见输出里的 `calibration`（D=1,2 约 ±5%，D=3,4 偏差可达 ±20%，且窗口太窄时
不可用 —— 每个区域都带 `quality` 字段，窗口 <1 个数量级的结果不要当真）。

**图的选择（可审计）**
    Γ_native = 𝒢_T（闭合类图：平衡项链 ＋ 循环相邻对换）：**从 Z0/Z1 导出的原生图**，
               T=8..24 → 10 … 112720 个顶点（GPU 枚举）。首选，没有手挑 Γ。
    对照     = ℤ^D 环面（均匀）：方法自检 —— 应给出 **1 个区域**且 d_s ≈ D。

**等级**：全部【数值证据】。D、T 是输入；Γ 的来源（L2-LATTICE-ORIGIN）仍开放。

用法：/usr/bin/python3 l2_regions.py [T ...]
输出：results/l2_regions.json
"""
from __future__ import annotations

import json
import os
import sys
import time
from collections import defaultdict

import numpy as np
from scipy.sparse import csr_matrix, diags, identity, kron
from scipy.sparse.linalg import eigsh

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

A_WIN, B_WIN = 3.0, 0.3          # d_s 拟合窗口 [A/μ_max, B/μ_min]（标定得来）


# ============================================================ 谱 / 分区
def kmeans(X, k, seed=0, restarts=8, iters=60):
    """朴素 k-means++（不依赖 sklearn）。"""
    X = np.asarray(X, float)
    n = len(X)
    if k >= n:
        return np.arange(n)
    rng = np.random.default_rng(seed)
    best, best_i = None, np.inf
    for _ in range(restarts):
        c = [X[rng.integers(n)]]
        for _ in range(k - 1):
            d2 = np.min(((X[:, None, :] - np.array(c)[None, :, :]) ** 2).sum(-1), axis=1)
            p = d2 / max(d2.sum(), 1e-300)
            c.append(X[rng.choice(n, p=p)])
        C = np.array(c)
        for _ in range(iters):
            d2 = ((X[:, None, :] - C[None, :, :]) ** 2).sum(-1)
            lab = d2.argmin(1)
            nC = np.array([X[lab == j].mean(0) if np.any(lab == j) else C[j] for j in range(k)])
            if np.allclose(nC, C):
                break
            C = nC
        inert = ((X - C[lab]) ** 2).sum()
        if inert < best_i:
            best_i, best = inert, lab.copy()
    return best


def slow_modes(A, m=20):
    """随机走动 P=D^{-1}A 的最大 m 个特征对（对称化 S=D^{-1/2}AD^{-1/2}）。"""
    n = A.shape[0]
    deg = np.asarray(A.sum(1)).ravel()
    deg[deg == 0] = 1.0
    Dm = diags(1.0 / np.sqrt(deg))
    S = (Dm @ A @ Dm).tocsr()
    k = int(min(m, max(2, n - 2)))
    ev, V = eigsh(S, k=k, which="LA")
    idx = np.argsort(-ev)
    ev, V = ev[idx], V[:, idx]
    V = Dm @ V
    V = V / np.maximum(np.linalg.norm(V, axis=0, keepdims=True), 1e-300)
    return ev, V, deg


def choose_k(ev, max_k=8, min_gap_ratio=1.15):
    """区域数由谱隙自己决定（不是手选的）。"""
    info = {"lambda": [round(float(x), 8) for x in ev[:max_k + 1]],
            "gap_definition": "gap_k = (1-λ_k)/(1-λ_{k+1})；取最大者，需 > %.2f" % min_gap_ratio}
    if len(ev) < 4:
        return 1, {**info, "reason": "慢模不足"}
    lam = ev[1:max_k + 1]
    if len(lam) < 2:
        return 1, {**info, "reason": "慢模不足"}
    gaps = (1.0 - lam[:-1]) / np.maximum(1.0 - lam[1:], 1e-14)
    j = int(np.argmax(gaps))
    info["gap_ratios"] = [round(float(g), 4) for g in gaps]
    if float(gaps[j]) < min_gap_ratio:
        return 1, {**info, "reason": f"无明显谱隙（最大间隙比 {gaps[j]:.3f} < {min_gap_ratio}）⇒ 单一区域"}
    return j + 1, {**info, "gap_ratio": round(float(gaps[j]), 4), "k_from_gap": j + 1}


# ============================================================ 维数估计
def lanczos_heat(L, ts, n_vec=8, n_lanc=150, seed=0):
    """
    P(t) = (1/m) tr e^{-tL}，随机向量 Lanczos 求积（只做 matvec）。
    返回 (P, ritz_min, ritz_max)  —— Ritz 值同时给出窗口端点。
    """
    m = L.shape[0]
    rng = np.random.default_rng(seed)
    acc = np.zeros(len(ts))
    all_ritz = []
    for _ in range(n_vec):
        v = rng.standard_normal(m)
        v /= np.linalg.norm(v)
        al, be = [], []
        q, qp, bp = v.copy(), np.zeros(m), 0.0
        nl = int(min(n_lanc, m))
        for _j in range(nl):
            w = L @ q
            a = float(q @ w)
            w = w - a * q - bp * qp
            b = float(np.linalg.norm(w))
            al.append(a)
            if b < 1e-10:
                break
            qp, q, bp = q, w / b, b
            be.append(b)
        bb = be[:len(al) - 1]
        T = np.diag(al) + (np.diag(bb, 1) + np.diag(bb, -1) if len(bb) else 0)
        th, S = np.linalg.eigh(T)
        all_ritz.append(th)
        w0 = S[0, :] ** 2
        acc += np.exp(-np.outer(ts, th)) @ w0
    ritz = np.concatenate(all_ritz)
    return acc / n_vec, float(ritz.min()), float(ritz.max())


def spectral_dim(L, n_t=240, n_vec=8, n_lanc=200, min_decades=1.0, wfrac=0.15):
    """
    内部谱维数：d_s = -2 dlnP/dlnt 的**平台**（P = (1/m)tr e^{-tL}，Lanczos 求积）。
    平台判据：窗口内 P 至少降 min_decades 个数量级，且段内 d_s 的相对起伏最小。
    **不需要 μ_min**（这正是上一版的病根）。
    quality['usable'] = 平台是否够宽 + 起伏是否够小。
    """
    m = L.shape[0]
    if m < 8:
        return None, {"usable": False, "reason": "区域太小"}, {}
    z = max(float(np.asarray(L.diagonal()).mean()), 1.0)   # 平均度 ~ 尺度
    t = np.logspace(np.log10(0.05 / z), np.log10(2000.0 / z), n_t)
    P, rmin, rmax = lanczos_heat(L, t, n_vec=n_vec, n_lanc=min(n_lanc, m))
    P = np.maximum(P, 1e-300)
    lt, lP = np.log(t), np.log(P)
    ds_all = -2 * np.gradient(lP, lt)
    w = max(8, int(wfrac * n_t))
    best = None
    for s0 in range(0, n_t - w):
        dec = float(lP[s0] - lP[s0 + w - 1])
        if dec < min_decades:
            continue
        seg = ds_all[s0:s0 + w]
        cv = float(seg.std() / max(abs(seg.mean()), 1e-9))
        if best is None or cv < best[0]:
            best = (cv, float(np.median(seg)), dec, float(t[s0]), float(t[s0 + w - 1]))
    if best is None:
        return None, {"usable": False, "reason": f"没有跨 {min_decades} 个数量级的窗口",
                      "ritz": [rmin, rmax]}, {}
    cv, d_s, dec, t0, t1 = best
    q = {"usable": bool(cv < 0.05), "plateau_cv": round(cv, 4),
         "P_decades": round(dec, 3), "t_window": [t0, t1],
         "mean_degree": round(z, 3)}
    return d_s, q, {"n_t": n_t, "n_vec": n_vec, "n_lanc": n_lanc}


def growth_dim(A_sub, r_max=24, frac=0.4, vol_frac_cap=0.4):
    """
    生长维数：vol(r) ~ r^{d_g}。**在 vol 超过 40% 节点数之前截断**（否则是有限尺寸饱和，
    这是上一版把 3 维环面测成 1.58 的原因）。只用后半段拟合。
    """
    n = A_sub.shape[0]
    if n < 8:
        return None, None
    deg = np.asarray(A_sub.sum(1)).ravel()
    src = int(np.argmax(deg))
    seen = np.zeros(n, bool)
    seen[src] = True
    frontier = [src]
    vol = []
    cap = vol_frac_cap * n
    for _ in range(r_max):
        nxt = []
        for u in frontier:
            for v in A_sub.indices[A_sub.indptr[u]:A_sub.indptr[u + 1]]:
                if not seen[v]:
                    seen[v] = True
                    nxt.append(int(v))
        frontier = nxt
        if len(seen) > cap:
            break
        vol.append(int(seen.sum()))
        if not frontier:
            break
    if len(vol) < 4:
        return None, None
    r = np.arange(1, len(vol) + 1, dtype=float)
    v = np.array(vol, float)
    s0 = max(1, int(frac * len(r)))
    sl = float(np.polyfit(np.log(r[s0:]), np.log(v[s0:]), 1)[0])
    return sl, {"r_range": [float(r[s0]), float(r[-1])],
                "vol_range": [float(v[s0]), float(v[-1])],
                "vol_frac_at_end": round(float(v[-1] / n), 3)}


def calibrate_ds(tori=((1, 2000), (2, 200), (3, 30), (4, 10), (2, 40), (3, 16))):
    """在已知答案的 ℤ^D 环面上标定 d_s 估计器 —— 误差直接写进结果。"""
    rows = []
    for D, L in tori:
        n = L ** D
        if n > 200_000:
            continue
        A = torus_adj(D, L)
        Lm = diags(np.asarray(A.sum(1)).ravel()) - A
        ds, q, info = spectral_dim(Lm.tocsr())
        rows.append({"D": D, "L": L, "N": int(n), "d_s_measured": None if ds is None else round(ds, 4),
                     "rel_err": None if ds is None else round((ds - D) / D, 4),
                     "quality": q})
    return rows


# ============================================================ 区域统计
def region_report(A, lab, k):
    n = A.shape[0]
    deg = np.asarray(A.sum(1)).ravel()
    rows = []
    for j in range(k):
        idx = np.nonzero(lab == j)[0]
        if len(idx) < 6:
            rows.append({"region": j, "size": int(len(idx)), "note": "太小"})
            continue
        sub = A[idx][:, idx].tocsc()
        d_in = np.asarray(sub.sum(1)).ravel()
        leak = float((deg[idx] - d_in).sum() / max(deg[idx].sum(), 1e-30))
        Lm = (diags(d_in) - sub).tocsr()
        ds, q, _ = spectral_dim(Lm)
        dg, gq = growth_dim(sub)
        rows.append({"region": j, "size": int(len(idx)),
                     "leak_frac": round(leak, 4),
                     "internal_deg_mean": round(float(d_in.mean()), 3),
                     "d_s_internal": None if ds is None else round(ds, 4),
                     "d_growth": None if dg is None else round(dg, 4),
                     "growth_info": gq, "quality": q})
    return rows


def walk_dynamics(A, lab, n_walk=600, steps=3000, seed=1):
    """真实随机轨迹：区域停留时间分布 ＋ 区域间转移矩阵（区域是动力学对象）。"""
    indptr, indices = A.indptr, A.indices
    deg = np.diff(indptr)
    n = A.shape[0]
    rng = np.random.default_rng(seed)
    cur = rng.integers(0, n, size=n_walk)
    lb = lab[cur].astype(int)
    run = np.zeros(n_walk, np.int64)
    dwell = defaultdict(list)
    trans = defaultdict(int)
    urand = rng.random((steps, n_walk))
    for s in range(steps):
        base = indptr[cur]
        pick = base + (urand[s] * deg[cur]).astype(np.int64)
        pick = np.minimum(pick, indptr[cur + 1] - 1)
        cur = indices[pick]
        nl = lab[cur].astype(int)
        run += 1
        ch = np.nonzero(nl != lb)[0]
        for i in ch:
            dwell[int(lb[i])].append(int(run[i]))
            trans[(int(lb[i]), int(nl[i]))] += 1
            run[i] = 0
        if len(ch):
            lb[ch] = nl[ch]
    for i in range(n_walk):
        dwell[int(lb[i])].append(int(run[i]))
    out = []
    for j in sorted(dwell):
        d = np.array(dwell[j], float)
        if len(d) < 10:
            out.append({"region": int(j), "n_dwell": int(len(d))})
            continue
        out.append({"region": int(j), "n_dwell": int(len(d)),
                    "mean_dwell_steps": round(float(d.mean()), 2),
                    "median_dwell_steps": round(float(np.median(d)), 2),
                    "max_dwell_steps": int(d.max()),
                    "std_over_mean": round(float(d.std() / max(d.mean(), 1e-9)), 3),
                    "exit_rate_per_step": round(float(1.0 / d.mean()), 6)})
    K = int(lab.max()) + 1
    M = np.zeros((K, K), int)
    for (u, v), c in trans.items():
        M[u, v] = c
    return out, M.tolist()


# ============================================================ 主流程
def torus_adj(D, L):
    r = np.arange(L)
    A1 = csr_matrix((np.ones(L), (r, (r + 1) % L)), shape=(L, L))
    A1 = A1 + A1.T
    A = A1
    for _ in range(D - 1):
        A = kron(A, identity(L, format="csr")) + kron(identity(A.shape[0], format="csr"), A1)
    return csr_matrix(A)


def analyse(A, tag, meta, m_modes=14, do_dynamics=True, max_k=8):
    t0 = time.time()
    n = A.shape[0]
    ncomp = _n_components(A)
    ev, V, deg = slow_modes(A, m=m_modes)
    k, why = choose_k(ev, max_k=max_k)
    res = {"tag": tag, "N": int(n), "edges": int(A.nnz // 2),
           "mean_degree": round(float(deg.mean()), 4), "components": int(ncomp),
           **meta,
           "slow_eigenvalues": [round(float(x), 8) for x in ev],
           "chosen_k": int(k), "why_k": why}
    if k > 1:
        X = V[:, 1:k]
        X = X / np.maximum(np.linalg.norm(X, axis=0, keepdims=True), 1e-300)
        lab = kmeans(X, k, seed=7)
    else:
        lab = np.zeros(n, int)
    res["region_sizes"] = [int((lab == j).sum()) for j in range(k)]
    res["regions"] = region_report(A, lab, k)
    if do_dynamics and n <= 300_000:
        dw, M = walk_dynamics(A, lab)
        res["dwell"] = dw
        res["transition_counts"] = M
    res["seconds"] = round(time.time() - t0, 2)
    return res, lab


def _n_components(A):
    n = A.shape[0]
    seen = np.zeros(n, bool)
    indptr, indices = A.indptr, A.indices
    c = 0
    for s in range(n):
        if seen[s]:
            continue
        c += 1
        stack = [s]
        seen[s] = True
        while stack:
            u = stack.pop()
            for v in indices[indptr[u]:indptr[u + 1]]:
                if not seen[v]:
                    seen[v] = True
                    stack.append(int(v))
    return c


def main(argv):
    Ts = [int(a) for a in argv if a.isdigit()] or [10, 12, 14, 16, 18, 20, 22, 24]
    out = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"), "runs": [],
           "method": "亚稳态/几乎不变集（P=D^{-1}A 的慢模）+ 谱嵌入 k-means；"
                     "d_s 用内部 Laplacian 热核的对数斜率",
           "window": {"a_over_mu_max": A_WIN, "b_over_mu_min": B_WIN}}

    print("=" * 84)
    print("① 估计器标定：均匀 ℤ^D 环面（已知答案）+ 方法自检（应 1 个区域）")
    out["calibration"] = calibrate_ds()
    for r in out["calibration"]:
        print(f"   D={r['D']} L={r['L']:>4} N={r['N']:>7}: d_s={r['d_s_measured']} "
              f"(相对误差 {r['rel_err']})  可用={r['quality'].get('usable')} "
              f"P跨度={r['quality'].get('P_decades_in_window')}")
    ctrl = []
    for D in (2, 3):
        L = {2: 40, 3: 16}[D]
        r, _ = analyse(torus_adj(D, L), f"torus_D{D}", {"D": D, "L": L},
                       m_modes=8, do_dynamics=False)
        ctrl.append(r)
        print(f"   对照 D={D}: 区域数 k={r['chosen_k']}  d_s={r['regions'][0].get('d_s_internal')}  "
              f"d_g={r['regions'][0].get('d_growth')}")
    out["control_torus"] = ctrl

    print("=" * 84)
    print("② 原生图 Γ = 𝒢_T（闭合类图，Z0/Z1 导出）—— 自发分域？")
    import l0_closure as L0
    from zcl import Engine
    import observable_sweep as OS
    eng = Engine()
    for T in Ts:
        reps, _st = L0.enumerate_necklaces(T, eng, chunk_bits=26, verbose=False)
        A = OS.build_sparse(T, reps, eng)
        r, lab = analyse(A, f"G{T}", {"T": T, "graph": "𝒢_T (native closure-class graph)"})
        out["runs"].append(r)
        ds = [x.get("d_s_internal") for x in r["regions"]]
        dg = [x.get("d_growth") for x in r["regions"]]
        lk = [x.get("leak_frac") for x in r["regions"]]
        print(f"   T={T:>2} N={r['N']:>7,} 分量={r['components']} k={r['chosen_k']} "
              f"({r['why_k'].get('reason', 'gap_ratio=' + str(r['why_k'].get('gap_ratio')))})")
        print(f"        区域大小 {r['region_sizes'][:8]}")
        print(f"        d_s {ds[:8]}")
        print(f"        d_g {dg[:8]}")
        print(f"        泄漏 {lk[:8]}   [{r['seconds']}s]")
        if r.get("dwell"):
            gw = r["dwell"][0]
            print(f"        动力学: 首个区域停留 mean={gw.get('mean_dwell_steps')} "
                  f"std/mean={gw.get('std_over_mean')} (≈1 ⇒ 指数分布)")
        np.save(os.path.join(OUT, f"l2_regions_labels_T{T}.npy"), lab)
        del A
    with open(os.path.join(OUT, "l2_regions.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("=" * 84)
    print("写出 results/l2_regions.json")


if __name__ == "__main__":
    main(sys.argv[1:])
