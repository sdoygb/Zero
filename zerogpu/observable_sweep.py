"""
observable_sweep.py --- 任务 1：中性可观测量普查（observable sweep）

**为什么要这样做**：此前我一直在"挑"可观测量去对 QFT（找极点、找谱密度），
那不是中性做法。本脚本用**机械的生成规则**枚举可观测量，只做**结构分类**，
物理名字留到 `OBSERVABLES.md` 里分开对。

**生成规则（可审计）**
    可观测量 = (对象 O) × (算子/构造 K) × (读出 R)
    O ∈ {Γ 上的行走占据 n_k(v)；𝒢_T 闭链图；词层 {±1}^n（无 Γ 对照）；缺陷算符；记录层账本}
    K ∈ {恒等(计数)、邻接 A、热核 e^{tA}、商/投影、扰动 A+V|0><0|、平移/旋转商}
    R ∈ {点值、谱、分布、熵、参与比、极值、支撑、矩、能隙、关联、标度}

**边界（诚实声明）**
  * 全部结论 = 【数值证据】，不是【导出】。
  * `defect` 扇区里的缺陷位势 V 是**外加输入**（不是 Z0 导出），单独标注。
  * Γ = ℤ^D 环面上的行走用到 D 作为输入（Z0 不给 D，见 G89 命题 1）。
  * 有限环面 L 的选择保证 k < L/2，故环面值与 ℤ^D 精确值一致（k ≤ K 时无绕回）。

**用法**
    /usr/bin/python3 observable_sweep.py walk closure word defect ledger
    /usr/bin/python3 observable_sweep.py all

输出：results/observable_sweep.json  +  OBSERVABLES.md
"""
from __future__ import annotations

import json
import math
import os
import sys
import time

import numpy as np

import structures as ST

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

RECORDS: list[dict] = []


def rec(sector, obj, readout, x, y, meta=None, note=None, maxpts=48, kind="curve"):
    """登记一条可观测量（去掉 NaN/inf，抽稀到 maxpts 点，但保留末点）。"""
    x = np.asarray(x, float).ravel()
    y = np.asarray(y, float).ravel()
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) == 0:
        return
    if len(x) > maxpts:
        idx = np.unique(np.r_[np.linspace(0, len(x) - 1, maxpts - 1).astype(int), len(x) - 1])
        x, y = x[idx], y[idx]
    r = {
        "id": f"{sector}:{obj}:{readout}",
        "sector": sector, "object": obj, "readout": readout, "kind": kind,
        "x": [round(float(v), 6) for v in x],
        "y": [float(f"{v:.10g}") for v in y],
        "meta": meta or {},
    }
    if note:
        r["note"] = note
    RECORDS.append(r)
    return r


def _sub(x, y, maxpts=48):
    x = np.asarray(x, float).ravel()
    y = np.asarray(y, float).ravel()
    if len(x) <= maxpts:
        return x, y
    idx = np.unique(np.r_[np.linspace(0, len(x) - 1, maxpts - 1).astype(int), len(x) - 1])
    return x[idx], y[idx]


# =====================================================================
# 扇区 A：Γ 上的行走占据（"传播"扇区；需要 Γ）
# =====================================================================
KS = {1: 96, 2: 64, 3: 48, 4: 24}


def walk_occupancy(D, K):
    """归一化占据 ñ_k = A^k δ_0 /(2D)^k，环面 L=2K+1（k<=K 时无绕回）。"""
    L = 2 * K + 1
    th = 2 * np.pi * np.fft.fftfreq(L)
    lam = np.zeros((L,) * D)
    for i in range(D):
        shape = [1] * D
        shape[i] = L
        lam = lam + 2 * np.cos(th).reshape(shape)
    lam = (lam / (2 * D)).astype(np.complex128)
    # 到原点的环面距离 r(v) = Σ min(v_i, L-v_i)
    ax = np.arange(L)
    d1 = np.minimum(ax, L - ax).astype(np.float64)
    r2 = np.zeros((L,) * D)
    for i in range(D):
        shape = [1] * D
        shape[i] = L
        r2 = r2 + (d1 ** 2).reshape(shape)
    r = np.sqrt(r2)
    return lam, r, L


def exact_return(D, K):
    """精确整数：n_{2m}(0) = Σ_{a_1+..+a_D=m} (2m)!/∏(a_i!)²。"""
    out = np.zeros(K + 1, dtype=object)
    for m in range(K // 2 + 1):
        tot = 0
        if D == 1:
            tot = math.factorial(2 * m) // (math.factorial(m) ** 2)
        else:
            for a1 in range(m + 1):
                f1 = math.factorial(a1) ** 2
                if D == 2:
                    tot += math.factorial(2 * m) // (f1 * math.factorial(m - a1) ** 2)
                else:
                    for a2 in range(m - a1 + 1):
                        f2 = f1 * math.factorial(a2) ** 2
                        if D == 3:
                            tot += math.factorial(2 * m) // (f2 * math.factorial(m - a1 - a2) ** 2)
                        else:
                            for a3 in range(m - a1 - a2 + 1):
                                a4 = m - a1 - a2 - a3
                                tot += math.factorial(2 * m) // (f2 * math.factorial(a3) ** 2
                                                                 * math.factorial(a4) ** 2)
        out[2 * m] = tot
    return out


def sector_walk(Ds=(1, 2, 3, 4), eps=1e-14):
    rows = []
    for D in Ds:
        K = KS[D]
        t0 = time.time()
        lam, r, L = walk_occupancy(D, K)
        cur = np.ones((L,) * D, dtype=np.complex128)
        R = np.zeros(K + 1)
        IPR = np.zeros(K + 1)
        H = np.zeros(K + 1)
        M2 = np.zeros(K + 1)
        M4 = np.zeros(K + 1)
        SUP = np.zeros(K + 1)
        MX = np.zeros(K + 1)
        FRONT = np.zeros(K + 1)
        for k in range(K + 1):
            n = np.fft.ifftn(cur).real
            np.clip(n, 0.0, None, out=n)
            R[k] = n[(0,) * D]
            IPR[k] = float((n ** 2).sum())
            nz = n[n > 0]
            H[k] = float(-(nz * np.log(nz)).sum())
            M2[k] = float((r ** 2 * n).sum())
            M4[k] = float((r ** 4 * n).sum())
            occ = n > eps
            SUP[k] = float(occ.sum())
            MX[k] = float(n.max())
            FRONT[k] = float(r[occ].max()) if occ.any() else 0.0
            cur = cur * lam
        kk = np.arange(K + 1)
        ks = kk[::2]
        meta = {"D": D, "K": K, "torus_L": L, "note_input": "D 是输入（G89 命题 1）"}
        rec("WALK", f"Z^{D} 行走占据 n_k(v)", "回返 R(k)=n_k(0)", ks, R[::2],
            {**meta, "parity": "odd k 恒为 0（二部图）"})
        rec("WALK", f"Z^{D} 行走占据 n_k(v)", "参与比 IPR=Σn²", ks, IPR[::2], meta)
        rec("WALK", f"Z^{D} 行走占据 n_k(v)", "占据熵 H=-Σp ln p", ks, H[::2], meta)
        rec("WALK", f"Z^{D} 行走占据 n_k(v)", "二阶矩 M2=Σr²p", ks, M2[::2], meta)
        rec("WALK", f"Z^{D} 行走占据 n_k(v)", "峰位概率 max_v p", ks, MX[::2], meta)
        rec("WALK", f"Z^{D} 行走占据 n_k(v)", f"可分辨支撑 #{'{'}p>{eps:g}{'}'}", ks, SUP[::2],
            {**meta, "note": "受浮点分辨限（前沿值 ~(2D)^-k）", "front_radius": FRONT[::2].tolist()})
        rec("WALK", f"Z^{D} 行走占据 n_k(v)", "超额峰度 M4/M2²", ks[2:], (M4 / M2 ** 2)[::2][2:],
            {**meta, "note": "高斯极限下应为常数"})
        # 精确整数核验
        ex = exact_return(D, K)
        exv = np.array([float(ex[i]) for i in ks])
        fftv = R[::2] * float(2 * D) ** ks
        rel = np.abs(fftv - exv) / np.maximum(exv, 1)
        rows.append({"D": D, "K": K,
                     "exact_return_check_max_rel_err": float(rel.max()),
                     "R_2m_exact_last": str(ex[2 * (K // 2)]),
                     "seconds": round(time.time() - t0, 2)})
        print(f"  WALK D={D} K={K}  精确核验 max_rel_err={rel.max():.2e}  "
              f"R({2*(K//2)})={float(ex[2*(K//2)]):.6g}  [{rows[-1]['seconds']}s]", flush=True)
        del lam, cur
    # 邻接算子谱密度（环面上）
    for D in (1, 2, 3, 4):
        L = {1: 1025, 2: 257, 3: 65, 4: 25}[D]
        th = 2 * np.pi * np.arange(L) / L
        lam = np.zeros((L,) * D)
        for i in range(D):
            shape = [1] * D
            shape[i] = L
            lam = lam + 2 * np.cos(th).reshape(shape)
        ev = lam.ravel()
        s = ST.classify_spectrum(ev, nbins=160)
        rec("WALK", f"Z^{D} 邻接算子 A 的谱", "谱密度形状", [0], [0],
            {"D": D, "spectrum_class": s, "L": L, "n_eig": int(ev.size)})
        print(f"  WALK D={D} 谱: {s['label']}  band={s.get('band')}", flush=True)
    return rows


# =====================================================================
# 扇区 B：闭链图 𝒢_T（L0 计数层）
# =====================================================================
def build_sparse(T, reps, eng, chunk=1 << 18):
    from scipy.sparse import csr_matrix
    import l0_closure as L0
    N = len(reps)
    mask = np.uint64((1 << T) - 1)
    kf = eng.kernel("l0_fused", L0.FUSED, "neighbors_canon")
    rs, cs = [], []
    for s in range(0, N, chunk):
        e = min(s + chunk, N)
        sub = np.ascontiguousarray(reps[s:e])
        m = e - s
        rb = eng.to_device(sub)
        cb = eng.empty(m * T, np.uint64)
        eng.run(kf, (m,), rb, cb, np.int32(T), mask)
        can = eng.from_device(cb, m * T, np.uint64).reshape(m, T)
        idx = np.searchsorted(reps, can)
        np.clip(idx, 0, N - 1, out=idx)
        ok = reps[idx] == can
        rows = np.arange(s, e, dtype=np.int64)[:, None]
        ok &= idx != rows
        ii = np.broadcast_to(rows, (m, T))[ok]
        jj = idx[ok]
        rs.append(ii)
        cs.append(jj)
    r = np.concatenate(rs)
    c = np.concatenate(cs)
    A = csr_matrix((np.ones(len(r)), (r, c)), shape=(N, N))
    A.sum_duplicates()
    A.data[:] = 1.0
    return A


def kpm_dos(A, n_mom=800, n_vec=8, seed=0, lam_lo=None, lam_hi=None):
    """KPM + Jackson 核的谱密度。返回 (网格, 密度)。"""
    from scipy.sparse.linalg import eigsh
    N = A.shape[0]
    if lam_hi is None:
        lam_hi = float(eigsh(A, k=1, which="LA", return_eigenvectors=False)[0])
    if lam_lo is None:
        lam_lo = float(eigsh(A, k=1, which="SA", return_eigenvectors=False)[0])
    # 6% 余量 + 内缩网格：避开 KPM 在 |x|=1 处 1/sqrt(1-x^2) 的假发散
    a = (lam_hi - lam_lo) / 2 * 1.06
    b = (lam_hi + lam_lo) / 2
    H = (A - b * __import__("scipy.sparse", fromlist=["eye"]).eye(N, format="csr")) / a
    rng = np.random.default_rng(seed)
    mu = np.zeros(n_mom + 1)
    for _ in range(n_vec):
        v = rng.standard_normal(N)
        v /= np.linalg.norm(v)
        t0 = v.copy()
        t1 = H @ v
        mu[0] += float(v @ t0)
        mu[1] += float(v @ t1)
        for n in range(2, n_mom + 1):
            t2 = 2 * (H @ t1) - t0
            mu[n] += float(v @ t2)
            t0, t1 = t1, t2
    mu /= n_vec
    # Jackson 核
    nm = n_mom
    g = ((nm - np.arange(nm + 1) + 1) * np.cos(np.pi * np.arange(nm + 1) / (nm + 1))
         + np.sin(np.pi * np.arange(nm + 1) / (nm + 1)) / np.tan(np.pi / (nm + 1)))
    g /= (nm + 1)
    x = np.linspace(-0.985, 0.985, 2001)
    Tn = np.cos(np.arange(nm + 1)[:, None] * np.arccos(np.clip(x, -1, 1))[None, :])
    dens = (g * mu)[:, None] * Tn
    dens = (dens[0] + 2 * dens[1:].sum(axis=0)) / (np.pi * np.sqrt(np.clip(1 - x ** 2, 1e-12, None)))
    ev = x * a + b
    return ev, dens / a


def spacing_ratio(ev, collapse: bool = False, ndigits: int = 8):
    """
    能级间距比 <r>（**不用展开**，这是 r 统计的标准定义）。
    Poisson 0.386 / GOE 0.531 / 刚性（等间距）→1。
    collapse=True 时先把简并本征值合并（只留不同值）。
    """
    ev = np.sort(np.asarray(ev, float))
    ev = ev[np.isfinite(ev)]
    if collapse:
        ev = np.unique(np.round(ev, ndigits))
    if len(ev) < 12:
        return None
    s = np.diff(ev)
    s = s[s > 0]
    if len(s) < 10:
        return None
    r = np.minimum(s[1:], s[:-1]) / np.maximum(s[1:], s[:-1])
    return float(r.mean())


def spectrum_shape_stats(ev, ndigits: int = 8):
    """简并/团簇/空隙的结构统计（比 gapped/smooth 更细一层）。"""
    ev = np.sort(np.asarray(ev, float))
    ev = ev[np.isfinite(ev)]
    N = len(ev)
    vals, cnt = np.unique(np.round(ev, ndigits), return_counts=True)
    w = float(ev.max() - ev.min())
    d = np.diff(vals)
    out = {"N": int(N), "distinct": int(len(vals)),
           "distinct_frac": round(float(len(vals) / max(1, N)), 4),
           "max_multiplicity": int(cnt.max()) if len(cnt) else 0,
           "mean_multiplicity": round(float(cnt.mean()), 4) if len(cnt) else 0,
           "max_gap_rel": round(float(d.max() / w), 4) if len(d) and w > 0 else 0.0,
           "bandwidth": round(w, 6)}
    if w > 0 and len(d):
        cut = 0.01 * w
        cl = 1 + int((d > cut).sum())
        out["clusters_at_1pct_bandwidth"] = cl
        out["cluster_frac"] = round(float(cl / max(1, len(vals))), 4)
    return out


def is_bipartite(A):
    """BFS 二着色。二部图 => 谱关于 0 对称（λmin = -λmax），可省一次 SA Lanczos。"""
    from collections import deque
    n = A.shape[0]
    indptr, indices = A.indptr, A.indices
    col = np.full(n, -1, np.int8)
    for s in range(n):
        if col[s] != -1:
            continue
        col[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for v in indices[indptr[u]:indptr[u + 1]]:
                if col[v] == -1:
                    col[v] = 1 - col[u]
                    q.append(int(v))
                elif col[v] == col[u]:
                    return False
    return True


def sector_closure(Ts=(8, 10, 12, 14, 16, 18, 20, 22, 24), dense_max=4200):
    from zcl import Engine
    import l0_closure as L0
    from scipy.sparse.linalg import eigsh
    eng = Engine()
    print("  device:", eng.info()["name"], flush=True)
    rows = []
    for T in Ts:
        t0 = time.time()
        reps, st = L0.enumerate_necklaces(T, eng, chunk_bits=26, verbose=False)
        N = len(reps)
        A = build_sparse(T, reps, eng)
        deg = np.asarray(A.sum(axis=1)).ravel()
        meta = {"T": T, "nodes": N, "edges": int(A.nnz // 2),
                "avg_degree": round(float(deg.mean()), 6)}
        # --- 度分布 ---
        hist = np.bincount(deg.astype(int))
        rec("CLOSURE", f"𝒢_{T}", "度分布 P(deg)", np.arange(len(hist)), hist / N, meta, kind="density")
        # --- 谱 ---
        bip = is_bipartite(A)
        if N <= dense_max:
            ev = np.linalg.eigvalsh(A.toarray())
            evs = ev
            lmax = float(evs.max())
            lmin = float(evs.min())
            shape = spectrum_shape_stats(evs)
            r_raw = spacing_ratio(evs, collapse=False)
            r_dis = spacing_ratio(evs, collapse=True)
            ipr_mean = None
            if N <= 1500:
                w, V = np.linalg.eigh(A.toarray())
                ipr_mean = float(((V ** 2) ** 2).sum(axis=0).mean())
        else:
            lmax = float(eigsh(A, k=1, which="LA", return_eigenvectors=False)[0])
            lmin = -lmax if bip else float(eigsh(A, k=1, which="SA", return_eigenvectors=False)[0])
            shape, r_raw, r_dis, ipr_mean = None, None, None, None
            evs = np.array([lmin, lmax])
        spec = ST.classify_spectrum(evs, nbins=200) if N <= dense_max else None
        # --- 谱密度 ---
        if N > dense_max:
            evk, dk = kpm_dos(A, n_mom=400, n_vec=6, lam_lo=lmin, lam_hi=lmax)
            spec = ST.classify_density(evk, dk, n_eff=N)
            rec("CLOSURE", f"𝒢_{T}", "KPM 谱密度 ρ(λ)", evk, dk, {**meta, "method": "KPM-Jackson"},
                kind="density")
        else:
            h, e = np.histogram(evs, bins=min(80, max(8, N // 4)))
            rec("CLOSURE", f"𝒢_{T}", "谱密度 ρ(λ)", 0.5 * (e[1:] + e[:-1]), h / max(1, N),
                {**meta, "method": "dense-hist"}, kind="density")
        rec("CLOSURE", f"𝒢_{T}", "谱密度形状", [0], [0],
            {**meta, "spectrum_class": spec, "shape_stats": shape,
             "spacing_ratio_r_raw": r_raw, "spacing_ratio_r_distinct": r_dis,
             "lambda_max": lmax, "lambda_min": lmin, "bipartite": bool(bip),
             "mean_IPR": ipr_mean, "dense": bool(N <= dense_max)}, kind="spectrum")
        # --- 热核 Z(t)=⟨e^{tλ}⟩ 的增长率（t→∞ 时应 →λmax）---
        ts = np.linspace(0, 6.0, 61)
        if N <= dense_max:
            Z = np.array([np.exp(np.clip(t * evs, -700, 700)).mean() for t in ts])
        else:
            Z = np.array([float(np.trapz(dk * np.exp(np.clip(t * evk, -700, 700)), evk))
                          for t in ts])
        rec("CLOSURE", f"𝒢_{T}", "热核 Z(t)=⟨e^{tλ}⟩", ts, Z, {**meta, "note": "增长率→λmax"})
        rows.append({**meta, "spectrum": spec["label"], "r_distinct": r_dis, "bipartite": bool(bip),
                     "lambda_max": float(lmax), "lambda_min": float(lmin),
                     "seconds": round(time.time() - t0, 2)})
        print(f"  CLOSURE T={T:>2} N={N:>7,} E={A.nnz//2:>8,}  {spec['label']:<26} "
              f"r_dis={('%.4f' % r_dis) if r_dis is not None else '  -  '}  "
              f"lmax={lmax:.3f}  bip={bip}  [{rows[-1]['seconds']}s]", flush=True)
        del A
    # --- 系综标度：本次 T=8..24 ＋ 复用的 l0_scan.json（T 到 32）---
    ens = {}
    for key, val in (("nodes", lambda r: r["nodes"]), ("edges", lambda r: r["edges"]),
                     ("avg_degree", lambda r: r["avg_degree"])):
        for r in rows:
            ens.setdefault(key, {})[r["T"]] = val(r)
    p = os.path.join(OUT, "l0_scan.json")
    if os.path.exists(p):
        for r in json.load(open(p))["rows"]:
            for key in ("nodes", "edges", "avg_degree"):
                ens.setdefault(key, {}).setdefault(r["T"], r[key])
    T = np.array(sorted(ens["nodes"]), float)
    ad = np.array([ens["avg_degree"][t] for t in T], float)
    nd = np.array([ens["nodes"][t] for t in T], float)
    rec("CLOSURE", "𝒢_T 系综", "平均度 vs T", T, ad, {"T_range": [float(T[0]), float(T[-1])],
        "source": "本次 T<=24 + l0_scan.json"})
    rec("CLOSURE", "𝒢_T 系综", "平均度 - T/2 vs T", T, ad - T / 2, {"source": "同上"})
    rec("CLOSURE", "𝒢_T 系综", "节点数 vs T", T, nd, {"source": "同上"})
    return rows


# =====================================================================
# 扇区 C：词层（无 Γ 对照）
# =====================================================================
def sector_word():
    """无 Γ 的最小对象：±1 词本身。看没有 Γ 时剩下什么结构。"""
    import itertools
    ns = list(range(2, 19, 2))
    # 1) 两点关联（对所有平衡词穷举）
    for n in (8, 10, 12, 14):
        words = np.array([w for w in itertools.product((1, -1), repeat=n) if sum(w) == 0], float)
        C = np.array([float((words[:, 0] * words[:, i]).mean()) for i in range(n)])
        rec("WORD", f"平衡词 {{±1}}^{n}", "两点关联 C(r)=⟨w_0 w_r⟩", np.arange(n), C,
            {"n": n, "n_words": int(len(words))})
        if n == 14:
            S = np.abs(np.fft.rfft(words, axis=1)) ** 2
            S = S.mean(axis=0)
            rec("WORD", f"平衡词 {{±1}}^{n}", "功率谱 S(q)", np.arange(len(S)), S, {"n": n})
    # 2) 平衡词计数占比（无 Γ 时唯一"增长"结构）
    nn = np.arange(2, 42, 2)
    tot = 2.0 ** nn
    bal = np.array([float(math.comb(int(k), int(k) // 2)) for k in nn])
    rec("WORD", "{±1}^n 全空间", "平衡词占比 C(n,n/2)/2^n", nn, bal / tot, {"note": "无 Γ 的对照"})
    rec("WORD", "{±1}^n 全空间", "平衡词计数增长", nn, bal, {"note": "无 Γ 的对照"})
    # 3) 一维零和随机游走回返（Γ 退化为线 = 无 Γ 时的隐式假设）
    n = np.arange(1, 61)
    R = np.array([float(math.comb(2 * int(m), int(m))) / 4.0 ** int(m) for m in n])
    rec("WORD", "1D 零和游走", "回返 R(2m)=C(2m,m)/4^m", 2 * n, R,
        {"note": "这就是 WALK D=1 —— '无 Γ' 的隐含假设其实是一条线"})
    print("  WORD 对照扇区完成", flush=True)


# =====================================================================
# 扇区 D：缺陷算符（带电/局域态的机器）—— V 是外加输入
# =====================================================================
def sector_defect():
    from scipy.sparse import csr_matrix, identity, kron
    from scipy.sparse.linalg import eigsh

    def torus_adj(D, L):
        # 1D 环的邻接
        r = np.arange(L)
        A1 = csr_matrix((np.ones(L), (r, (r + 1) % L)), shape=(L, L))
        A1 = A1 + A1.T
        A = A1
        for _ in range(D - 1):
            A = kron(A, identity(L, format="csr")) + kron(identity(A.shape[0], format="csr"), A1)
        return csr_matrix(A)

    LD = {1: (201,), 2: (101,), 3: (31,), 4: (13,)}
    Vs_pos = [0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0]
    Vs_neg = [-v for v in Vs_pos]
    rows = []
    for D in (1, 2, 3, 4):
        L = LD[D][0]
        N = L ** D
        A = torus_adj(D, L)
        idx0 = 0
        band = 2.0 * D
        meta = {"D": D, "L": L, "sites": int(N), "band_edge": band,
                "note_input": "缺陷位势 V 是外加输入（不是 Z0 导出）"}
        Epos, Eneg = [], []
        psi_store = {}
        for V in Vs_pos + Vs_neg:
            H = A + V * csr_matrix(([1.0], ([idx0], [idx0])), shape=(N, N))
            try:
                if V > 0:
                    w, v = eigsh(H, k=1, which="LA")
                else:
                    w, v = eigsh(H, k=1, which="SA")
            except Exception:
                continue
            E = float(w[0])
            (Epos if V > 0 else Eneg).append((V, E))
            if abs(V) in (1.0, 4.0):
                psi_store[V] = (E, np.asarray(v[:, 0]).ravel())
        # 束缚能 vs 耦合
        for tag, arr in (("V>0（上带边 +2D）", Epos), ("V<0（下带边 -2D）", Eneg)):
            if len(arr) < 3:
                continue
            arr = sorted(arr, key=lambda t: abs(t[0]))
            V = np.array([a[0] for a in arr], float)
            E = np.array([a[1] for a in arr], float)
            dE = E - np.sign(V) * band
            rec("DEFECT", f"Z^{D} + V|0><0|", f"束缚态能移 ΔE {tag}", V, dE, meta,
                note="ΔE<0 且 |ΔE| 单调增 = 束缚态；标度见结构分类")
        # 局域长度：沿径向的指数衰减
        for V, (E, psi) in psi_store.items():
            sh = psi.reshape((L,) * D)
            site0 = [0] * D
            vals = {}
            for off in range(1, L // 2):
                # 沿第一根轴取点
                ii = list(site0)
                ii[0] = off
                vals[off] = abs(float(sh[tuple(ii)]))
            off = np.array(sorted(vals), float)
            amp = np.array([vals[int(o)] for o in off], float)
            ok = amp > 0
            if ok.sum() >= 4:
                sl = np.polyfit(off[ok], np.log(amp[ok]), 1)
                xi = -1.0 / sl[0] if sl[0] < 0 else np.inf
                rec("DEFECT", f"Z^{D} + V|0><0|", f"束缚态径向衰减 |ψ(r)| (V={V})", off[ok], amp[ok],
                    {**meta, "E": E, "xi_fit": round(float(xi), 4)},
                    note="直线段 → 指数衰减 → 局域长度 ξ")
        rows.append({"D": D, "L": L, "band_edge": band,
                     "Epos": Epos, "Eneg": Eneg})
        print(f"  DEFECT D={D} L={L} sites={N}: "
              f"E(V=8)={[e for v,e in Epos if v==8][:1]}  "
              f"E(V=-8)={[e for v,e in Eneg if v==-8][:1]}", flush=True)
        del A
    # ---- 带边 DOS 指数：λ=2Σcos k_i 的分块随机采样，δ=2D-λ 对数分箱 ----
    rng = np.random.default_rng(0)
    expo = []
    for D in (1, 2, 3, 4):
        n_tot = 20_000_000 if D <= 3 else 12_000_000
        lo_b, hi_b = (1e-5, 0.3) if D <= 2 else (1e-3, 0.3)
        bins = np.logspace(np.log10(lo_b), np.log10(hi_b), 26)
        h = np.zeros(len(bins) - 1)
        edge = 2.0 * D
        done = 0
        while done < n_tot:
            m = min(2_000_000, n_tot - done)
            th = rng.uniform(-np.pi, np.pi, size=(m, D))
            lam = 2 * np.cos(th).sum(axis=1)
            d = edge - lam
            sel = (d > bins[0]) & (d < bins[-1])
            h += np.histogram(d[sel], bins=bins)[0]
            done += m
        w = np.diff(bins)
        c = 0.5 * (bins[1:] + bins[:-1])
        dens = h / (n_tot * w)
        ok = dens > 0
        e_exp = None
        if ok.sum() >= 6:
            sl = np.polyfit(np.log(c[ok]), np.log(dens[ok]), 1)
            e_exp = round(float(sl[0]), 4)
        expo.append({"D": D, "exponent": e_exp, "n_samples": n_tot,
                     "fit_points": int(ok.sum()), "theory_D_over_2_minus_1": D / 2 - 1})
        rec("DEFECT", f"Z^{D} 带边 DOS", "ρ(2D-δ) vs δ", c[ok], dens[ok],
            {"D": D, "edge": edge, "theory_exponent": D / 2 - 1,
             "note": "指数 = D/2-1：<0 发散(任意弱耦合成键)；=0 对数发散；>0 有限(有阈值耦合)"},
            kind="density")
        print(f"  DEFECT DOS 带边指数 D={D}: {e_exp} "
              f"(理论 D/2-1 = {D/2-1})", flush=True)
    rec("DEFECT", "带边 DOS", "带边指数汇总", [0], [0], {"edge_exponents": expo}, kind="spectrum")

    # ---- 临界耦合 Vc(N)：极值态的参与比 IPR(V) 的局域化转折 ----
    vc_rows = []
    for D in (1, 2, 3, 4):
        Ls = {1: [101, 201, 401, 801], 2: [21, 41, 81, 161],
              3: [11, 17, 25, 37], 4: [7, 9, 11, 13]}[D]
        Vgrid = np.logspace(-3, 1.3, 26)
        tab = []
        for L in Ls:
            N = L ** D
            A = torus_adj(D, L)
            ipr = np.full(len(Vgrid), np.nan)
            for j, V in enumerate(Vgrid):
                H = A + V * csr_matrix(([1.0], ([0], [0])), shape=(N, N))
                try:
                    w, v = eigsh(H, k=1, which="LA")
                except Exception:
                    continue
                p = np.asarray(v[:, 0]).ravel() ** 2
                ipr[j] = float((p ** 2).sum())
            thr = 20.0 / N
            vc = None
            idx = np.nonzero(np.isfinite(ipr) & (ipr > thr))[0]
            if len(idx):
                vc = float(Vgrid[idx[0]])
            rec("DEFECT", f"Z^{D} L={L}", "极值态 IPR vs V", Vgrid, ipr,
                {"D": D, "L": L, "sites": int(N), "IPR_extended=1/N": 1.0 / N,
                 "localized_threshold": thr,
                 "note": "IPR 从 1/N 跳到 O(1) = 从扩展态变成束缚态"})
            tab.append({"L": L, "sites": int(N), "Vc": vc,
                        "Vc_times_logN": None if vc is None else round(vc * math.log(N), 4),
                        "Vc_times_N": None if vc is None else round(vc * N, 4)})
            print(f"  DEFECT D={D} L={L:>4} N={N:>7}: Vc={vc}", flush=True)
        vc_rows.append({"D": D, "rows": tab})
        good = [t for t in tab if t["Vc"]]
        if len(good) >= 2:
            NN = np.array([t["sites"] for t in good], float)
            VV = np.array([t["Vc"] for t in good], float)
            rec("DEFECT", f"Z^{D} 缺陷", "临界耦合 Vc vs 格点数 N", NN, VV, {"D": D},
                note="Vc→0 = 任意弱耦合成键；Vc→常数 = 有阈值")
    rec("DEFECT", "缺陷算符", "临界耦合汇总", [0], [0], {"Vc_table": vc_rows})
    return rows


# =====================================================================
# 扇区 E：记录层账本（L1′）
# =====================================================================
def sector_ledger():
    p = os.path.join(OUT, "z0_record.json")
    if not os.path.exists(p):
        return
    d = json.load(open(p))
    rows = [r for r in d["rows"] if r.get("events", 0) > 0]
    n = np.array([r["n"] for r in rows], float)
    cov = np.array([r["coverage"] for r in rows], float)
    rec("LEDGER", "记录层 L1′", "闭合类覆盖率 vs n", n, cov,
        {"source": "z0_record.json", "note": "0.7657 = 文档值"})
    rec("LEDGER", "记录层 L1′", "对数缺口 ln(1-coverage) vs n", n, np.log1p(-cov),
        {"source": "z0_record.json"})
    rec("LEDGER", "记录层 L1′", "独立闭合类数 vs K_n", n,
        np.array([r["distinct_classes"] for r in rows], float),
        {"K_n": [r["K_n"] for r in rows], "source": "z0_record.json"})
    qp = os.path.join(OUT, "q_L.json")
    if os.path.exists(qp):
        q = json.load(open(qp))
        rr = q["rows"]
        L = np.array([r["L"] for r in rr], float)
        rec("LEDGER", "类权重 q_L", "q_L vs L", L, np.array([r["q_L"] for r in rr], float),
            {"source": "q_L.json"})
    print("  LEDGER 扇区完成", flush=True)


# =====================================================================
def main(argv):
    which = argv or ["all"]
    if "all" in which:
        which = ["walk", "closure", "word", "defect", "ledger"]
    info = {}
    t0 = time.time()
    if "walk" in which:
        print("[A] WALK 扇区", flush=True)
        info["walk"] = sector_walk()
    if "closure" in which:
        print("[B] CLOSURE 扇区", flush=True)
        info["closure"] = sector_closure()
    if "word" in which:
        print("[C] WORD 扇区（无 Γ 对照）", flush=True)
        sector_word()
    if "defect" in which:
        print("[D] DEFECT 扇区", flush=True)
        info["defect"] = sector_defect()
    if "ledger" in which:
        print("[E] LEDGER 扇区", flush=True)
        sector_ledger()

    # --- 分类 ---
    for r in RECORDS:
        k = r.get("kind", "curve")
        if k == "spectrum":
            sc = r["meta"].get("spectrum_class") or {}
            r["structure"] = {"label": sc.get("label", "declared"),
                              "params": sc.get("params", {}),
                              "detail": {kk: vv for kk, vv in r["meta"].items()
                                         if kk != "spectrum_class"}}
        elif k == "density":
            st = ST.classify_density(r["x"], r["y"])
            st["params"] = {**(st.get("params") or {}),
                            **{kk: vv for kk, vv in (r.get("meta") or {}).items()
                               if kk in ("method", "theory_exponent", "IPR_extended=1/N")}}
            r["structure"] = st
        elif len(r["x"]) >= 5 and not (len(set(r["y"])) == 1 and r["y"][0] == 0):
            r["structure"] = ST.classify_curve(r["x"], r["y"])
        else:
            r["structure"] = {"label": "declared", "params": {}}
    clusters = ST.summarize(RECORDS)
    out = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
           "seconds": round(time.time() - t0, 1),
           "classifier": "structures.classify_curve (对数空间 R²，复杂度和同等则取简)",
           "sectors": info, "clusters": clusters, "records": RECORDS}
    with open(os.path.join(OUT, "observable_sweep.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(f"\n登记 {len(RECORDS)} 条可观测量，用时 {out['seconds']}s")
    print("结构聚类：")
    for k, v in clusters.items():
        print(f"  {k:<28} {len(v):>3} 条")
    print("\n按结构分类的目录见 OBSERVABLES.md（用 make_catalogue.py 生成）")


if __name__ == "__main__":
    main(sys.argv[1:])
