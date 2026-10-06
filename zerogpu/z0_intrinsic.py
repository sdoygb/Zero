"""
z0_intrinsic.py --- 把 §12 的类标签**内禀化**：Weisfeiler–Lehman 着色是自同构不变的

判据：WL 颜色在 Aut(Γ) 下不变 ⟹ 同一 WL 色内的顶点**本应无法区分**。
    · 类数在**同色内仍有差异** ⟹ 那部分是**编号假象**（§12 的警告成立）
    · 类数**只随 WL 色变**       ⟹ 那部分是**内禀区分度**（D222 保留机制才有作用对象）

用法：/usr/bin/python3 z0_intrinsic.py → results/z0_intrinsic.json
"""
from __future__ import annotations
import json, os, time
import numpy as np
from scipy.sparse import csr_matrix, identity, kron

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results"); os.makedirs(OUT, exist_ok=True)


def wl_colors(A, iters=8):
    """Weisfeiler–Lehman 精化：颜色在 Aut(Γ) 下不变。"""
    V = A.shape[0]
    indptr, indices = A.indptr, A.indices
    col = np.asarray(A.sum(1)).ravel().astype(np.int64)
    for _ in range(iters):
        sigs = []
        for v in range(V):
            nb = np.sort(col[indices[indptr[v]:indptr[v + 1]]])
            sigs.append((int(col[v]), tuple(int(x) for x in nb)))
        uniq = {s: i for i, s in enumerate(sorted(set(sigs)))}
        new = np.array([uniq[s] for s in sigs], np.int64)
        if new.tolist() == col.tolist():
            break
        col = new
    return col


def torus_adj(D, L):
    r = np.arange(L)
    A1 = csr_matrix((np.ones(L), (r, (r + 1) % L)), shape=(L, L)); A1 = A1 + A1.T
    A = A1
    for _ in range(D - 1):
        A = kron(A, identity(L, format="csr")) + kron(identity(A.shape[0], format="csr"), A1)
    return csr_matrix(A)


def add_hubs(A, frac=0.25, extra=2, seed=0):
    rng = np.random.default_rng(seed); V = A.shape[0]; A = A.tolil()
    for u in rng.choice(V, size=max(1, int(frac * V)), replace=False):
        for _ in range(extra):
            v = int(rng.integers(0, V))
            if v != u:
                A[u, v] = 1.0; A[v, u] = 1.0
    return csr_matrix(A)


def min_rotation(d, k, base):
    best = cur = d; sh = base ** (k - 1)
    for _ in range(k - 1):
        cur = (cur % sh) * base + (cur // sh)
        if cur < best:
            best = cur
    return best


def class_counts(A, kmax=8, cap=6_000_000):
    V = A.shape[0]; indptr, indices = A.indptr, A.indices
    deg = np.diff(indptr); dmax = int(deg.max()); base = dmax
    slot = np.full((V, dmax), -1, np.int64)
    for v in range(V):
        slot[v, :deg[v]] = indices[indptr[v]:indptr[v + 1]]
    st0 = np.repeat(np.arange(V), dmax); sl = np.tile(np.arange(dmax), V)
    ok = slot.ravel() >= 0
    st0, sl = st0[ok], sl[ok]; cur = slot.ravel()[ok]; seq = sl.copy()
    sets = [set() for _ in range(V)]
    for k in range(1, kmax + 1):
        cl = cur == st0
        for s, q in zip(st0[cl], seq[cl]):
            sets[int(s)].add(min_rotation(int(q), k, base))
        if len(cur) > cap:
            break
        nxt = slot[cur]; m = nxt.shape[0]; valid = nxt >= 0
        dig = np.broadcast_to(np.arange(dmax)[None, :], (m, dmax))
        cur, seq, st0 = (nxt[valid], (np.repeat(seq[:, None], dmax, 1) * base + dig)[valid],
                         np.repeat(st0[:, None], dmax, 1)[valid])
    return np.array([len(s) for s in sets], float)


def main():
    t0 = time.time(); out = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"), "runs": []}
    L = 4; base = torus_adj(2, L)
    graphs = [("torus4x4_transitive", base), ("torus4x4_hubs", add_hubs(base, .25, 2, 1))]
    try:
        import l0_closure as L0, observable_sweep as OS
        from zcl import Engine
        eng = Engine(); reps, _ = L0.enumerate_necklaces(10, eng, chunk_bits=26, verbose=False)
        graphs.append(("G_10", OS.build_sparse(10, reps, eng)))
    except Exception as e:
        print("  [warn]", e)
    for name, A in graphs:
        V = A.shape[0]
        col = wl_colors(A); nc = class_counts(A)
        rows = []
        for c in np.unique(col):
            m = col == c
            rows.append({"wl_color": int(c), "n_vertices": int(m.sum()),
                         "class_count_min": int(nc[m].min()), "max": int(nc[m].max()),
                         "std": round(float(nc[m].std()), 2)})
        within = float(np.mean([r["std"] for r in rows]))
        between = float(np.std([nc[col == c].mean() for c in np.unique(col)]))
        print("=" * 92)
        print(f"### Γ={name}  V={V}   WL 色数 = {len(np.unique(col))}（自同构不变）")
        print(f"    {'WL色':>5} {'顶点数':>6} {'类数min':>9} {'类数max':>9} {'同色内std':>9}")
        for r in rows[:10]:
            print(f"    {r['wl_color']:>5} {r['n_vertices']:>6} {r['class_count_min']:>9} "
                  f"{r['class_count_max'] if False else r['max']:>9} {r['std']:>9}")
        print(f"    ⟹ 同色内平均 std = {within:.1f}（**编号假象**） vs 色间 std = {between:.1f}（**内禀区分度**）")
        out["runs"].append({"graph": name, "V": V, "n_wl_colors": int(len(np.unique(col))),
                            "per_color": rows, "within_color_std": within,
                            "between_color_std": between,
                            "ratio_intrinsic_over_artifact": round(between / max(within, 1e-9), 2)})
    with open(os.path.join(OUT, "z0_intrinsic.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("=" * 92)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_intrinsic.json")
    print("读法：色间 std ≫ 同色内 std ⟹ 区分度是内禀的（D222 保留机制可用）；")
    print("      色间 std ≲ 同色内 std ⟹ 之前测到的区分度基本是编号假象（硬 no-go）。")


if __name__ == "__main__":
    main()
