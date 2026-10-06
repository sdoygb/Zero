"""
z0_retain.py --- 用**内禀闭合类**做 D222 的保留竞争：存活集会碎裂吗？

规则（照 D222："只保留最高两层"、"被删去的低层不能再影响下一代"）：
    类人口 pop[c] = #{v : 闭合类 c 出现在 v 的类集合里}
    **保留类** = pop 最高的前 K 个
    **存活顶点** = {v : C_v ∩ 保留类 ≠ ∅}      （其余熄灭，不再被重播种）
测：存活集的连通块数（= 畴数）与块大小分布；与"长度标签"版（§11 恒 1 块）对照。

用法：/usr/bin/python3 z0_retain.py → results/z0_retain.json
"""
from __future__ import annotations
import json, os, time
from collections import deque
import numpy as np
from scipy.sparse import csr_matrix, identity, kron

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)


def torus_adj(D, L):
    r = np.arange(L); A1 = csr_matrix((np.ones(L), (r, (r+1) % L)), shape=(L, L)); A1 = A1 + A1.T
    A = A1
    for _ in range(D-1):
        A = kron(A, identity(L, format="csr")) + kron(identity(A.shape[0], format="csr"), A1)
    return csr_matrix(A)


def add_hubs(A, frac=.25, extra=2, seed=0):
    rng = np.random.default_rng(seed); V = A.shape[0]; A = A.tolil()
    for u in rng.choice(V, size=max(1, int(frac*V)), replace=False):
        for _ in range(extra):
            v = int(rng.integers(0, V))
            if v != u: A[u, v] = 1.0; A[v, u] = 1.0
    return csr_matrix(A)


def min_rot(d, k, base):
    best = cur = d; sh = base ** (k-1)
    for _ in range(k-1):
        cur = (cur % sh) * base + (cur // sh)
        if cur < best: best = cur
    return best


def class_sets(A, kmax=8):
    V = A.shape[0]; ip, ix = A.indptr, A.indices
    deg = np.diff(ip); dmax = int(deg.max()); base = dmax
    slot = np.full((V, dmax), -1, np.int64)
    for v in range(V): slot[v, :deg[v]] = ix[ip[v]:ip[v+1]]
    st0 = np.repeat(np.arange(V), dmax); sl = np.tile(np.arange(dmax), V)
    ok = slot.ravel() >= 0
    st0, cur, seq = st0[ok], slot.ravel()[ok], sl[ok]
    sets = [set() for _ in range(V)]
    for k in range(1, kmax+1):
        cl = cur == st0
        for s, q in zip(st0[cl], seq[cl]): sets[int(s)].add(min_rot(int(q), k, base))
        if len(cur) > 6_000_000: break
        nxt = slot[cur]; m = nxt.shape[0]; valid = nxt >= 0
        dig = np.broadcast_to(np.arange(dmax)[None, :], (m, dmax))
        cur, seq, st0 = (nxt[valid], (np.repeat(seq[:, None], dmax, 1)*base + dig)[valid],
                         np.repeat(st0[:, None], dmax, 1)[valid])
    return sets


def comps(mask, A):
    V = A.shape[0]; ip, ix = A.indptr, A.indices
    lab = np.full(V, -1, np.int64); c = 0
    for s in np.nonzero(mask)[0]:
        if lab[s] >= 0: continue
        lab[s] = c; q = deque([int(s)])
        while q:
            u = q.popleft()
            for v in ix[ip[u]:ip[u+1]]:
                if mask[v] and lab[v] < 0: lab[v] = c; q.append(int(v))
        c += 1
    return np.bincount(lab[lab >= 0]).tolist() if c else []


def main():
    t0 = time.time(); out = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"), "runs": []}
    L = 4; base = torus_adj(2, L)
    graphs = [("torus4x4_transitive", base), ("torus4x4_hubs", add_hubs(base, .25, 2, 1))]
    try:
        import l0_closure as L0, observable_sweep as OS
        from zcl import Engine
        eng = Engine(); reps, _ = L0.enumerate_necklaces(10, eng, chunk_bits=26, verbose=False)
        graphs.append(("G_10", OS.build_sparse(10, reps, eng)))
    except Exception as e: print("  [warn]", e)

    for name, A in graphs:
        V = A.shape[0]; sets = class_sets(A)
        pop = {}
        for s in sets:
            for c in s: pop[c] = pop.get(c, 0) + 1
        order = sorted(pop, key=lambda c: -pop[c])
        print("=" * 90)
        print(f"### Γ={name}  V={V}   类总数={len(pop)}   最热类的顶点数={pop[order[0]]}")
        print(f"{'K':>3} {'存活顶点':>8} {'畴数':>6} {'块大小分布(前6)':>34}")
        rows = []
        for K in (1, 2, 3, 5, 10, 20):
            keep = set(order[:K])
            alive = np.array([bool(sets[v] & keep) for v in range(V)])
            szs = comps(alive, A)
            rows.append({"K": K, "alive": int(alive.sum()), "clusters": len(szs),
                         "sizes": sorted(szs, reverse=True)[:8]})
            print(f"{K:>3} {int(alive.sum()):>8} {len(szs):>6} {str(sorted(szs, reverse=True)[:6]):>34}")
        out["runs"].append({"graph": name, "V": V, "n_classes": len(pop),
                            "hottest_class_pop": pop[order[0]], "rows": rows})
    with open(os.path.join(OUT, "z0_retain.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("=" * 90); print(f"用时 {round(time.time()-t0,1)}s → results/z0_retain.json")
    print("读法：畴数 > 1 ⟹ 真类标签让存活竞争**碎成了畴**（§11 用长度标签时恒为 1）。")


if __name__ == "__main__": main()
