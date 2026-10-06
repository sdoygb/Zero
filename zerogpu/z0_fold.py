"""
z0_fold.py --- 把单调饱和换成**折叠**（非单调）：能出现"每次重生都不同"吗？

三种重播种律（同一个逐顶点局部保留框架，只差非线性类型）：
  A 单调饱和  seed = min(a+b, K)              （= §15，预期收敛）
  B 秩截断    stack=[b,a,c] 取**最大的两个**，最小者清零（不连续）
  C 折叠      seed = (a+b) mod (K+1)          （锯齿/折叠，标准混沌配方）
判据：逐代距离 d_g 不衰减 ＋ 长时间不重复 ⟹ 非周期（"每次不一样"）。
"""
from __future__ import annotations
import json, os
import numpy as np
from scipy.sparse import csr_matrix, identity, kron
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "results")

def torus_adj(D, L):
    r = np.arange(L); A1 = csr_matrix((np.ones(L), (r, (r+1) % L)), shape=(L, L)); A1 = A1 + A1.T
    A = A1
    for _ in range(D-1):
        A = kron(A, identity(L, format="csr")) + kron(identity(A.shape[0], format="csr"), A1)
    return csr_matrix(A)

def add_hubs(A, frac=.25, extra=3, seed=0):
    rng = np.random.default_rng(seed); V = A.shape[0]; A = A.tolil()
    for u in rng.choice(V, size=max(1, int(frac*V)), replace=False):
        for _ in range(extra):
            v = int(rng.integers(0, V))
            if v != u: A[u, v] = 1.0; A[v, u] = 1.0
    return csr_matrix(A)

def run(A, kind="A", gens=80, K=2):
    V = A.shape[0]; deg = np.asarray(A.sum(1)).ravel()
    tau = np.maximum(1, np.round(0.6*deg.sum()/np.maximum(deg, 1))).astype(np.int64)
    N = np.eye(V)*(1.0/V); clock = np.zeros(V, np.int64)
    a = np.zeros(V); b = np.zeros(V); cur = np.zeros(V)
    trace = []
    for g in range(gens):
        for _ in range(int(np.median(tau))):
            N = N @ A
            mx = N.max()
            if mx > 0: N /= mx
            dg = np.diag(N).copy()
            if dg.sum() > 0:
                cur += (dg > 0); N[np.arange(V), np.arange(V)] = 0.0
            clock += 1
            for v in np.nonzero(clock >= tau)[0]:
                if N[:, v].sum() > 0: N[:, v] = 0.0
                if kind == "A":
                    b[v], a[v] = a[v], cur[v]
                    s = min(a[v] + b[v], K)
                elif kind == "B":                      # ★ 秩截断（不连续）
                    st = np.array([b[v], a[v], cur[v]])
                    st.sort(); st[0] = 0.0             # 丢掉最小者
                    b[v], a[v] = st[1], st[2]
                    s = min(a[v] + b[v], K)
                else:                                  # ★ 折叠
                    b[v], a[v] = a[v], cur[v]
                    s = (a[v] + b[v]) % (K + 1)
                cur[v] = 0
                N[v, v] += s
                clock[v] = 0
        trace.append(np.maximum(a + b, 0).copy())
    return np.array(trace)

res = {}
for gname, A in (("torus_L8_regular", torus_adj(2, 8)), ("torus_L8_hubs", add_hubs(torus_adj(2, 8)))):
    for kind in ("A", "B", "C"):
        S = run(A, kind=kind, gens=80)
        d = np.array([np.linalg.norm(S[i+1]-S[i])/max(np.linalg.norm(S[i]), 1e-12)
                      for i in range(len(S)-1)])
        seen = {}; rep = None
        for i, s in enumerate(S):
            k = tuple(np.round(s, 6).tolist())
            if k in seen: rep = (seen[k], i); break
            seen[k] = i
        res[f"{gname}|{kind}"] = {"d_first5": np.round(d[:5], 4).tolist(),
                                  "d_last10_mean": float(d[-10:].mean()),
                                  "repeat": rep, "n_distinct": len(seen)}
        print(f"{gname:<18} {kind}: d_g 首5={np.round(d[:5],3)}  末10均值={d[-10:].mean():.3e}  "
              f"重复={rep}  不同状态数={len(seen)}/80")
json.dump(res, open(os.path.join(OUT, "z0_fold.json"), "w"), ensure_ascii=False, indent=1)
print("→ results/z0_fold.json")
print("读法：末10均值不趋于 0 且无重复 ⟹ 非周期（每次重生不同）。")
