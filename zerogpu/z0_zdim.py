"""
z0_zdim.py --- 判定 `Z4-EVIDENCE-SCOPE`：**物理通道图 Γ 上谱维数稳不稳？**
语料：Z0 §2.6 逼出终端款所依据的"无稳定谱维数"是在**位形空间图 𝒢_T** 上算的，
对**物理通道图 Γ** 是否成立 —— `R3` §7.3 第 5 条登记为**未判定**。
本脚本用已标定的 d_s 估计器（l2_regions.spectral_dim，环面上误差 ±5%~20%）对两族图比较。
判据：d_s 随尺寸**漂移** ⟹ 无稳定谱维数；**收敛** ⟹ 有稳定谱维数。
"""
import json, os
import numpy as np
from scipy.sparse import diags, csr_matrix, identity, kron
import l2_regions as LR

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")

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

def ds_of(A):
    deg = np.asarray(A.sum(1)).ravel()
    L = (diags(deg) - A).tocsr()
    ds, q, info = LR.spectral_dim(L)
    return ds, q

res = {}
print("=" * 88)
print("① 物理通道图候选族：ℤ² 环面（对照）与 环面＋枢纽")
for tag, maker in (("torus", lambda L: torus_adj(2, L)),
                   ("torus_hubs", lambda L: add_hubs(torus_adj(2, L), .25, 3, 0))):
    vals = []
    for L in (8, 12, 16, 24, 32):
        A = maker(L); ds, q = ds_of(A)
        vals.append((L, None if ds is None else round(ds, 4)))
        print(f"  {tag:<12} L={L:>3} N={A.shape[0]:>5}  d_s={vals[-1][1]}  usable={q.get('usable')}")
    good = [(L, d) for L, d in vals if d]
    res[tag] = {"vals": vals,
                "drift": None if len(good) < 2 else round(good[-1][1] - good[0][1], 4)}
    print(f"   ⟹ {tag} 的漂移（末−首）= {res[tag]['drift']}")

print("=" * 88)
print("② 位形空间图 𝒢_T（语料说它无稳定谱维数）")
try:
    import l0_closure as L0, observable_sweep as OS
    from zcl import Engine
    eng = Engine(); vals = []
    for T in (14, 16, 18, 20, 22):
        reps, _ = L0.enumerate_necklaces(T, eng, chunk_bits=26, verbose=False)
        A = OS.build_sparse(T, reps, eng); ds, q = ds_of(A)
        vals.append((T, None if ds is None else round(ds, 4)))
        print(f"  G_{T:<3} N={A.shape[0]:>7}  d_s={vals[-1][1]}  usable={q.get('usable')}")
    good = [(t, d) for t, d in vals if d]
    res["G_T"] = {"vals": vals, "drift": None if len(good) < 2 else round(good[-1][1]-good[0][1], 4)}
    print(f"   ⟹ 𝒢_T 的漂移（末−首）= {res['G_T']['drift']}")
except Exception as e:
    print("  [warn]", e)
json.dump(res, open(os.path.join(OUT, "z0_zdim.json"), "w"), ensure_ascii=False, indent=1)
print("=" * 88)
print("读法：漂移≈0 ⟹ 该族**有**稳定谱维数（Z4 的论证在它上面不成立，需重审）；")
print("      漂移显著 ⟹ **无**稳定谱维数（Z4 的论证成立）。")
