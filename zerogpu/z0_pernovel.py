"""
z0_pernovel.py --- 「每次毁灭重长都不一样」吗？测**周期性 vs 非周期性**

用户的判据：重生的状态**每次应当不同**。在确定性系统里这等价于**非周期/混沌**（Z0③ 只禁概率，不禁混沌）。

★ 改正：D222「只保留最高两层」是**每个局部历史层 P_i 各留各的**（per-vertex），
  我此前实现成两种**全局**版本（全局最近两代 / 全局统一封顶），把局部非线性抹平了。

本脚本：每个顶点 $v$ 保留**自己**最近两代闭合计数 $(a_v,b_v)$，重播种 $=\\min(a_v+b_v,K)$（饱和非线性）。
测 **逐代距离** $d_g=\\|s_{g+1}-s_g\\|/\\|s_g\\|$：→0 ⟹ 收敛（每次重生相同）；不衰减 ⟹ 非周期（每次不同）。
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

def run(A, gens=60, K=2, seed0=None):
    """逐顶点局部保留（per-vertex top-2）＋ 饱和重播种。返回每代的 (seed 向量, 活动层总数)。"""
    V = A.shape[0]; deg = np.asarray(A.sum(1)).ravel()
    tau = np.maximum(1, np.round(0.6*deg.sum()/np.maximum(deg, 1))).astype(np.int64)
    N = np.eye(V)*(1.0/V); clock = np.zeros(V, np.int64)
    a = np.zeros(V); b = np.zeros(V)          # ★ 每个顶点自己的最近两代
    cur = np.zeros(V)                          # 本轮的闭合计数
    trace = []
    for g in range(gens):
        for t in range(int(np.median(tau))):
            N = N @ A
            mx = N.max()
            if mx > 0: N /= mx
            dg = np.diag(N).copy()
            if dg.sum() > 0:
                cur += (dg > 0)
                N[np.arange(V), np.arange(V)] = 0.0
            clock += 1
            for v in np.nonzero(clock >= tau)[0]:
                col = N[:, v].copy()
                if col.sum() > 0: N[:, v] = 0.0
                b[v], a[v] = a[v], cur[v]        # ★ 各顶点自己的 top-2
                cur[v] = 0
                N[v, v] += min(a[v] + b[v], K)   # 饱和非线性：由保留的两层播种
                clock[v] = 0
        seedv = np.maximum(a + b, 0)
        trace.append((float(N.sum()), seedv.copy()))
    return trace

out = {}
for name, A in (("torus_L8_regular", torus_adj(2, 8)), ("torus_L8_hubs", add_hubs(torus_adj(2, 8)))):
    tr = run(A, gens=60)
    S = np.array([s for _, s in tr]); A_ = np.array([x for x, _ in tr])
    d = np.array([np.linalg.norm(S[i+1]-S[i])/max(np.linalg.norm(S[i]), 1e-12) for i in range(len(S)-1)])
    # 是否出现重复状态（周期性）
    seen = {}; first_rep = None
    for i, s in enumerate(S):
        key = tuple(np.round(s, 6).tolist())
        if key in seen: first_rep = (seen[key], i); break
        seen[key] = i
    print("=" * 84)
    print(f"### {name}  V={A.shape[0]}")
    print(f"    逐代距离 d_g：g=1..5 {np.round(d[:5],4)}")
    print(f"                  g=20..24 {np.round(d[20:25],4)}   g=50..54 {np.round(d[50:55],4)}")
    print(f"    → 收敛? d 末值={d[-1]:.3e}  （→0 ⟹ 每次重生相同；不衰减 ⟹ 每次不同）")
    print(f"    状态重复：{'第 '+str(first_rep)+' 代重复' if first_rep else '60 代内无重复 ⟹ 非周期'}")
    print(f"    活动层 A(t)：首={A_[0]:.3e} 末={A_[-1]:.3e}  相对起伏={A_.std()/A_.mean():.4f}")
    out[name] = {"d_first5": d[:5].tolist(), "d_last": float(d[-1]),
                 "repeat": first_rep, "act_first": float(A_[0]), "act_last": float(A_[-1]),
                 "act_rel_std": float(A_.std()/A_.mean())}
json.dump(out, open(os.path.join(OUT, "z0_pernovel.json"), "w"), ensure_ascii=False, indent=1)
print("=" * 84); print("→ results/z0_pernovel.json")
