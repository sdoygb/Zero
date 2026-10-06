"""
z0_cycle.py --- L2 毁灭周期**观察到没有**？测时间结构。

正则 Γ ⟹ 所有 τ_v 相等 ⟹ **全局同步** ⟹ 活动应当以周期 τ **相干振荡**；
非正则 Γ ⟹ τ 有展布 ⟹ 时钟**失相** ⟹ 振荡应当**衰减**成平稳态。
（上一轮只测了时间平均，且采样间隔≈周期，会混淆 —— 这次每步采样。）
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

def run(A, steps=1200):
    V = A.shape[0]; deg = np.asarray(A.sum(1)).ravel()
    tau = np.maximum(1, np.round(0.6*deg.sum()/np.maximum(deg, 1))).astype(np.int64)
    N = np.eye(V)*(1.0/V); clock = np.zeros(V, np.int64); record = np.zeros(V); term = np.zeros(V)
    act = np.zeros(steps); nfire = np.zeros(steps)
    for t in range(steps):
        N = N @ A
        mx = N.max()
        if mx > 0: N /= mx
        dg = np.diag(N).copy()
        if dg.sum() > 0:
            record += dg; N[np.arange(V), np.arange(V)] = 0.0
        clock += 1
        fire = np.nonzero(clock >= tau)[0]
        nfire[t] = len(fire)
        for v in fire:
            col = N[:, v].copy()
            if col.sum() > 0:
                term[v] += col.sum(); term -= col; N[:, v] = 0.0
            N[v, v] += min(record[v], 2)
            clock[v] = 0
        act[t] = N.sum()
    return act, nfire, tau

def acf(x, maxlag):
    x = x - x.mean(); d = (x*x).sum()
    return np.array([1.0] + [float((x[:-k]*x[k:]).sum()/d) for k in range(1, maxlag)])

out = {}
for name, A in (("torus_L8_regular", torus_adj(2, 8)), ("torus_L8_hubs", add_hubs(torus_adj(2, 8)))):
    act, nfire, tau = run(A, steps=1200)
    a = acf(act, 300)
    # 主频：对去均值后的 act 做 FFT
    F = np.abs(np.fft.rfft(act - act.mean()))
    f = np.fft.rfftfreq(len(act))
    k = int(np.argmax(F[1:])) + 1
    print("=" * 84)
    print(f"### {name}  V={A.shape[0]}  τ∈[{tau.min()},{tau.max()}]  同步={tau.min()==tau.max()}")
    print(f"    活动 A(t)：均值={act.mean():.3e}  相对起伏 std/mean={act.std()/act.mean():.4f}")
    print(f"    自相关：lag1={a[1]:.4f} lag10={a[10]:.4f} lag50={a[50]:.4f} "
          f"lag100={a[100]:.4f} lag200={a[200]:.4f}")
    print(f"    自相关主峰（lag>5）：lag={int(np.argmax(a[5:]))+5}  值={a[5:].max():.4f}"
          f"   （若=τ≈{int(np.mean(tau))} ⟹ **看到周期**）")
    print(f"    FFT 主频={f[k]:.5f} → 周期≈{1/f[k]:.1f} 步")
    out[name] = {"tau_min": int(tau.min()), "tau_max": int(tau.max()),
                 "act_mean": float(act.mean()), "act_rel_std": float(act.std()/act.mean()),
                 "acf_lag1": float(a[1]), "acf_lag10": float(a[10]), "acf_lag50": float(a[50]),
                 "acf_lag100": float(a[100]), "acf_lag200": float(a[200]),
                 "acf_peak_lag": int(np.argmax(a[5:]))+5, "acf_peak": float(a[5:].max()),
                 "fft_period": float(1/f[k]), "mean_tau": float(np.mean(tau))}
json.dump(out, open(os.path.join(OUT, "z0_cycle.json"), "w"), ensure_ascii=False, indent=1)
print("=" * 84); print("→ results/z0_cycle.json")
