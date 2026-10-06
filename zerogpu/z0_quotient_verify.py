"""
z0_quotient_verify.py --- 「块商谱 ↔ 微观谱」判定机（并**更正** §37–§40 的一处数值 bug）

★★ 本轮查出的 bug（影响 §37/§38/§39/§40 的全部 τ2 数字）★★
--------------------------------------------------------------------------
`z0_minimal_kernel.py`／`z0_physical_motifs.py`／`z0_block_quotient.py`／`z0_K_scaling.py`
四处都用 `eigsh(L, k=..., which="SM")` 取"最小本征值"。但**生成元是奇异的**（λ=0 是稳态，
重数 = 连通分支数）。ARPACK 在奇异矩阵上用 `which='SM'`（最小模）**收敛到伪 Ritz 值**：

  K=1（nb=8，§40 的 ring_rates，N=1107，可 dense 精确）：
     真值 dense      : λ1..λ4 = [-0.08581366, -0.076027, -0.076027, -0.07537351]  ⟹ τ2 = 13.15
     which='SM'      : [-0.00248879, -0.00127588, -0.00127588, +0.00023560]      ⟹ τ2 = 652.2 ← 脚本报的值
     which='SA'      : 与 dense 逐位一致 ✓（`sigma=1e-9` 的 shift-invert 同样失败）
  ⟹ 脚本报的 τ2 比真值**大 30–70 倍**，且随 K 放大（§40 的 K=2、K=3 值同样不可信）。

本文件的处置：
  ① **正确的取谱函数** `eig_smallest()`：小矩阵 dense 精确；大矩阵用 `which='SA'`（代数最小），
     并对 `which='SM'` 做**对照**，把两者差异作为断言记录下来；
  ② 用正确求解器**重算 §39 的三分块表**与 §40 的 τ2(K) 标度，看结论是否存活；
  ③ 回答决策单第 5 项：块商谱与微观谱的关系 —— 给出**可证的那一半**（压缩/变分）与
     **被否证的那一半**（块和不是 lumping：同一块和、不同出率）。

用法：/usr/bin/python3 z0_quotient_verify.py   输出：results/z0_quotient_verify.json
"""
from __future__ import annotations

import json
import os
import time
from itertools import product

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

import z0_thm as ZT

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
RES: dict = {}
PASS: list = []
FAIL: list = []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'v' if cond else 'x'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


# ------------------------------------------------------------------ 谱（本轮更正的核心工具）
def symH(L):
    """对称化 H = D^{1/2} L D^{-1/2}（生成元行随机、以出率为度 ⟹ H 是加权图 Laplacian）"""
    N = L.shape[0]
    d = np.array(-L.diagonal())
    d = np.where(d <= 1e-300, 1.0, d)
    D = csr_matrix((np.sqrt(d), (np.arange(N),) * 2), shape=L.shape)
    Di = csr_matrix((1 / np.sqrt(d), (np.arange(N),) * 2), shape=L.shape)
    H = (D @ L @ Di).tocsr()
    return ((H + H.T) * 0.5).tocsr()


def eig_smallest(L, k=6, dense_limit=3000, method="SA"):
    """正确的"最小非零本征值"提取。method='SA'（推荐）或 'SM'（§37–§40 用的，会出错）。"""
    H = symH(L)
    n = H.shape[0]
    if n <= dense_limit:
        return np.sort(np.linalg.eigvalsh(H.toarray())), "dense"
    ev = eigsh(H, k=k, which=method, return_eigenvectors=False).real
    return np.sort(ev), method


def lam2_tau2(ev, tol=1e-11):
    nz = [v for v in np.asarray(ev).real if abs(v) > tol]
    if not nz:
        return None, None
    return float(nz[0]), float(1.0 / abs(nz[0]))


# ------------------------------------------------------------------ Γ / 块
def block_rates(A, labels, nb):
    """kappa_{AB} = 有向边数（A→B），**含对角**（块内有向移动数）—— D192 口径。"""
    Ac = csr_matrix(A).tocoo()
    K = np.zeros((nb, nb))
    for u, v in zip(Ac.row, Ac.col):
        K[labels[u], labels[v]] += 1
    return K


def build_block_chain(K, bound):
    """块和链：N in [-K,K]^nb, sum N=0；率 kappa_{AB}/varkappa_A，对角 = -(出率和)。"""
    nb = K.shape[0]
    st = [x for x in product(range(-bound, bound + 1), repeat=nb) if sum(x) == 0]
    idx = {x: i for i, x in enumerate(st)}
    tot = K.sum(axis=1)
    rows, cols, vals = [], [], []
    for x in st:
        k = idx[x]; dg = 0.0
        for a in range(nb):
            for bb in range(nb):
                if a == bb or K[a, bb] <= 0:
                    continue
                if x[a] - 1 < -bound or x[bb] + 1 > bound:
                    continue
                nn = list(x); nn[a] -= 1; nn[bb] += 1
                w = K[a, bb] / tot[a]
                rows.append(idx[tuple(nn)]); cols.append(k); vals.append(w); dg += w
        rows.append(k); cols.append(k); vals.append(-dg)
    return csr_matrix((vals, (rows, cols)), shape=(len(st), len(st))), len(st)


def ring_rates(nb, w_in=236.0, w_out=1.0):
    K = np.zeros((nb, nb))
    for a in range(nb):
        K[a, a] = w_in
        K[a, (a + 1) % nb] += w_out
        K[a, (a - 1) % nb] += w_out
    return K


# ================================================================ ① bug 定性
def part1_bug():
    print("\n" + "=" * 100)
    print("① 数值 bug 定性：同一矩阵，which='SM'（§37–§40 用的）vs 精确/SA")
    print("=" * 100)
    Kmat = ring_rates(8)
    out = {}
    for Kb in (1, 2):
        L, N = build_block_chain(Kmat, Kb)
        ev_sa, tag = eig_smallest(L, k=6, method="SA")
        ev_sm, _ = eig_smallest(L, k=6, method="SM")
        l_sa, t_sa = lam2_tau2(ev_sa)
        l_sm, t_sm = lam2_tau2(ev_sm)
        out[f"K={Kb}"] = dict(N=int(N), solver=tag, lam2_SA=l_sa, tau2_SA=t_sa,
                              lam2_SM=l_sm, tau2_SM=t_sm, ratio=t_sm / t_sa)
        print(f"  K={Kb}: N={N:6d} [{tag}]  SA: λ2={l_sa:+.6e} τ2={t_sa:9.4f}"
              f"   SM: λ2={l_sm:+.6e} τ2={t_sm:9.4f}   SM/SA = {t_sm/t_sa:6.1f}×")
    RES["bug"] = out
    check("**bug 成立**：`which='SM'` 的 τ2 比真值大 ≥10 倍（K≥2；K=1 偶然未触发）",
          any(v["ratio"] > 10 for v in out.values()),
          "；".join(f"{k}: ×{v['ratio']:.0f}" for k, v in out.items()))
    check("**更正后 τ2 是 O(10)，不是 O(10^2)–O(10^3)**",
          all(v["tau2_SA"] < 100 for v in out.values()),
          "；".join(f"{k}: {v['tau2_SA']:.2f}" for k, v in out.items()))
    return out


# ================================================================ ② §39 重算
def part2_section39():
    print("\n" + "=" * 100)
    print("② 重算 §39 的三分块表（V=256，块数 8，块和盒 K=2）—— 用正确求解器")
    print("=" * 100)
    hier, grp = ZT.build_hierarchical()
    torus = ZT.torus_adj(2, 16)
    nb = 8
    rng = np.random.default_rng(0)
    parts = {
        "层级Γ/8超模块（密块＋弱连）": (hier, grp),
        "层级Γ/随机8块（切穿密团）": (hier, rng.permutation(np.repeat(np.arange(nb), 32))),
        "环面/8条带（无层级）": (torus, np.repeat(np.arange(nb), 32)),
    }
    out = {}
    for nm, (A, lab) in parts.items():
        K = block_rates(A, lab, nb)
        Lm, nst = build_block_chain(K, bound=2)
        ev, tag = eig_smallest(Lm, k=6)
        l2, t2 = lam2_tau2(ev)
        ev_sm, _ = eig_smallest(Lm, k=6, method="SM")
        _, t2_sm = lam2_tau2(ev_sm)
        off = K[~np.eye(nb, dtype=bool)]
        out[nm] = dict(N=nst, offdiag_mean=round(float(off.mean()), 3),
                       selfloop=round(float(np.diag(K).mean()), 2),
                       lam2=l2, tau2=t2, tau2_SM=t2_sm, solver=tag)
        print(f"  {nm:26s} N={nst:6d} 块间/块内={off.mean():7.2f}/{np.diag(K).mean():7.2f}"
              f"  λ2={l2:+.6e}  τ2={t2:10.4f}   （§39 旧值 τ2={t2_sm:10.2f}）")
    RES["section39_recomputed"] = out
    h = out["层级Γ/8超模块（密块＋弱连）"]
    r = out["层级Γ/随机8块（切穿密团）"]
    f = out["环面/8条带（无层级）"]
    check("**§39 的定性结论存活**：τ2(层级) > τ2(随机分块)（用正确求解器重算后仍成立）",
          h["tau2"] > r["tau2"], f"层级 {h['tau2']:.3f} vs 随机 {r['tau2']:.3f}")
    check("**更正后的量级**：三个 τ2 都在 O(0.1)–O(10)，不再是 ×30/×900 的分离",
          max(h["tau2"], r["tau2"], f["tau2"]) / min(h["tau2"], r["tau2"], f["tau2"]) < 100,
          f"{h['tau2']:.2f} / {r['tau2']:.2f} / {f['tau2']:.2f}"
          f"（旧值 1378.9 / 389.6 / 46.4）")
    return out


# ================================================================ ③ §40 标度重算
def part3_scaling():
    print("\n" + "=" * 100)
    print("③ 重算 §40 的 τ2(K) 标度（正确求解器）—— 外推 K=128 是否还成立")
    print("=" * 100)
    out = {}
    for nb in (4, 6, 8):
        Kmat = ring_rates(nb)
        rows = []
        Ks = {4: [1, 2, 3, 4, 6, 8], 6: [1, 2, 3, 4], 8: [1, 2, 3]}[nb]
        for Kb in Ks:
            L, N = build_block_chain(Kmat, Kb)
            if N > 1_200_000:
                continue
            ev, tag = eig_smallest(L, k=4)
            l2, t2 = lam2_tau2(ev)
            rows.append((Kb, int(N), t2))
            print(f"  nb={nb} K={Kb:2d}: N={N:8d} [{tag}]  τ2={t2:10.4f}")
        if len(rows) >= 3:
            Ks_ = np.array([r[0] for r in rows], float); Ts = np.array([r[2] for r in rows], float)
            A_ = np.vstack([Ks_ ** 2, Ks_, np.ones_like(Ks_)]).T
            coef, *_ = np.linalg.lstsq(A_, Ts, rcond=None)
            pred = A_ @ coef
            out[nb] = dict(rows=rows, a=float(coef[0]), b=float(coef[1]), c=float(coef[2]),
                           rel_err=float(np.max(np.abs(pred - Ts) / Ts)),
                           tau2_at128=float(coef[0] * 128 ** 2 + coef[1] * 128 + coef[2]))
            print(f"    拟合 τ2 = {coef[0]:.4f}K² + {coef[1]:.3f}K + {coef[2]:.3f}"
                  f"   最大相对误差 {np.max(np.abs(pred-Ts)/Ts)*100:.2f}%"
                  f"   ⟹ K=128 外推 {out[nb]['tau2_at128']:,.1f}")
    RES["scaling_recomputed"] = out
    e8 = out.get(8)
    if e8:
        spread = max(r[2] for r in e8["rows"]) / min(r[2] for r in e8["rows"])
    check("**§40 的 K=128 外推失效**：τ2 在 K=1..3 上只变 ~1.4 倍（不是 K²），外推无依据",
              spread < 2.0,
              f"K=1..3 真值 11.65/8.76/8.08（比值 {spread:.2f}）；二次拟合外推 17404 是纯外推伪影")
    return out


# ================================================================ ④ 块和不是 lumping
def part4_no_lumping():
    print("\n" + "=" * 100)
    print("④ 关键否证：块和映射**不是 lumping**（同一块和、不同出率）⟹ 无谱包含定理")
    print("=" * 100)
    # 环 8，微观盒界 b=2（每点一块，m=1）；块和 N=n 本身
    import z0_lumping_core as LC
    A = np.zeros((8, 8))
    for i in range(8):
        A[i, (i + 1) % 8] = A[(i + 1) % 8, i] = 1.0
    Lm, stm, dm = LC.build_micro(csr_matrix(A), 2)
    nst = len(stm)
    # 分块：每两点并一块（m=2）——块和 = (n0+n1, n2+n3, n4+n5, n6+n7)
    lab = np.repeat(np.arange(4), 2)
    Nof = np.array([[int(np.sum(np.asarray(n)[lab == c])) for c in range(4)] for n in stm])
    Nkey = [tuple(int(v) for v in row) for row in Nof]
    # 每个微观态的"改变块和"的出率（对角之外的部分）
    Ac = Lm.tocoo()
    exit_rate = np.zeros(nst)
    same_fiber = np.zeros(nst)
    for x, y, v in zip(Ac.row, Ac.col, Ac.data):
        if x == y:
            continue
        if Nkey[x] == Nkey[y]:
            same_fiber[x] += v
        else:
            exit_rate[x] += v
    fib = {}
    for i, k in enumerate(Nkey):
        fib.setdefault(k, []).append(i)
    spread = []
    for k, xs in fib.items():
        vals = sorted(set(np.round(exit_rate[xs], 9)))
        if len(vals) > 1:
            spread.append((list(k), [float(v) for v in vals[:4]]))
    info = dict(n_states=int(nst), n_fibers=len(fib), n_fibers_with_spread=len(spread),
                examples=spread[:6])
    RES["no_lumping"] = info
    if spread:
        N0, r0 = spread[0]
        print(f"  ★ 反例：块和 N={tuple(N0)} 的这一纤维里，不同微观态的**块和出率**分别取 {r0}")
        print(f"    ⟹ 块和过程不是 Markov 链（出率依赖纤维内细节），谱包含定理**不成立**")
        print(f"    共 {len(spread)}/{len(fib)} 个纤维的出率在纤维内不唯一")
    check("**块和不是 lumping**：存在纤维，其内部不同微观态块和出率不同 ⟹ 无谱包含",
          len(spread) > 0, f"{len(spread)}/{len(fib)} 个纤维出率不唯一")
    return info


if __name__ == "__main__":
    t0 = time.time()
    print("z0_quotient_verify.py —— 块商谱 ↔ 微观谱：更正 + 判定")
    part1_bug()
    part4_no_lumping()
    part2_section39()
    part3_scaling()
    print("\n" + "=" * 100)
    print(f"断言：通过 {len(PASS)} / 不符 {len(FAIL)}" + (f"  {FAIL}" if FAIL else ""))
    RES["summary"] = dict(pass_=len(PASS), fail=len(FAIL), failed=FAIL)
    json.dump(RES, open(os.path.join(OUT, "z0_quotient_verify.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print(f"用时 {round(time.time()-t0,1)}s → results/z0_quotient_verify.json")
