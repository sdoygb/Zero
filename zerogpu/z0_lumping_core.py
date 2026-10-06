"""
z0_lumping_core.py --- 核心判定：块和链是微观链的 **lumping（可集计）** 还是只是 **compression（压缩）**？

分割线（本文件要判定的事）
--------------------------
块和映射 pi: n -> N=(N_B = sum_{v in B} n_v)。微观盒界为 |n_v| <= b（b = 物理的 L）。
两条链：

  L_micro  状态 = {n in [-b,b]^V : sum n = 0}，每个允许补偿移动率 1
  L_blk    状态 = {N in [-K,K]^{nb} : sum N = 0}，kappa_{AB} = 有向块间边数，率 kappa_{AB}/K_A

**关键**：块和的可达范围是 |N_B| <= b * m（m = 块内位点数），**不是** K。
于是：
  · b*m <= K  ⟹ pi 的像含于块和盒 ⟹ 块和函数空间**前向不变** ⟹ **lumping** ⟹ spec(L_blk) ⊆ spec(L_micro)
  · b*m >  K  ⟹ pi 的像**超出**块和盒 ⟹ 不封闭 ⟹ 只是 **compression** ⟹ 无谱包含

本文件对每个 (Gamma, 分块, b, K) 逐点验三件事：
  (i)   lumping 恒等式 L_micro pi*f = pi*(L_blk f)（对落在盒内的那些态）
  (ii)  微观态中块和落进盒内的比例（=1 ⟺ 前向不变）
  (iii) 两条链的谱（对称化后对角化）：τ2 的比较与块谱是否 ⊆ 微观谱
"""
from __future__ import annotations
import sys
from itertools import product

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

sys.path.insert(0, "/Users/oygb/Downloads/lh/zerogpu")


# ------------------------------------------------------------------ Γ 基元
def two_K4(weak_edges=((0, 4),)):
    A = np.zeros((8, 8))
    for i in range(4):
        for j in range(4):
            if i < j:
                A[i, j] = A[j, i] = 1.0
    for i in range(4, 8):
        for j in range(4, 8):
            if i < j:
                A[i, j] = A[j, i] = 1.0
    for u, v in weak_edges:
        A[u, v] = A[v, u] = 1.0
    return A, np.array([0] * 4 + [1] * 4)


def star(n=8):
    A = np.zeros((n, n))
    for j in range(1, n):
        A[0, j] = A[j, 0] = 1.0
    return A, np.array([0] + [1] * (n - 1))


def ring(n=8):
    A = np.zeros((n, n))
    for i in range(n):
        A[i, (i + 1) % n] = A[(i + 1) % n, i] = 1.0
    return A, np.arange(n)


# ------------------------------------------------------------------ 生成元
def build_micro(A, bound):
    V = A.shape[0]
    st = [n for n in product(range(-bound, bound + 1), repeat=V) if sum(n) == 0]
    idx = {n: k for k, n in enumerate(st)}
    nbr = [np.nonzero(A[:, v])[0] for v in range(V)]
    rows, cols, vals = [], [], []
    diag = np.zeros(len(st))
    for k, n in enumerate(st):
        for a in range(V):
            if n[a] <= -bound:
                continue
            na = n[a]
            for b in nbr[a]:
                if n[b] >= bound:
                    continue
                nn = list(n); nn[a] = na - 1; nn[b] += 1
                rows.append(idx[tuple(nn)]); cols.append(k); vals.append(1.0)
                diag[k] += 1.0
    rows += list(range(len(st))); cols += list(range(len(st))); vals += list(-diag)
    return csr_matrix((vals, (rows, cols)), shape=(len(st), len(st))), st, diag


def block_kappa(A, lab, nb):
    Ac = csr_matrix(A).tocoo()
    K = np.zeros((nb, nb))
    for u, v in zip(Ac.row, Ac.col):
        K[lab[u], lab[v]] += 1
    return K


def build_block(K, bound):
    nb = K.shape[0]
    st = [x for x in product(range(-bound, bound + 1), repeat=nb) if sum(x) == 0]
    idx = {x: k for k, x in enumerate(st)}
    tot = K.sum(axis=1)
    rows, cols, vals = [], [], []
    for x in st:
        k = idx[x]; dg = 0.0
        for a in range(nb):
            for b in range(nb):
                if a == b or K[a, b] == 0:
                    continue
                if x[a] - 1 < -bound or x[b] + 1 > bound:
                    continue
                nn = list(x); nn[a] -= 1; nn[b] += 1
                w = K[a, b] / tot[a]
                rows.append(idx[tuple(nn)]); cols.append(k); vals.append(w); dg += w
        rows.append(k); cols.append(k); vals.append(-dg)
    return csr_matrix((vals, (rows, cols)), shape=(len(st), len(st))), st


def sym(L):
    d = np.array(-L.diagonal())
    d = np.where(d <= 1e-300, 1.0, d)
    D = csr_matrix((np.sqrt(d), (np.arange(len(d)),) * 2), shape=L.shape)
    Di = csr_matrix((1 / np.sqrt(d), (np.arange(len(d)),) * 2), shape=L.shape)
    H = (D @ L @ Di).tocsr()
    return ((H + H.T) * 0.5).tocsr()


def low(H, k=6, dense_limit=30000):
    if H.shape[0] <= dense_limit:
        return np.sort(np.linalg.eigvalsh(H.toarray()))
    return np.sort(eigsh(H, k=k, which="SM", return_eigenvectors=False).real)


def tau2(ev, tol=1e-8):
    nz = [v for v in np.sort(np.asarray(ev).real) if abs(v) > tol]
    return None if not nz else float(1 / abs(nz[0]))


def run(tag, A, lab, b, Kbox, report=True):
    A = csr_matrix(A)
    nb = len(np.unique(lab))
    V = A.shape[0]
    counts = np.array([int(np.sum(lab == c)) for c in np.unique(lab)])
    m = int(counts.max())
    K = block_kappa(A, lab, nb)
    Lm, stm, _ = build_micro(A, b)
    Lb, stb = build_block(K, Kbox)
    bidx = {x: i for i, x in enumerate(stb)}
    N = np.array([[int(np.sum(np.asarray(n)[lab == c])) for c in np.unique(lab)] for n in stm])
    inb = np.array([tuple(x) in bidx for x in N])
    frac = float(inb.mean())
    # (i) lumping 恒等式（只在落进盒内的态上可比）
    rng = np.random.default_rng(0)
    f = rng.standard_normal(len(stb))
    sub = np.array([bidx[tuple(x)] for x in N[inb]])
    # π*f 的完整微观向量：盒外取 0（外面对 block 链无定义）
    full = np.zeros(len(stm)); full[inb] = f[sub]
    lhs = Lm @ full
    rhs = np.zeros(len(stm)); rhs[inb] = (Lb @ f)[sub]
    dev_in = float(np.max(np.abs(lhs[inb] - rhs[inb]))) if inb.any() else np.nan
    dev_out = float(np.max(np.abs(lhs[~inb]))) if (~inb).any() else 0.0
    evm = low(sym(Lm)); evb = low(sym(Lb))
    t2m, t2b = tau2(evm), tau2(evb)
    dev_spec = [float(np.min(np.abs(evm - x))) for x in evb if abs(x) > 1e-8]
    out = dict(tag=tag, V=V, nb=nb, m=m, b=b, Kbox=Kbox, Kmax_reach=b * m,
               N_micro=int(Lm.shape[0]), N_block=int(Lb.shape[0]),
               frac_in_box=frac, lumping_dev_in=dev_in, lumping_dev_out=dev_out,
               tau2_micro=t2m, tau2_block=t2b, ratio=(t2m / t2b if t2b else None),
               spec_incl_max=float(np.max(dev_spec)) if dev_spec else None,
               ev_micro=[round(float(x), 8) for x in evm[:6]],
               ev_block=[round(float(x), 8) for x in evb[:6]])
    if report:
        print(f"\n  [{tag}] V={V} nb={nb} m={m} b={b} Kbox={Kbox}（块和可达 {b*m}）")
        print(f"     N_micro={Lm.shape[0]:7d}  N_block={Lb.shape[0]:7d}  落进盒内比例={frac:.4f}")
        print(f"     (i) 盒内 lumping 偏差={dev_in:.3e}   箱外残差={dev_out:.3e}")
        print(f"     (iii) τ2_micro={t2m:10.4f}  τ2_block={t2b:10.4f}  比={t2m/t2b:.4f}"
              f"   块谱⊆微观谱最大偏差={np.max(dev_spec) if dev_spec else 0:.3e}")
    return out


if __name__ == "__main__":
    print("=" * 100)
    print("z0_lumping_core.py —— 块和链：lumping（谱包含）还是 compression（无包含）？")
    print("=" * 100)
    A, lab = two_K4()
    rows = []
    print("\n### 两个 K4 弱连（模块对）—— 同一 Γ，只换微观盒界 b 与块和盒界 K")
    for b in (1, 2):
        for Kbox in (1, 2, 4):
            rows.append(run(f"2K4 b={b} K={Kbox}", A, lab, b, Kbox))
    print("\n### 星 K1,7（每块 1 点 ⟹ m=1 ⟹ 任何 b,K 下块和 = 点密度）")
    As, labs = star(8)
    for b in (1, 2):
        rows.append(run(f"star b={b} K={b}", As, labs, b, b))
    print("\n### 环 8（每块 1 点，同上）")
    Ar, labr = ring(8)
    for b in (1, 2):
        rows.append(run(f"ring b={b} K={b}", Ar, labr, b, b))
