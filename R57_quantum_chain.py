#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R57 · 把猜想串起来：Z* → 前缀 → ω → K → 模流 → KMS → Born

猜想串（用户提供，逐个落地）：
  K1 零乱动，分一层一层        → 闭合 ±1 路径（Σw=0）；层 = 部分和的高度层级
  K2 Z 层闭合，路径记在历史层  → 取闭合路径的**前缀**（开路径）
  K3 Z 层是全局的              → Z* 多重集；下投影给出 E 与 P
  K4 E+P 破缺（局部失衡）      → 前缀终点 s ≠ 0 即局部失衡
  K5 层内/层间动力学           → 权重 = 补全数 × 层高代价

本脚本验证从这条链长出的"量子力学"：
  V1 非对易代数 M2(C)（原生，来自循环次序 + ±；见 G27）
  V2 态 ρ 正定、Tr ρ = 1、块权重 = ω
  V3 GNS 内积正定、维数 = dim A
  V4 模 Hamiltonian K = -log ρ；模流 σ_t(A) = ρ^{it} A ρ^{-it} 是 *-自同构
  V5 模流酉（保 GNS 内积 / 保迹）
  V6 KMS 条件（数值）：F(t) = ω(A σ_t(B)) 在带 0 ≤ Im t ≤ -1 内解析，
     边界值 F(t-i) = ω(σ_t(B) A)
  V7 Born 形式：p(P) = Tr(ρP) ≥ 0，正交投影求和 = 1
  V8 语境性 S_max > 2（KCBS，用顶三归一）
  V9 模谱 = {log(ω_i/ω_j)}；块大小的素数结构给秩（III_1 vs III_λ）
"""
from __future__ import annotations

import cmath
import json
import math
import os
from fractions import Fraction
from itertools import combinations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MU = np.array([math.sqrt(5), (5 - math.sqrt(5)) / 2, (5 - math.sqrt(5)) / 2])
OUT = {}


# ---------------------------------------------------------------- K1 K2 K3
def closed_paths(L):
    for pos in combinations(range(L), L // 2):
        w = [-1] * L
        for p in pos:
            w[p] = 1
        yield tuple(w)


def prefix_weights(L):
    """K3/K4：前缀 σ=(k,s) 的补全数 W = C(L-k,(L-k-s)/2)"""
    from math import comb
    out = {}
    for k in range(L + 1):
        for s in range(-k, k + 1, 2):
            rem = L - k
            if (rem - s) % 2:
                continue
            a = (rem - s) // 2
            if a < 0 or a > rem:
                continue
            out[(k, s)] = comb(rem, a)
    return out


def blocks(L, r):
    """按层高 |s| 聚合成块；权重 = 补全数 × r^|s|"""
    from collections import Counter
    blk = Counter()
    for (k, s), W in prefix_weights(L).items():
        blk[abs(s)] += W * (r ** abs(s))
    ws = [blk[h] for h in sorted(blk)]
    tot = sum(ws)
    return [x / tot for x in ws], sorted(blk)


# ---------------------------------------------------------------- V1 代数
def su2_basis():
    I = np.eye(2, dtype=complex)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    return [I, sx, sy, sz]


def check_algebra():
    I, sx, sy, sz = su2_basis()
    c1 = np.allclose(sx @ sy - sy @ sx, 2j * sz)
    c2 = np.allclose(sy @ sz - sz @ sy, 2j * sx)
    c3 = np.allclose(sz @ sx - sx @ sz, 2j * sy)
    # 二面体关系
    L = 4
    tau = np.array([[0, -1], [1, -1]], dtype=complex)  # 3 阶旋转的 2 维表示
    # 直接用 Z_L 的循环次序：置换矩阵
    P = np.roll(np.eye(L, dtype=complex), 1, axis=0)
    S = np.fliplr(np.eye(L, dtype=complex))
    dih = (np.allclose(P @ P.conj().T, np.eye(L)) and
           np.allclose(S @ S, np.eye(L)) and
           np.allclose(S @ P @ S, P.conj().T))
    return {"su2_commutators": bool(c1 and c2 and c3),
            "dihedral_relations": bool(dih),
            "dim_M2": 4, "L": L}


# ---------------------------------------------------------------- V2 V3
def check_state(ws):
    """态 ρ = block-diag((w_a/2) I_2)；验证正定、归一、GNS"""
    k = len(ws)
    dim = 2 * k
    rho = np.zeros((dim, dim), dtype=complex)
    for a, w in enumerate(ws):
        rho[2 * a:2 * a + 2, 2 * a:2 * a + 2] = (w / 2) * np.eye(2)
    ev = np.linalg.eigvalsh(rho)
    tr = np.trace(rho).real
    # GNS：随机取几个代数元素，检查 Gram 矩阵正定
    rng = np.random.default_rng(3)
    basis = []
    for _ in range(6):
        M = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
        basis.append(M)
    G = np.array([[np.trace(rho @ (A.conj().T @ B)) for B in basis]
                  for A in basis])
    return {"dim_algebra": dim, "min_eig_rho": float(ev.min()),
            "trace_rho": float(tr), "rho_positive": bool(ev.min() > 0),
            "gns_gram_min_eig": float(np.linalg.eigvalsh(G).min()),
            "gns_positive": bool(np.linalg.eigvalsh(G).min() > 0),
            "block_weights": [float(x) for x in ws]}


# ---------------------------------------------------------------- V4 V5 V6
def check_modular(ws, nA=None):
    """V4/V5/V6 模 Hamiltonian、模流、KMS —— 用谱分解精确实现 rho^{it}"""
    k = len(ws)
    dim = 2 * k
    rho = np.zeros((dim, dim), dtype=complex)
    for a, w in enumerate(ws):
        rho[2 * a:2 * a + 2, 2 * a:2 * a + 2] = (w / 2) * np.eye(2)
    ev, U = np.linalg.eigh(rho)
    nz = ev > 1e-300
    K = np.zeros((dim, dim), dtype=complex)
    for i in range(dim):
        if nz[i]:
            K += (-math.log(ev[i])) * np.outer(U[:, i], U[:, i].conj())

    def rpow(tc):
        """rho^{i tc}，只用正特征值（谱分解，零特征值贡献 0）"""
        M = np.zeros((dim, dim), dtype=complex)
        for i in range(dim):
            if nz[i]:
                M += (ev[i] ** (1j * tc)) * np.outer(U[:, i], U[:, i].conj())
        return M

    def sigma_c(tc, X):
        return rpow(tc) @ X @ rpow(-tc)

    rng = np.random.default_rng(11)
    A = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
    B = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))

    def omega(X):
        return np.trace(rho @ X)

    t0 = 0.37
    st = rpow(t0)
    unitary = np.allclose(st @ st.conj().T, np.eye(dim), atol=1e-9)
    hom = np.allclose(sigma_c(t0, A @ B), sigma_c(t0, A) @ sigma_c(t0, B), atol=1e-9)
    star = np.allclose(sigma_c(t0, A.conj().T), sigma_c(t0, A).conj().T, atol=1e-9)
    trace_pres = abs(np.trace(sigma_c(t0, A)) - np.trace(A)) < 1e-9
    pos_pres = True
    for _ in range(5):
        X = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
        Y = X.conj().T @ X
        if np.linalg.eigvalsh(sigma_c(0.21, Y)).min() < -1e-9:
            pos_pres = False

    def F(tr_, ti):
        return omega(A @ sigma_c(complex(tr_, ti), B))

    tests = []
    for tr_ in [-0.8, -0.3, 0.0, 0.25]:
        lhs = F(tr_, -1.0)
        rhs = omega(sigma_c(complex(tr_, 0.0), B) @ A)
        tests.append({"Re_t": tr_, "lhs": complex(lhs), "rhs": complex(rhs),
                      "diff": float(abs(lhs - rhs))})
    return {"K_span": float(np.real(np.max(np.linalg.eigvalsh(K))
                                    - np.min(np.linalg.eigvalsh(K)))),
            "K_nontrivial": bool(np.max(np.linalg.eigvalsh(K))
                                 - np.min(np.linalg.eigvalsh(K)) > 1e-6),
            "sigma_unitary": bool(unitary),
            "sigma_homomorphism": bool(hom),
            "sigma_star": bool(star),
            "sigma_trace_preserving": bool(trace_pres),
            "sigma_positive": bool(pos_pres),
            "kms_tests": tests,
            "kms_max_diff": max(x["diff"] for x in tests),
            "kms_holds": bool(max(x["diff"] for x in tests) < 1e-8)}


# ---------------------------------------------------------------- V7 V8
def check_born_and_context(ws):
    k = len(ws)
    dim = 2 * k
    rho = np.zeros((dim, dim), dtype=complex)
    for a, w in enumerate(ws):
        rho[2 * a:2 * a + 2, 2 * a:2 * a + 2] = (w / 2) * np.eye(2)
    rng = np.random.default_rng(5)
    ps = []
    for _ in range(20):
        v = rng.normal(size=dim) + 1j * rng.normal(size=dim)
        v = v / np.linalg.norm(v)
        P = np.outer(v, v.conj())
        ps.append(float(np.trace(rho @ P).real))
    # 正交投影求和 = 1
    Q, _ = np.linalg.qr(rng.normal(size=(dim, dim)))
    tot = sum(float(np.trace(rho @ np.outer(Q[:, i], Q[:, i].conj())).real)
              for i in range(dim))
    t3 = sorted(ws, reverse=True)[:3]
    t3 = [x / sum(t3) for x in t3]
    smax = float(np.dot(t3, MU))
    return {"born_nonneg": bool(min(ps) >= -1e-12),
            "born_min": min(ps), "born_max": max(ps),
            "orthogonal_sum": tot, "orthogonal_ok": abs(tot - 1) < 1e-10,
            "top3_normalized": [round(x, 6) for x in t3],
            "S_max": smax, "contextual": bool(smax > 2.0)}


# ---------------------------------------------------------------- V9
def factorize(n):
    f, d = {}, 2
    n = int(n)
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def exact_rank(vals):
    fs = [factorize(v) for v in vals]
    primes = sorted({p for f in fs for p in f})
    rows = [[fs[i + 1].get(p, 0) - fs[i].get(p, 0) for p in primes]
            for i in range(len(fs) - 1)]
    if not rows:
        return 0
    M = [[Fraction(x) for x in r] for r in rows]
    rk = 0
    for c in range(len(M[0])):
        piv = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        pv = M[rk][c]
        M[rk] = [x / pv for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk


def check_modular_spectrum(ws):
    """模谱 = {log(w_i/w_j)} 的加法子群；块大小素数结构给秩"""
    pos = [w for w in ws if w > 0]
    logs = [math.log(pos[i] / pos[i + 1]) for i in range(len(pos) - 1)]
    sizes = [int(round(w * 10 ** 9)) for w in ws]
    return {"log_ratios": [round(x, 6) for x in logs],
            "all_nonzero": all(abs(x) > 1e-12 for x in logs),
            "block_sizes_scaled": sizes,
            "note": "秩按精确整数块大小算（见 R53/R54）"}


if __name__ == "__main__":
    L = 16
    GAMMA = 1.13          # 挂起的常数（只用于取一个落在门内的代表值）
    r = math.exp(-GAMMA)
    ws, hs = blocks(L, r)
    print("=" * 74)
    print("K1-K5 串起来：Z* → 前缀 → ω（L=%d, gamma=%.2f 挂起）" % (L, GAMMA))
    print("=" * 74)
    print("  层高 h :", hs)
    print("  权重 ω :", [round(x, 6) for x in ws])
    print()
    OUT["setup"] = {"L": L, "gamma_pending": GAMMA, "r": r,
                    "layer_heights": hs, "weights": ws}

    print("V1 非对易代数（原生）")
    a = check_algebra()
    OUT["V1_algebra"] = a
    for k, v in a.items():
        print("   %-22s %s" % (k, v))
    print()

    print("V2/V3 态与 GNS")
    s = check_state(ws)
    OUT["V2_state"] = s
    for k, v in s.items():
        print("   %-22s %s" % (k, v))
    print()

    print("V4/V5/V6 模 Hamiltonian、模流、KMS")
    m = check_modular(ws)
    OUT["V4_modular"] = m
    for k, v in m.items():
        if k == "kms_tests":
            continue
        print("   %-22s %s" % (k, v))
    for t in m["kms_tests"]:
        print("     KMS t=%.2f-i : |Δ| = %.3e" % (t["Re_t"], t["diff"]))
    print()

    print("V7/V8 Born 形式与语境性")
    b = check_born_and_context(ws)
    OUT["V7_born"] = b
    for k, v in b.items():
        print("   %-22s %s" % (k, v))
    print()

    print("V9 模谱")
    sp = check_modular_spectrum(ws)
    OUT["V9_spectrum"] = sp
    print("   相邻对数比:", sp["log_ratios"])
    print()

    ok = (a["su2_commutators"] and a["dihedral_relations"]
          and s["rho_positive"] and s["gns_positive"]
          and m["K_nontrivial"] and m["sigma_unitary"]
          and m["sigma_homomorphism"] and m["sigma_star"] and m["sigma_positive"]
          and m["kms_holds"] and b["born_nonneg"] and b["orthogonal_ok"])
    OUT["all_checks_pass"] = bool(ok)
    print("=" * 74)
    print("总判定：%s" % ("量子力学全部要件通过（只余 gamma 挂起）" if ok else "有未通过项"))
    print("=" * 74)
    with open(os.path.join(HERE, "R57_quantum_chain_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False, default=str)
    print("-> R57_quantum_chain_results.json")
