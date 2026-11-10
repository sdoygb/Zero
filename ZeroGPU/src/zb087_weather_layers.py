#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Zb087 · 天气系统的三层账：守恒数分层、边界通量定均值、零和振荡的偏斜。

七组判据（对应正文 §2–§4）：
  S1 守恒数一维     随机连通图上 dim ker Bᵀ = 1（Zb004 定理 4.2 的数值对照）。
  S2 分层给第二个  两层块对角 → dim ker = 2；两层的守恒互不混合。
  S3 TOA 恒等式     大气层 A ＋ 海洋层 O，外部"太空"节点 s：Δ(全局均值) = -Φ_TOA/|S|（Zb086 定理 2.2 的实例）。
  S4 ENSO 统计      ONI 918 季度：均值、偏度、峰高比。
  S5 符号置换检验   随机翻符号重算偏度：检验正偏不是抽样偶然。
  S6 弱连接谱隙     两块 5+4、2 条权 ε 的层间边：λ1/ε → k(1/|A|+1/|B|)，斜率 → 1（定理 2.3）。
  S7 非线性读出     Q=C0φ(1-αφ)、Q=±q：φ(+q)+φ(-q)=2αq²/C0²，偏度 ≈ 6ασq/C0（定理 4.2）。

数据：ONI（NOAA CPC，geocore/examples/data/oni.txt，1950–2026）
运行：python3 ZeroGPU/src/zb087_weather_layers.py
输出：ZeroGPU/results/zb087_weather_layers.json
"""
from __future__ import annotations

import json
import random
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results" / "zb087_weather_layers.json"
ONI = Path("/Volumes/download/from_system_20261106/A_archives/Zero260401/geocore/examples/data/oni.txt")
SEASONS = ["DJF", "JFM", "FMA", "MAM", "AMJ", "MJJ", "JJA", "JAS", "ASO", "SON", "OND", "NDJ"]


def _graph_b(n, seed):
    rng = random.Random(seed)
    edges = [(i, i + 1) for i in range(n - 1)]
    edges += [(rng.randrange(n), rng.randrange(n)) for _ in range(2)]
    edges = [(i, j) for i, j in edges if i != j]
    B = np.zeros((n, len(edges)))
    for e, (i, j) in enumerate(edges):
        B[i, e], B[j, e] = -1.0, +1.0
    return edges, B


def _nullity(M, tol=1e-9):
    s = np.linalg.svd(M, compute_uv=False)
    return int((s < tol).sum()) + (M.shape[1] - len(s))


def s1_conservation_dim(n=7, seed=21):
    _, B = _graph_b(n, seed)
    return {"n": n, "dim_ker_Bt": _nullity(B.T), "rank": int(np.linalg.matrix_rank(B))}


def s2_two_layers(n1=6, n2=5, seed=22):
    _, B1 = _graph_b(n1, seed)
    _, B2 = _graph_b(n2, seed + 1)
    B = np.zeros((n1 + n2, B1.shape[1] + B2.shape[1]))
    B[:n1, :B1.shape[1]] = B1
    B[n1:, B1.shape[1]:] = B2
    return {"dim_ker_Bt_blockdiag": _nullity(B.T),
            "expected": 2,
            "single_layer": _nullity(B1.T)}


def s3_toa_identity(nA=5, nO=4, trials=300, seed=23):
    """A 层 nA 个格子、O 层 nO 个格子、另有 1 个太空节点。

    S = A ∪ O（地球系统）；太空节点 s 只与 A 的最后两个格子相连（TOA 的两条缝）。
    """
    n = nA + nO + 1
    s = n - 1
    rng = random.Random(seed)
    edges = [(i, i + 1) for i in range(nA - 1)]              # A 内部
    edges += [(nA + i, nA + i + 1) for i in range(nO - 1)]   # O 内部
    edges += [(0, nA), (nA - 1, nA + 1)]                     # A—O 耦合
    edges += [(nA - 2, s), (nA - 1, s)]                      # TOA
    edges = [(i, j) for i, j in edges if i != j]
    B = np.zeros((n, len(edges)))
    for e, (i, j) in enumerate(edges):
        B[i, e], B[j, e] = -1.0, +1.0
    S = set(range(nA + nO))
    worst = 0.0
    for _ in range(trials):
        x = np.array([rng.randint(-4, 4) for _ in range(n)], dtype=float)
        f = np.array([rng.randint(-4, 4) for _ in edges], dtype=float)
        xn = x + B @ f
        db = xn[list(S)].mean() - x[list(S)].mean()
        phi = 0.0
        for e, (i, j) in enumerate(edges):
            if i in S and j not in S:
                phi += f[e]
            elif j in S and i not in S:
                phi -= f[e]
        worst = max(worst, abs(db + phi / len(S)))
    return {"n_A": nA, "n_O": nO, "edges": len(edges), "trials": trials,
            "max_abs_residual": float(worst), "pass": bool(worst < 1e-12)}


def load_oni():
    v = []
    for line in ONI.read_text().splitlines():
        p = line.split()
        if len(p) == 4 and p[0] in SEASONS:
            v.append(float(p[3]))
    return np.array(v)


def s4_enso(a):
    n = len(a)
    sd = float(a.std(ddof=1))
    return {"n": n, "mean": float(a.mean()), "sd": sd,
            "skew": float(((a - a.mean()) ** 3).mean() / sd ** 3),
            "max": float(a.max()), "min": float(a.min()),
            "peak_ratio": float(abs(a.max() / a.min()))}


def s5_sign_permutation(a, trials=20000, seed=24):
    rng = np.random.default_rng(seed)
    def sk(v):
        m, s = v.mean(), v.std(ddof=1)
        return ((v - m) ** 3).mean() / s ** 3
    obs = abs(float(sk(a)))
    cnt = 0
    for _ in range(trials):
        v = a * rng.choice([-1.0, 1.0], size=a.size)
        if abs(float(sk(v))) >= obs:
            cnt += 1
    return {"observed_abs_skew": obs, "trials": trials,
            "p_value": (cnt + 1) / (trials + 1)}


def s6_spectral_gap(seed=25):
    """定理 2.3：两块 5+4 节点、2 条权 ε 的层间边；λ1/ε → k(1/|A|+1/|B|)。

    窄窗 (1e-4..1e-2) 看渐近斜率；宽窗 (1e-3..1e-1) 看 O(ε) 修正。
    """
    def lam1(eps):
        n1, n2 = 5, 4
        n = n1 + n2
        A = np.zeros((n, n))
        for i in range(n1 - 1):
            A[i, i + 1] = A[i + 1, i] = 1.0
        for i in range(n2 - 1):
            A[n1 + i, n1 + i + 1] = A[n1 + i + 1, n1 + i] = 1.0
        for i, j in [(2, n1), (4, n1)]:
            A[i, j] = A[j, i] = eps
        L = np.diag(A.sum(1)) - A
        return float(np.linalg.eigvalsh(L)[1])

    narrow = np.geomspace(1e-4, 1e-2, 5)
    wide = np.geomspace(1e-3, 1e-1, 5)
    ln = np.array([lam1(e) for e in narrow])
    lw = np.array([lam1(e) for e in wide])
    k, nA, nB = 2.0, 5, 4
    theory = k * (1.0 / nA + 1.0 / nB)
    return {"lambda1_narrow": ln.tolist(), "eps_narrow": narrow.tolist(),
            "slope_narrow": float(np.polyfit(np.log(narrow), np.log(ln), 1)[0]),
            "lambda1_wide": lw.tolist(), "eps_wide": wide.tolist(),
            "slope_wide": float(np.polyfit(np.log(wide), np.log(lw), 1)[0]),
            "lambda1_over_eps_small": float(ln[0] / narrow[0]),
            "theory_coefficient": float(theory),
            "pass": bool(abs(ln[0] / narrow[0] - theory) < 0.02)}


def s7_nonlinear_readout(C0=1.0, alpha=0.02, n=4_000_000, seed=26):
    """定理 4.2：Q=C0φ(1-αφ)，零和内容 Q=±q 给正偏读出 φ。

    核对：精确反解的回代、φ(+q)+φ(-q)、均值、偏度、峰比。
    """
    rng = np.random.default_rng(seed)
    q = rng.normal(0.0, 1.0, n)                 # σ_q = 1
    y = 4.0 * alpha * q / C0
    if np.max(np.abs(y)) >= 1.0:
        raise RuntimeError("采样越出 |4αq/C0|<1 的读出分支")
    phi = (1.0 - np.sqrt(1.0 - y)) / (2.0 * alpha)
    mean = float(phi.mean())
    sd = float(phi.std(ddof=1))
    skew = float(((phi - mean) ** 3).mean() / sd ** 3)
    qg = np.array([0.5, 1.0, 2.0, 2.59])
    xp = (1.0 - np.sqrt(1.0 - 4.0 * alpha * qg / C0)) / (2.0 * alpha)
    xm = (1.0 - np.sqrt(1.0 + 4.0 * alpha * qg / C0)) / (2.0 * alpha)
    theta = np.linspace(-2.0, 2.0, 9)
    qsum = C0 * (-theta * (1.0 + alpha * theta) + theta * (1.0 - alpha * theta))
    return {"C0": C0, "alpha": alpha, "samples": n,
            "sigma_q": 1.0,
            "max_abs_y": float(np.max(np.abs(y))),
            "mean": mean, "pred_mean": float(alpha),
            "sd": sd,
            "skew": skew, "pred_skew_6a_sigma_over_C0": float(6.0 * alpha),
            "peak_ratio": (xp / np.abs(xm)).tolist(),
            "peak_ratio_approx": (1.0 + 2.0 * alpha * qg / C0).tolist(),
            "phi_plus_plus_phi_minus": (xp + xm).tolist(),
            "pred_sum_2a_q2": (2.0 * alpha * qg ** 2 / C0 ** 2).tolist(),
            "max_abs_Qsum_plus_2C0_a_theta2": float(np.max(np.abs(qsum + 2.0 * C0 * alpha * theta ** 2))),
            "max_abs_Q_roundtrip": float(np.max(np.abs(C0 * phi * (1.0 - alpha * phi) - q))),
            "pass": bool(skew > 0 and abs(skew - 6.0 * alpha) < 0.02 and
                         abs(mean - alpha) < 0.01)}


def main():
    a = load_oni()
    res = {
        "S1_conservation_dim": s1_conservation_dim(),
        "S2_two_layers": s2_two_layers(),
        "S3_toa_identity": s3_toa_identity(),
        "S4_enso": s4_enso(a),
        "S5_sign_permutation": s5_sign_permutation(a),
        "S6_spectral_gap": s6_spectral_gap(),
        "S7_nonlinear_readout": s7_nonlinear_readout(),
    }
    res["PASS"] = all([res["S1_conservation_dim"]["dim_ker_Bt"] == 1,
                       res["S2_two_layers"]["dim_ker_Bt_blockdiag"] == 2,
                       res["S3_toa_identity"]["pass"],
                       res["S6_spectral_gap"]["pass"],
                       res["S7_nonlinear_readout"]["pass"]])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(res, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
