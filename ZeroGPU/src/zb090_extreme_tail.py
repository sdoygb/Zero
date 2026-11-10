#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Zb090 · 极端事件的尾部：零模不进配分函数。

四组判据（对应正文 §2–§4）：
  S1 零模不进能量   高斯链上 φ→φ+c1：能量不变、异常量 max(φ-mean) 的分布逐位不变。
  S2 判据方向       加一个非均匀扰动 ψ⊥1：同一个分布显著改变（KS 距离大）。
  S3 尾部形状       K=12 链上 P(X>x) 在 log 坐标里对 x 与 x² 拟合：x² 更直 ⟹ 不是幂律。
  S4 有限尺寸       把 K 从 6 加到 24：尾部尺度随 K 变（有限系统，无热力学极限）。

运行：python3 ZeroGPU/src/zb090_extreme_tail.py
输出：ZeroGPU/results/zb090_extreme_tail.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results" / "zb090_extreme_tail.json"
N_SAMP, SEED = 200_000, 20261107


def laplacian_path(K):
    L = 2 * np.eye(K) - np.eye(K, k=1) - np.eye(K, k=-1)
    L[0, 0] = L[-1, -1] = 1.0
    return L


def sample_phi(K, beta=1.0, rng=None):
    """φ ~ exp(-β/2 φᵀLφ) 的零和涨落（去掉零模）。"""
    L = laplacian_path(K)
    w, V = np.linalg.eigh(L)
    mask = w > 1e-9
    scale = 1.0 / np.sqrt(beta * w[mask])
    z = rng.standard_normal((N_SAMP, int(mask.sum())))
    return (V[:, mask] * scale) @ z.T          # K × N_SAMP


def ks(a, b):
    a, b = np.sort(a), np.sort(b)
    grid = np.union1d(a, b)
    ca = np.searchsorted(a, grid, side="right") / a.size
    cb = np.searchsorted(b, grid, side="right") / b.size
    return float(np.max(np.abs(ca - cb)))


def main():
    rng = np.random.default_rng(SEED)
    K = 12
    phi = sample_phi(K, 1.0, rng)
    x0 = phi.max(axis=0) - phi.mean(axis=0)

    # S1：均匀平移（零模）不改变异常量
    phi_shift = phi + 3.7
    x1 = phi_shift.max(axis=0) - phi_shift.mean(axis=0)
    ks_shift = ks(x0, x1)

    # S2：非均匀扰动
    j = np.arange(K)
    psi = np.cos(2 * np.pi * j / K)
    psi = psi - psi.mean()
    phi_pert = phi + 0.8 * psi[:, None]
    x2 = phi_pert.max(axis=0) - phi_pert.mean(axis=0)
    ks_pert = ks(x0, x2)

    # S3：尾部形状（P(X>x) 对 x 与 x² 的线性拟合）
    xs = np.quantile(x0, [0.5, 0.7, 0.8, 0.9, 0.95, 0.99])
    tail = 1.0 - np.array([(x0 <= x).mean() for x in xs])
    lx = np.log(tail)
    A1 = np.vstack([xs, np.ones_like(xs)]).T
    r1 = float(np.corrcoef(np.linalg.lstsq(A1, lx, rcond=None)[0] @ A1.T, lx)[0, 1] ** 2)
    A2 = np.vstack([xs ** 2, np.ones_like(xs)]).T
    r2 = float(np.corrcoef(np.linalg.lstsq(A2, lx, rcond=None)[0] @ A2.T, lx)[0, 1] ** 2)

    # S4：有限尺寸
    sizes = {}
    for k in (6, 12, 24):
        p = sample_phi(k, 1.0, np.random.default_rng(SEED + k))
        xx = p.max(axis=0) - p.mean(axis=0)
        sizes[k] = {"q99": float(np.quantile(xx, 0.99)), "sd": float(xx.std())}

    res = {
        "K": K,
        "S1_zero_mode": {"ks_shift": ks_shift, "pass": bool(ks_shift < 0.01)},
        "S2_direction": {"ks_perturb": ks_pert, "changes_as_expected": bool(ks_pert > 0.05)},
        "S3_tail": {"r2_linear_in_x": r1, "r2_linear_in_x2": r2,
                    "gaussian_like": bool(r2 > r1)},
        "S4_finite_size": sizes,
    }
    res["PASS"] = all([res["S1_zero_mode"]["pass"], res["S2_direction"]["changes_as_expected"],
                       res["S3_tail"]["gaussian_like"]])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(res, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
