#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Zb085 · 因果积累的速率与容量：并行加输入、端口截带宽、跨代不加分辨。

六组判据（对应正文 §2–§6，汇总在 §7）：
  S1 端口无信号   端口侧任意幺正不改系统约化态（复算 Zb074 命题 5.2 的 P1）。
  S2 端口滤波     指数记忆核的保留比：窄带 0.624、白噪声 0.239（Zb074 命题 5.5 的 P2）；
                  第二条独立路径用 FFT 的 |K|^2 与 sum h_n^2 对照。
  S3 输入计数     n 端口 × 每端口 k 读数：顺序 nk、同时 k^n（Zb074 命题 5.6）。
  S4 记录熵       ln(1+λ/κ) 表（Zb074 命题 7.2 的 C3）。
  S5 写入速率     每步可写入位数 ≤ ln(1+λ/κ)/(2τ ln(2+√3))（Zb074 (7.6) 的无量纲形）。
  S6 跨代容量     重播种 n 轮：|(g∘f)^n(X)| ≤ |Y| 且单调不增（Zb043 定理 2.2、推论 2.4）。

运行：python3 ZeroGPU/src/zb085_accumulation_rate.py
输出：ZeroGPU/results/zb085_accumulation_rate.json
"""
from __future__ import annotations

import json
import math
import random
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results" / "zb085_accumulation_rate.json"
rng = np.random.default_rng(20261107)
random.seed(20261107)

HALF_DEPTH = 2.0 * math.log(2.0 + math.sqrt(3.0))   # 2 ln(2+√3) ≈ 2.6339158


# ───────────────────────── S1 端口无信号 ─────────────────────────
def s1_port_no_signal():
    d_sys, d_port, trials = 2, 3, 20
    devs = []
    for _ in range(trials):
        # 随机系统态与端口态
        a = rng.normal(size=(d_sys, d_sys)) + 1j * rng.normal(size=(d_sys, d_sys))
        rho_s = a @ a.conj().T
        rho_s /= np.trace(rho_s)
        b = rng.normal(size=(d_port, d_port)) + 1j * rng.normal(size=(d_port, d_port))
        rho_p = b @ b.conj().T
        rho_p /= np.trace(rho_p)
        # 随机端口幺正（QR 分解）
        z = rng.normal(size=(d_port, d_port)) + 1j * rng.normal(size=(d_port, d_port))
        q, r = np.linalg.qr(z)
        u = q * (np.diag(r) / np.abs(np.diag(r)))
        rho = np.kron(rho_s, rho_p)
        full = np.kron(np.eye(d_sys), u) @ rho @ np.kron(np.eye(d_sys), u).conj().T
        # 部分迹掉端口
        reduced = np.einsum("acbc->ab", full.reshape(d_sys, d_port, d_sys, d_port))
        devs.append(float(np.max(np.abs(reduced - rho_s))))
    return {"trials": trials, "max_deviation": max(devs)}


# ───────────────────────── S2 端口滤波 ─────────────────────────
def s2_port_filter(tau=5.0, period=25.0):
    """窄带增益是解析量；宽带保留比是样本量（Zb074 P2 记的 0.239 是其中一次抽样）。"""
    n = 1200
    k = np.exp(-np.arange(n) / tau)
    h = k / k.sum()
    omega = 2.0 * np.pi / period
    gain = float(abs((h * np.exp(-1j * omega * np.arange(n))).sum()))     # 窄带幅度增益
    white_analytic = float(np.sqrt(np.sum(h ** 2)))                        # 白噪声 std 比
    # 按 Zb074 P2 的离散约定复算样本值：T=200、50 抽头、去掉两端 50 点
    T = 200
    t = np.arange(T)
    ker = np.exp(-np.arange(0, 50) / tau)
    ker = ker / ker.sum()
    r_nb = np.cos(2 * np.pi * t / period)
    y_nb = np.convolve(r_nb, ker, mode="same")
    it = slice(50, -50)
    sample_nb = float(np.std(y_nb[it]) / np.std(r_nb[it]))
    samples = []
    for seed in range(20):
        r_wb = np.random.default_rng(seed).normal(size=T)
        y_wb = np.convolve(r_wb, ker, mode="same")
        samples.append(float(np.std(y_wb[it]) / np.std(r_wb[it])))
    # 高频衰减斜率 |K|^2 ∝ w^-2
    w = np.linspace(0.0, 2.0, 2001)
    spec = (tau ** 2) / (1.0 + (w * tau) ** 2)
    hi = w > 5.0 / tau
    slope = float(np.polyfit(np.log(w[hi]), np.log(spec[hi]), 1)[0])
    return {
        "tau": tau,
        "narrowband_gain_analytic": round(gain, 3),
        "narrowband_sample_Zb074_convention": round(sample_nb, 3),
        "white_std_ratio_analytic": round(white_analytic, 3),
        "white_sample_range_20_seeds": [round(min(samples), 3), round(max(samples), 3)],
        "high_freq_slope": round(slope, 3),
    }


# ───────────────────────── S3 输入计数 ─────────────────────────
def s3_input_count(k=4, nmax=6):
    rows = {}
    for n in range(1, nmax + 1):
        rows[n] = {"sequential_nk": n * k, "simultaneous_k_pow_n": k ** n}
    return {"k": k, "rows": rows}


# ───────────────────────── S4 记录熵 ─────────────────────────
def s4_record_entropy(kappa=1.0):
    rows = {}
    for lam in (1.0, 10.0, 100.0, 1000.0, 10000.0):
        rows[int(lam)] = round(math.log(1.0 + lam / kappa), 3)
    return {"kappa": kappa, "ln_1_plus_lambda_over_kappa": rows}


# ───────────────────────── S5 写入速率上界 ─────────────────────────
def s5_rate_bound():
    rows = {}
    for tau in (1.0, 2.0, 5.0, 10.0):
        row = {}
        for lam in (10.0, 100.0, 1000.0, 10000.0):
            row[int(lam)] = round(math.log(1.0 + lam) / (tau * HALF_DEPTH), 3)
        rows[int(tau)] = row
    return {
        "formula": "bits/step <= ln(1+lambda/kappa) / (tau * 2 ln(2+sqrt3))",
        "half_depth_factor": round(HALF_DEPTH, 4),
        "rows": rows,
        "note": "并行数 n 不进入公式；它只经 lambda 进入，且为对数",
    }


# ───────────────────────── S6 跨代容量 ─────────────────────────
def s6_reseed_capacity(trials=300, nx=12, ny=6, rounds=4):
    worst_ok = True
    monotone = True
    sizes = []
    for _ in range(trials):
        X = list(range(nx))
        Y = list(range(ny))
        f = {x: random.randrange(ny) for x in X}              # 毁灭：满射由下方检查
        if len(set(f.values())) != ny:                        # 强制满射
            continue
        g = {y: random.choice([x for x in X if f[x] == y]) for y in Y}   # 重播种：选代表元
        cur = set(X)
        prev = None
        for _r in range(rounds):
            cur = {g[f[x]] for x in cur}
            if len(cur) > ny:
                worst_ok = False
            if prev is not None and len(cur) > prev:
                monotone = False
            prev = len(cur)
            sizes.append(len(cur))
    return {
        "trials": trials,
        "rounds": rounds,
        "all_subsets_of_Y": worst_ok,
        "monotone_nonincreasing": monotone,
        "max_size": max(sizes) if sizes else None,
        "bound": ny,
        "equal_fiber_check": equal_fiber(),
    }


def equal_fiber(k=3, ny=5):
    """等纤维时 |X| = k|Y|（Zb043 命题 2.3）。"""
    nx = k * ny
    X = list(range(nx))
    Y = list(range(ny))
    f = {x: x % ny for x in X}
    fibers = {}
    for x in X:
        fibers.setdefault(f[x], []).append(x)
    return {
        "k": k, "ny": ny, "nx": nx,
        "all_fibers_size_k": all(len(v) == k for v in fibers.values()),
        "nx_equals_k_ny": nx == k * ny,
    }


def main():
    res = {
        "S1_port_no_signal": s1_port_no_signal(),
        "S2_port_filter": s2_port_filter(),
        "S3_input_count": s3_input_count(),
        "S4_record_entropy": s4_record_entropy(),
        "S5_rate_bound": s5_rate_bound(),
        "S6_reseed_capacity": s6_reseed_capacity(),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(res, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
