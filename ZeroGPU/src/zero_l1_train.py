#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""L1 层分辨族的训练：原生面积相位 vs 高频景观。

两组实验
  T1 原生态制备：目标 |psi_theta*> 是原生面积相位族的一个随机成员；
     随机搜索找盆地 + Adam 精化，训练保真度 1-|〈psi*|psi_theta〉|^2。
  T2 原生 XY 能量：纯铁磁 XY 链（平衡子空间）的基态是均匀叠加；
     随机搜索 + Adam 最大化能量，E/E_max -> 1，保真度 -> 1。
  T5/T2 的梯度范数在随机初值处为 O(1)（BP-free 的梯度信号）。

实现：相位矩阵 M[i,h]=a(w_i) 若 h(w_i)=h，则 phi = M theta；
      所有代价/梯度用矩阵-向量运算，随机搜索 8000 点/重启。

输出：results/zero_l1_train.json
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve()
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
import zero_sim as zs  # noqa: E402

OUT = ROOT / "results" / "zero_l1_train.json"
OK, BAD, ROWS = [], [], []
RNG = np.random.default_rng(20261021)


def ck(name: str, cond: bool, detail: str = "") -> bool:
    (OK if cond else BAD).append(name)
    ROWS.append({"name": name, "ok": bool(cond), "detail": detail})
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def build(n):
    words = zs.balanced_words(n)
    N = len(words)
    H = n // 2
    a = np.array([zs.word_area(w) for w in words], dtype=float)
    h = np.array([zs.layer_of(w) for w in words], dtype=int)
    M = np.zeros((N, H))
    for i in range(N):
        M[i, h[i] - 1] = a[i]
    return words, N, a, h, M


def edges(n, words, index):
    ei, ej = [], []
    for i, w in enumerate(words):
        for k in range(n - 1):
            if w[k] != w[k + 1]:
                w2 = list(w)
                w2[k], w2[k + 1] = w2[k + 1], w2[k]
                ei.append(i); ej.append(index[tuple(w2)])
    return np.array(ei), np.array(ej)


def state_prep(n, steps=800, restarts=5, lr=0.05):
    words, N, a, h, M = build(n)
    H = M.shape[1]
    tgt = RNG.uniform(0, 2 * np.pi, H)
    phis = M @ tgt
    def cost_grad(theta):
        ph = M @ theta
        z = np.exp(1j * (ph - phis))
        S = z.mean()
        F = abs(S) ** 2
        g = 2.0 * np.real(np.conj(S) * (1j / N) * (M.T @ z))
        return 1.0 - F, -g
    best = (1e9, None)
    gnorm = []
    for r in range(restarts):
        th0 = RNG.uniform(0, 2 * np.pi, H)
        c0, g0 = cost_grad(th0)
        gnorm.append(float(np.linalg.norm(g0)))
        # 随机搜索：面积相位是高频景观，先找盆地
        cand = RNG.uniform(0, 2 * np.pi, (8000, H))
        ph = cand @ M.T
        z = np.exp(1j * (ph - phis[None, :]))
        costs = 1.0 - np.abs(z.mean(axis=1)) ** 2
        k = int(np.argmin(costs))
        th = cand[k].copy()
        m = np.zeros(H); v = np.zeros(H)
        for t in range(1, steps + 1):
            c, g = cost_grad(th)
            m = 0.9 * m + 0.1 * g; v = 0.999 * v + 0.001 * g * g
            th -= lr * (m / (1 - 0.9 ** t)) / (np.sqrt(v / (1 - 0.999 ** t)) + 1e-12)
            if t % 200 == 0:
                lr *= 0.7
        c, _ = cost_grad(th)
        if c < best[0]:
            best = (c, th.copy())
    return {"n": n, "H": H, "F_final": 1.0 - best[0],
            "grad_norm_init": float(np.mean(gnorm)), "theta_star": tgt.tolist()}


def xy_energy(n, steps=800, restarts=5, lr=0.05):
    words, N, a, h, M = build(n)
    H = M.shape[1]
    index = {w: i for i, w in enumerate(words)}
    ei, ej = edges(n, words, index)
    A = np.zeros((N, N))
    A[ei, ej] = 1.0
    emax = float(np.linalg.eigvalsh(A)[-1])  # 相邻交换 XY 的最大本征值
    def energy_grad(theta):
        ph = M @ theta
        d = ph[ei] - ph[ej]
        E = float(np.sum(np.cos(d)) / N)
        g = ((M[ei] - M[ej]).T @ np.sin(d)) / N
        return -E, g
    best = (1e9, None)
    gnorm = []
    for r in range(restarts):
        th0 = RNG.uniform(0, 2 * np.pi, H)
        c0, g0 = energy_grad(th0)
        gnorm.append(float(np.linalg.norm(g0)))
        cand = RNG.uniform(0, 2 * np.pi, (8000, H))
        ph = cand @ M.T
        d = ph[:, ei] - ph[:, ej]
        costs = -np.sum(np.cos(d), axis=1) / N
        k = int(np.argmin(costs))
        th = cand[k].copy()
        m = np.zeros(H); v = np.zeros(H)
        for t in range(1, steps + 1):
            c, g = energy_grad(th)
            m = 0.9 * m + 0.1 * g; v = 0.999 * v + 0.001 * g * g
            th -= lr * (m / (1 - 0.9 ** t)) / (np.sqrt(v / (1 - 0.999 ** t)) + 1e-12)
            if t % 200 == 0:
                lr *= 0.7
        c, _ = energy_grad(th)
        if c < best[0]:
            best = (c, th.copy())
    ph = M @ best[1]
    fid = float(abs(np.mean(np.exp(1j * ph))) ** 2)
    return {"n": n, "H": H, "E_final": -best[0], "E_max": emax,
            "ratio": float(-best[0] / emax), "fid_uniform": fid,
            "grad_norm_init": float(np.mean(gnorm))}


def main():
    print("=" * 72)
    print("zero_l1_train.py — L1 层分辨族的训练（随机搜索 + Adam）")
    print("=" * 72)
    prep, xy = [], []
    for n in [8, 10, 12]:
        p = state_prep(n)
        x = xy_energy(n)
        prep.append(p); xy.append(x)
        print(f"    n={n:2d} H={p['H']}  prep F={p['F_final']:.4f} grad0={p['grad_norm_init']:.3f}"
              f"  XY E/Emax={x['ratio']:.4f} fid={x['fid_uniform']:.4f} grad0={x['grad_norm_init']:.3f}")

    ck("T1 state-prep F > 0.85", all(p["F_final"] > 0.85 for p in prep),
       f"min F={min(p['F_final'] for p in prep):.4f}")
    ck("T2 state-prep grad O(1)", all(p["grad_norm_init"] > 0.05 for p in prep),
       f"min grad={min(p['grad_norm_init'] for p in prep):.3f}")
    ck("T3 XY phase-family ratio > 0.8", all(x["ratio"] > 0.8 for x in xy),
       f"min ratio={min(x['ratio'] for x in xy):.4f}")
    ck("T4 XY fidelity to uniform > 0.5", all(x["fid_uniform"] > 0.5 for x in xy),
       f"min fid={min(x['fid_uniform'] for x in xy):.4f}")
    ck("T5 XY grad O(1)", all(x["grad_norm_init"] > 0.1 for x in xy),
       f"min grad={min(x['grad_norm_init'] for x in xy):.3f}")

    out = {"script": "zero_l1_train.py", "checks": ROWS,
           "n_pass": len(OK), "n_fail": len(BAD),
           "state_prep": prep, "xy_energy": xy}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n{len(OK)} PASS / {len(BAD)} FAIL  ->  {OUT}")
    return 0 if not BAD else 1


if __name__ == "__main__":
    raise SystemExit(main())
