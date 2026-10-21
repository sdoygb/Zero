#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""L1 层内权重：runs × 零点穿越，把层振幅的 XY 缺口压到 <1%。

原生内层坐标
  r(w) = 1 + #{k: w_k != w_{k+1}}      （Z075 的线性段数）
  z(w) = #{1 <= k < L: prefix_k = 0}   （走者回到原点）

族
  |psi> = Σ_h sqrt(p_h) * (1/sqrt(Z_h)) Σ_{w∈h} exp((mu_h r(w)+nu_h z(w))/2) |w>
  参数 p_h（softmax）、mu_h、nu_h；H=L/2；共 3H 个参数。

对照
  层振幅：mu=nu=0（只有 p_h）。
  runs：nu=0。
  runs+zc：mu,nu 都开。
  全 simplex：精确 Perron（eigh），1.000000。

判据
  L1 runs+zc E/λ_max > 0.995（n=8/10/12）。
  L2 runs+zc 保真度 > 0.985。
  L3 runs+zc E/λ > runs E/λ + 0.005。
  L4 runs+zc E/λ > 层振幅 E/λ + 0.01。
  L5 参数数 3H < N。
  L6 全 simplex E/λ = 1（eigh 精确）。

输出：results/zero_l1_layerweights.json
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

OUT = ROOT / "results" / "zero_l1_layerweights.json"
OK, BAD, ROWS = [], [], []
RNG = np.random.default_rng(20261023)


def ck(name: str, cond: bool, detail: str = "") -> bool:
    (OK if cond else BAD).append(name)
    ROWS.append({"name": name, "ok": bool(cond), "detail": detail})
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def build(n):
    words = zs.balanced_words(n)
    N = len(words)
    idx = {w: i for i, w in enumerate(words)}
    a = np.array([zs.word_area(w) for w in words], float)
    h = np.array([zs.layer_of(w) for w in words], int)
    r = np.array([zs.word_runs(w) for w in words], float)
    z = np.array([sum(1 for k in range(1, n) if sum(w[:k]) == 0) for w in words], float)
    ei, ej = [], []
    for i, w in enumerate(words):
        for k in range(n - 1):
            if w[k] != w[k + 1]:
                w2 = list(w)
                w2[k], w2[k + 1] = w2[k + 1], w2[k]
                ei.append(i); ej.append(idx[tuple(w2)])
    ei, ej = np.array(ei), np.array(ej)
    A = np.zeros((N, N))
    A[ei, ej] = 1.0
    return words, N, a, h, r, z, ei, ej, A


def main():
    print("=" * 72)
    print("zero_l1_layerweights.py — 层内 runs × 零点穿越")
    print("=" * 72)
    rows = []
    for n in [8, 10, 12]:
        words, N, a, h, r, z, ei, ej, A = build(n)
        H = n // 2
        w, V = np.linalg.eigh(A)
        lam = float(w[-1])
        gs = V[:, -1]
        gs /= np.linalg.norm(gs)
        if gs.sum() < 0:
            gs = -gs
        perron_res = float(np.linalg.norm(A @ gs - lam * gs))
        perron_min = float(gs.min())
        bh = np.array([(h == k).sum() for k in range(1, H + 1)], float)
        # 层振幅：精确投影
        F_layer = 0.0
        for k in range(1, H + 1):
            m = h == k
            F_layer += (gs[m].sum() / np.sqrt(m.sum())) ** 2
        # 优化器
        def make_state(x, use_runs, use_zc):
            u = x[:H]
            off = H
            mu = x[off:off + H] if use_runs else np.zeros(H)
            nu = x[off + H:off + 2 * H] if use_zc else np.zeros(H)
            pu = np.exp(u - u.max()); pu /= pu.sum()
            amp = np.zeros(N)
            for k in range(1, H + 1):
                m = h == k
                wgt = mu[k - 1] * r[m] + nu[k - 1] * z[m]
                e = np.exp(wgt - wgt.max()); e /= e.sum()
                amp[m] = np.sqrt(pu[k - 1] * e)
            return amp

        def optimize(use_runs, use_zc, steps=600, restarts=3, lr=0.02):
            dim = H + (H if use_runs else 0) + (H if use_zc else 0)
            def E(x):
                p = make_state(x, use_runs, use_zc)
                return float(np.sum(p[ei] * p[ej]))
            def F(x):
                p = make_state(x, use_runs, use_zc)
                return abs(p @ gs) ** 2
            best = (-1e9, None)
            for rr in range(restarts):
                x = np.zeros(dim)
                x[:H] = np.log(bh)
                for _ in range(500):
                    cand = x.copy()
                    cand[:H] += RNG.normal(0, 0.2, H)
                    if dim > H:
                        cand[H:] += RNG.normal(0, 0.3, dim - H)
                    if E(cand) > E(x):
                        x = cand
                m = np.zeros_like(x); v = np.zeros_like(x)
                eps = 1e-5
                for t in range(1, steps + 1):
                    g = np.zeros_like(x)
                    for i in range(dim):
                        x1 = x.copy(); x1[i] += eps
                        x2 = x.copy(); x2[i] -= eps
                        g[i] = (E(x1) - E(x2)) / (2 * eps)
                    m = 0.9 * m + 0.1 * g
                    v = 0.999 * v + 0.001 * g * g
                    x += lr * (m / (1 - 0.9 ** t)) / (np.sqrt(v / (1 - 0.999 ** t)) + 1e-12)
                    if t % 200 == 0:
                        lr *= 0.7
                e = E(x)
                if e > best[0]:
                    best = (e, x.copy())
            return best[0] / lam, F(best[1]), dim

        ratio_r, F_r, dim_r = optimize(True, False)
        ratio_z, F_z, dim_z = optimize(True, True)
        # 层振幅能量比（用投影态）
        proj = np.zeros(N)
        for k in range(1, H + 1):
            m = h == k
            c = gs[m].sum() / np.sqrt(m.sum())
            proj[m] = c / np.sqrt(m.sum())
        proj /= np.linalg.norm(proj)
        E_layer_ratio = float(proj @ A @ proj) / lam
        rows.append({"n": n, "N": N, "H": H, "lambda_max": lam,
                     "F_layer": float(F_layer), "E_layer_ratio": E_layer_ratio,
                     "E_runs_ratio": float(ratio_r), "F_runs": float(F_r),
                     "E_runszc_ratio": float(ratio_z), "F_runszc": float(F_z),
                     "dim_runs": dim_r, "dim_runszc": dim_z,
                     "perron_res": perron_res, "perron_min": perron_min})
        print(f"    n={n:2d} N={N:4d} H={H}  layer E={E_layer_ratio:.4f} F={F_layer:.4f}"
              f"  runs E={ratio_r:.4f} F={F_r:.4f}  runs+zc E={ratio_z:.6f} F={F_z:.6f}")

    ck("L1 runs+zc E/lam > 0.995", all(r["E_runszc_ratio"] > 0.995 for r in rows),
       f"min={min(r['E_runszc_ratio'] for r in rows):.6f}")
    ck("L2 runs+zc fidelity > 0.985", all(r["F_runszc"] > 0.985 for r in rows),
       f"min={min(r['F_runszc'] for r in rows):.6f}")
    ck("L3 runs+zc beats runs by >0.005",
       all(r["E_runszc_ratio"] > r["E_runs_ratio"] + 0.005 for r in rows),
       f"min diff={min(r['E_runszc_ratio'] - r['E_runs_ratio'] for r in rows):.6f}")
    ck("L4 runs+zc beats layer by >0.01",
       all(r["E_runszc_ratio"] > r["E_layer_ratio"] + 0.01 for r in rows),
       f"min diff={min(r['E_runszc_ratio'] - r['E_layer_ratio'] for r in rows):.6f}")
    ck("L5 3H < N", all(3 * r["H"] < r["N"] for r in rows),
       f"max 3H={max(3 * r['H'] for r in rows)}, min N={min(r['N'] for r in rows)}")
    ck("L6 Perron residual/positivity",
       all(r["perron_res"] < 1e-10 and r["perron_min"] > 0 for r in rows),
       f"max res={max(r['perron_res'] for r in rows):.2e}, min gs={min(r['perron_min'] for r in rows):.3f}")

    out = {"script": "zero_l1_layerweights.py", "checks": ROWS,
           "n_pass": len(OK), "n_fail": len(BAD), "rows": rows,
           "model": "L1 word layers, XY adjacent-swap adjacency"}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n{len(OK)} PASS / {len(BAD)} FAIL  ->  {OUT}")
    return 0 if not BAD else 1


if __name__ == "__main__":
    raise SystemExit(main())
