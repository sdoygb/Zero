#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""L1 振幅族：层振幅与全 simplex 的 XY 优化；per-parameter BP vs gradient-norm。

对象
  L1 词基上的 XY 相邻交换矩阵 A（平衡子空间）；精确基态 = Perron 向量 gs，
  λ_max = 最大本征值。

两族
  ① 层振幅族：|psi> = Σ_h sqrt(p_h) |h>，|h> = 层 h 的均匀叠加（H=L/2 个参数 p_h）。
     最优 = gs 在层常量子空间的投影；保真度 F_layer = Σ_h (Σ_{w∈h} gs_w)^2 / b_h。
  ② 全 simplex 族：任意归一化 psi（N 个振幅）；XY 能量是二次型 psi^T A psi，
     梯度 g=2Apsi；用幂迭代/梯度上升找 gs。

判据
  A1 全 simplex 幂迭代到 λ_max（相对差 <1e-6）且保真度 >0.9999。
  A2 随机 psi 处 per-parameter 梯度方差 ~ 1/N（log-log 斜率 ≈ -1）。
  A3 同一处梯度范数 O(1)（>0.5）——per-parameter BP 不等于 gradient-norm BP。
  A4 层振幅族的保真度 <0.98（层常量子空间装不下 Perron 向量）。
  A5 层振幅族能量比 <1（层振幅有缺口）。

输出：results/zero_l1_amplitude.json
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

OUT = ROOT / "results" / "zero_l1_amplitude.json"
OK, BAD, ROWS = [], [], []
RNG = np.random.default_rng(20261022)


def ck(name: str, cond: bool, detail: str = "") -> bool:
    (OK if cond else BAD).append(name)
    ROWS.append({"name": name, "ok": bool(cond), "detail": detail})
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def build(n):
    words = zs.balanced_words(n)
    N = len(words)
    index = {w: i for i, w in enumerate(words)}
    h = np.array([zs.layer_of(w) for w in words], dtype=int)
    A = np.zeros((N, N))
    for i, w in enumerate(words):
        for k in range(n - 1):
            if w[k] != w[k + 1]:
                w2 = list(w)
                w2[k], w2[k + 1] = w2[k + 1], w2[k]
                A[i, index[tuple(w2)]] = 1.0
    return words, N, h, A


def main():
    print("=" * 72)
    print("zero_l1_amplitude.py — L1 振幅族：层振幅 vs 全 simplex")
    print("=" * 72)
    rows = []
    for n in [8, 10, 12]:
        words, N, h, A = build(n)
        w, V = np.linalg.eigh(A)
        lam = float(w[-1])
        gs = V[:, -1]
        gs = gs / np.linalg.norm(gs)
        if gs.sum() < 0:
            gs = -gs
        # ① 层常量子空间投影
        F_layer = 0.0
        proj = np.zeros(N)
        for hh in range(1, n // 2 + 1):
            m = h == hh
            c = gs[m].sum() / np.sqrt(m.sum())
            F_layer += c * c
            proj[m] = c / np.sqrt(m.sum())
        proj /= np.linalg.norm(proj)
        E_layer = float(proj @ A @ proj)
        # ② 全 simplex 幂迭代
        psi = RNG.normal(size=N)
        psi /= np.linalg.norm(psi)
        # 邻接图是二分的（反转奇偶），A 有 ±λ_max；用平移幂迭代隔离 Perron
        B = A + lam * np.eye(N)
        for _ in range(300):
            psi = B @ psi
            psi /= np.linalg.norm(psi)
        E_full = float(psi @ A @ psi)
        fid_full = float(abs(psi @ gs) ** 2)
        # per-parameter 方差 / 梯度范数（随机 psi）
        vg, gn = [], []
        for _ in range(200):
            p = RNG.normal(size=N)
            p /= np.linalg.norm(p)
            g = 2.0 * (A @ p)
            g -= 2.0 * p * (p @ (A @ p))  # 球面上的投影梯度
            vg.append(float(np.mean(g * g)))
            gn.append(float(np.linalg.norm(g)))
        rows.append({"n": n, "N": N, "lambda_max": lam,
                     "E_full": E_full, "fid_full": fid_full,
                     "F_layer": float(F_layer), "E_layer": E_layer,
                     "ratio_layer": E_layer / lam,
                     "var_per_param": float(np.mean(vg)),
                     "grad_norm": float(np.mean(gn))})
        print(f"    n={n:2d} N={N:4d}  lam={lam:.4f}  full: E={E_full:.6f} fid={fid_full:.6f}"
              f"  layer: F={F_layer:.4f} E/lam={E_layer/lam:.4f}"
              f"  var/param={np.mean(vg):.3e}  |grad|={np.mean(gn):.3f}")

    ck("A1 full simplex reaches lambda_max",
       all(abs(r["E_full"] - r["lambda_max"]) / r["lambda_max"] < 1e-6 and
           r["fid_full"] > 0.9999 for r in rows),
       f"max rel={max(abs(r['E_full']-r['lambda_max'])/r['lambda_max'] for r in rows):.2e}")
    ln = np.log([r["N"] for r in rows]); lv = np.log([r["var_per_param"] for r in rows])
    slope = float(np.polyfit(ln, lv, 1)[0])
    ck("A2 per-parameter variance ~ 1/N", -1.3 < slope < -0.7, f"slope={slope:.3f}")
    ck("A3 gradient norm O(1)", all(r["grad_norm"] > 0.5 for r in rows),
       f"min={min(r['grad_norm'] for r in rows):.3f}")
    ck("A4 layer-constant fidelity < 0.98", all(r["F_layer"] < 0.98 for r in rows),
       f"max F={max(r['F_layer'] for r in rows):.4f}")
    ck("A5 layer energy ratio < 1", all(r["ratio_layer"] < 1.0 for r in rows),
       f"max ratio={max(r['ratio_layer'] for r in rows):.4f}")

    out = {"script": "zero_l1_amplitude.py", "checks": ROWS,
           "n_pass": len(OK), "n_fail": len(BAD), "rows": rows,
           "slope_var_vs_N": slope,
           "model": "L1 word simplex, XY adjacent-swap adjacency"}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n{len(OK)} PASS / {len(BAD)} FAIL  ->  {OUT}")
    return 0 if not BAD else 1


if __name__ == "__main__":
    raise SystemExit(main())
