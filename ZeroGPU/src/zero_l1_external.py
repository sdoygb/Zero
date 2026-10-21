#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""L1 层分辨读出 vs 外部 XY 链的全局代价：BP 对照。

外部模型
  1D XY 链（n 比特，开边界），限制在平衡子空间 B_n = {b: Σ b_i = n/2}：
    H_XY = Σ_i J_i (X_i X_{i+1} + Y_i Y_{i+1})/2,  J_i = 1
  在计算基上，XX+YY 把 01 <-> 10（交换一个 1 与一个 0），矩阵元 1；
  Z 场是对角的，对面积相位的梯度没有贡献（Δa=0），故取纯 XY。

L1 嵌入
  平衡比特串 b 与 L1 词 w 一一对应：w_i = 2 b_i - 1（n=L）。
  w 的面积为 a(w)，高度为 h(w)；L1 层 b_h（Z161）。

代价与梯度
  |psi_theta> = N^{-1/2} Σ_w e^{i theta a(w)} |w>
  C(theta) = <psi_theta| H |psi_theta> = N^{-1} Σ_{w,w'} e^{i theta (a_w - a_{w'})} H_{w,w'}
  写 A_delta = N^{-1} Σ_{w,w': a_w-a_{w'}=delta} H_{w,w'}；
  Var_theta[dC/dtheta] = 2 Σ_{delta>0} delta^2 A_delta^2.

  全局代价：全部 H 的矩阵元。
  层分辨代价：P_h H P_h（两端都在层 h），层内归一化 N=b_h。

对照
  2-设计基线：同一 H_XY 在 n 比特全空间、Haar 态、局部参数 P=X_0；
    Var = ||[P,H_XY]||^2 / (d(d+1))。
  层投影：O_h = P_h H_XY P_h（P_h 是 L1 高度层 h 的 0/1 投影）；
    2-设计公式给 ||[P,O_h]||^2 / (d(d+1))。

判据
  E1 原生全局方差随 n 增长（不是 2^{-n}）。
  E2 2-设计全局方差按 2^{-n} 下降（log2 斜率 < -0.5）。
  E3 n=10 时原生全局 >> 2-设计全局（>1e2）。
  E4 n=10 时原生底层 >> 2-设计底层（>1e2）。
  E5 2-设计层投影仍是指数小（<1e-2，n=10）。
  E6 原生底层方差为正。

输出：results/zero_l1_external.json
"""
from __future__ import annotations

import json
import math
import pathlib
import sys
from collections import defaultdict

import numpy as np

HERE = pathlib.Path(__file__).resolve()
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
import zero_sim as zs  # noqa: E402

OUT = ROOT / "results" / "zero_l1_external.json"
OK, BAD, ROWS = [], [], []


def ck(name: str, cond: bool, detail: str = "") -> bool:
    (OK if cond else BAD).append(name)
    ROWS.append({"name": name, "ok": bool(cond), "detail": detail})
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def variance_from_H(H: np.ndarray, a: np.ndarray, N: int) -> float:
    """2 Σ_{δ>0} δ^2 A_δ^2，A_δ = N^{-1} Σ_{w,w': Δa=δ} H_{w,w'}."""
    acc = defaultdict(float)
    rows, cols = np.nonzero(np.triu(H, k=1))
    for i, j in zip(rows, cols):
        d = int(a[i] - a[j])
        if d == 0:
            continue
        acc[abs(d)] += H[i, j] / N
    return 2.0 * sum(d * d * v * v for d, v in acc.items())


def external_H(n: int, words, index):
    N = len(words)
    H = np.zeros((N, N))
    for i, w in enumerate(words):
        for k in range(n - 1):
            if w[k] != w[k + 1]:
                w2 = list(w)
                w2[k], w2[k + 1] = w2[k + 1], w2[k]
                j = index[tuple(w2)]
                H[i, j] += 1.0
    return H


def x0_full(n: int) -> np.ndarray:
    """X_0 在 2^n 维计算基上的矩阵（翻转第 0 位）。"""
    d = 1 << n
    P = np.zeros((d, d))
    for b in range(d):
        P[b, b ^ 1] = 1.0
    return P


def xy_full(n: int) -> np.ndarray:
    """1D XY 链（开边界，J=1）在 2^n 维计算基上的矩阵。"""
    d = 1 << n
    H = np.zeros((d, d))
    for b in range(d):
        for i in range(n - 1):
            bi = (b >> i) & 1
            bj = (b >> (i + 1)) & 1
            if bi != bj:
                bp = b ^ (1 << i) ^ (1 << (i + 1))
                H[b, bp] += 1.0
    return H


def comm_var(P: np.ndarray, O: np.ndarray, d: int) -> float:
    """2-设计纯态：Var[dC] = ||[P,O]||_HS^2 / (d(d+1))。"""
    C = P @ O - O @ P
    nrm = float(np.trace(C.conj().T @ C).real)
    return nrm / (d * (d + 1))


def layer_mask(n: int, h: int, words, a) -> np.ndarray:
    """L1 高度层 h 在 2^n 计算基上的 0/1 掩码。"""
    d = 1 << n
    mask = np.zeros(d)
    for w in words:
        if zs.layer_of(w) == h:
            b = sum((1 if x > 0 else 0) << i for i, x in enumerate(w))
            mask[b] = 1.0
    return mask

def main():
    print("=" * 72)
    print("zero_l1_external.py — L1 层分辨读出 vs 外部 XY 链的全局代价")
    print("=" * 72)
    table = []
    for n in range(4, 13, 2):
        words = zs.balanced_words(n)
        N = len(words)
        index = {w: i for i, w in enumerate(words)}
        a = np.array([zs.word_area(w) for w in words], dtype=np.int64)
        h = np.array([zs.layer_of(w) for w in words], dtype=np.int64)
        H = external_H(n, words, index)
        vg = variance_from_H(H, a, N)
        lay = []
        for hh in range(1, n // 2 + 1):
            m = h == hh
            bh = int(m.sum())
            if bh < 2:
                continue
            Hh = H[np.ix_(m, m)]
            vh = variance_from_H(Hh, a[m], bh)
            lay.append({"h": hh, "b_h": bh, "var": vh})
        table.append({"n": n, "N": N, "n_edges": int(H.sum() / 2),
                      "var_global": vg, "layers": lay,
                      "var_layer_sum": float(sum(x["var"] for x in lay))})
        print(f"    n={n:2d} N={N:4d} edges={int(H.sum()/2):4d}  Var_global={vg:.4e}"
              f"  layers: " + ", ".join(f"h{x['h']}({x['b_h']}):{x['var']:.3e}" for x in lay))

    # 2-设计基线：同一 XY 代价，Haar 2-设计态，局部参数 P=X_0
    design = []
    for n in range(4, 11, 2):
        d = 1 << n
        Hf = xy_full(n)
        P = x0_full(n)
        vg = comm_var(P, Hf, d)
        words = zs.balanced_words(n)
        a = np.array([zs.word_area(w) for w in words], dtype=np.int64)
        lay = []
        for hh in range(1, n // 2 + 1):
            m = layer_mask(n, hh, words, a)
            O = Hf * m[:, None] * m[None, :]
            lay.append({"h": hh, "b_h": int(m.sum()), "var": comm_var(P, O, d)})
        design.append({"n": n, "var_global": vg, "layers": lay})
        print(f"    [2d] n={n:2d} d={d:4d}  Var_global={vg:.4e}  "
              + ", ".join(f"h{x['h']}({x['b_h']}):{x['var']:.3e}" for x in lay))

    nat10 = [r for r in table if r["n"] == 10][0]
    des10 = [r for r in design if r["n"] == 10][0]
    ck("E1 native global grows with n", table[-1]["var_global"] > 2.0 * table[0]["var_global"],
       f"ratio={table[-1]['var_global']/table[0]['var_global']:.2f}")
    nn = np.array([r["n"] for r in design], float)
    lv = np.log2([r["var_global"] for r in design])
    slope = float(np.polyfit(nn, lv, 1)[0])
    ck("E2 2-design global decays (slope<-0.5)", slope < -0.5, f"slope={slope:.3f}")
    ratio_g = nat10["var_global"] / des10["var_global"]
    ck("E3 native global >> 2-design global (n=10)", ratio_g > 1e2, f"ratio={ratio_g:.2e}")
    ratio_l = nat10["layers"][0]["var"] / des10["layers"][0]["var"]
    ck("E4 native layer >> 2-design layer (n=10)", ratio_l > 1e2, f"ratio={ratio_l:.2e}")
    ck("E5 2-design layer still exponential (n=10)", des10["layers"][0]["var"] < 1e-2,
       f"{des10['layers'][0]['var']:.2e}")
    ck("E6 native layer positive", nat10["layers"][0]["var"] > 0,
       f"{nat10['layers'][0]['var']:.2e}")

    out = {"script": "zero_l1_external.py", "checks": ROWS,
           "n_pass": len(OK), "n_fail": len(BAD),
           "native": table, "design": design,
           "native_slope_log2": float(np.polyfit(
               np.array([r["n"] for r in table], float),
               np.log2([r["var_global"] for r in table]), 1)[0]),
           "design_slope_log2": slope, "ratio_global_n10": ratio_g,
           "ratio_layer_n10": ratio_l,
           "model": "1D XY chain; native = L1 area phase on balanced sector; "
                    "2-design = Haar state with local P=X_0"}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n{len(OK)} PASS / {len(BAD)} FAIL  ->  {OUT}")
    return 0 if not BAD else 1


if __name__ == "__main__":
    raise SystemExit(main())
