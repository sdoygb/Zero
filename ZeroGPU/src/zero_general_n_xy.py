#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""一般 N 的 XY 链：C4–C6（全局 3H vs 胞腔-wise）。

模型：固定 N 扇区的 XY 链 A（相邻 swap 的邻接矩阵）；Perron 基态 g>0。
振幅律：p_d ∝ exp(θ·f_d)
  3H：      层 one-hot + 每层 runs + 每层 zc（3H 参数）
  cell-wise：激发 rank one-hot（相对参考词：前 N 个 +1）+ 每层 runs + 每层 zc
参考词：w_ref = (+,...,+,-,...,-)，rank k = Hamming(w,w_ref)/2。

判据
  C4 一般 N（N≠L/2）：cell-wise 的 E/λ 或 F 优于 3H。
  C6 半满 N=L/2：cell-wise 不劣于 3H（容差 1e-3 E/λ）。
  C5 胞腔数 K=min(N,L-N)+1 随 L 多项式。

输出：results/zero_general_n_xy.json
"""
from __future__ import annotations

import itertools
import json
import pathlib
import sys

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

HERE = pathlib.Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import zero_sim as zs  # noqa: E402

OUT = HERE.parents[1] / "results" / "zero_general_n_xy.json"
OK, BAD, ROWS, TBL = [], [], [], []


def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    ROWS.append({"name": name, "ok": bool(cond), "detail": detail})
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(cond)


def words_general(L, N):
    out = []
    for pos in itertools.combinations(range(L), N):
        w = [-1] * L
        for i in pos:
            w[i] = 1
        out.append(tuple(w))
    return out


def feats(w):
    s = hi = lo = 0
    runs = 1
    zc = 0
    for i, x in enumerate(w):
        s += x
        hi = max(hi, s); lo = min(lo, s)
        if i > 0 and w[i] != w[i - 1]:
            runs += 1
        if 0 < i < len(w) - 1 and s == 0:
            zc += 1
    return max(hi, -lo), runs, zc


def build(L, N):
    words = words_general(L, N)
    idx = {w: i for i, w in enumerate(words)}
    n = len(words)
    A = np.zeros((n, n))
    for i, w in enumerate(words):
        for k in range(L - 1):
            if w[k] != w[k + 1]:
                w2 = list(w); w2[k], w2[k + 1] = w2[k + 1], w2[k]
                A[i, idx[tuple(w2)]] = 1.0
    val, vec = np.linalg.eigh(A)
    lam = float(val[-1]); g = vec[:, -1]
    g = g / np.linalg.norm(g)
    if g.sum() < 0:
        g = -g
    return words, A, lam, g


def mle(f, w, niter=3000):
    """max Σ w_d log p_d, p=softmax(f θ). 解析梯度。"""
    def obj(theta):
        z = f @ theta
        return logsumexp(z) - np.dot(w, z)
    def grad(theta):
        z = f @ theta
        p = np.exp(z - z.max()); p /= p.sum()
        return f.T @ (p - w)
    res = minimize(obj, np.zeros(f.shape[1]), jac=grad, method="L-BFGS-B",
                   options={"maxiter": niter, "ftol": 1e-15})
    z = f @ res.x
    p = np.exp(z - z.max()); p /= p.sum()
    return p, res.x


def main():
    print("=" * 74)
    print("zero_general_n_xy.py — 一般 N XY：3H vs 胞腔-wise（C4–C6）")
    print("=" * 74)
    cases = [(8, 3), (8, 4), (8, 5), (10, 4), (10, 5), (10, 6),
             (12, 5), (12, 6), (12, 7)]
    for L, N in cases:
        words, A, lam, g = build(L, N)
        n = len(words)
        w = g ** 2
        Hf = np.array([feats(x) for x in words])
        h = Hf[:, 0].astype(int); r = Hf[:, 1].astype(int); z = Hf[:, 2].astype(int)
        H = int(h.max()) + 1
        ref = (1,) * N + (-1,) * (L - N)
        k = np.array([sum(1 for a, b in zip(x, ref) if a != b) // 2 for x in words])
        K = int(k.max()) + 1
        # features
        f3 = np.zeros((n, 3 * H))
        fc = np.zeros((n, K + 2 * H))
        for i in range(n):
            f3[i, h[i]] = 1.0
            f3[i, H + h[i]] = r[i]
            f3[i, 2 * H + h[i]] = z[i]
            fc[i, k[i]] = 1.0
            fc[i, K + h[i]] = r[i]
            fc[i, K + H + h[i]] = z[i]
        flr = np.zeros((n, 3 * H + K))
        for i in range(n):
            flr[i, h[i]] = 1.0
            flr[i, H + k[i]] = 1.0
            flr[i, H + K + h[i]] = r[i]
            flr[i, H + K + H + h[i]] = z[i]
        base = 3 * H + K
        fsp = np.zeros((n, base + 2 * (L - 1)))
        fsp[:, :base] = flr
        for i, wd in enumerate(words):
            walls = [j for j in range(L - 1) if wd[j] != wd[j + 1]]
            fsp[i, base + walls[0]] = 1.0
            fsp[i, base + (L - 1) + walls[-1]] = 1.0
        fm = np.zeros((n, base + 2 * (L - 1) + 2 * L))
        fm[:, :base + 2 * (L - 1)] = fsp
        for i, wd in enumerate(words):
            walls = [j for j in range(L - 1) if wd[j] != wd[j + 1]]
            jbar = sum(walls) / len(walls)
            var = sum((j - jbar) ** 2 for j in walls) / len(walls)
            b1 = max(0, min(L - 2, int(round(jbar))))
            b2 = max(0, min(L - 1, int(round(var / max(1e-12, ((L - 2) / 2) ** 2) * (L - 1)))))
            fm[i, base + 2 * (L - 1) + b1] = 1.0
            fm[i, base + 2 * (L - 1) + L + b2] = 1.0
        p3, _ = mle(f3, w)
        pc, _ = mle(fc, w)
        plr, _ = mle(flr, w)
        psp, _ = mle(fsp, w)
        pm, _ = mle(fm, w)
        def energy(p):
            amp = np.sqrt(p)
            return float(amp @ A @ amp)
        def fid(p): return float(np.sum(g * np.sqrt(p)) ** 2)
        E3, Ec, Elr, Esp, Em = energy(p3), energy(pc), energy(plr), energy(psp), energy(pm)
        F3, Fc, Flr, Fsp, Fm = fid(p3), fid(pc), fid(plr), fid(psp), fid(pm)
        pmin = float(g.min()); pres = float(np.linalg.norm(A @ g - lam * g))
        TBL.append({"L": L, "N": N, "M": 2 * N - L, "n": n, "lam": lam,
                    "pmin": pmin, "pres": pres,
                    "H": H, "K": K, "params3": 3 * H, "paramsc": K + 2 * H,
                    "Eratio3": E3 / lam, "Eratioc": Ec / lam,
                    "F3": F3, "Fc": Fc,
                    "Elr": Elr / lam, "Flr": Flr, "paramslr": 3 * H + K,
                    "Esp": Esp / lam, "Fsp": Fsp, "paramssp": base + 2 * (L - 1),
                    "Em": Em / lam, "Fm": Fm, "paramsm": base + 2 * (L - 1) + 2 * L})
        print(f"  L={L} N={N} n={n:3d} 3H {3*H:2d}p E/λ={E3/lam:.5f} F={F3:.5f} | "
              f"cell {K+2*H:2d}p E/λ={Ec/lam:.5f} F={Fc:.5f} | "
              f"L+R {3*H+K:2d}p E/λ={Elr/lam:.5f} F={Flr:.5f} | "
              f"sp {base+2*(L-1):2d}p E/λ={Esp/lam:.5f} F={Fsp:.5f} | "
              f"mom {base+2*(L-1)+2*L:2d}p E/λ={Em/lam:.5f} F={Fm:.5f}  K={K}")
    # checks
    ck("X1 Perron positive/residual",
       all(x["pmin"] > 0 and x["pres"] < 1e-10 for x in TBL),
       f"min g={min(x['pmin'] for x in TBL):.3e}, max res={max(x['pres'] for x in TBL):.2e}")
    ck("C5 K=min(N,L-N)+1, poly", all(x["K"] == min(x["N"], x["L"] - x["N"]) + 1 for x in TBL),
       f"max K={max(x['K'] for x in TBL)}")
    gen = [x for x in TBL if x["N"] != x["L"] // 2]
    half = [x for x in TBL if x["N"] == x["L"] // 2]
    ck("C4 cell-wise beats 3H at N!=L/2 (E/λ)",
       all(x["Eratioc"] > x["Eratio3"] + 1e-4 for x in gen),
       f"min diff={min(x['Eratioc'] - x['Eratio3'] for x in gen):+.5f}")
    ck("C6 half-filled not worse (E/λ)",
       all(x["Eratioc"] >= x["Eratio3"] - 1e-3 for x in half),
       f"min diff={min(x['Eratioc'] - x['Eratio3'] for x in half):+.5f}")
    ck("C4b L+R beats 3H at N!=L/2 (E/λ and F)",
       all(x["Elr"] > x["Eratio3"] + 1e-4 and x["Flr"] > x["F3"] + 1e-4 for x in gen),
       f"min dE={min(x['Elr']-x['Eratio3'] for x in gen):+.5f}, min dF={min(x['Flr']-x['F3'] for x in gen):+.5f}")
    ck("C6b L+R half-filled not worse (E/λ)",
       all(x["Elr"] >= x["Eratio3"] - 1e-4 for x in half),
       f"min dE={min(x['Elr']-x['Eratio3'] for x in half):+.6f}")
    ck("C5b L+R params 3H+K, poly",
       all(x["paramslr"] == 3 * x["H"] + x["K"] for x in TBL),
       f"max={max(x['paramslr'] for x in TBL)}")
    ck("C4c spatial beats L+R at N!=L/2 (F)",
       all(x["Fsp"] > x["Flr"] + 1e-4 for x in gen),
       f"min dF={min(x['Fsp']-x['Flr'] for x in gen):+.5f}")
    ck("C6c spatial half-filled not worse (E/λ)",
       all(x["Esp"] >= x["Eratio3"] - 1e-4 for x in half),
       f"min dE={min(x['Esp']-x['Eratio3'] for x in half):+.6f}")
    ck("C5c spatial params base+2(L-1), poly",
       all(x["paramssp"] == 3 * x["H"] + x["K"] + 2 * (x["L"] - 1) for x in TBL),
       f"max={max(x['paramssp'] for x in TBL)}")
    ck("C4d moments beat first/last at N!=L/2 (F)",
       all(x["Fm"] > x["Fsp"] + 1e-4 for x in gen),
       f"min dF={min(x['Fm']-x['Fsp'] for x in gen):+.5f}")
    ck("C6d moments half-filled not worse (E/λ)",
       all(x["Em"] >= x["Eratio3"] - 1e-4 for x in half),
       f"min dE={min(x['Em']-x['Eratio3'] for x in half):+.6f}")
    ck("C5d moments params poly",
       all(x["paramsm"] == 3 * x["H"] + x["K"] + 2 * (x["L"] - 1) + 2 * x["L"] for x in TBL),
       f"max={max(x['paramsm'] for x in TBL)}")
    ck("T1 general-N F >= 0.99 (moments)",
       all(x["Fm"] >= 0.99 for x in gen),
       f"min={min(x['Fm'] for x in gen):.4f}")
    out = {"script": "zero_general_n_xy.py", "checks": ROWS, "table": TBL,
           "n_pass": len(OK), "n_fail": len(BAD)}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n{len(OK)} PASS / {len(BAD)} FAIL  ->  {OUT}")
    return 0 if not BAD else 1


if __name__ == "__main__":
    raise SystemExit(main())
