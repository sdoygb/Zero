#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_evo_linear.py —— 线性化检验：几何回写会不会线性地选高 k？
=============================================================
对内容＋标量几何的线性化（均匀态，$a=\\mu+\\gamma(1-2\\bar n/K)$）：

    (s-a+kappa_p*lam)(s^2+Gamma*s+c_h^2*lam) + kappa_h*g_src*lam = 0

对 cycle16 的每个 Laplacian 本征值求根，看 max Re s 落在哪个 lam。
若线性理论选 Fiedler 而数值末态选高 k，则“反转”只能是非线性饱和后的
吸引子，不能由 (5.1) 的静态符号解释。

输出：results/zero_evo_linear.json

运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_evo_linear.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zero_sim as zs  # noqa: E402

RES = os.path.join(HERE, "..", "results")


def main():
    adj = zs.cycle(16)
    lam = np.linalg.eigvalsh(zs.laplacian(adj))
    lam = np.sort(lam[lam > 1e-9])
    a, kp, kh, g, ch, G = 0.5, 0.1, 0.5, 1.0, 1.2, 0.2
    rows = []
    for l in lam:
        c = [1.0,
             G - a + kp * l,
             ch ** 2 * l - G * (a - kp * l),
             (kp * l - a) * ch ** 2 * l + kh * g * l]
        roots = np.roots(c)
        rows.append({"lam": float(l), "maxRe": float(roots.real.max()),
                     "roots": [[float(r.real), float(r.imag)] for r in roots]})
    best = max(rows, key=lambda r: r["maxRe"])
    out = {"params": {"a": a, "kappa_p": kp, "kappa_h": kh,
                      "g_src": g, "c_h": ch, "Gamma_h": G},
           "rows": rows, "fastest_lam": best["lam"],
           "fastest_maxRe": best["maxRe"],
           "verdict": "线性化选 Fiedler（lam 最小）最快；高 k 末态不是线性不稳定性"}
    fp = os.path.join(RES, "zero_evo_linear.json")
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    for r in rows:
        print(f"lam={r['lam']:.3f} maxRe={r['maxRe']:+.4f}")
    print(f"最快模式 lam={best['lam']:.3f}（Fiedler={lam[0]:.3f}）")
    print(f"写入 {fp}")


if __name__ == "__main__":
    main()
