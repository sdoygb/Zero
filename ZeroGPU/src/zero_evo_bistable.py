#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_evo_bistable.py —— c_h=1.55 的 FID 支单独复核（10^6 步）
=============================================================
`zero_evo_long.py` 的四条长跑里，c_h=1.55 只跑了 random 初值（行波支）。
本脚本补同参数的 fiedler 初值（FID 支），确认双稳的两个支在 10^6 步内
都不衰减。输出：results/zero_evo_bistable_155.json

运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_evo_bistable.py
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from zero_evo import ZeroEvo, lap_matrix  # noqa: E402
from zero_evo_phase import base_params, classify  # noqa: E402
from zero_evo_long import lyap_local  # noqa: E402
import zero_sim as zs  # noqa: E402

RES = os.path.join(HERE, "..", "results")


def main():
    os.makedirs(RES, exist_ok=True)
    adj = zs.cycle(16)
    w, Q = np.linalg.eigh(lap_matrix(adj))
    f = Q[:, 1]
    f = f / np.linalg.norm(f) * 0.3
    e = ZeroEvo(adj, base_params(c_h=1.55), n0=1.0 + f)
    t0 = time.time()
    e.run(1_000_000, every=500)
    cls = classify(e, adj)
    lam = lyap_local(e, nsteps=20000, eps=1e-10)
    out = {"c_h": 1.55, "init": "fiedler", "steps": 10 ** 6,
           "class": cls, "lyapunov": lam,
           "windows": [
               {"mean": float(np.mean(e.hist["n_std"][i::5])),
                "std": float(np.std(e.hist["n_std"][i::5]))}
               for i in range(5)],
           "wall_s": round(time.time() - t0, 1)}
    fp = os.path.join(RES, "zero_evo_bistable_155.json")
    with open(fp, "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print(f"c_h=1.55 fiedler: {cls['tag']} lam={cls['lam']:.3f} "
          f"high={cls['high']:.3f} temp={cls['temporal']:.2e} lyap={lam:+.4f}")
    print(f"写入 {fp}（{out['wall_s']}s）")


if __name__ == "__main__":
    main()
