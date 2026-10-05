#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R35_type_iii_probe.py -- P2″ 探针：极限因子的类型由"年龄块权重的差集所生成的加法群"决定。

判据（本文核验并数值实现）：
  设块权重 w_0..w_T > 0，模谱 = {log(w_i/w_j)}，其生成的加法子群 G ⊆ R：
    · G = {0}                 -> 模流（渐近）平凡（R12 的常数剖面）
    · G ≅ cZ（循环/离散）      -> type III_lambda，lambda = e^{-c}
    · G 在 R 中稠密            -> type III_1（BW 的几何模流所需）
  数值判据：取非零差集 D，令 c* = min|D|，若 max|d/c* - round(d/c*)| 很小 -> 循环；否则 -> 稠密。

输出：R35_type_iii_results.json
"""

from __future__ import annotations

import io
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def differences(logw: np.ndarray, window: float = 30.0, cap: int = 4000):
    vals = []
    n = len(logw)
    for i in range(n):
        for j in range(i + 1, n):
            d = abs(logw[i] - logw[j])
            if 1e-12 < d <= window:
                vals.append(d)
    vals = sorted(set(round(v, 12) for v in vals))
    if len(vals) > cap:
        idx = np.linspace(0, len(vals) - 1, cap).astype(int)
        vals = [vals[k] for k in idx]
    return np.array(vals)


def classify(name: str, w: np.ndarray, note: str = ""):
    w = np.asarray(w, dtype=float)
    w = w / w.max()
    logw = np.log(w)
    logw = logw - logw[0]
    span = float(np.ptp(logw))
    if span < 1e-9:
        return {"name": name, "note": note, "n_blocks": int(len(w)),
                "log_span": span, "type": "trivial ({0})", "lambda": 1.0,
                "cyclic_residual": 0.0, "n_diffs": 0}
    D = differences(logw)
    if len(D) == 0:
        return {"name": name, "note": note, "n_blocks": int(len(w)),
                "log_span": span, "type": "undecided", "lambda": None,
                "cyclic_residual": None, "n_diffs": 0}
    c = float(D.min())
    resid = float(np.abs(D / c - np.round(D / c)).max())
    cyclic = resid < 1e-3
    lam = math.exp(-c) if cyclic else None
    return {"name": name, "note": note, "n_blocks": int(len(w)),
            "log_span": span, "c_min": c, "cyclic_residual": resid,
            "n_diffs": int(len(D)),
            "type": ("III_lambda" if cyclic else "III_1 (dense)"),
            "lambda": lam}


T = 48
profiles = [
    ("uniform", np.ones(T + 1), "常数剖面：R12 的机制"),
    ("native_branching_2^a", 2.0 ** np.arange(T + 1), "原生满分支：w_a ∝ 2^a"),
    ("geometric_0.99^a", 0.99 ** np.arange(T + 1), "缓慢几何"),
    ("polynomial_(a+1)^2", (np.arange(T + 1) + 1.0) ** 2, "幂律"),
    ("factorial_(a+1)!", np.array([math.factorial(a + 1) for a in range(T + 1)], dtype=float), "超指数"),
    ("stairs_2^a_times_(a+1)", (2.0 ** np.arange(T + 1)) * (np.arange(T + 1) + 1.0), "几何×多项式"),
]
results = [classify(n, w, note) for n, w, note in profiles]

# 结论判定
native = [r for r in results if r["name"] == "native_branching_2^a"][0]
poly = [r for r in results if r["name"].startswith("polynomial")][0]
fact = [r for r in results if r["name"].startswith("factorial")][0]
uni = [r for r in results if r["name"] == "uniform"][0]
geomx = [r for r in results if r["name"].startswith("stairs")][0]

out = {
    "criterion": "G=<log(w_i/w_j)>：{0}->平凡；cZ->III_lambda；稠密->III_1",
    "profiles": results,
    "verdict": {
        "uniform_is_trivial": uni["type"].startswith("trivial"),
        "native_full_branching_type": native["type"],
        "native_lambda": native["lambda"],
        "polynomial_type": poly["type"],
        "factorial_type": fact["type"],
        "geometric_times_polynomial_type": geomx["type"],
        "bw_requires": "III_1",
        "consequence": (
            "若块权重按原生满分支指数增长，极限为 III_lambda(=1/2) 而非 III_1，"
            "几何模流不可得；几何模流反过来要求粗粒化 pi 的渐近轮廓非指数（幂律/超指数），"
            "即 P2' 变成对缺失输入 pi 的一个新约束"
        ),
    },
}

with io.open(os.path.join(HERE, "R35_type_iii_results.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=1)

for r in results:
    print(json.dumps(r, ensure_ascii=False))
print()
print(json.dumps(out["verdict"], ensure_ascii=False, indent=1))
