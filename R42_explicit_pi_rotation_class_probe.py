#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R42_explicit_pi_rotation_class_probe.py -- 显式 pi（旋转类）的块结构，及两项二值判据。

pi = 旋转类（R31 路线 A）：块 = 平衡词的 Z_L 旋转轨道，块大小 o_c，N_L = C(L,L/2)。
  权重 omega_c = o_c / N_L  =>  原生态的谱 lambda = {omega_c}
校验：q_L = sum_c omega_c^2 应与 R31 的表一致（5/9, 7/25, 19/175, ...）
判据：
  R37（单体语境性）：S_max = sum_{k<=3} lambda_k mu_k > 2 ?   mu = (sqrt5, 1.3819660, 1.3819660)
  R35（模论类型）：差集 <log(w_i/w_j)> 循环 -> III_lambda；稠密 -> III_1

输出：R42_explicit_pi_rotation_class_results.json
"""
from __future__ import annotations
import io, itertools, json, math, os
from collections import Counter, defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MU = np.array([math.sqrt(5), 1.381966011250105, 1.381966011250105])


def balanced_words(L):
    return [w for w in itertools.product((1, -1), repeat=L) if sum(w) == 0]


def canon(w):
    return min(tuple(w[i:] + w[:i]) for i in range(len(w)))


def orbits(L):
    words = balanced_words(L)
    groups = defaultdict(list)
    for w in words:
        groups[canon(w)].append(w)
    # 轨道大小 = 不同旋转的个数
    sizes = []
    for c, ws in groups.items():
        rots = {tuple(c[i:] + c[:i]) for i in range(L)}
        sizes.append(len(rots))
    # 校验：sum(sizes) == N_L
    return sorted(sizes, reverse=True), len(words)


def s_max_native(spectrum):
    lam = np.sort(np.asarray(spectrum, dtype=float))[::-1]
    if lam.size < 3:                      # 不足 3 个本征值时补零（3 维子空间可用空能级）
        lam = np.concatenate([lam, np.zeros(3 - lam.size)])
    return float(np.dot(lam[:3], MU))


def r35_classify(spectrum, window=30.0):
    w = np.sort(np.asarray(spectrum, dtype=float))[::-1]
    logw = np.log(w)
    if np.ptp(logw) < 1e-12:
        return "trivial", None, 1.0
    vals = sorted({round(abs(logw[i] - logw[j]), 12)
                   for i in range(len(logw)) for j in range(i + 1, len(logw))
                   if 1e-12 < abs(logw[i] - logw[j]) <= window})
    if not vals:
        return "undecided", None, None
    D = np.array(vals); c = float(D.min())
    rho = float(np.abs(D / c - np.round(D / c)).max())
    kind = "III_lambda" if rho < 1e-3 else "III_1"
    lam = math.exp(-c) if rho < 1e-3 else None
    return kind, lam, rho


rows = []
r31 = {2: 1.0, 4: 5 / 9, 6: 7 / 25, 8: 19 / 175}
for L in range(2, 19, 2):
    sizes, N = orbits(L)
    assert sum(sizes) == N, (L, sum(sizes), N)
    omega = [s / N for s in sizes]
    q = sum(o * o for o in omega)
    Smax = s_max_native(omega)
    kind, lam, rho = r35_classify(omega)
    rows.append({
        "L": L, "N_L": N, "n_classes": len(sizes),
        "orbit_sizes": sizes[:8], "q_L": q,
        "q_L_matches_R31": (abs(q - r31[L]) < 1e-9) if L in r31 else None,
        "lambda_1": max(omega),
        "S_max_top3": Smax, "contextual": bool(Smax > 2.0),
        "r35_type": kind, "r35_lambda": lam, "r35_cyclic_residual": rho,
    })

# 对照：年龄分割（满分支 2^a）
T = 12
age = 2.0 ** np.arange(T + 1); age = age / age.sum()
age_S = s_max_native(age)
age_kind, age_lam, age_rho = r35_classify(age)

out = {
    "pi_candidate": "旋转类（R31 路线 A）：块 = 平衡词的 Z_L 旋转轨道",
    "cross_check_R31": {str(k): v for k, v in r31.items()},
    "rows": rows,
    "alternative_age_partition": {
        "profile": "2^a (满分支)", "lambda_1": float(max(age)),
        "S_max_top3": age_S, "contextual": bool(age_S > 2.0),
        "r35_type": age_kind, "r35_lambda": age_lam,
    },
    "verdict": {
        "rotation_class_all_contextual": all(r["contextual"] for r in rows),
        "rotation_class_types": sorted({r["r35_type"] for r in rows}),
        "rotation_class_lambda1_range": [min(r["lambda_1"] for r in rows), max(r["lambda_1"] for r in rows)],
        "age_partition_contextual": bool(age_S > 2.0),
        "age_partition_type": age_kind,
    },
}
with io.open(os.path.join(HERE, "R42_explicit_pi_rotation_class_results.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=1)


# ---------------------------------------------------------------- 渐近类型：有理独立对数比
def rational_rank(spectrum, tol=1e-9):
    """差集 {log(w_i/w_j)} 中"有理独立"的生成元个数（>=2 且比值无理 => 稠密 => III_1）"""
    w = np.sort(np.asarray(spectrum, dtype=float))[::-1]
    D = sorted({round(abs(math.log(w[i]/w[j])), 12) for i in range(len(w)) for j in range(i+1, len(w))
                if abs(math.log(w[i]/w[j])) > 1e-12})
    if not D:
        return 0, []
    basis = []
    for d in D:
        v = [d]
        for b in basis:
            v.append(round(d / b))
        frac = [abs(d / b - round(d / b)) for b in basis]
        if all(f > 1e-6 for f in frac):
            basis.append(d)
    return len(basis), [round(b, 6) for b in basis]

rank_rows = []
for L in range(4, 19, 2):
    sz, N = orbits(L)
    om = [s / N for s in sz]
    r, basis = rational_rank(om)
    rank_rows.append({"L": L, "rank": r, "basis_logs": basis, "type": "III_1" if r >= 2 else "III_lambda"})
out["asymptotic_type"] = {
    "method": "差集生成元的秩：秩>=2 且含无理比 => 稠密 => III_1",
    "rows": rank_rows,
    "final_type": "III_1" if rank_rows[-1]["rank"] >= 2 else "III_lambda",
    "note": "轨道大小形如 L/d (d|L)，故比值为整数比 log p；L 含 >=2 个素因子时秩>=2 => III_1",
}
with io.open(os.path.join(HERE, "R42_explicit_pi_rotation_class_results.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=1)
print("rank rows:", [(r["L"], r["rank"], r["type"]) for r in rank_rows])

print("L   N_L 类数  q_L        q(L)=R31?  λ1      S_max   语境?  R35 类型")
for r in rows:
    print("%-3d %-4d %-4d %.6f  %-8s %.4f  %.4f  %-5s %s" % (
        r["L"], r["N_L"], r["n_classes"], r["q_L"],
        str(r["q_L_matches_R31"]), r["lambda_1"], r["S_max_top3"],
        "是" if r["contextual"] else "否", r["r35_type"]))
print()
print(json.dumps(out["alternative_age_partition"], ensure_ascii=False))
print()
print(json.dumps(out["verdict"], ensure_ascii=False, indent=1))
