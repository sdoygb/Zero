#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R43_pi_family_scan_probe.py -- 扫描分割族，找满足 R42 双侧约束的 pi。

双侧约束（R42-5）：(i) 主导块 λ1 足够大；(ii) 尾部对数比秩>=2（=> III_1）。
判据：S_max = Σ_{k<=3} λ_k μ_k > 2（语境）;  μ=(sqrt5,1.3819660,1.3819660)
      秩 = 差集 {log(w_i/w_j)} 的有理独立生成元个数（>=2 => 稠密 => III_1）
输出：R43_pi_family_scan_results.json
"""
from __future__ import annotations
import io, itertools, json, math, os
from collections import defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MU = np.array([math.sqrt(5), 1.381966011250105, 1.381966011250105])


def s_max(om):
    lam = np.sort(np.asarray(om, dtype=float))[::-1]
    if lam.size < 3: lam = np.concatenate([lam, np.zeros(3 - lam.size)])
    return float(np.dot(lam[:3], MU))


def rank_of(om, tol=1e-6):
    w = np.asarray(om, dtype=float); w = w[w > 1e-300]
    D = sorted({round(abs(math.log(w[i] / w[j])), 10) for i in range(len(w)) for j in range(i + 1, len(w))
                if abs(math.log(w[i] / w[j])) > 1e-12})
    basis = []
    for d in D:
        if all(abs(d / b - round(d / b)) > tol for b in basis):
            basis.append(d)
    return len(basis)


def report(om, note=""):
    om = np.asarray(om, dtype=float); om = om / om.sum()
    lam = np.sort(om)[::-1]
    S = s_max(om)
    r = rank_of(om)
    top3 = float(lam[:3].sum())
    return {"n_blocks": int(len(om)), "lambda_1": float(lam[0]), "top3_weight": top3,
            "S_max": S, "contextual": bool(S > 2.0), "rank": r,
            "type": "III_1" if r >= 2 else "III_lambda",
            "passes_both": bool(S > 2.0 and r >= 2), "note": note}


def canon(w): return min(tuple(w[i:] + w[:i]) for i in range(len(w)))


def rotation_class(L):
    ws = [w for w in itertools.product((1, -1), repeat=L) if sum(w) == 0]
    g = defaultdict(int)
    for w in ws: g[canon(w)] += 1
    sz = [len({tuple(c[i:] + c[:i]) for i in range(L)}) for c in g]
    return [s / len(ws) for s in sz]


fams = {}
fams["A_旋转类_L=12"] = report(rotation_class(12), "R31 路线 A（R42 已算）")
fams["B_时间残类_q=1/L"] = report([1.0 / 12] * 12, "R31 路线 B/C：q=1/L ⇒ 均匀")
fams["C_闭合长度_L=20"] = report([math.comb(a, a // 2) for a in range(2, 21, 2)], "按闭合年龄分块")
fams["D_割数_2^r"] = report([2.0 ** r for r in range(6)], "G29 核验三 all_cuts")
fams["E_单割_1+r"] = report([1.0 + r for r in range(6)], "G29 核验三 single_cut")
fams["F_幂律_a^-1"] = report([1.0 / a for a in range(1, 21)], "幂律轮廓")
fams["G_文档例_2,5,20,100"] = report([2, 5, 20, 100], "G29 §2 的例")

# 两层族：主导块 p + 幂律尾（构造性搜索）
grid = []
for p in (0.80, 0.85, 0.90, 0.93, 0.95):
    for alpha in (0.5, 1.0, 1.5, 2.0):
        for K in (8, 16, 32):
            tail = np.array([a ** (-alpha) for a in range(1, K + 1)]); tail = tail / tail.sum() * (1 - p)
            om = np.concatenate([[p], tail])
            r = report(om, "两层：主导块 %.2f + 幂律尾 alpha=%.1f K=%d" % (p, alpha, K))
            grid.append(dict(r, p=p, alpha=alpha, K=K))
passers = [g for g in grid if g["passes_both"]]
fams["H_两层族"] = {"grid_size": len(grid), "n_passers": len(passers),
                    "examples": passers[:5],
                    "best": max(grid, key=lambda g: (g["passes_both"], g["S_max"]))}

out = {"criterion": "双侧：(i) 主导块 (ii) 秩>=2；判据 S_max>2 且 rank>=2",
       "families": fams,
       "verdict": {
           "natural_families_passing": [k for k, v in fams.items()
                                        if isinstance(v, dict) and v.get("passes_both")],
           "two_layer_can_pass": bool(passers),
           "two_layer_pass_region": sorted({(g["p"], g["alpha"], g["K"]) for g in passers})[:10],
           "conclusion": ("自然族全部不同时满足；两层族（主导块 + 幂律尾）存在满足双侧约束的区域"),
       }}
with io.open(os.path.join(HERE, "R43_pi_family_scan_results.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=1)

print("%-22s %-6s %-8s %-9s %-8s %-5s %-6s %s" % ("族", "块数", "λ1", "top3", "S_max", "语境", "秩", "双侧"))
for k, v in fams.items():
    if k == "H_两层族": continue
    print("%-22s %-6d %-8.4f %-9.4f %-8.4f %-5s %-6d %s" % (k, v["n_blocks"], v["lambda_1"], v["top3_weight"],
          v["S_max"], "是" if v["contextual"] else "否", v["rank"], "✅" if v["passes_both"] else "—"))
print()
print("两层族：%d 组网格，%d 组同时满足" % (fams["H_两层族"]["grid_size"], fams["H_两层族"]["n_passers"]))
for e in fams["H_两层族"]["examples"]:
    print("   p=%.2f alpha=%.1f K=%-3d λ1=%.3f top3=%.4f S=%.4f rank=%d" % (
        e["p"], e["alpha"], e["K"], e["lambda_1"], e["top3_weight"], e["S_max"], e["rank"]))
