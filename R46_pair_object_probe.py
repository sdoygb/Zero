#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R46_pair_object_probe.py -- 对账本的"对象"是什么？导出或排除 C(D+1,2)。

思路：C(k,2) 的 k = 对象的个数。
  · 若对象 = D 个独立方向（洛伦兹平面）  => C(D,2) = dim so(D-1,1)   [R31 标准]
  · 若对象 = D+1 个通道/顶点（D-单纯形的顶点） => C(D+1,2) = 单纯形边数  [欧氏字典]
检验：仓库自身模型里的宏类计数是否恰为单纯形顶点数。
  G29 核验三：single_cut: M = 1+r ; all_cuts: M = 2^r ; sterile: M = 1
  单纯形：r-单纯形顶点数 = r+1  => M = 1+r 恰为 r-单纯形顶点数 ✓
                边数 = C(r+1,2)
输出：R46_pair_object_results.json
"""
from __future__ import annotations
import io, json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))


def simplex_vertices(r):  # r-单纯形
    return r + 1


def simplex_edges(r):
    return math.comb(r + 1, 2)


rows = []
for r in range(1, 8):
    V = simplex_vertices(r)
    E = simplex_edges(r)
    rows.append({
        "r": r,
        "M_single_cut": 1 + r,          # G29 核验三
        "M_all_cuts": 2 ** r,
        "simplex_vertices_r": V,
        "single_cut_equals_simplex_vertices": (1 + r) == V,
        "all_cuts_equals_simplex_vertices": (2 ** r) == V,
        "simplex_edges": E,
        "C(r,2)": math.comb(r, 2),
        "C(r+1,2)": math.comb(r + 1, 2),
    })

# 两种字典的对账
face_off = {
    "lorentzian": {
        "objects": "D 个独立方向（旋转平面）",
        "pairs": "C(D,2) = dim so(D-1,1)",
        "signature": "洛伦兹内建 ✓",
        "window": [0.5, 0.6],
        "contextuality_compatible": False,
    },
    "euclidean_simplex": {
        "objects": "D+1 个通道（D-单纯形顶点）",
        "pairs": "C(D+1,2) = 单纯形边数 = dim so(D+1)",
        "signature": "欧氏/无签名 ✗（须外供）",
        "window": [0.6, 2 / 3],
        "contextuality_compatible": True,
    },
}

out = {
    "question": "对账本的『对象』是 D 个方向，还是 D+1 个通道？",
    "simplex_table": rows,
    "face_off": face_off,
    "verdict": {
        "single_cut_matches_simplex_vertices": all(r["single_cut_equals_simplex_vertices"] for r in rows),
        "all_cuts_matches_simplex_vertices": all(r["all_cuts_equals_simplex_vertices"] for r in rows),
        "derivation_candidate": ("仓库自身模型 single_cut: M=1+r 恰为 r-单纯形顶点数 ⇒ 若对账本数的是通道对（单纯形边），"
                                 "则多重度 = C(r+1,2) = 欧氏字典"),
        "trade_off": "洛伦兹字典签名内建但与语境性互斥；单纯形字典与语境性相容但签名须外供（G59/G33 的因果结构）",
        "next": "判据：因果结构能否供出签名 —— 能则 C(D+1,2) 合法（D=4 与量子性同时到手）；不能则 R44 的互斥重新生效",
    },
}
with io.open(os.path.join(HERE, "R46_pair_object_results.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=1)

print("r   M_single_cut=1+r  单纯形顶点r+1  相等?  M_all_cuts=2^r  边数C(r+1,2)  C(r,2)")
for r in rows:
    print("%-3d %-16d %-13d %-5s %-15d %-13d %d" % (
        r["r"], r["M_single_cut"], r["simplex_vertices_r"],
        "✓" if r["single_cut_equals_simplex_vertices"] else "✗",
        r["M_all_cuts"], r["simplex_edges"], r["C(r,2)"]))
print()
print(json.dumps(out["verdict"], ensure_ascii=False, indent=1))
