#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R45_ledger_form_scan_probe.py -- 扫描账本形式 F_D = M(D) * q^{E(D)}，找"D=4 峰窗口"落在 q>3/5 一侧者。

(R44) 标准形式 M=C(D,2), E=D 的窗口 = (1/2, 3/5)，与语境性区 (3/5,1] 恰好互补 -> no-go。
本探针：扫描 M 与 E 的一般族，看是否存在窗口与 (3/5,1] 相交者。
判据：argmax_D F_D(q) == 4 且 q > 3/5  =>  选维与语境性可共存。
输出：R45_ledger_form_scan_results.json
"""
from __future__ import annotations
import io, json, math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MU2 = (5.0 - math.sqrt(5.0)) / 2.0
MU = (math.sqrt(5.0), MU2, MU2)
QS = np.round(np.arange(0.30, 1.0001, 0.002), 6)
Q_CTX = 0.6


def M_forms():
    return {
        "C(D,2)": lambda D: math.comb(D, 2),
        "C(D+1,2)": lambda D: math.comb(D + 1, 2),
        "C(D,3)": lambda D: math.comb(D, 3) if D >= 3 else 0,
        "D": lambda D: D,
        "D^2": lambda D: D * D,
        "D^3": lambda D: D ** 3,
        "2^D": lambda D: 2 ** D,
        "(D-1)!": lambda D: math.factorial(D - 1),
    }


def E_forms():
    return {
        "D": lambda D: D,
        "D+1": lambda D: D + 1,
        "D-1": lambda D: D - 1,
        "C(D,2)": lambda D: math.comb(D, 2),
        "C(D+1,2)": lambda D: math.comb(D + 1, 2),
        "D^2": lambda D: D * D,
        "2^D": lambda D: 2 ** D,
        "1": lambda D: 1,
    }


def peak_D(M, E, q, Dmax=30):
    vals = [(D, M(D) * q ** E(D)) for D in range(2, Dmax + 1)]
    return max(vals, key=lambda t: t[1])[0]


def window_of(M, E, Dmax=30):
    """使 argmax = 4 的 q 区间"""
    good = [q for q in QS if peak_D(M, E, q, Dmax) == 4]
    if not good:
        return None
    return (float(min(good)), float(max(good)))


rows = []
for mn, M in M_forms().items():
    for en, E in E_forms().items():
        w = window_of(M, E)
        if w is None:
            rows.append({"M": mn, "E": en, "window": None, "overlap_ctx": False})
            continue
        lo, hi = w
        # 与语境性区 (0.6, 1] 的交
        over_lo, over_hi = max(lo, Q_CTX), hi
        overlap = over_hi > over_lo + 1e-9
        rows.append({"M": mn, "E": en, "window": [lo, hi], "overlap_ctx": bool(overlap),
                     "overlap_interval": [over_lo, over_hi] if overlap else None})

# 具体载体
carriers = {
    "documented_2_5_20_100": [0.0157, 0.0394, 0.1575, 0.7874],
    "routeA_L4": [4 / 6, 2 / 6],
    "two_layer_p0.8_a2_K8": None,
}
tail = np.array([a ** (-2.0) for a in range(1, 9)]); tail = tail / tail.sum() * 0.2
carriers["two_layer_p0.8_a2_K8"] = np.concatenate([[0.8], tail]).tolist()

def s_spec(lam):
    lam = np.sort(np.asarray(lam, float))[::-1]
    if lam.size < 3: lam = np.concatenate([lam, np.zeros(3 - lam.size)])
    return float(np.dot(lam[:3], MU))

car = {}
for k, v in carriers.items():
    q = float(sum(w * w for w in v))
    car[k] = {"q": q, "S_max": s_spec(v), "contextual": bool(s_spec(v) > 2),
              "peak_D_standard": peak_D(M_forms()["C(D,2)"], E_forms()["D"], q),
              "peak_D_Cplus": peak_D(M_forms()["C(D+1,2)"], E_forms()["D"], q)}

out = {
    "family": "F_D = M(D) * q^{E(D)}",
    "contextuality_region": "q > 3/5",
    "standard_form": "M=C(D,2), E=D  ->  window (1/2, 3/5)  (R44 no-go)",
    "scan": rows,
    "carriers": car,
    "verdict": {
        "n_forms": len(rows),
        "forms_with_overlap": [{"M": r["M"], "E": r["E"], "window": r["window"],
                                "overlap": r["overlap_interval"]}
                               for r in rows if r["overlap_ctx"]],
        "escape_exists": any(r["overlap_ctx"] for r in rows),
    },
}
with io.open(os.path.join(HERE, "R45_ledger_form_scan_results.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=1)

print("%-12s %-10s %-22s %s" % ("M(D)", "E(D)", "峰在D=4的 q 窗口", "与 q>0.6 相交"))
for r in rows:
    w = r["window"]
    print("%-12s %-10s %-22s %s" % (r["M"], r["E"],
          ("[%.4f, %.4f]" % (w[0], w[1])) if w else "无",
          ("✅ [%.4f, %.4f]" % tuple(r["overlap_interval"])) if r["overlap_ctx"] else "—"))
print()
for k, v in car.items():
    print("%-24s q=%.4f S=%.4f 语境=%s  标准字典峰D=%d  C(D+1,2)字典峰D=%d" % (
        k, v["q"], v["S_max"], v["contextual"], v["peak_D_standard"], v["peak_D_Cplus"]))
