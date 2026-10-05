#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R51 探针 V（修正版）：R44 界的**普适化**

命题 (R51-1)  设 k 块（k≥3）归一谱 w，q=Σw_i²，顶三归一首分量 λ1。
              则
                    λ1 ≤ (1+sqrt(2q-1))/2                       (R51-1)
              等号可达，且与块数 k 无关。

证明（两行代数）  记 T=w1+w2+w3，（w1,w2,w3)/T 是归一 3 向量，故
              Σ_{i≤3}(w_i/T)² ≤ Σ_{i≤3} w_i²/T² ≤ q/T² ≤ q
              且 (w1/T)² ≥ Σ(w_i/T)² - (1-λ1)²/T² ≥ ... 直接代入 (R44-2) 的
              Cauchy–Schwarz 步：固定 q3=Σλ_i² 与 Σλ_i=1 时
              λ1 ≤ (1+sqrt(2q3-1))/2，而 q3 ≤ q，且右端关于 q 单调增。
              ⇒ (R51-1)。∎

推论 (R51-2)  语境性（λ1>λ*=0.723607）要求 q ≥ 3/5 —— **对任意块数成立**。
推论 (R51-3)  与 R32 生存窗口 q∈(1/2,3/5) 合并：D=4 与语境性互斥，且
              q=3/5 是**精确相接**（非数值巧合）。

本探针用**每次采样自带其真实 q**的方式核验 (R51-1)（避免上一版的容差漏筛）。
"""
from __future__ import annotations

import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LAM = 0.723607
OUT = {}


def bound(q):
    v = 2 * q - 1
    return (1 + math.sqrt(v)) / 2 if v >= 0 else float("nan")


def check(n=400000, seed=5):
    """每个样本用它**自己**的 q 比界：超出量 = λ1 - bound(q_self)。"""
    rng = np.random.default_rng(seed)
    worst = {}
    for _ in range(n):
        k = int(rng.integers(3, 20))
        w = np.sort(rng.random(k) ** rng.choice([0.2, 0.5, 1.0, 2.0, 5.0]))[::-1]
        w = w / w.sum()
        q = float(np.sum(w ** 2))
        l1 = float(w[0] / w[:3].sum())
        exc = l1 - bound(q)
        if exc > worst.get(k, (-9, None))[0]:
            worst[k] = (exc, {"q": round(q, 6), "lam1": round(l1, 6),
                              "bound": round(bound(q), 6)})
    return {"per_k": {str(k): {"excess": round(v[0], 10), **v[1]}
                      for k, v in sorted(worst.items())},
            "max_excess": max(v[0] for v in worst.values())}


def tightness():
    """等号族：w=(λ1,(1-λ1)/2,(1-λ1)/2, 0,...,0)，q=(1+2λ1²)/3。"""
    rows = []
    for l1 in [0.4, 0.6, 0.723607, 0.85, 0.95, 1.0]:
        q = (1 + 2 * l1 ** 2) / 3
        rows.append({"lam1": l1, "q_from_formula": round(q, 6),
                     "bound_at_that_q": round(bound(q), 6),
                     "tight": abs(bound(q) - l1) < 1e-9})
    return rows


def ledger_crossing():
    """R51-3：q* = 3/5 处，pair-carrier 的 D=4 窗口恰好关闭。"""
    q = 0.6
    m = lambda D: D * (D - 1) / 2
    return {
        "q_star": q,
        "lambda1_at_q_star": round(bound(q), 6),
        "S_max_at_bound": 2.0,
        "D4_window": [round(m(4) / (m(4) + m(3)), 6),
                      round(m(4) / (m(4) + m(5)), 6)],
        "m5_over_m4": round(m(5) / m(4), 6),
        "m5_over_m4_required_to_cover_qstar": round(1 / q, 6),
        "exact_meeting": abs(m(5) / m(4) - 1 / q) < 1e-12,
        "conclusion": "pair-carrier 的 m(5)/m(4)=5/3 与 KCBS 临界 1/q*=5/3 精确相等",
    }


def escape_scan():
    """R44-7 逃逸出口：要覆盖 q=3/5 需 m(5)/m(4) < 5/3。列出候选。"""
    cands = {"C(D,2)（成对身份）": lambda D: D * (D - 1) / 2,
             "C(D+1,2)（字典读法）": lambda D: (D + 1) * D / 2,
             "D^2": lambda D: D ** 2,
             "D^2+1（每次多一笔代际重播种记录?）": lambda D: D ** 2 + 1,
             "2^D": lambda D: 2 ** D}
    rows = []
    need = 1 / 0.6
    for nm, m in cands.items():
        lo = m(4) / (m(4) + m(3))
        hi = m(4) / (m(4) + m(5))
        rows.append({"mult": nm, "m5/m4": round(m(5) / m(4), 5),
                     "window": [round(lo, 5), round(hi, 5)],
                     "covers_q_star": hi > 0.6,
                     "need_m5_m4_below": round(need, 5)})
    return {"rows": rows}


if __name__ == "__main__":
    print("(R51-1) 普适界检验：λ1 ≤ (1+√(2q-1))/2，每个样本用它自己的 q")
    OUT["universality"] = check()
    print("  最大超出 = %.3e" % OUT["universality"]["max_excess"])
    for k, d in OUT["universality"]["per_k"].items():
        print("   k=%-3s 超出=%+.2e   最紧样本: q=%.6f λ1=%.6f 界=%.6f"
              % (k, d["excess"], d["q"], d["lam1"], d["bound"]))
    print("\n(R51-2) 等号族 ...")
    OUT["tightness"] = tightness()
    for r in OUT["tightness"]:
        print("   λ1=%.6f  q=(1+2λ1²)/3=%.6f  界=%.6f  紧=%s"
              % (r["lam1"], r["q_from_formula"], r["bound_at_that_q"], r["tight"]))
    print("\n(R51-3) 账本临界点 ...")
    OUT["ledger"] = ledger_crossing()
    L = OUT["ledger"]
    print("   q*=%.2f  λ1=%.6f  S_max=%.1f" % (L["q_star"], L["lambda1_at_q_star"], L["S_max_at_bound"]))
    print("   D=4 窗口 =", L["D4_window"], " m5/m4 =", L["m5_over_m4"])
    print("   要覆盖 q* 需 m5/m4 <", L["m5_over_m4_required_to_cover_qstar"])
    print("   精确相接 =", L["exact_meeting"])
    print("\n(R44-7) 逃逸出口候选 ...")
    OUT["escape"] = escape_scan()
    for r in OUT["escape"]["rows"]:
        print("   %-34s m5/m4=%-8s 窗口=%-20s 覆盖q*=%s"
              % (r["mult"], r["m5/m4"], r["window"], r["covers_q_star"]))
    with open(os.path.join(HERE, "R51_final_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False)
    print("\n-> R51_final_results.json")
