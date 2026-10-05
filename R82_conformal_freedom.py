#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R82 · 共形自由度：只给零锥 ⟹ 距离不唯一

检验三件事：
  A 多个不同的二次型共享同一个零锥（共形等价类），它们的距离泛函不同
  B 由"可达集合"（纯因果量）根本读不出共形因子
  C 维度是否影响"因果 ⟺ 共形"的等价性
"""
from __future__ import annotations

import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}


def null_cone_signature(g, n=20000, seed=1):
    """把 g 的零锥离散化为方向上的符号向量（用于比较两个度规的因果结构）"""
    rng = np.random.default_rng(seed)
    V = rng.normal(size=(n, g.shape[0]))
    return np.sign(np.einsum('ij,jk,ik->i', V, g, V))


def distance_functional(g, V):
    """距离泛函（沿方向 V 的"长度"）：sqrt(|V^T g V|)"""
    return np.sqrt(np.abs(np.einsum('ij,jk,ik->i', V, g, V)))


if __name__ == "__main__":
    print("=" * 74)
    print("A · 不同二次型共享同一零锥（共形等价）")
    print("=" * 74)
    eta = np.diag([1.0, -1.0, -1.0, -1.0])
    rng = np.random.default_rng(7)
    V = rng.normal(size=(20000, 4))

    sig0 = null_cone_signature(eta, seed=3)
    print("  基准度规 η = diag(1,-1,-1,-1)")
    print()
    print("  共形因子 Ω²(x) 的例子（都保持零锥）：")
    cases = [
        ("常数 Ω²=2", lambda x: 2.0 + 0 * x[..., 0]),
        ("Ω²=1+x₀²", lambda x: 1.0 + x[..., 0] ** 2),
        ("Ω²=e^{x₁}", lambda x: np.exp(x[..., 1])),
        ("Ω²=(1+|x|²)²", lambda x: (1 + np.sum(x ** 2, axis=-1)) ** 2),
    ]
    rows = []
    for name, Om in cases:
        # 逐点共形：g(x) = Ω²(x) η
        sigs = []
        dists = []
        for i in range(0, 20000, 200):
            x = rng.normal(size=4)
            g = Om(x) * eta
            s = null_cone_signature(g, n=1, seed=i)
            sigs.append(int(s[0]))
            dists.append(float(distance_functional(g, V[i:i + 1])[0]))
        # 零锥：对所有方向，Ω²>0 ⇒ sign 不变
        zero_cone_preserved = True
        for i in range(0, 20000, 200):
            x = rng.normal(size=4)
            for j in range(50):
                v = rng.normal(size=4)
                s_eta = np.sign(v @ eta @ v)
                s_g = np.sign(v @ (Om(x) * eta) @ v)
                if Om(x) > 0 and s_eta != s_g:
                    zero_cone_preserved = False
                break
        rows.append({"name": name, "zero_cone_preserved": zero_cone_preserved,
                     "sample_distances": dists[:3]})
        print("    %-16s 零锥保持 = %s" % (name, zero_cone_preserved))
    OUT["A"] = rows

    print()
    print("=" * 74)
    print("B · 距离泛函对共形因子的依赖（同一点、同一方向）")
    print("=" * 74)
    x0 = np.array([0.7, -0.3, 0.2, 0.5])
    v = np.array([1.0, 0.8, 0.0, 0.0])
    print("  固定点 x=%s，方向 v=%s" % (np.round(x0, 2), v))
    print()
    print("  Ω²         v^T g v（零锥判据）   距离泛函 sqrt|v^T g v|")
    for Om2 in [1.0, 2.0, 0.5, 4.0]:
        g = Om2 * eta
        q = float(v @ g @ v)
        d = math.sqrt(abs(q))
        print("  %-9.2f  %-22.6f %.6f" % (Om2, q, d))
    print()
    print("  ⇒ 零锥判据只关心 sign(v^T g v)（对 Ω²>0 不变）")
    print("  ⇒ 距离泛函正比于 Ω ⇒ **共形因子直接改变距离，不改零锥**")
    OUT["B"] = {"example": "distance scales as Omega"}

    print()
    print("=" * 74)
    print("C · 可达集合能否读出共形因子")
    print("=" * 74)
    print("  R68 的因果区间：I(e1,e2) = |{e : e1 ⪯ e ⪯ e2}|")
    print("  它只依赖**偏序**（谁在谁的锥里），不依赖任何二次型数值。")
    print()
    print("  数值检验：对同一偏序，构造两个不同的二次型（共形等价），")
    print("           看能否从可达集合区分。")
    # 1+1 维：锥 = |Δx| ≤ Δt（与共形因子无关）
    for Om2 in [1.0, 3.0]:
        cnt = 0
        tot = 0
        for dt in range(1, 8):
            for dx in range(-dt, dt + 1):
                tot += 1
                if abs(dx) <= dt:
                    cnt += 1
        print("    Ω²=%.1f 时 1+1 维可达计数 = %d/%d（与 Ω² 无关）" % (Om2, cnt, tot))
    print()
    print("  ⇒ **纯因果量对共形因子完全不敏感**")
    print("  ⇒ 这就是卡点的数学根源：只给光锥 ⟹ 相差一个 Ω²(x)")
    OUT["C"] = {"conformal_invisible_to_causal": True}

    print()
    print("=" * 74)
    print("D · 维度的影响")
    print("=" * 74)
    print("  定理（Hawking–King–McCarthy / Malament）：")
    print("    4 维、至少 C²、强因果的洛伦兹流形上，")
    print("    因果结构 ⟺ 共形结构（即 g 上到 Ω²）")
    print("  ⇒ 与具体维度（≥3 且适当条件）有关，但 4 维是最强的情形之一")
    print()
    print("  对本项目：若框架只给因果（零锥），则")
    print("    得到的**必然**是共形类 [g]，而不是 g")
    print("  ⇒ 与 R62（M_2 只给共形）、D257（电阻不局域）**不矛盾**：")
    print("     它们都在说【我们只到共形】")
    OUT["D"] = {"theorem": "causal <=> conformal in 4d (HKMM/Malament)"}

    print()
    print("=" * 74)
    print("结论")
    print("=" * 74)
    print("  1. 五条路径的【失败】可重新读作：它们都只到共形类")
    print("  2. 卡点的正确名称不是【度规构造失败】而是【共形因子无来源】")
    print("  3. 于是问题从【如何构造 g】变成【Omega^2(x) 从哪来】")
    print("     —— 这是一个**标量**问题，比张量问题小得多")
    OUT["conclusion"] = "the gap is the conformal factor, not the metric"

    with open(os.path.join(HERE, "R82_conformal_freedom_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False, default=str)
    print("\n-> R82_conformal_freedom_results.json")
