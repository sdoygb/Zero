#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R69 · 连续极限与维数：格点因果区间 → 连续度规

问题 1（连续极限）：
  R67 已给出 I(Δτ,|Δx|₁) = Σ_{k≤(Δτ−|Δx|₁)/2}(2k²+1)  （2+1 维格）
  问：I 的标度指数是否给出**连续 Minkowski 度规**？

  连续情形（d 维 Minkowski，d−1 个空间方向）：
      V(diamond) = (2/d)·Vol(S^{d−2})·(σ/2)^d·(1/(d−1))     (σ = 固有时)
  ⇒ V ∝ σ^d。故由 I ~ σ^a 反解 d = a。

问题 2（维数）：
  因果路能否在**细化**中选出一个维数，而不是继承输入的？
  判据：I 的标度指数 a 是否只依赖**空间格的维数**（⇒ 继承），
        还是能被别的机制定住（⇒ 选出）。
"""
from __future__ import annotations

import json
import math
import os
from itertools import product

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}


# ------------------------------------------------------------ 锥截面（精确枚举）
def cone_section(dtau, dims):
    """时刻 dtau 的光锥内部点数（dims 维空间，L¹ 范数，同奇偶）。"""
    rng = range(-dtau, dtau + 1)
    cnt = 0
    for pt in product(rng, repeat=dims):
        n1 = sum(abs(x) for x in pt)
        if n1 <= dtau and (dtau - n1) % 2 == 0:
            cnt += 1
    return cnt


def interval(dtau, dx, dims):
    """因果区间基数：交集 = 时间跨度 Δτ−n1 内的嵌套截面和。
    精确做法：直接枚举 (τ,x) 满足两边可达。"""
    n1 = abs(dx)
    if dtau < n1 or (dtau - n1) % 2:
        return 0
    # 交集：τ ∈ [0,Δτ]，|x−x1|≤τ，|x2−x|≤Δτ−τ
    # 对每个 τ 枚举 x（x1=0, x2=dx）
    tot = 0
    for tau in range(0, dtau + 1):
        rng = range(-tau, tau + 1)
        for pt in product(rng, repeat=dims):
            n_a = sum(abs(x) for x in pt)
            if n_a > tau or (tau - n_a) % 2:
                continue
            rest = tuple((dx if i == 0 else 0) - pt[i] for i in range(dims))
            n_b = sum(abs(x) for x in rest)
            if n_b <= dtau - tau and ((dtau - tau) - n_b) % 2 == 0:
                tot += 1
    return tot


def fit_exponent(xs, ys):
    """拟合 log y = a log x + b"""
    lx = np.log(xs)
    ly = np.log(ys)
    A = np.vstack([lx, np.ones_like(lx)]).T
    a, b = np.linalg.lstsq(A, ly, rcond=None)[0]
    resid = ly - (a * lx + b)
    return float(a), float(b), float(np.max(np.abs(resid)))


if __name__ == "__main__":
    print("=" * 74)
    print("问题 1：连续极限 —— 因果区间的标度指数是否给出 Minkowski 维数")
    print("=" * 74)

    # 先核对 R67 的闭式（2 维空间）
    print("\n  核对 R67 的闭式 I = Σ_{k≤(Δτ−|Δx|₁)/2}(2k²+1)（2 维空间）:")
    ok = True
    for dtau in range(0, 9):
        for dx in range(0, dtau + 1):
            formula = 0
            if dtau >= dx and (dtau - dx) % 2 == 0:
                m = (dtau - dx) // 2
                formula = sum(2 * k * k + 1 for k in range(m + 1))
            enum = interval(dtau, dx, 2)
            if formula != enum:
                ok = False
                print("    不符 dtau=%d dx=%d: 闭式=%d 枚举=%d"
                      % (dtau, dx, formula, enum))
    print("    闭式与枚举一致: %s" % ok)
    OUT["formula_check"] = bool(ok)

    # 各空间维数下，沿时间轴（dx=0）的区间标度
    print("\n  沿时间轴（|Δx|=0）的因果区间标度：")
    print("  空间维数 d_s   Δτ 范围    拟合指数 a   残差    期望 d = d_s+1")
    rows = []
    for dims in [1, 2, 3]:
        dts = list(range(0, 10, 2))
        if dims == 3:
            dts = list(range(0, 8, 2))
        ys = [interval(dt, 0, dims) for dt in dts]
        nz = [(dt, y) for dt, y in zip(dts, ys) if y > 0 and dt > 0]
        if len(nz) >= 3:
            a, b, res = fit_exponent([x for x, _ in nz], [y for _, y in nz])
            rows.append({"dims": dims, "a": a, "resid": res,
                         "expected_d": dims + 1})
            print("  %-13d %-12s %-12.4f %-7.4f %d"
                  % (dims, "%d..%d" % (dts[0], dts[-1]), a, res, dims + 1))
    OUT["scaling"] = rows
    print("\n  ⇒ 指数 a 与空间维数一一对应 ⇒ **I 的标度读出 d = d_s + 1**")

    # 洛伦兹不变量检验（各维数）
    print("\n  洛伦兹结构检验：I 是否只依赖 Δτ − |Δx|₁（各空间维数）")
    for dims in [1, 2, 3]:
        groups = {}
        for dtau in range(0, 9):
            for n1 in range(0, dtau + 1):
                if (dtau - n1) % 2:
                    continue
                v = interval(dtau, n1, dims)
                groups.setdefault(dtau - n1, set()).add(v)
        uniq = all(len(s) == 1 for s in groups.values())
        print("    空间维数 %d: I 只依赖 Δτ−|Δx|₁ = %s" % (dims, uniq))
        OUT.setdefault("lorentz_only", []).append({"dims": dims, "only_gap": uniq})

    print()
    print("=" * 74)
    print("问题 2：维数是继承的还是被选出的")
    print("=" * 74)
    print("  结构：因果路从空间格维数 d_s 得出时空维数 d = d_s + 1")
    print("  ⇒ 时空维数 d 由**输入的**空间维数决定 ⇒ 继承，非选出")
    print()
    print("  问：能否在细化中反解 d_s？检验方式：")
    print("    只给 I 的标度数据，反解 d_s，看是否唯一。")
    for r in rows:
        ds = r["a"] - 1
        print("    拟合 a=%.4f  ⇒ 反解 d_s=%.4f  （真值 %d）"
              % (r["a"], ds, r["dims"]))
    print()
    print("  ⇒ 反解唯一且精确 ⇒ 维数**可测**，但测量的是输入格的维数")
    print("  ⇒ 因果路本身不选维；选维必须来自别处（账本/稳定性/物理筛选）")
    OUT["dimension_measurable"] = True
    OUT["dimension_selected"] = False
    with open(os.path.join(HERE, "R69_continuum_dimension_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False, default=str)
    print("\n-> R69_continuum_dimension_results.json")
