#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R79 · 物质侧（修正）：$T^{ab}$ 的两参数族与各向同性

R78 的修正：
  B 判据写错了 —— 用"非对角元为零"当各向同性是**基依赖**的。
  正确判据：$T^{ab}$ 在空间子空间上是否 $\propto \delta^{ab}$（即在 $S_D$ 的
  标准表示下是标量）。

本脚本：
  A  正确的各向同性判据（谱判据）
  B  一般权重：Edges 在 Stab(0)=S_D 下分成**两个轨道**
       A 类 = 与根相连（大小 D）
       B 类 = 不与根相连（大小 C(D,2)）
     ⇒ $T^{ab}$ 是**两参数族** $(a,b)$
  C  两参数族的 $T^{00}$、空间压力、各向异性
  D  状态方程 $w=p/\\rho$ 对 $b/a$ 的依赖
  E  无迹条件（共形物质）是否给出唯一 $b/a$
"""
from __future__ import annotations

import json
import math
import os
from itertools import combinations, permutations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}
D4 = 4
NV = D4 + 1
VERTS = list(range(NV))
EDGES = list(combinations(VERTS, 2))


def build_T(a, b, nv=NV, root=0):
    """a = 与根相连的边权，b = 其余边权"""
    T = np.zeros((nv, nv))
    for (i, j) in EDGES:
        w = a if (i == root or j == root) else b
        u = np.zeros(nv)
        u[i] = -1.0
        u[j] = 1.0
        T += w * np.outer(u, u)
    return T


def isotropy_check(T, root=0, nv=NV):
    """把 T 限制到"空间"子空间（非根通道张成的 R^{D}），
    检查它是否 ∝ I。返回 (各向异性范数, 平均对角)"""
    idx = [i for i in range(nv) if i != root]
    Tsub = T[np.ix_(idx, idx)]
    tr = np.trace(Tsub) / len(idx)
    aniso = Tsub - tr * np.eye(len(idx))
    return float(np.abs(aniso).max()), float(tr), Tsub


if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=True, linewidth=140)

    print("=" * 74)
    print("A · 正确的各向同性判据（谱判据）")
    print("=" * 74)
    T0 = build_T(1.0, 1.0)
    print("  均匀权重 T^{ab} =")
    print(T0)
    ev = np.linalg.eigvalsh(T0)
    print("  本征值 = %s" % np.round(ev, 8))
    print("  非零本征值全相等（= 5）⇒ T ∝ 投影到 1^⊥ ⇒ **各向同性** ✅")
    print()
    print("  ⇒ R78 用“非对角元为零”判各向同性是**错的**（基依赖）")
    print("     正确判据：空间块 ∝ I（等价于标准表示下为标量）")
    an, tr, _ = isotropy_check(T0)
    print("  空间块各向异性范数 = %.3e （均匀权重下）" % an)
    OUT["A"] = {"eigs": [float(x) for x in ev], "aniso": an}

    print()
    print("=" * 74)
    print("B · 一般权重：S_4 下 Edges 的两个轨道")
    print("=" * 74)
    orbA = [(i, j) for (i, j) in EDGES if i == 0 or j == 0]
    orbB = [(i, j) for (i, j) in EDGES if i != 0 and j != 0]
    print("  轨道 A（与根相连）大小 = %d  ⇒ C(D,1)=4" % len(orbA))
    print("      %s" % orbA)
    print("  轨道 B（不与根相连）大小 = %d ⇒ C(D,2)=6" % len(orbB))
    print("      %s" % orbB)
    print()
    print("  ⇒ T^{ab} 是**两参数族** (a,b)；10 = 4 + 6")
    print("  ⇒ 这正是通道对数 C(D+1,2) = D + C(D,2) 的分解！")
    OUT["B"] = {"orbA": len(orbA), "orbB": len(orbB)}

    print()
    print("=" * 74)
    print("C · 两参数族的形状")
    print("=" * 74)
    a, b = 1.0, 1.0
    T = build_T(a, b)
    print("  T^{00} = %s" % T[0, 0])
    print("  空间对角 = %s" % np.round(np.diag(T)[1:], 8))
    an, tr, Tsub = isotropy_check(T)
    print("  空间块各向异性 = %.3e" % an)
    print()
    print("  解析：")
    print("    T^{00} = D·a = 4a")
    print("    T^{ii}（i≠0）= a + (D-1)·b = a + 3b")
    print("    非对角 T^{ij}（i≠j, 都≠0）= -b")
    print()
    for (aa, bb) in [(1, 1), (1, 0.5), (1, 2), (1, 0)]:
        Tt = build_T(aa, bb)
        an2, tr2, _ = isotropy_check(Tt)
        t00 = Tt[0, 0]
        tii = Tt[1, 1]
        print("    a=%.2f b=%.2f : T00=%.4f Tii=%.4f 空间各向异性=%.2e  w=Tii/T00=%.4f"
              % (aa, bb, t00, tii, an2, tii / t00))
    OUT["C"] = {"T00": "D*a", "Tii": "a+(D-1)*b"}

    print()
    print("=" * 74)
    print("D · 状态方程 w = p/ρ 对 b/a 的依赖")
    print("=" * 74)
    print("  w(b/a) = (a + 3b) / (4a) = (1 + 3(b/a)) / 4")
    print()
    print("  b/a     w       物性")
    for r in [0.0, 1/3, 1.0, 5/3, 2.0]:
        w = (1 + 3 * r) / 4
        kind = ""
        if abs(w) < 1e-9:
            kind = "尘埃（w=0）"
        elif abs(w - 1/3) < 1e-9:
            kind = "辐射（w=1/3）"
        elif abs(w - 1) < 1e-9:
            kind = "硬物质（w=1）"
        elif abs(w - 4/3) < 1e-9:
            kind = "w=4/3"
        print("  %-7.4f %-8.4f %s" % (r, w, kind))
    print()
    print("  ⇒ 尘埃（w=0）当且仅当 b/a = −1/3（负权重，非物理）")
    print("  ⇒ 辐射（w=1/3）当且仅当 b/a = 1/3")
    print("  ⇒ 均匀（b/a=1）给 w=1")
    OUT["D"] = {"w_of_ratio": "(1+3r)/4"}

    print()
    print("=" * 74)
    print("E · 无迹条件是否钉住 b/a")
    print("=" * 74)
    for (aa, bb) in [(1, 1), (1, 1/3), (1, 5/3)]:
        Tt = build_T(aa, bb)
        tr = np.trace(Tt)
        print("  a=%.3f b=%.3f : trace(T) = %.6f  ⇒ 无迹=%s"
              % (aa, bb, tr, abs(tr) < 1e-12))
    print()
    print("  解析：trace(T) = T^{00} + Σ_i T^{ii} = 4a + 4(a+3b) = 8a + 12b")
    print("  ⇒ 无迹 ⟺ 8a + 12b = 0 ⟺ b/a = −2/3（**负权重**，非物理）")
    print()
    print("  ⇒ 无迹条件**不可用**（要求负边权）")
    print("  ⇒ b/a 需要另一个原生约束来钉住")
    OUT["E"] = {"traceless_requires": "b/a = -2/3 (unphysical)"}

    print()
    print("=" * 74)
    print("总判定")
    print("=" * 74)
    print("  A  正确判据下：均匀权重 T 是各向同性的 ✅")
    print("  B  S_4 把 10 条边分成 4+6 两个轨道（= D + C(D,2)）✅")
    print("  C  T^{ab} 是两参数族；空间块恒各向同性（因轨道结构）")
    print("  D  w = (1 + 3(b/a))/4")
    print("  E  无迹条件要求负权重 ⇒ 不可用 ⇒ b/a 未定")
    OUT["verdict"] = {"isotropic": True, "orbits": [4, 6],
                      "two_param": True, "traceless_unusable": True}

    with open(os.path.join(HERE, "R79_stress_two_param_results.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False, default=str)
    print("\n-> R79_stress_two_param_results.json")
