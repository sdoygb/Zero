#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
area_rank_parity_check.py —— 面积秩 A_h 的奇偶定理：核验
============================================================
对象：Z048 §8 第 2 条「A_h 奇偶 ↔ 层内存在零面积词」，其定理形式。

定义（Z048 (2.1)）：
  w ∈ {±1}^L, Σw=0;  s_k = Σ_{i≤k} w_i,  h(w)=max_k |s_k|;
  a(w)=Σ_{k<L} s_k;  A_h = #{a(w): h(w)=h}.

定理：A_h 为奇数 ⟺ 层 h 内存在 a(w)=0 的词。
证明：a(-w)=-a(w)（Z048 §1.1）、h(-w)=h(w)（Z047 §1.2），
      故层内面积谱在 a↦-a 下不变；有限对称集含 0 时元素数为奇数，否则为偶数。

核验：
  A. Z₂：a(-w)=-a(w)、h(-w)=h(w)（逐词，L=4..16）
  B. 每层：A_h 奇 ⟺ 0 在谱中（L=4..20）
  C. L≡2 (mod 4) 每层 A_h 偶；L≡0 (mod 4) 低层可奇

判据方向固定：PASS = 逐层判据全成立；FAIL = 该篇要报的判定。
"""
from __future__ import annotations
import os
import json
from itertools import combinations
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []


def P(*a):
    s = " ".join(str(x) for x in a)
    LOG.append(s)
    print(s)


def words(L):
    """零和词：选 L/2 个位置为 +1，其余 -1。"""
    for pos in combinations(range(L), L // 2):
        w = [-1] * L
        for i in pos:
            w[i] = 1
        yield w


def layer_area(w):
    s = 0
    h = 0
    a = 0
    for k, x in enumerate(w):
        s += x
        if abs(s) > h:
            h = abs(s)
        a += s
    return h, a


def main():
    P("=" * 80)
    P("面积秩 A_h 的奇偶定理：核验 (area_rank_parity_check)")
    P("=" * 80)

    checks = []
    bad_z2 = 0
    bad_parity = 0
    table = []
    for L in range(4, 21, 2):
        # 逐词 Z2 + 分层
        spec = {}  # h -> set of areas
        for w in words(L):
            h, a = layer_area(w)
            h2, a2 = layer_area([-x for x in w])
            if h2 != h or a2 != -a:
                bad_z2 += 1
            spec.setdefault(h, set()).add(a)
        rows = []
        for h in sorted(spec):
            S = spec[h]
            A = len(S)
            has0 = 0 in S
            if (A % 2 == 1) != has0:
                bad_parity += 1
            rows.append((h, A, has0, A % 2 == 1 == has0))
        table.append((L, rows))
        P(f"L={L}: " + "  ".join(f"h{h}:A={A},0={'Y' if z else 'N'}" for h, A, z, _ in rows))

    P(f"\nZ2 违例 = {bad_z2}；奇偶判据违例 = {bad_parity}")

    def chk(name, cond, detail=""):
        checks.append({"name": name, "pass": bool(cond), "detail": detail})
        P(("  [v] " if cond else "  [x] ") + name + (f"   {detail}" if detail else ""))

    P("\n[A] Z₂ 对称")
    chk("A  a(-w)=-a(w)、h(-w)=h(w)（L=4..16 逐词）", bad_z2 == 0, f"{bad_z2}")
    P("\n[B] 奇偶定理")
    chk("B  A_h 奇 ⟺ 0 在谱中（L=4..20 逐层）", bad_parity == 0, f"{bad_parity}")
    P("\n[C] 模 4")
    mod4 = True
    for L, rows in table:
        if L % 4 == 2:
            if any(A % 2 == 1 for _, A, _, _ in rows):
                mod4 = False
    chk("C  L≡2 (mod 4) 全部 A_h 偶", mod4)

    npass = sum(1 for c in checks if c["pass"])
    P("\n" + "=" * 80)
    P(f"核验：通过 {npass} / {len(checks)}")
    P("=" * 80)

    out = {"bad_z2": bad_z2, "bad_parity": bad_parity,
           "table": [[L, [[h, A, z] for h, A, z, _ in rows]] for L, rows in table],
           "checks": checks, "pass": npass, "total": len(checks)}
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "area_rank_parity_check.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "area_rank_parity_check.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    P("\n已存 area_rank_parity_check.json / .txt")


if __name__ == "__main__":
    main()
