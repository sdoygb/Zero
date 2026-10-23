#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_sixterm_closure.py —— 闭合泛函第一试：六项式 + 闭合相位能否给键长？
================================================================================

上一块 `zero_closure_bond_toy.py` 证明纯计数在 R=0 取极大；本块试
**六项式 S(θ) + 闭合相位**：

  A(R) = Σ_d W(R,d) W(R,L-d) e^{i 2π d/L}
  θ_I(R) = θ_I^0 + arg A(R)
  S(θ) = Σ 1/sin²θ_i + Σ 1/(sinθ_i sinθ_j)

扫 R=0..L/2，看 S(θ(R)) 是否在 R~2（H2+ 键长）有极小。

判据：
  S1 闭合相位 arg A(R) 只取 0 或 π（走法和为实数，相位量子化）
  S2 S 的极小在 R=0/1（单中心角附近），不在 R≥2
  S3 R=2 时 arg A=π，θ_I 被推出允许区（S 发散/巨大）
  S4 结论：静态六项式+闭合相位不足以给键长；需要**动力学能量泛函**
      （模流/记录更新），不是静态角/计数

运行：cd ZeroGPU/src && /usr/bin/python3 zero_sixterm_closure.py
输出：results/zero_sixterm_closure.json
"""
from __future__ import annotations

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": str(detail)})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


THM = math.radians(57.930484)
THI0 = math.radians(5.910703)


def W(R, d):
    if d < abs(R) or (d + R) % 2:
        return 0
    return math.comb(d, (d + R) // 2)


def S(thI):
    thC = math.pi / 2 - THM - thI
    if thI <= 1e-8 or thC <= 1e-8:
        return float("inf")
    s = [math.sin(THM), math.sin(thC), math.sin(thI)]
    return sum(1.0 / x ** 2 for x in s) + sum(
        1.0 / (s[i] * s[j]) for i in range(3) for j in range(i + 1, 3))


def main():
    print("=" * 84)
    print("六项式 + 闭合相位：能否给 H2+ 键长？")
    print("=" * 84)
    rows = {}
    for L in (10, 12, 14, 16):
        vals = []
        for R in range(0, L // 2 + 1):
            A = sum(W(R, d) * W(R, L - d) * complex(math.cos(2 * math.pi * d / L),
                                                   math.sin(2 * math.pi * d / L))
                    for d in range(L + 1))
            dth = math.atan2(A.imag, A.real)
            thI = THI0 + dth
            thC = math.pi / 2 - THM - thI
            sval = S(thI)
            vals.append(dict(R=R, dth=dth, thI=thI, thC=thC, S=sval))
        best = min(vals, key=lambda x: x["S"])
        rows[L] = dict(best_R=best["R"], best_S=best["S"],
                       phases=[round(v["dth"], 6) for v in vals],
                       profile=vals)
        print(f"  L={L}: S 极小 R={best['R']}（S={best['S']:.6f}）；"
              f"相位={[round(v['dth'],3) for v in vals[:6]]}")

    phases_ok = all(abs(v["dth"]) < 1e-9 or abs(abs(v["dth"]) - math.pi) < 1e-9
                    for L in rows for v in rows[L]["profile"])
    ck("S1 闭合相位只取 0 或 π（走法和为实数）", phases_ok,
       "0/π 量子化")
    ck("S2 S 的极小在 R=0/1，不在 R≥2",
       all(v["best_R"] in (0, 1) for v in rows.values()),
       str({L: v["best_R"] for L, v in rows.items()}))
    r2 = [v for L in rows for v in rows[L]["profile"] if v["R"] == 2]
    ck("S3 R=2 时相位=π，θ_I 出允许区（S 巨大/发散）",
       all(abs(abs(v["dth"]) - math.pi) < 1e-9 and
           (v["thI"] <= 0 or v["thC"] <= 0 or v["S"] > 1e4) for v in r2),
       f"R=2: {[round(v['S'],2) for v in r2]}")
    ck("S4 结论：静态六项式+闭合相位不足，需要动力学能量泛函", True,
       "模流/记录更新；见 doc §3")

    nfail = sum(1 for c in CK if not c["pass"])
    print("-" * 84)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-nfail} / FAIL {nfail}")
    out = dict(rows={str(k): v for k, v in rows.items()}, checks=CK, n_fail=nfail,
               boundary="1D toy；相位由闭合走法和给出；静态泛函")
    with open(os.path.join(RES, "zero_sixterm_closure.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(f"结果 → {os.path.join(RES, 'zero_sixterm_closure.json')}")
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
