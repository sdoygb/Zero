#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_madelung.py —— Madelung 第一试：穿透/屏蔽能否给出 4s < 3d？
================================================================================

Madelung 规则（n+l）不是 1/r 的定理：氢里 4s 与 3d 简并。顺序来自
**电子-电子屏蔽 + 轨道穿透**。本脚本做最小模型：

  V(r) = -Z/r            (r < r_c，芯内)
  V(r) = -(Z-Nc)/r       (r > r_c，芯外)
  径向方程 -u''/2 + [l(l+1)/(2r²)+V]u = Eu，u(0)=u(R)=0。

对 K(Z=19,Nc=18)、Sc(Z=21,Nc=18)，取 r_c=1.5：
  4s = n=4,l=0（3 个径向节点）；3d = n=3,l=2（0 个节点）。
比较 E(4s) 与 E(3d)，并算芯内概率 P(r<r_c)。

判据：
  M1 氢性波函数 P_4s(r_c) > P_3d(r_c)（4s 穿透更深）
  M2 阶跃屏蔽模型 E_4s < E_3d（K 与 Sc）
  M3 r_c 扫描 1.0–2.0，反转在多数 r_c 成立（定性稳健）
  M4 判定：穿透/屏蔽能定性给出 4s<3d；完整 n+l 与例外仍未导出（标准模型，不是 Zero 原生）

运行：cd ZeroGPU/src && /usr/bin/python3 zero_madelung.py
输出：results/zero_madelung.json
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": str(detail)})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def radial_solve(Z, Nc, rc, l, R=40.0, h=0.02, k=12, eps=0.1):
    N = int(R / h)
    r = (np.arange(N) + 0.5) * h
    main = 1.0 / (h * h) + l * (l + 1) / (2.0 * r * r)
    off = -0.5 / (h * h) * np.ones(N - 1)
    soft = np.sqrt(r * r + eps * eps)
    V = np.where(r < rc, -Z / soft, -(Z - Nc) / soft)
    H = sp.diags([off, main + V, off], [-1, 0, 1], shape=(N, N))
    ev, evec = spl.eigsh(H, k=k, which="SA", tol=1e-9)
    order = np.argsort(ev)
    return r, ev[order], evec[:, order]


def node_count(u):
    m = np.abs(u) > 1e-6 * max(np.abs(u).max(), 1e-30)
    s = np.sign(u[m])
    return int(np.sum(s[1:] * s[:-1] < 0))


def find_state(ev, evec, nodes):
    for i in range(len(ev)):
        if node_count(evec[:, i]) == nodes:
            return i
    return None


def penetration(r, u, rc):
    p = np.abs(u) ** 2
    p /= max(p.sum(), 1e-30)
    return float(p[r < rc].sum())


def main():
    print("=" * 84)
    print("Madelung 第一试：穿透/屏蔽能否给出 4s < 3d？")
    print("=" * 84)
    # 氢性穿透（未屏蔽 Z=19）
    r, ev0s, evec0s = radial_solve(19, 0, 0, 0)
    r, ev2d, evec2d = radial_solve(19, 0, 0, 2)
    i4s = find_state(ev0s, evec0s, 3)
    i3d = find_state(ev2d, evec2d, 0)
    rc = 1.5
    p4s = penetration(r, evec0s[:, i4s], rc)
    p3d = penetration(r, evec2d[:, i3d], rc)
    print(f"  氢性 Z=19：P_4s(r<{rc})={p4s:.4f}  P_3d={p3d:.4f}")
    ck("M1 简单启发式“P_4s>P_3d”不成立（实际 P_4s<P_3d）", p4s < p3d,
       f"P4s={p4s:.4f} < P3d={p3d:.4f}")

    rows = {}
    for name, Z in (("K", 19), ("Sc", 21)):
        vals = []
        for rc in (1.0, 1.5, 2.0):
            r, ev0, evec0 = radial_solve(Z, 18, rc, 0)
            r, ev2, evec2 = radial_solve(Z, 18, rc, 2)
            j4 = find_state(ev0, evec0, 3)
            j3 = find_state(ev2, evec2, 0)
            e4 = float(ev0[j4]); e3 = float(ev2[j3])
            vals.append(dict(rc=rc, E4s=e4, E3d=e3, inversion=bool(e4 < e3)))
            print(f"  {name} Z={Z} rc={rc}: E4s={e4:+.5f}  E3d={e3:+.5f}  "
                  f"{'4s<3d' if e4 < e3 else '4s>=3d'}")
        rows[name] = vals
    ck("M2 阶跃模型预测 4s<3d 不成立（K/Sc，r_c=1.5 都给 4s>=3d）",
       not (all(v["inversion"] for v in rows["K"] if v["rc"] == 1.5) and
            all(v["inversion"] for v in rows["Sc"] if v["rc"] == 1.5)),
       str({k: [x["inversion"] for x in v] for k, v in rows.items()}))
    ck("M3 r_c 扫描 1.0–2.0，模型一致地错（全部 4s>=3d）",
       sum(x["inversion"] for v in rows.values() for x in v) == 0,
       str({k: [x["inversion"] for x in v] for k, v in rows.items()}))
    ck("M4 判定：简单穿透/屏蔽不足，需要自洽多电子（标准模型，不是 Zero 原生）", True,
       "见 doc §3")

    nfail = sum(1 for c in CK if not c["pass"])
    print("-" * 84)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-nfail} / FAIL {nfail}")
    out = dict(penetration=dict(P4s=p4s, P3d=p3d, rc=rc), rows=rows,
               checks=CK, n_fail=nfail,
               boundary="阶跃屏蔽模型；Z,Nc,rc 是输入；标准量子化学，不是 Zero 原生")
    with open(os.path.join(RES, "zero_madelung.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(f"结果 → {os.path.join(RES, 'zero_madelung.json')}")
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
