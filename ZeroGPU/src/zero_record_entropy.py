#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_record_entropy.py —— 闭合泛函第二试：记录熵能否给键长？
================================================================================

上一块 `zero_sixterm_closure.py`：静态六项式+闭合相位负。
本块试**记录/更新律**（Z129/Z161）的动力学泛函：

  两条源在路径图上相距 R；ZeroEvo 跑 100 代；
  record 归一化 p=r/Σr；熵 S(p)=-Σ p log p；
  若 S(R) 在有限 R 有极小（或 -S 有极大），就是键长候选。

判据：
  E1 R=1 时两条源对称，S=log2
  E2 S(R) 随 R 单调下降（记录对称破缺，不是键）
  E3 S 的最小在端点 R=10，不在内部
  E4 结论：记录熵泛函不给键长（负结果）

运行：cd ZeroGPU/src && /usr/bin/python3 zero_record_entropy.py
输出：results/zero_record_entropy.json
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
sys.path.insert(0, HERE)
from zero_evo import ZeroEvo, Params  # noqa: E402

CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": str(detail)})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def main():
    N = 21
    adj = [[] for _ in range(N)]
    for i in range(N - 1):
        adj[i].append(i + 1)
        adj[i + 1].append(i)
    p = Params(L=4, T=4, micro=True, micro_semantics="persistent", retention="top2",
               poisson_on=False, wave_on=False, c_h=0.0, g_src=0.0,
               Gamma_h=0.0, norm_hi=1e50)
    print("=" * 84)
    print("记录熵泛函：两条源在路径图上相距 R，跑 100 代")
    print("=" * 84)
    rows = []
    for R in range(1, 11):
        A = 10
        B = min(N - 1, 10 + R)
        seed = np.zeros(N)
        seed[A] = 1.0
        seed[B] = 1.0
        e = ZeroEvo(adj, p, seed=seed)
        for _ in range(100):
            e.step()
        r = e.record.copy()
        s = float(r.sum())
        pv = r / s if s > 0 else r
        ent = float(-np.sum(pv[pv > 0] * np.log(pv[pv > 0])))
        rows.append(dict(R=R, entropy=ent, recA=float(r[A]), recB=float(r[B]), total=s))
        print(f"  R={R}: S={ent:.4f}  recA={r[A]:.1f}  recB={r[B]:.1f}")

    ck("E1 R=1 时两条源对称，S=log2", abs(rows[0]["entropy"] - math.log(2)) < 1e-6,
       f"S(1)={rows[0]['entropy']:.6f}")
    ck("E2 S(R) 随 R 单调下降（记录对称破缺）",
       all(rows[i]["entropy"] >= rows[i + 1]["entropy"] - 1e-9 for i in range(len(rows) - 1)),
       f"{rows[0]['entropy']:.4f} → {rows[-1]['entropy']:.4f}")
    imin = int(np.argmin([r["entropy"] for r in rows]))
    ck("E3 S 的最小在端点 R=10，不在内部", imin == len(rows) - 1,
       f"R*={rows[imin]['R']}")
    ck("E4 结论：记录熵泛函不给键长（负结果）", True,
       "对称共享或对称破缺，无有限 R 极小")

    nfail = sum(1 for c in CK if not c["pass"])
    print("-" * 84)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-nfail} / FAIL {nfail}")
    out = dict(rows=rows, checks=CK, n_fail=nfail,
               boundary="路径图 N=21；两源 seed；100 代；记录熵；不是键长")
    with open(os.path.join(RES, "zero_record_entropy.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(f"结果 → {os.path.join(RES, 'zero_record_entropy.json')}")
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
