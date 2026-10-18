#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""16/8/4 纤维塔：top-down 核与 3 个 Z2xZ2 字符。

16 = 4 base x (Z2 x Z2)；顶4 -> 纤维4 -> 底16。
比较两种 top-down 核：
  global q：所有 base 用同一个纤维分布 -> 3 字符秩1；
  base-dependent Q：每个 base 的纤维分布不同 -> 3 字符秩3。
输出：results/zero_fiber_topdown.json
"""
from __future__ import annotations
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
H2 = np.array([[1,1],[1,-1]], float)
H4 = np.kron(H2, H2) / 2.0
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")


def mode_rank(X):
    """X: 4 base x 4 fiber; 返回 3 个非平凡字符在 base 上的 4x3 矩阵的秩。"""
    Y = X @ H4.T          # base x character
    M = Y[:, 1:]          # 3 nontrivial
    s = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(s > 1e-10 * max(s.max(), 1e-30))), s


def main():
    t = np.array([1.0, 0.8, 0.6, 0.4])
    q = np.array([0.4, 0.3, 0.2, 0.1])
    Q = np.array([
        [0.4, 0.3, 0.2, 0.1],
        [0.1, 0.4, 0.3, 0.2],
        [0.2, 0.1, 0.4, 0.3],
        [0.3, 0.2, 0.1, 0.4],
    ])
    Xg = np.outer(t, q)
    Xd = np.diag(t) @ Q
    rg, sg = mode_rank(Xg)
    rd, sd = mode_rank(Xd)
    print("global q  char singular values", np.round(sg, 6), "rank", rg)
    print("base-dep Q char singular values", np.round(sd, 6), "rank", rd)
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    idx = np.arange(3); w = 0.38
    ax.bar(idx - w/2, sg, w, label="global q (rank 1)")
    ax.bar(idx + w/2, sd, w, label="base-dependent Q (rank 3)")
    ax.set_xticks(idx); ax.set_xticklabels(["char 1", "char 2", "char 3"])
    ax.set_ylabel("singular value"); ax.legend(fontsize=8); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(RES, "zero_fiber_topdown.png"), dpi=130)
    ck = [
        {"name": "global q 的 3 字符秩=1", "pass": rg == 1, "detail": str(rg)},
        {"name": "base-dependent Q 的 3 字符秩=3", "pass": rd == 3, "detail": str(rd)},
        {"name": "Q 每行都是全支撑", "pass": bool(np.all(Q > 0)), "detail": ""},
    ]
    for c in ck:
        print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']}" + (f" —— {c['detail']}" if c['detail'] else ""))
    out = {"global_q_rank": rg, "global_q_s": sg.tolist(),
           "basedep_Q_rank": rd, "basedep_Q_s": sd.tolist(),
           "q": q.tolist(), "Q": Q.tolist(), "checks": ck,
           "n_fail": sum(1 for c in ck if not c["pass"]),
           "note": "structural: base-dependent connection needed for rank3"}
    with open(os.path.join(RES, "zero_fiber_topdown.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(ck)} 条：PASS {len(ck)-out['n_fail']} / FAIL {out['n_fail']}")


if __name__ == "__main__":
    main()
