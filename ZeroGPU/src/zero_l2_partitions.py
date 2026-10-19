#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_l2_partitions.py —— L2 的不同划分是否给出不同的"世界"？
================================================================
假设（用户）：L2 可以从不同角度再划分，划出不同的世界。

两种读法：
  (A) **同一状态的基变换/读出**（相邻比、Helmert、Laplacian 本征基、C3 基）；
  (B) **不同的粗粒化/分组**（16/8/4、6/18/4、8/4/2 …）。

测试：
  1. 同一形状算子 M 在四种基下的谱：是否不变？
  2. C3 基本身是否把 M 谱分裂成 3 个？
  3. 换成手征算子 P·M 后是否出现 3 个不同本征值？
  4. 不同的粗粒化是否给出不同的层数/周期？它们是不是输入？

结果：
  · M 在四种基下都是 (1/81)·I —— **谱不变，1 个不同本征值**；
  · C3 基不分裂 M（基变换不产生新动力学）；
  · 手征算子 P·M 给 1/81·(1, ω, ω²) —— **3 个不同本征值，但这是换算子，不是换划分**；
  · 16/8/4、6/18/4、8/4/2 给不同层数/周期，但全部是**输入**，没有选择原则。

判定：**划分给"有效世界/读出"，不给"基本世界"；要不同世界，必须换定律（耦合/
手征）或加选择原则（seed/对称性/观测）。**

输出：results/zero_l2_partitions.json、zero_l2_partitions.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_l2_partitions.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zero_sim as zs  # noqa: E402

RES = os.path.join(HERE, "..", "results")


def distinct(ev, tol=1e-9):
    cl = []
    for z in ev:
        if not any(abs(z - c) < tol for c in cl):
            cl.append(z)
    return len(cl)


def main() -> None:
    A = zs.adj_matrix(zs.k4()).astype(float)
    AL = np.linalg.matrix_power(A, 4)
    H = np.array([[1, -1, 0, 0], [1, 1, -2, 0], [1, 1, 1, -3]], float).T
    Q, _ = np.linalg.qr(H)
    Q = Q[:, :3]
    B = np.zeros((3, 4))
    for i in range(3):
        B[i, i] = 1.0
        B[i, i + 1] = -1.0
    u = 0.25

    M = Q.T @ AL @ Q / 81.0                     # Helmert 基
    Sf = B @ Q / u                               # 相邻比基
    L = 3.0 * np.eye(4) - A
    w, V = np.linalg.eigh(L)
    EL = V[:, np.isclose(w, 4.0)]                # K4 的 -1 空间（3 维）
    P = np.zeros((4, 4))
    P[0, 0] = 1.0
    P[2, 1] = P[3, 2] = P[1, 3] = 1.0

    S = EL.T @ Q                                  # Helmert -> Laplacian 坐标
    PL = EL.T @ P @ EL
    wc, Vc = np.linalg.eig(PL)
    SC = Vc.conj().T @ S                          # Helmert -> C3 坐标
    bases = {
        "Helmert": M,
        "adjacent_ratio": Sf @ M @ np.linalg.inv(Sf),
        "laplacian_-1": S @ M @ np.linalg.inv(S),
        "C3_irrep": SC @ M @ np.linalg.inv(SC),   # 同一 M，C3 基下仍是标量
    }
    specs = {k: [complex(z) for z in np.linalg.eigvals(v)] for k, v in bases.items()}

    chiral = EL.T @ (AL @ P) @ EL / 81.0
    chiral_spec = [complex(z) for z in np.linalg.eigvals(chiral)]

    # 粗粒化（来自既有文档的登记）
    coarse = [
        dict(name="16/8/4", layers=[16, 8, 4], period=4, status="输入（zero_sublayers）"),
        dict(name="6/18/4", layers=[6, 18, 4], period=4, status="计数测试/反例"),
        dict(name="8/4/2", layers=[8, 4, 2], period=4, status="L=4 下的输入"),
        dict(name="3 slots", layers=[3], period=4, status="运动学坐标（无亚亚层）"),
    ]

    checks = [
        dict(name="同一 M 在四种基下谱不变（1 个不同本征值）",
             pass_=all(distinct(specs[k]) == 1 for k in specs),
             detail={k: [round(abs(z), 6) for z in specs[k]] for k in specs}),
        dict(name="C3 基不分裂 M（基变换不产生新动力学）",
             pass_=distinct(specs["C3_irrep"]) == 1,
             detail=f"distinct={distinct(specs['C3_irrep'])}"),
        dict(name="手征算子 P·M 出现 3 个不同本征值",
             pass_=distinct(chiral_spec) == 3,
             detail=[complex(round(z.real, 6), round(z.imag, 6)) for z in chiral_spec]),
        dict(name="不同粗粒化给不同层数",
             pass_=len({tuple(c["layers"]) for c in coarse}) == len(coarse),
             detail={c["name"]: c["layers"] for c in coarse}),
        dict(name="粗粒化全部是输入（无选择原则）",
             pass_=all("输入" in c["status"] or "反例" in c["status"]
                       or "坐标" in c["status"] for c in coarse),
             detail={c["name"]: c["status"] for c in coarse}),
    ]

    res = dict(
        M=[[complex(round(z.real, 9), round(z.imag, 9)) for z in row] for row in M],
        specs={k: [[float(z.real), float(z.imag)] for z in v] for k, v in specs.items()},
        chiral=[[float(z.real), float(z.imag)] for z in chiral_spec],
        coarse=coarse,
        checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                 "pass": bool(c["pass_"])} for c in checks],
    )

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0))
    names = list(specs.keys()) + ["chiral P·M"]
    counts = [distinct(specs[k]) for k in specs] + [distinct(chiral_spec)]
    axes[0].bar(range(len(names)), counts,
                color=["#777"] * len(specs) + ["#c33"])
    axes[0].set_xticks(range(len(names)))
    axes[0].set_xticklabels(names, rotation=20, ha="right", fontsize=8)
    axes[0].set_ylabel("distinct eigenvalues")
    axes[0].set_title("partition does not change the spectrum")
    axes[0].grid(alpha=0.3, axis="y")
    for i, c in enumerate(coarse):
        axes[1].bar([i], [max(c["layers"])], width=0.5)
        axes[1].text(i, max(c["layers"]) + 0.4, str(c["layers"]),
                     ha="center", fontsize=7)
    axes[1].set_xticks(range(len(coarse)))
    axes[1].set_xticklabels([c["name"] for c in coarse], fontsize=8)
    axes[1].set_ylabel("max layers")
    axes[1].set_title("different partitions = different layer counts (all inputs)")
    axes[1].grid(alpha=0.3, axis="y")
    fig.tight_layout()
    plot = os.path.join(RES, "zero_l2_partitions.png")
    fig.savefig(plot, dpi=140)
    plt.close(fig)

    def _j(o):
        if isinstance(o, complex):
            return [float(o.real), float(o.imag)]
        return str(o)

    with open(os.path.join(RES, "zero_l2_partitions.json"), "w") as f:
        json.dump(res, f, indent=1, default=_j)

    print("L2 划分测试")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    print(f"  谱：{[distinct(specs[k]) for k in specs]}；手征 distinct={distinct(chiral_spec)}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
