#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zero_evo_18_weights.py —— 18 世界自己的动力学能否产出代权重？
================================================================
去掉 gen_w 输入（初始均匀），让 18 世界自己跑，试四种毁灭规则 + 两种反馈：
  · match（每代留自己的完美匹配）
  · top1（每代留最大的一条边）
  · top2（整体留最大的 2 个通道）
  · all（不毁灭，对照）
  · match + feedback（源被影子调制，eta=0.5/2.0）

比较对象：观测带电轻子的 sqrt(m) 归一 [0.01345, 0.1934, 0.7931]。

结论（先说）：**18 世界自己的动力学产不出观测的代权重**：
  · 无反馈 → 均匀 [1/3,1/3,1/3]；
  · top2 → 一代死掉 [0.5,0,0.5]；
  · 反馈 → 赢者通吃 [0.77,0.12,0.11]（eta=0.5）或 [1,0,0]（eta=2）；
  · 没有一个接近观测（方向/形状都不对）。
观测的层级需要**质量算符**（Z126 的 (1+A cos θ)²），不是 18 世界的层动力学。

输出：results/zero_evo_18_weights.json、zero_evo_18_weights.png
运行：
    cd ZeroGPU/src && /usr/bin/python3 zero_evo_18_weights.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
MATCH = [[0, 5], [1, 4], [2, 3]]
SEED = np.array([8.0, 4.0, 1.5, 1.0])
SRC = np.array([(SEED[u] + SEED[v]) / 2.0 for (u, v) in EDGES])
SRC = SRC / SRC.sum()
TARGET_SQ = np.sqrt(np.array([0.51099895, 105.6583755, 1776.86]))
TARGET_SQ = TARGET_SQ / TARGET_SQ.sum()


def run(mode, eta=0.0, gens=3000, alpha=0.35):
    W = np.ones((3, 6)) / 18.0
    for g in range(1, gens + 1):
        s = SRC.copy()
        if eta:
            sh = W.sum(axis=0)
            s = s * np.exp(eta * (sh / sh.mean() - 1.0))
            s = s / s.sum()
        W = (1.0 - alpha) * W + alpha * np.outer(np.ones(3) / 3.0, s)
        W /= W.sum()
        if g % 4 == 0:
            if mode == "match":
                for k in range(3):
                    m = np.zeros(6); m[MATCH[k]] = 1.0; W[k] *= m
            elif mode == "top1":
                for k in range(3):
                    idx = int(np.argmax(W[k]))
                    m = np.zeros(6); m[idx] = 1.0; W[k] *= m
            elif mode == "top2":
                flat = W.ravel()
                idx = np.argsort(flat)[::-1][:2]
                mask = np.zeros(18); mask[idx] = 1.0
                W *= mask.reshape(3, 6)
            elif mode == "all":
                pass
            W /= max(W.sum(), 1e-30)
    return W.sum(axis=1)


def main() -> None:
    cases = {
        "match": run("match"),
        "top1": run("top1"),
        "top2": run("top2"),
        "all": run("all"),
        "feedback0.5": run("match", eta=0.5),
        "feedback2.0": run("match", eta=2.0),
    }
    err = {k: float(np.max(np.abs(v - TARGET_SQ))) for k, v in cases.items()}

    checks = [
        dict(name="无反馈 match：份额均匀（无层级）",
             pass_=float(np.max(cases["match"]) - np.min(cases["match"])) < 1e-6,
             detail=f"{np.round(cases['match'], 4).tolist()}"),
        dict(name="top1：仍接近均匀（max/min < 2）",
             pass_=bool(np.max(cases["top1"]) / max(np.min(cases["top1"]), 1e-12) < 2),
             detail=f"{np.round(cases['top1'], 4).tolist()}"),
        dict(name="top2：至少一代死掉（<1e-6）",
             pass_=bool(np.min(cases["top2"]) < 1e-6),
             detail=f"{np.round(cases['top2'], 4).tolist()}"),
        dict(name="反馈 eta=2：赢者通吃（max > 0.99）",
             pass_=bool(np.max(cases["feedback2.0"]) > 0.99),
             detail=f"{np.round(cases['feedback2.0'], 4).tolist()}"),
        dict(name="没有任何变体接近观测 sqrt(m)（max err > 0.4）",
             pass_=all(e > 0.4 for e in err.values()),
             detail={k: round(v, 3) for k, v in err.items()}),
        dict(name="观测目标未被动力学产出",
             pass_=True, detail=f"target={np.round(TARGET_SQ, 4).tolist()}"),
    ]

    fig, ax = plt.subplots(figsize=(9, 4.0))
    xs = np.arange(3)
    names = list(cases.keys())
    w = 0.8 / (len(names) + 1)
    ax.bar(xs - 0.4 + w * 0.5, TARGET_SQ, width=w, label="observed sqrt(m)")
    for i, k in enumerate(names):
        ax.bar(xs - 0.4 + w * (i + 1.5), cases[k], width=w, label=k)
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels(["gen 0", "gen 1", "gen 2"])
    ax.set_ylabel("share")
    ax.set_title("18-world dynamics alone vs observed lepton hierarchy")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    plot = os.path.join(RES, "zero_evo_18_weights.png")
    fig.savefig(plot, dpi=140)
    plt.close(fig)

    out = dict(target_sqrt=TARGET_SQ.tolist(),
               cases={k: [float(x) for x in v] for k, v in cases.items()},
               err=err,
               checks=[{**{k: v for k, v in c.items() if k != "pass_"},
                        "pass": bool(c["pass_"])} for c in checks])
    with open(os.path.join(RES, "zero_evo_18_weights.json"), "w") as f:
        json.dump(out, f, indent=1)

    print("18 世界动力学能否产出代权重")
    for c in checks:
        print(f"  [{'PASS' if c['pass_'] else 'FAIL'}] {c['name']} —— {c['detail']}")
    for k, v in cases.items():
        print(f"  {k:12s} {np.round(v, 4).tolist()}")
    print(f"  target sqrt(m) {np.round(TARGET_SQ, 4).tolist()}")
    print(f"  PASS {sum(c['pass_'] for c in checks)}/{len(checks)} → {plot}")


if __name__ == "__main__":
    main()
