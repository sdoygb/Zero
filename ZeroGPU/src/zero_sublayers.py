#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""L2 亚亚层分级毁灭 prototype。

结构（用户定）：L2 三个亚层，亚亚层数 16/8/4（底/中/顶）。
毁灭（每 T=4）：底 16 个全灭；中只保留最上面 2 个（灭 6/8）；顶 4 个不动。

这是一个**层级 prototype**，不是忠实 4 顶点走者引擎：
  · 每层用有序亚亚层通道表示，通道权重之和 = 该层总量；
  · 层间有下行流 顶→中→底，顶部有 L1 seed，底部是汇；
  · 归一化总量到 1；毁灭前后分别记录，检验周期 4 的“全灭/半灭/不动”。
输出：results/zero_sublayers.json、results/zero_sublayers.png
"""
from __future__ import annotations
import argparse, json, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
COUNTS = (16, 8, 4)          # 底/中/顶
TOP_KEEP = (0, 0, 2)         # 底灭 0/0（全灭），中留 2，顶留 4（不灭）
CK = []


def ck(name, cond, detail=""):
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))


def profile(n, r):
    if abs(r - 1.0) < 1e-12:
        return np.ones(n) / n
    p = r ** np.arange(n)
    return p / p.sum()


def init_state(profs, totals):
    return [profs[g] * totals[g] for g in range(3)]


def totals(w):
    return np.array([float(x.sum()) for x in w])


def active(w):
    return np.array([int((x > 1e-14).sum()) for x in w])


def destroy(w):
    """底全灭；中只留最上面 2 个；顶不动。"""
    w[0][:] = 0.0
    w[1][2:] = 0.0
    # 顶不动
    return w


def flow(w, profs, k, seed):
    """顶→中→底的下行流；顶部 seed；归一化。"""
    P = totals(w)
    # 顶 → 中
    j32 = k * P[2]
    if P[2] > 0:
        take = w[2] * (j32 / P[2])
        w[2] -= take
        w[1] += j32 * profs[1]
    # 中 → 底
    P = totals(w)
    j21 = k * P[1]
    if P[1] > 0:
        take = w[1] * (j21 / P[1])
        w[1] -= take
        w[0] += j21 * profs[0]
    # L1 seed 到顶
    w[2] += seed * profs[2]
    # 归一化总量到 1（保持比例，避免总量漂移）
    tot = sum(float(x.sum()) for x in w)
    if tot > 0:
        for g in range(3):
            w[g] /= tot
    return w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", type=int, default=200)
    ap.add_argument("--T", type=int, default=4)
    ap.add_argument("--r", type=float, default=0.8)
    ap.add_argument("--k", type=float, default=0.2)
    ap.add_argument("--seed", type=float, default=0.05)
    args = ap.parse_args()

    profs = [profile(n, args.r) for n in COUNTS]
    # 初值：观测 L2 的粗分组（底 p2+p3，中 p1，顶 p0），归一到 1
    totals0 = np.array([0.1034 + 0.0690, 0.2759, 0.5517])
    totals0 /= totals0.sum()
    w = init_state(profs, totals0)

    rows = []
    for gen in range(1, args.steps + 1):
        pre = totals(w)
        destroyed = (gen % args.T == 0)
        if destroyed:
            w = destroy(w)
            post_d = totals(w)
            act_d = active(w)
        else:
            post_d = pre.copy()
            act_d = active(w)
        w = flow(w, profs, args.k, args.seed)
        post_f = totals(w)
        act_f = active(w)
        rows.append({"gen": gen, "pre": pre.tolist(),
                     "post_destroy": post_d.tolist(),
                     "post_flow": post_f.tolist(),
                     "active_destroy": act_d.tolist(),
                     "active_flow": act_f.tolist(),
                     "destroyed": bool(destroyed)})
        if gen <= 12 or destroyed:
            print(f"gen={gen:>3} pre={np.round(pre,4)} "
                  f"post_d={np.round(post_d,4)} post_f={np.round(post_f,4)} "
                  f"act_d={act_d.tolist()}")

    print("=" * 70)
    print("判据")
    print("=" * 70)
    ck("亚亚层数 = 16/8/4", COUNTS == (16, 8, 4), str(COUNTS))
    dd = [r for r in rows if r["destroyed"]]
    ck("每次毁灭：底全灭（总=0）", all(abs(r["post_destroy"][0]) < 1e-15 for r in dd),
       f"n={len(dd)}")
    ck("每次毁灭：中只留最上 2 个", all(r["active_destroy"][1] == 2 for r in dd),
       str([r["active_destroy"][1] for r in dd[:5]]))
    ck("每次毁灭：顶 4 个不动",
       all(r["active_destroy"][2] == 4 and abs(r["post_destroy"][2] - r["pre"][2]) < 1e-12
           for r in dd), "top unchanged")
    ck("毁灭后流/seed 恢复，无 NaN",
       all(np.all(np.isfinite(r["post_flow"])) and min(r["post_flow"]) >= 0 for r in rows),
       "finite/nonnegative")
    ck("周期 4 的毁灭事件", all(r["destroyed"] == (r["gen"] % 4 == 0) for r in rows),
       "every T=4")

    # 尾段均值（最后 20 代）
    tail = rows[-20:]
    mean_post = np.mean([r["post_flow"] for r in tail], axis=0)
    print("尾段 post_flow 均值（底,中,顶）:", np.round(mean_post, 6))

    # 图
    gens = np.array([r["gen"] for r in rows])
    pd = np.array([r["post_destroy"] for r in rows])
    pf = np.array([r["post_flow"] for r in rows])
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].plot(gens, pd[:, 0], label="bottom 16 (post-destroy)", lw=1.0)
    ax[0].plot(gens, pd[:, 1], label="middle 8 (post-destroy)", lw=1.0)
    ax[0].plot(gens, pd[:, 2], label="top 4 (post-destroy)", lw=1.0)
    ax[0].set_xlabel("generation"); ax[0].set_ylabel("share"); ax[0].set_title("(a) just after destruction")
    ax[0].legend(fontsize=8)
    ax[1].plot(gens, pf[:, 0], label="bottom", lw=1.0)
    ax[1].plot(gens, pf[:, 1], label="middle", lw=1.0)
    ax[1].plot(gens, pf[:, 2], label="top", lw=1.0)
    ax[1].set_xlabel("generation"); ax[1].set_ylabel("share"); ax[1].set_title("(b) after flow/seed")
    ax[1].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout()
    fp = os.path.join(RES, "zero_sublayers.png")
    fig.savefig(fp, dpi=130)

    out = {"counts": list(COUNTS), "top_keep": list(TOP_KEEP),
           "T": args.T, "r": args.r, "k": args.k, "seed": args.seed,
           "rows": rows, "checks": CK,
           "n_fail": sum(1 for c in CK if not c["pass"]),
           "tail_mean_post_flow": mean_post.tolist(),
           "anchor": "integer counts; no dimensionful anchor"}
    with open(os.path.join(RES, "zero_sublayers.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print("写入", os.path.join(RES, "zero_sublayers.json"), "和", fp)


if __name__ == "__main__":
    main()
