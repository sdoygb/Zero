"""
self_interaction_scaling.py --- 自相互作用的标度：D=4 是**临界维数**吗？（任务 3 候选）

**为什么测这个**
    AUDIT_missing_graph.md §4：补回 Γ 之后，"**合并就是相互作用**"（不同词落到同一顶点
    ⟹ 整数重数 ⟹ 相互作用）。所以"相互作用有多强"= **世界线自交/相撞有多少**：
        I(n) = Σ_{顶点 s} C(k_s, 2),   k_s = 该世界线访问 s 的次数
    随机游走的经典标度：I(n) ~ n² / R(n)^D ~ n^{2-D/2}。
        D=1: n^{3/2} ↑   D=2: n¹ ↑   D=3: n^{1/2} ↑   **D=4: n⁰ = log n（临界）**   D≥5: → 常数 ↓
    ⟹ **D=4 是"相互作用相关／无关"的分界**（自避行走与 Edwards 模型的上临界维数）。

    这正是任务 3 要问的："最小系统里有没有任何 Z0 原生的泛函在 D=4 处取峰/变号？"
    —— 峰值没有，但**变号**有：指数 α(D)=2-D/2 在 D=4 处过零。

**测量**
    · α(D)：I(n) ~ n^{α} 的拟合指数（D=1..6）
    · 访问点数 S(n) ~ n^{νD}（对照：D=1 时 1/2，D≥2 时 1）
    · 相交概率 P(I>0) 随 n 的行为
    · 每条世界线的 I(n) 分布（是否重尾）

**边界**：D 是输入；随机游走在 ℤ^D 上（无限格，坐标不取模，越界样本剔除并计数）。
等级：【数值证据】。

用法：/usr/bin/python3 self_interaction_scaling.py
输出：results/self_interaction.json
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)


def collisions(D, n, W, seed=0):
    """
    返回 (I_mean, I_std, S_mean, out_of_bounds_frac)：
      I = Σ_s C(k_s,2)（自相交对数），S = 访问过的不同顶点数。
    坐标不取模；|coord| 超过 R=6*sqrt(n)+20 的样本剔除。
    """
    rng = np.random.default_rng(seed)
    step = rng.integers(0, 2 * D, size=(W, n))
    axis = step // 2
    sign = np.where(step % 2 == 0, 1, -1).astype(np.int64)
    delta = np.zeros((W, n, D), np.int64)
    rows = np.arange(W)[:, None]
    delta[rows, np.arange(n)[None, :], axis] = sign
    coords = np.zeros((W, n + 1, D), np.int64)
    np.cumsum(delta, axis=1, out=coords[:, 1:, :])

    R = int(6 * np.sqrt(n)) + 20
    bad = (np.abs(coords) > R).any(axis=(1, 2))
    if bad.all():
        return None
    coords = coords[~bad]
    W2 = coords.shape[0]
    base = 2 * R + 1
    idx = coords[..., 0] + R
    for d in range(1, D):
        idx = idx * base + (coords[..., d] + R)
    idx.sort(axis=1)
    # 每行的 run 长度 → Σ C(k,2) 与不同点数
    same = np.diff(idx, axis=1) == 0
    I = np.zeros(W2, np.int64)
    S = np.ones(W2, np.int64)
    # 用"变化点"计数：k_s = run 长度
    newrun = np.ones((W2, idx.shape[1]), bool)
    newrun[:, 1:] = ~same
    S = newrun.sum(axis=1)
    # run 长度：对每个 run 起点累计
    run_len = np.zeros((W2, idx.shape[1]), np.int64)
    # 逐列扫描（列数 = n+1，可接受）
    cur = np.zeros(W2, np.int64)
    for j in range(idx.shape[1]):
        cur = np.where(newrun[:, j], 1, cur + 1)
        run_len[:, j] = cur
    # 每个 run 的最终长度在下一个 run 起点前一列
    endmask = np.zeros_like(newrun)
    endmask[:, :-1] = newrun[:, 1:]
    endmask[:, -1] = True
    kk = run_len[endmask]
    rows_idx = np.repeat(np.arange(W2), endmask.sum(axis=1))
    I = np.zeros(W2, np.int64)
    np.add.at(I, rows_idx, kk * (kk - 1) // 2)
    return (float(I.mean()), float(I.std()), float(S.mean()), float(bad.mean()), W2)


def two_walk_intersection(D, ns, W=2000, seed=0):
    """
    两条独立世界线（同起点，长度各 n）在 t≥1 是否相交。
    经典判据：P→1 ⟺ D≤4；D≥5 有正概率永不相交；D=4 临界（1-P ~ c/log n）。
    """
    rng = np.random.default_rng(seed)
    rows = []
    for n in ns:
        step = rng.integers(0, 2 * D, size=(2 * W, n))
        axis = step // 2
        sign = np.where(step % 2 == 0, 1, -1).astype(np.int64)
        delta = np.zeros((2 * W, n, D), np.int64)
        rr = np.arange(2 * W)[:, None]
        delta[rr, np.arange(n)[None, :], axis] = sign
        coords = np.zeros((2 * W, n + 1, D), np.int64)
        np.cumsum(delta, axis=1, out=coords[:, 1:, :])
        R = int(6 * np.sqrt(n)) + 20
        bad = (np.abs(coords) > R).any(axis=(1, 2))
        base = 2 * R + 1
        idx = coords[..., 0] + R
        for d in range(1, D):
            idx = idx * base + (coords[..., d] + R)
        A = idx[0::2, 1:]      # 世界线 1（t>=1）
        B = idx[1::2, 1:]      # 世界线 2（t>=1）
        hit = np.zeros(W, bool)
        for i in range(W):
            if bad[2 * i] or bad[2 * i + 1]:
                continue
            hit[i] = np.intersect1d(A[i], B[i], assume_unique=False).size > 0
        rows.append({"n": n, "P_intersect": round(float(hit.mean()), 5),
                     "one_minus_P": round(float(1 - hit.mean()), 5), "pairs": int(W)})
    return rows


def main(argv):
    Ds = [1, 2, 3, 4, 5, 6]
    ns = [25, 50, 100, 200, 400, 800, 1600, 3200]
    W = 4000
    out = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
           "definition": "I(n)=Σ_s C(k_s,2)（同一条世界线的自相交对数）；S(n)=访问顶点数",
           "prediction": ("①自交强度 I(n)：瞬态行走（D≥3）每步相遇概率可和 ⟹ I(n)~c_D·n，"
                          "**领头项对所有 D≥2 都是线性**；差别在**次领头修正**：D<4 是幂律 n^{2-D/2}，"
                          "**D=4 是对数（临界）**，D>4 是常数。"
                          "②两条独立世界线是否**必相交**：D≤4 概率→1，D≥5 概率<1（D=4 临界，1-P~c/log n）"),
           "runs": []}
    print("=" * 96)
    print(f"{'D':>2} " + " ".join(f"{('n=%d' % n):>10}" for n in ns) + f"{'   α(I)':>9}{'  α(S)':>8}")
    for D in Ds:
        Is, Ss, alphas = [], [], []
        row = {"D": D, "points": []}
        for n in ns:
            r = collisions(D, n, W, seed=D * 1000 + n)
            if r is None:
                continue
            I_m, I_s, S_m, badfrac, W2 = r
            row["points"].append({"n": n, "I_mean": I_m, "I_std": I_s, "S_mean": S_m,
                                  "kept_walkers": W2, "out_of_bounds_frac": round(badfrac, 6)})
            Is.append(I_m)
            Ss.append(S_m)
        nn = np.array([p["n"] for p in row["points"]], float)
        Is = np.array(Is, float)
        Ss = np.array(Ss, float)
        good = Is > 0
        aI = float(np.polyfit(np.log(nn[good]), np.log(Is[good]), 1)[0]) if good.sum() >= 4 else None
        aS = float(np.polyfit(np.log(nn), np.log(Ss), 1)[0]) if len(nn) >= 4 else None
        row["alpha_I"] = None if aI is None else round(aI, 4)
        row["alpha_I_predicted"] = round(2 - D / 2, 4)
        row["alpha_S"] = None if aS is None else round(aS, 4)
        out["runs"].append(row)
        print(f"{D:>2} " + " ".join(f"{v:>10.3g}" for v in Is) +
              f"{('%.4f' % aI) if aI is not None else '   -  ':>9}{aS:>8.4f}", flush=True)
    # ---- ② 两条世界线必相交？----
    print("=" * 96)
    print("② 两条独立世界线（同起点）在 t≥1 是否相交：经典判据 D≤4 必交 / D≥5 未必（D=4 临界）")
    ns2 = [50, 100, 200, 400, 800, 1600, 3200]
    out["two_walk"] = []
    print(f"{'D':>2} " + " ".join(f"{('n=%d' % n):>9}" for n in ns2))
    for D in Ds:
        rows2 = two_walk_intersection(D, ns2, W=1500, seed=77 + D)
        out["two_walk"].append({"D": D, "rows": rows2})
        print(f"{D:>2} " + " ".join(f"{r['P_intersect']:>9.4f}" for r in rows2), flush=True)
    print("  读法：P 单调升向 1 ⟹ 相互作用**必然**（D≤4）；P 明显停在 <1 ⟹ 可以**永不相交**（D≥5）")
    print("        D=4 的逼近最慢（1-P ~ c/log n）—— 它正是临界维数。")

    with open(os.path.join(OUT, "self_interaction.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("=" * 96)
    print(f"{'D':>2} {'α(I) 实测':>10} {'α=2-D/2':>9}  读法")
    for r in out["runs"]:
        a = r["alpha_I"]
        tag = ("相互作用**相关**（自交随 n 增长）" if (a or 0) > 0.05 else
               "**临界**（自交 ~ log n，相互作用边缘）" if abs(a or 0) <= 0.05 else
               "相互作用**无关**（自交饱和）")
        print(f"{r['D']:>2} {a:>10} {r['alpha_I_predicted']:>9}  {tag}")
    print("\n若 α 在 D=4 过零 ⟹ **D=4 是自相互作用的相关性边界**（Z0 原生泛函，非人为选择）。")


if __name__ == "__main__":
    main(sys.argv[1:])
