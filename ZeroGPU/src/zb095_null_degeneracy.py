#!/usr/bin/env python3
"""P2 的对照：随机密码表下，简并度是否也集中在 5 档？

检验思路：如果"简并度只出现 5 档 {1,2,3,4,6}"是密码表的结构约束（不是巧合），
则随机把 64 个密码子分成 20 个氨基酸组 + 3 个终止（保持组大小分布不变）时，
简并度档数应同样少；反之随机表的档数更多。

两种随机化：
  · 置换：保持 20 个氨基酸的简并度集合 {1,2,3,4,6} 不变，随机指派到氨基酸——档数必然相同，只检验"与非氨基酸的耦合"；
  · 完全随机分配：64 个密码子随机分给 21 个类别（20 氨基酸 + 1 终止），组大小随机——档数会变多。
"""
import json
import random
from collections import Counter
from pathlib import Path

OUT = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb095_null_degeneracy.json")
N_SIM = 20000


def levels(assign):
    d = Counter(assign.values())
    d.pop("*", None)
    return len(set(d.values()))


def main():
    ext = json.load(open("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb095_codon_degeneracy.json"))
    dist = {int(k): v for k, v in ext["degeneracy_distribution"].items()}
    obs_levels = ext["n_levels"]
    print(f"实测：简并度 {obs_levels} 档 {ext['distinct_levels']}")

    rng = random.Random(20261110)
    # (A) 完全随机：64 个密码子随机指派到 21 类
    counts = Counter()
    for _ in range(N_SIM):
        assign = {i: rng.choice(list(range(20)) + ["*"]) for i in range(64)}
        counts[levels(assign)] += 1
    p_ge = sum(v for k, v in counts.items() if k <= obs_levels) / N_SIM
    print(f"(A) 完全随机 64→21 类：档数分布 {dict(sorted(counts.items())[:8])} …")
    print(f"    档数 ≤ {obs_levels} 的比例：{p_ge:.4%}  ⟹ 实测 5 档"
          f"{'不太可能由随机产生' if p_ge < 0.05 else '在随机范围内'}")

    # (B) 随机化后的期望档数
    exp_levels = sum(k * v for k, v in counts.items()) / N_SIM
    print(f"(B) 随机表的平均档数：{exp_levels:.2f}；实测 {obs_levels}")

    out = {
        "observed_levels": obs_levels,
        "observed_distribution": dist,
        "null_full_random": {str(k): v / N_SIM for k, v in sorted(counts.items())},
        "p_levels_le_observed": p_ge,
        "null_mean_levels": exp_levels,
        "n_sim": N_SIM,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print("落盘", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
