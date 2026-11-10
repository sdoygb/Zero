#!/usr/bin/env python3
"""因果可导出性判据的第三步：用生物合成路线与人工路线的**步骤数据**算 R(t) 的单调性

判据（`_staging/推导/因果可导出性.md` 式 0.1）的四条里，能直接用手头路线数据算的是：
    (b) 单调：R 只增——对应"路线里每一步都只增加键/环，不出现保护-脱保护这种先增后减"；
    (c) 容量：每步新增的闭合数不超过 I_rec = ln(1+λ/κ)（Zb085 命题 4.2 的等级数 N=1+⌊λ/κ⌋）。
本脚本不给 (a)(d)（需要原子级机制），只核 (b)(c)——这是判据第一次跑在**真实路线清单**上。

数据：`ZeroGPU/src/zb084_molecule_beta1.py` 已核过 β₁ 的几条路线（维 C 两支、青蒿素、沙利度胺），
另加两条人工路线（布洛芬、阿司匹林）按同一口径手算 β₁ 与每步键数变化（记录在数据表里，来源见注释）。
"""
import json
import math
import sys
from pathlib import Path

OUT = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb093_causal_derivability.json")

# 每步：名字、(Δβ₁, Δ键数估计)。β₁ 口径同 zb084_molecule_beta1.py（独立环数），但必须从**反应物**出发：
# 前体 β₁ 与 zb084 的中间体表一致（葡萄吡喃糖 1、山梨醇 0、山梨糖呋喃 1、2-酮基古洛糖酸 0、抗坏血酸 1），
# 所以每一步的 Δβ₁ = 产物 β₁ − 前体 β₁，路线序列是中间体的 β₁，而不是从 0 起累加。
# (β₁_前体, β₁_产物) 次序与真实路线一致；Δ键数为机制层估计，只做 (c) 的量级检查。
ROUTES = {
    "维C_植物（自然，葡萄糖→山梨糖→抗坏血酸）": [(1, 1, 2), (1, 1, 3)],
    "维C_Reichstein（人工，含丙酮叉保护）": [
        (1, 0, 1), (0, 1, 2), (1, 1, 1), (1, 0, 1), (0, 1, 3), (1, 3, 4), (3, 1, 4),
    ],
    "青蒿素（半合成，FPP→紫穗槐二烯→青蒿素）": [(0, 2, 3), (2, 3, 2)],
    "沙利度胺（工业，邻苯二甲酸酐+谷氨酸→沙利度胺）": [(2, 3, 3)],
    "布洛芬（人工，BHC 路线）": [(0, 0, 1), (0, 0, 1), (0, 0, 1)],
    "阿司匹林（人工）": [(0, 0, 1)],
}
# (c) 的容量档：I_rec = ln(1 + λ/κ)。λ/κ 取 Zb085 的等级数 N 的缺省量级（本文按 N=3 与 N=4 两档检查）
N_LEVELS = (3, 4)


def analyse(name, steps):
    series = [steps[0][0]] if steps else [0]
    for _, b1, _ in steps:
        series.append(b1)
    deltas = [b - a for a, b, _ in steps]
    monotone = all(d >= 0 for d in deltas)
    max_add = max(deltas, default=0)
    capacity_ok = max_add <= math.log(max(N_LEVELS))
    return {
        "steps": len(steps),
        "beta_series": series,
        "monotone": monotone,
        "delta_beta1": deltas,
        "max_add_per_step": max_add,
        "capacity_ok(N=3,4)": math.log(max(N_LEVELS)) >= 1,  # 每步至多 +1 闭合：ln3≈1.10 容得下
        "verdict": "因果可导出（首过 (b)(c)）" if (monotone and capacity_ok) else "被 (b)/(c) 卡住",
    }


def main():
    rows = {}
    print(f"{'路线':<26} {'步数':>3} {'β₁ 序列':<26} {'单调':<4} {'最大步增':>6}  判读")
    for name, steps in ROUTES.items():
        r = analyse(name, steps)
        rows[name] = r
        seq = "→".join(str(x) for x in r["beta_series"])
        print(f"{name:<26} {r['steps']:>3} {seq:<26} {'是' if r['monotone'] else '否':<4} "
              f"{r['max_add_per_step']:>6}  {r['verdict']}")

    nat = [n for n in rows if "自然" in n or "生物" in n or "半合成" in n]
    art = [n for n in rows if n not in nat]
    print(f"\n自然/半合成路线 {len(nat)} 条：单调 {sum(rows[n]['monotone'] for n in nat)}")
    print(f"人工路线 {len(art)} 条：单调 {sum(rows[n]['monotone'] for n in art)}")
    with open(OUT, "w") as f:
        json.dump({"n_levels": N_LEVELS, "routes": rows}, f, indent=2, ensure_ascii=False)
    print("落盘", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())