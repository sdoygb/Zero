#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模平衡 · 探索模拟：零和模式的生殖选择
================================================================================
问题：
  不是“谁存活得久”，而是“谁能产生更多后代谱系”。

模型：
  取五种零和模式。每种模式有：
    cost        单次循环作用量成本
    return_time 回零所需时间
    fidelity    复制规整度 / 保真度
    base_rate   整体生殖标度

  每个模式每一代产生 Poisson 个后代：

      birth_rate = base_rate * fidelity
                   * min_cost / cost
                   * min_return_time / return_time

      N_i(t+1) = N_i(t) + Poisson(birth_rate_i * N_i(t)).

  所有模式的荷向量都满足 q_i 的和为零，因此总荷为零严格保持。

  对照模式 long_lived_sterile 永远不死亡，但生殖率为零。它说明：
  “存活”与“留下后代”是两件不同的事。

这不是任何外部公理的推论，只是恢复层候选机制的数值试验。

运行：
  python3 simulations/zero_generative_selection.py
  python3 simulations/zero_generative_selection.py --plot simulations/zero_generative_selection.png
  python3 simulations/zero_generative_selection.py --json simulations/zero_generative_selection_results.json
"""
from __future__ import annotations

import argparse
import json
import os
import statistics
import sys
from dataclasses import asdict, dataclass

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import numpy as np


@dataclass(frozen=True)
class Mode:
    name: str
    cost: float
    return_time: float
    fidelity: float
    base_rate: float
    charge: tuple[int, ...]
    note: str


@dataclass
class ModeSummary:
    name: str
    note: str
    cost: float
    return_time: float
    fidelity: float
    birth_rate: float
    mean_final_population: float
    mean_descendants: float
    final_share: float
    rank: int


MODES = (
    Mode(
        name="cheap_fast_regular",
        cost=1.0,
        return_time=1.0,
        fidelity=1.0,
        base_rate=0.36,
        charge=(1, -1, 0, 0, 0),
        note="省、回零快、复制规整",
    ),
    Mode(
        name="cheap_fast_irregular",
        cost=1.0,
        return_time=1.0,
        fidelity=0.58,
        base_rate=0.36,
        charge=(1, 0, -1, 0, 0),
        note="省、回零快，但复制不规整",
    ),
    Mode(
        name="expensive_fast_regular",
        cost=4.0,
        return_time=1.0,
        fidelity=1.0,
        base_rate=0.36,
        charge=(0, 1, 0, -1, 0),
        note="回零快且规整，但每次生殖成本高",
    ),
    Mode(
        name="cheap_slow_regular",
        cost=1.0,
        return_time=4.0,
        fidelity=1.0,
        base_rate=0.36,
        charge=(0, 0, 1, 0, -1),
        note="省且规整，但回零慢",
    ),
    Mode(
        name="long_lived_sterile",
        cost=1.0,
        return_time=1.0,
        fidelity=1.0,
        base_rate=0.0,
        charge=(1, 1, -1, -1, 0),
        note="永不死亡，但生殖率为零",
    ),
)


def reproduction_rate(mode: Mode, min_cost: float, min_return_time: float) -> float:
    """保零生殖核的期望后代速率。"""
    return (
        mode.base_rate
        * mode.fidelity
        * (min_cost / mode.cost)
        * (min_return_time / mode.return_time)
    )


def total_charge(populations: np.ndarray, modes: tuple[Mode, ...]) -> int:
    charges = np.array([mode.charge for mode in modes], dtype=int)
    return int(np.sum(populations @ charges))


def simulate_once(
    modes: tuple[Mode, ...],
    rates: np.ndarray,
    generations: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    n_modes = len(modes)
    populations = np.ones(n_modes, dtype=np.int64)
    history = [populations.copy()]

    for _ in range(generations):
        births = rng.poisson(rates * populations)
        populations = populations + births.astype(np.int64)
        history.append(populations.copy())

        if np.sum(populations) > 50_000_000:
            raise RuntimeError("population exceeded simulation safety cap")

    return np.array(history, dtype=np.int64), np.array(
        [total_charge(row, modes) for row in history],
        dtype=np.int64,
    )


def summarize(
    modes: tuple[Mode, ...],
    rates: np.ndarray,
    histories: np.ndarray,
) -> list[ModeSummary]:
    final_populations = histories[:, -1, :]
    mean_final = final_populations.mean(axis=0)
    mean_descendants = mean_final - 1.0
    total_mean = float(np.sum(mean_final))
    order = np.argsort(-mean_final)
    rank = np.empty(len(modes), dtype=int)
    for position, index in enumerate(order, start=1):
        rank[index] = position

    summaries: list[ModeSummary] = []
    for index, mode in enumerate(modes):
        summaries.append(
            ModeSummary(
                name=mode.name,
                note=mode.note,
                cost=mode.cost,
                return_time=mode.return_time,
                fidelity=mode.fidelity,
                birth_rate=float(rates[index]),
                mean_final_population=float(mean_final[index]),
                mean_descendants=float(mean_descendants[index]),
                final_share=float(mean_final[index] / total_mean),
                rank=int(rank[index]),
            )
        )
    return summaries


def run_experiment(
    runs: int,
    generations: int,
    seed: int,
) -> tuple[list[ModeSummary], np.ndarray, np.ndarray, dict[str, object]]:
    min_cost = min(mode.cost for mode in MODES)
    min_return_time = min(mode.return_time for mode in MODES)
    rates = np.array(
        [
            reproduction_rate(mode, min_cost, min_return_time)
            for mode in MODES
        ],
        dtype=float,
    )

    all_histories = []
    max_charge_error = 0
    for run_index in range(runs):
        history, charge_history = simulate_once(
            MODES,
            rates,
            generations,
            seed + 104729 * run_index,
        )
        all_histories.append(history)
        max_charge_error = max(
            max_charge_error,
            int(np.max(np.abs(charge_history))),
        )

    histories = np.array(all_histories, dtype=np.int64)
    summaries = summarize(MODES, rates, histories)
    exemplar = histories[0].astype(float)
    diagnostics = {
        "modes": len(MODES),
        "generations": generations,
        "runs": runs,
        "seed": seed,
        "max_total_charge_error": max_charge_error,
        "population_safety_cap": 50_000_000,
    }
    return summaries, exemplar, rates, diagnostics


def print_report(
    summaries: list[ModeSummary],
    diagnostics: dict[str, object],
) -> None:
    print("零和模式的生殖选择模拟")
    print("=" * 108)
    print(
        f"runs={diagnostics['runs']}  "
        f"generations={diagnostics['generations']}  "
        f"max total-charge error={diagnostics['max_total_charge_error']}"
    )
    print("-" * 108)
    print(
        f"{'rank':>4} {'mode':<26} {'cost':>6} {'return':>7} "
        f"{'fidelity':>9} {'birth':>8} {'desc. mean':>14} {'share':>9}"
    )
    for row in sorted(summaries, key=lambda item: item.rank):
        print(
            f"{row.rank:>4} "
            f"{row.name:<26} "
            f"{row.cost:>6.2f} "
            f"{row.return_time:>7.2f} "
            f"{row.fidelity:>9.2f} "
            f"{row.birth_rate:>8.4f} "
            f"{row.mean_descendants:>14.2f} "
            f"{row.final_share:>8.2%}"
        )
    print("-" * 108)
    print("关键对照：long_lived_sterile 的存活数始终为 1，但后代始终为 0。")
    print("因此本图的排序来自生殖，不来自寿命或存活率。")


def save_json(
    path: str,
    summaries: list[ModeSummary],
    diagnostics: dict[str, object],
) -> None:
    payload = {
        "diagnostics": diagnostics,
        "summaries": [asdict(row) for row in summaries],
        "model": {
            "birth_rate": (
                "base_rate * fidelity * min_cost/cost "
                "* min_return_time/return_time"
            ),
            "population_update": "N(t+1)=N(t)+Poisson(birth_rate*N(t))",
            "zero_constraint": "all mode charge vectors sum to zero",
        },
    }
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def save_plot(
    path: str,
    summaries: list[ModeSummary],
    exemplar: np.ndarray,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.6), constrained_layout=True)
    generations = np.arange(exemplar.shape[0])

    for index, mode in enumerate(MODES):
        axes[0].plot(
            generations,
            np.maximum(exemplar[:, index], 1.0),
            linewidth=2.0,
            label=mode.name,
        )
    axes[0].set_yscale("log")
    axes[0].set_xlabel("generation")
    axes[0].set_ylabel("population / descendants")
    axes[0].set_title("single-run descendant growth")
    axes[0].grid(alpha=0.25)
    axes[0].legend(fontsize=8)

    ordered = sorted(summaries, key=lambda row: row.mean_final_population)
    names = [row.name.replace("_", "\n") for row in ordered]
    values = [row.mean_final_population for row in ordered]
    colors = [
        "#c23b22" if row.name == "long_lived_sterile" else "#2f6f9f"
        for row in ordered
    ]
    axes[1].barh(names, values, color=colors)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("mean final population (log scale)")
    axes[1].set_title("reproductive selection")
    axes[1].grid(axis="x", alpha=0.25)

    fig.suptitle(
        "Zero-sum generative selection: descendants, not mere persistence",
        fontsize=14,
    )
    fig.savefig(path, dpi=180)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Simulate zero-sum mode reproduction and lineage selection."
    )
    parser.add_argument("--runs", type=int, default=80)
    parser.add_argument("--generations", type=int, default=34)
    parser.add_argument("--seed", type=int, default=210)
    parser.add_argument("--json", type=str, default="")
    parser.add_argument("--plot", type=str, default="")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.runs <= 0 or args.generations <= 0:
        print("runs and generations must be positive", file=sys.stderr)
        return 2

    summaries, exemplar, _rates, diagnostics = run_experiment(
        runs=args.runs,
        generations=args.generations,
        seed=args.seed,
    )
    print_report(summaries, diagnostics)

    if args.json:
        save_json(args.json, summaries, diagnostics)
        print(f"\nJSON 已写入：{args.json}")
    if args.plot:
        save_plot(args.plot, summaries, exemplar)
        print(f"图像已写入：{args.plot}")

    winner = max(summaries, key=lambda row: row.mean_final_population)
    sterile = next(
        row for row in summaries if row.name == "long_lived_sterile"
    )
    if winner.name != "cheap_fast_regular":
        print(
            "预期最省、最快、最规整的模式应排名第一。",
            file=sys.stderr,
        )
        return 1
    if sterile.mean_final_population != 1.0 or sterile.mean_descendants != 0.0:
        print("对照模式 long_lived_sterile 不应产生后代。", file=sys.stderr)
        return 1

    print("\n判决：生殖选择排名与“省、回零快、规整”一致。")
    print("对照模式长期存在但没有后代，说明存活本身不构成谱系扩张。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
