#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Zero-sum three-layer universe
=============================

This sandbox implements the R/P/E/D design:

    R: result layer
       closed histories that have become stable records

    P: history layer
       exact closed histories, including how they were assembled

    E: evolution layer
       finite-lifetime layers containing open and closed continuations

    D: dead paths
       paths that disappear when a layer is destroyed

No probabilities, mutation rates, or fertility weights are used. Every
zero-sum continuation is generated.

A layer is seeded by one closed history w. During its lifetime L it generates
all extensions

    e = one balanced history of even length tau <= L.

The resulting history is w + e. It is recorded in P, its cyclic result is
recorded in R, and it seeds a new layer. When the original layer reaches age
L, all remaining paths in that layer go to D.

The initial seed is the empty history. The empty seed simply means a layer
starting from the zero state.

Run:
    python3 simulations/zero_sum_tri_layer_universe.py
"""
from __future__ import annotations

import json
import math
import os
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    os.path.join(tempfile.gettempdir(), "mpl-zero-sum"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    from simulations.zero_sum_cycle_evolution import canonical_cycle
except ModuleNotFoundError:
    from zero_sum_cycle_evolution import canonical_cycle


ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "simulations" / "zero_sum_tri_layer_universe_results.json"
OUT_PNG = ROOT / "simulations" / "zero_sum_tri_layer_universe.png"

STEPS = 10
SCENARIOS = (0, 4, 6)
COUNT_SCALE_LIMIT = 1.0e12


@dataclass(frozen=True)
class Scenario:
    lifetime: int
    label: str


SCENARIO_DATA = tuple(
    Scenario(
        lifetime=lifetime,
        label=(
            "no evolution layer"
            if lifetime == 0
            else f"finite layer L={lifetime}"
        ),
    )
    for lifetime in SCENARIOS
)


def balanced_words(period: int) -> tuple[tuple[int, ...], ...]:
    words = []
    positive = period // 2
    for positions in combinations(range(period), positive):
        positive_positions = set(positions)
        words.append(
            tuple(
                1 if index in positive_positions else -1
                for index in range(period)
            )
        )
    return tuple(words)


def word_to_string(word: tuple[int, ...]) -> str:
    if not word:
        return "0"
    return "".join("+" if value == 1 else "-" for value in word)


def shannon(values) -> float:
    counts = np.asarray(list(values), dtype=float)
    total = float(np.sum(counts))
    if total <= 0.0:
        return 0.0
    probabilities = counts[counts > 0.0] / total
    result = float(-np.sum(probabilities * np.log(probabilities)))
    return 0.0 if abs(result) < 1.0e-15 else result


def top_entries(counter: Counter, limit: int = 12) -> list[dict[str, object]]:
    total = float(sum(counter.values()))
    if total <= 0.0:
        return []
    result = []
    for key, value in counter.most_common(limit):
        result.append(
            {
                "word": word_to_string(key),
                "length": len(key),
                "count": float(value),
                "fraction": float(value / total),
            }
        )
    return result


def rescale_state(
    future: list[Counter],
    layer_deaths: list[float],
    path_deaths: list[float],
    future_closure_time: list[float],
    result: Counter,
    history: Counter,
    scale: float,
) -> tuple[
    list[Counter],
    list[float],
    list[float],
    list[float],
    Counter,
    Counter,
]:
    scaled_future = [
        Counter({key: value * scale for key, value in bucket.items()})
        for bucket in future
    ]
    scaled_deaths = [value * scale for value in layer_deaths]
    scaled_path_deaths = [value * scale for value in path_deaths]
    scaled_closure_times = [value * scale for value in future_closure_time]
    return (
        scaled_future,
        scaled_deaths,
        scaled_path_deaths,
        scaled_closure_times,
        Counter({key: value * scale for key, value in result.items()}),
        Counter({key: value * scale for key, value in history.items()}),
    )


def simulate(
    scenario: Scenario,
    extensions: dict[int, tuple[tuple[int, ...], ...]],
) -> dict[str, object]:
    lifetime = scenario.lifetime
    horizon = STEPS + lifetime + 1
    future: list[Counter] = [Counter() for _ in range(horizon + 1)]
    layer_deaths = [0.0 for _ in range(horizon + 1)]
    path_deaths = [0.0 for _ in range(horizon + 1)]
    future_closure_time = [0.0 for _ in range(horizon + 1)]
    future[0][()] = 1.0

    result: Counter = Counter()
    history: Counter = Counter()
    active_layers = 0.0
    total_layer_deaths = 0.0
    total_path_deaths = 0.0
    total_closures = 0.0
    closure_time_mass = 0.0

    snapshots = []
    rescale_steps = 0

    for time in range(STEPS + 1):
        born = future[time]
        born_count = float(sum(born.values()))
        active_layers += born_count - layer_deaths[time]
        total_layer_deaths += layer_deaths[time]
        total_path_deaths += path_deaths[time]
        closure_time_mass += future_closure_time[time]

        for closed_history, layer_count in born.items():
            if not closed_history:
                continue
            result[canonical_cycle(closed_history)] += layer_count
            history[closed_history] += layer_count
            total_closures += layer_count

        total_counts = (
            sum(float(sum(bucket.values())) for bucket in future)
            + float(sum(result.values()))
            + float(sum(history.values()))
            + total_path_deaths
        )
        if total_counts > COUNT_SCALE_LIMIT:
            scale = COUNT_SCALE_LIMIT / total_counts
            (
                future,
                layer_deaths,
                path_deaths,
                result,
                history,
            ) = rescale_state(
                future,
                layer_deaths,
                path_deaths,
                future_closure_time,
                result,
                history,
                scale,
            )
            active_layers *= scale
            total_layer_deaths *= scale
            total_path_deaths *= scale
            total_closures *= scale
            closure_time_mass *= scale
            born = future[time]
            rescale_steps += 1

        snapshots.append(
            {
                "time": time,
                "active_layers": active_layers,
                "result_records": float(sum(result.values())),
                "result_modes": int(len(result)),
                "result_diversity": shannon(result.values()),
                "history_records": float(sum(history.values())),
                "history_paths": int(len(history)),
                "history_diversity": shannon(history.values()),
                "dead_paths": total_path_deaths,
                "closed_events": total_closures,
                "mean_closure_time": (
                    closure_time_mass / total_closures
                    if total_closures > 0.0
                    else 0.0
                ),
            }
        )

        if time == STEPS:
            break

        if lifetime == 0:
            continue

        for seed, layer_count in list(born.items()):
            if layer_count <= 0.0:
                continue

            for closure_time, extension_words in extensions.items():
                if closure_time > lifetime:
                    continue
                target_time = time + closure_time
                if target_time > STEPS:
                    continue
                for extension in extension_words:
                    closed_history = seed + extension
                    future[target_time][closed_history] += layer_count
                    future_closure_time[target_time] += layer_count * closure_time

            layer_deaths[time + lifetime] += layer_count
            path_deaths[time + lifetime] += layer_count * (2**lifetime)

    final = snapshots[-1]
    return {
        "scenario": scenario.label,
        "lifetime": lifetime,
        "snapshots": snapshots,
        "final_active_layers": final["active_layers"],
        "final_result_records": final["result_records"],
        "final_result_modes": final["result_modes"],
        "final_history_records": final["history_records"],
        "final_history_paths": final["history_paths"],
        "final_result_diversity": final["result_diversity"],
        "final_history_diversity": final["history_diversity"],
        "final_dead_paths": final["dead_paths"],
        "final_mean_closure_time": final["mean_closure_time"],
        "top_results": top_entries(result),
        "top_histories": top_entries(history),
        "rescale_steps": rescale_steps,
    }


def make_plot(runs: list[dict[str, object]]) -> None:
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "legend.fontsize": 8,
            "axes.unicode_minus": False,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.2))
    colors = {
        0: "#777777",
        4: "#C77C28",
        6: "#2B7A78",
    }

    for run in runs:
        times = [snapshot["time"] for snapshot in run["snapshots"]]
        axes[0, 0].plot(
            times,
            [max(snapshot["result_records"], 0.5) for snapshot in run["snapshots"]],
            linewidth=1.8,
            color=colors[run["lifetime"]],
            label=run["scenario"],
        )
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_title("Result layer R")
    axes[0, 0].set_xlabel("global time")
    axes[0, 0].set_ylabel("closed records")
    axes[0, 0].grid(alpha=0.22)
    axes[0, 0].legend(frameon=False)

    for run in runs:
        times = [snapshot["time"] for snapshot in run["snapshots"]]
        axes[0, 1].plot(
            times,
            [max(snapshot["history_records"], 0.5) for snapshot in run["snapshots"]],
            linewidth=1.8,
            color=colors[run["lifetime"]],
            label=run["scenario"],
        )
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_title("History layer P")
    axes[0, 1].set_xlabel("global time")
    axes[0, 1].set_ylabel("recorded histories")
    axes[0, 1].grid(alpha=0.22)
    axes[0, 1].legend(frameon=False)

    for run in runs:
        times = [snapshot["time"] for snapshot in run["snapshots"]]
        axes[1, 0].plot(
            times,
            [snapshot["active_layers"] for snapshot in run["snapshots"]],
            linewidth=1.8,
            color=colors[run["lifetime"]],
            label=run["scenario"],
        )
    axes[1, 0].set_yscale("symlog", linthresh=1.0)
    axes[1, 0].set_title("Evolution layers E")
    axes[1, 0].set_xlabel("global time")
    axes[1, 0].set_ylabel("active layers")
    axes[1, 0].grid(alpha=0.22)
    axes[1, 0].legend(frameon=False)

    for run in runs:
        times = [snapshot["time"] for snapshot in run["snapshots"]]
        axes[1, 1].plot(
            times,
            [max(snapshot["dead_paths"], 0.5) for snapshot in run["snapshots"]],
            linewidth=1.8,
            color=colors[run["lifetime"]],
            label=run["scenario"],
        )
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_title("Destroyed open paths D")
    axes[1, 1].set_xlabel("global time")
    axes[1, 1].set_ylabel("dead path records")
    axes[1, 1].grid(alpha=0.22)
    axes[1, 1].legend(frameon=False)

    fig.suptitle("Zero-sum three-layer universe", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUT_PNG, dpi=180)
    plt.close(fig)


def main() -> None:
    max_closure_time = max(SCENARIOS)
    extensions = {
        period: balanced_words(period)
        for period in range(2, max_closure_time + 1, 2)
    }
    runs = [
        simulate(scenario, extensions)
        for scenario in SCENARIO_DATA
    ]
    make_plot(runs)

    payload = {
        "model": {
            "R": "cyclic canonical results of exact closed histories",
            "P": "exact closed histories",
            "E": "finite-lifetime layers",
            "D": "paths destroyed when a layer reaches its lifetime",
            "initial_seed": "empty zero history",
            "probabilities": "none",
            "mutations": "none",
            "scenarios": {
                scenario.label: scenario.lifetime
                for scenario in SCENARIO_DATA
            },
            "steps": STEPS,
            "count_scale_limit": COUNT_SCALE_LIMIT,
            "note": (
                "All zero-sum continuations are generated. Scaling is "
                "numerical only and does not change proportions."
            ),
        },
        "runs": runs,
    }
    OUT_JSON.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print("three-layer universe summary")
    for run in runs:
        print(
            f"L={run['lifetime']:>2}  "
            f"R={run['final_result_records']:.1f}  "
            f"Rmodes={run['final_result_modes']}  "
            f"P={run['final_history_records']:.1f}  "
            f"Ppaths={run['final_history_paths']}  "
            f"E={run['final_active_layers']:.1f}  "
            f"D={run['final_dead_paths']:.1f}  "
            f"meanT={run['final_mean_closure_time']:.3f}"
        )
    print(f"\nplot: {OUT_PNG}")
    print(f"json: {OUT_JSON}")


if __name__ == "__main__":
    main()
