#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Persistent open reservoir with rare closure exits
=================================================

This sandbox tests the proposal:

    open vibrations continue to evolve;
    only the small subset that closes exits the active layer.

There is no finite layer lifetime and no reseeding in this model. A path is
an integer vector q in Z^d. Each step changes one coordinate by +1 or -1.
The path closes when q becomes the zero vector.

All successors are instantiated. No probabilities or mutation rates are used.

Run:
    python3 simulations/zero_sum_open_reservoir.py
"""
from __future__ import annotations

import json
import math
import os
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    os.path.join(tempfile.gettempdir(), "mpl-zero-sum"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "simulations" / "zero_sum_open_reservoir_results.json"
OUT_PNG = ROOT / "simulations" / "zero_sum_open_reservoir.png"

DIMENSIONS = (1, 2, 4)
EPOCHS = 10


@dataclass(frozen=True)
class Scenario:
    dimension: int

    @property
    def name(self) -> str:
        return f"d{self.dimension}"

    @property
    def label(self) -> str:
        return f"dimension {self.dimension}"


SCENARIOS = tuple(Scenario(dimension) for dimension in DIMENSIONS)


def initial_states(dimension: int) -> tuple[tuple[int, ...], ...]:
    plus = [0] * dimension
    minus = [0] * dimension
    plus[0] = 1
    minus[0] = -1
    return tuple(plus), tuple(minus)


def is_closed(state: tuple[int, ...]) -> bool:
    return all(value == 0 for value in state)


def successors(state: tuple[int, ...]) -> list[tuple[int, ...]]:
    result = []
    for index in range(len(state)):
        for delta in (1, -1):
            child = list(state)
            child[index] += delta
            result.append(tuple(child))
    return result


def total_balance(active: Counter[tuple[int, ...]]) -> int:
    return sum(
        count * sum(state)
        for state, count in active.items()
    )


def simulate(scenario: Scenario) -> dict[str, object]:
    seeds = initial_states(scenario.dimension)
    active: Counter[tuple[int, ...]] = Counter({seed: 1 for seed in seeds})
    cumulative_closures = 0
    cumulative_transitions = 0

    snapshots: list[dict[str, object]] = [
        {
            "epoch": 0,
            "active_paths": sum(active.values()),
            "distinct_active_states": len(active),
            "closures_this_epoch": 0,
            "closure_fraction": 0.0,
            "cumulative_closures": 0,
            "cumulative_closure_fraction": 0.0,
            "total_balance": total_balance(active),
        }
    ]

    for epoch in range(EPOCHS):
        next_active: Counter[tuple[int, ...]] = Counter()
        closures = 0
        for state, count in active.items():
            for child in successors(state):
                if is_closed(child):
                    closures += count
                else:
                    next_active[child] += count

        cumulative_closures += closures
        open_successors = sum(next_active.values())
        total_successors = open_successors + closures
        cumulative_transitions += total_successors
        closure_fraction = (
            closures / total_successors
            if total_successors
            else 0.0
        )
        active = next_active
        snapshots.append(
            {
                "epoch": epoch + 1,
                "active_paths": open_successors,
                "distinct_active_states": len(active),
                "closures_this_epoch": closures,
                "closure_fraction": closure_fraction,
                "cumulative_closures": cumulative_closures,
                "cumulative_closure_fraction": (
                    cumulative_closures / cumulative_transitions
                    if cumulative_transitions
                    else 0.0
                ),
                "total_balance": total_balance(active),
            }
        )

    final = snapshots[-1]
    return {
        "scenario": scenario.name,
        "label": scenario.label,
        "dimension": scenario.dimension,
        "epochs": EPOCHS,
        "snapshots": snapshots,
        "final_active_paths": final["active_paths"],
        "final_closure_fraction": final["closure_fraction"],
        "final_cumulative_closure_fraction": final[
            "cumulative_closure_fraction"
        ],
        "final_cumulative_closures": final["cumulative_closures"],
        "final_total_balance": final["total_balance"],
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
    colors = ["#2B7A78", "#C77C28", "#B44A5A", "#3A6EA5"]
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8))

    for run, color in zip(runs, colors):
        epochs = [snapshot["epoch"] for snapshot in run["snapshots"]]
        axes[0].plot(
            epochs,
            [snapshot["active_paths"] for snapshot in run["snapshots"]],
            color=color,
            marker="o",
            markersize=3,
            linewidth=1.8,
            label=run["label"],
        )
        axes[1].plot(
            epochs,
            [
                snapshot["cumulative_closure_fraction"]
                for snapshot in run["snapshots"]
            ],
            color=color,
            marker="o",
            markersize=3,
            linewidth=1.8,
            label=run["label"],
        )
        axes[2].plot(
            epochs,
            [
                snapshot["cumulative_closures"]
                for snapshot in run["snapshots"]
            ],
            color=color,
            marker="o",
            markersize=3,
            linewidth=1.8,
            label=run["label"],
        )

    axes[0].set_yscale("log")
    axes[0].set_title("Active open paths")
    axes[0].set_xlabel("epoch")
    axes[0].set_ylabel("count (log scale)")
    axes[0].grid(alpha=0.22)
    axes[0].legend(frameon=False)

    axes[1].set_title("Cumulative closure fraction")
    axes[1].set_xlabel("epoch")
    axes[1].set_ylabel("closures / all successors")
    axes[1].grid(alpha=0.22)
    axes[1].legend(frameon=False)

    axes[2].set_yscale("log")
    axes[2].set_title("Cumulative closures")
    axes[2].set_xlabel("epoch")
    axes[2].set_ylabel("count (log scale)")
    axes[2].grid(alpha=0.22)
    axes[2].legend(frameon=False)

    fig.suptitle("Persistent open reservoir with closure exits", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT_PNG, dpi=180)
    plt.close(fig)


def main() -> None:
    runs = [simulate(scenario) for scenario in SCENARIOS]
    make_plot(runs)

    verification = {
        "active_layer_never_empties": all(
            snapshot["active_paths"] > 0
            for run in runs
            for snapshot in run["snapshots"]
        ),
        "open_paths_grow_without_reseeding": all(
            run["final_active_paths"] > run["snapshots"][0]["active_paths"]
            for run in runs
        ),
        "closure_fraction_is_less_than_one": all(
            snapshot["cumulative_closure_fraction"] < 1.0
            for run in runs
            for snapshot in run["snapshots"]
        ),
        "global_balance_stays_zero": all(
            snapshot["total_balance"] == 0
            for run in runs
            for snapshot in run["snapshots"]
        ),
        "larger_dimension_lowers_cumulative_closure_fraction": (
            runs[-1]["final_cumulative_closure_fraction"]
            < runs[0]["final_cumulative_closure_fraction"]
        ),
    }

    payload = {
        "model": {
            "state_space": "integer vectors in Z^d",
            "closure_rule": "state becomes the zero vector",
            "step_alphabet": "one coordinate changes by +1 or -1",
            "layer_lifetime": None,
            "reseeding": None,
            "probabilities": "none",
            "mutations": "none",
            "dimensions": list(DIMENSIONS),
            "epochs": EPOCHS,
        },
        "runs": runs,
        "verification": verification,
    }
    OUT_JSON.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print("persistent open reservoir experiment")
    for key, value in verification.items():
        print(f"  [{'v' if value else 'x'}] {key}")
    for run in runs:
        print(
            f"{run['scenario']:>4}  "
            f"active={run['final_active_paths']:>12}  "
            f"cumulative_closure_fraction="
            f"{run['final_cumulative_closure_fraction']:.6f}  "
            f"cumulative_closures={run['final_cumulative_closures']:>12}  "
            f"balance={run['final_total_balance']}"
        )
    print(f"\nplot: {OUT_PNG}")
    print(f"json: {OUT_JSON}")


if __name__ == "__main__":
    main()
