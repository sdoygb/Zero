#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Closure exit from the evolution layer
=====================================

This sandbox tests the proposed rule:

    once a path closes, it exits the evolution layer;
    only its exact history and closed zero mode remain.

An active path is a word of plus/minus steps. It closes when its cumulative
balance returns to zero. We compare two rules:

    continue_after_closure:
        closed paths remain active and can continue.

    exit_after_closure:
        closed paths are removed from the active layer, written to the local
        history layer P_i, and their cyclic class is written to the global
        closure-zero layer Z.

Open paths that exceed the finite layer lifetime go to D_i.

No probabilities or mutation rates are used. Every allowed next step is
instantiated.

Run:
    python3 simulations/zero_sum_closure_exit.py
"""
from __future__ import annotations

import json
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

try:
    from simulations.zero_sum_cycle_evolution import canonical_cycle
except ModuleNotFoundError:
    from zero_sum_cycle_evolution import canonical_cycle


ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "simulations" / "zero_sum_closure_exit_results.json"
OUT_PNG = ROOT / "simulations" / "zero_sum_closure_exit.png"

LAYER_LIFETIME = 4
EPOCHS = 6
INITIAL_WORDS = (
    "+",
    "-",
    "++",
    "--",
    "+++",
    "---",
    "+-+",
    "-+-",
)


@dataclass(frozen=True)
class Scenario:
    name: str
    keep_closed_active: bool
    label: str


SCENARIOS = (
    Scenario(
        name="continue_after_closure",
        keep_closed_active=True,
        label="closed paths remain active",
    ),
    Scenario(
        name="exit_after_closure",
        keep_closed_active=False,
        label="closed paths exit E",
    ),
)


def word_from_string(text: str) -> tuple[int, ...]:
    return tuple(1 if char == "+" else -1 for char in text)


def word_to_string(word: tuple[int, ...]) -> str:
    if not word:
        return "0"
    return "".join("+" if value == 1 else "-" for value in word)


def is_closed(word: tuple[int, ...]) -> bool:
    return bool(word) and sum(word) == 0


def closed_class(word: tuple[int, ...]) -> tuple[int, ...]:
    return canonical_cycle(word)


def active_state(
    active: list[set[tuple[int, ...]]],
) -> dict[str, int]:
    all_paths = [word for paths in active for word in paths]
    return {
        "active_paths": len(all_paths),
        "active_closed_paths": sum(
            1 for word in all_paths if is_closed(word)
        ),
        "active_open_paths": sum(
            1 for word in all_paths if not is_closed(word)
        ),
    }


def simulate(scenario: Scenario) -> dict[str, object]:
    active = [
        {word_from_string(text)}
        for text in INITIAL_WORDS
    ]
    histories = [set() for _ in INITIAL_WORDS]
    dead = [set() for _ in INITIAL_WORDS]
    zero_layer: Counter[tuple[int, ...]] = Counter()

    snapshots: list[dict[str, object]] = []
    for epoch in range(EPOCHS + 1):
        state = active_state(active)
        snapshots.append(
            {
                "epoch": epoch,
                **state,
                "history_records": sum(len(book) for book in histories),
                "zero_events": sum(zero_layer.values()),
                "zero_modes": len(zero_layer),
                "dead_paths": sum(len(book) for book in dead),
                "per_site_active": [
                    len(paths) for paths in active
                ],
                "per_site_closed_active": [
                    sum(1 for word in paths if is_closed(word))
                    for paths in active
                ],
                "per_site_history_records": [
                    len(book) for book in histories
                ],
                "per_site_dead_paths": [
                    len(book) for book in dead
                ],
            }
        )

        if epoch == EPOCHS:
            break

        next_active = [set() for _ in INITIAL_WORDS]
        exited_this_epoch = 0
        died_this_epoch = 0

        for site, paths in enumerate(active):
            for word in paths:
                if len(word) >= LAYER_LIFETIME:
                    dead[site].add(word)
                    died_this_epoch += 1
                    continue

                for step in (1, -1):
                    child = word + (step,)
                    if is_closed(child):
                        histories[site].add(child)
                        zero_layer[closed_class(child)] += 1
                        exited_this_epoch += 1
                        if scenario.keep_closed_active:
                            next_active[site].add(child)
                    else:
                        next_active[site].add(child)

        active = next_active
        snapshots[-1]["exits_this_epoch"] = exited_this_epoch
        snapshots[-1]["dead_this_epoch"] = died_this_epoch

    final = snapshots[-1]
    return {
        "scenario": scenario.name,
        "label": scenario.label,
        "layer_lifetime": LAYER_LIFETIME,
        "epochs": EPOCHS,
        "snapshots": snapshots,
        "final_active_paths": final["active_paths"],
        "final_active_closed_paths": final["active_closed_paths"],
        "final_history_records": final["history_records"],
        "final_zero_events": final["zero_events"],
        "final_zero_modes": final["zero_modes"],
        "final_dead_paths": final["dead_paths"],
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
    colors = {
        "continue_after_closure": "#C77C28",
        "exit_after_closure": "#2B7A78",
    }
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8))

    for run in runs:
        epochs = [snapshot["epoch"] for snapshot in run["snapshots"]]
        axes[0].plot(
            epochs,
            [snapshot["active_paths"] for snapshot in run["snapshots"]],
            color=colors[run["scenario"]],
            marker="o",
            markersize=3,
            linewidth=1.8,
            label=run["label"],
        )
        axes[0].plot(
            epochs,
            [
                snapshot["active_closed_paths"]
                for snapshot in run["snapshots"]
            ],
            color=colors[run["scenario"]],
            linestyle="--",
            marker="s",
            markersize=2.5,
            linewidth=1.3,
        )
    axes[0].set_title("Active paths (solid) and closed active paths (dashed)")
    axes[0].set_xlabel("epoch")
    axes[0].set_ylabel("count")
    axes[0].grid(alpha=0.22)
    axes[0].legend(frameon=False)

    for run in runs:
        epochs = [snapshot["epoch"] for snapshot in run["snapshots"]]
        axes[1].plot(
            epochs,
            [snapshot["history_records"] for snapshot in run["snapshots"]],
            color=colors[run["scenario"]],
            marker="o",
            markersize=3,
            linewidth=1.8,
            label=run["label"],
        )
        axes[1].plot(
            epochs,
            [snapshot["zero_modes"] for snapshot in run["snapshots"]],
            color=colors[run["scenario"]],
            linestyle="--",
            marker="s",
            markersize=2.5,
            linewidth=1.3,
        )
    axes[1].set_title("History records (solid) and zero modes (dashed)")
    axes[1].set_xlabel("epoch")
    axes[1].set_ylabel("count")
    axes[1].grid(alpha=0.22)
    axes[1].legend(frameon=False)

    labels = []
    zero_events = []
    dead_paths = []
    for run in runs:
        labels.append(run["scenario"].replace("_", "\n"))
        zero_events.append(run["final_zero_events"])
        dead_paths.append(run["final_dead_paths"])
    x = np.arange(len(labels))
    width = 0.38
    axes[2].bar(
        x - width / 2,
        zero_events,
        width=width,
        color="#3A6EA5",
        label="zero events",
    )
    axes[2].bar(
        x + width / 2,
        dead_paths,
        width=width,
        color="#B44A5A",
        label="dead paths",
    )
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(labels, fontsize=7)
    axes[2].set_title("Final closure and death counts")
    axes[2].set_ylabel("count")
    axes[2].grid(axis="y", alpha=0.22)
    axes[2].legend(frameon=False)

    fig.suptitle("Closure exit from the evolution layer", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT_PNG, dpi=180)
    plt.close(fig)


def main() -> None:
    runs = [simulate(scenario) for scenario in SCENARIOS]
    make_plot(runs)

    run_by_name = {run["scenario"]: run for run in runs}
    exit_run = run_by_name["exit_after_closure"]
    continue_run = run_by_name["continue_after_closure"]
    verification = {
        "exit_removes_all_closed_paths": all(
            snapshot["active_closed_paths"] == 0
            for snapshot in exit_run["snapshots"]
        ),
        "control_retains_closed_paths": any(
            snapshot["active_closed_paths"] > 0
            for snapshot in continue_run["snapshots"]
        ),
        "exit_writes_history_records": exit_run["final_history_records"] > 0,
        "exit_writes_zero_modes": exit_run["final_zero_modes"] > 0,
        "exit_leaves_no_active_paths_after_lifetime": (
            exit_run["final_active_paths"] == 0
        ),
    }

    payload = {
        "model": {
            "sites": len(INITIAL_WORDS),
            "initial_words": list(INITIAL_WORDS),
            "layer_lifetime": LAYER_LIFETIME,
            "epochs": EPOCHS,
            "step_alphabet": ["+", "-"],
            "closure_rule": "cumulative balance returns to zero",
            "probabilities": "none",
            "mutations": "none",
            "rules": {
                scenario.name: scenario.label
                for scenario in SCENARIOS
            },
        },
        "runs": runs,
        "verification": verification,
    }
    OUT_JSON.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print("closure-exit experiment")
    for key, value in verification.items():
        print(f"  [{'v' if value else 'x'}] {key}")
    for run in runs:
        print(
            f"{run['scenario']:>24}  "
            f"active={run['final_active_paths']:>4}  "
            f"closed_active={run['final_active_closed_paths']:>3}  "
            f"history={run['final_history_records']:>4}  "
            f"zero_events={run['final_zero_events']:>5}  "
            f"zero_modes={run['final_zero_modes']:>3}  "
            f"dead={run['final_dead_paths']:>4}"
        )
    print(f"\nplot: {OUT_PNG}")
    print(f"json: {OUT_JSON}")


if __name__ == "__main__":
    main()
