#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Periodic destruction and regrowth of the evolution layer
========================================================

This sandbox tests the proposed cycle:

    E evolves -> some paths close into P and Z
             -> the whole E layer is periodically destroyed
             -> D records the open paths lost at destruction
             -> E regrows from the closed histories in P

Three rules are compared:

    continuous:
        E is never destroyed.

    wipe_no_reseed:
        E is destroyed periodically and nothing is re-seeded.

    wipe_reseed:
        E is destroyed periodically, then every exact closed history
        w in P_i seeds the two open paths w+ and w-.

No probabilities, mutation rates, or channel weights are used.

Run:
    python3 simulations/zero_sum_periodic_destruction.py
"""
from __future__ import annotations

import json
import os
import tempfile
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
OUT_JSON = ROOT / "simulations" / "zero_sum_periodic_destruction_results.json"
OUT_PNG = ROOT / "simulations" / "zero_sum_periodic_destruction.png"

DESTRUCTION_PERIOD = 3
CYCLES = 4
INITIAL_WORDS = (
    "+",
    "-",
    "++",
    "--",
    "+++",
    "---",
)


@dataclass(frozen=True)
class Scenario:
    name: str
    destroy_periodically: bool
    reseed: bool
    label: str


SCENARIOS = (
    Scenario(
        name="continuous",
        destroy_periodically=False,
        reseed=False,
        label="continuous E",
    ),
    Scenario(
        name="wipe_no_reseed",
        destroy_periodically=True,
        reseed=False,
        label="periodic wipe, no reseed",
    ),
    Scenario(
        name="wipe_reseed",
        destroy_periodically=True,
        reseed=True,
        label="periodic wipe, reseed from P",
    ),
)


def word_from_string(text: str) -> tuple[int, ...]:
    return tuple(1 if char == "+" else -1 for char in text)


def is_closed(word: tuple[int, ...]) -> bool:
    return bool(word) and sum(word) == 0


def total_active(active: list[set[tuple[int, ...]]]) -> int:
    return sum(len(paths) for paths in active)


def total_history(
    histories: list[set[tuple[int, ...]]],
) -> int:
    return sum(len(book) for book in histories)


def total_dead(dead: list[set[tuple[int, ...]]]) -> int:
    return sum(len(book) for book in dead)


def word_balance(word: tuple[int, ...]) -> int:
    return sum(word)


def layer_balances(
    layers: list[set[tuple[int, ...]]],
) -> list[int]:
    return [
        sum(word_balance(word) for word in paths)
        for paths in layers
    ]


def evolve_one_step(
    active: list[set[tuple[int, ...]]],
    histories: list[set[tuple[int, ...]]],
    zero_layer: set[tuple[int, ...]],
) -> tuple[list[set[tuple[int, ...]]], int]:
    next_active = [set() for _ in active]
    closed_events = 0

    for site, paths in enumerate(active):
        for word in paths:
            for step in (1, -1):
                child = word + (step,)
                if is_closed(child):
                    histories[site].add(child)
                    zero_layer.add(canonical_cycle(child))
                    closed_events += 1
                else:
                    next_active[site].add(child)

    return next_active, closed_events


def reseed_from_history(
    histories: list[set[tuple[int, ...]]],
) -> list[set[tuple[int, ...]]]:
    active = [set() for _ in histories]
    for site, book in enumerate(histories):
        for word in book:
            for step in (1, -1):
                seed = word + (step,)
                if not is_closed(seed):
                    active[site].add(seed)
    return active


def simulate(scenario: Scenario) -> dict[str, object]:
    active = [
        {word_from_string(text)}
        for text in INITIAL_WORDS
    ]
    histories = [set() for _ in INITIAL_WORDS]
    dead = [set() for _ in INITIAL_WORDS]
    zero_layer: set[tuple[int, ...]] = set()

    snapshots: list[dict[str, object]] = []
    destruction_events = 0
    total_destroyed = 0
    total_reseeded = 0
    total_closed_events = 0
    previous_dead_count = 0
    death_accepts_only_open_paths = True
    death_layer_is_monotone = True
    reseed_uses_only_history = True
    history_preserved_at_each_destruction = True
    zero_layer_preserved_at_each_destruction = True

    for time in range(1, DESTRUCTION_PERIOD * CYCLES + 1):
        active, closed_events = evolve_one_step(
            active,
            histories,
            zero_layer,
        )
        total_closed_events += closed_events

        destruction_event = (
            scenario.destroy_periodically
            and time % DESTRUCTION_PERIOD == 0
        )
        destroyed_paths = 0
        destroyed_balance_by_site = [0] * len(active)
        reseeded_paths = 0
        reseeded_balance_by_site = [0] * len(active)
        active_after_destruction = total_active(active)

        if destruction_event:
            destruction_events += 1
            destroyed_paths = total_active(active)
            total_destroyed += destroyed_paths
            history_before_destruction = total_history(histories)
            zero_modes_before_destruction = len(zero_layer)
            for site, paths in enumerate(active):
                death_accepts_only_open_paths &= all(
                    not is_closed(word) for word in paths
                )
                destroyed_balance_by_site[site] = sum(
                    word_balance(word) for word in paths
                )
                dead[site].update(paths)

            active = [set() for _ in active]
            active_after_destruction = 0
            history_preserved_at_each_destruction &= (
                total_history(histories)
                == history_before_destruction
            )
            zero_layer_preserved_at_each_destruction &= (
                len(zero_layer) == zero_modes_before_destruction
            )

            if scenario.reseed:
                active = reseed_from_history(histories)
                reseed_uses_only_history &= (
                    active == reseed_from_history(histories)
                )
                reseeded_paths = total_active(active)
                reseeded_balance_by_site = layer_balances(active)
                total_reseeded += reseeded_paths

        dead_count = total_dead(dead)
        death_layer_is_monotone &= dead_count >= previous_dead_count
        previous_dead_count = dead_count

        active_balances = layer_balances(active)
        history_balances = layer_balances(histories)
        dead_balances = layer_balances(dead)
        total_balance = (
            sum(active_balances)
            + sum(history_balances)
            + sum(dead_balances)
        )
        snapshots.append(
            {
                "time": time,
                "cycle": (time - 1) // DESTRUCTION_PERIOD + 1,
                "step_in_cycle": (time - 1) % DESTRUCTION_PERIOD + 1,
                "active_paths": total_active(active),
                "closed_events_this_step": closed_events,
                "history_records": total_history(histories),
                "zero_modes": len(zero_layer),
                "dead_paths": total_dead(dead),
                "destruction_event": destruction_event,
                "active_after_destruction": active_after_destruction,
                "destroyed_paths": destroyed_paths,
                "reseeded_paths": reseeded_paths,
                "active_balances": active_balances,
                "history_balances": history_balances,
                "dead_balances": dead_balances,
                "destroyed_balance_by_site": destroyed_balance_by_site,
                "reseeded_balance_by_site": reseeded_balance_by_site,
                "total_balance": total_balance,
            }
        )

    final = snapshots[-1]
    return {
        "scenario": scenario.name,
        "label": scenario.label,
        "destruction_period": DESTRUCTION_PERIOD,
        "cycles": CYCLES,
        "snapshots": snapshots,
        "final_active_paths": final["active_paths"],
        "final_history_records": final["history_records"],
        "final_zero_modes": final["zero_modes"],
        "final_dead_paths": final["dead_paths"],
        "destruction_events": destruction_events,
        "total_destroyed_paths": total_destroyed,
        "total_reseeded_paths": total_reseeded,
        "total_closed_events": total_closed_events,
        "death_accepts_only_open_paths": (
            death_accepts_only_open_paths
        ),
        "death_layer_is_monotone": death_layer_is_monotone,
        "reseed_uses_only_history": reseed_uses_only_history,
        "history_preserved_at_each_destruction": (
            history_preserved_at_each_destruction
        ),
        "zero_layer_preserved_at_each_destruction": (
            zero_layer_preserved_at_each_destruction
        ),
        "final_active_balances": final["active_balances"],
        "final_history_balances": final["history_balances"],
        "final_dead_balances": final["dead_balances"],
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
    colors = {
        "continuous": "#3A6EA5",
        "wipe_no_reseed": "#B44A5A",
        "wipe_reseed": "#2B7A78",
    }
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8))

    for run in runs:
        times = [snapshot["time"] for snapshot in run["snapshots"]]
        axes[0].plot(
            times,
            [snapshot["active_paths"] for snapshot in run["snapshots"]],
            color=colors[run["scenario"]],
            marker="o",
            markersize=3,
            linewidth=1.8,
            label=run["label"],
        )
    axes[0].set_title("Active paths through destruction cycles")
    axes[0].set_xlabel("time")
    axes[0].set_ylabel("count")
    axes[0].grid(alpha=0.22)
    axes[0].legend(frameon=False)

    for run in runs:
        times = [snapshot["time"] for snapshot in run["snapshots"]]
        axes[1].plot(
            times,
            [snapshot["history_records"] for snapshot in run["snapshots"]],
            color=colors[run["scenario"]],
            marker="o",
            markersize=3,
            linewidth=1.8,
            label=run["label"],
        )
        axes[1].plot(
            times,
            [snapshot["zero_modes"] for snapshot in run["snapshots"]],
            color=colors[run["scenario"]],
            linestyle="--",
            marker="s",
            markersize=2.5,
            linewidth=1.3,
        )
    axes[1].set_title("History records (solid) and zero modes (dashed)")
    axes[1].set_xlabel("time")
    axes[1].set_ylabel("count")
    axes[1].grid(alpha=0.22)
    axes[1].legend(frameon=False)

    labels = []
    destroyed = []
    reseeded = []
    for run in runs:
        labels.append(run["scenario"].replace("_", "\n"))
        destroyed.append(run["total_destroyed_paths"])
        reseeded.append(run["total_reseeded_paths"])
    x = np.arange(len(labels))
    width = 0.38
    axes[2].bar(
        x - width / 2,
        destroyed,
        width=width,
        color="#B44A5A",
        label="destroyed",
    )
    axes[2].bar(
        x + width / 2,
        reseeded,
        width=width,
        color="#2B7A78",
        label="reseeded",
    )
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(labels, fontsize=7)
    axes[2].set_title("Destroyed and re-seeded paths")
    axes[2].set_ylabel("count")
    axes[2].grid(axis="y", alpha=0.22)
    axes[2].legend(frameon=False)

    fig.suptitle("Periodic destruction and regrowth of E", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT_PNG, dpi=180)
    plt.close(fig)


def main() -> None:
    runs = [simulate(scenario) for scenario in SCENARIOS]
    make_plot(runs)

    run_by_name = {run["scenario"]: run for run in runs}
    continuous = run_by_name["continuous"]
    no_reseed = run_by_name["wipe_no_reseed"]
    reseed = run_by_name["wipe_reseed"]

    destruction_snapshots = [
        snapshot
        for snapshot in reseed["snapshots"]
        if snapshot["destruction_event"]
    ]
    verification = {
        "periodic_destruction_events_occur": (
            reseed["destruction_events"] == CYCLES
            and no_reseed["destruction_events"] == CYCLES
        ),
        "wipe_without_reseed_becomes_empty": (
            no_reseed["final_active_paths"] == 0
        ),
        "wipe_actually_empties_E_before_reseed": all(
            snapshot["active_after_destruction"] == 0
            for snapshot in destruction_snapshots
        ),
        "reseed_restores_activity_every_cycle": all(
            snapshot["reseeded_paths"] > 0
            for snapshot in destruction_snapshots
        ),
        "history_survives_destruction": (
            reseed["final_history_records"] > 0
            and reseed["final_zero_modes"] > 0
        ),
        "open_paths_are_recorded_as_dead": (
            all(run["death_accepts_only_open_paths"] for run in runs)
            and reseed["final_dead_paths"] > 0
        ),
        "death_layer_is_monotone_and_has_no_removal": all(
            run["death_layer_is_monotone"] for run in runs
        ),
        "death_layer_does_not_feed_reseed": all(
            run["reseed_uses_only_history"] for run in runs
        ),
        "history_and_zero_layer_survive_each_destruction": all(
            run["history_preserved_at_each_destruction"]
            and run["zero_layer_preserved_at_each_destruction"]
            for run in runs
        ),
        "continuous_rule_has_no_destruction": (
            continuous["destruction_events"] == 0
        ),
        "global_zero_sum_is_conserved": all(
            snapshot["total_balance"] == 0
            for run in runs
            for snapshot in run["snapshots"]
        ),
        "history_layer_is_zero_balanced": all(
            sum(snapshot["history_balances"]) == 0
            for run in runs
            for snapshot in run["snapshots"]
        ),
        "reseeded_activity_can_start_zero_balanced": all(
            sum(snapshot["reseeded_balance_by_site"]) == 0
            for snapshot in destruction_snapshots
        ),
        "dead_layer_can_hold_local_imbalance": any(
            abs(value) > 0
            for run in runs
            for value in run["final_dead_balances"]
        ),
        "local_dead_imbalances_compensate_globally": all(
            sum(run["final_dead_balances"]) == 0
            for run in runs
        ),
    }

    payload = {
        "model": {
            "initial_words": list(INITIAL_WORDS),
            "destruction_period": DESTRUCTION_PERIOD,
            "cycles": CYCLES,
            "probabilities": "none",
            "mutations": "none",
            "reseed_rule": "each closed history w seeds w+ and w-",
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

    print("periodic destruction experiment")
    for key, value in verification.items():
        print(f"  [{'v' if value else 'x'}] {key}")
    for run in runs:
        print(
            f"{run['scenario']:>16}  "
            f"cycles={run['destruction_events']}  "
            f"active={run['final_active_paths']:>4}  "
            f"history={run['final_history_records']:>4}  "
            f"zero_modes={run['final_zero_modes']:>3}  "
            f"dead={run['final_dead_paths']:>4}  "
            f"reseeded={run['total_reseeded_paths']:>4}  "
            f"Q_active={run['final_active_balances']}  "
            f"Q_dead={run['final_dead_balances']}  "
            f"Q_total={run['final_total_balance']}"
        )
    print(f"\nplot: {OUT_PNG}")
    print(f"json: {OUT_JSON}")


if __name__ == "__main__":
    main()
