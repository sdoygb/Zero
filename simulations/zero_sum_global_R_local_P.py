#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Global result layer R and local history layer P
==============================================

This sandbox tests the conjecture:

    R is global
    P is local

There are several local sites. Each site has its own history library P_i.
The result layer R is always global in the main scenarios.

At every epoch, each active history h in a site produces all balanced
extensions e up to period L. The new history h+e is recorded locally in P_i,
and its cyclic result is recorded in the global R.

The scenarios differ only in how a site rebuilds its active library:

    local_P:
        rebuild from its own local P_i

    global_P:
        rebuild from the union of all histories

    global_R_reconstruct:
        rebuild from one global canonical history for each result in R

No probabilities or mutation rates are used.

Run:
    python3 simulations/zero_sum_global_R_local_P.py
"""
from __future__ import annotations

import json
import math
import os
import tempfile
from collections import Counter
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
OUT_JSON = ROOT / "simulations" / "zero_sum_global_R_local_P_results.json"
OUT_PNG = ROOT / "simulations" / "zero_sum_global_R_local_P.png"

LAYER_LIFETIME = 4
EPOCHS = 4
INITIAL_HISTORIES = (
    "+-",
    "++--",
    "+-+-+-",
    "++-+--",
    "+++---",
    "++--++--",
)

INITIALIZATION_MODES = {
    "distinct": INITIAL_HISTORIES,
    "shared": ("+-",) * len(INITIAL_HISTORIES),
    "one_local_variant": ("+-",) * (len(INITIAL_HISTORIES) - 1)
    + ("++--",),
}


@dataclass(frozen=True)
class Scenario:
    name: str
    rebuild_rule: str
    label: str


SCENARIOS = (
    Scenario(
        name="local_P",
        rebuild_rule="local_P",
        label="global R, local P",
    ),
    Scenario(
        name="global_P",
        rebuild_rule="global_P",
        label="global R, global P control",
    ),
    Scenario(
        name="global_R_reconstruct",
        rebuild_rule="global_R_reconstruct",
        label="global R reconstructs history",
    ),
)


def word_from_string(text: str) -> tuple[int, ...]:
    word = tuple(1 if char == "+" else -1 for char in text)
    if sum(word) != 0:
        raise ValueError(f"history is not zero-sum: {text}")
    return word


def word_to_string(word: tuple[int, ...]) -> str:
    if not word:
        return "0"
    return "".join("+" if value == 1 else "-" for value in word)


def balanced_words(period: int) -> tuple[tuple[int, ...], ...]:
    words = []
    for positions in combinations(range(period), period // 2):
        positive = set(positions)
        words.append(
            tuple(
                1 if index in positive else -1
                for index in range(period)
            )
        )
    return tuple(words)


def jaccard(left: set[tuple[int, ...]], right: set[tuple[int, ...]]) -> float:
    if not left and not right:
        return 1.0
    union = left | right
    if not union:
        return 1.0
    return len(left & right) / len(union)


def mean_pairwise_overlap(libraries: list[set[tuple[int, ...]]]) -> float:
    values = []
    for left in range(len(libraries)):
        for right in range(left + 1, len(libraries)):
            values.append(jaccard(libraries[left], libraries[right]))
    return float(np.mean(values)) if values else 1.0


def shannon_counts(counter: Counter) -> float:
    total = float(sum(counter.values()))
    if total <= 0.0:
        return 0.0
    probabilities = np.asarray(list(counter.values()), dtype=float) / total
    probabilities = probabilities[probabilities > 0.0]
    value = float(-np.sum(probabilities * np.log(probabilities)))
    return 0.0 if abs(value) < 1.0e-15 else value


def top_result(
    histories: set[tuple[int, ...]],
) -> tuple[str, int]:
    counter = Counter(canonical_cycle(history) for history in histories)
    if not counter:
        return "", 0
    result, _ = counter.most_common(1)[0]
    return word_to_string(result), len(counter)


def simulate(
    scenario: Scenario,
    extensions: tuple[tuple[int, ...], ...],
    initial_histories: tuple[str, ...] = INITIAL_HISTORIES,
) -> dict[str, object]:
    local_histories = [
        {word_from_string(initial)}
        for initial in initial_histories
    ]
    active_histories = [
        set(local_histories[index])
        for index in range(len(initial_histories))
    ]

    global_results: set[tuple[int, ...]] = set()
    global_histories: set[tuple[int, ...]] = set()
    representative_history: dict[tuple[int, ...], tuple[int, ...]] = {}

    snapshots = []
    for epoch in range(EPOCHS + 1):
        local_overlap = mean_pairwise_overlap(local_histories)
        active_overlap = mean_pairwise_overlap(active_histories)
        top_results = [
            top_result(histories)
            for histories in local_histories
        ]
        distinct_top_results = len({item[0] for item in top_results if item[0]})

        snapshots.append(
            {
                "epoch": epoch,
                "global_results": len(global_results),
                "global_histories": len(global_histories),
                "local_history_overlap": local_overlap,
                "active_history_overlap": active_overlap,
                "distinct_top_results": distinct_top_results,
                "local_history_sizes": [
                    len(histories) for histories in local_histories
                ],
                "top_results": [item[0] for item in top_results],
            }
        )

        if epoch == EPOCHS:
            break

        produced_per_site = []
        for site, active in enumerate(active_histories):
            produced = set()
            for history in active:
                for extension in extensions:
                    child = history + extension
                    produced.add(child)
                    result = canonical_cycle(child)
                    global_results.add(result)
                    representative_history.setdefault(result, child)
            local_histories[site].update(produced)
            global_histories.update(produced)
            produced_per_site.append(produced)

        if scenario.rebuild_rule == "local_P":
            active_histories = [
                set(local_histories[index])
                for index in range(len(local_histories))
            ]
        elif scenario.rebuild_rule == "global_P":
            active_histories = [
                set(global_histories)
                for _ in range(len(local_histories))
            ]
        elif scenario.rebuild_rule == "global_R_reconstruct":
            reconstruction = set(representative_history.values())
            active_histories = [
                set(reconstruction)
                for _ in range(len(local_histories))
            ]
        else:
            raise ValueError(f"unknown rebuild rule: {scenario.rebuild_rule}")

    final = snapshots[-1]
    return {
        "scenario": scenario.name,
        "label": scenario.label,
        "epochs": EPOCHS,
        "layer_lifetime": LAYER_LIFETIME,
        "initial_histories": list(initial_histories),
        "snapshots": snapshots,
        "final_global_results": final["global_results"],
        "final_global_histories": final["global_histories"],
        "final_local_history_overlap": final["local_history_overlap"],
        "final_active_overlap": final["active_history_overlap"],
        "final_distinct_top_results": final["distinct_top_results"],
        "final_local_history_sizes": final["local_history_sizes"],
        "final_top_results": final["top_results"],
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
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8))
    colors = {
        "local_P": "#2B7A78",
        "global_P": "#C77C28",
        "global_R_reconstruct": "#B44A5A",
    }

    for run in runs:
        epochs = [snapshot["epoch"] for snapshot in run["snapshots"]]
        axes[0].plot(
            epochs,
            [snapshot["local_history_overlap"] for snapshot in run["snapshots"]],
            linewidth=1.9,
            color=colors[run["scenario"]],
            marker="o",
            markersize=3,
            label=run["label"],
        )
    axes[0].set_ylim(0.0, 1.02)
    axes[0].set_title("History overlap between sites")
    axes[0].set_xlabel("epoch")
    axes[0].set_ylabel("mean pairwise Jaccard")
    axes[0].grid(alpha=0.22)
    axes[0].legend(frameon=False)

    for run in runs:
        epochs = [snapshot["epoch"] for snapshot in run["snapshots"]]
        axes[1].plot(
            epochs,
            [snapshot["distinct_top_results"] for snapshot in run["snapshots"]],
            linewidth=1.9,
            color=colors[run["scenario"]],
            marker="o",
            markersize=3,
            label=run["label"],
        )
    axes[1].set_title("Distinct local leading results")
    axes[1].set_xlabel("epoch")
    axes[1].set_ylabel("number of distinct top results")
    axes[1].grid(alpha=0.22)
    axes[1].legend(frameon=False)

    labels = []
    global_results = []
    global_histories = []
    for run in runs:
        labels.append(run["scenario"].replace("_", "\n"))
        global_results.append(run["final_global_results"])
        global_histories.append(run["final_global_histories"])
    x = np.arange(len(labels))
    width = 0.38
    axes[2].bar(
        x - width / 2,
        global_results,
        width=width,
        label="global R results",
        color="#3A6EA5",
    )
    axes[2].bar(
        x + width / 2,
        global_histories,
        width=width,
        label="global P histories",
        color="#C77C28",
    )
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(labels, fontsize=7)
    axes[2].set_title("Final record counts")
    axes[2].set_ylabel("count")
    axes[2].grid(axis="y", alpha=0.22)
    axes[2].legend(frameon=False)

    fig.suptitle("Global R versus local P experiment", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT_PNG, dpi=180)
    plt.close(fig)


def main() -> None:
    extensions = tuple(
        word
        for period in range(2, LAYER_LIFETIME + 1, 2)
        for word in balanced_words(period)
    )
    runs = [simulate(scenario, extensions) for scenario in SCENARIOS]
    initialization_controls = {
        mode: [
            simulate(scenario, extensions, starts)
            for scenario in SCENARIOS
        ]
        for mode, starts in INITIALIZATION_MODES.items()
    }
    make_plot(runs)

    payload = {
        "model": {
            "sites": len(INITIAL_HISTORIES),
            "initial_histories": list(INITIAL_HISTORIES),
            "layer_lifetime": LAYER_LIFETIME,
            "epochs": EPOCHS,
            "extensions": len(extensions),
            "probabilities": "none",
            "mutations": "none",
            "rules": {
                scenario.name: scenario.label for scenario in SCENARIOS
            },
            "note": (
                "R is global in all main scenarios. The only difference is "
                "whether a site rebuilds from local P, global P, or a global "
                "canonical history reconstructed from R."
            ),
        },
        "runs": runs,
        "initialization_controls": initialization_controls,
    }
    OUT_JSON.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print("global R / local P experiment")
    for run in runs:
        print(
            f"{run['scenario']:>22}  "
            f"R={run['final_global_results']:>4}  "
            f"P={run['final_global_histories']:>5}  "
            f"history_overlap={run['final_local_history_overlap']:.3f}  "
            f"top_result_choices={run['final_distinct_top_results']}  "
            f"tops={run['final_top_results']}"
        )
    print("\ninitialization controls")
    for mode, control_runs in initialization_controls.items():
        for run in control_runs:
            print(
                f"{mode:>18}  {run['scenario']:>22}  "
                f"history_overlap={run['final_local_history_overlap']:.3f}  "
                f"top_result_choices={run['final_distinct_top_results']}  "
                f"tops={run['final_top_results']}"
            )
    print(f"\nplot: {OUT_PNG}")
    print(f"json: {OUT_JSON}")


if __name__ == "__main__":
    main()
