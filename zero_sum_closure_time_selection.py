#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Zero-sum closure-time selection
===============================

This sandbox follows the rule:

    all zero-sum fluctuation modes are present from the beginning
    faster closure means more completed cycles per internal time
    a completed closure may leave a descendant

There is no mutation probability, no cost parameter, and no split
multiplicity.

The primitive object is a balanced cyclic word whose first positive return
to zero is at its end. Its closure time T is the length of that primitive
return. Longer words that immediately repeat a shorter primitive return are
not counted as new modes; they are repeated cycles of the same mode.

Two rules are compared:

    persist_only:
        w -> w after one closure

    closure_daughter:
        w -> w + w after one closure

The second rule is the minimal explicit assumption that closure produces a
descendant. It is an input, not a consequence of zero-sum alone.

Run:
    python3 simulations/zero_sum_closure_time_selection.py
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
OUT_JSON = ROOT / "simulations" / "zero_sum_closure_time_selection_results.json"
OUT_PNG = ROOT / "simulations" / "zero_sum_closure_time_selection.png"

MAX_PERIOD = 18
STEPS = 120
POPULATION_SCALE_LIMIT = 1.0e9


@dataclass(frozen=True)
class Scenario:
    name: str
    daughters_per_closure: int
    label: str


SCENARIOS = (
    Scenario(
        name="persist_only",
        daughters_per_closure=0,
        label="persist only",
    ),
    Scenario(
        name="closure_daughter",
        daughters_per_closure=1,
        label="one daughter per closure",
    ),
)


def first_return_time(word: tuple[int, ...]) -> int:
    total = 0
    for index, value in enumerate(word, start=1):
        total += value
        if total == 0:
            return index
    raise ValueError("word never returns to zero")


def primitive_closure(word: tuple[int, ...]) -> tuple[int, ...]:
    return canonical_cycle(word[: first_return_time(word)])


def primitive_modes(max_period: int = MAX_PERIOD) -> tuple[tuple[int, ...], ...]:
    modes = set()
    for period in range(2, max_period + 1, 2):
        for positions in combinations(range(period), period // 2):
            positive = set(positions)
            word = tuple(
                1 if index in positive else -1
                for index in range(period)
            )
            modes.add(primitive_closure(word))
    return tuple(sorted(modes, key=lambda item: (len(item), item)))


def word_to_string(word: tuple[int, ...]) -> str:
    return "".join("+" if value == 1 else "-" for value in word)


def shannon(counts: np.ndarray) -> float:
    total = float(np.sum(counts))
    if total <= 0.0:
        return 0.0
    probabilities = counts[counts > 0.0] / total
    value = float(-np.sum(probabilities * np.log(probabilities)))
    return 0.0 if abs(value) < 1.0e-15 else value


def simulate(
    scenario: Scenario,
    modes: tuple[tuple[int, ...], ...],
) -> dict[str, object]:
    cohorts: dict[int, np.ndarray] = {}
    for index, word in enumerate(modes):
        cohorts[index] = np.zeros(len(word), dtype=float)
        cohorts[index][0] = 1.0

    snapshots = []
    rescale_steps = 0
    closure_time_histories: dict[int, list[dict[str, float]]] = defaultdict(list)

    for step in range(STEPS + 1):
        population = np.array(
            [float(np.sum(cohorts[index])) for index in range(len(modes))],
            dtype=float,
        )
        total = float(np.sum(population))
        if total > 0.0:
            shares = population / total
            top_index = int(np.argmax(population))
            mean_closure_time = float(
                sum(shares[index] * len(modes[index]) for index in range(len(modes)))
            )
            mean_rate = float(
                sum(
                    shares[index]
                    * (
                        math.log(1.0 + scenario.daughters_per_closure)
                        / len(modes[index])
                    )
                    for index in range(len(modes))
                )
            )
            top_fraction = float(population[top_index] / total)
            top_word = word_to_string(modes[top_index])
            lineages = int(np.count_nonzero(population > 1.0e-12 * total))
        else:
            top_index = -1
            mean_closure_time = 0.0
            mean_rate = 0.0
            top_fraction = 0.0
            top_word = ""
            lineages = 0

        by_period = Counter()
        if total > 0.0:
            for index, word in enumerate(modes):
                by_period[len(word)] += float(population[index] / total)
        for closure_time in sorted({len(word) for word in modes}):
            closure_time_histories[closure_time].append(
                {
                    "time": step,
                    "share": by_period.get(closure_time, 0.0),
                }
            )

        snapshots.append(
            {
                "time": step,
                "population": total,
                "lineages": lineages,
                "diversity": shannon(population),
                "mean_closure_time": mean_closure_time,
                "mean_rate": mean_rate,
                "top_word": top_word,
                "top_fraction": top_fraction,
            }
        )

        if step == STEPS:
            break

        next_cohorts: dict[int, np.ndarray] = {}
        for index, word in enumerate(modes):
            next_cohorts[index] = np.zeros(len(word), dtype=float)

        for index, word in enumerate(modes):
            current = cohorts[index]
            period = len(word)
            if period > 1:
                next_cohorts[index][1:] += current[:-1]
            finishers = float(current[-1])
            if finishers <= 0.0:
                continue
            next_cohorts[index][0] += finishers
            next_cohorts[index][0] += (
                scenario.daughters_per_closure * finishers
            )

        cohorts = next_cohorts
        total = float(sum(float(np.sum(cohort)) for cohort in cohorts.values()))
        if total > POPULATION_SCALE_LIMIT:
            scale = POPULATION_SCALE_LIMIT / total
            for index in cohorts:
                cohorts[index] *= scale
            rescale_steps += 1

    final_population = np.array(
        [float(np.sum(cohorts[index])) for index in range(len(modes))],
        dtype=float,
    )
    total = float(np.sum(final_population))
    top_entries = []
    if total > 0.0:
        order = np.argsort(final_population)[::-1]
        for index in order[:15]:
            if final_population[index] <= 1.0e-12 * total:
                continue
            top_entries.append(
                {
                    "word": word_to_string(modes[int(index)]),
                    "closure_time": len(modes[int(index)]),
                    "fraction": float(final_population[index] / total),
                }
            )

    return {
        "scenario": scenario.name,
        "label": scenario.label,
        "daughters_per_closure": scenario.daughters_per_closure,
        "snapshots": snapshots,
        "closure_time_histories": closure_time_histories,
        "final_population": total,
        "final_lineages": int(np.count_nonzero(final_population > 0.0)),
        "final_diversity": shannon(final_population),
        "top_entries": top_entries,
        "rescale_steps": rescale_steps,
    }


def summarize(runs: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    summary = {}
    for run in runs:
        final = run["snapshots"][-1]
        summary[run["scenario"]] = {
            "label": run["label"],
            "daughters_per_closure": run["daughters_per_closure"],
            "final_population": run["final_population"],
            "final_lineages": run["final_lineages"],
            "final_diversity": run["final_diversity"],
            "mean_closure_time": final["mean_closure_time"],
            "mean_rate": final["mean_rate"],
            "top_word": final["top_word"],
            "top_fraction": final["top_fraction"],
            "rescale_steps": run["rescale_steps"],
        }
    return summary


def make_plot(
    runs: list[dict[str, object]],
    modes: tuple[tuple[int, ...], ...],
) -> None:
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
        "persist_only": "#777777",
        "closure_daughter": "#2B7A78",
    }

    for run in runs:
        times = [snapshot["time"] for snapshot in run["snapshots"]]
        populations = [
            max(float(snapshot["population"]), 0.5)
            for snapshot in run["snapshots"]
        ]
        axes[0, 0].plot(
            times,
            populations,
            linewidth=1.8,
            color=colors[run["scenario"]],
            label=run["label"],
        )
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_title("Total population")
    axes[0, 0].set_xlabel("internal time")
    axes[0, 0].set_ylabel("closed lineages")
    axes[0, 0].grid(alpha=0.22)
    axes[0, 0].legend(frameon=False)

    for run in runs:
        times = [snapshot["time"] for snapshot in run["snapshots"]]
        values = [snapshot["mean_closure_time"] for snapshot in run["snapshots"]]
        axes[0, 1].plot(
            times,
            values,
            linewidth=1.8,
            color=colors[run["scenario"]],
            label=run["label"],
        )
    axes[0, 1].set_title("Mean closure time")
    axes[0, 1].set_xlabel("internal time")
    axes[0, 1].set_ylabel("mean T")
    axes[0, 1].grid(alpha=0.22)
    axes[0, 1].legend(frameon=False)

    selected_run = next(
        run for run in runs if run["scenario"] == "closure_daughter"
    )
    times = [0, 20, 60, 120]
    periods = sorted({len(word) for word in modes})
    width = 0.18
    x = np.arange(len(periods))
    for offset, time in enumerate(times):
        snapshot_index = min(
            time,
            len(selected_run["closure_time_histories"][periods[0]]) - 1,
        )
        shares = [
            selected_run["closure_time_histories"][period][snapshot_index]["share"]
            for period in periods
        ]
        axes[1, 0].bar(
            x + (offset - 1.5) * width,
            shares,
            width=width,
            label=f"t={time}",
        )
    axes[1, 0].set_xticks(x)
    axes[1, 0].set_xticklabels([str(period) for period in periods])
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_title("Closure-time distribution")
    axes[1, 0].set_xlabel("closure time T")
    axes[1, 0].set_ylabel("population share")
    axes[1, 0].grid(axis="y", alpha=0.22)
    axes[1, 0].legend(frameon=False, fontsize=7)

    top = selected_run["top_entries"][:12]
    labels = [entry["word"] for entry in top]
    values = [entry["fraction"] for entry in top]
    y = np.arange(len(top))
    axes[1, 1].barh(y, values, color="#C77C28", alpha=0.9)
    axes[1, 1].set_yticks(y)
    axes[1, 1].set_yticklabels(labels, fontsize=7)
    axes[1, 1].invert_yaxis()
    axes[1, 1].set_title("Leading closure modes")
    axes[1, 1].set_xlabel("final population share")
    axes[1, 1].grid(axis="x", alpha=0.22)

    fig.suptitle("Zero-sum closure-time selection", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUT_PNG, dpi=180)
    plt.close(fig)


def main() -> None:
    modes = primitive_modes()
    runs = [simulate(scenario, modes) for scenario in SCENARIOS]
    summary = summarize(runs)
    make_plot(runs, modes)

    payload = {
        "model": {
            "zero_constraint": "every mode is a balanced primitive return",
            "initial_condition": "one lineage of every mode",
            "closure_time": "first positive return to zero",
            "repeated_words": "reduced to their primitive closure mode",
            "probabilities": "none",
            "mutation": "none",
            "rules": {
                scenario.name: scenario.label for scenario in SCENARIOS
            },
            "steps": STEPS,
            "max_period": MAX_PERIOD,
            "population_scale_limit": POPULATION_SCALE_LIMIT,
            "note": (
                "The daughter rule is an explicit minimal reproduction "
                "assumption. Population rescaling is numerical only."
            ),
        },
        "summary": summary,
        "runs": runs,
    }
    OUT_JSON.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"\nmodes: {len(modes)}")
    print(f"plot: {OUT_PNG}")
    print(f"json: {OUT_JSON}")


if __name__ == "__main__":
    main()
