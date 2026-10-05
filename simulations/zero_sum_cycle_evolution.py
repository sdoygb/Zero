#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Zero-sum cycle sandbox
======================

This is an exploratory model. It does not use any D* derivation.

The model keeps one microscopic rule exact:

    sum(charges) = 0.

A closed fluctuation is represented by a cyclic word of +1 and -1 steps.
Open words have no lineage. A closed word returns to zero. Every internal
return to zero is a possible split point.

Three closure algebras are compared:

    sterile       M = 1
    single_cut    M = 1 + r
    all_cuts      M = 2**r

Here r is the number of internal returns to zero and M is the number of
inherited closure programs produced after one period. No free fertility
parameter is introduced. The lineage exponent is

    lambda = log(M) / T.

Important limitation:

M counts abstract branch programs or descriptions. Treating each program as
a full new loop assumes an additional reproduction or copy operation. That
operation is not derived from the zero-sum condition. The later audit script
zero_sum_reproduction_audit.py separates branch count, child-object count,
and actual population growth.

Mutation is an optional deterministic kernel. It only uses charge-preserving
operations:

    swap adjacent +- or -+
    insert an adjacent +- or -+ pair
    delete an adjacent +- or -+ pair

The last two operations change the cycle length. Thus mutation can explore
the space of zero-sum cycles rather than a pre-labelled fitness landscape.

Run:
    python3 simulations/zero_sum_cycle_evolution.py
"""
from __future__ import annotations

import json
import math
import os
import tempfile
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


ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "simulations" / "zero_sum_cycle_evolution_results.json"
OUT_PNG = ROOT / "simulations" / "zero_sum_cycle_evolution.png"

MAX_PERIOD = 12
MIN_PERIOD = 4
MUTATION_PROBABILITY = 0.03
STEPS = 300
POPULATION_SCALE_LIMIT = 1.0e6


def word_from_string(text: str) -> tuple[int, ...]:
    values = tuple(1 if char == "+" else -1 for char in text)
    if sum(values) != 0:
        raise ValueError(f"word is not zero-sum: {text}")
    return values


def canonical_cycle(word: tuple[int, ...]) -> tuple[int, ...]:
    rotations = [word[index:] + word[:index] for index in range(len(word))]
    return min(rotations)


def internal_return_count(word: tuple[int, ...]) -> int:
    total = 0
    returns = 0
    for value in word[:-1]:
        total += value
        if total == 0:
            returns += 1
    return returns


def closure_multiplicity(rule: str, returns: int) -> int:
    if rule == "sterile":
        return 1
    if rule == "single_cut":
        return 1 + returns
    if rule == "all_cuts":
        return 2**returns
    raise ValueError(f"unknown closure rule: {rule}")


def build_modes() -> tuple[tuple[int, ...], ...]:
    modes = set()
    for period in range(2, MAX_PERIOD + 1, 2):
        positive = period // 2
        for positions in combinations(range(period), positive):
            positive_positions = set(positions)
            word = tuple(
                1 if index in positive_positions else -1
                for index in range(period)
            )
            modes.add(canonical_cycle(word))
    return tuple(sorted(modes, key=lambda item: (len(item), item)))


def raw_zero_sum_neighbors(word: tuple[int, ...]) -> set[tuple[int, ...]]:
    neighbors: set[tuple[int, ...]] = set()
    period = len(word)

    # A local swap does not change the number of + and - steps.
    for index in range(period):
        next_index = (index + 1) % period
        if word[index] == word[next_index]:
            continue
        changed = list(word)
        changed[index], changed[next_index] = changed[next_index], changed[index]
        neighbors.add(tuple(changed))

    # Delete one adjacent compensating pair.
    if period - 2 >= MIN_PERIOD:
        for index in range(period):
            next_index = (index + 1) % period
            if word[index] == word[next_index]:
                continue
            if index == period - 1:
                reduced = word[1:-1]
            else:
                reduced = word[:index] + word[index + 2 :]
            if len(reduced) >= MIN_PERIOD and sum(reduced) == 0:
                neighbors.add(reduced)

    # Insert one adjacent compensating pair.
    if period + 2 <= MAX_PERIOD:
        for index in range(period + 1):
            for pair in ((1, -1), (-1, 1)):
                expanded = word[:index] + pair + word[index:]
                neighbors.add(expanded)

    return neighbors


def mutation_kernel(
    modes: tuple[tuple[int, ...], ...],
) -> dict[int, tuple[int, ...]]:
    mode_index = {word: index for index, word in enumerate(modes)}
    kernel: dict[int, tuple[int, ...]] = {}
    for index, word in enumerate(modes):
        targets = set()
        for neighbor in raw_zero_sum_neighbors(word):
            canonical = canonical_cycle(neighbor)
            target = mode_index.get(canonical)
            if target is not None and target != index:
                targets.add(target)
        kernel[index] = tuple(sorted(targets))
    return kernel


@dataclass(frozen=True)
class Scenario:
    name: str
    rule: str
    mutation: float
    seed_word: str
    label: str


SCENARIOS = (
    Scenario(
        name="sterile_regular",
        rule="sterile",
        mutation=0.0,
        seed_word="+-+-+-",
        label="regular but sterile",
    ),
    Scenario(
        name="one_cut_clone",
        rule="single_cut",
        mutation=0.0,
        seed_word="+-+-+-",
        label="one-cut branching",
    ),
    Scenario(
        name="one_cut_evolution",
        rule="single_cut",
        mutation=MUTATION_PROBABILITY,
        seed_word="++--",
        label="one-cut with mutation",
    ),
    Scenario(
        name="all_cuts_evolution",
        rule="all_cuts",
        mutation=MUTATION_PROBABILITY,
        seed_word="++--",
        label="all-cut with mutation",
    ),
)


def mode_table(
    modes: tuple[tuple[int, ...], ...],
    rule: str,
) -> dict[int, dict[str, float]]:
    table = {}
    for index, word in enumerate(modes):
        returns = internal_return_count(word)
        multiplicity = closure_multiplicity(rule, returns)
        period = len(word)
        table[index] = {
            "period": period,
            "returns": returns,
            "return_density": returns / period if period else 0.0,
            "multiplicity": multiplicity,
            "lambda": math.log(multiplicity) / period if multiplicity > 1 else 0.0,
        }
    return table


def shannon(counts: np.ndarray) -> float:
    total = float(np.sum(counts))
    if total <= 0.0:
        return 0.0
    positive = counts[counts > 0.0] / total
    value = float(-np.sum(positive * np.log(positive)))
    return 0.0 if abs(value) < 1.0e-15 else value


def simulate(
    scenario: Scenario,
    modes: tuple[tuple[int, ...], ...],
    kernel: dict[int, tuple[int, ...]],
) -> dict[str, object]:
    table = mode_table(modes, scenario.rule)
    mode_index = {word: index for index, word in enumerate(modes)}
    seed = canonical_cycle(word_from_string(scenario.seed_word))
    if seed not in mode_index:
        raise ValueError(f"seed is not a represented mode: {scenario.seed_word}")

    cohorts: dict[int, np.ndarray] = {}
    for index, word in enumerate(modes):
        period = len(word)
        cohorts[index] = np.zeros(period, dtype=float)
    cohorts[mode_index[seed]][0] = 1.0

    snapshots = []
    rescale_steps = 0

    for step in range(STEPS + 1):
        population = np.array(
            [float(np.sum(cohorts[index])) for index in range(len(modes))],
            dtype=float,
        )
        total = float(np.sum(population))
        if total > 0.0:
            shares = population / total
            top_index = int(np.argmax(population))
            top_word = modes[top_index]
            lineage_count = int(np.count_nonzero(population > 1.0e-12 * total))
            mean_period = float(
                sum(shares[index] * table[index]["period"] for index in range(len(modes)))
            )
            mean_returns = float(
                sum(shares[index] * table[index]["returns"] for index in range(len(modes)))
            )
            mean_return_density = float(
                sum(
                    shares[index] * table[index]["return_density"]
                    for index in range(len(modes))
                )
            )
            mean_lambda = float(
                sum(shares[index] * table[index]["lambda"] for index in range(len(modes)))
            )
            top_fraction = float(population[top_index] / total)
            top_lambda = float(table[top_index]["lambda"])
        else:
            top_index = -1
            top_word = ()
            lineage_count = 0
            mean_period = 0.0
            mean_returns = 0.0
            mean_return_density = 0.0
            mean_lambda = 0.0
            top_fraction = 0.0
            top_lambda = 0.0

        snapshots.append(
            {
                "time": step,
                "population": total,
                "lineages": lineage_count,
                "diversity": shannon(population),
                "mean_period": mean_period,
                "mean_returns": mean_returns,
                "mean_return_density": mean_return_density,
                "mean_lambda": mean_lambda,
                "top_word": "".join("+" if value == 1 else "-" for value in top_word),
                "top_fraction": top_fraction,
                "top_lambda": top_lambda,
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

            multiplicity = float(table[index]["multiplicity"])
            offspring = finishers * multiplicity
            targets = kernel[index]
            if scenario.mutation <= 0.0 or not targets:
                next_cohorts[index][0] += offspring
            else:
                faithful = offspring * (1.0 - scenario.mutation)
                mutated = offspring * scenario.mutation / len(targets)
                next_cohorts[index][0] += faithful
                for target in targets:
                    next_cohorts[target][0] += mutated

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
        for index in order[:12]:
            if final_population[index] <= 1.0e-12 * total:
                continue
            word = modes[int(index)]
            top_entries.append(
                {
                    "word": "".join("+" if value == 1 else "-" for value in word),
                    "period": int(table[int(index)]["period"]),
                    "returns": int(table[int(index)]["returns"]),
                    "multiplicity": int(table[int(index)]["multiplicity"]),
                    "lambda": float(table[int(index)]["lambda"]),
                    "fraction": float(final_population[index] / total),
                }
            )

    return {
        "scenario": scenario.name,
        "label": scenario.label,
        "rule": scenario.rule,
        "mutation_probability": scenario.mutation,
        "seed_word": scenario.seed_word,
        "snapshots": snapshots,
        "rescale_steps": rescale_steps,
        "final_population": total,
        "final_lineages": int(np.count_nonzero(final_population > 0.0)),
        "final_diversity": shannon(final_population),
        "top_entries": top_entries,
    }


def summarize(runs: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    summary = {}
    for run in runs:
        final_snapshot = run["snapshots"][-1]
        summary[run["scenario"]] = {
            "label": run["label"],
            "rule": run["rule"],
            "mutation_probability": run["mutation_probability"],
            "seed_word": run["seed_word"],
            "final_population": run["final_population"],
            "final_lineages": run["final_lineages"],
            "final_diversity": run["final_diversity"],
            "mean_return_density": final_snapshot["mean_return_density"],
            "mean_lambda": final_snapshot["mean_lambda"],
            "top_word": final_snapshot["top_word"],
            "top_fraction": final_snapshot["top_fraction"],
            "top_lambda": final_snapshot["top_lambda"],
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
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.0))
    colors = {
        "sterile_regular": "#777777",
        "one_cut_clone": "#2B7A78",
        "one_cut_evolution": "#C77C28",
        "all_cuts_evolution": "#B44A5A",
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
            color=colors[run["scenario"]],
            linewidth=1.8,
            label=run["label"],
        )

    axes[0, 0].set_yscale("log")
    axes[0, 0].set_title("Inherited closed cycles")
    axes[0, 0].set_xlabel("internal time")
    axes[0, 0].set_ylabel("population (log scale)")
    axes[0, 0].grid(alpha=0.22)
    axes[0, 0].legend(frameon=False)

    evolution = next(
        run for run in runs if run["scenario"] == "one_cut_evolution"
    )
    times = [snapshot["time"] for snapshot in evolution["snapshots"]]
    axes[0, 1].plot(
        times,
        [snapshot["mean_return_density"] for snapshot in evolution["snapshots"]],
        color="#C77C28",
        linewidth=2.0,
        label="mean zero-return density",
    )
    axes[0, 1].plot(
        times,
        [snapshot["top_fraction"] for snapshot in evolution["snapshots"]],
        color="#2B7A78",
        linewidth=1.7,
        label="leading lineage share",
    )
    axes[0, 1].set_ylim(0.0, 1.02)
    axes[0, 1].set_title("Differentiation without a fertility landscape")
    axes[0, 1].set_xlabel("internal time")
    axes[0, 1].set_ylabel("fraction")
    axes[0, 1].grid(alpha=0.22)
    axes[0, 1].legend(frameon=False)

    rule = "single_cut"
    table = mode_table(modes, rule)
    periods = np.array([table[index]["period"] for index in range(len(modes))])
    returns = np.array([table[index]["returns"] for index in range(len(modes))])
    lambdas = np.array([table[index]["lambda"] for index in range(len(modes))])
    axes[1, 0].scatter(
        returns / periods,
        lambdas,
        c=periods,
        cmap="viridis",
        s=28,
        alpha=0.8,
    )
    axes[1, 0].set_title("Closure algebra sets the growth exponent")
    axes[1, 0].set_xlabel("return density r/T")
    axes[1, 0].set_ylabel("lambda = log(1+r)/T")
    axes[1, 0].grid(alpha=0.22)

    top = evolution["top_entries"][:10]
    labels = [
        f'{entry["word"]}\nT={entry["period"]} r={entry["returns"]}'
        for entry in top
    ]
    fractions = [entry["fraction"] for entry in top]
    y = np.arange(len(top))
    axes[1, 1].barh(y, fractions, color="#C77C28", alpha=0.9)
    axes[1, 1].set_yticks(y)
    axes[1, 1].set_yticklabels(labels, fontsize=7)
    axes[1, 1].invert_yaxis()
    axes[1, 1].set_xlabel("final share")
    axes[1, 1].set_title("Leading surviving zero-sum cycles")
    axes[1, 1].grid(axis="x", alpha=0.22)

    fig.suptitle("Zero-sum cycle evolution sandbox", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUT_PNG, dpi=180)
    plt.close(fig)


def main() -> None:
    modes = build_modes()
    kernel = mutation_kernel(modes)
    max_charge_error = max(abs(sum(word)) for word in modes)
    mutation_charge_error = 0
    for index, targets in kernel.items():
        source_charge = sum(modes[index])
        for target in targets:
            mutation_charge_error = max(
                mutation_charge_error,
                abs(source_charge - sum(modes[target])),
            )
    runs = [simulate(scenario, modes, kernel) for scenario in SCENARIOS]
    summary = summarize(runs)
    make_plot(runs, modes)

    payload = {
        "model": {
            "zero_constraint": "sum(charges)=0 for every cycle and every mutation",
            "period": "number of signed steps in one closed cycle",
            "internal_returns": "number of proper times the cycle returns to zero",
            "rules": {
                "sterile": "M=1",
                "single_cut": "M=1+r",
                "all_cuts": "M=2**r",
            },
            "lineage_rate": "lambda=log(M)/T",
            "mutation_kernel": "charge-preserving swaps, insertions, deletions",
            "mutation_probability": MUTATION_PROBABILITY,
            "steps": STEPS,
            "min_period": MIN_PERIOD,
            "max_period": MAX_PERIOD,
            "population_scale_limit": POPULATION_SCALE_LIMIT,
            "max_mode_charge_error": max_charge_error,
            "max_mutation_charge_error": mutation_charge_error,
            "note": (
                "Population rescaling is numerical only. The closure algebra "
                "and the initial seed are explicit inputs, not derived claims. "
                "M is a branch-program multiplicity; interpreting it as a "
                "physical copy number assumes an extra reproduction rule."
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
    print(f"\nplot: {OUT_PNG}")
    print(f"json: {OUT_JSON}")


if __name__ == "__main__":
    main()
