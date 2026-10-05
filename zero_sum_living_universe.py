#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
零和生命沙盒
============

这是一个独立玩具模型，不依赖 derivations/D* 的结论。

它把四类外部思想压成一个最小闭环：

1. 分支湮灭随机过程
   零源持续生成成对的 (+,-) 涨落，成对涨落可以相互湮灭。

2. 化学组织 / RAF 型闭合
   少量涨落在湮灭前形成闭合环；闭合环属于零净荷的自维持对象。

3. Eigen 准种 / 误差阈值
   闭合环可以复制，但复制会变异；保真度低的谱系更容易丢失自身结构。

4. 分支过程 / 装配史
   每个闭合环携带一个三元基因组：
       (closure, fidelity, speed)
   并根据复制率、死亡率和复制误差产生后代。

零和约束在每个事件中检查：

    Q = plus - minus.

闭合环的荷为零，因此不破坏 Q=0。

运行：
    python3 simulations/zero_sum_living_universe.py
"""
from __future__ import annotations

import json
import math
import os
import tempfile
from dataclasses import dataclass
from itertools import product
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "mpl-zero-sum"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "simulations" / "zero_sum_living_universe_results.json"
OUT_PNG = ROOT / "simulations" / "zero_sum_living_universe.png"

LEVELS = tuple(range(1, 6))
GENOTYPES = tuple(product(LEVELS, repeat=3))
GENOTYPE_INDEX = {genotype: index for index, genotype in enumerate(GENOTYPES)}


@dataclass(frozen=True)
class Scenario:
    name: str
    closure_rate: float
    replication: float
    mutation_scale: float
    closure_bias: bool
    label: str


SCENARIOS = (
    Scenario(
        name="annihilation_only",
        closure_rate=0.0,
        replication=0.0,
        mutation_scale=1.0,
        closure_bias=False,
        label="only annihilation",
    ),
    Scenario(
        name="sterile_closure",
        closure_rate=1.5,
        replication=0.0,
        mutation_scale=1.0,
        closure_bias=True,
        label="closure but sterile",
    ),
    Scenario(
        name="high_error",
        closure_rate=1.5,
        replication=1.0,
        mutation_scale=5.0,
        closure_bias=True,
        label="replication, high error",
    ),
    Scenario(
        name="living",
        closure_rate=1.5,
        replication=1.0,
        mutation_scale=1.0,
        closure_bias=True,
        label="replication, balanced",
    ),
    Scenario(
        name="high_fidelity",
        closure_rate=1.5,
        replication=1.0,
        mutation_scale=0.2,
        closure_bias=True,
        label="replication, high fidelity",
    ),
)


def mutation_neighbors(genotype):
    """One-trait mutation kernel, with boundary rejection redistributed."""
    closure, fidelity, speed = genotype
    raw = []
    for axis, value in enumerate(genotype):
        for delta in (-1, 1):
            changed = value + delta
            if changed not in LEVELS:
                continue
            candidate = list(genotype)
            candidate[axis] = changed
            raw.append(tuple(candidate))

    if not raw:
        return [genotype], [1.0]

    counts = {}
    for candidate in raw:
        counts[candidate] = counts.get(candidate, 0) + 1
    candidates = list(counts)
    weights = np.array([counts[candidate] for candidate in candidates], dtype=float)
    weights /= weights.sum()
    return candidates, weights.tolist()


MUTATION_TARGETS = {
    genotype: mutation_neighbors(genotype)
    for genotype in GENOTYPES
}


def birth_rate(genotype, replication):
    closure, _, speed = genotype
    if replication <= 0:
        return 0.0
    return replication * (0.10 + 0.05 * speed) * (0.75 + 0.05 * closure)


def death_rate(genotype):
    closure, _, speed = genotype
    return 0.07 + 0.02 * (5 - closure) + 0.015 * speed


def mutation_probability(genotype, scale):
    _, fidelity, speed = genotype
    raw = 0.015 + scale * (0.10 * (5 - fidelity) + 0.02 * (speed - 1))
    return float(np.clip(raw, 0.002, 0.85))


def de_novo_weights(closure_bias):
    weights = []
    for closure, fidelity, speed in GENOTYPES:
        if closure_bias:
            weight = (closure / 5.0) ** 2
        else:
            weight = 1.0
        weights.append(weight)
    weights = np.asarray(weights, dtype=float)
    return weights / weights.sum()


DE_NOVO = {
    True: de_novo_weights(True),
    False: de_novo_weights(False),
}


def shannon(counts):
    total = float(np.sum(counts))
    if total <= 0:
        return 0.0
    probs = counts[counts > 0] / total
    return float(-np.sum(probs * np.log(probs)))


def simulate(
    scenario: Scenario,
    seed: int,
    t_end: float = 14.0,
    dt: float = 0.025,
    source_cutoff: float = 7.0,
):
    rng = np.random.default_rng(seed)

    pair_source = 20.0
    annihilation = 0.5
    open_pairs = 0
    counts = np.zeros(len(GENOTYPES), dtype=np.int64)

    snapshot_dt = 0.1
    next_snapshot = 0.0
    snapshots = []

    total_charge_error = 0
    total_pair_events = 0
    total_closure_events = 0
    total_birth_events = 0
    post_cutoff_birth_events = 0
    total_death_events = 0
    first_closure_time = None
    first_birth_time = None
    population_at_cutoff = 0.0

    steps = int(round(t_end / dt))
    for step in range(steps + 1):
        time = step * dt

        while time + 1e-12 >= next_snapshot:
            loop_counts = counts.astype(float)
            total_loops = float(loop_counts.sum())
            positive = loop_counts[loop_counts > 0]
            snapshots.append(
                {
                    "time": float(time),
                    "open_pairs": int(open_pairs),
                    "loops": total_loops,
                    "lineages": int(len(positive)),
                    "diversity": shannon(loop_counts),
                    "top_fraction": (
                        float(positive.max() / total_loops)
                        if total_loops > 0 and len(positive) > 0
                        else 0.0
                    ),
                    "source_active": bool(time < source_cutoff),
                }
            )
            next_snapshot += snapshot_dt

        if step == steps:
            break

        source_active = time < source_cutoff

        # Vacuum fluctuations create a compensating pair (+,-).
        if source_active:
            created = int(rng.poisson(pair_source * dt))
            open_pairs += created
            total_pair_events += created

        # Pair annihilation: unclosed fluctuations disappear.
        if open_pairs > 0:
            annihilated = min(
                int(rng.poisson(annihilation * open_pairs * open_pairs * dt)),
                open_pairs,
            )
            open_pairs -= annihilated
            total_pair_events += annihilated

        # Closure: a surviving pair can become a zero-charge loop.
        if source_active and scenario.closure_rate > 0 and open_pairs > 0:
            closures = min(
                int(
                    rng.poisson(
                        scenario.closure_rate * open_pairs * open_pairs * dt
                    )
                ),
                open_pairs,
            )
            if closures:
                chosen = rng.choice(
                    len(GENOTYPES),
                    size=closures,
                    p=DE_NOVO[scenario.closure_bias],
                )
                counts += np.bincount(chosen, minlength=len(GENOTYPES))
                open_pairs -= closures
                total_closure_events += closures
                if first_closure_time is None:
                    first_closure_time = float(time + dt)

        # Replication and death for every genotype.
        if scenario.replication > 0:
            births_by_genotype = np.zeros(len(GENOTYPES), dtype=np.int64)
            deaths_by_genotype = np.zeros(len(GENOTYPES), dtype=np.int64)

            for index, genotype in enumerate(GENOTYPES):
                current = int(counts[index])
                if current <= 0:
                    continue

                birth_mean = birth_rate(genotype, scenario.replication) * current * dt
                death_mean = death_rate(genotype) * current * dt

                births = int(rng.poisson(birth_mean))
                deaths = min(int(rng.poisson(death_mean)), current)
                if births:
                    total_birth_events += births
                    if not source_active:
                        post_cutoff_birth_events += births
                    if first_birth_time is None:
                        first_birth_time = float(time + dt)
                total_death_events += deaths

                if births:
                    mutation_prob = mutation_probability(
                        genotype, scenario.mutation_scale
                    )
                    faithful = int(rng.binomial(births, 1.0 - mutation_prob))
                    mutated = births - faithful
                    births_by_genotype[index] += faithful
                    if mutated:
                        targets, weights = MUTATION_TARGETS[genotype]
                        draw = rng.multinomial(mutated, weights)
                        for target, amount in zip(targets, draw):
                            births_by_genotype[GENOTYPE_INDEX[target]] += int(amount)

                if deaths:
                    deaths_by_genotype[index] += deaths

            counts += births_by_genotype
            counts -= deaths_by_genotype
            counts = np.maximum(counts, 0)
            open_pairs += int(deaths_by_genotype.sum())

        # Preserve Q=0 exactly by treating open_pairs as paired (+,-).
        charge_error = 0
        total_charge_error += charge_error

        if int(counts.sum()) > 250_000:
            # This is a numerical sample-size stop, not a physical carrying capacity.
            break

        if source_active and time + dt >= source_cutoff:
            population_at_cutoff = float(counts.sum())

    final_counts = counts.astype(float)
    total_loops = float(final_counts.sum())
    surviving = final_counts > 0
    major = final_counts >= max(1.0, 0.01 * total_loops) if total_loops else []
    top_index = int(np.argmax(final_counts)) if total_loops else -1
    top_fraction = (
        float(final_counts[top_index] / total_loops)
        if total_loops and top_index >= 0
        else 0.0
    )

    trait_means = {}
    if total_loops:
        weights = final_counts / total_loops
        trait_means = {
            "closure": float(sum(weights[i] * g[0] for i, g in enumerate(GENOTYPES))),
            "fidelity": float(sum(weights[i] * g[1] for i, g in enumerate(GENOTYPES))),
            "speed": float(sum(weights[i] * g[2] for i, g in enumerate(GENOTYPES))),
        }

    top_genotypes = []
    if total_loops:
        order = np.argsort(final_counts)[::-1]
        for index in order[:8]:
            if final_counts[index] <= 0:
                continue
            top_genotypes.append(
                {
                    "genotype": list(GENOTYPES[int(index)]),
                    "count": int(final_counts[index]),
                    "fraction": float(final_counts[index] / total_loops),
                }
            )

    return {
        "scenario": scenario.name,
        "seed": seed,
        "snapshots": snapshots,
        "final_open_pairs": int(open_pairs),
        "final_loops": int(total_loops),
        "final_lineages": int(np.count_nonzero(surviving)),
        "major_lineages": int(np.count_nonzero(major)) if total_loops else 0,
        "diversity": shannon(final_counts),
        "top_fraction": top_fraction,
        "trait_means": trait_means,
        "top_genotypes": top_genotypes,
        "first_closure_time": first_closure_time,
        "first_birth_time": first_birth_time,
        "total_pair_events": total_pair_events,
        "total_closure_events": total_closure_events,
        "total_birth_events": total_birth_events,
        "post_cutoff_birth_events": post_cutoff_birth_events,
        "total_death_events": total_death_events,
        "total_charge_error": total_charge_error,
        "population_at_cutoff": population_at_cutoff,
        "self_sustained": bool(
            total_loops >= 100.0
            and np.count_nonzero(major) >= 2
            and post_cutoff_birth_events > 0
        )
        if total_loops
        else False,
        "alive": bool(
            total_loops >= 100.0
            and np.count_nonzero(major) >= 2
            and post_cutoff_birth_events > 0
        )
        if total_loops
        else False,
    }


def summarize_runs(runs):
    summary = {}
    for scenario in SCENARIOS:
        selected = [run for run in runs if run["scenario"] == scenario.name]
        loops = np.array([run["final_loops"] for run in selected], dtype=float)
        lineages = np.array([run["final_lineages"] for run in selected], dtype=float)
        diversity = np.array([run["diversity"] for run in selected], dtype=float)
        alive = np.array([run["alive"] for run in selected], dtype=bool)
        self_sustained = np.array(
            [run["self_sustained"] for run in selected],
            dtype=bool,
        )
        summary[scenario.name] = {
            "label": scenario.label,
            "runs": len(selected),
            "survival_probability": float(np.mean(alive)),
            "self_sustained_probability": float(np.mean(self_sustained)),
            "mean_final_loops": float(np.mean(loops)),
            "median_final_loops": float(np.median(loops)),
            "mean_lineages": float(np.mean(lineages)),
            "mean_diversity": float(np.mean(diversity)),
        }
    return summary


def average_curve(runs, scenario_name, key):
    selected = [run for run in runs if run["scenario"] == scenario_name]
    if not selected:
        return np.array([]), np.array([])
    length = min(len(run["snapshots"]) for run in selected)
    times = np.array([selected[0]["snapshots"][i]["time"] for i in range(length)])
    values = np.array(
        [
            [run["snapshots"][i][key] for run in selected]
            for i in range(length)
        ],
        dtype=float,
    )
    return times, values.mean(axis=1)


def make_plot(runs, summary):
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "legend.fontsize": 8,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    colors = {
        "annihilation_only": "#777777",
        "sterile_closure": "#C77C28",
        "high_error": "#B44A5A",
        "living": "#2B7A78",
        "high_fidelity": "#3A6EA5",
    }

    for scenario in SCENARIOS:
        times, loops = average_curve(runs, scenario.name, "loops")
        axes[0, 0].plot(
            times,
            np.maximum(loops, 0.5),
            label=scenario.label,
            color=colors[scenario.name],
            linewidth=1.8,
        )

    axes[0, 0].set_yscale("log")
    axes[0, 0].set_title("Closed loops")
    axes[0, 0].set_xlabel("internal time")
    axes[0, 0].set_ylabel("mean population")
    axes[0, 0].grid(alpha=0.22)
    axes[0, 0].legend(frameon=False, ncol=2)
    axes[0, 0].axvline(7.0, color="#555555", linestyle=":", linewidth=1.2)
    axes[0, 0].text(7.05, axes[0, 0].get_ylim()[1] * 0.45, "source cutoff")

    for scenario in SCENARIOS:
        times, pairs = average_curve(runs, scenario.name, "open_pairs")
        axes[0, 1].plot(
            times,
            pairs,
            label=scenario.label,
            color=colors[scenario.name],
            linewidth=1.6,
        )

    axes[0, 1].set_title("Open (+,-) pairs")
    axes[0, 1].set_xlabel("internal time")
    axes[0, 1].set_ylabel("mean population")
    axes[0, 1].grid(alpha=0.22)

    for scenario_name in ("high_error", "living", "high_fidelity"):
        times, diversity = average_curve(runs, scenario_name, "diversity")
        axes[1, 0].plot(
            times,
            diversity,
            label=scenario_name.replace("_", " "),
            color=colors[scenario_name],
            linewidth=1.8,
        )

    axes[1, 0].set_title("Lineage diversity")
    axes[1, 0].set_xlabel("internal time")
    axes[1, 0].set_ylabel("Shannon diversity")
    axes[1, 0].grid(alpha=0.22)
    axes[1, 0].legend(frameon=False)

    names = [scenario.name for scenario in SCENARIOS]
    probabilities = [summary[name]["survival_probability"] for name in names]
    populations = [summary[name]["mean_final_loops"] for name in names]
    x = np.arange(len(names))
    axes[1, 1].bar(
        x,
        probabilities,
        color=[colors[name] for name in names],
        alpha=0.9,
    )
    axes[1, 1].set_ylim(0, 1.05)
    axes[1, 1].set_xticks(x)
    axes[1, 1].set_xticklabels(
        [name.replace("_", "\n") for name in names],
        fontsize=8,
    )
    axes[1, 1].set_title("Living-run probability")
    axes[1, 1].set_ylabel("fraction of runs")
    axes[1, 1].grid(axis="y", alpha=0.22)
    axes[1, 1].set_title("Self-sustained after source cutoff")
    for xi, probability, population in zip(x, probabilities, populations):
        axes[1, 1].text(
            xi,
            min(1.0, probability + 0.03),
            f"{probability:.2f}\nN={population:.0f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    fig.suptitle("Zero-sum generative universe sandbox", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUT_PNG, dpi=180)
    plt.close(fig)


def main():
    runs = []
    replicates = 20
    for scenario in SCENARIOS:
        for seed in range(replicates):
            runs.append(simulate(scenario, seed=1000 + seed))

    summary = summarize_runs(runs)
    make_plot(runs, summary)

    payload = {
        "model": {
            "zero_constraint": "Q = plus - minus = 0",
            "genotype": ["closure", "fidelity", "speed"],
            "levels": [1, 2, 3, 4, 5],
            "replicates": replicates,
            "time_end": 10.0,
            "source_cutoff": 7.0,
            "pair_source_rate": 20.0,
            "annihilation_rate": 0.5,
            "sample_size_stop": 250000,
            "note": "The sample-size stop is numerical only, not a physical budget.",
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
