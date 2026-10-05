#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Zero-sum closure versus reproduction audit
==========================================

The earlier cycle-evolution sandbox used multipliers such as

    M = 1 + r
    M = 2**r

where r was the number of internal returns to zero. Those numbers counted
possible descriptions or branch programs. They did not automatically count
new physical loops.

This audit separates four objects that had been conflated:

    closure
    persistence
    splitting
    copying

All rules are deterministic and have no channel weights:

    persist       w -> [w]
    copy          w -> [w, w]
    split_parent  w -> [w] plus one split for every internal zero-return cut
    split_only    w -> one split for every internal zero-return cut
    split_all     w -> all decompositions at subsets of zero-return cuts

The exact transition matrix is built over all cyclic zero-sum words up to
period 12. There are no probabilities and no fitted fertility parameter.

Run:
    python3 simulations/zero_sum_reproduction_audit.py
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

try:
    from simulations.zero_sum_cycle_evolution import canonical_cycle
except ModuleNotFoundError:
    from zero_sum_cycle_evolution import canonical_cycle


ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "simulations" / "zero_sum_reproduction_audit_results.json"
OUT_PNG = ROOT / "simulations" / "zero_sum_reproduction_audit.png"

MAX_PERIOD = 12
GENERATIONS = 18
SEEDS = (
    "++--",
    "+-+-",
    "+-+-+-",
    "+-+-+-+-+-+-",
)


@dataclass(frozen=True)
class Rule:
    name: str
    label: str


RULES = (
    Rule("persist", "persist only"),
    Rule("copy", "copy only"),
    Rule("split_parent", "persist + one-cut splits"),
    Rule("split_only", "splits without parent"),
    Rule("split_all", "all zero-return decompositions"),
)


def word_from_string(text: str) -> tuple[int, ...]:
    value = tuple(1 if char == "+" else -1 for char in text)
    if sum(value) != 0:
        raise ValueError(f"word is not zero-sum: {text}")
    return value


def word_to_string(word: tuple[int, ...]) -> str:
    return "".join("+" if value == 1 else "-" for value in word)


def all_modes(max_period: int = MAX_PERIOD) -> tuple[tuple[int, ...], ...]:
    modes = set()
    for period in range(2, max_period + 1, 2):
        for positions in combinations(range(period), period // 2):
            positive = set(positions)
            word = tuple(
                1 if index in positive else -1
                for index in range(period)
            )
            modes.add(canonical_cycle(word))
    return tuple(sorted(modes, key=lambda item: (len(item), item)))


def internal_cuts(word: tuple[int, ...]) -> tuple[int, ...]:
    total = 0
    cuts = []
    for index, value in enumerate(word[:-1], start=1):
        total += value
        if total == 0:
            cuts.append(index)
    return tuple(cuts)


def split_at(
    word: tuple[int, ...],
    cut: int,
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    return (
        canonical_cycle(word[:cut]),
        canonical_cycle(word[cut:]),
    )


def child_branches(
    word: tuple[int, ...],
    rule: str,
) -> list[tuple[tuple[int, ...], ...]]:
    cuts = internal_cuts(word)

    if rule == "persist":
        return [(word,)]

    if rule == "copy":
        return [(word, word)]

    if rule == "split_parent":
        branches = [(word,)]
        branches.extend(split_at(word, cut) for cut in cuts)
        return branches

    if rule == "split_only":
        return [split_at(word, cut) for cut in cuts]

    if rule == "split_all":
        branches = []
        cut_count = len(cuts)
        for mask in range(1 << cut_count):
            positions = [0]
            positions.extend(
                cuts[index]
                for index in range(cut_count)
                if (mask >> index) & 1
            )
            positions.append(len(word))
            branches.append(
                tuple(
                    canonical_cycle(word[positions[index] : positions[index + 1]])
                    for index in range(len(positions) - 1)
                )
            )
        return branches

    raise ValueError(f"unknown rule: {rule}")


def transition_matrix(
    modes: tuple[tuple[int, ...], ...],
    rule: str,
) -> np.ndarray:
    mode_index = {word: index for index, word in enumerate(modes)}
    matrix = np.zeros((len(modes), len(modes)), dtype=np.int64)
    for parent_index, word in enumerate(modes):
        for branch in child_branches(word, rule):
            for child in branch:
                matrix[mode_index[child], parent_index] += 1
    return matrix


def sequence(
    matrix: np.ndarray,
    modes: tuple[tuple[int, ...], ...],
    seed_word: str,
    generations: int,
) -> list[int]:
    mode_index = {word: index for index, word in enumerate(modes)}
    seed = canonical_cycle(word_from_string(seed_word))
    if seed not in mode_index:
        raise ValueError(f"seed is not in finite mode set: {seed_word}")

    population = np.zeros(len(modes), dtype=np.int64)
    population[mode_index[seed]] = 1
    values = []
    for _ in range(generations + 1):
        values.append(int(np.sum(population)))
        population = matrix @ population
    return values


def nilpotent_index(matrix: np.ndarray) -> int | None:
    if not np.any(matrix):
        return 1
    power = np.eye(matrix.shape[0], dtype=np.int64)
    for index in range(1, matrix.shape[0] + 1):
        power = power @ matrix
        if not np.any(power):
            return index
    return None


def classify_sequence(values: list[int]) -> str:
    if values[-1] == 0:
        return "finite extinction"
    if all(value == values[0] for value in values):
        return "constant"
    if len(values) >= 3 and values[-2] > 0:
        ratio = values[-1] / values[-2]
        if ratio > 1.5:
            return "exponential"
    return "polynomial/subexponential"


def matrix_audit(matrix: np.ndarray) -> dict[str, object]:
    float_matrix = matrix.astype(float)
    eigenvalues = np.linalg.eigvals(float_matrix)
    spectral_radius = float(np.max(np.abs(eigenvalues)))

    result: dict[str, object] = {
        "spectral_radius": spectral_radius,
        "zero_eigenvalues": int(np.count_nonzero(np.abs(eigenvalues) < 1.0e-10)),
    }
    if math.isclose(spectral_radius, 1.0, abs_tol=1.0e-10):
        unipotent = matrix - np.eye(matrix.shape[0], dtype=np.int64)
        result["unipotent_nilpotent_index"] = nilpotent_index(unipotent)
    elif math.isclose(spectral_radius, 0.0, abs_tol=1.0e-10):
        result["nilpotent_index"] = nilpotent_index(matrix)
    return result


def multiplicity_table(modes: tuple[tuple[int, ...], ...]) -> list[dict[str, object]]:
    rows = []
    seen_returns = set()
    for word in sorted(modes, key=lambda item: (len(item), internal_cuts(item))):
        period = len(word)
        returns = len(internal_cuts(word))
        key = (period, returns)
        if key in seen_returns:
            continue
        seen_returns.add(key)

        split_parent = child_branches(word, "split_parent")
        split_only = child_branches(word, "split_only")
        split_all = child_branches(word, "split_all")
        rows.append(
            {
                "word": word_to_string(word),
                "period": period,
                "returns": returns,
                "split_parent_paths": len(split_parent),
                "split_parent_objects": sum(len(branch) for branch in split_parent),
                "split_only_paths": len(split_only),
                "split_only_objects": sum(len(branch) for branch in split_only),
                "split_all_paths": len(split_all),
                "split_all_objects": sum(len(branch) for branch in split_all),
            }
        )
    return rows


def make_plot(
    audits: dict[str, dict[str, object]],
    sequences: dict[str, dict[str, list[int]]],
    multiplicity: list[dict[str, object]],
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
    seed = "+-+-+-+-+-+-"
    colors = {
        "persist": "#777777",
        "copy": "#2B7A78",
        "split_parent": "#C77C28",
        "split_only": "#B44A5A",
        "split_all": "#3A6EA5",
    }
    generations = np.arange(GENERATIONS + 1)

    for rule in RULES:
        values = np.maximum(sequences[rule.name][seed], 0.5)
        axes[0, 0].plot(
            generations,
            values,
            marker="o",
            markersize=3,
            linewidth=1.7,
            color=colors[rule.name],
            label=rule.label,
        )
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_title("Actual loop population after each cycle")
    axes[0, 0].set_xlabel("cycle generation")
    axes[0, 0].set_ylabel("number of loop objects")
    axes[0, 0].grid(alpha=0.22)
    axes[0, 0].legend(frameon=False, fontsize=7)

    for rule in RULES:
        values = audits[rule.name]
        axes[0, 1].bar(
            rule.label,
            values["spectral_radius"],
            color=colors[rule.name],
            alpha=0.9,
        )
    axes[0, 1].set_title("Spectral radius alone does not distinguish growth")
    axes[0, 1].set_ylabel("rho(A)")
    axes[0, 1].tick_params(axis="x", rotation=25, labelsize=7)
    axes[0, 1].grid(axis="y", alpha=0.22)

    returns = np.array([row["returns"] for row in multiplicity], dtype=float)
    parent_objects = np.array(
        [row["split_parent_objects"] for row in multiplicity],
        dtype=float,
    )
    all_paths = np.array(
        [row["split_all_paths"] for row in multiplicity],
        dtype=float,
    )
    axes[1, 0].scatter(
        returns,
        parent_objects,
        label="split with parent: child objects",
        color="#C77C28",
        s=30,
    )
    axes[1, 0].scatter(
        returns,
        all_paths,
        label="all decompositions: branch programs",
        color="#3A6EA5",
        s=30,
    )
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_title("Internal zero returns count branches, not copies")
    axes[1, 0].set_xlabel("r = internal zero returns")
    axes[1, 0].set_ylabel("count")
    axes[1, 0].grid(alpha=0.22)
    axes[1, 0].legend(frameon=False, fontsize=7)

    labels = []
    degrees = []
    for rule in RULES:
        audit = audits[rule.name]
        index = audit.get("unipotent_nilpotent_index")
        if index is None:
            index = audit.get("nilpotent_index", 0)
        labels.append(rule.label)
        degrees.append(index)
    positions = np.arange(len(labels))
    axes[1, 1].bar(positions, degrees, color="#6B5B95", alpha=0.9)
    axes[1, 1].set_xticks(positions)
    axes[1, 1].set_xticklabels(labels, rotation=25, ha="right", fontsize=7)
    axes[1, 1].set_title("Jordan/nilpotent depth")
    axes[1, 1].set_ylabel("smallest k for the relevant power")
    axes[1, 1].grid(axis="y", alpha=0.22)

    fig.suptitle("Zero-sum reproduction formula audit", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUT_PNG, dpi=180)
    plt.close(fig)


def main() -> None:
    modes = all_modes()
    audits = {}
    sequences = {}
    for rule in RULES:
        matrix = transition_matrix(modes, rule.name)
        audits[rule.name] = matrix_audit(matrix)
        sequences[rule.name] = {
            seed: sequence(matrix, modes, seed, GENERATIONS)
            for seed in SEEDS
        }

    for rule in RULES:
        for seed, values in sequences[rule.name].items():
            audits[rule.name].setdefault("sequence_classes", {})[seed] = (
                classify_sequence(values)
            )

    multiplicity = multiplicity_table(modes)
    make_plot(audits, sequences, multiplicity)

    payload = {
        "model": {
            "zero_constraint": "each cyclic word has total charge zero",
            "modes": "all cyclic balanced words up to period 12",
            "weights": "none",
            "probabilities": "none",
            "rules": {rule.name: rule.label for rule in RULES},
            "generations": GENERATIONS,
        },
        "matrix_audit": audits,
        "multiplicity_table": multiplicity,
        "sequences": sequences,
        "interpretation": {
            "copying": (
                "among the five tested rules, only the explicit copy rule "
                "gives exponential growth; multi-type autocatalytic cycles "
                "remain a separate possibility"
            ),
            "closure": (
                "closure alone gives persistence, not reproduction"
            ),
            "splitting": (
                "splitting into non-reproducing subloops gives at most "
                "polynomial accumulation; if the parent is removed, the "
                "process can die in finite generations"
            ),
            "formula_warning": (
                "M=1+r and M=2**r are branch-program counts. They are not "
                "automatically physical offspring counts and do not by "
                "themselves imply exponential lineage growth"
            ),
        },
    }
    OUT_JSON.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(payload["interpretation"], ensure_ascii=False, indent=2))
    print("\nmatrix audit")
    for rule in RULES:
        audit = audits[rule.name]
        print(
            f"{rule.name:>12}  rho={audit['spectral_radius']:.3f}  "
            f"zero_eigs={audit['zero_eigenvalues']:>3}  "
            f"unipotent_depth={audit.get('unipotent_nilpotent_index')}  "
            f"nilpotent_depth={audit.get('nilpotent_index')}"
        )
    print(f"\nplot: {OUT_PNG}")
    print(f"json: {OUT_JSON}")


if __name__ == "__main__":
    main()
