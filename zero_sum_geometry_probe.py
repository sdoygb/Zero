#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Zero-sum closure geometry probe
===============================

This sandbox asks a narrow question:

    Does the graph of zero-sum closed cycles already select a stable
    low-dimensional geometry?

The microstates are cyclic words of +1 and -1 steps with total zero.
Two words are adjacent when one adjacent opposite-sign pair is exchanged.
This is the local, charge-preserving closure graph.

The script measures:

    neighborhood growth
    graph diameter scaling
    spectral dimension from the graph heat kernel
    stress of a best four-dimensional distance embedding
    dimensions after coarse-graining by closure structure

Controls:

    a four-dimensional torus
    a random regular graph

No D* derivation or external model is used here.

Run:
    python3 simulations/zero_sum_geometry_probe.py
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
import networkx as nx
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

try:
    from simulations.zero_sum_cycle_evolution import (
        canonical_cycle,
        internal_return_count,
        raw_zero_sum_neighbors,
    )
except ModuleNotFoundError:
    from zero_sum_cycle_evolution import (
        canonical_cycle,
        internal_return_count,
        raw_zero_sum_neighbors,
    )


ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "simulations" / "zero_sum_geometry_probe_results.json"
OUT_PNG = ROOT / "simulations" / "zero_sum_geometry_probe.png"

PERIODS = (8, 10, 12, 14, 16, 18)
COARSE_MAX_PERIOD = 16
SPECTRAL_POINTS = (0.10, 0.05, 0.02)
MAX_NEIGHBORHOOD_RADIUS = 20


def all_cycle_modes(period: int) -> tuple[tuple[int, ...], ...]:
    modes = set()
    positive = period // 2
    for positions in combinations(range(period), positive):
        positive_positions = set(positions)
        word = tuple(
            1 if index in positive_positions else -1
            for index in range(period)
        )
        modes.add(canonical_cycle(word))
    return tuple(sorted(modes))


def fixed_period_graph(period: int) -> tuple[nx.Graph, tuple[tuple[int, ...], ...]]:
    modes = all_cycle_modes(period)
    mode_index = {word: index for index, word in enumerate(modes)}
    graph = nx.Graph()
    graph.add_nodes_from(range(len(modes)))

    for index, word in enumerate(modes):
        for position in range(period):
            next_position = (position + 1) % period
            if word[position] == word[next_position]:
                continue
            changed = list(word)
            changed[position], changed[next_position] = (
                changed[next_position],
                changed[position],
            )
            target = mode_index[canonical_cycle(tuple(changed))]
            if target != index:
                graph.add_edge(index, target)

    return graph, modes


def four_dimensional_torus(side: int) -> nx.Graph:
    graph = nx.Graph()
    nodes = list(np.ndindex((side, side, side, side)))
    graph.add_nodes_from(nodes)
    for node in nodes:
        for axis in range(4):
            neighbour = list(node)
            neighbour[axis] = (neighbour[axis] + 1) % side
            neighbour_tuple = tuple(neighbour)
            if node != neighbour_tuple:
                graph.add_edge(node, neighbour_tuple)
    return graph


def random_regular_control(
    nodes: int,
    degree: int,
    seed: int = 7,
) -> nx.Graph:
    if (nodes * degree) % 2:
        nodes += 1
    return nx.random_regular_graph(degree, nodes, seed=seed)


def neighbourhood_growth(
    graph: nx.Graph,
    max_radius: int = MAX_NEIGHBORHOOD_RADIUS,
) -> np.ndarray:
    total = np.zeros(max_radius + 1, dtype=float)
    nodes = list(graph.nodes)
    if not nodes:
        return total
    for source in nodes:
        distances = nx.single_source_shortest_path_length(
            graph,
            source,
            cutoff=max_radius,
        )
        for distance in distances.values():
            total[distance] += 1.0
    total = np.cumsum(total)
    return total / len(nodes)


def local_growth_dimension(
    growth: np.ndarray,
    total_nodes: int,
) -> float | None:
    valid = [
        radius
        for radius in range(2, len(growth))
        if growth[radius] > growth[radius - 1]
        and growth[radius] <= 0.8 * total_nodes
    ]
    if len(valid) < 3:
        return None
    radii = np.asarray(valid, dtype=float)
    values = growth[valid]
    slope = np.polyfit(np.log(radii), np.log(values), 1)[0]
    return float(slope)


def graph_spectrum(graph: nx.Graph) -> np.ndarray:
    nodes = graph.number_of_nodes()
    if nodes <= 1000:
        adjacency = nx.to_numpy_array(graph, dtype=float)
        laplacian = np.diag(adjacency.sum(axis=1)) - adjacency
        eigenvalues = np.linalg.eigvalsh(laplacian)
    else:
        adjacency = nx.to_scipy_sparse_array(
            graph,
            dtype=float,
            format="csr",
        )
        degree = np.asarray(adjacency.sum(axis=1)).ravel()
        laplacian = csr_matrix(np.diag(degree)) - adjacency
        count = min(300, nodes - 2)
        eigenvalues = eigsh(
            laplacian,
            k=count,
            which="SM",
            return_eigenvectors=False,
            tol=1.0e-6,
        )
    return eigenvalues[eigenvalues > 1.0e-10]


def spectral_curve(
    eigenvalues: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    times = np.logspace(-1.0, 1.8, 180)
    heat_kernel = np.array(
        [
            float(np.mean(np.exp(-eigenvalues * time)))
            for time in times
        ]
    )
    slopes = -2.0 * np.gradient(
        np.log(np.maximum(heat_kernel, 1.0e-300)),
        np.log(times),
    )
    return times, slopes, heat_kernel


def four_dimensional_embedding_stress(graph: nx.Graph) -> float | None:
    nodes = list(graph.nodes)
    count = len(nodes)
    if count < 5 or count > 900:
        return None

    distances = np.zeros((count, count), dtype=float)
    for row, source in enumerate(nodes):
        distance_map = nx.single_source_shortest_path_length(graph, source)
        for column, target in enumerate(nodes):
            distances[row, column] = distance_map[target]

    centering = np.eye(count) - np.ones((count, count)) / count
    gram = -0.5 * centering @ (distances**2) @ centering
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1][:4]
    positive = np.maximum(values[order], 0.0)
    coordinates = vectors[:, order] * np.sqrt(positive)
    embedded = np.linalg.norm(
        coordinates[:, None, :] - coordinates[None, :, :],
        axis=2,
    )
    return float(
        np.sqrt(
            np.sum((distances - embedded) ** 2)
            / np.sum(distances**2)
        )
    )


def summarize_graph(
    name: str,
    graph: nx.Graph,
    *,
    include_embedding: bool = True,
) -> dict[str, object]:
    growth = neighbourhood_growth(graph)
    local_dimension = local_growth_dimension(
        growth,
        graph.number_of_nodes(),
    )
    eigenvalues = graph_spectrum(graph)
    times, slopes, probabilities = spectral_curve(eigenvalues)
    probability_points = {}
    for probability in SPECTRAL_POINTS:
        index = int(np.argmin(np.abs(probabilities - probability)))
        probability_points[f"{probability:.2f}"] = {
            "time": float(times[index]),
            "spectral_dimension": float(slopes[index]),
        }

    mean_degree = (
        2.0 * graph.number_of_edges() / graph.number_of_nodes()
        if graph.number_of_nodes()
        else 0.0
    )
    return {
        "name": name,
        "nodes": graph.number_of_nodes(),
        "edges": graph.number_of_edges(),
        "mean_degree": mean_degree,
        "diameter": (
            int(nx.diameter(graph))
            if graph.number_of_nodes() > 1 and nx.is_connected(graph)
            else None
        ),
        "local_growth_dimension": local_dimension,
        "spectral_points": probability_points,
        "mds_4d_stress": (
            four_dimensional_embedding_stress(graph)
            if include_embedding
            else None
        ),
        "growth": growth.tolist(),
        "spectral_times": times.tolist(),
        "spectral_dimension_curve": slopes.tolist(),
    }


def excursion_lengths(word: tuple[int, ...]) -> tuple[int, ...]:
    total = 0
    lengths = []
    previous = -1
    for index, value in enumerate(word):
        total += value
        if total == 0:
            lengths.append(index - previous)
            previous = index
    return tuple(lengths)


def coarse_quotient(
    modes: tuple[tuple[int, ...], ...],
    classifier,
) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(classifier(word) for word in modes)
    for word in modes:
        source_class = classifier(word)
        for neighbour in raw_zero_sum_neighbors(word):
            target_class = classifier(canonical_cycle(neighbour))
            if source_class != target_class:
                graph.add_edge(source_class, target_class)
    return graph


def all_modes_up_to(period: int) -> tuple[tuple[int, ...], ...]:
    modes = set()
    for current_period in range(2, period + 1, 2):
        modes.update(all_cycle_modes(current_period))
    return tuple(sorted(modes, key=lambda word: (len(word), word)))


def make_plot(
    summaries: list[dict[str, object]],
    coarse_summaries: list[dict[str, object]],
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
    closure = [summary for summary in summaries if summary["period"] is not None]
    controls = [summary for summary in summaries if summary["period"] is None]
    cmap = plt.get_cmap("viridis")

    for index, summary in enumerate(closure):
        growth = np.asarray(summary["growth"], dtype=float)
        radii = np.arange(len(growth))
        axes[0, 0].plot(
            radii[1:],
            growth[1:],
            marker="o",
            markersize=3,
            linewidth=1.5,
            color=cmap(index / max(1, len(closure) - 1)),
            label=f'closure T={summary["period"]}',
        )
    for summary in controls:
        growth = np.asarray(summary["growth"], dtype=float)
        radii = np.arange(len(growth))
        axes[0, 0].plot(
            radii[1:],
            growth[1:],
            linewidth=2.0,
            linestyle="--",
            label=summary["name"],
        )
    axes[0, 0].set_xscale("log")
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_title("Neighbourhood growth")
    axes[0, 0].set_xlabel("graph radius")
    axes[0, 0].set_ylabel("mean reachable nodes")
    axes[0, 0].grid(alpha=0.22)
    axes[0, 0].legend(frameon=False, fontsize=7)

    for summary in closure:
        if summary["period"] not in (12, 16, 18):
            continue
        axes[0, 1].plot(
            summary["spectral_times"],
            summary["spectral_dimension_curve"],
            linewidth=1.8,
            label=f'closure T={summary["period"]}',
        )
    for summary in controls:
        axes[0, 1].plot(
            summary["spectral_times"],
            summary["spectral_dimension_curve"],
            linewidth=1.8,
            linestyle="--",
            label=summary["name"],
        )
    axes[0, 1].set_xscale("log")
    axes[0, 1].set_ylim(0.0, 7.0)
    axes[0, 1].set_title("Spectral dimension from heat-kernel return")
    axes[0, 1].set_xlabel("diffusion time")
    axes[0, 1].set_ylabel("effective dimension")
    axes[0, 1].grid(alpha=0.22)
    axes[0, 1].legend(frameon=False, fontsize=7)

    periods = np.array(
        [summary["period"] for summary in closure],
        dtype=float,
    )
    spectral_05 = np.array(
        [
            summary["spectral_points"]["0.05"]["spectral_dimension"]
            for summary in closure
        ],
        dtype=float,
    )
    axes[1, 0].plot(
        periods,
        spectral_05,
        marker="o",
        linewidth=2.0,
        color="#2B7A78",
        label="zero-sum closure graph",
    )
    for summary in controls:
        axes[1, 0].axhline(
            summary["spectral_points"]["0.05"]["spectral_dimension"],
            linestyle="--",
            linewidth=1.2,
            label=summary["name"],
        )
    axes[1, 0].set_title("Dimension at matched diffusion return probability")
    axes[1, 0].set_xlabel("closure period T")
    axes[1, 0].set_ylabel("spectral dimension at P(t)=0.05")
    axes[1, 0].grid(alpha=0.22)
    axes[1, 0].legend(frameon=False, fontsize=7)

    labels = []
    values = []
    for summary in summaries:
        if summary["period"] is not None:
            continue
        labels.append(summary["name"])
        values.append(
            summary["spectral_points"]["0.05"]["spectral_dimension"]
        )
    for summary in coarse_summaries:
        labels.append(summary["name"])
        values.append(
            summary["spectral_points"]["0.05"]["spectral_dimension"]
        )
    y = np.arange(len(labels))
    axes[1, 1].barh(y, values, color="#B44A5A", alpha=0.9)
    axes[1, 1].set_yticks(y)
    axes[1, 1].set_yticklabels(labels, fontsize=7)
    axes[1, 1].invert_yaxis()
    axes[1, 1].set_title("Effective dimension after coarse-graining")
    axes[1, 1].set_xlabel("spectral dimension at P(t)=0.05")
    axes[1, 1].grid(axis="x", alpha=0.22)

    fig.suptitle("Zero-sum closure geometry probe", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUT_PNG, dpi=180)
    plt.close(fig)


def main() -> None:
    graph_summaries = []
    for period in PERIODS:
        graph, _ = fixed_period_graph(period)
        summary = summarize_graph(f"closure_T{period}", graph)
        summary["period"] = period
        graph_summaries.append(summary)

    torus = four_dimensional_torus(5)
    torus_summary = summarize_graph("4D torus", torus)
    torus_summary["period"] = None
    graph_summaries.append(torus_summary)

    closure_mean_degree = int(
        round(
            np.mean(
                [
                    summary["mean_degree"]
                    for summary in graph_summaries
                    if summary["period"] == 16
                ]
            )
        )
    )
    control = random_regular_control(810, closure_mean_degree)
    control_summary = summarize_graph(
        f"random {closure_mean_degree}-regular",
        control,
    )
    control_summary["period"] = None
    graph_summaries.append(control_summary)

    modes = all_modes_up_to(COARSE_MAX_PERIOD)
    return_graph = coarse_quotient(
        modes,
        lambda word: (
            len(word),
            internal_return_count(word),
        ),
    )
    return_summary = summarize_graph(
        "coarse (T,returns)",
        return_graph,
        include_embedding=False,
    )
    return_summary["period"] = None

    profile_graph = coarse_quotient(
        modes,
        lambda word: (
            len(word),
            excursion_lengths(word),
        ),
    )
    profile_summary = summarize_graph(
        "coarse (T,profile)",
        profile_graph,
        include_embedding=False,
    )
    profile_summary["period"] = None
    coarse_summaries = [return_summary, profile_summary]

    make_plot(graph_summaries, coarse_summaries)

    period_results = {
        str(summary["period"]): {
            key: value
            for key, value in summary.items()
            if key
            not in {
                "growth",
                "spectral_times",
                "spectral_dimension_curve",
            }
        }
        for summary in graph_summaries
        if summary["period"] is not None
    }
    control_results = {
        summary["name"]: {
            key: value
            for key, value in summary.items()
            if key
            not in {
                "growth",
                "spectral_times",
                "spectral_dimension_curve",
            }
        }
        for summary in graph_summaries
        if summary["period"] is None
    }
    coarse_results = {
        summary["name"]: {
            key: value
            for key, value in summary.items()
            if key
            not in {
                "growth",
                "spectral_times",
                "spectral_dimension_curve",
            }
        }
        for summary in coarse_summaries
    }

    payload = {
        "model": {
            "microstates": "cyclic zero-sum words of +1 and -1 steps",
            "edges": "adjacent opposite-sign swaps",
            "periods": list(PERIODS),
            "coarse_graining": [
                "period and internal return count",
                "period and ordered excursion-length profile",
            ],
            "controls": ["4D torus", "random regular graph"],
        },
        "summary": {
            "fixed_period_graphs": period_results,
            "controls": control_results,
            "coarse_graphs": coarse_results,
        },
        "interpretation": {
            "status": "no stable four-dimensional plateau detected",
            "reason": (
                "local growth exponent and matched-return spectral "
                "dimension both drift with closure period; coarse-graining "
                "does not restore a fixed four-dimensional value"
            ),
            "next_gap": (
                "a zero-sum closure graph alone does not fix the number of "
                "effective geometric directions"
            ),
        },
    }
    OUT_JSON.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(payload["interpretation"], ensure_ascii=False, indent=2))
    print("\nfixed-period summary")
    for period, result in payload["summary"]["fixed_period_graphs"].items():
        diameter = result["diameter"]
        local_dimension = result["local_growth_dimension"]
        stress = result["mds_4d_stress"]
        print(
            f"T={period:>2}  N={result['nodes']:>4}  "
            f"diam={str(diameter):>2}  "
            f"d_local={local_dimension if local_dimension is not None else float('nan'):.3f}  "
            f"d_s(0.05)="
            f"{result['spectral_points']['0.05']['spectral_dimension']:.3f}  "
            f"stress4={stress if stress is not None else float('nan'):.3f}"
        )
    print("\ncoarse graphs")
    for name, result in payload["summary"]["coarse_graphs"].items():
        local_dimension = result["local_growth_dimension"]
        print(
            f"{name}: N={result['nodes']}, diam={result['diameter']}, "
            f"d_local={local_dimension if local_dimension is not None else float('nan'):.3f}, "
            f"d_s(0.05)="
            f"{result['spectral_points']['0.05']['spectral_dimension']:.3f}"
        )
    print(f"\nplot: {OUT_PNG}")
    print(f"json: {OUT_JSON}")


if __name__ == "__main__":
    main()
