#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R23_dim_desc_probe.py
=====================

Numerical probe for R23_dimension_descendant_selection.md.

The probe does not claim a Zero-native dimension sector. It checks the
algebra of three model classes:

1. product sectors with no cross-axis loss:
       F_D = a**D
   The log growth is linear in D, so there is no interior D=4 peak.

2. coherence-descendant model:
       F_D = B * D * q**D
   A unique integer maximum at D=4 occurs exactly for 3/4 < q < 4/5.

3. lifetime-linked model:
       F_D(L) = D * exp(-D/L)
   The unique integer maximum is D=L.

4. free-lifetime normalization:
       P(L) = L/e
   grows without a finite maximum per macrocycle; the per-step rate
   r(L) = (log L - 1)/L peaks at L=e^2 in the continuum and at the
   even closure length L=8 among integer candidates.

5. GR-compatible product survival:
       S_D = a**D, 0 < a < 1
   Strictly decreases with D. Since GR-compatible sectors require D>=4,
   D=4 is the unique conditional-survival maximum. Absolute surviving
   branch counts increase with D, so upgrading this fraction to a
   long-term evolution-layer share requires the open bridge EVO-NORM.

6. External survival multipliers (secondary):
       F_D^GR = s_D * B * D * q**D
   If one abandons product survival, s_4/s_D must exceed the exact
   base-fitness ratio.

Run:
    python3 R23_dim_desc_probe.py
    python3 R23_dim_desc_probe.py --json R23_dim_desc_results.json
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_JSON = ROOT / "R23_dim_desc_results.json"
DIMS = tuple(range(1, 25))


def product_weight(D: int, a: float) -> float:
    return a**D


def descendant_weight(D: int, q: float, B: float = 1.0) -> float:
    return B * D * q**D


def lifetime_weight(D: int, L: float) -> float:
    return D * math.exp(-D / L)


def argmax_unique(values: dict[int, float], rel_tol: float = 1.0e-12) -> dict[str, object]:
    maximum = max(values.values())
    winners = [
        D
        for D, value in values.items()
        if math.isclose(value, maximum, rel_tol=rel_tol, abs_tol=0.0)
    ]
    return {
        "argmax": winners[0],
        "winners": winners,
        "unique": len(winners) == 1,
        "maximum": maximum,
    }


def product_no_go(a: float) -> dict[str, object]:
    values = {D: product_weight(D, a) for D in DIMS}
    slopes = [
        values[D + 1] / values[D]
        for D in DIMS[:-1]
    ]
    return {
        "a": a,
        "log_growth": {str(D): math.log(values[D]) for D in DIMS},
        "adjacent_ratios": slopes,
        "ratio_spread": max(slopes) - min(slopes),
        "argmax": argmax_unique(values),
    }


def descendant_case(q: float) -> dict[str, object]:
    values = {D: descendant_weight(D, q) for D in DIMS}
    result = argmax_unique(values)
    result.update(
        {
            "q": q,
            "F4_over_F3": descendant_weight(4, q) / descendant_weight(3, q),
            "F5_over_F4": descendant_weight(5, q) / descendant_weight(4, q),
            "weights": {str(D): values[D] for D in DIMS},
        }
    )
    return result


def lifetime_case(L: int) -> dict[str, object]:
    values = {D: lifetime_weight(D, L) for D in DIMS}
    result = argmax_unique(values)
    result.update(
        {
            "L": L,
            "weights": {str(D): values[D] for D in DIMS},
        }
    )
    return result


def lifetime_family_case(L: int) -> dict[str, object]:
    return {
        "L": L,
        "peak_dimension": L,
        "per_cycle_peak": L * math.exp(-1.0),
        "per_step_log_growth": (math.log(L) - 1.0) / L,
    }


def gr_survival_case(
    q: float,
    s_four: float,
    label: str,
) -> dict[str, object]:
    base = {D: descendant_weight(D, q) for D in DIMS}
    thresholds = {
        D: base[D] / base[4]
        for D in DIMS
        if D != 4
    }
    threshold_dimension = max(thresholds, key=thresholds.get)
    boosted = {
        D: base[D] * (s_four if D == 4 else 1.0)
        for D in DIMS
    }
    return {
        "label": label,
        "q": q,
        "s_four": s_four,
        "base_winners": argmax_unique(base)["winners"],
        "boosted_winners": argmax_unique(boosted)["winners"],
        "required_s4_ratio": thresholds[threshold_dimension],
        "threshold_dimension": threshold_dimension,
    }


def exact_direction_survival(L: int) -> float:
    return math.comb(L, L // 2) / 2**L


def product_survival_case(L: int) -> dict[str, object]:
    gr_dims = tuple(D for D in DIMS if D >= 4)
    a = exact_direction_survival(L)
    values = {D: a**D for D in gr_dims}
    result = argmax_unique(values)
    return {
        "L": L,
        "a": a,
        "gr_dims": list(gr_dims),
        "winners": result["winners"],
        "maximum": result["maximum"],
        "S4": values[4],
    }


def q_scan(
    q_min: float = 0.50,
    q_max: float = 0.90,
    steps: int = 401,
) -> list[dict[str, object]]:
    rows = []
    for index in range(steps + 1):
        q = q_min + (q_max - q_min) * index / steps
        result = descendant_case(q)
        rows.append(
            {
                "q": q,
                "argmax": result["argmax"],
                "unique": result["unique"],
                "winners": result["winners"],
            }
        )
    return rows


def scan_unique_four(
    q_min: float = 0.50,
    q_max: float = 0.90,
    steps: int = 40000,
) -> list[dict[str, float]]:
    intervals: list[dict[str, float]] = []
    start: float | None = None
    previous = 0.50
    for index in range(steps + 1):
        q = q_min + (q_max - q_min) * index / steps
        result = descendant_case(q)
        in_interval = result["unique"] and result["winners"] == [4]
        if in_interval and start is None:
            start = q
        if not in_interval and start is not None:
            intervals.append({"start": start, "end": previous})
            start = None
        previous = q
    if start is not None:
        intervals.append({"start": start, "end": q_max})
    return intervals


def build_payload() -> dict[str, object]:
    product_cases = [product_no_go(a) for a in (0.80, 1.00, 1.20)]
    descendant_cases = [
        descendant_case(q)
        for q in (0.70, 0.75, math.exp(-0.25), 0.80, 0.82)
    ]
    lifetime_cases = [lifetime_case(L) for L in range(1, 17)]
    lifetime_family = [
        lifetime_family_case(L)
        for L in range(4, 25, 2)
    ]
    gr_survival_cases = [
        gr_survival_case(0.70, 1.08, "above the 3D threshold"),
        gr_survival_case(0.7788, 1.00, "base model already selects 4D"),
        gr_survival_case(math.exp(-0.25), 1.00, "exact 4D window point"),
        gr_survival_case(0.82, 1.01, "higher but below the 5D threshold"),
        gr_survival_case(0.82, 1.03, "above the 5D threshold"),
    ]
    product_survival_cases = [
        product_survival_case(L)
        for L in (2, 4, 6, 8)
    ]
    per_step_winner = max(
        lifetime_family,
        key=lambda row: row["per_step_log_growth"],
    )
    return {
        "model": {
            "dims": list(DIMS),
            "product": "F_D = a**D",
            "descendant": "F_D = B*D*q**D",
            "lifetime": "F_D(L) = D*exp(-D/L)",
            "status": (
                "conditional algebraic model; DIM-SECTOR and DIM-COST "
                "are not derived from Z0; EVO-NORM is open"
            ),
        },
        "product_no_go": product_cases,
        "descendant_cases": descendant_cases,
        "unique_four_intervals": scan_unique_four(),
        "q_scan": q_scan(),
        "lifetime_cases": lifetime_cases,
        "lifetime_family": lifetime_family,
        "product_survival_cases": product_survival_cases,
        "gr_survival_cases": gr_survival_cases,
        "summary": {
            "unique_four_window": [0.75, 0.8],
            "unique_four_window_open": "3/4 < q < 4/5",
            "lifetime_identity": "argmax_D D*exp(-D/L) = L",
            "free_lifetime_per_cycle_optimum": "unbounded",
            "free_lifetime_per_step_continuous_optimum": math.exp(2.0),
            "free_lifetime_per_step_even_integer_optimum": per_step_winner["L"],
            "survival_ordering_proof": (
                "a_L=binom(L,L/2)/2**L; S_D=a_L**D decreases; "
                "D>=4 gives argmax 4 in conditional survival rate"
            ),
            "evo_norm_open": True,
            "gr_surv_threshold": "s4/sD > F_D^0/F_4^0",
            "no_unconditional_zero_to_four": True,
        },
    }


def print_report(payload: dict[str, object]) -> None:
    print("R23 DIM-DESC probe")
    print("=" * 72)
    print("product-sector no-go")
    for row in payload["product_no_go"]:
        print(
            "  a=%.2f  argmax=%s  unique=%-5s  ratio-spread=%.3e"
            % (
                row["a"],
                row["argmax"]["argmax"],
                row["argmax"]["unique"],
                row["ratio_spread"],
            )
        )

    print("")
    print("coherence-descendant model")
    for row in payload["descendant_cases"]:
        print(
            "  q=%.8f  argmax=%-2d  unique=%-5s  winners=%s"
            % (
                row["q"],
                row["argmax"],
                row["unique"],
                row["winners"],
            )
        )

    print("")
    print("lifetime-linked model")
    for row in payload["lifetime_cases"]:
        print(
            "  L=%-2d  argmax=%-2d  unique=%s"
            % (row["L"], row["argmax"], row["unique"])
        )

    print("")
    print("free-lifetime normalization on even closure lengths")
    for row in payload["lifetime_family"]:
        print(
            "  L=%-2d  per-cycle-peak=%.8f  per-step-rate=%.8f"
            % (
                row["L"],
                row["per_cycle_peak"],
                row["per_step_log_growth"],
            )
        )

    print("")
    print("exact GR-compatible product survival")
    for row in payload["product_survival_cases"]:
        print(
            "  L=%-2d  a_L=%.8f  S4=%.8f  argmax-on-D>=4=%s"
            % (row["L"], row["a"], row["S4"], row["winners"])
        )

    print("")
    print("external survival multiplier thresholds (secondary)")
    for row in payload["gr_survival_cases"]:
        print(
            "  %-38s base=%-8s boosted=%-8s threshold=%.8f at D=%d"
            % (
                row["label"],
                row["base_winners"],
                row["boosted_winners"],
                row["required_s4_ratio"],
                row["threshold_dimension"],
            )
        )

    intervals = payload["unique_four_intervals"]
    print("")
    print("unique D=4 intervals on q in [0.50, 0.90]:")
    for interval in intervals:
        print("  [%.8f, %.8f]" % (interval["start"], interval["end"]))
    print("")
    print("honest boundary: DIM-SECTOR and DIM-COST are named inputs; EVO-NORM is open.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", type=Path, default=None)
    args = parser.parse_args()

    payload = build_payload()
    print_report(payload)
    if args.json is not None:
        args.json.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        print("json: %s" % args.json)


if __name__ == "__main__":
    main()
