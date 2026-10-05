#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R25_native_pair_cost_probe.py

Independent probe for:
  1. rotation-class ledger purity q_L;
  2. the single-carrier model F_D = D q^D;
  3. the pair-carrier model F_D = C(D,2) q^D;
  4. the exact windows and the L=4 peak.

The probe does not import any Zero program.
"""

from __future__ import annotations

import io
import json
import os
from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb, exp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "R25_native_pair_cost_results.json")


def canonical(word: tuple[int, ...]) -> tuple[int, ...]:
    return min(word[i:] + word[:i] for i in range(len(word)))


def rotation_spectrum(L: int) -> tuple[Counter[int], int]:
    words = [w for w in product((1, -1), repeat=L) if sum(w) == 0]
    classes = Counter(canonical(w) for w in words)
    return classes, len(words)


def ledger_purity(L: int) -> tuple[Fraction, dict[int, int], int]:
    classes, total = rotation_spectrum(L)
    sizes = Counter(classes.values())
    q = sum(Fraction(size, total) ** 2 for size in classes.values())
    return q, dict(sorted(sizes.items())), total


def argmax(values: dict[int, Fraction]) -> tuple[list[int], Fraction]:
    maximum = max(values.values())
    return [key for key, value in values.items() if value == maximum], maximum


def single_weights(q: Fraction, dmax: int = 40) -> dict[int, Fraction]:
    return {D: D * q**D for D in range(1, dmax + 1)}


def pair_weights(q: Fraction, dmax: int = 40) -> dict[int, Fraction]:
    return {D: comb(D, 2) * q**D for D in range(2, dmax + 1)}


def main() -> None:
    rows = []
    for L in range(2, 21, 2):
        q, sizes, total = ledger_purity(L)
        one_winners, _ = argmax(single_weights(q))
        pair_winners, _ = argmax(pair_weights(q))
        bound = Fraction(L, total)
        rows.append(
            {
                "L": L,
                "N": total,
                "class_size_distribution": sizes,
                "q": str(q),
                "q_float": float(q),
                "q_le_L_over_N": q <= bound,
                "L_over_N": str(bound),
                "L_over_N_lt_3_4": bound < Fraction(3, 4),
                "single_argmax": one_winners,
                "pair_argmax": pair_winners,
            }
        )

    q4 = Fraction(5, 9)
    q4_single_winners, _ = argmax(single_weights(q4))
    q4_pair_winners, _ = argmax(pair_weights(q4))

    low_boundary_single, _ = argmax(single_weights(Fraction(3, 4)))
    high_boundary_single, _ = argmax(single_weights(Fraction(4, 5)))
    low_boundary_pair, _ = argmax(pair_weights(Fraction(1, 2)))
    high_boundary_pair, _ = argmax(pair_weights(Fraction(3, 5)))

    e_quarter = Fraction(0)
    # Keep the exp(-1/L) comparison in floating point intentionally: it is not the
    # native ledger value and is used only as a negative-control comparison.
    q_cont = exp(-0.25)
    cont_winners = []
    cont_max = -1.0
    for D in range(2, 41):
        value = comb(D, 2) * q_cont**D
        if value > cont_max:
            cont_max = value
            cont_winners = [D]
        elif abs(value - cont_max) < 1.0e-15:
            cont_winners.append(D)

    payload = {
        "model": "rotation-class purity and pair-carrier dimension selection",
        "rows": rows,
        "q4": str(q4),
        "q4_single_argmax": q4_single_winners,
        "q4_pair_argmax": q4_pair_winners,
        "single_boundaries": {
            "q=3/4": low_boundary_single,
            "q=4/5": high_boundary_single,
        },
        "pair_boundaries": {
            "q=1/2": low_boundary_pair,
            "q=3/5": high_boundary_pair,
        },
        "continuous_q_exp_minus_1_over_4_pair_argmax": cont_winners,
        "continuous_q": q_cont,
        "continuous_argmax_note": "negative control only; not a native ledger value",
    }

    with io.open(OUT, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")

    print("R25 probe")
    for row in rows:
        print(
            "  L=%-2d N=%-6d q=%-12s single=%s pair=%s"
            % (
                row["L"],
                row["N"],
                row["q"],
                row["single_argmax"][:5],
                row["pair_argmax"][:5],
            )
        )
    print("  exact q4 =", q4)
    print("  q4 single argmax =", q4_single_winners)
    print("  q4 pair argmax   =", q4_pair_winners)
    print("  pair boundaries  =", payload["pair_boundaries"])
    print("  q=exp(-1/4) pair argmax (control) =", cont_winners)
    print("  wrote", OUT)


if __name__ == "__main__":
    main()
