#!/usr/bin/env python3
"""Independent numerical check of the H5 dictionary-error proposition.

This script deliberately avoids sparse matrix powers.  For every examined
edge it enumerates all closed walks of length <= k by dynamic programming on
the universal cover, multiplying the scalar c at every traversed edge
midpoint.  The result is compared with

    P_{k,d}(c) = sum_{m even, 2 <= m <= k} m W_d(m-1) c^m.

The smooth and kink fields depend only on x_0.  In that case the full torus
edge scan is reduced to exact edge equivalence classes: transverse translations
do not change any edge weight.  A separate all-edge comparison on small tori
checks this reduction.
"""

import itertools
import math
import sys


AMPLITUDE = 0.2
JUMP_HIGH = 1.2
JUMP_LOW = 0.8
MAIN_K = 4
FAILURES = 0
CHECKS = 0


def check(name, condition, detail=""):
    global FAILURES, CHECKS
    CHECKS += 1
    ok = bool(condition)
    if not ok:
        FAILURES += 1
    tag = "PASS" if ok else "FAIL"
    suffix = ("  " + detail) if detail else ""
    print("  [%s] %s%s" % (tag, name, suffix))


def head(title):
    print("\n" + "=" * 76)
    print(title)
    print("=" * 76)


def direction_steps(d):
    steps = []
    for axis in range(d):
        for sign in (1, -1):
            step = [0] * d
            step[axis] = sign
            steps.append(tuple(step))
    return tuple(steps)


def walk_count_to_neighbor(length, d):
    """Exact number of length-walks from 0 to e_0 on Z^d."""
    current = {(0,) * d: 1}
    steps = direction_steps(d)
    for _ in range(length):
        nxt = {}
        for state, count in current.items():
            for step in steps:
                target = tuple(state[i] + step[i] for i in range(d))
                nxt[target] = nxt.get(target, 0) + count
        current = nxt
    target = tuple(1 if i == 0 else 0 for i in range(d))
    return current.get(target, 0)


def dictionary_coefficients(d, k):
    return tuple(
        (m, walk_count_to_neighbor(m - 1, d))
        for m in range(2, k + 1, 2)
    )


def dictionary_value(c, coefficients):
    return sum(m * count * c ** m for m, count in coefficients)


def smooth_field(x):
    return 1.0 + AMPLITUDE * math.cos(2.0 * math.pi * x)


def kink_field(x):
    periodic = x - math.floor(x)
    triangular = 1.0 - 2.0 * abs(periodic - 0.5)
    return 1.0 + AMPLITUDE * triangular


def jump_field(x):
    periodic = x - math.floor(x)
    return JUMP_HIGH if periodic < 0.5 else JUMP_LOW


def uniform_field(value):
    return lambda _x: value


def edge_closed_walk_weight(n, d, midpoint, first_step, k, field):
    """Sum m times all fixed-first-edge closed walks of lengths 2..k.

    midpoint is the physical midpoint of the oriented edge U -> U + first_step.
    States are lattice displacements measured from U.  Every traversed edge is
    weighted by field evaluated at that edge midpoint.
    """
    steps = direction_steps(d)
    zero = (0,) * d
    current = {first_step: field(midpoint[0])}
    total = 0.0

    # Before this loop current contains walks with the fixed first edge only.
    for length in range(2, k + 1):
        nxt = {}
        for state, weight in current.items():
            state_0 = state[0]
            for step in steps:
                target = tuple(state[i] + step[i] for i in range(d))
                edge_offset = (2 * state_0 + step[0] - first_step[0]) / (2.0 * n)
                edge_weight = field(midpoint[0] + edge_offset)
                nxt[target] = nxt.get(target, 0.0) + weight * edge_weight
        current = nxt
        if length % 2 == 0:
            total += length * current.get(zero, 0.0)

    return total


def brute_force_edge_weight(n, d, midpoint, first_step, k, field):
    """Independent literal enumeration of every walk sequence."""
    steps = direction_steps(d)
    zero = (0,) * d
    total = 0.0

    for length in range(2, k + 1, 2):
        for tail in itertools.product(steps, repeat=length - 1):
            state = first_step
            weight = field(midpoint[0])
            for step in tail:
                edge_offset = (
                    2 * state[0] + step[0] - first_step[0]
                ) / (2.0 * n)
                weight *= field(midpoint[0] + edge_offset)
                state = tuple(state[i] + step[i] for i in range(d))
            if state == zero:
                total += length * weight
    return total


def midpoint_for_class(n, d, q0, axis):
    midpoint = [0.0] * d
    midpoint[0] = (q0 + (0.5 if axis == 0 else 0.0)) / n
    first_step = tuple(1 if i == axis else 0 for i in range(d))
    return tuple(midpoint), first_step


def aligned_edge_scan(n, d, k, field):
    """Exact max scan for fields depending only on x_0."""
    coefficients = dictionary_coefficients(d, k)
    maximum = 0.0
    for q0 in range(n):
        for axis in range(d):
            midpoint, first_step = midpoint_for_class(n, d, q0, axis)
            actual = edge_closed_walk_weight(
                n, d, midpoint, first_step, k, field
            )
            predicted = dictionary_value(field(midpoint[0]), coefficients)
            rel = abs(actual / predicted - 1.0)
            if rel > maximum:
                maximum = rel
    return maximum


def all_edge_scan(n, d, k, field):
    """Direct scan over every vertex and positive edge direction."""
    coefficients = dictionary_coefficients(d, k)
    maximum = 0.0
    for vertex in itertools.product(range(n), repeat=d):
        for axis in range(d):
            midpoint = [
                (vertex[i] + (0.5 if i == axis else 0.0)) / n
                for i in range(d)
            ]
            first_step = tuple(1 if i == axis else 0 for i in range(d))
            actual = edge_closed_walk_weight(
                n, d, tuple(midpoint), first_step, k, field
            )
            predicted = dictionary_value(field(midpoint[0]), coefficients)
            rel = abs(actual / predicted - 1.0)
            if rel > maximum:
                maximum = rel
    return maximum


def metric_series(ns, d, k, field):
    return [(n, aligned_edge_scan(n, d, k, field)) for n in ns]


def print_metric_series(label, series):
    print("  %s" % label)
    previous = None
    for n, error in series:
        ratio = None if previous is None else previous / error
        constant = error * n * n
        ratio_text = "  ratio=----" if ratio is None else "  ratio=%.6f" % ratio
        print(
            "    N=%4d  rel=%.12e  C=rel*N^2=%.9f%s"
            % (n, error, constant, ratio_text)
        )
        previous = error


def loglog_slope(xs, ys):
    log_x = [math.log(x) for x in xs]
    log_y = [math.log(y) for y in ys]
    mean_x = sum(log_x) / len(log_x)
    mean_y = sum(log_y) / len(log_y)
    numerator = sum(
        (x - mean_x) * (y - mean_y) for x, y in zip(log_x, log_y)
    )
    denominator = sum((x - mean_x) ** 2 for x in log_x)
    return numerator / denominator


def relative_change(a, b):
    return abs(b / a - 1.0)


head("F0  Independent implementation checks")

for d in (1, 2, 3):
    for k in (4, 6):
        coefficients = dictionary_coefficients(d, k)
        for c in (0.8, 1.0, 1.2):
            n = 25
            midpoint = tuple([0.13] + [0.0] * (d - 1))
            first_step = tuple(1 if i == 0 else 0 for i in range(d))
            actual = edge_closed_walk_weight(
                n, d, midpoint, first_step, k, uniform_field(c)
            )
            predicted = dictionary_value(c, coefficients)
            check(
                "uniform c: direct walk DP matches dictionary (d=%d,k=%d,c=%.1f)"
                % (d, k, c),
                abs(actual / predicted - 1.0) < 2e-14,
                "rel=%.3e" % abs(actual / predicted - 1.0),
            )

for d in (2, 3):
    n = 7
    k = 4
    worst = 0.0
    for q0 in range(n):
        for axis in range(d):
            midpoint, first_step = midpoint_for_class(n, d, q0, axis)
            dp = edge_closed_walk_weight(
                n, d, midpoint, first_step, k, smooth_field
            )
            brute = brute_force_edge_weight(
                n, d, midpoint, first_step, k, smooth_field
            )
            worst = max(worst, abs(dp / brute - 1.0))
    check(
        "literal walk enumeration equals edge DP (d=%d,k=%d)" % (d, k),
        worst < 3e-13,
        "worst rel=%.3e" % worst,
    )

for d in (2, 3):
    n = 5
    k = 4
    full = all_edge_scan(n, d, k, kink_field)
    classes = aligned_edge_scan(n, d, k, kink_field)
    check(
        "all-edge scan equals aligned equivalence classes (d=%d,N=%d)"
        % (d, n),
        abs(full - classes) < 2e-14,
        "diff=%.3e" % abs(full - classes),
    )


head("F1  Smooth field: max relative error and adjacent-N ratios")

main_sizes = {
    1: (16, 32, 64),
    2: (8, 16, 32),
    3: (8, 16, 32),
}

main_results = {}
for d in (1, 2, 3):
    ns = main_sizes[d]
    check(
        "no wrapping in main smooth runs (d=%d)" % d,
        min(ns) > MAIN_K,
        "min N=%d > k=%d" % (min(ns), MAIN_K),
    )
    series = metric_series(ns, d, MAIN_K, smooth_field)
    main_results[d] = series
    print_metric_series("d=%d, k=%d" % (d, MAIN_K), series)

    errors = [error for _n, error in series]
    constants = [error * n * n for n, error in series]
    ratios = [errors[i] / errors[i + 1] for i in range(len(errors) - 1)]
    check(
        "d=%d smooth errors decrease" % d,
        all(errors[i] > errors[i + 1] for i in range(len(errors) - 1)),
    )
    check(
        "d=%d last adjacent-N ratio approaches 4" % d,
        ratios[-1] > 3.85 and abs(ratios[-1] / 4.0 - 1.0) < 0.04,
        "ratios=%s" % ["%.6f" % r for r in ratios],
    )
    check(
        "d=%d C=rel*N^2 stabilizes" % d,
        relative_change(constants[-2], constants[-1]) < 0.03
        and relative_change(constants[0], constants[1])
        > relative_change(constants[-2], constants[-1]),
        "C=%s" % ["%.9f" % c for c in constants],
    )


head("F2  Fixed-k constant growth in d=2")

growth_k = (4, 8, 12, 16, 20, 24)
growth_n_factor = 16
growth_constants = []
for k in growth_k:
    n = growth_n_factor * k
    error = aligned_edge_scan(n, 2, k, smooth_field)
    constant = error * n * n
    growth_constants.append(constant)
    print(
        "    k=%2d  N=%3d  rel=%.12e  C=%.9f  C/k^3=%.9f"
        % (k, n, error, constant, constant / k ** 3)
    )

last_four_slope = loglog_slope(growth_k[-4:], growth_constants[-4:])
last_interval_slope = math.log(
    growth_constants[-1] / growth_constants[-2]
) / math.log(growth_k[-1] / growth_k[-2])
check(
    "C(k)/k^3 reaches a plateau",
    relative_change(growth_constants[-2] / growth_k[-2] ** 3,
                    growth_constants[-1] / growth_k[-1] ** 3) < 0.002,
    "last two C/k^3=%.9f, %.9f"
    % (
        growth_constants[-2] / growth_k[-2] ** 3,
        growth_constants[-1] / growth_k[-1] ** 3,
    ),
)
check(
    "empirical C(k) exponent is close to 3",
    abs(last_four_slope - 3.0) < 0.08
    and abs(last_interval_slope - 3.0) < 0.02,
    "last-four fit=%.6f, last interval=%.6f"
    % (last_four_slope, last_interval_slope),
)


head("F3  Boundary cases: Lipschitz kink versus jump")

kink_sizes = (16, 32, 64, 128)
kink_series = metric_series(kink_sizes, 2, MAIN_K, kink_field)
print_metric_series("Lipschitz triangular kink, d=2, k=4", kink_series)
kink_errors = [error for _n, error in kink_series]
kink_ratios = [
    kink_errors[i] / kink_errors[i + 1]
    for i in range(len(kink_errors) - 1)
]
kink_constants = [error * n for n, error in kink_series]
check(
    "kink errors follow O(a): adjacent-N ratio approaches 2",
    kink_ratios[-1] > 1.98 and kink_ratios[-1] < 2.02,
    "ratios=%s" % ["%.6f" % r for r in kink_ratios],
)
check(
    "kink C1=rel*N stabilizes",
    relative_change(kink_constants[-2], kink_constants[-1]) < 0.01,
    "C1=%s" % ["%.9f" % c for c in kink_constants],
)

jump_sizes = (16, 32, 64, 128)
jump_series = metric_series(jump_sizes, 2, MAIN_K, jump_field)
print_metric_series("Jump c, d=2, k=4", jump_series)
jump_errors = [error for _n, error in jump_series]
jump_ratios = [
    jump_errors[i] / jump_errors[i + 1]
    for i in range(len(jump_errors) - 1)
]
check(
    "jump error is O(1): adjacent-N ratios remain 1",
    all(abs(ratio - 1.0) < 2e-12 for ratio in jump_ratios),
    "ratios=%s" % ["%.12f" % r for r in jump_ratios],
)
check(
    "jump error has no convergence",
    max(jump_errors) - min(jump_errors) < 2e-13,
    "range=%.3e" % (max(jump_errors) - min(jump_errors)),
)


head("Summary")
for d in (1, 2, 3):
    errors = [error for _n, error in main_results[d]]
    ratios = [errors[i] / errors[i + 1] for i in range(len(errors) - 1)]
    constants = [
        error * n * n for n, error in main_results[d]
    ]
    print(
        "  d=%d: rel=%s; ratios=%s; C_last_two=%.9f, %.9f"
        % (
            d,
            ["%.6e" % e for e in errors],
            ["%.6f" % r for r in ratios],
            constants[-2],
            constants[-1],
        )
    )
print(
    "  k growth (d=2): C=%s; exponent_last_four=%.6f; "
    "exponent_last_interval=%.6f"
    % (
        ["%.6g" % c for c in growth_constants],
        last_four_slope,
        last_interval_slope,
    )
)
print(
    "  kink (d=2): ratios=%s; jump ratios=%s"
    % (
        ["%.6f" % r for r in kink_ratios],
        ["%.6f" % r for r in jump_ratios],
    )
)

if FAILURES:
    print("\nFAILED: %d of %d checks." % (FAILURES, CHECKS))
    sys.exit(1)

print("\nALL %d CHECKS PASSED." % CHECKS)
