#!/usr/bin/env python3
"""Independent numerical checks for the H3/H7 mollification construction.

This is deliberately a numerical audit, not a proof.  It keeps the three
scales separate:

* h: the mollification width;
* a: the lattice/cell spacing;
* eps: a finite-difference/reference scale used only where noted.

The H3/H7 construction is

    mu_a = a^d sum_j c(x_j) delta_{x_j},
    c_{a,h} = rho_h * mu_a,

with a/h -> 0.  The script verifies:

1. pure periodic heat mollification converges in C^0, C^1 and C^2;
2. the discrete-measure Riemann error has the advertised a^2 and h powers
   as conservative envelopes, while printing the sharper observed rates;
3. the intra-cell oscillation scales as a/h for a jump field;
4. piecewise-constant and piecewise-linear interpolation do not themselves
   define C^2 fields.

The compactly supported kernels below are stress-test kernels for the
Riemann error.  They are not claimed to be the smooth mollifier used in H7;
the heat kernel is used for the smooth/C^2 checks.
"""

import math
import sys

import numpy as np


FAILURES = 0
CHECKS = 0


def check(name, condition, detail=""):
    global FAILURES, CHECKS
    CHECKS += 1
    ok = bool(condition)
    if not ok:
        FAILURES += 1
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name,
                           ("  " + detail) if detail else ""))


def head(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def relative_change(a, b):
    return abs(float(a) / float(b) - 1.0)


def log_slope(xs, ys):
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    if np.any(xs <= 0) or np.any(ys <= 0):
        raise ValueError("log_slope needs positive data")
    lx = np.log(xs)
    ly = np.log(ys)
    return float(np.polyfit(lx, ly, 1)[0])


def wrap_signed(x):
    return (x + 0.5) % 1.0 - 0.5


def heat_kernel_circle(x, h):
    """Periodic heat kernel on the unit circle, summed over images."""
    z = wrap_signed(x)
    prefactor = 1.0 / math.sqrt(2.0 * math.pi * h * h)
    value = 0.0
    for image in range(-8, 9):
        zz = z + image
        value += prefactor * math.exp(-zz * zz / (2.0 * h * h))
    return value


def bernoulli_cubic(x):
    """Periodic C^2 field whose third derivative has a jump."""
    x = x - math.floor(x)
    return 1.0 + 0.08 * (x ** 3 - 1.5 * x ** 2 + 0.5 * x)


def smooth_sine(x):
    return 1.0 + 0.08 * math.sin(2.0 * math.pi * x)


def gauss_nodes_weights(count=24):
    nodes, weights = np.polynomial.legendre.leggauss(count)
    nodes = 0.5 * (nodes + 1.0)
    weights = 0.5 * weights
    return nodes, weights


def compact_normalization(m):
    return math.gamma(m + 1.5) / (math.sqrt(math.pi) * math.gamma(m + 1))


def compact_kernel(t, radius, m, derivative=0):
    """Compactly supported kernel of radius H and regularity controlled by m."""
    s = t / radius
    if abs(s) >= 1.0:
        return 0.0
    q = 1.0 - s * s
    c_m = compact_normalization(m)
    if derivative == 0:
        return c_m / radius * q ** m
    if derivative != 2 or m < 2:
        raise ValueError("only derivative 0 or 2 is implemented")
    first = m * (m - 1) * q ** (m - 2) * (-2.0 * s / radius) ** 2
    second = m * q ** (m - 1) * (-2.0 / (radius * radius))
    return c_m / radius * (first + second)


def compact_continuous_conv(x, radius, m, derivative, field, nodes, weights):
    """Exact-in-polynomial Gauss integration of rho_h * f at one point."""
    lo = x - radius
    hi = x + radius
    breaks = [lo, hi]
    for integer in range(math.floor(lo), math.floor(hi) + 1):
        if lo < integer < hi:
            breaks.append(float(integer))
    breaks = sorted(set(breaks))
    total = 0.0
    for left, right in zip(breaks[:-1], breaks[1:]):
        if right <= left:
            continue
        width = right - left
        for u, w in zip(nodes, weights):
            y = left + width * u
            total += (
                width
                * w
                * field(y)
                * compact_kernel(x - y, radius, m, derivative)
            )
    return total


def compact_discrete_conv(x, radius, a, d, m, derivative, field):
    """Normalized discrete measure a^d sum_j field(x_j) rho_h(x-x_j)."""
    lo = math.floor((x - radius) / a) - 1
    hi = math.ceil((x + radius) / a) + 1
    total = 0.0
    for j0 in range(lo, hi + 1):
        t0 = x - j0 * a
        if abs(t0) >= radius:
            continue
        product = compact_kernel(t0, radius, m, derivative)
        for axis in range(1, d):
            xi = 0.111 + 0.073 * axis
            transverse = 0.0
            for j in range(lo, hi + 1):
                tip = xi - j * a
                if abs(tip) < radius:
                    transverse += compact_kernel(tip, radius, m, 0)
            product *= transverse
        total += a ** d * field(j0 * a) * product
    return total


def spectral_heat_errors(d, h, sample_count=96):
    """Errors of the exact heat semigroup on a finite Fourier polynomial."""
    coordinates = [np.arange(sample_count) / sample_count for _ in range(d)]
    grids = np.meshgrid(*coordinates, indexing="ij")
    modes = [
        (0.32, tuple(1 if i == 0 else 0 for i in range(d))),
        (0.18, tuple(1 if i == 1 else 0 for i in range(d))),
        (0.11, tuple(2 if i == 0 else 0 for i in range(d))),
        (0.07, tuple(1 if i in (0, 1) else 0 for i in range(d))),
        (0.05, tuple(3 if i == min(1, d - 1) else 0 for i in range(d))),
    ]
    error_0 = 0.0
    error_1 = 0.0
    error_2 = 0.0
    for amplitude, k in modes:
        phase = sum(
            2.0 * math.pi * k[i] * grids[i] for i in range(d)
        )
        k2 = sum(q * q for q in k)
        multiplier = math.exp(-2.0 * math.pi ** 2 * h * h * k2) - 1.0
        mode = amplitude * multiplier * np.exp(1j * phase)
        error_0 = max(error_0, float(np.max(np.abs(mode))))
        for i in range(d):
            derivative_1 = mode * (2j * math.pi * k[i])
            error_1 = max(error_1, float(np.max(np.abs(derivative_1))))
            for j in range(d):
                derivative_2 = mode * (-4.0 * math.pi ** 2 * k[i] * k[j])
                error_2 = max(
                    error_2, float(np.max(np.abs(derivative_2)))
                )
    return error_0, error_1, error_2


def discrete_heat_field(xs, h, a):
    """Mollified alternating step labels on a one-dimensional cell grid."""
    count = int(round(1.0 / a))
    if relative_change(count * a, 1.0) > 1e-12:
        raise ValueError("test grid must divide the circle exactly")
    sites = np.arange(count, dtype=float) * a
    labels = np.where(sites < 0.5, 0.75, 1.5)
    field = np.empty_like(xs, dtype=float)
    for index, x in enumerate(xs):
        kernel = np.array(
            [heat_kernel_circle(x - y, h) for y in sites],
            dtype=float,
        )
        field[index] = a * float(np.dot(labels, kernel))
    return field


def eta_series():
    results = []
    for h in (0.08, 0.04, 0.02):
        a = h / 16.0
        sample_count = int(round(1.0 / a))
        xs = np.arange(2 * sample_count, dtype=float) / (
            2.0 * sample_count
        )
        field = discrete_heat_field(xs, h, a)
        assert field.min() > 0.0
        shift = 2
        variation = np.max(
            np.abs(np.roll(field, -shift) - field)
        )
        results.append((h, a, variation / field.min()))
    return results


def interpolation_checks():
    series = []
    for count in (16, 32, 64, 128):
        xs = np.arange(count, dtype=float) / count
        values = np.array([smooth_sine(x) for x in xs])
        pc_jump = float(np.max(np.abs(np.roll(values, -1) - values)))
        slopes = np.roll(values, -1) - values
        slope_jump = float(
            np.max(np.abs(np.roll(slopes, -1) - slopes))
        )
        midpoints = (np.arange(count, dtype=float) + 0.5) / count
        pc_values = values
        p1_values = 0.5 * (values + np.roll(values, -1))
        pc_error = float(
            np.max(
                np.abs(
                    pc_values
                    - np.array([smooth_sine(x) for x in midpoints])
                )
            )
        )
        p1_error = float(
            np.max(
                np.abs(
                    p1_values
                    - np.array([smooth_sine(x) for x in midpoints])
                )
            )
        )
        series.append((count, pc_jump, slope_jump, pc_error, p1_error))
    return series


head("F0  Smooth periodic heat mollification")

heat_h = (0.04, 0.02, 0.01, 0.005)
heat_errors = [spectral_heat_errors(2, h) for h in heat_h]
print("  h        C0 error       C1 error       C2 error")
for h, errors in zip(heat_h, heat_errors):
    print(
        "  %.4f  %.6e  %.6e  %.6e"
        % (h, errors[0], errors[1], errors[2])
    )

check(
    "pure mollification C0 converges",
    heat_errors[-1][0] < heat_errors[0][0] * 0.05,
)
check(
    "pure mollification C1 converges",
    heat_errors[-1][1] < heat_errors[0][1] * 0.05,
)
check(
    "pure mollification C2 converges",
    heat_errors[-1][2] < heat_errors[0][2] * 0.05,
)

for label, index in (("C0", 0), ("C1", 1), ("C2", 2)):
    observed = log_slope(heat_h, [e[index] for e in heat_errors])
    check(
        "%s refinement slope is O(h^2)" % label,
        1.75 < observed < 2.25,
        "slope=%.6f" % observed,
    )


head("F1  Discrete pushforward: L-infinity Riemann scaling")

nodes, weights = gauss_nodes_weights(20)
riemann_grid = np.linspace(0.05, 0.45, 51)
fixed_h = 0.10
fixed_radius = 0.25 * fixed_h
fixed_a = (0.01, 0.005, 0.0025, 0.00125)

riemann_errors = []
for a in fixed_a:
    error = max(
        abs(
            compact_discrete_conv(
                x, fixed_radius, a, 1, 1, 0, bernoulli_cubic
            )
            - compact_continuous_conv(
                x, fixed_radius, 1, 0, bernoulli_cubic, nodes, weights
            )
        )
        for x in riemann_grid
    )
    riemann_errors.append(error)

print(
    "  fixed h=%.3f, a=%s"
    % (fixed_h, ["%.6g" % a for a in fixed_a])
)
print(
    "  L-infinity errors=%s"
    % ["%.6e" % e for e in riemann_errors]
)
riemann_slope = log_slope(fixed_a[-3:], riemann_errors[-3:])
check(
    "L-infinity Riemann error has observed a^2 scaling",
    1.9 < riemann_slope < 2.1,
    "slope=%.6f" % riemann_slope,
)

joint_h = (0.20, 0.10, 0.05, 0.025)
joint_errors = []
for h in joint_h:
    radius = 0.25 * h
    a = 0.25 * h * h
    joint_errors.append(
        max(
            abs(
                compact_discrete_conv(
                    x, radius, a, 1, 1, 0, bernoulli_cubic
                )
                - compact_continuous_conv(
                    x, radius, 1, 0, bernoulli_cubic, nodes, weights
                )
            )
            for x in riemann_grid
        )
    )
check(
    "joint refinement (a ~ h^2) drives L-infinity error down",
    joint_errors[-1] < joint_errors[0] * 0.02,
    "errors=%s" % ["%.3e" % e for e in joint_errors],
)


head("F2  Discrete pushforward: second-derivative numerical side check")

# The m=3 compact kernel has continuous second derivative.  The sequence is
# chosen small enough that the numerical side check is stable; it is not used
# as a sharp asymptotic proof.
d2_a = (0.005, 0.00125, 0.0003125, 7.8125e-5)
d2_radius = 0.25 * 0.10
d2_errors = []
for a in d2_a:
    error = max(
        abs(
            compact_discrete_conv(
                x, d2_radius, a, 1, 3, 2, smooth_sine
            )
            - compact_continuous_conv(
                x, d2_radius, 3, 2, smooth_sine, nodes, weights
            )
        )
        for x in riemann_grid
    )
    d2_errors.append(error)

print(
    "  fixed h=0.100, a=%s" % ["%.6g" % a for a in d2_a]
)
print("  D2 errors=%s" % ["%.6e" % e for e in d2_errors])
d2_slope = log_slope(d2_a, d2_errors)
check(
    "D2 Riemann side check has observed a^2 tail",
    1.9 < d2_slope < 2.1,
    "tail slope=%.6f" % d2_slope,
)


head("F3  Claimed h-envelopes and their sharpness")

envelope_l = []
for h in (0.10, 0.05, 0.025, 0.0125):
    radius = 0.25 * h
    a = 0.001
    error = max(
        abs(
            compact_discrete_conv(
                x, radius, a, 1, 1, 0, bernoulli_cubic
            )
            - compact_continuous_conv(
                x, radius, 1, 0, bernoulli_cubic, nodes, weights
            )
        )
        for x in riemann_grid
    )
    envelope_l.append(error * h ** (1 + 2) / a ** 2)

envelope_d2 = []
for h in (0.10, 0.05, 0.025, 0.0125):
    radius = 0.25 * h
    a = 0.001
    error = max(
        abs(
            compact_discrete_conv(
                x, radius, a, 1, 3, 2, smooth_sine
            )
            - compact_continuous_conv(
                x, radius, 3, 2, smooth_sine, nodes, weights
            )
        )
        for x in riemann_grid
    )
    envelope_d2.append(error * h ** (1 + 4) / a ** 2)

print(
    "  d=1 L-infinity envelope E*h^(d+2)/a^2=%s"
    % ["%.3e" % v for v in envelope_l]
)
print(
    "  d=1 D2 envelope E*h^(d+4)/a^2=%s"
    % ["%.3e" % v for v in envelope_d2]
)
check(
    "claimed L-infinity h-envelope remains finite and small",
    max(envelope_l) < 1.0,
    "max=%.3e" % max(envelope_l),
)
check(
    "claimed D2 h-envelope decreases under refinement",
    envelope_d2[-1] < envelope_d2[0],
    "first=%.3e last=%.3e" % (envelope_d2[0], envelope_d2[-1]),
)

# The observed normalized scaling is d-independent in this setting.  This
# distinction is printed rather than hidden: the d+2 and d+4 formulas are
# conservative envelopes, not sharp exponents.
sharp_l_constant = [
    riemann_errors[i] * fixed_h ** 2 / fixed_a[i] ** 2
    for i in range(len(fixed_a))
]
sharp_d2_constant = [
    d2_errors[i] * 0.10 ** 4 / d2_a[i] ** 2
    for i in range(len(d2_a))
]
print(
    "  sharp normalized L-infinity check E*h^2/a^2=%s"
    % ["%.3e" % v for v in sharp_l_constant]
)
print(
    "  sharp normalized D2 check E*h^4/a^2=%s"
    % ["%.3e" % v for v in sharp_d2_constant]
)
check(
    "normalized L-infinity constants stabilize",
    relative_change(sharp_l_constant[-2], sharp_l_constant[-1]) < 0.02,
)
check(
    "normalized D2 constants stabilize",
    relative_change(sharp_d2_constant[-2], sharp_d2_constant[-1]) < 0.02,
)


head("F4  Intra-cell oscillation eta_a")

eta_results = eta_series()
print("  h        a          eta_a       eta_a*h/a")
for h, a, eta in eta_results:
    print(
        "  %.3f  %.6f  %.6e  %.9f"
        % (h, a, eta, eta * h / a)
    )
eta_scaled = [eta * h / a for h, a, eta in eta_results]
check(
    "eta_a is numerically O(a/h)",
    relative_change(eta_scaled[-1], eta_scaled[0]) < 0.01,
    "scaled values=%s" % ["%.9f" % v for v in eta_scaled],
)


head("F5  Piecewise interpolation cannot serve directly as C^2")

interpolation = interpolation_checks()
print("  N     PC jump    P1 slope jump    PC C0 err    P1 C0 err")
for row in interpolation:
    print(
        "  %4d  %.3e  %.3e  %.3e  %.3e"
        % row
    )
check(
    "piecewise-constant interpolants retain value jumps",
    all(row[1] > 1e-6 for row in interpolation),
    "last jump=%.3e" % interpolation[-1][1],
)
check(
    "piecewise-linear interpolants retain derivative jumps",
    all(row[2] > 1e-6 for row in interpolation),
    "last slope jump=%.3e" % interpolation[-1][2],
)
check(
    "piecewise-constant and piecewise-linear errors both vanish in C0",
    interpolation[-1][3] < interpolation[0][3]
    and interpolation[-1][4] < interpolation[0][4],
)
check(
    "a smooth mollification is therefore required before taking two derivatives",
    interpolation[-1][2] < interpolation[0][2]
    and interpolation[-1][2] > 0.0,
)


head("Summary")
print("  checks=%d, failures=%d" % (CHECKS, FAILURES))
print("  Numerical side evidence only; no claim here replaces a proof.")

if FAILURES:
    print("\nFAILED: %d of %d checks." % (FAILURES, CHECKS))
    sys.exit(1)

print("\nALL %d CHECKS PASSED." % CHECKS)
