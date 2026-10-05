#!/usr/bin/env python3
"""Numerical probes for the critical free-fermion interval modular Hamiltonian.

This script is intentionally offline and self-contained.  It constructs the
restricted correlation matrix

    C_ij = sin(pi (i-j) / 2) / (pi (i-j)),   i != j,
    C_ii = 1/2,

and diagonalizes it in mpmath at high precision.  Everything below is a
finite-N numerical check, not a proof of strong resolvent convergence.
"""

import sys

import mpmath as mp


# 80 digits is comfortably above the requested ~50 digits and keeps spurious
# even-distance matrix elements below 1e-40 through N=20.
mp.mp.dps = 80

NDIMS = (8, 12, 16, 20)
CHECKERBOARD_TOL = mp.mpf("1e-40")
HERMITIAN_TOL = mp.mpf("1e-60")
PHS_TOL = mp.mpf("1e-50")

ASSERTIONS = 0
FAILURES = 0


def check(name, condition, detail=""):
    """Print one [v]/[x] assertion and update the summary counters."""
    global ASSERTIONS, FAILURES
    ASSERTIONS += 1
    ok = bool(condition)
    if not ok:
        FAILURES += 1
    suffix = ("   " + detail) if detail else ""
    print("  [%s] %s%s" % ("v" if ok else "x", name, suffix))
    return ok


def fmt(value, digits=12):
    return mp.nstr(value, digits)


def build_restricted_matrix(N):
    """Build C, diagonalize it, and form h = log((1-C)/C)."""
    C = mp.matrix(N, N)
    for i in range(N):
        for j in range(N):
            distance = i - j
            if distance == 0:
                C[i, j] = mp.mpf(1) / 2
            else:
                C[i, j] = mp.sin(mp.pi * distance / 2) / (mp.pi * distance)

    eigenvalues, eigenvectors = mp.eigsy(C)
    modular_eigenvalues = [
        mp.log((1 - eigenvalues[k]) / eigenvalues[k]) for k in range(N)
    ]
    h = (
        eigenvectors
        * mp.diag(modular_eigenvalues)
        * eigenvectors.T
    )
    return C, eigenvalues, eigenvectors, h, modular_eigenvalues


def linear_fit(xs, ys):
    """Least-squares fit y = slope*x + intercept, in mp precision."""
    n = len(xs)
    mean_x = mp.fsum(xs) / n
    mean_y = mp.fsum(ys) / n
    denominator = mp.fsum((x - mean_x) ** 2 for x in xs)
    slope = mp.fsum(
        (xs[i] - mean_x) * (ys[i] - mean_y) for i in range(n)
    ) / denominator
    intercept = mean_y - slope * mean_x

    ss_tot = mp.fsum((y - mean_y) ** 2 for y in ys)
    ss_res = mp.fsum(
        (ys[i] - (slope * xs[i] + intercept)) ** 2 for i in range(n)
    )
    r_squared = 1 - ss_res / ss_tot if ss_tot != 0 else mp.mpf(1)
    correlation = (
        mp.fsum(
            (xs[i] - mean_x) * (ys[i] - mean_y) for i in range(n)
        )
        / mp.sqrt(
            mp.fsum((x - mean_x) ** 2 for x in xs)
            * mp.fsum((y - mean_y) ** 2 for y in ys)
        )
    )
    return slope, intercept, r_squared, correlation


def correlation_and_scaled_residual(xs, ys):
    """Pearson correlation and min_A ||y-A*x||_2 / ||y||_2."""
    n = len(xs)
    mean_x = mp.fsum(xs) / n
    mean_y = mp.fsum(ys) / n
    covariance = mp.fsum(
        (xs[i] - mean_x) * (ys[i] - mean_y) for i in range(n)
    )
    variance_x = mp.fsum((x - mean_x) ** 2 for x in xs)
    variance_y = mp.fsum((y - mean_y) ** 2 for y in ys)
    correlation = covariance / mp.sqrt(variance_x * variance_y)

    scale = mp.fsum(xs[i] * ys[i] for i in range(n)) / mp.fsum(
        x * x for x in xs
    )
    residual = mp.sqrt(
        mp.fsum((ys[i] - scale * xs[i]) ** 2 for i in range(n))
        / mp.fsum(y * y for y in ys)
    )
    return correlation, scale, residual


def quadratic_form(vector, matrix):
    N = len(vector)
    return mp.fsum(
        vector[i] * matrix[i, j] * vector[j]
        for i in range(N)
        for j in range(N)
    )


def make_bw_hopping_matrix(N, beta):
    """Nearest-neighbor hopping matrix with bond weights beta[b]."""
    B = mp.matrix(N, N)
    for bond in range(N - 1):
        B[bond, bond + 1] = beta[bond]
        B[bond + 1, bond] = beta[bond]
    return B


def smooth_test_vectors(N):
    """Smooth boundary-vanishing probes for the quadratic-form comparison."""
    length = mp.mpf(N)
    positions = [mp.mpf(i) + mp.mpf(1) / 2 for i in range(N)]
    parabola = [
        position * (length - position) / length for position in positions
    ]
    return {
        "sin(pi x/l)": [
            mp.sin(mp.pi * position / length) for position in positions
        ],
        "parabola": parabola,
        "parabola^2": [value * value for value in parabola],
        "parabola modulation": [
            parabola[i]
            * (
                1
                + mp.mpf("0.2")
                * mp.cos(2 * mp.pi * positions[i] / length)
            )
            for i in range(N)
        ],
    }


def main():
    print("=" * 72)
    print("R13 numeric probe: critical free-fermion interval")
    print("=" * 72)
    print("mpmath precision: %d decimal digits" % mp.mp.dps)
    print("interval sizes N: %s" % ", ".join(str(N) for N in NDIMS))

    data = {}
    for N in NDIMS:
        data[N] = build_restricted_matrix(N)

    # ------------------------------------------------------------------
    print("\n[1] Restricted correlation matrix")
    for N in NDIMS:
        C, eigenvalues, _, _, _ = data[N]
        hermitian_error = max(
            abs(C[i, j] - C[j, i])
            for i in range(N)
            for j in range(N)
        )
        pairing_error = max(
            abs(eigenvalues[i] + eigenvalues[N - 1 - i] - 1)
            for i in range(N)
        )
        print(
            "    N=%2d: ||C-C^T||_max=%s, "
            "max|lambda_i+lambda_{N-1-i}-1|=%s, "
            "lambda_min=%s, lambda_max=%s"
            % (
                N,
                fmt(hermitian_error),
                fmt(pairing_error),
                fmt(eigenvalues[0], 20),
                fmt(eigenvalues[N - 1], 20),
            )
        )
        check(
            "N=%d: C is Hermitian" % N,
            hermitian_error <= HERMITIAN_TOL,
            "max residual %s" % fmt(hermitian_error),
        )
        check(
            "N=%d: spec(C)=1-spec(C) by pairing" % N,
            pairing_error <= PHS_TOL,
            "max pairing residual %s" % fmt(pairing_error),
        )
        check(
            "N=%d: all C eigenvalues lie strictly in (0,1)" % N,
            eigenvalues[0] > 0 and eigenvalues[N - 1] < 1,
            "lambda_min=%s, lambda_max=%s"
            % (fmt(eigenvalues[0], 20), fmt(eigenvalues[N - 1], 20)),
        )

    # ------------------------------------------------------------------
    print("\n[2] Checkerboard structure of h=log((1-C)/C)")
    for N in NDIMS:
        h = data[N][3]
        diagonal_max = max(abs(h[i, i]) for i in range(N))
        even_max = max(
            abs(h[i, j])
            for i in range(N)
            for j in range(N)
            if i != j and (i - j) % 2 == 0
        )
        odd_min = min(
            abs(h[i, j])
            for i in range(N)
            for j in range(N)
            if (i - j) % 2 == 1
        )
        print(
            "    N=%2d: max|diag|=%s, max even-distance |h|=%s, "
            "min odd-distance |h|=%s"
            % (
                N,
                fmt(diagonal_max),
                fmt(even_max),
                fmt(odd_min),
            )
        )
        check(
            "N=%d: even-distance h matrix elements vanish to 1e-40" % N,
            diagonal_max <= CHECKERBOARD_TOL
            and even_max <= CHECKERBOARD_TOL,
            "diag=%s, even=%s" % (fmt(diagonal_max), fmt(even_max)),
        )
        check(
            "N=%d: odd-distance h matrix elements are nonzero" % N,
            odd_min > CHECKERBOARD_TOL,
            "min odd-distance |h|=%s" % fmt(odd_min),
        )

    # ------------------------------------------------------------------
    print("\n[3] Center nearest-neighbor growth")
    center_signed = []
    center_magnitude = []
    for N in NDIMS:
        h = data[N][3]
        center = N // 2
        signed = h[center, center + 1]
        magnitude = abs(signed)
        center_signed.append(signed)
        center_magnitude.append(magnitude)
        print(
            "    N=%2d: c=%2d, h_{c,c+1}=%s, |h_{c,c+1}|=%s"
            % (N, center, fmt(signed), fmt(magnitude))
        )

    magnitudes = [center_magnitude[i] for i in range(len(NDIMS))]
    slope, intercept, r_squared, _ = linear_fit(
        [mp.mpf(N) for N in NDIMS], magnitudes
    )
    signed_slope, signed_intercept, _, _ = linear_fit(
        [mp.mpf(N) for N in NDIMS], center_signed
    )
    print(
        "    magnitude fit:  |h_c| = %s N + (%s), R^2=%s"
        % (fmt(slope), fmt(intercept), fmt(r_squared))
    )
    print(
        "    signed fit:     h_c  = %s N + (%s)"
        % (fmt(signed_slope), fmt(signed_intercept))
    )
    check(
        "center nearest-neighbor magnitude grows linearly in N",
        mp.mpf("0.75") < slope < mp.mpf("0.95")
        and mp.mpf("-1.2") < intercept < mp.mpf("-0.4")
        and r_squared > mp.mpf("0.99"),
        "slope=%s, intercept=%s, R^2=%s"
        % (fmt(slope), fmt(intercept), fmt(r_squared)),
    )

    # ------------------------------------------------------------------
    print("\n[4] Nearest-neighbor profile versus beta(x)=x(l-x)/l")
    trim = 1  # remove one endpoint bond on each side
    profile_correlations = []
    profile_residuals = []
    for N in NDIMS:
        h = data[N][3]
        length = mp.mpf(N)
        weights = [abs(h[i, i + 1]) for i in range(N - 1)]
        beta = [
            (mp.mpf(i) + mp.mpf(1) / 2)
            * (length - mp.mpf(i) - mp.mpf(1) / 2)
            / length
            for i in range(N - 1)
        ]
        used_weights = weights[trim : len(weights) - trim]
        used_beta = beta[trim : len(beta) - trim]
        correlation, scale, residual = correlation_and_scaled_residual(
            used_beta, used_weights
        )
        profile_correlations.append(correlation)
        profile_residuals.append(residual)
        print(
            "    N=%2d: corr=%s, best-scale A=%s, "
            "||w-A beta||_2/||w||_2=%s"
            % (
                N,
                fmt(correlation),
                fmt(scale),
                fmt(residual),
            )
        )
    check(
        "nearest-neighbor profile correlates with beta after endpoint trim",
        min(profile_correlations) > mp.mpf("0.70"),
        "minimum correlation=%s" % fmt(min(profile_correlations)),
    )
    check(
        "normalized profile residual decreases with N",
        all(
            profile_residuals[i + 1] < profile_residuals[i]
            for i in range(len(profile_residuals) - 1)
        ),
        "residuals %s"
        % " -> ".join(fmt(value, 8) for value in profile_residuals),
    )

    # ------------------------------------------------------------------
    print("\n[5] Nearest-neighbor BW quadratic-form probe")
    print(
        "    B_N is the hopping matrix with "
        "B_{b,b+1}=B_{b+1,b}=beta_b and no diagonal."
    )
    form_errors = {}
    for N in NDIMS:
        h = data[N][3]
        W = -h
        length = mp.mpf(N)
        c = N // 2
        beta = [
            (mp.mpf(i) + mp.mpf(1) / 2)
            * (length - mp.mpf(i) - mp.mpf(1) / 2)
            / length
            for i in range(N - 1)
        ]
        B = make_bw_hopping_matrix(N, beta)
        alpha = W[c, c + 1] / B[c, c + 1]
        print("    N=%2d: central-bond scale alpha_N=%s" % (N, fmt(alpha)))

        for name, vector in smooth_test_vectors(N).items():
            q_h = quadratic_form(vector, W)
            q_b = quadratic_form(vector, B)
            relative = abs(q_h - alpha * q_b) / abs(q_h)
            form_errors.setdefault(name, []).append(relative)
            print(
                "        %-22s q_h=%s, q_B=%s, relative diff=%s"
                % (name, fmt(q_h, 8), fmt(q_b, 8), fmt(relative, 8))
            )

    monotone_errors = all(
        all(
            values[i + 1] < values[i]
            for i in range(len(values) - 1)
        )
        for values in form_errors.values()
    )
    final_max_error = max(values[-1] for values in form_errors.values())
    check(
        "quadratic-form relative differences decrease with N",
        monotone_errors,
        "max final relative difference=%s" % fmt(final_max_error),
    )
    check(
        "all N=20 quadratic-form relative differences are below 3%",
        final_max_error < mp.mpf("0.03"),
        "max N=20 relative difference=%s" % fmt(final_max_error),
    )

    # ------------------------------------------------------------------
    print("\n[6] Spectral radius scaling")
    spectral_radii = []
    for N in NDIMS:
        modular_eigenvalues = data[N][4]
        spectral_radius = max(abs(value) for value in modular_eigenvalues)
        spectral_radii.append(spectral_radius)
        print(
            "    N=%2d: ||h||_spec=%s, ||h||_spec/N=%s"
            % (
                N,
                fmt(spectral_radius),
                fmt(spectral_radius / N),
            )
        )
    spectral_slope, spectral_intercept, _, spectral_correlation = linear_fit(
        [mp.mpf(N) for N in NDIMS], spectral_radii
    )
    print(
        "    fit: ||h||_spec = %s N + (%s), correlation=%s"
        % (
            fmt(spectral_slope),
            fmt(spectral_intercept),
            fmt(spectral_correlation),
        )
    )
    check(
        "spectral radius scales linearly with N",
        spectral_slope > 1 and spectral_correlation > mp.mpf("0.99"),
        "slope=%s, correlation=%s"
        % (fmt(spectral_slope), fmt(spectral_correlation)),
    )

    # ------------------------------------------------------------------
    print("\n[7] Third-neighbor / nearest-neighbor ratio (report only)")
    for N in NDIMS:
        h = data[N][3]
        center = N // 2 - 1
        ratio = abs(h[center, center + 3]) / abs(h[center, center + 1])
        print(
            "    N=%2d: c=%2d, |h_{c,c+3}|/|h_{c,c+1}|=%s"
            % (N, center, fmt(ratio))
        )
    print(
        "    This is a fixed-spacing lattice ratio.  Whether the tail "
        "contributes in the continuum limit is a document-level question."
    )

    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("ASSERTIONS: %d" % ASSERTIONS)
    print("FAILURES: %d" % FAILURES)
    print("RESULT: %s" % ("PASS" if FAILURES == 0 else "FAIL"))
    print("=" * 72)
    return 0 if FAILURES == 0 else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        FAILURES += 1
        print("  [x] uncaught exception: %r" % (exc,))
        print("\nASSERTIONS: %d" % (ASSERTIONS + 1))
        print("FAILURES: %d" % FAILURES)
        print("RESULT: FAIL")
        sys.exit(1)
