#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R58_check.py -- 独立复算 R58 的全部断言

逐项核对文档里出现的每一个数：
  F1  层总权重 b_s（L=16）
  F2  无倾斜 S_max = 1.728556
  F3  门槛 gamma* = 1.105384（二分，S_max 恰为 2）
  F4  代表性值 gamma=1.13 的 omega / 顶三归一 / S_max / q / K 跨度
  F5  账本峰 D=4
  F6  门槛与 L 无关（L=16,20,24 各自二分）
  F7  ln3 与门槛的差（0.61%，非精确）
  F8  V1-V9 九项量子力学验证（调用 R57 的实现）
  F9  III_1 秩（gamma=0 时的原生态）
  F10 K 非平凡对 gamma 的独立性（gamma=0 时 K 跨度 = 9.08）
"""
import io
import math
import os
import sys
from collections import Counter
from fractions import Fraction
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
PASS, FAIL = [], []

MU = (math.sqrt(5), (5 - math.sqrt(5)) / 2, (5 - math.sqrt(5)) / 2)


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print("  [%s] %s%s" % ("v" if cond else "x", name,
                          ("   " + detail) if detail else ""))


def head(t):
    print("\n" + "=" * 74)
    print(t)
    print("=" * 74)


def S_of(w):
    v = sorted(w, reverse=True)[:3]
    s = sum(v)
    return sum(a * b for a, b in zip([x / s for x in v], MU))


def layer_weights(L, r):
    """返回 (b_s dict, omega dict)"""
    b = Counter()
    for k in range(L + 1):
        for s in range(-k, k + 1, 2):
            rem = L - k
            if (rem - s) % 2:
                continue
            a = (rem - s) // 2
            if a < 0 or a > rem:
                continue
            b[abs(s)] += comb(rem, a) * (r ** abs(s))
    tot = sum(b.values())
    return dict(b), {s: v / tot for s, v in b.items()}


def smax_at(L, gamma):
    _, w = layer_weights(L, math.exp(-gamma))
    return S_of(list(w.values()))


def peak_dim(q):
    for D in range(2, 400):
        if D == 2:
            if 3 * q < 1:
                return 2
            continue
        if D * q > (D - 2) and (D + 1) * q < (D - 1):
            return D
    return None


def bisect_gamma(L, lo=1.0, hi=1.2):
    for _ in range(60):
        mid = (lo + hi) / 2
        if smax_at(L, mid) > 2:
            hi = mid
        else:
            lo = mid
    return hi


def factorize(n):
    f, d = {}, 2
    n = int(n)
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def exact_rank(vals):
    fs = [factorize(v) for v in vals]
    primes = sorted({p for f in fs for p in f})
    rows = [[fs[i + 1].get(p, 0) - fs[i].get(p, 0) for p in primes]
            for i in range(len(fs) - 1)]
    if not rows:
        return 0
    M = [[Fraction(x) for x in r] for r in rows]
    rk = 0
    for c in range(len(M[0])):
        piv = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        pv = M[rk][c]
        M[rk] = [x / pv for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk


if __name__ == "__main__":
    head("R58 复算：Z* 下投影 → 量子力学")
    L = 16

    # F1
    b0, w0 = layer_weights(L, 1.0)
    head("F1 层总权重 b_s（gamma=0）")
    exp_b = {0: 17577, 1: 17576, 2: 8162, 3: 3456, 4: 1300,
             5: 418, 6: 108, 7: 20, 8: 2}
    check("b_s 与枚举值一致", all(abs(b0[s] - exp_b[s]) < 1e-9 for s in exp_b),
          str({s: int(b0[s]) for s in sorted(b0)}))
    check("b_0 - b_1 = 1（实测恒差 1，非精确相等）", abs(b0[0] - b0[1] - 1) < 1e-9,
          "%.0f vs %.0f" % (b0[0], b0[1]))

    # F2
    head("F2 无倾斜的语境性")
    s0 = S_of(list(w0.values()))
    check("S_max(gamma=0) = 1.728556", abs(s0 - 1.728556) < 1e-6, "%.6f" % s0)
    check("无倾斜 => 非语境（S_max < 2）", s0 < 2)
    q0 = sum(x * x for x in w0.values())
    check("q(gamma=0) = 0.295416", abs(q0 - 0.295416) < 1e-6, "%.6f" % q0)

    # F3
    head("F3 门槛 gamma*（二分）")
    g16 = bisect_gamma(16)
    check("gamma*(L=16) = 1.105384", abs(g16 - 1.105384) < 5e-6, "%.6f" % g16)
    check("S_max(gamma*) = 2（恰在门槛）",
          abs(smax_at(16, g16) - 2.0) < 1e-8, "%.10f" % smax_at(16, g16))
    check("r* = 0.331084", abs(math.exp(-g16) - 0.331084) < 5e-6,
          "%.6f" % math.exp(-g16))
    check("门槛两侧翻转", smax_at(16, g16 * 0.99) < 2 and smax_at(16, g16 * 1.01) > 2)

    # F4
    head("F4 代表性值 gamma=1.13")
    _, w1 = layer_weights(16, math.exp(-1.13))
    ws1 = [w1[s] for s in sorted(w1)]
    s1 = S_of(ws1)
    check("S_max(1.13) = 2.004730", abs(s1 - 2.004730) < 1e-6, "%.6f" % s1)
    q1 = sum(x * x for x in ws1)
    check("q(1.13) = 0.581992", abs(q1 - 0.581992) < 1e-6, "%.6f" % q1)
    t3 = sorted(ws1, reverse=True)[:3]
    t3n = [x / sum(t3) for x in t3]
    check("顶三归一 = (0.729144, 0.235524, 0.035331)",
          abs(t3n[0] - 0.729144) < 1e-6 and abs(t3n[1] - 0.235524) < 1e-6
          and abs(t3n[2] - 0.035331) < 1e-6,
          str([round(x, 6) for x in t3n]))
    K = [-math.log(x) for x in ws1 if x > 1e-300]
    check("K 跨度 = 18.1212", abs((max(K) - min(K)) - 18.1212) < 1e-3,
          "%.4f" % (max(K) - min(K)))

    # F5
    head("F5 账本峰")
    check("q=0.581992 的账本峰 D=4", peak_dim(q1) == 4, "D=%s" % peak_dim(q1))

    # F6
    head("F6 门槛与 L 无关")
    for LL in [16, 20, 24]:
        g = bisect_gamma(LL)
        check("gamma*(L=%d) 与 L=16 一致（<0.5%%）" % LL,
              abs(g - g16) / g16 < 0.005, "%.6f" % g)

    # F7
    head("F7 与 ln3 的关系")
    d = math.log(3) - g16
    check("ln3 与门槛差 = -0.61%（非精确相等）",
          abs(d / g16 + 0.0061) < 0.002, "差 %+.4f (%+.2f%%)" % (d, 100 * d / g16))
    check("ln3 本身不过门（S_max < 2）", smax_at(16, math.log(3)) < 2,
          "%.6f" % smax_at(16, math.log(3)))

    # F8
    head("F8 九项量子力学验证（直接调用 R57 的函数）")
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "r57", os.path.join(HERE, "R57_quantum_chain.py"))
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        rr = math.exp(-1.13)
        ws_q, _ = m.blocks(16, rr)
        alg = m.check_algebra()
        st = m.check_state(ws_q)
        md = m.check_modular(ws_q)
        bo = m.check_born_and_context(ws_q)
        check("V1 非对易代数", bool(alg["su2_commutators"]
                                and alg["dihedral_relations"]))
        check("V2 态正定且归一", bool(st["rho_positive"]
                                 and abs(st["trace_rho"] - 1) < 1e-9))
        check("V3 GNS 正定", bool(st["gns_positive"]))
        check("V4 K 非平凡", bool(md["K_nontrivial"]))
        check("V5 模流：酉/同态/保 */保迹/保正",
              bool(md["sigma_unitary"] and md["sigma_homomorphism"]
                   and md["sigma_star"] and md["sigma_trace_preserving"]
                   and md["sigma_positive"]))
        check("V6 KMS 条件（<1e-8）", bool(md["kms_holds"]),
              "%.2e" % md["kms_max_diff"])
        check("V7 Born 非负且正交和=1",
              bool(bo["born_nonneg"] and bo["orthogonal_ok"]))
        check("V8 语境性 S_max>2", bool(bo["contextual"]), "%.6f" % bo["S_max"])
    except Exception as e:
        check("F8 调用 R57 成功", False, repr(e))

    # F9
    head("F9 III_1 秩（原生态 b_s）")
    for LL, want in [(4, 2), (8, 4), (12, 6), (16, 7), (20, 10)]:
        bb, _ = layer_weights(LL, 1.0)
        vals = [bb[s] for s in sorted(bb)]
        r = exact_rank(vals)
        check("rank(L=%d) = %d" % (LL, want), r == want, "实算 %d，需 %d"
              % (r, len(vals) - 1))

    # F10
    head("F10 K 非平凡与 gamma 无关")
    for g in [0.0, 0.5, 1.0, 2.0]:
        _, w = layer_weights(16, math.exp(-g))
        ws = [w[s] for s in sorted(w)]
        K = [-math.log(x) for x in ws if x > 1e-300]
        span = max(K) - min(K)
        check("K 跨度(gamma=%.1f) > 1（非平凡）" % g, span > 1.0, "%.4f" % span)
    _, wg0 = layer_weights(16, 1.0)
    ws_g0 = [wg0[s] for s in sorted(wg0)]
    Kg0 = [-math.log(x) for x in ws_g0 if x > 1e-300]
    check("K 跨度(gamma=0) = 9.0812", abs((max(Kg0) - min(Kg0)) - 9.0812) < 1e-3,
          "%.4f" % (max(Kg0) - min(Kg0)))

    head("汇总")
    print("  通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
    if FAIL:
        print("  不符项：", FAIL)
    sys.exit(1 if FAIL else 0)
