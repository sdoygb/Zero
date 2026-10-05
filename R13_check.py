#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R13_check.py -- L1 正路由（临界自由费米链 -> 强预解收敛）的独立核验。

检查对象:
  - R13 主线文档与三份配套文件的边界表述
  - 有限维准自由模 Hamiltonian 的两条结构事实（棋盘、粒子-空穴配对）
  - 有限尺寸定量结构（近邻权重线性增长、谱半径线性增长、抛物线剖面）
  - G79 的 0.865 相关对正则化截断的依赖（本轮更正）

不访问网络，不修改任何项目文件，不依赖其他 R*_check.py。
"""

from __future__ import annotations

import io
import math
import os
import sys

import numpy as np
import mpmath as mp


HERE = os.path.dirname(os.path.abspath(__file__))
MAIN = "R13_L1_strong_resolvent_attempt.md"
EXTRA = ["R13_external_limit_lemmas.md", "R13_refutation_attempt.md", "R13_numeric_probe.py"]

checks: list[tuple[str, bool, str]] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    checks.append((name, bool(condition), detail))


def read_local(name: str) -> str:
    with io.open(os.path.join(HERE, name), encoding="utf-8") as handle:
        return handle.read()


def all_tokens(text: str, tokens: list[str]) -> bool:
    return all(token in text for token in tokens)


def _corr(N: int) -> mp.matrix:
    C = mp.matrix(N, N)
    for a in range(N):
        for b in range(N):
            d = a - b
            C[a, b] = mp.mpf(1) / 2 if d == 0 else mp.sin(mp.pi * d / 2) / (mp.pi * d)
    return C


def corr_eigs(N: int) -> list:
    E = mp.eigsy(_corr(N))[0]
    return [E[i] for i in range(N)]


def modular_h(N: int, w: list) -> mp.matrix:
    Q = mp.eigsy(_corr(N))[1]
    return Q * mp.diag([mp.log((1 - w[i]) / w[i]) for i in range(N)]) * Q.T


def linfit(xs: list, ys: list) -> tuple:
    n = len(xs)
    sx = sum(xs)
    sy = sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(xs[i] * ys[i] for i in range(n))
    syy = sum(y * y for y in ys)
    denom = n * sxx - sx * sx
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    r2 = (n * sxy - sx * sy) ** 2 / (denom * (n * syy - sy * sy))
    return slope, intercept, r2


def pearson(xs: list, ys: list) -> float:
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    cov = sum((xs[i] - mx) * (ys[i] - my) for i in range(n))
    vx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    vy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return cov / (vx * vy)


def g79_ring(corr_region: int, clip: float) -> tuple:
    """N=96 周期环，区间 [32,64)（与 G79_check.py F1 同配置），返回 (中心键, 相关)。"""
    N = 96
    H = np.zeros((N, N))
    for i in range(N):
        j = (i + 1) % N
        H[i, j] = -1.0
        H[j, i] = -1.0
    ev, U = np.linalg.eigh(H)
    C = U[:, ev < 0] @ U[:, ev < 0].conj().T
    sel = list(range(32, 32 + corr_region))
    CA = C[np.ix_(sel, sel)].real
    w, V = np.linalg.eigh((CA + CA.conj().T) / 2)
    wc = np.clip(w, clip, 1 - clip)
    h = V @ np.diag(np.log((1 - wc) / wc)) @ V.conj().T
    bond = np.array([abs(h[i, i + 1]) for i in range(corr_region - 1)])
    xx = np.arange(corr_region - 1, dtype=float) + 0.5
    para = xx * (corr_region - xx)
    cc = float(np.corrcoef(bond, para)[0, 1])
    return float(bond[len(bond) // 2]), cc


def main() -> int:
    check("R13 主线文档存在", os.path.isfile(os.path.join(HERE, MAIN)))
    if not os.path.isfile(os.path.join(HERE, MAIN)):
        for name, ok, detail in checks:
            print("[%s] %s%s" % ("v" if ok else "x", name, ("  " + detail) if detail else ""))
        print("\nASSERTIONS=%d" % len(checks))
        print("FAILURES=%d" % sum(not ok for _, ok, _ in checks))
        print("RESULT=FAIL")
        return 1

    # ------------------------------------------------------------------
    # A. 文件与边界表述
    # ------------------------------------------------------------------
    text = read_local(MAIN)
    check("配套文件齐全", all(os.path.isfile(os.path.join(HERE, n)) for n in EXTRA))
    check("R13 明确 J1 仍未关闭",
          all_tokens(text, ["J1 仍未关闭", "L1 保持", "开放"]))
    check("R13 明确原样命题被排除",
          all_tokens(text, ["原样", "排除", "谱发散"]))
    check("R13 不冒充 Zero 定理",
          "不是 Zero 定理" in text and "识别" in text)
    check("R13 不写已关闭 J1",
          ("已关闭 J1" not in text) and ("L1 已证" not in text)
          and ("关闭了 L1" not in text))
    check("七个缺口标签齐全",
          all_tokens(text, ["Z-CRIT-DER", "Z-SCALE", "Z-HILB", "Z-CORE",
                            "Z-TAIL", "Z-STRESS", "Z-CONF"]))
    check("引用外部三件工具",
          all_tokens(text, ["Eisler–Tonni–Peschel", "Kato", "Brunetti–Guido–Longo"]))
    check("G79 截断更正被登记",
          all_tokens(text, ["0.865", "截断", "0.964"]))
    check("对抗审计登记谱半径线性发散",
          all_tokens(read_local("R13_refutation_attempt.md"),
                     ["谱半径", "线性", "排除"]))
    check("外部引理文件登记抛物线权重归属表",
          all_tokens(read_local("R13_external_limit_lemmas.md"),
                     ["谁真正给出抛物线权重", "Eisler–Tonni–Peschel"]))

    # ------------------------------------------------------------------
    # B. 有限维结构：棋盘 + 粒子-空穴配对
    # ------------------------------------------------------------------
    mp.mp.dps = 50
    N = 12
    w = corr_eigs(N)
    h = modular_h(N, w)
    even_off = max(abs(h[i, j]) for i in range(N) for j in range(N) if (i - j) % 2 == 0)
    odd_off = min(abs(h[i, j]) for i in range(N) for j in range(N) if (i - j) % 2 == 1)
    check("棋盘结构：偶数距离矩阵元严格为零", even_off < mp.mpf("1e-40"),
          "max_even=%.2e" % float(even_off))
    check("棋盘结构：奇数距离非零", odd_off > 0,
          "min_odd=%.2e" % float(odd_off))

    ws = sorted(w, key=lambda x: float(x))
    ph_pair = max(abs(ws[i] + ws[N - 1 - i] - 1) for i in range(N))
    check("粒子-空穴配对 spec(C)=1-spec(C)", ph_pair < mp.mpf("1e-40"),
          "max|lambda+lambda'-1|=%.2e" % float(ph_pair))

    # ------------------------------------------------------------------
    # C. 有限尺寸定量结构
    # ------------------------------------------------------------------
    Ns = [8, 12, 16, 20]
    centers, spec_ratio, shape_corrs = [], [], []
    for n in Ns:
        wn = corr_eigs(n)
        hn = modular_h(n, wn)
        c = n // 2
        centers.append(float(abs(hn[c, c + 1])))
        spec_ratio.append(float(max(abs(mp.log((1 - x) / x)) for x in wn)) / n)
        prof = [float(abs(hn[j, j + 1])) for j in range(n - 1)]
        beta = [float(mp.mpf(j + 1) * (n - 1 - j) / n) for j in range(n - 1)]
        shape_corrs.append(pearson(prof, beta))

    slope, intercept, r2 = linfit([float(n) for n in Ns], centers)
    check("中心近邻权重随 N 线性增长", abs(slope - 0.8691) < 0.02 and r2 > 0.999,
          "slope=%.4f intercept=%.3f R2=%.6f" % (slope, intercept, r2))
    check("谱半径/N 随 N 单调上升（线性发散）",
          all(spec_ratio[i] < spec_ratio[i + 1] for i in range(len(spec_ratio) - 1))
          and spec_ratio[-1] > 1.5,
          "||h||/N: " + " -> ".join("%.3f" % x for x in spec_ratio))
    check("近邻剖面与抛物线剖面高度相关（>0.99）",
          min(shape_corrs) > 0.99,
          "corr: " + ", ".join("%.5f" % c for c in shape_corrs))

    # ------------------------------------------------------------------
    # D. G79 的 0.865 是截断依赖量（本轮更正）
    # ------------------------------------------------------------------
    clips = [1e-9, 1e-12, 1e-14, 1e-15]
    bond_c, corrs = [], []
    for c in clips:
        b, cc = g79_ring(32, c)
        bond_c.append(b)
        corrs.append(cc)
    check("复现 G79 报告值 0.865（clip=1e-9）",
          abs(corrs[0] - 0.8654) < 0.01,
          "cc(1e-9)=%.5f" % corrs[0])
    check("相关随截断放松单调上升",
          all(corrs[i] <= corrs[i + 1] + 1e-6 for i in range(len(corrs) - 1)),
          " -> ".join("%.4f" % c for c in corrs))
    check("截断依赖显著（1e-15 比 1e-9 高 >0.05）",
          corrs[-1] > corrs[0] + 0.05,
          "cc(1e-9)=%.4f cc(1e-15)=%.4f" % (corrs[0], corrs[-1]))
    check("中心键随截断放松显著增大",
          bond_c[-1] > bond_c[0] * 1.4,
          "bond: " + ", ".join("%.3f" % b for b in bond_c))

    failures = sum(not ok for _, ok, _ in checks)
    for name, ok, detail in checks:
        suffix = ("  " + detail) if detail else ""
        print("[%s] %s%s" % ("v" if ok else "x", name, suffix))
    print()
    print("ASSERTIONS=%d" % len(checks))
    print("FAILURES=%d" % failures)
    print("RESULT=%s" % ("PASS" if failures == 0 else "FAIL"))
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
