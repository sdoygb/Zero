#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R37_check.py -- 核验 KCBS 探针：闭式判据、对照、原生态语境性、阈值。
"""

from __future__ import annotations

import io
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


def head(title: str) -> None:
    print("")
    print("=" * 72)
    print(title)
    print("=" * 72)


def check(name: str, cond: bool, detail: str = "") -> None:
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def read(name: str) -> str:
    with io.open(os.path.join(HERE, name), encoding="utf-8") as handle:
        return handle.read()


DOC = read("R37_kcbs_contextuality.md")
R36 = read("R36_chsh_bell_locality.md")
G29 = read("G29_probability_as_derived_not_postulated.md")
G72 = read("G72_kappa1_from_the_ledger.md")
STATUS = read("STATUS.md")
INDEX = read("INDEX.md")
RES = json.load(io.open(os.path.join(HERE, "R37_kcbs_results.json"), encoding="utf-8"))

NCTX = 2.0
c45 = math.cos(4 * math.pi / 5)
th0 = math.atan(math.sqrt(-1.0 / c45))
vs = [np.array([math.cos(th0),
                math.sin(th0) * math.cos(4 * math.pi * i / 5),
                math.sin(th0) * math.sin(4 * math.pi * i / 5)], dtype=complex) for i in range(5)]
A = sum(np.outer(v, v.conj()) for v in vs)
mu = np.sort(np.linalg.eigvalsh(A))[::-1]


def s_max(spec):
    return float(np.dot(np.sort(np.asarray(spec, dtype=float))[::-1], mu))


# ======================================================================
head("F1  文档结构：判据、正面结论、三个边界")

check("标题与性质为 KCBS ＋ 正面结论",
      "KCBS 探针" in DOC and "是语境的" in DOC and "量子性第一个硬证据" in DOC)
check("闭式判据 (R37-1)(R37-2) 与阈值 (R37-3) 在位",
      "(R37-1)" in DOC and "(R37-2)" in DOC and "(R37-3)" in DOC
      and "von Neumann" in DOC)
check("与 R36 的分工 (R37-4) 在位",
      "(R37-4)" in DOC and "单体统计" in DOC and "多体（Bell）层面" in DOC)
check("三个边界与二值可证伪入口 (R37-5) 在位",
      "(R37-5)" in DOC and "state-dependent" in DOC
      and "0.7236" in DOC and "翻转" in DOC)
check("没有夸大：写明未算出原生谱",
      "没有算出**原生**谱" in DOC)

# ======================================================================
head("F2  独立复算：机制与闭式")

check("相邻正交（误差 < 1e-12）",
      max(abs(np.vdot(vs[i], vs[(i + 1) % 5])) for i in range(5)) < 1e-12)
check("mu(A) = (sqrt5, 1.3819660, 1.3819660) 且和 = 5",
      abs(mu[0] - math.sqrt(5)) < 1e-12
      and abs(mu[1] - 1.381966011250105) < 1e-9
      and abs(mu.sum() - 5.0) < 1e-12)
check("完全混合态给 5/3 < 2", abs(s_max([1 / 3, 1 / 3, 1 / 3]) - 5 / 3) < 1e-12)
check("最优纯态给 sqrt5", abs(s_max([1, 0, 0]) - math.sqrt(5)) < 1e-12)
check("闭式是上确界（3000 随机态抽查不超界）",
      RES["P2_closed_form_check"]["closed_form_upper_bounds_all_samples"] is True)
rng = np.random.default_rng(11)
over = 0.0
for _ in range(400):
    Z = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    Q, _ = np.linalg.qr(Z)
    lam = rng.dirichlet([1.0, 1.0, 1.0])
    rho = Q @ np.diag(lam) @ Q.conj().T
    S = float(sum(np.real(np.vdot(vs[i], rho @ vs[i])) for i in range(5)))
    over = max(over, S - s_max(lam))
check("独立抽查 400 例：随机取向的 S 不超过闭式", over <= 1e-12, "超出=%.2e" % over)

# ======================================================================
head("F3  原生态：语境性判定与阈值")

w = np.array([2.0, 5.0, 20.0, 100.0]); w = w / w.sum()
top3 = np.sort(w[:3] / w[:3].sum())[::-1]
drop = np.sort(w[1:] / w[1:].sum())[::-1]
S_top3, S_drop = s_max(top3), s_max(drop)
check("原生态 top-3 归约：S_max > 2",
      S_top3 > NCTX and abs(S_top3 - 2.01463413439802) < 1e-9, "S=%.4f" % S_top3)
check("原生态去最小块：S_max > 2",
      S_drop > NCTX and abs(S_drop - 2.0652475842498537) < 1e-9, "S=%.4f" % S_drop)
lam_star = (NCTX - mu[1]) / (mu[0] - mu[1])
check("阈值 λ* = (2-μ2)/(μ1-μ2) ≈ 0.7236",
      abs(lam_star - 0.7236067977499785) < 1e-9, "λ*=%.6f" % lam_star)
check("原生 λ1 超过阈值（余量为正）", top3[0] > lam_star,
      "λ1=%.4f 余量=%.4f" % (top3[0], top3[0] - lam_star))
check("λ1 低于阈值时结论翻转为非语境（可证伪性）",
      s_max([0.70, 0.15, 0.15]) < NCTX)

# ======================================================================
head("F4  与结果 JSON 及账本一致")

check("JSON 判定原生态为语境",
      RES["verdict"]["native_state_contextual"] is True
      and all(v > NCTX for v in RES["verdict"]["native_S_values"].values()))
check("JSON 的 μ 与阈值写明", abs(RES["P1_control"]["mu_eigenvalues"][0] - math.sqrt(5)) < 1e-12
      and abs(RES["P4_threshold"]["lambda_star"] - lam_star) < 1e-9)
check("JSON 登记三条 caveat", len(RES["verdict"]["caveats"]) == 3)
check("G29／G72 的权重来源在位", "0.787402" in G29 and "3.912" in G72)
check("R36 的结论未被改动（多体仍无场地）",
      "Bell 局域" in R36)
check("STATUS 已登记 R37",
      "### 2.37" in STATUS
      and "R37_kcbs_contextuality.md" in STATUS
      and "R37_check.py" in STATUS
      and "语境" in STATUS)
check("INDEX 已收录 R37",
      "R37_kcbs_contextuality.md" in INDEX and "R37_check.py" in INDEX)

print("")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
