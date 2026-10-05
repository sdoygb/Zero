#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R20_check.py —— 【A1 判定：形状成立，常数不成立】的独立核验
============================================================
不访问网络，不修改任何项目文件。独立实现（不导入 R20_A1_probe.py）。

独立实断言：
  F1  文档锚点：A1 目标、结论框、L1 仍开放
  F2  §1  A1 的判定式与 B 的键中点定义
  F3  §2  数值表与跨 f 散布
  F4  §3  外推表、2pi 排除、R15 估计器依赖
  F5  §4  命题 R20.1 三条
  F6  §5  R21／R22 对原文 Z-STRESS-2π 处置的校正
  F7  独立复算：三个估计器、2pi 排除、形状残差、跨 f 散布
  F8  跨文档状态一致：R0 / STATUS / INDEX / R15 / R19
"""
import io
import os
import sys

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 50

FAIL = 0
N = 0


def check(name, cond, detail="", level="ind"):
    global FAIL, N
    N += 1
    ok = bool(cond)
    if not ok:
        FAIL += 1
    tag = "v" if ok else "x"
    if level == "sup" and ok:
        tag = "i"
    print("  [%s] %s%s" % (tag, name, ("   " + detail) if detail else ""))


def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def read(name):
    return io.open(os.path.join(HERE, name), encoding="utf-8").read()


DOC = read("R20_A1_verdict_shape_holds_constant_fails.md")
R0 = read("R0_publication_theorem.md")
R15 = read("R15_zcar_double_cover_and_zstress_scale.md")
R19 = read("R19_L1_upstream_probability_phase_and_missing_boost.md")
STATUS = read("STATUS.md")
IDX = read("INDEX.md")


# ======================================================================
head("F1  文档锚点")

check("标题含‘形状成立，常数不成立’", "形状成立，常数不成立" in DOC)
check("A1 判定式在位", "\\sigma_{\\theta/2\\pi}=\\alpha_\\theta" in DOC)
check("结论框含 2pi=6.283", "2\\pi=6.283" in DOC)
check("L1 仍标为未关闭", "**仍未关闭**" in DOC)


# ======================================================================
head("F2  §1 判定式与 B 的定义")

check("h = log((1-C)/C) 在位", "h=\\log\\frac{1-C}{C}" in DOC)
check("几何权 xi(x)=x(l-x)/l 在位", "\\xi(x)=\\frac{x(l-x)}{l}" in DOC)
check("B 的键中点权在位", "B_{b,b+1}=B_{b+1,b}=\\xi(x_b)" in DOC)
check("两个检验层（近邻 + 二次型）在位",
      "近邻权重" in DOC and "二次型" in DOC)
check("约定歧义约 8% 明写", "8\\%" in DOC or "8 %" in DOC or "约 $8\\%$" in DOC)


# ======================================================================
head("F3  §2 数值表")

check("N=16 的 A_N 记录在位", "3.6074" in DOC)
check("N=80 的 A_N 记录在位", "3.3649" in DOC)
check("N=80 的中心键 rho 记录在位", "3.4182" in DOC)
check("N=80 的二次型记录在位", "3.4613" in DOC)
check("跨光滑测试函数散布 1.03% 在位", "1.03\\%" in DOC or "1.03 %" in DOC)


# ======================================================================
head("F4  §3 外推与估计器")

check("三条外推常数在位",
      all(x in DOC for x in ["3.3048", "3.3802", "3.3943"]))
check("2pi 排除比 1.851 在位", "1.851" in DOC)
check("R15 估计器依赖明写", "估计器依赖" in DOC)
check("候选常数表含 pi 与 pi^2/3", "3.1416" in DOC and "3.2899" in DOC)


# ======================================================================
head("F5  §4 命题 R20.1")

check("命题 R20.1 在位", "命题 R20.1" in DOC)
check("形状成立与常数不成立两条都写",
      "形状成立" in DOC and "常数不成立" in DOC)
check("缺口写成纯常数", "纯常数" in DOC)
check("标为数值证成非定理", "数值证成，非定理" in DOC)


# ======================================================================
head("F6  §5 处置")

check("Z-STRESS-2π 判定式在位",
      "Z-STRESS-2π" in DOC and "1.851\\pm0.05" in DOC)
check("后续校正指出原文处置已被改写",
      "原文的 `Z-STRESS-2π` 判定式" in DOC and "校正后的处置" in DOC)
check("符号层与算子层分开",
      "符号层" in DOC and "算子层" in DOC)
check("排除用 1.851 反推 kappa",
      "排除用 $1.851$ 反推" in DOC)


# ======================================================================
head("F7  独立复算（N=16,24,32）")


def build_h(n):
    C = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            d = i - j
            C[i, j] = mp.mpf(1) / 2 if d == 0 else mp.sin(mp.pi * d / 2) / (mp.pi * d)
    ev, evec = mp.eigsy(C)
    h = mp.matrix(n, n)
    for k in range(n):
        val = mp.log((1 - ev[k]) / ev[k])
        for i in range(n):
            for j in range(n):
                h[i, j] += val * evec[i, k] * evec[j, k]
    return h


def beta_bond(n):
    length = mp.mpf(n - 1)
    return [
        (mp.mpf(b) + mp.mpf(1) / 2) * (length - mp.mpf(b) - mp.mpf(1) / 2) / length
        for b in range(n - 1)
    ]


def quad(M, vec, n):
    return mp.fsum(vec[i] * M[i, j] * vec[j] for i in range(n) for j in range(n))


NS = (16, 24, 32)
A16 = None
last = {}
resid = []
lam_spr = []
for n in NS:
    h = build_h(n)
    w = [abs(h[i, i + 1]) for i in range(n - 1)]
    beta = beta_bond(n)
    uw, ub = w[1:-1], beta[1:-1]
    A = mp.fsum(ub[i] * uw[i] for i in range(len(uw))) / mp.fsum(b * b for b in ub)
    R = mp.sqrt(mp.fsum((uw[i] - A * ub[i]) ** 2 for i in range(len(uw)))) / mp.sqrt(
        mp.fsum(u * u for u in uw))
    rho = w[n // 2] / beta[n // 2]
    l = mp.mpf(n - 1)
    us = [mp.mpf(i) / l for i in range(n)]
    funcs = [[mp.sin(mp.pi * u) for u in us], [u * (1 - u) for u in us], [(u * (1 - u)) ** 2 for u in us]]
    lam = []
    for f in funcs:
        B = mp.matrix(n, n)
        for b in range(n - 1):
            B[b, b + 1] = beta[b]
            B[b + 1, b] = beta[b]
        lam.append(quad(h, f, n) / quad(B, f, n))
    mean_lam = mp.fsum(lam) / len(lam)
    spr = mp.sqrt(mp.fsum((x - mean_lam) ** 2 for x in lam) / len(lam)) / abs(mean_lam)
    resid.append(R)
    lam_spr.append(spr)
    if n == 16:
        A16 = A
    last = {"A": A, "rho": rho, "lam": abs(mean_lam), "R": R}

check("复算 A_N(N=16) 与文档记录 3.6074 一致",
      abs(A16 - mp.mpf("3.6074")) < mp.mpf("0.001"),
      "复算=%s" % mp.nstr(A16, 8))
check("N=32 三个估计器都在 [3.2, 3.7]",
      all(mp.mpf("3.2") < last[k] < mp.mpf("3.7") for k in ("A", "rho", "lam")),
      "A=%s rho=%s |lam|=%s" % (mp.nstr(last["A"], 7), mp.nstr(last["rho"], 7), mp.nstr(last["lam"], 7)))
check("N=32 时 2pi 被排除逾 1.6 倍",
      (2 * mp.pi) / last["lam"] > mp.mpf("1.6"),
      "2pi/|lam|=%s" % mp.nstr((2 * mp.pi) / last["lam"], 6))
check("形状残差随 N 下降",
      resid[2] < resid[1] < resid[0],
      "R: %s" % " -> ".join(mp.nstr(r, 6) for r in resid))
check("跨光滑测试函数散布 < 1.5%",
      max(lam_spr) < mp.mpf("0.015"),
      "spread: %s" % " -> ".join(mp.nstr(s, 6) for s in lam_spr))


# ======================================================================
head("F8  跨文档状态一致")

check("R0 范围已扩到 R20",
      any(("R8`–`R%d" % n) in R0 for n in range(20, 40)) or ("R8-R20" in R0))
check("STATUS 已登记 R20", "R20" in STATUS and "2.20" in STATUS)
check("INDEX 已链接 R20 文档",
      "](R20_A1_verdict_shape_holds_constant_fails.md)" in IDX)
check("INDEX 已链接 R20 核验脚本", "](R20_check.py)" in IDX)
check("R15 已加估计器口径修正",
      "最小二乘" in R15 and ("估计器" in R15 or "口径" in R15))
check("R19 已更新为形状成立／常数差因子",
      "形状成立" in R19 and "1.85" in R19)
check("R19 的 A1 处置仍指向 Z-STRESS-2π", "Z-STRESS-2π" in R19)


# ======================================================================
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
