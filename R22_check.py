#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R22_check.py —— 【A1 主符号与 R20 比值口径的分离】独立核验
=========================================================
不访问网络，不修改任何项目文件。

独立实断言：
  F1  文档锚点：R21 的 pi 保留、R20 估计器被排除、L1 仍开放
  F2  ETP 的 h=-H/N、T_N 与精确级数结构
  F3  主符号只保留 m=0，R20 B_N 与 T_N 不同
  F4  独立数值：N=24 近零模比值趋于 -pi，光滑向量不检验 pi
  F5  跨文档状态一致：R0 / STATUS / INDEX / R20 / R21 / R19
  F6  反过度主张：不关闭 L1、不替代 Z-WICK / Z-WEDGE
"""
import io
import os
import sys

import mpmath as mp


HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 50

FAIL = 0
NCHK = 0


def check(name, cond, detail="", level="ind"):
    global FAIL, NCHK
    NCHK += 1
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


DOC = read("R22_principal_symbol_vs_r20_estimator.md")
R0 = read("R0_publication_theorem.md")
R19 = read("R19_L1_upstream_probability_phase_and_missing_boost.md")
R20 = read("R20_A1_verdict_shape_holds_constant_fails.md")
R21 = read("R21_vf_normalization_resolves_the_r20_factor.md")
STATUS = read("STATUS.md")
IDX = read("INDEX.md")


# ======================================================================
head("F1  文档锚点")

check("标题含主符号与 R20 估计器", "主符号" in DOC and "R20" in DOC)
check("保留 R21 的 pi", "\\pi=2\\pi/v_F" in DOC)
check("排除 R20 有限尺寸解释", "不是待补的有限尺寸常数" in DOC)
check("L1 仍标为未关闭", "J1 仍未关闭" in DOC)
check("不新增 Zero 参数", "不新增 Zero 参数" in DOC)


# ======================================================================
head("F2  ETP 的精确结构")

check("h=-H/N 在位", "h=-\\frac{H}{N}" in DOC)
check("T_N 的 ETP 定义在位",
      "(T_N)_{b,b+1}" in DOC and "\\frac{b+1}{N}" in DOC)
check("c0=pi 在位", "\\alpha_m\\beta_m\\big|_{m=0}=\\pi" in DOC)
check("K_N=-N sum 在位移位", "K_N=-N\\sum_{m\\ge0}" in DOC)


# ======================================================================
head("F3  主符号与估计器分离")

check("m=0 的唯一导数贡献在位",
      "S_m=\\delta_{m,0}" in DOC and "m\\ge1" in DOC)
check("R20 B_N 定义在位",
      "(B_N^{\\rm R20})_{b,b+1}" in DOC and "N-1" in DOC)
check("逐点差 O(1/N) 在位", "O(1/N)" in DOC)
check("命题 R22.1 在位", "命题 R22.1" in DOC and "条件证成" in DOC)
check("命题 R22.2 在位", "命题 R22.2" in DOC and "【排除】" in DOC)
check("0.9255 无定理地位", "没有定理地位" in DOC)
check("条件强预解骨架在位", "条件强预解骨架（未证）" in DOC)
check("余项因子分解与 O(N) 在位",
      "R_N=T_N^2q(T_N)B_N" in DOC and "\\lVert R_N\\rVert=O(N)" in DOC)
check("四条具名输入齐备",
      all(x in DOC for x in ["R22-STAG", "R22-UV", "R22-SA", "R22-FOCK"]))
check("禁止误用下半有界 Kato", "不能直接套“下半有界二次型的 Kato 收敛”" in DOC)


# ======================================================================
head("F4  独立数值复核")

N = 24
C = mp.matrix(N, N)
for i in range(N):
    for j in range(N):
        d = i - j
        C[i, j] = mp.mpf(1) / 2 if d == 0 else mp.sin(mp.pi * d / 2) / (mp.pi * d)

ev, evec = mp.eigsy(C)
K = mp.matrix(N, N)
for k in range(N):
    val = mp.log((1 - ev[k]) / ev[k])
    for i in range(N):
        for j in range(N):
            K[i, j] += val * evec[i, k] * evec[j, k]


def qform(M, f):
    return mp.fsum(f[i] * M[i, j] * f[j] for i in range(N) for j in range(N))


T = mp.matrix(N, N)
B = mp.matrix(N, N)
length = mp.mpf(N - 1)
for b in range(N - 1):
    zt = mp.mpf(b + 1) / N * (1 - mp.mpf(b + 1) / N)
    zb = (mp.mpf(b) + mp.mpf(1) / 2) * (length - mp.mpf(b) - mp.mpf(1) / 2) / length
    T[b, b + 1] = T[b + 1, b] = zt
    B[b, b + 1] = B[b + 1, b] = zb

wt, vt = mp.eigsy(T)
idx = min(range(N), key=lambda k: abs(wt[k]))
f_low = [vt[i, idx] for i in range(N)]
ratio_low = qform(K, f_low) / (N * qform(T, f_low))

wb, vb = mp.eigsy(B)
idx_b = min(range(N), key=lambda k: abs(wb[k]))
f_blow = [vb[i, idx_b] for i in range(N)]
ratio_blow = qform(K, f_blow) / qform(B, f_blow)

u = [mp.mpf(i) / (N - 1) for i in range(N)]
f_smooth = [mp.sin(mp.pi * x) for x in u]
ratio_ts = qform(K, f_smooth) / (N * qform(T, f_smooth))
ratio_bs = qform(K, f_smooth) / qform(B, f_smooth)

check("ETP T 的最近零模本征值小", abs(wt[idx]) < mp.mpf("0.02"),
      "lambda=%s" % mp.nstr(wt[idx], 12))
check("ETP T 近零模比值趋近 -pi", abs(ratio_low + mp.pi) < mp.mpf("2e-3"),
      "ratio=%s  err=%s" % (mp.nstr(ratio_low, 12), mp.nstr(abs(ratio_low + mp.pi), 8)))
check("光滑向量不是主符号检验", abs(ratio_ts + mp.pi) > mp.mpf("0.2"),
      "ratio=%s" % mp.nstr(ratio_ts, 12))
check("R20 B 的最近零模给出明显不同比值", abs(ratio_blow + mp.pi) > mp.mpf("3"),
      "ratio=%s" % mp.nstr(ratio_blow, 12))
check("R20 B 的光滑向量比值也在非 pi 区", abs(ratio_bs + mp.pi) > mp.mpf("0.4"),
      "ratio=%s" % mp.nstr(ratio_bs, 12))


# ======================================================================
head("F5  跨文档状态一致")

check("R0 范围已扩到 R22",
      any(("R8`–`R%d" % n) in R0 for n in range(22, 40)) or ("R8-R22" in R0))
check("STATUS 已登记 R22", "R22" in STATUS and "2.22" in STATUS)
check("STATUS 已登记条件强预解骨架",
      "R22-STAG" in STATUS and "R22-UV" in STATUS and "R_N=T_N^2q(T_N)B_N" in STATUS)
check("INDEX 已链接 R22 文档",
      "](R22_principal_symbol_vs_r20_estimator.md)" in IDX)
check("INDEX 已链接 R22 核验脚本", "](R22_check.py)" in IDX)
check("R21 已指向 R22 校正",
      "R22" in R21 and ("估计器" in R21 or "口径" in R21))
check("R20 已标出估计器偏差", "R22" in R20 and ("估计器" in R20 or "口径" in R20))
check("R19 的 A1 不再只指向裸数值", "R22" in R19 or "主符号" in R19)


# ======================================================================
head("F6  反过度主张")

check("不关闭 L1", "没有证明四维 L1" in DOC)
check("保留 Z-WICK / Z-WEDGE", "`Z-WICK`" in DOC and "`Z-WEDGE`" in DOC)
check("保留 Z-CORE / Z-TAIL", "`Z-CORE`／`Z-TAIL`" in DOC)
check("条件恢复判定不变", "条件恢复" in DOC or "不关闭" in DOC)
check("四维主符号不可替代楔形", "不能替代楔形几何" in DOC)


# ======================================================================
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (NCHK, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
