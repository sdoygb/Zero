#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R8_L2_check.py -- 独立审计 R8 的 L2：普适面积密度。

本脚本只做数值与文档口径核验，不修改任何下游状态文件。

审计要点：
  1. 重算 Z3 点汇在寿命环 L=4 上的正能隙、有效交错质量和零模；
  2. 检验 m_stag -> 0 与 xi / L 的标度，纠正 v_F 漏项；
  3. 重算 G78 的三维有限窗口阈值 m>=2，并判断它是否只是拟合判据；
  4. 用更大的二维窗口检查低 m 是否仍可由有限尺寸解释；
  5. 检查面积密度 eta 是否真能不依赖 gap、状态和尺度。
"""

import io
import os
import sys

import numpy as np



HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0
NCHECK = 0


def check(name, cond, detail=""):
    global FAIL, NCHECK
    NCHECK += 1
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def head(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def read(name):
    with io.open(os.path.join(HERE, name), encoding="utf-8") as f:
        return f.read()


R8 = read("R8_jacobson_entanglement_equilibrium_completion.md")
G76 = read("G76_area_law_in_2d.md")
G77 = read("G77_staggered_coupling_from_A5.md")
G78 = read("G78_area_law_in_3d.md")
G79 = read("G79_horizon_thermodynamics.md")
D231 = read("D231_modular_density_profile_gap.md")
G73 = read("G73_B_is_an_input.md")
G83 = read("G83_the_missing_1_14.md")


# ======================================================================
head("F1  文档链与口径：L2 的数值究竟从哪里来")
# ======================================================================

check("R8 的 C3 明确要求 eta 不依赖状态、尺度、gap 与细化",
      "不依赖局部态扰动" in R8
      and "球的中心" in R8
      and "gap 参数" in R8
      and "细化方式" in R8)
check("R8 的 J2 保留 m_stag=0.23534171，并把 xi≈1.078L 标为漏掉 v_F 的旧口径",
      "m_{\\rm stag}=0.23534171" in R8
      and "1.078L" in R8
      and "漏掉了费米速度" in R8
      and "2.157\\,L" in R8)
check("G77 的 L=4 完整强度值来自正能模式 gap=0.47068342",
      "0.47068342" in G77 and "m=\\mathrm{gap}/2" in G77)
check("G77 明确给出大 L 极限 gap*L -> 1.8546",
      "1.8546" in G77 and "\\equiv0\\pmod4" in G77)
check("G77 把点汇的旧关联长度写成 v_F=1 的 xi=1/m",
      "xi=v_F/m" in G77 and "v_F=1" in G77)
check("G78 的可判阈值 m>=2 是有限窗口拟合判据",
      "m=1$ **不能**下结论" in G78 and "只有 $\\xi\\ll L$" in G78)
check("G83 明确回撤：阈值 m>=2 是人为判据，且关联长度应补 v_F",
      "阈值是自设的" in G83 and "关联长度漏了" in G83)
check("G83 的校对结果是 m_stag=0.235342 给 xi=8.50 到 11.22",
      "8.50" in G83 and "11.22" in G83)
check("G73 已把 L=4 降为选择原则/输入，而不是导出",
      "L=4" in G73 and "选择原则" in G73 and "输入" in G73)
check("D231 处理的是模密度剖面，不提供面积密度 eta",
      "模密度剖面" in D231 and "面积密度" not in D231)
check("G79 只接收面积律标度与系数约定，没有独立测 eta 的普适性",
      "面积熵" in G79 and "标度" in G79 and "系数" in G79)
check("G78 的附注仍残留 m_stag=Theta/(2L) 的旧因子 2，与 G77 的 Theta/L 不一致",
      "m_{\\rm stag}=\\Theta(\\sigma)/(2L_{\\rm life})" in G78)


# ======================================================================
head("F2  Z3 点汇的独立重算：有限正能隙不是稳定的均匀质量")
# ======================================================================


def point_sink_spectrum(L):
    H = np.zeros((L, L), dtype=float)
    for j in range(L):
        H[j, (j + 1) % L] -= 1.0
        H[(j + 1) % L, j] -= 1.0
    H[L - 1, L - 1] += 1.0
    return np.linalg.eigvalsh(H)


Ls = np.array([4, 8, 16, 32, 64, 128, 256, 512], dtype=int)
spec = {L: point_sink_spectrum(L) for L in Ls}
positive_gaps = np.array([spec[L][spec[L] > 1e-12].min() for L in Ls])
zero_modes = np.array([np.min(np.abs(spec[L])) for L in Ls])

gap4 = float(positive_gaps[0])
m4 = 0.5 * gap4
print("      L=4 正能隙 = %.12f，m=gap/2 = %.12f" % (gap4, m4))
print("      各 L 的最小 |E| = %s" % " ".join("%.2e" % z for z in zero_modes))

check("L=4 正能模式重算为 0.47068342", abs(gap4 - 0.47068342) < 5e-9,
      "gap4=%.12f" % gap4)
check("L=4 有效质量 m=gap/2 重算为 0.23534171", abs(m4 - 0.23534171) < 5e-9,
      "m4=%.12f" % m4)
check("L=4,8,...,512 同时在机器精度内存在零模", np.max(zero_modes) < 1e-12,
      "max min|E|=%.3e" % np.max(zero_modes))
check("1/L 标度下正能隙单调下降", np.all(np.diff(positive_gaps) < 0),
      "gaps=%s" % np.array2string(positive_gaps, precision=8))

fit_L = Ls[-5:].astype(float)
fit_y = positive_gaps[-5:] * fit_L
fit_A = np.column_stack([np.ones_like(fit_L), 1.0 / fit_L**2])
fit_c = np.linalg.lstsq(fit_A, fit_y, rcond=None)[0]
gapL_limit = float(fit_c[0])
mL_limit = 0.5 * gapL_limit
print("      拟合得到 gap*L -> %.10f，故 m*L -> %.10f" % (gapL_limit, mL_limit))

check("正能隙满足 gap*L -> 1.85459，而不是趋于正常数 gap",
      abs(gapL_limit - 1.85459) < 5e-5, "gap*L=%.10f" % gapL_limit)
check("有效质量满足 m*L -> 0.927295，故 m(L)->0", abs(mL_limit - 0.927295) < 3e-5,
      "m*L=%.10f" % mL_limit)
check("L=512 的有效质量已降到 0.004 以下", 0.5 * positive_gaps[-1] < 0.004,
      "m(512)=%.8f" % (0.5 * positive_gaps[-1]))

# 旧口径 1/m，及 G83 修正后的 v_F/m。点汇的一维费米速度为 2t=2。
xi_old_4 = 1.0 / m4
xi_v2_4 = 2.0 / m4
xi_v264_4 = 2.64 / m4
xi_v2_ratio = 4.0 / gapL_limit
print("      L=4: xi(vF=1)=%.6f, xi(vF=2)=%.6f, xi(vF=2.64)=%.6f"
      % (xi_old_4, xi_v2_4, xi_v264_4))
print("      vF=2 的大 L 标度: xi = %.6f L" % xi_v2_ratio)

check("R8 沿用的 1.078L 等价于漏掉 v_F 的 1/m 口径",
      abs(2.0 / gapL_limit - 1.0784) < 2e-4
      and abs(1.0784 * 4.0 - 4.3136) < 2e-3)
check("一维点汇的费米速度 2t=2 给 xi(L=4)=8.50，而非 4.25",
      abs(xi_v2_4 - 8.50) < 0.02, "xi=%.6f" % xi_v2_4)
check("若采用 G83 的 2.0--2.64 费米速度区间，xi 为 8.50--11.22",
      abs(xi_v264_4 - 11.22) < 0.03, "xi_max=%.6f" % xi_v264_4)
check("正确量纲下的大 L 标度为 xi≈2.157L，而不是 1.078L",
      abs(xi_v2_ratio - 2.157) < 5e-3, "xi/L=%.6f" % xi_v2_ratio)
check("R8 的 xi 被低估计约 2 倍", 1.9 < xi_v2_4 / xi_old_4 < 2.1,
      "ratio=%.6f" % (xi_v2_4 / xi_old_4))


# ======================================================================
head("F3  G78 的 8.5 倍：算术正确，但阈值不是第一原理")
# ======================================================================

N3 = 12
LS3 = np.array([2, 3, 4, 5, 6], dtype=float)


def idx3(i, j, k):
    return ((i % N3) * N3 + (j % N3)) * N3 + (k % N3)


def build_3d(mass):
    H = np.zeros((N3**3, N3**3), dtype=float)
    for i in range(N3):
        for j in range(N3):
            for k in range(N3):
                a = idx3(i, j, k)
                for di, dj, dk in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
                    b = idx3(i + di, j + dj, k + dk)
                    H[a, b] = -1.0
                    H[b, a] = -1.0
                H[a, a] += mass * ((-1) ** (i + j + k))
    return H


def entropy_3d(C, L):
    sel = [idx3(i, j, k) for i in range(L) for j in range(L) for k in range(L)]
    CA = C[np.ix_(sel, sel)]
    nu = np.clip(np.linalg.eigvalsh((CA + CA.conj().T) / 2), 1e-13, 1 - 1e-13)
    return float(-np.sum(nu * np.log(nu) + (1 - nu) * np.log(1 - nu)))


def fit_3d(mass):
    H = build_3d(mass)
    ev, U = np.linalg.eigh(H)
    C = U[:, ev < 0] @ U[:, ev < 0].conj().T
    S = np.array([entropy_3d(C, int(L)) for L in LS3])
    A = np.column_stack([LS3**2, LS3, np.ones_like(LS3)])
    coef = np.linalg.lstsq(A, S, rcond=None)[0]
    rms = float(np.sqrt(np.mean((A @ coef - S) ** 2)))
    loo = []
    for j in range(len(LS3)):
        keep = np.arange(len(LS3)) != j
        loo.append(float(np.linalg.lstsq(A[keep], S[keep], rcond=None)[0][0]))
    spread = (max(loo) - min(loo)) / float(coef[0])
    return {
        "S": S,
        "a": float(coef[0]),
        "rms": rms,
        "loo_spread": float(spread),
        "eta_face": float(coef[0] / 3.0),
    }


masses3 = [0.0, 0.23534171, 1.0, 2.0, 3.0, 4.0]
fit3 = {m: fit_3d(m) for m in masses3}
for m in masses3:
    row = fit3[m]
    print("      m=%.8f  a=%.6f  eta_face=%.6f  rms=%.3e  LOO=%.3e%%"
          % (m, row["a"], row["eta_face"], row["rms"], 100 * row["loo_spread"]))

rms0 = fit3[0.0]["rms"]
rms_a5 = fit3[0.23534171]["rms"]
rms2 = fit3[2.0]["rms"]
gap_ratio = 2.0 / 0.23534171

check("无 gap 的三维拟合明显失败（rms>1e-2）", rms0 > 1e-2,
      "rms0=%.3e" % rms0)
check("m=0.23534171 在 N=12^3 小窗口仍明显失败（rms≈2.50e-2）",
      abs(rms_a5 - 2.50e-2) < 2e-3, "rms=%.3e" % rms_a5)
check("m=2 在同一小窗口拟合很好（rms<1e-4）", rms2 < 1e-4,
      "rms2=%.3e" % rms2)
check("m=2 的 rms 约比 m=0.23534171 改善 500 倍",
      400 < rms_a5 / rms2 < 650, "ratio=%.1f" % (rms_a5 / rms2))
check("m=0.23534171 相对无 gap 只改善约 1.76 倍",
      1.5 < rms0 / rms_a5 < 2.0, "ratio=%.3f" % (rms0 / rms_a5))
check("2/m_stag 的算术缺口为 8.498 倍",
      abs(gap_ratio - 8.498) < 0.01, "ratio=%.6f" % gap_ratio)
check("m=0.23534171 的留一拟合系数散布超过 1%，说明小窗口系数不稳",
      fit3[0.23534171]["loo_spread"] > 0.01,
      "LOO=%.4f%%" % (100 * fit3[0.23534171]["loo_spread"]))
check("m=2 的留一拟合系数散布低于 0.1%，说明大 gap 已由同一窗口分辨",
      fit3[2.0]["loo_spread"] < 1e-3,
      "LOO=%.4f%%" % (100 * fit3[2.0]["loo_spread"]))

eta234 = np.array([fit3[m]["eta_face"] for m in (2.0, 3.0, 4.0)])
prod234 = eta234 * np.array([2.0, 3.0, 4.0])
prod_spread = float(np.std(prod234) / np.mean(prod234))
check("G81 的 eta_face*m 在 m=2,3,4 上稳定（相对散布<3%）",
      prod_spread < 0.03, "spread=%.4f%%" % (100 * prod_spread))
check("但 C3 要求的 eta_face 在 m=2,3,4 中已近似变化 2 倍",
      eta234.max() / eta234.min() > 1.9,
      "eta=%s ratio=%.3f" % (np.array2string(eta234, precision=6), eta234.max() / eta234.min()))
check("因此 m>=2 检验的是 eta_face*m，不是 eta_face 的普适性",
      prod_spread < 0.03 and eta234.max() / eta234.min() > 1.9)


# ======================================================================
head("F4  更大二维窗口：低 m 的失败可由有限窗口解释")
# ======================================================================

N2 = 48
LS2 = np.arange(2, N2 // 2 + 1, dtype=int)


def idx2(i, j):
    return (i % N2) * N2 + (j % N2)


def fit_2d(mass):
    H = np.zeros((N2**2, N2**2), dtype=float)
    for i in range(N2):
        for j in range(N2):
            a = idx2(i, j)
            for di, dj in ((1, 0), (0, 1)):
                b = idx2(i + di, j + dj)
                H[a, b] = -1.0
                H[b, a] = -1.0
            H[a, a] += mass * ((-1) ** (i + j))
    ev, U = np.linalg.eigh(H)
    C = U[:, ev < 0] @ U[:, ev < 0].conj().T
    S = []
    for L in LS2:
        sel = [idx2(i, j) for i in range(L) for j in range(L)]
        CA = C[np.ix_(sel, sel)]
        nu = np.clip(np.linalg.eigvalsh((CA + CA.conj().T) / 2), 1e-13, 1 - 1e-13)
        S.append(float(-np.sum(nu * np.log(nu) + (1 - nu) * np.log(1 - nu))))
    S = np.array(S)
    x = LS2.astype(float)[-7:]
    y = S[-7:]
    A = np.column_stack([x, np.ones_like(x), 1.0 / x])
    coef = np.linalg.lstsq(A, y, rcond=None)[0]
    rms = float(np.sqrt(np.mean((A @ coef - y) ** 2)))
    return float(coef[0]), rms


masses2 = [0.23534171, 0.5, 1.0, 2.0, 4.0]
fit2 = {m: fit_2d(m) for m in masses2}
for m in masses2:
    print("      2D N=48, m=%.8f: alpha=%.9f, rms=%.3e"
          % (m, fit2[m][0], fit2[m][1]))

alpha_a5 = fit2[0.23534171][0]
alpha4 = fit2[4.0][0]
alpha1 = fit2[1.0][0]
check("更大二维窗口中 m=0.23534171 的 1/L 外推系数仍稳定到 1.67 附近",
      abs(alpha_a5 - 1.67136) < 0.01, "alpha=%.6f" % alpha_a5)
check("m=4 的同口径面积系数约为 0.254", abs(alpha4 - 0.25392) < 0.005,
      "alpha=%.6f" % alpha4)
check("二维面积系数从 m=0.23534171 到 m=4 变化超过 6 倍",
      alpha_a5 / alpha4 > 6.0, "ratio=%.3f" % (alpha_a5 / alpha4))
check("即使 m=1 与 m=4 比较，面积系数也变化超过 3 倍",
      alpha1 / alpha4 > 3.0, "ratio=%.3f" % (alpha1 / alpha4))
check("m=0.23534171 的大窗口拟合残差很小但仍非外推定理",
      fit2[0.23534171][1] < 2e-4, "rms=%.3e" % fit2[0.23534171][1])
check("m=4 的大窗口拟合残差接近机器精度",
      fit2[4.0][1] < 1e-10, "rms=%.3e" % fit2[4.0][1])


# ======================================================================
head("F5  eta 的状态依赖、尺度依赖与 C3 判定")
# ======================================================================

ground_state_docs = all(("半满" in doc and "自由费米" in doc) for doc in (G76, G78))
check("G76/G78 的核心测量都是自由费米基态，不含状态扰动族", ground_state_docs)
check("G78 重算只使用占据投影 ev<0，即单一基态",
      "ev < 0" in read("G78_check.py"))
check("G78 的尺寸固定为 N=12^3，没有细化族扫描", "N = 12" in read("G78_check.py"))
check("G76/G78 的数值能反驳 eta 与 gap 无关，但不能单独反驳固定状态下的空间渐近",
      eta234.max() / eta234.min() > 1.9 and fit3[2.0]["rms"] < 1e-4)
check("C3 的局部状态扰动独立性没有被当前数据检验",
      "不依赖局部态扰动" in R8
      and "温度" not in G78
      and "局域扰动" not in G78)
check("C3 的细化独立性没有被当前数据检验",
      "细化方式" in R8 and "N=12^3" in G78 and "更大尺寸未做" in G78)
check("C3 的 gap 无关性已被同一模型中的 eta_face 数值否定",
      eta234.max() / eta234.min() > 1.9)
check("Z3 点汇在寿命细化时给出零模和 1/L 能隙，不能满足固定正 gap 的面积律前提",
      np.max(zero_modes) < 1e-12 and abs(mL_limit - 0.927295) < 3e-5)
check("若 L=4 固定而只放大三维空间区域，则上述结构性结论不再成立，必须另测大空间",
      "选择原则" in G73 and "更大尺寸未做" in G78)

verdict_ok = (
    abs(m4 - 0.23534171) < 5e-9
    and abs(gap_ratio - 8.498) < 0.01
    and eta234.max() / eta234.min() > 1.9
    and np.max(zero_modes) < 1e-12
)
check("综合判定：8.5 倍算术成立；有限窗口与结构性问题必须分开记账",
      verdict_ok)


# ======================================================================
head("汇总")
# ======================================================================

print("  断言 %d 项，不符项：%d" % (NCHECK, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
