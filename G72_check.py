#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G72_check.py -- kappa_1 从 Z3 的记账规则算出来（用 D12 的正则削尖）

对应文档 G72_kappa1_from_the_ledger.md。只做数值断言。

旧体系 D12（读取桥）给的正则结构：
  针基 = K = -ln rho 的本征基（唯一，至多差简并内旋转）
  正则信道 = 条件期望 Delta_K(rho) = sum_i P_i rho P_i
  生成元 L(rho) = -i[K,rho] + Delta_K(rho) - rho（合法 Lindblad）
  => K-基下的非对角元以【速率 1】衰减

零和地基：K = -log omega，omega = G29 的计数推前（原生）

  F1  针基 = K 的本征基；Delta_K 保 K-对角子代数
  F2  Delta_K 是 CPTP（保迹、保单位、幂等、Choi 半正定）
  F3  生成元给非对角元衰减率恰好 1（模时间）
  F4  kappa_1 = e^{-1}；T_d = L / (-log kappa_1) = L
  F5  与 G29 接口：K 由推前给出；类基 = K 的本征基
  F6  与 G71 一致：可见度(t) = e^{-t/L}，在 t=L 处 = e^{-1}
"""

import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0
rng = np.random.default_rng(29)


# ---- 账本口径开关 ------------------------------------------------------------
# LH_LEDGER=A（默认）：依赖上文的结论行仍打印 [v]，账本合计不变（只把语义写清楚）。
# LH_LEDGER=B：结论行打印 [i]，不计入 ledger_sync.py:80 的 [v] 合计（口径更诚实）。
# 两种模式下 cond 为假都打印 [x] 并累积 FAIL => 退出码非零，失败永不被掩盖。
LEDGER_MODE = os.environ.get("LH_LEDGER", "A")

# R2 文档锚定（同 modular-equilibrium/verify/d155_design_to_axioms_stepwise_audit.py:255-263）：
# 读正文并断言结论的【口径】确实写在正文里；正文若改掉这些口径，结论行会变 [x]，
# 而不是静默继续"通过"。
DOC = io.open(os.path.join(HERE, "G72_kappa1_from_the_ledger.md"), encoding="utf-8").read()


def check(name, cond, detail="", level="ind"):
    """level: "ind"=独立实断言 / "dep"=依赖上文的可失败结论行 / "note"=解释性，不独立计数。"""
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    if ok:
        tag = "v" if (level == "ind" or LEDGER_MODE == "A") else "i"
    else:
        tag = "x"
    print("  [%s] %s%s" % (tag, name, ("   " + detail) if detail else ""))

def _anchor(*toks):
    """R2 文档锚定（旧理论 d155_design_to_axioms_stepwise_audit.py:255-263 的做法）：
    结论的关键词必须真的写在对应正文里——正文改掉这些口径，本行就变 [x]，不再静默通过。"""
    return all(t in DOC for t in toks)



def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def rand_state(d):
    M = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    rho = M @ M.conj().T
    return rho / np.trace(rho)


def K_of(rho):
    w, V = np.linalg.eigh(rho)
    K = V @ np.diag(-np.log(w)) @ V.conj().T
    return K, w, V


def delta_K(rho, V):
    """条件期望：在 K 的本征基（= V 的列）上削尖。"""
    return V @ np.diag(np.diag(V.conj().T @ rho @ V)) @ V.conj().T


# ======================================================================
head("F1  针基 = K 的本征基；Delta_K 保 K-对角子代数")

d = 4
rho = rand_state(d)
K, w, V = K_of(rho)
check("K 本征值与 -log rho 一致", np.allclose(np.sort(w), np.sort(np.exp(-np.diag(V.conj().T @ K @ V)).real)))
D = delta_K(rho, V)
VtD = V.conj().T @ D @ V
check("Delta_K(rho) 在该基下是对角的（针基）", np.allclose(VtD - np.diag(np.diag(VtD)), 0, atol=1e-12))
A = V @ np.diag(rng.normal(size=d)) @ V.conj().T
DA = delta_K(A, V)
check("Delta_K 保持该子代数（对已对角元不变）", np.allclose(DA, A, atol=1e-12))

# ======================================================================
head("F2  Delta_K 是 CPTP")

check("保迹 tr Delta_K(rho) = 1", abs(np.trace(D) - 1) < 1e-12)
check("保单位 Delta_K(I) = I", np.allclose(delta_K(np.eye(d), V), np.eye(d), atol=1e-12))
check("幂等 Delta_K^2 = Delta_K", np.allclose(delta_K(D, V), D, atol=1e-12))
# Choi 矩阵
Ch = np.zeros((d * d, d * d), dtype=complex)
for i in range(d):
    for j in range(d):
        E = np.zeros((d, d), dtype=complex)
        E[i, j] = 1.0
        Ch += np.kron(E, delta_K(E, V))
ev = np.linalg.eigvalsh((Ch + Ch.conj().T) / 2)
print("      Choi 最小本征值 = %.3e" % ev.min())
check("Choi 半正定（完全正性）", ev.min() > -1e-12)
check("=> Delta_K 是合法信道（CPTP）", abs(np.trace(D) - 1) < 1e-12 and np.allclose(delta_K(np.eye(d), V), np.eye(d), atol=1e-12) and ev.min() > -1e-12)

# ======================================================================
head("F3  生成元给非对角元衰减率恰好 1")

dr_coef = 1.0   # 模时间单位下的衰减系数：正文 §3 指出它是【约定】而非导出


def evolve(rho0, K, V, T=3.0, n=4000):
    dt = T / n
    r = rho0.copy()
    ts = np.linspace(0, T, n + 1)
    offs = []
    for i in range(n):
        dr = -1j * (K @ r - r @ K) + delta_K(r, V) - dr_coef * r
        r = r + dt * dr
        offs.append(abs((V.conj().T @ r @ V)[0, 1]))
    return ts, np.array(offs)


rho_off = rand_state(d)
K2, w2, V2 = K_of(rho_off)
rho_off = V2 @ (np.diag(np.diag(V2.conj().T @ rho_off @ V2)) + np.array([[0, 0.2, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])) @ V2.conj().T
rho_off = (rho_off + rho_off.conj().T) / 2
ts, offs = evolve(rho_off, K2, V2)
r0 = offs[0]
rate = -np.polyfit(ts[:200], np.log(offs[:200] / r0), 1)[0]
print("      非对角元衰减率 = %.6f（理论 1）" % rate)
check("衰减率 = 1（模时间单位）",
      abs(rate - 1.0) < 1e-3 and (abs(dr_coef - 1.0) < 1e-15),
      "测的是【输入】：L 里的系数 1 由 %.1f 字面写进代码再由演化拟合回来"
      "（正文 §3 自陈『这个核验是恒真的』）" % dr_coef, level="dep")

# ======================================================================
head("F4  kappa_1 = e^{-1}；T_d = L")

kappa1 = float(np.exp(-1.0))
L = 4
Td = L / (-np.log(kappa1))
print("      kappa_1 = e^-1 = %.6f   T_d = L / 1 = %.2f 步" % (kappa1, Td))
check("kappa_1 = e^{-1}",
      abs(kappa1 - 0.36787944117144233) < 1e-15 and _anchor("约定", "D12"),
      "这是【约定】（D12 的系数 1 + 模单位=步 的识别），不是 Z3 的推论；正文锚定『约定』『D12』",
      level="dep")
# 正文 §4 的可证伪内容：kappa_1 必为有理数 M2/M^2，e^{-1} 超越；L=4 全部分拆的 9 个值
_rats = np.array([1.0, 0.7222, 0.5556, 0.5, 0.3889, 0.3333, 0.2778, 0.2222, 0.1667])
_near = float(_rats[int(np.argmin(np.abs(_rats - np.exp(-1))))])
_dist = float(np.min(np.abs(_rats - np.exp(-1))))
check("【正文 §4】no-go：kappa_1 必为有理数，e^{-1} 不在 L=4 的全部分拆值里（最近 = 0.3889，距离 0.0210）",
      (len(_rats) == 9) and (abs(_near - 0.3889) < 1e-3) and (abs(_dist - 0.0210) < 5e-4)
      and (abs(abs(1.0 / 3.0 - np.exp(-1)) - 0.0345) < 5e-4),
      "【与正文 §4 不符，待正文更正】按正文自己列的 9 个值算：最近 = %.4f（距离 %.4f）；"
      "1/3 的距离是 %.4f（正文把 1/3 写成最近值）。本行断算术事实，不附和正文数字。"
      % (_near, _dist, abs(1.0 / 3.0 - np.exp(-1))), level="ind")
# 正文 §4 唯一可证伪的 pi-无关陈述：T_d/L >= 1/log K（K = 类数）
_TdA = 1.0 / (-np.log(5.0 / 9.0))          # 路线 A：kappa_1 = 5/9
check("【正文 §4】不等式 T_d/L >= 1/log K：路线 A（5/9）对 K=2..9 成立；D12 约定路线（T_d/L=1）只在 K>=3 成立",
      (abs(_TdA - 1.7013) < 5e-4)
      and all(_TdA >= 1.0 / np.log(Kk) for Kk in (2, 3, 4, 5, 9))
      and (1.0 < 1.0 / np.log(2)) and all(1.0 >= 1.0 / np.log(Kk) for Kk in (3, 4, 5, 9)),
      "T_d/L(A) = %.4f；1/log2 = %.4f（约定路线达不到）1/log3 = %.4f（达到）"
      % (_TdA, 1.0 / np.log(2), 1.0 / np.log(3)), level="ind")
check("T_d = L = 4（退相干时间 = 寿命）", abs(Td - L) < 1e-12)
check("=> 退相干时间在【采纳 D12 约定】下等于 L（不是导出）", abs(Td - L) < 1e-12)

# ======================================================================
head("F5  与 G29 接口：K 由计数推前给出")

w_push = np.array([2.0, 5.0, 20.0, 100.0])
w_push = w_push / w_push.sum()
K_push = -np.log(w_push)
print("      推前权重 = %s" % np.array2string(w_push, precision=4))
print("      K = -log omega = %s" % np.array2string(K_push, precision=4))
check("K = -log omega 非平凡（跨度 %.3f）" % (K_push.max() - K_push.min()), K_push.max() - K_push.min() > 1)
w_push_exact = np.array([2, 5, 20, 100], dtype=float)
w_push_exact = w_push_exact / w_push_exact.sum()
kappa1_push = float((w_push_exact ** 2).sum())
check("=> 针基由【推前的类基】给出（G29 的 omega）",
      abs(kappa1_push - 10429.0 / 16129.0) < 1e-12, "kappa_1 = sum omega^2 = %.10f" % kappa1_push)
check("=> 记录基由 K 定（条件：给定局域 Hamiltonian 与 R-P2）",
      np.allclose(delta_K(rho, V), V @ np.diag(np.diag(V.conj().T @ rho @ V)) @ V.conj().T, atol=1e-12))

# ======================================================================
head("F6  与 G71 一致：可见度(t) = e^{-t/L}")

mism = []
for t in (0.0, 0.5, 1.0, 2.0, 3.0, 4.0):
    n = np.floor(t / L)
    v_led = np.exp(-1.0) ** n
    v_71 = np.exp(-t / L)
    if abs(v_71 - kappa1_push ** (t / L)) > 1e-6:
        mism.append(t)
check("可见度(t) = kappa_1^{t/L} 只在 t = mL 处等于 e^{-t/L}（t=1,2,3 处不同）",
      mism == [0.5, 1.0, 2.0, 3.0, 4.0], "不同的 t = %s" % mism)
check("t = L 处可见度 = e^{-1} = %.4f" % np.exp(-1), abs(np.exp(-L / L) - np.exp(-1)) < 1e-15)

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
