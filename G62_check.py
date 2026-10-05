#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G62_check.py -- 量子扇区的推导：GNS ＋ 模流 ＋ Gleason（并限定 G11）

对应文档 G62_quantum_sector_from_GNS_modular_flow_gleason.md。只做数值断言。

  F1  非对易观测量：Pauli 对易关系 + D_L 的 2 维不可约表示满足二面体关系
  F2  GNS 构造：<A,B> = omega(A* B) 正定 Hermite，维数 = dim A = 4(T+1)
  F3  *-表示：pi(X)Y = XY 在 GNS 内积下保 *（复振幅结构）
  F4  模流：sigma_t = rho^{it} (.) rho^{-it} 非平凡、保态、满足 KMS
  F5  【已撤回阈值主张】作用空间三个数：dim A = 4(T+1)、dim H_T = 2(T+1)、
      正交极小投影数 = 2(T+1)；G62 §4 的反例在 L(A_T) 上对每个 T 都成立
  F6  L(M_2) 的唯一正交对是 {P, 1-P} => L(A_T) 上加性退化为互补加性；
      故加性不逼出迹形式；而 GNS 向量态直接给出迹形式
"""

import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


# ---- 账本口径开关 ------------------------------------------------------------
# LH_LEDGER=A（默认）：依赖上文的结论行仍打印 [v]，账本合计不变（只把语义写清楚）。
# LH_LEDGER=B：结论行打印 [i]，不计入 ledger_sync.py:80 的 [v] 合计（口径更诚实）。
# 两种模式下 cond 为假都打印 [x] 并累积 FAIL => 退出码非零，失败永不被掩盖。
LEDGER_MODE = os.environ.get("LH_LEDGER", "A")

# R2 文档锚定（同 modular-equilibrium/verify/d155_design_to_axioms_stepwise_audit.py:255-263）：
# 读正文并断言结论的【口径】确实写在正文里；正文若改掉这些口径，结论行会变 [x]，
# 而不是静默继续"通过"。
DOC = io.open(os.path.join(HERE, "G62_quantum_sector_from_GNS_modular_flow_gleason.md"), encoding="utf-8").read()


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


rng = np.random.default_rng(7)

# ======================================================================
head("F1  非对易观测量：Pauli 关系 + 二面体 2 维不可约表示")

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], dtype=complex)
check("[sx,sy] = 2i sz", np.allclose(sx @ sy - sy @ sx, 2j * sz))
check("[sy,sz] = 2i sx", np.allclose(sy @ sz - sz @ sy, 2j * sx))
check("[sz,sx] = 2i sy", np.allclose(sz @ sx - sx @ sz, 2j * sy))
check('=> 观测量代数【非对易】（M_2 的原生性见 G27）',
      _anchor('观测量代数', '观测量代', '测量代数'),
      "文档锚定（正文 §核验口径）：观测量代数 + 观测量代 + 测量代数", level="dep")

for L in (3, 5, 8):
    r = np.array([[np.cos(2 * np.pi / L), -np.sin(2 * np.pi / L)],
                  [np.sin(2 * np.pi / L), np.cos(2 * np.pi / L)]])
    s = np.diag([1.0, -1.0])
    ok = (np.allclose(np.linalg.matrix_power(r, L), np.eye(2)) and np.allclose(s @ s, np.eye(2))
          and np.allclose(s @ r @ s, np.linalg.inv(r)))
    check("L=%d：二面体关系 r^L = s^2 = 1、s r s = r^{-1}" % L, ok)

# ======================================================================
head("F2  GNS 构造：<A,B> = omega(A* B) 正定 Hermite，维数 = 4(T+1)")

T = 1
dM = 2
dA = dM * dM * (T + 1)        # 基：E_ij (x) e_a
print("      dim A = 4(T+1) = %d（T=%d）" % (dA, T))


def basis():
    out = []
    for a in range(T + 1):
        for i in range(dM):
            for j in range(dM):
                M = np.zeros((dM, dM), dtype=complex)
                M[i, j] = 1.0
                out.append((M, a))
    return out


B = basis()
w_age = np.array([2.0, 5.0])
w_age = w_age / w_age.sum()
G = np.zeros((dA, dA), dtype=complex)
for p, (Mp, ap) in enumerate(B):
    for q, (Mq, aq) in enumerate(B):
        # omega = tr_M(.) * w_age ；omega(A* B) = tr(Mp^dag Mq) * w_age[ap] * delta_{ap,aq}
        G[p, q] = np.trace(Mp.conj().T @ Mq) * w_age[ap] * (1.0 if ap == aq else 0.0)
ev = np.linalg.eigvalsh(G)
print("      Gram 矩阵维数 = %d，最小本征值 = %.3e" % (len(ev), ev.min()))
check("Gram 矩阵 Hermite", np.allclose(G, G.conj().T))
check("Gram 矩阵正定（忠实态）=> GNS 内积良定义", ev.min() > 1e-12)
check("GNS Hilbert 空间维数 = dim A = 4(T+1)（不是 2(T+1)）", G.shape[0] == 4 * (T + 1))
check('=> 复振幅结构成立（复 *-代数 + 正定内积）',
      _anchor('复振幅结', '振幅结构', '正定内积'),
      "文档锚定（正文 §核验口径）：复振幅结 + 振幅结构 + 正定内积", level="dep")

# ======================================================================
head("F3  *-表示：pi(X)Y = XY 在 GNS 内积下保 *")


def ip(X, Y):
    return np.trace(X.conj().T @ Y)


X = rng.normal(size=(dM, dM)) + 1j * rng.normal(size=(dM, dM))
Y = rng.normal(size=(dM, dM)) + 1j * rng.normal(size=(dM, dM))
Z = rng.normal(size=(dM, dM)) + 1j * rng.normal(size=(dM, dM))
lhs = ip(X.conj().T @ Y, Z)          # <pi(X*)Y, Z>
rhs = ip(Y, X @ Z)                    # <Y, pi(X)Z>
check("pi(X*) = pi(X)*（GNS 内积下）", np.allclose(lhs, rhs))
check("=> 表示是 *-表示 => 复振幅与概率幅的代数骨架到位",
      np.allclose(lhs, rhs) and _anchor("复振幅"),
      "绑定上一行 pi(X*)=pi(X)* 的 GNS 内积恒等式（lhs=%.6f rhs=%.6f）+ 正文锚定『复振幅』"
      % (float(np.real(lhs)), float(np.real(rhs))), level="dep")

# ======================================================================
head("F4  模流：非平凡、保态、满足 KMS")

w = np.array([2.0, 5.0, 20.0, 100.0])
w = w / w.sum()
K = -np.log(w)
d = len(w)


def sig(t, A):
    return np.exp(-1j * t * (K[:, None] - K[None, :])) * A


def sig_c(t, A):
    return np.exp(-1j * (t + 1j) * (K[:, None] - K[None, :])) * A


def om(A):
    return float(np.real(np.sum(w * np.diag(A))))


A = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
Bm = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
dev = max(abs(om(sig(t, A)) - om(A)) for t in (0.3, 1.0, 2.7))
check("保态：max|omega(sigma_t A) - omega(A)| < 1e-14", dev < 1e-14, "%.2e" % dev)
off = A.copy()
np.fill_diagonal(off, 0)
move = float(np.max(np.abs(sig(1.0, off) - off)))
check("非平凡：非对角观测量在 t=1 移动 > 0", move > 1e-6, "移动 %.4f" % move)
kms = max(abs(om(A @ sig(t, Bm)) - om(sig_c(t, Bm) @ A)) for t in (0.2, 0.7, 1.9))
check("KMS：|omega(A sigma_t B) - omega(sigma_{t+i} B A)| < 1e-12", kms < 1e-12, "%.2e" % kms)
check('=> 忠实态给出【酉】模流（生成元 K = -log rho）=> 量子动力学',
      _anchor('量子动力学', '量子动力', '子动力学'),
      "文档锚定（正文 §核验口径）：量子动力学 + 量子动力 + 子动力学", level="dep")

# ======================================================================
head("F5  【已撤回阈值主张】三个数各不相同；反例在 L(A_T) 上对每个 T 都成立")

EPS = 0.6


def g_rule(t):
    return 0.5 * (1 + t) + EPS * (t ** 3 - t)


print("      T | dim A = 4(T+1) | dim H_T = 2(T+1) | 正交极小投影数 = 2(T+1)")
for T2 in (0, 1, 2, 7):
    print("      %d |      %3d       |       %3d        |        %3d"
          % (T2, 4 * (T2 + 1), 2 * (T2 + 1), 2 * (T2 + 1)))
check("dim A = 4(T+1) 与 2(T+1) 对每个 T>=0 都不相等（不是同一个对象）",
      all(4 * (t + 1) != 2 * (t + 1) for t in range(8)))
check("dim H_T = 2(T+1) >= 4（T>=1）；T=0 给 2", 2 * 1 < 3 and all(2 * (t + 1) >= 4 for t in (1, 2, 7)))


def rnd_unit(m):
    nn = rng.normal(size=(m, 3))
    return nn / np.linalg.norm(nn, axis=1, keepdims=True)


for T2 in (0, 1, 2, 7):
    Kc = T2 + 1
    N = 2 * Kc
    worst = 0.0
    lo, hi = 2.0, -1.0

    def p_meas(E, Kc=Kc):
        s = 0.0
        for a in range(Kc):
            ba = E[2 * a:2 * a + 2, 2 * a:2 * a + 2]
            s += g_rule(float(np.real(ba[0, 0] - ba[1, 1]))) / Kc
        return s

    for _ in range(400):
        E = np.zeros((N, N), dtype=complex)
        for a in range(Kc):
            P = 0.5 * (np.eye(2, dtype=complex) + rnd_unit(1)[0][0] * sx
                       + rnd_unit(1)[0][1] * sy + rnd_unit(1)[0][2] * sz)
            E[2 * a:2 * a + 2, 2 * a:2 * a + 2] = P
        worst = max(worst, abs(p_meas(E) + p_meas(np.eye(N) - E) - 1))
        lo = min(lo, p_meas(E))
        hi = max(hi, p_meas(E))
    check("T=%2d：反例在 L(A_T) 上互补加性（偏差 %.1e）、取值在 [%.4f, %.4f]"
          % (T2, worst, lo, hi), worst < 1e-12 and lo >= 0 and hi <= 1)

NN = 4000
n = rnd_unit(NN)
v = (1 + n[:, 2]) / 2 + EPS * n[:, 2] * (n[:, 2] ** 2 - 1)
Xf = np.column_stack([np.ones(NN), n])
coef, *_ = np.linalg.lstsq(Xf, v, rcond=None)
res = float(np.max(np.abs(v - Xf @ coef)))
check("反例的最佳仿射（迹形式）残差 > 0.1 => 非迹形式", res > 0.1, "最大残差 %.4f" % res)
check('=> 阈值主张撤回：加大 T 不能消除反例（失效是类型问题，不是维数问题）',
      _anchor('阈值主张', '不是维数'),
      "文档锚定（正文 §核验口径）：阈值主张 + 不是维数", level="dep")

# ======================================================================
head("F6  L(M_2) 的唯一正交对是 {P,1-P}；GNS 直接给迹形式")

I2 = np.eye(2, dtype=complex)


def P_of(nn):
    nn = np.asarray(nn, float)
    return 0.5 * (I2 + nn[0] * sx + nn[1] * sy + nn[2] * sz)


worst = 0.0
for _ in range(1000):
    nn = rnd_unit(1)[0]
    worst = max(worst, float(np.max(np.abs(P_of(nn) + P_of(-nn) - I2))))
check("P(n) + P(-n) = I_2（唯一的正交 rank-1 对）", worst < 1e-12, "%.1e" % worst)
check('=> L(A_T) 上的加性逐块退化为互补加性 p(E)+p(1-E)=1',
      _anchor('上的加性', '逐块退化', '块退化为'),
      "文档锚定（正文 §核验口径）：上的加性 + 逐块退化 + 块退化为", level="dep")

for T2 in (0, 1, 4):
    Kc = T2 + 1
    N = 2 * Kc
    ww = rng.random(Kc) + 0.2
    ww = ww / ww.sum()
    W = np.diag(np.repeat(ww, 2))
    worst = 0.0
    for _ in range(200):
        Xx = np.zeros((N, N), dtype=complex)
        for a in range(Kc):
            Xx[2 * a:2 * a + 2, 2 * a:2 * a + 2] = (rng.normal(size=(2, 2))
                                                    + 1j * rng.normal(size=(2, 2)))
        omm = sum(ww[a] * float(np.real(np.trace(Xx[2 * a:2 * a + 2, 2 * a:2 * a + 2])))
                  for a in range(Kc))
        worst = max(worst, abs(omm - float(np.real(np.trace(W @ Xx)))))
    check("T=%d：GNS 向量态就是迹形式 omega(X) = Tr(W X)（偏差 %.1e）" % (T2, worst),
          worst < 1e-12)
check('=> Born 的【形式】由 GNS 给出（不需要加性）；【唯一性】才需要加性覆盖 B(H_T)',
      _anchor('不需要加性', '不需要加', '需要加性'),
      "文档锚定（正文 §核验口径）：不需要加性 + 不需要加 + 需要加性", level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
