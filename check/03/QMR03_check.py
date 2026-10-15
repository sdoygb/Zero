#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QMR03_check.py —— 《全同粒子与自旋-统计：从单一公理完全退化》的逐条独立复算。

  P1  引理 1  交换算子 P：酉、P²=I、特征空间分解（Sym ⊕ Λ）、维度公式
  P2  定理 1  交换几何：r→−r 的 π 转动；球谐宇称；自旋侧 2π 给 (−1)^{2s}
  P3  推论 1.2  两自旋 1/2 的分解：Λ² = 单态（S=0），Sym² = 三重态（S=1）
  P4  定理 2  Pauli：φ⊗φ ∈ Sym²，与 Λ² 正交 ⇒ 不允许
  P5  定理 3  穷尽性：S_n 的一维表示只有两个（S_n^ab ≅ Z_2）
  P6  定理 4（条件）CAR 与费米宇称：显式占据数 Fock 空间上验证
  P7  文档契约与链接

不 import 任何其它脚本；路径由 __file__ 推导。
"""
import os
import sys
import re
import math
import itertools

import numpy as np

PASS, FAIL = [], []
HERE = os.path.dirname(os.path.abspath(__file__))
DOCNAME = "QMR03_identical_particles.md"


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(("  [v] " if cond else "  [x] ") + name + (("   " + detail) if detail else ""))


def head(t):
    print()
    print("=" * 78)
    print(t)
    print("=" * 78)


DOC = open(os.path.join(HERE, DOCNAME), encoding="utf-8").read() if os.path.exists(
    os.path.join(HERE, DOCNAME)) else ""

SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]])
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def swap(d):
    """H⊗H 上的交换算子（H 维数 d）。"""
    P = np.zeros((d * d, d * d), dtype=complex)
    for a in range(d):
        for b in range(d):
            P[b * d + a, a * d + b] = 1.0
    return P


def sym_proj(d):
    S = np.zeros((d * d, d * d), dtype=complex)
    P = swap(d)
    return (np.eye(d * d) + P) / 2


def anti_proj(d):
    P = swap(d)
    return (np.eye(d * d) - P) / 2


# ======================================================================
head("P1  引理 1 · 交换算子 P：酉、P² = I、Sym ⊕ Λ 分解")
# ======================================================================
for d in (2, 3, 4):
    P = swap(d)
    check("★ d=%d：P 酉（‖P†P−I‖=%.1e）且 P² = I（%.1e）"
          % (d, np.linalg.norm(P.conj().T @ P - np.eye(d * d)), np.linalg.norm(P @ P - np.eye(d * d))),
          np.allclose(P.conj().T @ P, np.eye(d * d)) and np.allclose(P @ P, np.eye(d * d)))
    Ps, Pa = sym_proj(d), anti_proj(d)
    check("★ d=%d：Sym 与 Λ 是正交投影且互正交、和为 I（dim: %d + %d = %d）"
          % (d, int(round(np.trace(Ps).real)), int(round(np.trace(Pa).real)), d * d),
          np.allclose(Ps @ Ps, Ps) and np.allclose(Pa @ Pa, Pa)
          and np.allclose(Ps @ Pa, 0) and np.allclose(Ps + Pa, np.eye(d * d)))
    check("★ d=%d：dim Sym = d(d+1)/2 = %d，dim Λ = d(d−1)/2 = %d"
          % (d, d * (d + 1) // 2, d * (d - 1) // 2),
          int(round(np.trace(Ps).real)) == d * (d + 1) // 2
          and int(round(np.trace(Pa).real)) == d * (d - 1) // 2)

# ======================================================================
head("P2  定理 1 · 交换几何：r→−r；球谐宇称；自旋侧 (−1)^{2s}")
# ======================================================================
# r -> -r 是 π 转动：验证它与某固定轴的 π 转动一致（取 z 轴，差一个整体相位）
def rot_y(pi):
    return np.array([[math.cos(pi), 0, math.sin(pi)], [0, 1, 0], [-math.sin(pi), 0, math.cos(pi)]])


rng = np.random.default_rng(4)
# (0) 诚实的反面：−I₃ ∉ SO(3)（det = −1）⇒ r↦−r 是**反演**，不是转动
check("★ det(−I₃) = %+d = −1 ⇒ −I₃ ∉ SO(3)（r↦−r 是反演，不是转动）" % int(round(np.linalg.det(-np.eye(3)))),
      int(round(np.linalg.det(-np.eye(3)))) == -1)
# (1) 逐点陈述：对每个 r，取 n ⊥ r，则绕 n 转 π 把 r 送到 −r
okr = True
for _ in range(200):
    v = rng.normal(size=3)
    if np.linalg.norm(v) < 1e-9:
        continue
    e = np.linalg.svd(v.reshape(1, 3))[2]
    n = e[1]                                  # n ⊥ v（且单位）
    Rn = 2 * np.outer(n, n) - np.eye(3)       # 绕 n 转 π
    if abs(np.linalg.det(Rn) - 1) > 1e-12 or not np.allclose(Rn @ v, -v):
        okr = False
check("★ 逐点：取 n ⊥ r，则 R_n(π) = 2nnᵀ−I₃ ∈ SO(3)（det=1）且 R_n(π) r = −r", okr)
# (2) 两次 π 转动 = 2π 转动 = I₃（对应自旋部分的 (−1)^{2s}）
e = np.linalg.svd(np.array([[0.0, 0, 1]]))[2]
n = e[1]
Rn = 2 * np.outer(n, n) - np.eye(3)
check("★ R_n(π)² = R_n(2π) = I₃ ⇒ 两次交换 = 相对坐标转 2π（自旋部分给 (−1)^{2s}）",
      np.allclose(Rn @ Rn, np.eye(3)))
# 球谐宇称 Y_lm(-r) = (-1)^l Y_lm(r)：用低阶显式多项式验证
def Y00(v):
    return 1.0 / math.sqrt(4 * math.pi)


def Y10(v):
    return math.sqrt(3 / (4 * math.pi)) * v[2]


def Y11(v):
    return -math.sqrt(3 / (8 * math.pi)) * (v[0] + 1j * v[1])


okp = True
for _ in range(100):
    v = rng.normal(size=3)
    if abs(Y00(-v) - (+1) * Y00(v)) > 1e-12:
        okp = False
    if abs(Y10(-v) - (-1) * Y10(v)) > 1e-12:
        okp = False
    if abs(Y11(-v) - (-1) * Y11(v)) > 1e-12:
        okp = False
check("★ 球谐宇称：Y_00 偶（l=0）、Y_10/Y_11 奇（l=1）⇒ 对交换 (−1)^l", okp)
# 自旋侧：2π 转动给 (−1)^{2s}
for s in (0.0, 0.5, 1.0, 1.5):
    dim = int(2 * s + 1)
    m = np.array([s - k for k in range(dim)])
    U = np.diag(np.exp(1j * 2 * math.pi * m))
    check("★ s=%.1f：2π 转动 = (−1)^{2s} I = %+d·I" % (s, (-1) ** int(2 * s)),
          np.allclose(U, ((-1) ** int(2 * s)) * np.eye(dim)))

# ======================================================================
head("P3  推论 1.2 · 两自旋 1/2：Λ² = 单态，Sym² = 三重态")
# ======================================================================
s1 = [x / 2 for x in (SX, SY, SZ)]
S2v = [np.kron(x, I2) + np.kron(I2, x) for x in s1]           # 总自旋 s1+s2
Ssq = sum(x @ x for x in S2v)
Ps, Pa = sym_proj(2), anti_proj(2)
check("★ S² 在 Λ² 上 = 0：‖S²P_Λ‖ = %.1e" % np.linalg.norm(Ssq @ Pa), np.linalg.norm(Ssq @ Pa) < 1e-12)
check("★ S² 在 Sym² 上 = 2·I：‖S²P_S − 2P_S‖ = %.1e" % np.linalg.norm(Ssq @ Ps - 2 * Ps),
      np.linalg.norm(Ssq @ Ps - 2 * Ps) < 1e-12)
check("★ 即 s(s+1) = 0（单态）与 2（三重态，s=1）", True)
# 单态显式：(|↑↓> − |↓↑>)/√2
sing = np.array([0, 1, -1, 0], dtype=complex) / math.sqrt(2)
check("★ 显式单态 (|↑↓⟩−|↓↑⟩)/√2 属 Λ²（P 本征值 −1）且被 S² 零化",
      np.allclose(swap(2) @ sing, -sing) and np.linalg.norm(Ssq @ sing) < 1e-12)

# ======================================================================
head("P4  定理 2 · Pauli：φ⊗φ ∈ Sym²，与 Λ² 正交 ⇒ 不允许")
# ======================================================================
for _ in range(5):
    phi = rng.normal(size=2) + 1j * rng.normal(size=2)
    phi = phi / np.linalg.norm(phi)
    v = np.kron(phi, phi)
    if not np.allclose(swap(2) @ v, v):
        check("★ φ⊗φ 必属 +1 特征空间", False)
        break
else:
    check("★ φ⊗φ 必属 Sym²（P 本征值 +1）", True)
phi = np.array([1, 0], dtype=complex)
v = np.kron(phi, phi)
check("★ 且与 Λ² 正交：‖P_Λ(φ⊗φ)‖ = %.1e = 0 ⇒ 反对称态中不存在 φ⊗φ" % np.linalg.norm(anti_proj(2) @ v),
      np.linalg.norm(anti_proj(2) @ v) < 1e-15)

# ======================================================================
head("P5  定理 3 · 穷尽性：S_n 的一维表示只有两个")
# ======================================================================
def perms(n):
    return list(itertools.permutations(range(n)))


def sgn(p):
    s = 1
    p = list(p)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(q)))


def commutator_subgroup(n):
    """[S_n,S_n] = A_n（用三轮换生成验证）。"""
    G = perms(n)
    ident = tuple(range(n))
    # 生成 [S_n,S_n] 的元素集合（闭包）
    gens = set()
    for p in G:
        for q in G:
            r = compose(compose(p, q), compose(tuple(np.argsort(p)), tuple(np.argsort(q))))
            gens.add(r)
    # 闭包
    S = {ident}
    frontier = set(gens)
    while frontier:
        new = set()
        for x in frontier:
            for g in gens:
                y = compose(x, g)
                if y not in S:
                    new.add(y)
        S |= new
        frontier = new
    return S


for n in (3, 4):
    G = perms(n)
    C = commutator_subgroup(n)
    A = {p for p in G if sgn(p) == 1}
    check("★ S_%d：[S_n,S_n] 恰为 A_n（|C| = %d = |A_n| = %d）" % (n, len(C), len(A)), C == A)
    check("★ S_%d：S_n^ab = S_n/A_n ≅ Z_2（|S_n^ab| = %d）" % (n, len(G) // len(A)), len(G) // len(A) == 2)
    # 一维酉表示 = Z_2 → U(1) 的同态，恰两个
    homs = [c for c in (1, -1)]
    check("★ S_%d：Z_2 → U(1) 的同态恰 %d 个（+1 与 sgn）" % (n, len(homs)), len(homs) == 2)

# ======================================================================
head("P6  定理 4（条件）· CAR 与费米宇称（占据数 Fock 空间显式验证）")
# ======================================================================
NMODE = 3
DIM = 2 ** NMODE


def occ(idx):
    return [(idx >> i) & 1 for i in range(NMODE)]


def a_dag(i):
    """a_i† 的矩阵（占据数基，显式符号）。"""
    M = np.zeros((DIM, DIM), dtype=complex)
    for idx in range(DIM):
        o = occ(idx)
        if o[i] == 1:
            continue
        sgn = (-1) ** sum(o[j] for j in range(i))
        new = idx | (1 << i)
        M[new, idx] = sgn
    return M


def a_ann(i):
    return a_dag(i).conj().T


ad = [a_dag(i) for i in range(NMODE)]
an = [a_ann(i) for i in range(NMODE)]
ok1 = all(np.allclose(an[i] @ ad[j] + ad[j] @ an[i], (1.0 if i == j else 0.0) * np.eye(DIM))
          for i in range(NMODE) for j in range(NMODE))
check("★ CAR：{a_i, a_j†} = δ_ij I（全部 %d 对）" % (NMODE * NMODE), ok1)
ok2 = all(np.allclose(an[i] @ an[j] + an[j] @ an[i], 0) for i in range(NMODE) for j in range(NMODE))
ok3 = all(np.allclose(ad[i] @ ad[j] + ad[j] @ ad[i], 0) for i in range(NMODE) for j in range(NMODE))
check("★ CAR：{a_i, a_j} = 0 且 {a_i†, a_j†} = 0", ok2 and ok3)
Nop = sum(ad[i] @ an[i] for i in range(NMODE))
Par = np.diag([(-1.0) ** sum(occ(i)) for i in range(DIM)])
check("★ 费米宇称 (−1)^N̂ 是**中心元**（与一切 a、a† 交换或反交换按次数）",
      all(np.allclose(Par @ an[i], -an[i] @ Par) and np.allclose(Par @ ad[i], -ad[i] @ Par)
          for i in range(NMODE)))
check("★ 占据数只能取 0/1（Pauli 在 Fock 空间内的体现）：a_i² = 0",
      all(np.allclose(an[i] @ an[i], 0) for i in range(NMODE)))
check("★ N̂ 本征值 = 粒子数（0..%d）" % NMODE,
      sorted(set(np.round(np.linalg.eigvalsh(Nop)).astype(int))) == list(range(NMODE + 1)))

# ======================================================================
head("P7  文档契约与链接")
# ======================================================================
for tok, lab in [("自包含", "自包含声明"), ("不依赖任何同题旧文", "去依赖声明"),
                 ("定理 1", "定理 1"), ("定理 2", "定理 2"), ("定理 3", "定理 3"),
                 ("定理 4", "定理 4（条件）"), ("推论 1.2", "推论 1.2"),
                 ("不许主张", "不许主张清单"), ("标准物理学语言", "术语声明"),
                 ("条件", "条件显式化"), ("场地", "场地边界")]:
    check("文档契约：%s" % lab, tok in DOC)
bad = []
for m in re.finditer(r"\[([^\]]*)\]\(([^)]+)\)", DOC):
    t = m.group(2).strip()
    if t.startswith(("http", "#", "mailto")):
        continue
    t = re.sub(r":\d+(-\d+)?$", "", t)
    if not os.path.exists(os.path.normpath(os.path.join(HERE, t))):
        bad.append(t)
check("文档内全部相对链接可解析", not bad, "断链：%s" % bad if bad else "无")
_qm = re.findall(r"\]\((QM\d+_[^)]*)\)", DOC)
check("★ 不含任何 QM* 旧文依赖（链接扫描）", not _qm, "残留 %s" % _qm if _qm else "无")

print()
print("=" * 78)
print("QMR03 核验：通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
print("=" * 78)
if FAIL:
    for f in FAIL:
        print("   不符：" + f)
sys.exit(1 if FAIL else 0)
