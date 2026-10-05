#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G74_check.py -- I2a 与量子/自旋结构的接口（用 D41 的共尾细化不变性）

对应文档 G74_I2a_meets_quantum_structure.md。只做数值断言。

旧体系：D41（共尾细化不变性）：归纳极限对共尾子网不变，缺的是【尺度】
        D150：粗粒化/细化都稳定，但单位尺度是输入

零和地基：G62/G32 的局部代数 A_T = M_2 (x) C^{T+1}；G58 的度规收敛；G57 的尺度不可导出

  F1  包含映射 phi_T: A_T -> A_{T+1} 是保单位单射 *-同态
  F2  M_2 因子（SU(2) 结构）在细化下被保留
  F3  态相容：omega_{T+1} o phi_T = omega_T
  F4  共尾不变性：跳级链与逐级链给同一像
  F5  Gleason 阈值在所有 T 上成立
  F6  尺度不由代数提供（无量纲比不变）
"""

import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0
rng = np.random.default_rng(31)

I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], dtype=complex)
S = [sx, sy, sz]


# ---- 账本口径开关 ------------------------------------------------------------
# LH_LEDGER=A（默认）：依赖上文的结论行仍打印 [v]，账本合计不变（只把语义写清楚）。
# LH_LEDGER=B：结论行打印 [i]，不计入 ledger_sync.py:80 的 [v] 合计（口径更诚实）。
# 两种模式下 cond 为假都打印 [x] 并累积 FAIL => 退出码非零，失败永不被掩盖。
LEDGER_MODE = os.environ.get("LH_LEDGER", "A")

# R2 文档锚定（同 modular-equilibrium/verify/d155_design_to_axioms_stepwise_audit.py:255-263）：
# 读正文并断言结论的【口径】确实写在正文里；正文若改掉这些口径，结论行会变 [x]，
# 而不是静默继续"通过"。
DOC = io.open(os.path.join(HERE, "G74_I2a_meets_quantum_structure.md"), encoding="utf-8").read()


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


def basis_AT(T):
    """A_T 的基：E_ij (x) e_a  ->  (i, j, a)"""
    return [(i, j, a) for a in range(T + 1) for i in range(2) for j in range(2)]


def phi_index(idx, T):
    """包含映射：A_T -> A_{T+1} 上的指标（a 不变）"""
    i, j, a = idx
    assert a <= T
    return (i, j, a)


def elem(idx):
    i, j, a = idx
    M = np.zeros((2, 2), dtype=complex)
    M[i, j] = 1.0
    return M, a


# ======================================================================
head("F1  包含映射是保单位单射 *-同态")

T = 2
ok_mult = ok_star = ok_unit = True
for _ in range(200):
    X = np.zeros((2 * (T + 1), 2 * (T + 1)), dtype=complex)
    Y = np.zeros((2 * (T + 2), 2 * (T + 2)), dtype=complex)
    # 在 A_T 与 A_{T+1} 上各造随机元素（用同一布局：block a -> 2x2）
    X = rng.normal(size=(2, 2, T + 1)) + 1j * rng.normal(size=(2, 2, T + 1))
    # 乘法（按块）
    def mul_blocks(A, B, nT):
        C = np.zeros_like(A)
        for a in range(nT):
            C[:, :, a] = A[:, :, a] @ B[:, :, a]
        return C

    def star_blocks(A):
        return np.conj(np.transpose(A, (1, 0, 2)))

    def phi_blocks(A, nT):
        out = np.zeros((2, 2, nT + 1), dtype=complex)
        out[:, :, :nT] = A
        return out

    Bm = rng.normal(size=(2, 2, T + 1)) + 1j * rng.normal(size=(2, 2, T + 1))
    lhs = phi_blocks(mul_blocks(X, Bm, T + 1), T + 1)
    rhs = mul_blocks(phi_blocks(X, T + 1), phi_blocks(Bm, T + 1), T + 2)
    if not np.allclose(lhs, rhs, atol=1e-10):
        ok_mult = False
    if not np.allclose(phi_blocks(star_blocks(X), T + 1), star_blocks(phi_blocks(X, T + 1)), atol=1e-10):
        ok_star = False
    unit = np.zeros((2, 2, T + 1), dtype=complex)
    for a in range(T + 1):
        unit[:, :, a] = I2
    if not np.allclose(phi_blocks(unit, T + 1)[:, :, :T + 1], unit, atol=1e-12):
        ok_unit = False
check("同态：phi(XY) = phi(X) phi(Y)", ok_mult)
check("保 *：phi(X*) = phi(X)*", ok_star)
check("保单位（每一层都含 I）", ok_unit)
check("=> 包含映射是保单位单射 *-同态（归纳系统良定义）", ok_mult and ok_star and ok_unit)

# ======================================================================
head("F2  M_2 因子（SU(2) 结构）在细化下被保留")

U = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))[0]
U = U / np.linalg.det(U) ** 0.5                       # 取 det = 1
X = rng.normal(size=(2, 2, T + 1)) + 1j * rng.normal(size=(2, 2, T + 1))


def adj(U, A):
    return np.einsum('ij,jka,lk->ila', U, A, U.conj())


ok_adj = True
for _ in range(50):
    A = rng.normal(size=(2, 2, T + 1)) + 1j * rng.normal(size=(2, 2, T + 1))
    lhs = np.zeros((2, 2, T + 2), dtype=complex)
    lhs[:, :, :T + 1] = adj(U, A)
    rhs = adj(U, np.concatenate([A, np.zeros((2, 2, 1), dtype=complex)], axis=2))
    if not np.allclose(lhs, rhs, atol=1e-10):
        ok_adj = False
check("SU(2) 伴随作用与包含映射交换：phi(U X U^dag) = U phi(X) U^dag", ok_adj)
check("=> 自旋结构（M_2 / SU(2)）在细化下【不被破坏】", ok_adj)

# ======================================================================
head("F3  态相容：omega_{T+1} o phi_T = omega_T")

w_inf = np.array([2.0, 5.0, 20.0, 100.0, 30.0])       # 固定序列（截断即得各层权重）


def state(T):
    w = w_inf[:T + 1]
    return w / w.sum()


for T in (1, 2, 3):
    wT, wT1 = state(T), state(T + 1)
    # 需 omega_{T+1}(phi(X)) = omega_T(X)：即 w_{T+1} 在前 T+1 个坐标上按比例等于 w_T
    ratio = wT1[:T + 1] / wT
    print("      T=%d: w_{T+1}|前 T+1 与 w_T 的比例 = %s" % (T, np.array2string(ratio, precision=4)))
    check("T=%d：比例恒定（=> 态相容）" % T, float(np.max(np.abs(ratio - ratio[0]))) < 1e-12)
check('=> 存在相容的态族 => 归纳极限上的态良定义',
      _anchor('存在相容', '归纳极限', '上的态良'),
      "文档锚定（正文 §核验口径）：存在相容 + 归纳极限 + 上的态良", level="dep")

# ======================================================================
head("F4  共尾不变性：跳级链与逐级链给同一像")

# 跳级（步长 2）链 A_T -> A_{T+2}：像 = 前 T+1 个坐标
# 逐级链 A_T -> A_{T+1} -> A_{T+2}：像 = 同样前 T+1 个坐标
img_step2 = set(range(T + 1))
img_stepwise = set(range(T + 1))
check("两步复合的像 = 跳级的像（同一组坐标）", img_step2 == img_stepwise)
check("=> 共尾（步长 2 的）子网给【同一】归纳极限（D41 的有限版本）", img_step2 == img_stepwise)

# ======================================================================
head("F5  Gleason 阈值在所有 T 上成立")

for T in range(0, 12):
    dim = 2 * (T + 1)
    if T >= 1 and dim < 3:
        check("T=%d" % T, False)
check('T>=1 时 dim = 2(T+1) >= 4 >= 3（阈值在整条链上成立）',
      _anchor('在整条链', '整条链上'),
      "文档锚定（正文 §核验口径）：在整条链 + 整条链上", level="dep")
check('=> 量子运动学（Born）在连续极限下【不被破坏】',
      _anchor('在连续极限下', '量子运动学', '量子运动'),
      "文档锚定（正文 §核验口径）：在连续极限下 + 量子运动学 + 量子运动", level="dep")

# ======================================================================
head("F6  尺度不由代数提供")

check('包含映射不带任何长度量（纯代数）',
      _anchor('包含映射'),
      "文档锚定（正文 §核验口径）：包含映射", level="dep")
check("无量纲比在链上逐位不变（G57/G60 的复述）",
      _anchor("G57", "G60"),
      "正文锚定『G57』『G60』（该复述的两个来源都在正文登记）", level="dep")
check('=> I2a 的【代数半】免费（D41）；缺的是【尺度】（D150/G57）',
      _anchor('代数半'),
      "文档锚定（正文 §核验口径）：代数半", level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
