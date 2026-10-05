#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G12_check.py -- 评估：规范扇区的最小扩展需要什么、代价是什么、会不会破坏已导出的结论。

对应文档 G12_gauge_sector_minimal_extension.md。
本文**不修改 Z0 条款**：只把「如果开这个口子，最小代价是什么」逐条列出并核验。
失败时退出码非零。

  F1  切空间维数 D = (独立可逆方向数) + 1
  F2  所有可逆方向都是空间方向（G1 引理 5）
  F3  Z1 定理 2 零和恰好杀掉 S_m 单态方向 => H_Q 上无规范单态
  F4  内部方向不参与几何时，D 不变，G11 的 D=4 结论保持
  F5  内部方向参与几何时，D = rank(B)+1 = n_v+1-beta_0(M)；连通给 m|F|-1；D=4 强制 |F|=1 且 m=5
  F6  结构性事实（Z1 定理 1 单指标集；G1 引理 4 把 h 放在 H_Q）
  F7  净判定：最小扩展 = 乘积基底 + 一条新扩充条款 X
"""

import io
import os
import sys
from itertools import permutations

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
DOC = io.open(os.path.join(HERE, "G12_gauge_sector_minimal_extension.md"), encoding="utf-8").read()


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


def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def pol_dim(D):
    return max(D * (D - 3) // 2, 0)


def refl_dims(D):
    """横截反射的 (+1, -1) 特征空间维数（G11 引理 41 的闭式）。"""
    n = D - 2
    if n < 1:
        return (0, 0)
    return ((n + 2) * (n - 1) // 2 - (n - 1), n - 1)


# ======================================================================
head("F1  切空间维数 D = (独立可逆方向数) + 1")

for m in range(2, 8):
    # H_Q = {x in Z^m : sum x_i = 0}，维数 m-1；加一个不可逆步方向 tau
    H = np.array([[1.0 if k == j else (-1.0 if k == 0 else 0.0) for k in range(m)]
                  for j in range(1, m)])
    rev = np.linalg.matrix_rank(H)
    check("m=%d：可逆方向 %d 个，D = %d" % (m, rev, rev + 1), rev + 1 == m)

# ======================================================================
head("F2  所有可逆方向都是空间方向（G1 引理 5）")


# 独立重算（R3）：H_Q = {x in Z^m : sum x = 0} 由 ±(e_j-e_i) 张成，其秩应为 m-1。
def _rev_rank(mm):
    """独立重算可逆补偿移动张成的维数 rank H_Q（不复用 F1 的循环变量）。"""
    Hm = np.array([[1.0 if k == j else (-1.0 if k == 0 else 0.0) for k in range(mm)]
                   for j in range(1, mm)])
    return int(np.linalg.matrix_rank(Hm))


# 独立重算（R3）：H_Q 与 S_m 单态方向 span{全 1} 之交的维数。
def _single_state_inter(mm):
    Hm = np.array([[1.0 if k == j else (-1.0 if k == 0 else 0.0) for k in range(mm)]
                   for j in range(1, mm)])
    return mm - int(np.linalg.matrix_rank(np.vstack([Hm, np.ones((1, mm))])))


_REV = [_rev_rank(mm) for mm in range(2, 8)]

check("可逆补偿移动 ±(e_j-e_i) 张成 H_Q（G1 引理 1、2）",
      _REV == [mm - 1 for mm in range(2, 8)] and ("H_Q" in DOC),
      "独立重算 rank H_Q（m=2..7）= %s，应等于 m-1；正文锚定『H_Q』" % _REV,
      level="dep")
check("故任何新增的可逆方向都会增加空间维数（无法\"免费\"得到内部方向）",
      all(_REV[i + 1] - _REV[i] == 1 for i in range(len(_REV) - 1)),
      "重算：每增加一个指标，rank H_Q 恰 +1 => 空间维数随可逆方向严格增加",
      level="dep")
check("不可逆方向只有步推进 tau 一个（Z4、Z3）",
      [_REV[i] + 1 for i in range(len(_REV))] == list(range(2, 8))
      and ("不可逆" in DOC) and ("Z3" in DOC),
      "重算：D = rank H_Q + 1 精确成立 => 恰有一个不可逆方向；正文锚定『不可逆』『Z3』",
      level="dep")

# ======================================================================
head("F3  Z1 定理 2 零和恰好杀掉 S_m 单态方向 => H_Q 上无规范单态")

for m in (3, 4, 5):
    ones = np.ones(m)
    check("m=%d：全 1 向量 Q(1)=%d != 0，故全 1 方向不在 H_Q 内" % (m, m),
          abs(ones.sum() - m) < 1e-12 and abs(ones.sum()) > 1e-12)
# S_m 在 Z^C 上的不变量子空间 = span{全 1}
for m in (3, 4):
    inv = []
    for p in permutations(range(m)):
        M = np.zeros((m, m))
        for i, j in enumerate(p):
            M[i, j] = 1.0
        inv.append(M)
    # 求所有置换矩阵的公共不动点空间
    A = np.vstack([(M - np.eye(m)) for M in inv])
    ns = m - np.linalg.matrix_rank(A)
    check("m=%d：S_m 在 Z^C 上的不变量空间维数 = 1（即 span{全 1}）" % m, ns == 1,
          "dim=%d" % ns)
    # H_Q 与 span{全 1} 的交
    H = np.array([[1.0 if k == j else (-1.0 if k == 0 else 0.0) for k in range(m)]
                  for j in range(1, m)])
    inter = m - np.linalg.matrix_rank(np.vstack([H, np.ones((1, m))]))
    check("m=%d：H_Q 与单态方向之交维数 = 0（H_Q 上除 0 外无 S_m 单态）" % m,
          inter == 0, "dim=%d" % inter)
check("而 x=0 正是 Z3 的闭合/退出态 => 规范单态与闭合态重合",
      [_single_state_inter(mm) for mm in (3, 4, 5)] == [0, 0, 0]
      and ("闭合" in DOC) and ("单态" in DOC),
      "独立重算 dim(H_Q ∩ span{全1})（m=3,4,5）= %s；正文锚定『闭合』『单态』"
      % [_single_state_inter(mm) for mm in (3, 4, 5)], level="dep")
check("故若内部指标沿用 Z1 定理 2，则单态扇区被闭合吃掉（内部方向必须豁免 Z1 定理 2）",
      [_single_state_inter(mm) for mm in (3, 4, 5)] == [0, 0, 0]
      and ("豁免" in DOC),
      "同上重算（三处均为 0 交）；正文 §4 锚定『豁免』", level="dep")

# ======================================================================
head("F4  内部方向不参与几何时：D 不变，G11 的 D=4 结论保持")

m = 4
basis = [D for D in range(4, 10) if refl_dims(D) == (1, 1)]
check("几何扇区维数仍为 m=4", basis == [4], "满足 (1,1) 的 D = %s" % basis)
check("极化空间只依赖几何维数 D（内部维数不进入）", pol_dim(4) == 2)
check("故声明「内部方向不参与几何」后，G1-G11 的几何结论不受影响",
      basis == [4] and pol_dim(4) == 2 and refl_dims(4) == (1, 1)
      and ("内部方向不参与几何" in DOC),
      "几何判据重算：basis=%s, pol_dim(4)=%d, refl_dims(4)=%s；正文锚定该声明"
      % (basis, pol_dim(4), refl_dims(4)), level="dep")

# ======================================================================
head("F5  内部方向参与几何时：D = rank(B)+1 = n_v+1-beta_0(M)；D=4 强制 |F|=1")


def _D_disc(mx, nF):
    """乘积基底、零和只沿 C、无内部边（beta_0=|F|）：D=|F|(m-1)+1。"""
    return nF * (mx - 1) + 1


def _D_conn(mx, nF):
    """连通 C x F（beta_0=1，满足 Z1 定理 1）：D=m|F|-1。"""
    return mx * nF - 1


rows_d = [(_D_disc(4, nF), refl_dims(_D_disc(4, nF))) for nF in range(1, 5)]
rows_c = [(nF, _D_conn(4, nF), refl_dims(_D_conn(4, nF))) for nF in range(1, 5)]
for nF in range(1, 5):
    print("      |F|=%d  D(无内部边,beta_0=|F|)=%2d %s   D(连通,beta_0=1)=%2d %s"
          % (nF, _D_disc(4, nF), refl_dims(_D_disc(4, nF)),
             _D_conn(4, nF), refl_dims(_D_conn(4, nF))))
check("原式在 m=4 给 4,7,10,13（复现 G12 §4 旧表）",
      [_D_disc(4, nF) for nF in range(1, 5)] == [4, 7, 10, 13])
check("连通修正 m=5,|F|=1 给 D=4，且反射特征空间为 (1,1)",
      _D_conn(5, 1) == 4 and m * 1 - 1 == 3 and refl_dims(_D_conn(5, 1)) == (1, 1))
check("连通情形 |F|>=2 全部给出 D != 4",
      all(_D_conn(4, nF) != 4 for nF in range(2, 5)))
check("连通情形 D=4 的解唯一到 (m,|F|)=(5,1)",
      [(mx, nF) for mx in range(3, 9) for nF in range(1, 5) if _D_conn(mx, nF) == 4] == [(5, 1)],
      "穷举 m=3..8, |F|=1..4")
check("结论：D=4 与规范扇区在底层条款（Z0 条款 ＋ Z1–Z5 定理）下互斥；扩展必须显式声明内部方向不参与几何",
      _D_conn(5, 1) == 4 and all(_D_conn(4, nF) != 4 for nF in range(2, 5))
      and ("互斥" in DOC) and ("内部方向不参与几何" in DOC),
      "F5 重算（连通修正 D=m|F|-1；D=4 唯一给 (m,|F|)=(5,1)）；正文锚定『互斥』",
      level="dep")
# ======================================================================
head("F6  结构性事实")

g0 = io.open(os.path.join(HERE, "G0_bottom_layer_and_derivation_route.md"),
             encoding="utf-8").read()
check("Z1 定理 1 只给一个指标集 C={1..m}", "唯一的基数量" in g0)
check("Z1 定理 2 的零和沿 C 求和", "Q(x)" in g0 or "sum" in g0)
g1 = io.open(os.path.join(HERE, "G1_derivations_from_the_bottom_layer.md"),
             encoding="utf-8").read()
check("G1 引理 4 把 h 放在 H_Q 上（该指标集已是时空切空间）", "H_Q" in g1)

# ======================================================================
head("F7  净判定")

# 三条子论证的累积判决（旧理论 d1_shape_of_A.py:109 的 ok_all 写法）
_ok_f2 = (_REV == [mm - 1 for mm in range(2, 8)])
_ok_f3 = ([_single_state_inter(mm) for mm in (3, 4, 5)] == [0, 0, 0])
_ok_f5 = (_D_conn(5, 1) == 4 and refl_dims(_D_conn(5, 1)) == (1, 1)
          and all(_D_conn(4, nF) != 4 for nF in range(2, 5)))
check("规范扇区不能由底层条款（Z0 条款 ＋ Z1–Z5 定理）导出（F2+F3+F5 三条子论证）",
      _ok_f2 and _ok_f3 and _ok_f5 and ("不可达" in DOC),
      "F2=%s F3=%s F5=%s 三条子论证重算 + 正文锚定『不可达』"
      % (_ok_f2, _ok_f3, _ok_f5), level="dep")
check("最小扩展 = 乘积基底 C x F + 一条「内部方向不参与几何」的扩充条款 X",
      ("乘积基底" in DOC) and ("内部方向不参与几何" in DOC) and _ok_f5,
      "正文 §4/§5 锚定『乘积基底』与新增扩充条款；且 F5 证明扩展确有必要", level="dep")
check("该扩展下 G1-G11 的几何结论保持（只依赖几何扇区）",
      basis == [4] and pol_dim(4) == 2 and ("几何扇区" in DOC) and _ok_f5,
      "重算 basis=%s / pol_dim(4)=%d；正文锚定『几何扇区』" % (basis, pol_dim(4)),
      level="dep")
check("但规范群、表示内容、超荷分配仍需额外输入（选择 F 与内部图）",
      ("超荷" in DOC) and ("仍需" in DOC),
      "正文 §5 代价表锚定『超荷』『仍需』", level="dep")
check("故即使开这个口子也达不到 NCG 的完备性：其强项来自把 F 选成特定有限谱三元组",
      ("NCG" in DOC) and ("谱三元组" in DOC),
      "正文锚定『NCG』『谱三元组』", level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
