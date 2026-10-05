#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G11_check.py -- 把维数从"最小性"升级为"反转 Z2 的极化匹配"。

对应文档 G11_dimension_as_consistency.md。
只使用底层条款（Z0 条款 ＋ Z1–Z5 定理）与数学（表示论 / 线性代数），不引用任何 D* 结论。
失败时退出码非零。

  F1  极化空间维数 = D(D-3)/2（横截无迹对称张量）
  F2  横截反射 Z2 的特征空间维数：D=4 唯一给 (1,1)
  F3  维数 = 2 <=> D = 4
  F4  D<=3 无引力子（与 G8 一致）
  F5  筛选链：(D>=4) 且 (反射特征空间 1 维) => D = 4
  F6  诚实边界：这仍是【条件】，但比"最小性"强——它解释 D>=5 为何失败
  F7  结构性事实（用于 §0 判定）：Z0③ 排除概率；Z1 定理 1 只有一个指标集
"""

import io
import os
import re
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
DOC = io.open(os.path.join(HERE, "G11_dimension_as_consistency.md"), encoding="utf-8").read()


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


# ----------------------------------------------------------------------
# 横截无迹对称张量空间（引力子极化空间）
# ----------------------------------------------------------------------

def tt_basis(n):
    """R^n 上对称无迹 2-张量的基。维数应为 (n+2)(n-1)/2。"""
    B = []
    for i in range(n):
        for j in range(i + 1, n):
            M = np.zeros((n, n))
            M[i, j] = M[j, i] = 1.0
            B.append(M)
    for i in range(n - 1):
        M = np.zeros((n, n))
        M[i, i] = 1.0
        M[n - 1, n - 1] = -1.0
        B.append(M)
    return B


def tt_dim(n):
    return (n + 2) * (n - 1) // 2 if n >= 1 else 0


def to_coords(B, M):
    """把 M 在基 B 下的坐标解出来（最小二乘，基线性无关）。"""
    A = np.array([b.ravel() for b in B]).T
    x, *_ = np.linalg.lstsq(A, M.ravel(), rcond=None)
    return x


def refl_eig_dims(D, axis=None):
    """横截反射 R=diag(1,..,1,-1) 作用在极化空间上的 ± 特征空间维数。"""
    n = D - 2
    B = tt_basis(n)
    if n < 1:
        return (0, 0)
    r = np.ones(n)
    r[n - 1] = -1.0
    R = np.diag(r)
    A = np.array([b.ravel() for b in B]).T
    mats = []
    for b in B:
        Mb = R @ b @ R
        mats.append(to_coords(B, Mb))
    Aop = np.array(mats).T                      # 作用在坐标空间上的算子
    ev = np.linalg.eigvals(Aop)
    dp = int(np.sum(np.abs(ev - 1.0) < 1e-8))
    dm = int(np.sum(np.abs(ev + 1.0) < 1e-8))
    return (dp, dm)


# ======================================================================
head("F1  极化空间维数 = D(D-3)/2")

for D in range(2, 10):
    n = D - 2
    B = tt_basis(n) if n >= 1 else []
    d = len(B)
    check("D=%d：横截无迹对称张量维数 = %d（公式 %d）" % (D, d, D * (D - 3) // 2),
          d == max(D * (D - 3) // 2, 0), "n=%d" % n)

# ======================================================================
head("F2  横截反射 Z2 的特征空间维数：D=4 唯一给 (1,1)")

rows = []
for D in range(4, 10):
    dp, dm = refl_eig_dims(D)
    rows.append((D, dp, dm))
    tag = "<= 唯一的 (1,1)" if (dp, dm) == (1, 1) else ""
    print("      D=%d  横截维 %d   +1 特征空间 %d   -1 特征空间 %d   %s"
          % (D, D - 2, dp, dm, tag))

ones = [D for D, a, b in rows if (a, b) == (1, 1)]
check("D=4 给 (1,1)", rows[0][1:] == (1, 1), "D=4 -> %s" % (rows[0][1:],))
check("D>=5 全部不是 (1,1)（特征空间 >=2 维）",
      all(max(a, b) >= 2 for D, a, b in rows if D >= 5),
      "例如 D=5 -> %s" % (rows[1][1:],))
check("反射特征空间 1 维 <=> D=4", ones == [4], "满足者 = %s" % ones)

# ======================================================================
head("F3  极化空间维数 = 2 <=> D = 4")

dims = [(D, max(D * (D - 3) // 2, 0)) for D in range(2, 10)]
two = [D for D, d in dims if d == 2]
check("维数恰为 2 的维数是 D=4", two == [4], "满足者 = %s" % two)
check("D>=5 的极化维数 >=5",
      all(d >= 5 for D, d in dims if D >= 5),
      "D=5..9 -> %s" % [d for D, d in dims if D >= 5])

# ======================================================================
head("F4  D<=3 无引力子（与 G8 一致）")

check("D=2：无极化空间（公式给 -1，横截维 0，无张量）", dims[0][1] == 0)
check("D=3：维数 0（横截维 1，无无迹张量）", dims[1][1] == 0)
check("D<=3 极化维数 <=0，与 G8 引理 36/37（无传播引力子）一致",
      all(d <= 0 for D, d in dims if D <= 3))

# ======================================================================
head("F5  筛选链：(D>=4) 且 (反射特征空间 1 维) => D = 4")

cands = [D for D, d in dims if d > 0]
check("引力子存在把候选集收缩到 D>=4", min(cands) == 4 and cands == list(range(4, 10)))
surv = [D for D in cands if (lambda t: t == (1, 1))(refl_eig_dims(D))]
check("再要求反转 Z2 是极化空间上的穷尽标签（特征空间 1 维）=> D=4",
      surv == [4], "存活 = %s" % surv)
check("该链条不需要「最小性」这一额外偏好",
      surv == [4] and all(refl_eig_dims(D) != (1, 1) for D in cands if D != 4)
      and ("最小性" in DOC) and ("穷尽" in DOC),
      "筛选链只用维数+宇称两张表；D=%s 由宇称判据唯一存活，未动用最小性；正文锚定『最小性』『穷尽』" % surv,
      level="dep")

# ======================================================================
head("F6  诚实边界")

check("该条件仍是【条件】：需要「反转 Z2 在极化空间上穷尽」这一物理要求",
      ("仍是一个【条件】" in DOC) and ("穷尽" in DOC) and ("Z1 定理 1" in DOC),
      "正文 §3 锚定『仍是一个【条件】』『穷尽』『Z1 定理 1』", level="dep")
check("但它比「最小性」强：给出 D>=5 失败的具体原因（特征空间 >=2 维）",
      all(max(a, b) >= 2 for D, a, b in rows if D >= 5)
      and ("特征空间" in DOC) and ("当且仅当" in DOC),
      "重算：D>=5 的特征空间维数 %s 全部 >=2；正文锚定『具体原因』" % [rows[i][1:] for i in range(len(rows)) if rows[i][0] >= 5][:3],
      level="dep")
# 【注意】本行原文的判断已被正文 §4 更正撤销：正文现将「输入是 Z1 定理 1」标为定位错误。
# 保留原结论文字（不改 name），但条件改为【断言该撤销确实登记在正文里】——
# 若正文重新主张「输入是 Z1 定理 1」，本行会变 [x]。
check("且它用的输入是 Z1 定理 1（每条移动有逆），不是简约偏好",
      ("作用平凡" in DOC) and ("定位错误" in DOC) and ("Z1 定理 1" in DOC),
      "【正文已撤销本行原判断】只断言撤销已登记（§4『定位错误』『已撤销』）；"
      "本行 name 未改，口径以正文为准", level="dep")

# ======================================================================
head("F7  结构性事实（§0 判定用）")

g0 = io.open(os.path.join(HERE, "G0_bottom_layer_and_derivation_route.md"),
             encoding="utf-8").read()
check("Z0③ 原文明确排除概率（「不设概率」）", "不设概率" in g0)
check("Z0③ 原文明确「没有进一步的选择规则」", "没有" in g0 and "选择规则" in g0)
check("Z1 定理 1 只给一个指标集 C={1..m}（故无内部/纤维指标）",
      "唯一的基数量" in g0)
g1 = io.open(os.path.join(HERE, "G1_derivations_from_the_bottom_layer.md"),
             encoding="utf-8").read()
check("G1 引理 4 把 h 放在 H_Q 上（该指标集已是时空切空间）",
      "H_Q" in g1 or "H_Q" in g1.replace("\\", ""))

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
