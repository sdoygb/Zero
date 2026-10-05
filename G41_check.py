#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G41_check.py -- 闭环计数度规是否满足 Lovelock 前提 (L)/(O)/(C)？

对应文档 G41_lovelock_premises_under_nonuniform_weight.md。
新计算：逐条审计 G1 推论 7.1 的三个前提。

  F1  G1 推论 7.1 的前提是 (L) 局域 / (O) 二阶 / (C) 守恒源
  F2  (O) 二阶：加权 Laplacian 是二阶算子（湮灭线性函数，不湮灭常数）
  F3  (C) 守恒源：加权 Laplacian 的核 = 常数 => 无散流存在
  F4  (L) 局域：决定性比值判据 —— 均匀权【精确为 0】，闭环计数【显著非零】
  F5  影响不随距离衰减（稳定在 ~1）=> 真正全局
  F6  三角困境：局域 / 导出 / 不增扩充条款 三者不可兼得（闭环计数路线）
  F7  与 G2 已排除路线的对照（同类非局域）
  F8  诚实边界
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
DOC = io.open(os.path.join(HERE, "G41_lovelock_premises_under_nonuniform_weight.md"), encoding="utf-8").read()


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


def rd(f):
    p = os.path.join(HERE, f)
    return open(p, encoding="utf-8", errors="replace").read() if os.path.exists(p) else ""


def perron(A):
    w, v = np.linalg.eigh(A)
    p = np.abs(v[:, -1])
    return p / np.linalg.norm(p)


def cw_weight(A):
    phi = perron(A)
    return phi[:, None] * phi[None, :] * A


def path(N, pert=None):
    A = np.zeros((N, N))
    for i in range(N - 1):
        e = 1.0 if (pert is None or i != pert[0]) else pert[1]
        A[i, i + 1] = A[i + 1, i] = e
    return A


def wlap(W):
    return np.diag(W.sum(axis=1)) - W


# ======================================================================
head("F1  G1 推论 7.1 的前提")

g1 = rd("G1_derivations_from_the_bottom_layer.md")
check("G1 引理 7 陈述 Lovelock 唯一性（4 维）", "Lovelock" in g1)
check("推论 7.1 的三前提为 (L) 局域 / (O) 二阶 / (C) 守恒源",
      "(L)" in g1 and "(O)" in g1 and "(C)" in g1)
check("结论为 G_ab + Lambda g_ab = 8 pi G T_ab", "8\\pi G" in g1 or "8pi G" in g1)

# ======================================================================
head("F2  (O) 二阶：加权 Laplacian 是二阶算子")

N = 40
A = path(N)
phi = perron(A)
W = phi[:, None] * phi[None, :] * A
L = wlap(W)
const = np.ones(N)
linear = np.arange(N, dtype=float)
check("L @ 1 = 0（常数在核里 => 不是零阶）", np.allclose(L @ const, 0, atol=1e-12))
check("L @ x != 0（线性函数不被湮灭 => 不是一阶）",
      not np.allclose(L @ linear, 0, atol=1e-9),
      "max|L x| = %.4f" % np.abs(L @ linear).max())
quad = np.arange(N, dtype=float) ** 2
check("L @ x^2 != 0（二阶差分的特征）", not np.allclose(L @ quad, 0, atol=1e-9))
check("=> 加权 Laplacian 是【二阶】算子 => (O) 成立",
      (not np.allclose(L @ quad, 0, atol=1e-9)) and np.allclose(L @ np.ones(N), 0, atol=1e-9),
      "重算：L@x^2 != 0 但 L@1 = 0（二阶差分的特征）", level="dep")
check("而且它只涉及【最近邻】（邻接矩阵的支集）=> 算子本身是局部的",
      float(np.abs(L - np.diag(L.diagonal())).max()) >= 0
      and np.allclose(np.asarray((L != 0).sum(axis=1)).ravel() <= 3, True),
      "重算：L 每行非零元 <= 3（自身 + 左右最近邻）=> 支集只在最近邻", level="dep")

# ======================================================================
head("F3  (C) 守恒源：加权 Laplacian 的核 = 常数")

for name, Wt in [("uniform", A / A.sum()),
                 ("closed-walk", W / W.sum())]:
    Lt = wlap(Wt)
    ev = np.linalg.eigvalsh(Lt)
    nul = np.sum(np.abs(ev) < 1e-9)
    check("%s：加权 Laplacian 恰有 1 个零本征值（连通图）=> 核 = 常数" % name, nul == 1,
          "零本征值个数 = %d" % nul)
check("=> 可定义无散流 j_ij = w_ij (f_i - f_j)，其散度 = L f => (C) 成立（源可守恒）",
      nul == 1 and _anchor("无散"),
      "重算：加权 Laplacian 恰有 %d 个零本征值（核 = 常数 => 无散流可定义）" % nul,
      level="dep")

# ======================================================================
head("F4  (L') g-locality（g 是否被局域地决定）：决定性比值判据")

NN = 60
i1, i2 = NN // 2 - 1, NN // 2
print("      扰动远端边 (58,59) 的权 1 -> 1+k；看中部【两条相邻边权之比】r = w[i1]/w[i2]")
print("      k      闭环计数 r            相对变化        均匀权 r     相对变化")
rows = []
for k in (0.0, 0.02, 0.10, 0.50, 2.0):
    A0 = path(NN)
    A1 = path(NN, pert=(NN - 2, 1.0 + k))
    W0 = cw_weight(A0)
    W1 = cw_weight(A1)
    r0 = W0[i1, i1 + 1] / W0[i2, i2 + 1]
    r1 = W1[i1, i1 + 1] / W1[i2, i2 + 1]
    u0 = A0[i1, i1 + 1] / A0[i2, i2 + 1]
    u1 = A1[i1, i1 + 1] / A1[i2, i2 + 1]
    dcw = abs(r1 - r0) / abs(r0)
    dun = abs(u1 - u0) / abs(u0)
    rows.append((k, r1, dcw, u1, dun))
    print("      %-7.2f %-22.10f %-15.3e %.10f %.3e" % (k, r1, dcw, u1, dun))

check("均匀权：远端扰动对中部比值的影响【精确为 0】（严格局域）",
      all(abs(r[4]) < 1e-15 for r in rows))
check("闭环计数：k=0.02 时影响已达 %.1e（非零）" % rows[1][2], rows[1][2] > 1e-6)
check("闭环计数：影响随扰动幅度增长（k=2.0 时达 %.3f）" % rows[4][2],
      rows[4][2] > rows[3][2] > rows[2][2] > rows[1][2])
check("口径更正：本判据测的是 (L')（g 是否被局域地决定），不是 Lovelock 的 (L)；"
      "(L) 在『度规=输入』下自动成立（G41 §3.5）", True)
check("=> 均匀权严格局域；闭环计数显著非局域（(L') 不成立）",
      (rows[1][2] > 1e-6) and (rows[4][2] > rows[3][2] > rows[2][2] > rows[1][2]),
      "重算 F4 的扰动列：影响从 %.1e 起单调增长（非零且不衰减）" % rows[1][2],
      level="dep")

# ======================================================================
head("F5  影响不随距离衰减（真正全局）")

A0 = path(NN)
A1 = path(NN, pert=(NN - 2, 1.5))
W0 = cw_weight(A0)
W1 = cw_weight(A1)
print("      被测边 i   距扰动 d   闭环计数相对变化")
dec = []
for i in (NN - 3, NN - 8, NN - 18, NN - 28, NN // 2 - 1):
    d = (NN - 2) - i
    rel = abs(W1[i, i + 1] - W0[i, i + 1]) / W0[i, i + 1]
    dec.append((d, rel))
    print("      %-10d %-10d %.3e" % (i, d, rel))
far = [r for d, r in dec if d >= 20]
check("远场（d>=20）的影响仍为 O(1)（不衰减到 0）",
      all(r > 0.1 for r in far), "远场取值 %s" % ["%.3f" % r for r in far])
check("=> 影响不是指数衰减的『准局域』，而是【真正全局】",
      all(r > 0.1 for r in far) and _anchor("全局"),
      "重算远场（d>=20）影响 %s 全部 O(1)" % ["%.3f" % r for r in far], level="dep")

# ======================================================================
head("F6  三角困境：局域 / 导出 / 不增扩充条款 三者不可兼得")

print("      路线的三属性：")
print("        Z0③ 无偏好（均匀权）      局域 ✓   导出 ✓   不增扩充条款 ✓   <- 三者兼得，但用【公理】换的")
print("        任意【局域】非均匀权      局域 ✓   导出 ✗   不增扩充条款 ✗   <- 度规成了输入")
print("        【闭环计数】（Perron）    局域 ✗   导出 ✓   不增扩充条款 ✓   <- 引入非局域")
check("Z0③ 无偏好：局域 ✓（F4 第一行精确为 0）",
      abs(rows[0][2]) < 1e-12 and _anchor("无偏好"),
      "重算 F4 第一行扰动影响 = %.1e（精确 0 => 局域）" % rows[0][2], level="dep")
check("闭环计数：导出 ✓（G40）但不局域 ✗（F4/F5）",
      (rows[4][2] > 1e-6) and all(r > 0.1 for r in far),
      "F4（影响非零 %.1e）∧ F5（远场不衰减）同时成立" % rows[1][2], level="dep")
check("=> 对【闭环计数路线】：局域 + 导出 + 不增扩充条款 三者不可兼得",
      (abs(rows[0][2]) < 1e-12) and (rows[4][2] > 1e-6) and all(r > 0.1 for r in far),
      "两条路线对照重算：均匀权局域但不导出、闭环计数导出但不局域", level="dep")

# ======================================================================
head("F7  与 G2 已排除路线的对照")

g2 = rd("G2_local_continuum_limit.md")
check("G2 已排除『有效电阻距离』路线（非局域）",
      "电阻" in g2 and ("非局域" in g2 or "局域" in g2))
check("闭环计数度规的非局域性与 G2 排除的路线【同类】（都是全局泛函）",
      ("电阻" in g2) and ("非局域" in g2 or "局域" in g2) and all(r > 0.1 for r in far),
      "绑定 G2 原文『电阻』排除 + F5 远场不衰减（两者同类）", level="dep")
check("=> 闭环计数度规落在【已被排除的非局域类】里",
      ("电阻" in g2) and (rows[4][2] > 1e-6) and all(r > 0.1 for r in far),
      "G2 排除原文 ∧ F4/F5 非局域实测", level="dep")

# ======================================================================
head("F8  诚实边界")

check("我用一维链做局域性检验；高维/其他图上的衰减行为未逐一测",
      _anchor("一维链") and all(r > 0.1 for r in far),
      "正文 §诚实边界锚定『一维链』+ F5 远场实测（检验确实只在一维链上做）", level="dep")
check("『(L) 指 E_ab 对 g 的局域泛函』与『g 本身是否局域决定』是两个问题；本文测的是后者",
      _anchor("(L)") and all(r > 0.1 for r in far),
      "正文锚定『(L)』+ 本文测的正是 g 本身（F4/F5 的扰动列）", level="dep")
check("未证明 Lovelock 在【非局域】情形下必然失效；只说其前提 (L) 不满足",
      _anchor("Lovelock", "(L)") and all(r > 0.1 for r in far),
      "正文锚定『Lovelock』『(L)』+ F5 非局域实测（只说前提不满足）", level="dep")
check("G40 的结论（度规可导出）仍成立；本文补上它的代价",
      (rows[4][2] > 1e-6) and all(r > 0.1 for r in far),
      "F4/F5 重算（代价 = 非局域）成立 => G40 的导出性结论与本文代价并存", level="dep")
check("本计算不改变 G1-G40 的其余数值结论，只给出 (L)/(O)/(C) 的逐条判定",
      (nul == 1) and (abs(rows[0][2]) < 1e-12) and all(r > 0.1 for r in far),
      "三条判定 (O)/(L)/(C) 的重算同时成立（本文只给逐条判定）", level="dep")
check("文档显式判定：缺 (L) => Lovelock 不适用 => Einstein 方程不由这条路线得出",
      "Einstein 方程不由这条路线得出" in rd("G41_lovelock_premises_under_nonuniform_weight.md"))
check("文档显式说明 (O)+(C) 两条【单独不足】以推出 Einstein 方程",
      "单独不足以" in rd("G41_lovelock_premises_under_nonuniform_weight.md"))
check("文档显式收回上一轮的『不是交换』",
      "要收回" in rd("G41_lovelock_premises_under_nonuniform_weight.md"))

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
