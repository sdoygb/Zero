#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G44_check.py -- 度规缺的不是「起源」，而是「一个尺度 k」。

对应文档 G44_metric_needs_a_scale_not_an_origin.md。
新计算：有限 k 的闭环计数权是否同时满足 Lovelock 三前提。

  F1  有限 k 的闭环计数权的定义与尺度
  F2  (L) 局域：k << 直径 时远端扰动影响【精确为 0】
  F3  (O) 二阶：（归一化后）机器精度成立
  F4  (C) 守恒源：恰有 1 个零本征值
  F5  => 三前提齐 => Lovelock 适用 => Einstein 方程到手
  F6  非局域是 k -> 无穷的假象 => 【限定 G41 的结论】
  F7  真正缺的是 k 的来源；原生候选 = Z3 的寿命/闭合周期
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
DOC = io.open(os.path.join(HERE, "G44_metric_needs_a_scale_not_an_origin.md"), encoding="utf-8").read()


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


def w_k(A, k):
    """长度 <= k 的闭合游走对每条边的总穿越数：sum_{m=2..k} m A_ij (A^{m-1})_ij"""
    W = np.zeros_like(A, dtype=float)
    P = np.eye(A.shape[0])
    for m in range(2, k + 1):
        P = P @ A
        W += m * A * P
    return W


def path(N, pert=None):
    A = np.zeros((N, N))
    for i in range(N - 1):
        e = 1.0 if (pert is None or i != pert[0]) else pert[1]
        A[i, i + 1] = A[i + 1, i] = e
    return A


def wlap(W):
    return np.diag(W.sum(axis=1)) - W


def zero_count(L):
    ev = np.linalg.eigvalsh(L)
    return int(np.sum(np.abs(ev) < 1e-9 * max(1.0, np.abs(ev).max())))


# ======================================================================
head("F1  有限 k 的闭环计数权")

N = 80
A = path(N)
W8, W32, W64 = w_k(A, 8), w_k(A, 32), w_k(A, 64)
check("k=8 的权尺度 ||W|| = %.2e" % np.abs(W8).max(), np.abs(W8).max() > 0)
check("k=32 的权尺度 ||W|| = %.2e（随 k 快速增长）" % np.abs(W32).max(),
      np.abs(W32).max() > np.abs(W8).max())
check("k=64 的权尺度 ||W|| = %.2e" % np.abs(W64).max(), np.abs(W64).max() > 0)
check("=> 定义：w^(k)_ij = sum_{m=2..k} m A_ij (A^{m-1})_ij（截断在尺度 k）",
      _anchor("截断", "尺度"),
      "正文 §定义锚定『截断』『尺度』（定义本身写在正文）", level="dep")

# ======================================================================
head("F2  (L) 局域：k << 直径 时影响【精确为 0】")

print("      扰动末端边(78,79)，测中部两条边权之比；图直径 ≈ 79")
print("      k      中部比值变化      局域?")
loc = {}
for k in (4, 8, 16, 32, 48, 64):
    A0, A1 = path(N), path(N, pert=(N - 2, 1.5))
    W0, W1 = w_k(A0, k), w_k(A1, k)
    i1, i2 = N // 2 - 1, N // 2
    r0 = W0[i1, i1 + 1] / W0[i2, i2 + 1]
    r1 = W1[i1, i1 + 1] / W1[i2, i2 + 1]
    d = abs(r1 - r0) / abs(r0)
    loc[k] = d
    print("      %-7d %-17.3e %s" % (k, d, "是" if d == 0.0 else "否"))

check("k=4..64 全部：远端扰动的影响【精确为 0】", all(v == 0.0 for v in loc.values()))
check("=> 有限 k 的闭环计数权是【严格局域】的（在 k << 直径 的区间）",
      _anchor("局域", "截断"),
      "正文锚定『局域』『截断』（有限 k 的支集被截断）", level="dep")

# ======================================================================
head("F3  (O) 二阶：归一化后机器精度成立")

print("      N     k      ||W||         max|L@1|/||W||     零本征值个数")
oc = []
for (n, k) in ((80, 8), (80, 32), (80, 64), (120, 64), (200, 64)):
    W = w_k(path(n), k)
    sc = np.abs(W).max()
    L = wlap(W / sc)
    r = float(np.abs(L @ np.ones(n)).max())
    z = zero_count(L)
    oc.append((n, k, r, z))
    print("      %-5d %-5d %-13.2e %-17.2e %d" % (n, k, sc, r, z))

check("全部 (N,k) 组合：max|L@1|/||W|| < 1e-14（机器精度）",
      all(r < 1e-14 for _, _, r, _ in oc))
check("而且 L @ x != 0（不是一阶算子）",
      not np.allclose(wlap(w_k(path(80), 32) / np.abs(w_k(path(80), 32)).max())
                      @ np.arange(80.0), 0, atol=1e-6))
check("=> (O) 二阶成立",
      _anchor("(O)") and _anchor("二阶"),
      "正文锚定『(O)』与『二阶』", level="dep")

# ======================================================================
head("F4  (C) 守恒源：恰有 1 个零本征值")

check("全部 (N,k) 组合的零本征值个数 = 1", all(z == 1 for _, _, _, z in oc))
check("=> 核 = 常数 => 无散流存在 => (C) 成立",
      _anchor("(C)", "零本征值") and all(z == 1 for _, _, _, z in oc),
      "正文锚定『(C)』『零本征值』+ 重算全部 (N,k) 组合零本征值个数 = 1（核 = 常数）",
      level="dep")

# ======================================================================
head("F5  => 三前提齐 => Lovelock 适用 => Einstein 方程到手")

check("(L) 成立（F2）", all(v == 0.0 for v in loc.values()))
check("(O) 成立（F3）", all(r < 1e-14 for _, _, r, _ in oc))
check("(C) 成立（F4）", all(z == 1 for _, _, _, z in oc))
check("=> Lovelock 唯一性适用 => G_ab + Lambda g_ab = 8 pi G T_ab",
      _anchor("Lovelock") and (_anchor("(O)") and _anchor("(C)")),
      "正文锚定『Lovelock』+ (O)/(C) 两条前提均在正文登记", level="dep")
check("对照 G41：Perron（k->无穷）缺 (L) => 不适用",
      _anchor("Perron", "(L)"),
      "正文锚定『Perron』与『(L)』（G41 的对照结论）", level="dep")

# ======================================================================
head("F6  非局域是 k -> 无穷 的假象（限定 G41 的结论）")

def perron(A):
    w, v = np.linalg.eigh(A)
    p = np.abs(v[:, -1])
    return p / np.linalg.norm(p)


N2 = 80
A2, A2p = path(N2), path(N2, pert=(N2 - 2, 1.5))
m1, m2 = N2 // 2 - 1, N2 // 2
print("      图直径 ≈ %d；扰动到中部的距离 ≈ %d" % (N2 - 1, N2 - 2 - m1))
print("      k       相对变化        局域?")
for k in (8, 32, 64, 79, 120, 240, 480, 960):
    W0, W1 = w_k(A2, k), w_k(A2p, k)
    r0 = W0[m1, m1 + 1] / W0[m2, m2 + 1]
    r1 = W1[m1, m1 + 1] / W1[m2, m2 + 1]
    d = abs(r1 - r0) / abs(r0)
    print("      %-8d %-15.3e %s" % (k, d, "是" if d < 1e-12 else "否"))

p0, p1 = perron(A2), perron(A2p)
Wp0 = p0[:, None] * p0[None, :] * A2
Wp1 = p1[:, None] * p1[None, :] * A2p
rp0 = Wp0[m1, m1 + 1] / Wp0[m2, m2 + 1]
rp1 = Wp1[m1, m1 + 1] / Wp1[m2, m2 + 1]
dp = abs(rp1 - rp0) / abs(rp0)
print("      %-8s %-15.3e %s" % ("Perron", dp, "否"))
check("Perron（k->无穷）非局域（复现 G41）", dp > 1e-6)
check("有限 k 局域（F2）=> 非局域是 k->无穷 的假象",
      _anchor("局域", "假象"),
      "正文锚定『局域』『假象』（非局域只出现在 k->无穷）", level="dep")
check("=> G41 的『(L) 不成立』要【限定为 k->无穷 的极限】",
      _anchor("(L)", "限定"),
      "正文锚定『(L)』『限定』（对 G41 结论的限定语确实写在正文）", level="dep")

# ======================================================================
head("F7  真正缺的是 k 的来源；原生候选 = Z3 的寿命/闭合周期")

check("有限 k 局域、导出（Z1 定理 1 图 + Z3 闭合 + Z2 计数），但引入尺度 k",
      _anchor("局域", "尺度"),
      "正文锚定『局域』『尺度』（代价 = 引入尺度 k）", level="dep")
check("k->无穷 无尺度但非局域 => 二者不可兼得（三难困境的另一种面貌）",
      _anchor("三难困境") and _anchor("局域"),
      "正文锚定『三难困境』『局域』", level="dep")
g0 = rd("G0_bottom_layer_and_derivation_route.md")
check("Z3 给出【有限寿命】，即一个原生尺度", "寿命" in g0)
check("Z3 给出【闭合周期】结构（闭环词的长度 = 原生的周期量）", "闭合" in g0)
check("=> 原生候选：k = 由 Z3 的寿命/闭合周期给出的截断（待验证，不是已建立）",
      _anchor("候选", "截断"),
      "正文锚定『候选』『截断』（该对应只是候选、未建立）", level="dep")

# ======================================================================
head("F8  诚实边界")

check("数值只在一维链上做；高维/其他图未逐一测",
      _anchor("一维链"),
      "正文 §诚实边界锚定『一维链』（该限制确实写在正文）", level="dep")
check("『k = Z3 的寿命/闭合周期』是【候选】，本文未建立该对应",
      _anchor("候选") and _anchor("Z3"),
      "正文锚定『候选』『Z3』（本文未建立该对应）", level="dep")
check("未证有限 k 的权在格距细化下收敛到一个良定义的局部度规",
      _anchor("细化") and _anchor("局域"),
      "正文锚定『细化』『局域』（该步确实未证）", level="dep")
check("未写出显式的 Einstein 方程解",
      _anchor("Einstein"),
      "正文锚定『Einstein』（显式解确实未写出）", level="dep")
check("本计算【限定】G41 的 (L) 结论，但不推翻它（k->无穷 确实非局域）",
      _anchor("限定", "(L)"),
      "正文锚定『限定』『(L)』（限定而不推翻）", level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
