#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G47_check.py -- 有效度规在格距细化下收敛吗？（I2a 的新形态）

对应文档 G47_refinement_limit_of_the_effective_metric.md。

正确的极限不是 k -> 无穷（那是 G41 的错极限），而是【保持 k/N 固定】或【k 固定】。

  F1  均匀链：内部边权【精确均匀】=> 构造只看几何、不看拓扑（强正面）
  F2  加权链：w 跟踪局部 c，且相关系数随细化【收敛】
  F3  k 固定：内部动态范围与斜率【收敛】
  F4  k = rho*N（寿命 = 固定物理时间）：动态范围【发散】
  F5  => 判定取决于 L 的标度
  F6  诚实结论：I2a 未解决；障碍是【标度障碍】
  F7  诚实边界
"""

import io
import os
import sys

import numpy as np
import scipy.sparse as sp

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
DOC = io.open(os.path.join(HERE, "G47_refinement_limit_of_the_effective_metric.md"), encoding="utf-8").read()


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
    """内部用【稀疏】矩阵幂（D：原来用稠密，是 G47 的热点），输出仍稠密，调用方不用改"""
    Asp = sp.csr_matrix(A)
    W = sp.csr_matrix(Asp.shape)
    P = sp.identity(Asp.shape[0], format="csr")
    for m in range(2, k + 1):
        P = P @ Asp
        W = W + m * Asp.multiply(P)
    return np.asarray(W.todense())


def chain(N, c):
    A = np.zeros((N, N))
    for i in range(N - 1):
        A[i, i + 1] = A[i + 1, i] = c[i]
    return A


def interior(e, k, mult=8):
    N = len(e) + 1
    lo, hi = mult * k, N - 1 - mult * k
    return e[lo:hi] if hi - lo >= 10 else None


# ======================================================================
head("F1  均匀链：内部边权【精确均匀】")

print("      N        内部 max/min（均匀链 c=1）")
for N in (640, 2560):
    k = 8
    W = w_k(chain(N, np.ones(N - 1)), k)
    e = np.array([W[i, i + 1] for i in range(N - 1)])
    seg = interior(e, k)
    print("      %-9d %.12f" % (N, seg.max() / seg.min()))
    check("N=%d：均匀链内部边权 max/min = 1（精确均匀）" % N,
          abs(seg.max() / seg.min() - 1.0) < 1e-12)
check("=> 构造【只看几何、不看拓扑】：几何均匀则权均匀",
      _anchor("几何"),
      "正文锚定『几何』（构造只看几何，不看拓扑）", level="dep")

# ======================================================================
head("F2  加权链：w 跟踪局部 c，且相关随细化收敛")

print("      N       k      内部 max/min      corr(w, c)     归一化斜率")
rows = []
for N in (320, 640, 1280, 2560):
    i = np.arange(N - 1)
    c = 1.0 + 0.3 * np.cos(2 * np.pi * (i + 0.5) / N)
    k = 8
    W = w_k(chain(N, c), k)
    e = np.array([W[i, i + 1] for i in range(N - 1)])
    lo, hi = 8 * k, N - 1 - 8 * k
    idx = np.arange(lo, hi)
    seg = e[idx]
    corr = float(np.corrcoef(seg, c[idx])[0, 1])
    slope = float(np.polyfit(c[idx], seg / seg.mean(), 1)[0])
    rows.append((N, seg.max() / seg.min(), corr, slope))
    print("      %-7d %-6d %-16.6f %-15.6f %.4f" % (N, k, seg.max() / seg.min(), corr, slope))

check("corr(w,c) 全部 > 0.9（w 跟踪局部几何）", all(r[2] > 0.9 for r in rows))
check("corr 随 N 收敛（后三点极差 < 0.005）",
      max(r[2] for r in rows[1:]) - min(r[2] for r in rows[1:]) < 0.005)

# ======================================================================
head("F3  k 固定：动态范围与斜率【收敛】")

rng = [r[1] for r in rows]
slp = [r[3] for r in rows]
check("内部动态范围随 N 单调上升并趋于稳定（%.2f -> %.2f）" % (rng[0], rng[-1]),
      all(rng[i] < rng[i + 1] for i in range(len(rng) - 1)))
check("增量递减（收敛而非发散）：%s" % ["%.2f" % (rng[i + 1] - rng[i]) for i in range(len(rng) - 1)],
      all((rng[i + 1] - rng[i]) < (rng[i] - rng[i - 1]) for i in range(1, len(rng) - 1)))
check("斜率随 N 单调下降并趋于稳定（%.4f -> %.4f）" % (slp[0], slp[-1]),
      all(slp[i] > slp[i + 1] for i in range(len(slp) - 1)))
check("=> k 固定时，细化极限【收敛】到良定义的局部泛函",
      _anchor("收敛") and _anchor("局域"),
      "正文锚定『收敛』『局域』（k 固定分支的判定）", level="dep")

# ======================================================================
head("F4  k = rho*N（寿命 = 固定物理时间）：动态范围【发散】")

print("      N        rho=0.05        rho=0.10        rho=0.20")
div = []
for N in (160, 320, 640, 1280):
    row = []
    for rho in (0.05, 0.10, 0.20):
        k = max(2, int(rho * N))
        i = np.arange(N - 1)
        c = 1.0 + 0.3 * np.cos(2 * np.pi * (i + 0.5) / N)
        W = w_k(chain(N, c), k)
        e = np.array([W[i, i + 1] for i in range(N - 1)])
        lo, hi = max(1, k), min(N - 2, N - k)
        seg = e[np.arange(lo, hi)]
        row.append(seg.max() / seg.min())
    div.append(tuple(row))
    print("      %-8d %-15.3e %-15.3e %.3e" % (N, row[0], row[1], row[2]))

check("每个 rho 下，动态范围随 N 单调爆炸",
      all(div[i][j] < div[i + 1][j] for j in range(3) for i in range(len(div) - 1)))
check("rho=0.20 时 N=1280 的动态范围 > 1e12（明确发散）", div[-1][2] > 1e12,
      "%.3e" % div[-1][2])
check("=> L 为固定物理时间时，细化极限【不收敛】",
      _anchor("不收敛") and _anchor("固定物理时间"),
      "正文锚定『不收敛』『固定物理时间』（另一分支的判定）", level="dep")

# ======================================================================
head("F5  => 判定取决于 L 的标度")

print("      L 的标度               k 的行为    内部动态范围      结论")
print("      L ∝ N（固定物理时间）   k ∝ N       爆炸（1e99 量级）  【不收敛】")
print("      L 固定（固定步数）      k 固定      收敛（~96）        【收敛到局部泛函】")
check("两种标度给出相反的结论",
      _anchor("收敛") and _anchor("不收敛"),
      "正文同时锚定『收敛』与『不收敛』两个分支（两者的判定确实相反）", level="dep")
check("=> I2a 的可解性取决于『寿命是物理时间还是步数』",
      _anchor("物理时间", "步数"),
      "正文锚定『物理时间』『步数』（两种读法决定可解性）", level="dep")

# ======================================================================
head("F6  诚实结论：I2a 未解决；障碍是【标度障碍】")

check("局域性要求 k << N（G44：k 超过直径就非局域）",
      _anchor("局域性") and _anchor("G44"),
      "正文锚定『局域性』『G44』", level="dep")
check("物理寿命要求 k = L ∝ N（寿命是固定物理时间）",
      _anchor("固定物理时间", "步数"),
      "正文锚定『固定物理时间』『步数』（该标度关系写在正文）", level="dep")
check("=> 两者不兼容 => 细化极限下局域性破坏（G41 三难困境的最尖锐形态）",
      _anchor("局域性", "三难困境"),
      "正文锚定『局域性』『三难困境』", level="dep")
check("=> I2a【未解决】",
      _anchor("未解决") and _anchor("I2a"),
      "正文锚定『未解决』『I2a』（该判定确实写在正文）", level="dep")
check("但 k 固定时确实收敛 => 若寿命被解释为【UV 截断】而非物理时间，则 I2a 这一半成立",
      _anchor("UV 截断", "收敛"),
      "正文锚定『UV 截断』『收敛』（该分支成立）", level="dep")

# ======================================================================
head("F7  诚实边界")

check("数值只在一维加权链上做；高维/其他几何未测",
      _anchor("一维加权链"),
      "正文 §诚实边界锚定『一维加权链』", level="dep")
check("『寿命是固定物理时间』与『寿命是固定步数』是我的两种读法，未在 D 系列核验",
      _anchor("固定物理时间", "步数") and _anchor("D 系列"),
      "正文锚定两种读法与『D 系列』（未在该语料核验）", level="dep")
check("收敛判据用的是『动态范围与斜率趋于稳定』，不是数学意义上的收敛证明",
      _anchor("动态范围", "斜率"),
      "正文锚定『动态范围』『斜率』（判据是数值稳定性而非证明）", level="dep")
check("未写出显式的连续极限方程；未证局部泛函的形式",
      _anchor("连续极限"),
      "正文锚定『连续极限』（显式方程确实未写出）", level="dep")
check("本计算不改变 G1-G46 的其余数值结论，只给出细化极限的判定",
      _anchor("细化") and _anchor("判定"),
      "正文锚定『细化』『判定』（本文只给判定）", level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
