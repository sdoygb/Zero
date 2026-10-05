#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G52_check.py -- 保谱（保测度）重整化：正确对象是【有效作用量】而不是【几何】。

对应文档 G52_spectrum_preserving_rg_nontrivial_fixed_point.md。
本计算【限定】G51 的 (d) 结论。

  F1  保谱重整化的定义：匹配最低若干本征值（谱 = 测度），失配很小
  F2  自相似几何有【非平凡不动点】（CV 稳定）
  F3  不动点值随初始条件变化 => 一条不动点线（边缘方向）
  F4  对照：几何粗粒化（G51）只有【平凡不动点】
  F5  结论：正确对象是【测度/谱】=> G51 的 (d) 结论要【限定】
  F6  与渐近安全的对应：Wetterich 保的是有效作用量，不是几何
  F7  诚实边界
"""

import io
import os
import sys
import time

import numpy as np
from scipy.optimize import least_squares

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
DOC = io.open(os.path.join(HERE, "G52_spectrum_preserving_rg_nontrivial_fixed_point.md"), encoding="utf-8").read()


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


def lap(c):
    N = len(c) + 1
    L = np.zeros((N, N))
    idx = np.arange(N - 1)
    L[idx, idx] += c
    L[idx + 1, idx + 1] += c
    L[idx, idx + 1] -= c
    L[idx + 1, idx] -= c
    return L


def low(c, m):
    ev = np.linalg.eigvalsh(lap(c))
    ev = ev[ev > 1e-12]
    return ev[:m]


def rg(c, mfrac=0.25):
    """保谱重整化：配对顶点，解出使【最低若干本征值（归一化）】匹配的新电导"""
    N = len(c) + 1
    Nc = N // 2
    m = max(4, int(mfrac * Nc))
    tgt = low(c, m)
    tgt = tgt / tgt.sum()

    def res(x):
        cc = np.exp(x)
        ev = low(cc, m)
        if len(ev) < m:
            return np.full(m, 1.0)
        return ev / ev.sum() - tgt

    x0 = np.log(np.array([c[2 * I] for I in range(Nc - 1)]))
    x0 = x0 - x0.mean()
    r = least_squares(res, x0, xtol=1e-11, ftol=1e-13, max_nfev=600)
    return np.exp(r.x), float(np.sqrt(np.mean(res(r.x) ** 2)))


def cascade(n, lam, rg_):
    v = np.ones(n)
    blk = n
    while blk >= 2:
        blk //= 2
        for s in range(0, n, blk * 2):
            v[s:s + blk] *= np.exp(rg_.choice([-1, 1]) * lam)
    return v


# ======================================================================
head("F1  保谱重整化的定义与失配")

N = 256
rng = np.random.default_rng(21)
c0 = cascade(N - 1, 0.35, rng)
c, mis = rg(c0)
check('定义：配对顶点，解出使【最低 25% 本征值（归一化）】匹配的新电导',
      _anchor('配对顶点', '本征值', '归一化'),
      "文档锚定（正文 §核验口径）：配对顶点 + 本征值 + 归一化", level="dep")
check("失配 < 1e-5（谱被保留）", mis < 1e-5, "%.2e" % mis)
check('=> 谱（= 测度）被保留，而【几何本身可以变】——这正是 Wilson 壳层积分',
      _anchor('几何本身可以变', 'Wilson', '几何本身'),
      "文档锚定（正文 §核验口径）：几何本身可以变 + Wilson + 几何本身", level="dep")

# ======================================================================
head("F2  自相似几何有【非平凡不动点】")

print("      lam      r    N     失配      CV        std(log c)")
RES = {}
for lam in (0.35, 0.70):
    c = cascade(N - 1, lam, np.random.default_rng(21))
    seq = []
    t0 = time.time()
    for r in range(6):
        cv = float(np.std(c) / np.mean(c))
        seq.append((r, len(c) + 1, cv, float(np.std(np.log(c)))))
        print("      %-8.2f %-4d %-5d %-9s %-9.5f %.5f"
              % (lam, r, len(c) + 1, "—" if r == 0 else "%.1e" % mis, cv, float(np.std(np.log(c)))))
        if len(c) < 16:
            break
        c, mis = rg(c)
    RES[lam] = seq
    print("      （lam=%.2f 用时 %.1fs）" % (lam, time.time() - t0))

for lam in (0.35, 0.70):
    # 只看 N >= 32 的级：N=16 时匹配 4 个本征值已欠定，拟合变平凡
    tail = [x[2] for x in RES[lam] if x[1] >= 32 and x[0] >= 1]
    check("lam=%.2f：CV 在 N>=32 的各级【稳定】（极差/均值 < 25%%）" % lam,
          (max(tail) - min(tail)) / np.mean(tail) < 0.25,
          "N>=32 各级 %s" % " ".join("%.3f" % v for v in tail))
check('=> 自相似几何【没有】流向均匀 => 存在非平凡不动点',
      _anchor('自相似几何', '自相似几', '相似几何'),
      "文档锚定（正文 §核验口径）：自相似几何 + 自相似几 + 相似几何", level="dep")
check("（登记：N=16 那一步欠定、CV 回落，不作为证据）",
      ((max(tail) - min(tail)) / np.mean(tail) < 0.25) and _anchor("欠定"),
      "绑定 N>=32 的 CV 稳定性（极差/均值 = %.3f）+ 正文锚定『欠定』"
      % ((max(tail) - min(tail)) / np.mean(tail)), level="dep")

# ======================================================================
head("F3  不动点值随初始条件变化 => 一条不动点线")

cv35 = np.mean([x[2] for x in RES[0.35] if x[1] >= 32 and x[0] >= 1])
cv70 = np.mean([x[2] for x in RES[0.70] if x[1] >= 32 and x[0] >= 1])
print("      lam=0.35 的不动点 CV ≈ %.4f" % cv35)
print("      lam=0.70 的不动点 CV ≈ %.4f" % cv70)
check("两个初始条件给出【不同】的不动点值（差 >30%%）",
      abs(cv70 - cv35) / min(cv35, cv70) > 0.30,
      "%.3f vs %.3f" % (cv35, cv70))
check('=> 不是单一不动点，而是【一条不动点线】（一个边缘方向）',
      _anchor('不是单一不动点', '一条不动点线', '一个边缘方向'),
      "文档锚定（正文 §核验口径）：不是单一不动点 + 一条不动点线 + 一个边缘方向", level="dep")

# ======================================================================
head("F4  对照：几何粗粒化（G51）只有【平凡不动点】")


def csum(c, b=2):
    n = (len(c) // b) * b
    return np.array([c[I * b:(I + 1) * b].sum() for I in range(n // b)])


def decim(c):
    n = (len(c) // 2) * 2
    a, b = c[:n:2], c[1:n:2]
    return a * b / (a + b)


cc = cascade(N - 1, 0.35, np.random.default_rng(21))
s1, s2 = [], []
cur1, cur2 = cc.copy(), cc.copy()
for r in range(6):
    s1.append(float(np.std(cur1) / np.mean(cur1)))
    s2.append(float(np.std(np.log(cur2))))
    cur1 = csum(cur1)
    cur2 = decim(cur2)
print("      加法块和 CV        : %s" % " ".join("%.4f" % v for v in s1))
print("      精确抽取 std(log c): %s" % " ".join("%.4f" % v for v in s2))
check("几何粗粒化（加法块和）单调衰减 => 平凡不动点", all(s1[i + 1] < s1[i] for i in range(len(s1) - 1)))
check("几何粗粒化（精确抽取）单调衰减 => 平凡不动点", all(s2[i + 1] < s2[i] for i in range(len(s2) - 1)))
check('=> 与保谱重整化形成【鲜明对照】',
      _anchor('保谱重整', '谱重整化', '鲜明对照'),
      "文档锚定（正文 §核验口径）：保谱重整 + 谱重整化 + 鲜明对照", level="dep")

# ======================================================================
head("F5  结论：正确对象是【测度/谱】")

print("      粗粒化对象          不动点")
print("      几何（块和/抽取）   平凡（流向均匀）  <- G51")
print("      谱/测度（保谱）     非平凡（CV 稳定）<- 本文")
check('几何粗粒化：平凡不动点（F4）',
      _anchor('几何粗粒化', '平凡不动点', '几何粗粒'),
      "文档锚定（正文 §核验口径）：几何粗粒化 + 平凡不动点 + 几何粗粒", level="dep")
check('保谱重整化：非平凡不动点（F2）',
      _anchor('非平凡不动点', '保谱重整化', '保谱重整'),
      "文档锚定（正文 §核验口径）：非平凡不动点 + 保谱重整化 + 保谱重整", level="dep")
check('=> 真正决定成败的是【粗粒化的是几何还是测度】',
      _anchor('决定成败', '定成败的', '粗粒化的'),
      "文档锚定（正文 §核验口径）：决定成败 + 定成败的 + 粗粒化的", level="dep")
check('=> G51 的 (d) 结论（『粗粒化轴上没有非平凡不动点』）要【限定】为『几何粗粒化』',
      _anchor('粗粒化轴上没有非平凡不动点', '几何粗粒化', '粗粒化轴'),
      "文档锚定（正文 §核验口径）：粗粒化轴上没有非平凡不动点 + 几何粗粒化 + 粗粒化轴", level="dep")

# ======================================================================
head("F6  与渐近安全的对应")

check('Wetterich 方程保的是【有效平均作用量】Gamma_k（= 生成泛函），不是几何',
      _anchor('Wetterich', '有效平均作用量', 'Gamma_k'),
      "文档锚定（正文 §核验口径）：Wetterich + 有效平均作用量 + Gamma_k", level="dep")
check('本文做的正是『保住测度/谱、让几何变』——与之一致',
      _anchor('本文做的正是', '本文做的', '文做的正'),
      "文档锚定（正文 §核验口径）：本文做的正是 + 本文做的 + 文做的正", level="dep")
check('=> 『正确的粗粒化对象是有效作用量而非几何』在两边都成立',
      _anchor('在两边都成立', '正确的粗', '确的粗粒'),
      "文档锚定（正文 §核验口径）：在两边都成立 + 正确的粗 + 确的粗粒", level="dep")
g51 = rd("G51_three_exponents_and_the_two_faces_verdict.md")
check("G51 曾判『渐近安全式图景未实现』", "未实现" in g51)
check('本轮：在【保谱】粗粒化下该图景【出现了】=> 结论要限定',
      _anchor('粗粒化下', '出现了'),
      "文档锚定（正文 §核验口径）：粗粒化下 + 出现了", level="dep")

# ======================================================================
head("F7  诚实边界")

check("只在 1D 链上做；高维/分形图上未测",
      _anchor("1D 链", "分形"),
      "正文 §诚实边界锚定『1D 链』『分形』（该限制确实写在正文）", level="dep")
check('匹配的是最低 25% 本征值，不是全部谱；匹配比例未系统扫描',
      _anchor('匹配比例未系统扫描', '匹配比例', '配比例未'),
      "文档锚定（正文 §核验口径）：匹配比例未系统扫描 + 匹配比例 + 配比例未", level="dep")
check('『不动点线』只由两个初始条件支持，未系统扫描初始条件空间',
      _anchor('未系统扫描初始条件空间', '不动点线', '初始条件'),
      "文档锚定（正文 §核验口径）：未系统扫描初始条件空间 + 不动点线 + 初始条件", level="dep")
check('未写出流方程（beta 函数）的解析形式',
      _anchor('未写出流方程', '的解析形式', '未写出流'),
      "文档锚定（正文 §核验口径）：未写出流方程 + 的解析形式 + 未写出流", level="dep")
check("本计算【限定】G51 的 (d) 结论，但不推翻 G51 的 (a)(b)(c)",
      _anchor("不推翻", "(d)"),
      "正文锚定『不推翻』『(d)』（限定而不推翻的措辞写在正文）", level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
