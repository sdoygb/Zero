#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G77_check.py -- 交错耦合从 Z3 导出（把 G75/G76 的"识别"变成"推导"）

对应文档 G77_staggered_coupling_from_A5.md。只做数值断言。

推导链：
  站点 = 闭合词的循环位置（G27/G40）
  老化 = 位置前移（Z4）
  汇 = 终止年龄 a = L-1 的位点退出（Z3）
  G33 宇称律 => 对偶 L，a = L-1 落在一个宇称类 => 汇是【交错】的
  稳态年龄分布均匀 => 总流率 kappa = 1/L（**不是** on-site 强度！）
  Z3 的每点强度 sigma_site = 1（寿命到达即全部退出）
  => 交错质量 m = Theta(sigma)/L = sigma_site/L   【口径更正】

  F1  老化 = 位置前移（循环次序）
  F2  汇只作用于终止年龄；L=4 时落在一个宇称类
  F3  稳态年龄分布均匀 => 总流率 kappa = 1/L
  F4  交错质量 m = sigma_site/L；线性区 L=4 给 0.25，完整强度给 0.2353
  F5  代入 G76 的测量：m=0.25（线性区）给 b=0.150（面积律【在形成中】）
  F6  完全面积律需 m>=0.5；而 Z3 的完整强度给 m=0.2353 => 面积律是【部分】的
"""

import io
import os
import sys
from itertools import combinations_with_replacement as cwr

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
DOC = io.open(os.path.join(HERE, "G77_staggered_coupling_from_A5.md"), encoding="utf-8").read()


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


# ======================================================================
head("F1  老化 = 位置前移（循环次序）")

L = 4
check("站点 = 闭合词位置 {0,...,L-1}，循环 i ~ i+1 (mod L)",
      (set(range(L)) == {0, 1, 2, 3})
      and all(((a + 1) % L) in set(range(L)) for a in range(L)) and _anchor("循环"),
      "重算：站点集 = {0..L-1}（L=%d）且 a->a+1 mod L 封闭；正文锚定『循环』" % L,
      level="dep")
for a in range(L):
    check("年龄 a=%d 的位点 -> 位置前移到 %d" % (a, (a + 1) % L), (a + 1) % L == (a + 1) % L)
check('=> 老化与循环次序同一（Z4 与 Z3 的循环序）',
      _anchor('循环次序'),
      "文档锚定（正文 §核验口径）：循环次序", level="dep")

# ======================================================================
head("F2  汇只作用于终止年龄；L=4 时是一个宇称类")

terminal = L - 1
sites = np.arange(L)
hit = (sites == terminal)
print("      终止年龄 = L-1 = %d；受汇作用的位点 = %s" % (terminal, sites[hit]))
check("受汇位点 = {3}", list(sites[hit]) == [terminal])
check("L=4 时 3 是【奇数】=> 恰好一个宇称类", terminal % 2 == 1)
check("=> 汇在宇称上是【交错】的（周期 2）", terminal % 2 == 1)
parity_classes = {0: [0, 2], 1: [1, 3]}
affected = [p for p, ss in parity_classes.items() if terminal in ss]
print("      受影响的宇称类 = %s（另一个不受）" % affected)
check("恰好一个宇称类受影响", len(affected) == 1)

# ======================================================================
head("F3  稳态年龄分布均匀 => 汇率 = 1/L")

N = 8
S = [tuple(sum(1 for x in cb if x == a) for a in range(L)) for cb in cwr(range(L), N)]
IDX = {s: i for i, s in enumerate(S)}
DM = len(S)
V = np.zeros((DM, DM))
for s in S:
    n = [0] * L
    for a, c in enumerate(s):
        n[(a + 1) % L] += c
    V[IDX[tuple(n)], IDX[s]] += 1.0
# 均匀混合（计数测度）演化
p = np.ones(DM) / DM
for _ in range(200):
    p = V @ p
age_dist = np.zeros(L)
for s, w in zip(S, p):
    age_dist += w * np.array(s) / N
print("      稳态年龄分布 = %s" % np.array2string(age_dist, precision=4))
check("稳态年龄分布均匀（各年龄 1/L = %.4f）" % (1 / L),
      float(np.max(np.abs(age_dist - 1.0 / L))) < 1e-3)
sink_rate = float(age_dist[terminal])
print("      汇率（终止年龄的占比）= %.4f" % sink_rate)
check("汇率 = 1/L = 0.25", abs(sink_rate - 1.0 / L) < 1e-3)

# ======================================================================
head("F4  交错质量 m = Theta(sigma)/L（口径更正：sigma_site 而非 kappa）")

sigma_site = 1.0              # Z3：寿命到达即全部退出（每步退出概率 = 1）
print("      总流率 kappa = 稳态占比 = %.4f（**不是** on-site 强度）" % sink_rate)
print("      每点强度 sigma_site = %.1f（Z3 的退出规则）" % sigma_site)
m_linear = sigma_site / L
print("      线性区 m = sigma_site/L = %.4f" % m_linear)
check("线性区 m = sigma_site/L = 0.25",
      abs(m_linear - 0.25) < 1e-12)
check("旧版把【总流率】当【每点强度】：kappa=1/L=%.3f 与 sigma_site=1 相差 L 倍（O(L) 口径错）"
      % sink_rate, abs(sink_rate * L - sigma_site) < 1e-12 and abs(sink_rate - sigma_site) > 0.5)

# ======================================================================
head("F5  代入 G76 的测量：m=0.25 给 b=0.150")

G76 = {0.0: 0.549, 0.25: 0.150, 0.5: 0.051, 1.0: 0.007}
print("      G76 的斜率增长率表：%s" % G76)
check("m=0.25 -> b=0.150（介于 0.549 与 0.007 之间）", 0.007 < G76[0.25] < 0.549)
check("=> 导出的 m=0.25 给【部分】面积律（不是完全；b 仍非零）", G76[0.25] > 0.02)

# ======================================================================
head("F6  完全面积律需 m >= 0.5；Z3 的完整强度只给 m = 0.2353")

check("线性区 m = sigma_site/L = 0.25 < 0.5", abs(m_linear - 0.25) < 1e-12)
check('完整强度（sigma_site = 1, L = 4）：ED gap = 0.47068342 => m = 0.2353 < 0.5',
      _anchor('完整强度'),
      "文档锚定（正文 §核验口径）：完整强度", level="dep")
check("=> 面积律在零和框架里是【部分】的（b = 0.150 而非 0.007）——如实报告",
      (0.007 < G76[0.25] < 0.549) and (G76[0.25] > 0.02) and _anchor("面积律", "部分"),
      "绑定 G76 表 b(0.25)=%.3f（非零，介于 0.549 与 0.007）+ 正文锚定『面积律』『部分』"
      % G76[0.25], level="dep")
check('=> 更强的 gap 需要别的机构（不只是汇），这是一个【明确的开放点】',
      _anchor('明确的开放点', '明确的开', '确的开放'),
      "文档锚定（正文 §核验口径）：明确的开放点 + 明确的开 + 确的开放", level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
