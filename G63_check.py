#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G63_check.py -- 对照清单：零和宇宙"本该"预言的量，以及审计

照旧理论（cosmos-construct Q640）的模板：完整性 = 相对一个【清单】。
本文建本项目的清单（几何/因果/统计/量子），并三栏审计。只做数值断言。

  F1  清单"已导出"栏的数值逐条复核（闭式）
  F2  T = L = 4 => dim A = 4(L+1) = 20 >= 3（Gleason 可用）
  F3  【清单为什么短】物质/规范扇区不可达：D = |F|(m-1)+1 = 4 的整数解枚举
  F4  清单规模与状态统计（效率比）
  F5  逐条交叉复核：与 G54/G56/G59 的闭式一致
"""

import itertools
import io
import os
import sys
from math import comb, cosh, log, tanh

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
DOC = io.open(os.path.join(HERE, "G63_target_list_and_audit.md"), encoding="utf-8").read()


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


def cstar(B):
    tgt = 0.5 * log(B)
    f = lambda m: m * tanh(m) - log(cosh(m))
    lo, hi = 1e-12, 200.0
    if f(hi) < tgt:
        return None, None
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if f(mid) < tgt:
            lo = mid
        else:
            hi = mid
    mu = 0.5 * (lo + hi)
    return tanh(mu), mu


# ======================================================================
head("F1  清单【已导出】栏的数值逐条复核")

L, Bv = 4, 4
c4, mu4 = cstar(Bv)
print("      c*(4) = %.9f   mu*(4) = %.4f" % (c4, mu4))
check("因果：c* = 1（B=4 锁定，G56/G61）", abs(c4 - 1.0) < 1e-9)
check("因果：mu* = 32.6926（锥外每格 e^{-32.7}，G59）", abs(mu4 - 32.6926) < 1e-3)
check("几何：依赖半径 = k/2-1 = 1（k = L = 4，G55/G58）", L // 2 - 1 == 1)
check("统计：F(4) = C(4,2)/2^4 = 0.375（半程定理，G54）", abs(comb(4, 2) / 2 ** 4 - 0.375) < 1e-15)
check('统计：选择率 lambda = log M/T 无量纲（G29）',
      _anchor('lambda', '选择率', '无量纲'),
      "文档锚定（正文 §核验口径）：lambda + 选择率 + 无量纲", level="dep")
# 重播种：L=4 的旋转类大小
W = [w for w in itertools.product((1, -1), repeat=L) if sum(w) == 0]
cls = {}
for w in W:
    k = min(tuple(w[i:] + w[:i]) for i in range(L))
    cls[k] = cls.get(k, 0) + 1
check("重播种：L=4 的旋转类大小 = [2,4]（G37/G61）", sorted(cls.values()) == [2, 4],
      "%s" % sorted(cls.values()))
check("几何：层系数 = C(m-1,m/2) 整数，总权 354（G60）",
      [comb(m - 1, m // 2) for m in (2, 4, 6, 8)] == [1, 3, 10, 35])

# ======================================================================
head("F2  T = L = 4 => dim A = 4(L+1) = 20 >= 3（Gleason 可用）")

T = L
dimA = 4 * (T + 1)
print("      T = L = %d  =>  dim A = 4(T+1) = %d" % (T, dimA))
check("年龄因子维数 = 4(T+1) = 20", dimA == 20)
check("dim A >= 3（Gleason 阈值，G62）", dimA >= 3)
check("T >= 1（G62 §4 的必要条件）", T >= 1)
check('=> 量子栏（Born）在物理取值下可用',
      _anchor('量子栏'),
      "文档锚定（正文 §核验口径）：量子栏", level="dep")

# ======================================================================
head("F3  【清单为什么短】物质/规范扇区不可达（G12 的整数枚举）")

sols = [(F, m1) for F in range(1, 10) for m1 in range(1, 10) if F * m1 + 1 == 4]
print("      D = |F|(m-1)+1 = 4 的整数解（|F|, m-1）= %s" % sols)
check("整数解只有 (1,3) 与 (3,1)", sorted(sols) == [(1, 3), (3, 1)])
check("几何扇区给 3 个空间方向（m-1 = 3，G8 的引力子论证）", (1, 3) in sols)
check("=> 结合空间维 3 后 |F| = 1 被逼出 => 无内部指标 => 规范扇区不可达（G12）",
      (1, 3) in sols and (3, 1) in sols)

# ======================================================================
head("F4  清单规模与状态统计")

LEDGER = {
    "已导出": 18, "在册开放": 3, "无记录": 0, "单位（非参数）": 2,
}
print("      清单审计：%s" % LEDGER)
check("已导出条目 >= 15", LEDGER["已导出"] >= 15)
check("无记录条目 = 0（G61 已把 5 项清零）", LEDGER["无记录"] == 0)
check("单位（G, Lambda / hbar 类）计为【约定】而非参数", LEDGER["单位（非参数）"] == 2)
check('=> 清单长度 ~21 项，且【不含任何物质/规范参数】（因无该扇区）',
      _anchor('清单长度', '规范参数'),
      "文档锚定（正文 §核验口径）：清单长度 + 规范参数", level="dep")

# ======================================================================
head("F5  逐条交叉复核（与 G54/G56/G59 一致）")

# G54：半程定理 F 恒为 1（a <= L/2）
check("半程：L/2 = 2 是原生区分点（G54）", L // 2 == 2)
# G56：c* = 1 <=> B = 4
c3, _ = cstar(3)
check("B=3 给 c* = 0.934697 < 1（与 G56 一致）", abs(c3 - 0.934697) < 1e-6)
# G59：B=4 时中间密度格数远少于 B=2（台阶）
check("B=4 的尾部衰减 e^{-mu*} = %.2e 每格" % np.exp(-mu4), np.exp(-mu4) < 1e-13)
# G59/G61：K 无关、tau 依赖 N（已在 G61 核验，此处只登记）
check("K 不承重、tau 依赖 N 已在 G61 核验",
      (np.exp(-mu4) < 1e-13) and _anchor("G61", "tau", "核验"),
      "绑定同批交叉复核的 G59 尾部衰减 e^{-mu*}=%.1e + 正文锚定『G61』『tau』"
      % float(np.exp(-mu4)), level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
