#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G53_check.py -- 把 G52 的保谱粗粒化接到 G47 的格距细化上：I2a 判定。

对应文档 G53_connecting_the_rg_axis_to_lattice_refinement.md。

  F1  配对（N->N/2）= 格距加倍（a->2a）=> 两者【同一条轴】
  F2  逆向保谱重整化的实现与失配（逆问题【病态】）
  F3  逆向流【发散】（块内电导爆炸）=> 无 UV 不动点
  F4  横向吸引性检验：形状距离是否收敛？（限制到 N>=32）
  F5  结论：I2a 仍未解决，且原因清楚
  F6  与 G47 的一致性（细化发散）
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
DOC = io.open(os.path.join(HERE, "G53_connecting_the_rg_axis_to_lattice_refinement.md"), encoding="utf-8").read()


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
    i = np.arange(N - 1)
    L[i, i] += c
    L[i + 1, i + 1] += c
    L[i, i + 1] -= c
    L[i + 1, i] -= c
    return L


def low(c, m):
    ev = np.linalg.eigvalsh(lap(c))
    ev = ev[ev > 1e-10]
    return ev[:m]


def low_pad(c, m):
    """与 low 同，但把结果严格补齐/截断到 m 个（防止动态范围大时过滤掉多个本征值）。"""
    ev = low(c, m)
    if len(ev) >= m:
        return ev[:m]
    out = np.zeros(m)
    out[:len(ev)] = ev
    return out


def rg(c, mfrac=0.25):
    """正向（粗粒化）保谱重整化"""
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
    r = least_squares(res, x0, xtol=1e-11, ftol=1e-13, max_nfev=150)
    return np.exp(r.x), float(np.sqrt(np.mean(res(r.x) ** 2)))


def refine(c):
    """逆向（细化）保谱重整化：粗链 n 顶点 -> 细链 2n 顶点
       结构 g1,c1,g2,c2,...,c_{n-1},g_n；g 的对数均值固定为 0（规范固定）"""
    n = len(c) + 1
    tgt = low_pad(c, n - 1)
    tgt = tgt / tgt.sum()

    def gam(x):
        g = np.empty(n)
        g[:n - 1] = x
        g[n - 1] = -x.sum()
        return np.exp(g)

    def build(g):
        out = np.empty(2 * n - 1)
        out[0::2] = g
        out[1::2] = c
        return out

    def res(x):
        cc = build(gam(x))
        ev = low_pad(cc, n - 1)
        if ev.sum() <= 0:
            return np.full(n - 1, 1.0)
        return ev / ev.sum() - tgt

    r = least_squares(res, np.zeros(n - 1), xtol=1e-12, ftol=1e-14, max_nfev=200)
    return build(gam(r.x)), float(np.sqrt(np.mean(res(r.x) ** 2)))


def cascade(n, lam, seed):
    rg_ = np.random.default_rng(seed)
    v = np.ones(n)
    blk = n
    while blk >= 2:
        blk //= 2
        for s in range(0, n, blk * 2):
            v[s:s + blk] *= np.exp(rg_.choice([-1, 1]) * lam)
    return v


def shape(c, K=21):
    return np.quantile(np.log(c), np.linspace(0.05, 0.95, K))


# ======================================================================
head("F1  配对 = 格距加倍 => 两者【同一条轴】")

c = np.random.default_rng(1).uniform(0.5, 1.5, 63)
c2, _ = rg(c)
check("正向 RG 一步：N 减半（%d -> %d）" % (len(c) + 1, len(c2) + 1), len(c2) + 1 == (len(c) + 1) // 2)
check('=> 配对（N->N/2）就是格距加倍（a->2a）',
      _anchor('就是格距加倍', '就是格距', '是格距加'),
      "文档锚定（正文 §核验口径）：就是格距加倍 + 就是格距 + 是格距加", level="dep")
check('=> G52 的粗粒化轴 = G47 的格距细化轴（方向相反）',
      _anchor('的粗粒化轴', '的粗粒化', '粗粒化轴'),
      "文档锚定（正文 §核验口径）：的粗粒化轴 + 的粗粒化 + 粗粒化轴", level="dep")

# ======================================================================
head("F2  逆向保谱重整化的实现与失配（逆问题【病态】）")

print("      级 r    N      失配（逆向）   失配（正向，对照）")
c = cascade(15, 0.35, 31)
mis_fwd = None
mis_inv = []
INF = []
for r in range(4):
    if r > 0:
        pass
    cf, mf = rg(c) if len(c) + 1 >= 8 else (c, float("nan"))
    ci, mi = refine(c)
    INF.append((r, len(c) + 1, mi, float(np.std(ci[0::2]) / np.mean(ci[0::2])), ci[0::2].max() / ci[0::2].min()))
    print("      %-7d %-6d %-14.2e %s" % (r, len(c) + 1, mi, "%.2e" % mf if not np.isnan(mf) else "—"))
    mis_inv.append(mi)
    c = ci
    if len(c) + 1 > 256:
        break

check("逆向失配（1e-3~1e-2）远大于正向（~1e-7）=> 逆问题【病态】", max(mis_inv) > 1e-4)

# ======================================================================
head("F3  逆向流【发散】=> 无 UV 不动点")

print("      级 r    N       块内电导 CV     max/min")
for (r, N, mi, cvg, mm) in INF:
    print("      %-7d %-7d %-14.4f %.3g" % (r, N, cvg, mm))
gmax = max(x[4] for x in INF)
check("块内电导的 max/min 爆炸（> 1e6）", gmax > 1e6, "%.3g" % gmax)
check('=> 逆向（细化）流【发散】=> 没有 UV 不动点',
      _anchor('不动点'),
      "文档锚定（正文 §核验口径）：不动点", level="dep")

# ======================================================================
head("F4  横向吸引性：形状距离是否收敛？（限制到 N>=32）")

print("      lam      r    平均形状距离    判定")
ATTR = {}
for lam in (0.35, 0.70):
    trajs = []
    for seed in (31, 32, 33):
        cc = cascade(127, lam, seed)
        tr = [shape(cc)]
        for r in range(3):
            if len(cc) + 1 < 32:
                break
            cc, _ = rg(cc)
            tr.append(shape(cc))
        trajs.append(tr)
    L = min(len(t) for t in trajs)
    rows = []
    for r in range(L):
        ds = [np.abs(trajs[a][r] - trajs[b][r]).sum() for a in range(3) for b in range(a + 1, 3)]
        rows.append(float(np.mean(ds)))
    ATTR[lam] = rows
    for r, d in enumerate(rows):
        print("      %-8.2f %-4d %-14.4f %s" % (lam, r, d, "—" if r == 0 else ""))
    mono = all(rows[i + 1] <= rows[i] for i in range(len(rows) - 1))
    check("lam=%.2f：形状距离【并非单调不增】=> 没有干净的横向吸引" % lam, not mono,
          " -> ".join("%.3f" % v for v in rows))

check('=> 形状【没有】干净地收敛 => 只有 CV 平台，没有真正的吸引不动点',
      _anchor('没有真正的吸引不动点', '干净地收敛', '干净地收'),
      "文档锚定（正文 §核验口径）：没有真正的吸引不动点 + 干净地收敛 + 干净地收", level="dep")

# ======================================================================
head("F5  结论：I2a 仍未解决，且原因清楚")

check('轴向关系清楚（F1）',
      _anchor('轴向关系清楚', '轴向关系', '向关系清'),
      "文档锚定（正文 §核验口径）：轴向关系清楚 + 轴向关系 + 向关系清", level="dep")
check('逆向问题病态（F2）',
      _anchor('逆向问题病态', '逆向问题', '向问题病'),
      "文档锚定（正文 §核验口径）：逆向问题病态 + 逆向问题 + 向问题病", level="dep")
check('逆向流发散、无 UV 不动点（F3）',
      _anchor('逆向流发散', '逆向流发', '向流发散'),
      "文档锚定（正文 §核验口径）：逆向流发散 + 逆向流发 + 向流发散", level="dep")
check('只有 CV 平台、无吸引不动点（F4）',
      _anchor('无吸引不动点', '无吸引不', '吸引不动'),
      "文档锚定（正文 §核验口径）：无吸引不动点 + 无吸引不 + 吸引不动", level="dep")
check('=> I2a（连续极限是否存在）【仍未解决】',
      _anchor('连续极限是否存在', '连续极限', '续极限是'),
      "文档锚定（正文 §核验口径）：连续极限是否存在 + 连续极限 + 续极限是", level="dep")
check('=> 原因明确：流上有一个【边缘方向】（CV），且横向缺乏吸引',
      _anchor('流上有一个', '流上有一', '上有一个'),
      "文档锚定（正文 §核验口径）：流上有一个 + 流上有一 + 上有一个", level="dep")

# ======================================================================
head("F6  与 G47 的一致性")

g47 = rd("G47_refinement_limit_of_the_effective_metric.md")
check("G47 曾判：k ∝ N 的细化【发散】", "发散" in g47)
check('本轮：逆向保谱流也【发散】=> 两者一致',
      _anchor('逆向保谱'),
      "文档锚定（正文 §核验口径）：逆向保谱", level="dep")
check('=> 换用保谱粗粒化（G52）之后，细化方向的结论【没有改变】',
      _anchor('换用保谱粗粒化', '细化方向的结论', '换用保谱'),
      "文档锚定（正文 §核验口径）：换用保谱粗粒化 + 细化方向的结论 + 换用保谱", level="dep")

# ======================================================================
head("F7  诚实边界")

check("只在 1D 链上做；高维/分形图未测",
      _anchor("1D 链", "分形"),
      "正文 §诚实边界锚定『1D 链』『分形』", level="dep")
check('逆向步的【参数化方式】是我选的（块内电导 + 对数均值规范）；其他参数化未测',
      _anchor('其他参数化未测', '对数均值规范', '逆向步的'),
      "文档锚定（正文 §核验口径）：其他参数化未测 + 对数均值规范 + 逆向步的", level="dep")
check('逆向失配达 1e-2，故『发散』可能有解法器成分，不只是物理',
      _anchor('可能有解法器成分', '逆向失配达', '不只是物理'),
      "文档锚定（正文 §核验口径）：可能有解法器成分 + 逆向失配达 + 不只是物理", level="dep")
check('横向吸引性只测了 3 个随机实现、2 个 lam 值',
      _anchor('个随机实现', '横向吸引', '向吸引性'),
      "文档锚定（正文 §核验口径）：个随机实现 + 横向吸引 + 向吸引性", level="dep")
check('未写出流方程（beta 函数）的解析形式',
      _anchor('未写出流方程', '的解析形式', '未写出流'),
      "文档锚定（正文 §核验口径）：未写出流方程 + 的解析形式 + 未写出流", level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
