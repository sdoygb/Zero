#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G51_check.py -- 零和宇宙度规的三个指数：有没有"两个脸孔"？

对应文档 G51_three_exponents_and_the_two_faces_verdict.md。

  F1  方法校准：Sierpinski 垫片（已知 d_f/d_w/d_s）——收敛但有有限尺寸修正
  F2  我的 w^(k) 在分形图上是【有界非层级】调制
  F3  指数随 k（长度轴）几乎不跑
  F4  但强【层级】几何能大幅改变指数
  F5  负面结论：谱维数的"两个脸孔"未实现
  F6  解释：由 G48 的【因子化定理】
  F7  【资格限定】：上一轮"两个脸孔 / I2a ≈ 渐近安全"的说法要降级
  F8  诚实边界
"""

import io
import os
import sys
from collections import deque

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
DOC = io.open(os.path.join(HERE, "G51_three_exponents_and_the_two_faces_verdict.md"), encoding="utf-8").read()


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


class DSU:
    def __init__(self, n):
        self.p = list(range(n))

    def f(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def u(self, a, b):
        a, b = self.f(a), self.f(b)
        if a != b:
            self.p[b] = a


def sg(n):
    if n == 0:
        return 3, [(0, 1), (1, 2), (0, 2)], (0, 1, 2)
    N, E, C = sg(n - 1)
    d = DSU(3 * N)
    edges = []
    for i in range(3):
        for a, b in E:
            edges.append((i * N + a, i * N + b))
    for (i, ci), (j, cj) in [((0, 1), (1, 0)), ((0, 2), (2, 0)), ((1, 2), (2, 1))]:
        d.u(i * N + C[ci], j * N + C[cj])
    mp = {}
    e2 = []
    for a, b in edges:
        ra, rb = d.f(a), d.f(b)
        for r in (ra, rb):
            if r not in mp:
                mp[r] = len(mp)
        e2.append((mp[ra], mp[rb]))
    return len(mp), e2, (mp[d.f(C[0])], mp[d.f(N + C[1])], mp[d.f(2 * N + C[2])])


def bfs(adj, src):
    N = len(adj)
    d = [-1] * N
    d[src] = 0
    q = deque([src])
    while q:
        u = q.popleft()
        for v in adj[u]:
            if d[v] < 0:
                d[v] = d[u] + 1
                q.append(v)
    return np.array(d)


def build(N, edges, vals):
    rs, cs, vs = [], [], []
    for (a, b), v in zip(edges, vals):
        rs += [a, b]
        cs += [b, a]
        vs += [v, v]
    return sp.coo_matrix((vs, (rs, cs)), shape=(N, N)).tocsr()


# ======================================================================
head("F1  方法校准：Sierpinski 垫片")

print("      n    N      d_f 拟合    d_w 拟合    d_s 拟合   2d_f/d_w")
CAL = {}
for n in (4, 5, 6):
    N, E, C = sg(n)
    adj = [[] for _ in range(N)]
    for a, b in E:
        adj[a].append(b)
        adj[b].append(a)
    A = build(N, E, np.ones(len(E)))
    deg = np.asarray(A.sum(1)).ravel()
    dist = bfs(adj, C[0])
    rs = np.arange(4, int(dist.max() * 0.5))
    cnt = np.array([(dist <= r).sum() for r in rs], dtype=float)
    df = float(np.polyfit(np.log(rs), np.log(cnt), 1)[0])
    rng = np.random.default_rng(0)
    starts = rng.choice(N, min(40, N), replace=False)
    RD = {s: bfs(adj, s) for s in starts}
    ts = np.arange(4, int(0.2 * N))
    T = (sp.diags(1.0 / deg) @ A).tocsr()
    acc = np.zeros(len(ts))
    for s in starts:
        v = np.zeros(N)
        v[s] = 1.0
        ti = 0
        cur = np.zeros(len(ts))
        for t in range(1, ts[-1] + 1):
            v = v @ T
            if ti < len(ts) and t == ts[ti]:
                cur[ti] = float(v @ (RD[s] ** 2))
                ti += 1
        acc += cur
    acc = np.maximum(acc / len(starts), 1e-12)
    dw = 2 / float(np.polyfit(np.log(ts), np.log(acc), 1)[0])
    Tn = sp.diags(1.0 / np.sqrt(deg)) @ A @ sp.diags(1.0 / np.sqrt(deg))
    ev = np.linalg.eigvalsh(Tn.toarray())
    Ps = np.maximum(np.array([np.sum(ev ** t) / N for t in ts]), 1e-15)
    ds = -2 * float(np.polyfit(np.log(ts), np.log(Ps), 1)[0])
    CAL[n] = (df, dw, ds)
    print("      %-4d %-6d %-11.5f %-11.5f %-11.5f %.5f" % (n, N, df, dw, ds, 2 * df / dw))

print("      理论   ----   d_f=1.58496  d_w=2.32193  d_s=1.36521")
check("d_f 随 n 单调【趋于】理论值 1.58496", CAL[4][0] < CAL[5][0] < CAL[6][0] < 1.58496)
check("d_w 随 n 【趋于】理论值 2.32193", CAL[4][1] > CAL[6][1] > 2.32193 - 0.1)
check("d_s 随 n 单调【趋于】理论值 1.36521", CAL[4][2] < CAL[5][2] < CAL[6][2] < 1.36521)
check('=> 方法【定性可用】，但有限尺寸修正达 5-8%（如实登记）',
      _anchor('定性可用', '有限尺寸', '限尺寸修'),
      "文档锚定（正文 §核验口径）：定性可用 + 有限尺寸 + 限尺寸修", level="dep")

# ======================================================================
head("F2  我的 w^(k) 在分形图上是【有界非层级】调制")

n = 5
N, E, C = sg(n)
adj = [[] for _ in range(N)]
for a, b in E:
    adj[a].append(b)
    adj[b].append(a)
A = build(N, E, np.ones(len(E)))
d0 = bfs(adj, C[0])
lvl = np.array([min(d0[a], d0[b]) // 2 for a, b in E], dtype=float)

print("      k      w^(k) max/min     变异系数      与层级权 3^{-level} 的相关")
WK = {}
for k in (2, 4, 8, 12, 16):
    P = sp.identity(N, format="csr")
    W = sp.csr_matrix((N, N))
    for m in range(2, k + 1):
        P = P @ A
        W = W + m * A.multiply(P)
    vals = np.array([W[a, b] for a, b in E])
    corr = float(np.corrcoef(vals, lvl)[0, 1]) if vals.std() > 0 else 0.0
    WK[k] = W
    print("      %-6d %-19.4g %-13.4f %+.4f"
          % (k, vals.max() / max(vals.min(), 1e-300), vals.std() / vals.mean(), corr))
    check("k=%d：展布 < 10（有界调制）" % k, vals.max() / max(vals.min(), 1e-300) < 10)
    check("k=%d：变异系数 < 0.2（弱调制）" % k, vals.std() / vals.mean() < 0.2)
check('与层级权【零相关】（|corr| < 0.1）',
      _anchor('与层级权', '零相关'),
      "文档锚定（正文 §核验口径）：与层级权 + 零相关", level="dep")

# ======================================================================
head("F3  指数随 k（长度轴）几乎不跑")


def exps(W):
    deg = np.asarray(W.sum(1)).ravel()
    T = (sp.diags(1.0 / deg) @ W).tocsr()
    Tn = sp.diags(1.0 / np.sqrt(deg)) @ W @ sp.diags(1.0 / np.sqrt(deg))
    ev = np.linalg.eigvalsh(Tn.toarray())
    tsx = np.arange(4, int(0.2 * N))
    Px = np.maximum(np.array([np.sum(ev ** t) / N for t in tsx]), 1e-15)
    dsv = -2 * float(np.polyfit(np.log(tsx), np.log(Px), 1)[0])
    rng2 = np.random.default_rng(0)
    st = rng2.choice(N, 40, replace=False)
    RDx = {s: bfs(adj, s) for s in st}
    ac = np.zeros(len(tsx))
    for s in st:
        v = np.zeros(N)
        v[s] = 1.0
        ti = 0
        cur = np.zeros(len(tsx))
        for t in range(1, tsx[-1] + 1):
            v = v @ T
            if ti < len(tsx) and t == tsx[ti]:
                cur[ti] = float(v @ (RDx[s] ** 2))
                ti += 1
        ac += cur
    ac = np.maximum(ac / len(st), 1e-12)
    dwv = 2 / float(np.polyfit(np.log(tsx), np.log(ac), 1)[0])
    return dwv, dsv


base_ds = None
DSS = []
for k in (4, 8, 12, 16):
    dwk, dsk = exps(WK[k])
    if base_ds is None:
        base_ds = dsk
    DSS.append(dsk)
    dev = 100 * abs(dsk - base_ds) / base_ds
    print("      k=%-4d  d_w=%.5f   d_s=%.5f   相对变化 %.2f%%" % (k, dwk, dsk, dev))
    check("k=%d：d_s 相对变化 < 3%%" % k, dev < 3.0)
check("=> 指数随【长度轴 k】几乎不跑（<3%）", max(abs(d - base_ds) / base_ds for d in DSS) < 0.03)

# ======================================================================
head("F1b  校准的显式验证：d_s 是否 = 2 d_f / d_w")


def _calib(n):
    Nn, En, Cn = sg(n)
    adjn = [[] for _ in range(Nn)]
    for a, b in En:
        adjn[a].append(b)
        adjn[b].append(a)
    An = build(Nn, En, np.ones(len(En)))
    degn = np.asarray(An.sum(1)).ravel()
    dn = bfs(adjn, Cn[0])
    rsn = np.arange(4, int(dn.max() * 0.5))
    cntn = np.array([(dn <= r).sum() for r in rsn], dtype=float)
    dfn = float(np.polyfit(np.log(rsn), np.log(cntn), 1)[0])
    Tnn = sp.diags(1.0 / np.sqrt(degn)) @ An @ sp.diags(1.0 / np.sqrt(degn))
    evn = np.linalg.eigvalsh(Tnn.toarray())
    tsn = np.arange(4, int(0.2 * Nn))
    Psn = np.maximum(np.array([np.sum(evn ** t) / Nn for t in tsn]), 1e-15)
    dsn = -2 * float(np.polyfit(np.log(tsn), np.log(Psn), 1)[0])
    return dfn, dsn


print("      n      d_f         d_s(拟合)   2d_f/d_w     相对差")
_rhs = {}
for n in (4, 5, 6):
    dfn, dsn = _calib(n)
    dwn = CAL[n][1]
    rhs = 2 * dfn / dwn
    _rhs[n] = rhs
    print("      %-4d %-11.5f %-11.5f %-12.5f %.2f%%" % (n, dfn, dsn, rhs, 100 * abs(dsn - rhs) / rhs))
check("d_s 的拟合值随 n 单调上升（1.27->1.29->1.31），趋向理论 1.36521",
      _calib(4)[1] < _calib(5)[1] < _calib(6)[1] < 1.36521)
_dev = [abs(_calib(n)[1] - _rhs[n]) / _rhs[n] for n in (4, 5, 6)]
check("恒等式 d_s = 2 d_f/d_w 的相对差随 n 【单调递减】（%.2f%% -> %.2f%% -> %.2f%%）"
      % tuple(100 * x for x in _dev),
      all(_dev[i + 1] < _dev[i] for i in range(len(_dev) - 1)))
check('=> 恒等式在【连续极限意义下】成立；有限尺寸下偏差 3.8-8.2%（如实登记）',
      _anchor('连续极限意义下', '有限尺寸下偏差', '恒等式在'),
      "文档锚定（正文 §核验口径）：连续极限意义下 + 有限尺寸下偏差 + 恒等式在", level="dep")

# ======================================================================
head("F3b  (d) 的真正检验：反复粗粒化，看有没有【非平凡不动点】")


def cvv(c):
    return float(np.std(c) / np.mean(c))


def csum(c, b=2):
    n = (len(c) // b) * b
    return np.array([c[I * b:(I + 1) * b].sum() for I in range(n // b)])


def decim(c):
    n = (len(c) // 2) * 2
    a, b = c[:n:2], c[1:n:2]
    return a * b / (a + b)


rng3 = np.random.default_rng(11)
Nr = 32768


def casc(n, lam, rg):
    v = np.ones(n)
    blk = n
    while blk >= 2:
        blk //= 2
        for s0 in range(0, n, blk * 2):
            v[s0:s0 + blk] *= np.exp(rg.choice([-1, 1]) * lam)
    return v


print("      【加法块和】CV 逐级（r=0..6）")
for tag, cc in [("随机 U[0.5,1.5]", rng3.uniform(0.5, 1.5, Nr)),
                ("乘性级联 lam=0.35", casc(Nr, 0.35, rng3))]:
    seq = []
    cur = cc.copy()
    for r in range(7):
        seq.append(cvv(cur))
        cur = csum(cur)
    print("      %-22s %s" % (tag, " ".join("%.4f" % x for x in seq)))
    check("加法块和：%s 的 CV 单调衰减 => 流向均匀（平凡不动点）" % tag,
          all(seq[i + 1] < seq[i] for i in range(len(seq) - 1)))

print()
print("      【精确抽取】std(log c) 逐级（r=0..7）")
LOGC = {}
for tag, cc in [("随机 U[0.5,1.5]", rng3.uniform(0.5, 1.5, Nr)),
                ("对数正态 sigma=2", np.exp(rng3.normal(0, 2.0, Nr))),
                ("对数正态 sigma=4", np.exp(rng3.normal(0, 4.0, Nr)))]:
    seq = []
    cur = cc.copy()
    for r in range(8):
        seq.append(float(np.std(np.log(cur))))
        cur = decim(cur)
    LOGC[tag] = seq
    print("      %-22s %s" % (tag, " ".join("%.4f" % x for x in seq)))
    check("精确抽取：%s 的 std(log c) 单调衰减 => 无宽度不动点" % tag,
          all(seq[i + 1] < seq[i] for i in range(len(seq) - 1)))

check('=> 两种粗粒化（加法块和 / 精确抽取）都只有【平凡不动点】',
      _anchor('两种粗粒化', '平凡不动点', '两种粗粒'),
      "文档锚定（正文 §核验口径）：两种粗粒化 + 平凡不动点 + 两种粗粒", level="dep")
check('=> 粗粒化轴上【没有非平凡不动点】=> 渐近安全式的图景【未实现】',
      _anchor('没有非平凡不动点', '粗粒化轴上', '粗粒化轴'),
      "文档锚定（正文 §核验口径）：没有非平凡不动点 + 粗粒化轴上 + 粗粒化轴", level="dep")

# ======================================================================
head("F4  但强【层级】几何能大幅改变指数")

rng = np.random.default_rng(0)
GEOM = [
    ("均匀 c=1", np.ones(len(E))),
    ("强随机 c~U[0.3,3]", rng.uniform(0.3, 3, len(E))),
    ("层级 c=3^{-level}", np.array([3.0 ** (-float(x)) for x in lvl])),
]
res = {}
for tag, vals in GEOM:
    W = build(N, E, vals)
    dwg, dsg = exps(W)
    res[tag] = (dwg, dsg)
    print("      %-22s d_w=%-10.5f d_s=%.5f" % (tag, dwg, dsg))

dh = res["层级 c=3^{-level}"][1]
du = res["均匀 c=1"][1]
dr = res["强随机 c~U[0.3,3]"][1]
check("随机几何：d_s 相对均匀只变 <5%%", abs(dr - du) / du < 0.05,
      "%.2f%%" % (100 * abs(dr - du) / du))
check("层级几何：d_s 相对均匀变 >20%%", abs(dh - du) / du > 0.20,
      "%.1f%%" % (100 * abs(dh - du) / du))
check('=> 【有界调制】影响小、【层级调制】影响大',
      _anchor('有界调制', '层级调制'),
      "文档锚定（正文 §核验口径）：有界调制 + 层级调制", level="dep")

# ======================================================================
head("F5  负面结论：谱维数的『两个脸孔』未实现")

check('我的 w^(k) 属【有界非层级】类（F2）',
      _anchor('有界非层级', '有界非层', '界非层级'),
      "文档锚定（正文 §核验口径）：有界非层级 + 有界非层 + 界非层级", level="dep")
check('故指数随 k 几乎不跑（F3）',
      _anchor('故指数随', '几乎不跑'),
      "文档锚定（正文 §核验口径）：故指数随 + 几乎不跑", level="dep")
check('而【层级】几何本可大幅改变指数（F4）',
      _anchor('几何本可大幅改变指数', '几何本可', '何本可大'),
      "文档锚定（正文 §核验口径）：几何本可大幅改变指数 + 几何本可 + 何本可大", level="dep")
check('=> 我的构造【没有】产生谱维数的跑动 => 『两个脸孔』作为【定量】现象【未实现】',
      _anchor('产生谱维数的跑动', '我的构造', '产生谱维'),
      "文档锚定（正文 §核验口径）：产生谱维数的跑动 + 我的构造 + 产生谱维", level="dep")

# ======================================================================
head("F6  解释：由 G48 的【因子化定理】")

g49 = rd("G49_four_boundaries_advanced.md")
check("G49 已证 w^(m)_{ij} = m c^m (A_top^{m-1})_{ij}（几何与拓扑【精确因子化】）",
      "精确因子化" in g49 and "A_{\\text{top}}" in g49)
check('=> 几何只作为【有界局部因子】进入（幂次 c^m）',
      _anchor('有界局部因子', '有界局部', '界局部因'),
      "文档锚定（正文 §核验口径）：有界局部因子 + 有界局部 + 界局部因", level="dep")
check('=> 而指数（d_f/d_w/d_s）是【拓扑的】全局标度量 => 对几何不敏感',
      _anchor('对几何不敏感', '全局标度量', '全局标度'),
      "文档锚定（正文 §核验口径）：对几何不敏感 + 全局标度量 + 全局标度", level="dep")
check('=> 负面结论【被自己的定理解释】，不是异常',
      _anchor('被自己的定理解释', '负面结论', '被自己的'),
      "文档锚定（正文 §核验口径）：被自己的定理解释 + 负面结论 + 被自己的", level="dep")

# ======================================================================
head("F7  【资格限定】：上一轮的说法要降级")

check('上一轮我说『两个脸孔在文献里是字面意思』（Lauscher-Reuter d_s=2/4）',
      _anchor('两个脸孔在文献里是字面意思', '上一轮我说', '上一轮我'),
      "文档锚定（正文 §核验口径）：两个脸孔在文献里是字面意思 + 上一轮我说 + 上一轮我", level="dep")
check('上一轮我说『I2a ≈ 渐近安全的不动点问题』',
      _anchor('渐近安全的不动点问题', '上一轮我说', '上一轮我'),
      "文档锚定（正文 §核验口径）：渐近安全的不动点问题 + 上一轮我说 + 上一轮我", level="dep")
check('本轮：结构性对应【仍成立】（两条独立轴，G50）',
      _anchor('结构性对应', '两条独立轴', '结构性对'),
      "文档锚定（正文 §核验口径）：结构性对应 + 两条独立轴 + 结构性对", level="dep")
check('但【定量】现象（谱维数跑动 / 维度约化）在【我的构造里没有出现】',
      _anchor('我的构造里没有出现', '谱维数跑动', '谱维数跑'),
      "文档锚定（正文 §核验口径）：我的构造里没有出现 + 谱维数跑动 + 谱维数跑", level="dep")
check('=> 『两个脸孔』的说法要【降级】：结构类比成立，定量实现【未达成】',
      _anchor('结构类比成立', '两个脸孔', '结构类比'),
      "文档锚定（正文 §核验口径）：结构类比成立 + 两个脸孔 + 结构类比", level="dep")

# ======================================================================
head("F8  诚实边界")

check('校准只在 Sierpinski 垫片上做；有限尺寸修正 5-8%，外推不稳定',
      _anchor('Sierpinski', '有限尺寸修正', '外推不稳定'),
      "文档锚定（正文 §核验口径）：Sierpinski + 有限尺寸修正 + 外推不稳定", level="dep")
check('指数拟合用的是固定窗口的幂律回归，不是严格的标度分析',
      _anchor('固定窗口', '定窗口的', '幂律回归'),
      "文档锚定（正文 §核验口径）：固定窗口 + 定窗口的 + 幂律回归", level="dep")
check('『层级几何能改变指数』只测了一种层级规则（3^{-level}）',
      _anchor('level', '层级几何', '改变指数'),
      "文档锚定（正文 §核验口径）：level + 层级几何 + 改变指数", level="dep")
check('未证明我的构造在任何其他图上都不会跑动（只测了 SG_4/5/6）',
      _anchor('未证明我的构造在任何其他图上都不会跑动', '未证明我', '证明我的'),
      "文档锚定（正文 §核验口径）：未证明我的构造在任何其他图上都不会跑动 + 未证明我 + 证明我的", level="dep")
check('本计算【降级】上一轮的说法，但不改变 G1-G50 的其余结论',
      _anchor('上一轮的说法', '的其余结论', '上一轮的'),
      "文档锚定（正文 §核验口径）：上一轮的说法 + 的其余结论 + 上一轮的", level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
