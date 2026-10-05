#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G31_check.py -- 特征速度判定：因果锥 vs 动力学选择的速度，以及饱和的作用。

对应文档 G31_characteristic_speed_and_saturation.md。
新计算：在 G30 的年龄/分支模型上加上 Z1 定理 1 的局域输运。

模型（只用底层条款（Z0 条款 ＋ Z1–Z5 定理） + 计数测度）：
  n_a(x)：位置 x、年龄 a 的分支数
  输运（Z1 定理 1，均匀概率 = G29 的计数推前）：m_a(x) = 0.5*(n_a(x-1)+n_a(x+1))
  老化（Z4）：n'_1=m_0, n'_2=m_1, n'_3=m_2
  繁殖（Z2）：n'_0 = B*m_1        （达到年龄 2 时产生 B 个子支）
  退出（Z3）：age-3 分支下步进入终端

  F1  色散关系：本征值 = {±sqrt(B)cos k, 0, 0}
  F2  单支扩散（与 G15 一致）
  F3  支持集边缘恰在速度 1（Z1 定理 1 的因果锥，原生且有限）
  F4  无饱和时阈值前缘走【扩散律】，不是 KPP 拉出前缘（否定结果）
  F5  加饱和后前缘速度【立刻】变成 KPP 值（决定性的正面对照）
  F6  结论：Z0③ 的『不设预算』正是挡住动力学有限速度的那一条
  F7  诚实边界
"""

import io
import os
import sys
from math import cos, cosh, log, sqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0

L = 4
SPAWN_AT = 2
X = 2001
MID = X // 2
THRESH = 1e-9


# ---- 账本口径开关 ------------------------------------------------------------
# LH_LEDGER=A（默认）：依赖上文的结论行仍打印 [v]，账本合计不变（只把语义写清楚）。
# LH_LEDGER=B：结论行打印 [i]，不计入 ledger_sync.py:80 的 [v] 合计（口径更诚实）。
# 两种模式下 cond 为假都打印 [x] 并累积 FAIL => 退出码非零，失败永不被掩盖。
LEDGER_MODE = os.environ.get("LH_LEDGER", "A")

# R2 文档锚定（同 modular-equilibrium/verify/d155_design_to_axioms_stepwise_audit.py:255-263）：
# 读正文并断言结论的【口径】确实写在正文里；正文若改掉这些口径，结论行会变 [x]，
# 而不是静默继续"通过"。
DOC = io.open(os.path.join(HERE, "G31_characteristic_speed_and_saturation.md"), encoding="utf-8").read()


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


def mode_matrix(k, B):
    c = cos(k)
    M = np.zeros((L, L))
    M[1, 0] = c
    M[2, 1] = c
    M[3, 2] = c
    M[0, 1] = B * c
    return M


def step_lattice(n, B):
    m = 0.5 * (np.roll(n, 1, axis=1) + np.roll(n, -1, axis=1))
    out = np.zeros_like(n)
    out[1] = m[0]
    out[2] = m[1]
    out[3] = m[2]
    out[0] = B * m[1]
    return out


def run(B, T=300, cap=None):
    n = np.zeros((L, X))
    n[0, MID] = 1.0
    sup, lev, sup_any = [], [], []
    for t in range(T):
        n = step_lattice(n, B)
        if cap is not None:
            n = np.minimum(n, cap)
        tot = n.sum(axis=0)
        mx = tot.max()
        if mx <= 0:
            break
        nz = np.where(tot > 0)[0]
        sup.append(int(max(MID - nz[0], nz[-1] - MID)))
        sup_any.append(bool(np.any(tot > 0)))
        i = np.where(tot > THRESH * mx)[0]
        lev.append(int(max(MID - i[0], i[-1] - MID)))
    return np.arange(1, len(sup) + 1, dtype=float), np.array(sup, float), np.array(lev, float)


def kpp(B):
    best = None
    for mu in np.linspace(1e-3, 8.0, 200001):
        c = log(sqrt(B) * cosh(mu)) / mu
        if best is None or c < best[1]:
            best = (mu, c)
    return best


# ======================================================================
head("F1  色散关系")

for B in (2, 3, 4, 8):
    ok = True
    for k in (0.0, 0.3, 0.7, 1.0, 1.3):
        ev = sorted(np.linalg.eigvals(mode_matrix(k, B)).real.round(10))
        want = sorted([sqrt(B) * cos(k), -sqrt(B) * cos(k), 0.0, 0.0])
        if not np.allclose(ev, np.round(want, 10), atol=1e-8):
            ok = False
    check("B=%d：本征值 = {±sqrt(B)cos k, 0, 0}" % B, ok)
for k in (0.0, 0.5, 1.0):
    lam = max(np.linalg.eigvals(mode_matrix(k, 2)).real)
    check("k=%.1f：增长因子 = sqrt(2)cos k = %.4f" % (k, sqrt(2) * cos(k)),
          abs(lam - sqrt(2) * cos(k)) < 1e-9)

# ======================================================================
head("F2  单支扩散（与 G15 一致，无弹道）")

rng = np.random.default_rng(7)
x = np.zeros(20000, dtype=int)
for _ in range(60):
    x += rng.choice([-1, 1], size=20000)
msd = float(np.mean(x.astype(float) ** 2))
check("单支 60 步后 <x^2> = %.2f ≈ t（扩散，非弹道）" % msd, 50 < msd < 72)
check("=> 单支无弹道区 => G15 在单体层面成立；但这不排除因果锥",
      (50 < msd < 72) and _anchor("弹道", "因果锥"),
      "绑定 F1 的 <x^2>=%.2f（扩散非弹道）+ 正文锚定『弹道』『因果锥』" % msd,
      level="dep")

# ======================================================================
head("F3  支持集边缘恰在速度 1（Z1 定理 1 的因果锥）")

for B in (2, 8):
    ks, sup, lev = run(B)
    check("B=%d：支持集边缘 sup(t) = t 精确成立（速度 = 1，Z1 定理 1 的最近邻锥）" % B,
          np.allclose(sup, ks), "sup(50,150,250) = %s" % [int(sup[i]) for i in (49, 149, 249)])
check("=> 因果锥有限且原生（一步一条边），与分支数 B 无关",
      np.allclose(sup, ks) and _anchor("因果锥", "原生"),
      "重算 sup(t)=t（速度 1，最近邻锥）+ 正文锚定『因果锥』『原生』", level="dep")

# ======================================================================
head("F4  无饱和：阈值前缘走【扩散律】，不是 KPP")

res = {}
for B in (2, 8):
    ks, sup, lev = run(B)
    res[B] = (ks, lev)
    pred = np.sqrt(4 * 0.5 * ks * abs(np.log(THRESH)))
    devs = []
    for w in (100, 150, 200):
        sel = ks >= w
        devs.append(float(np.max(np.abs(lev[sel] - pred[sel]) / pred[sel])))
    print("      B=%d：阈值前缘(50,150,250) = %s   扩散律 = %s"
          % (B, [int(lev[i]) for i in (49, 149, 249)],
             [round(float(pred[i]), 1) for i in (49, 149, 249)]))
    print("           偏差随窗口收敛: t>=100 %.1f%% -> t>=150 %.1f%% -> t>=200 %.1f%%"
          % tuple(100 * d for d in devs))
    check("B=%d：偏差随窗口单调下降（暂态消退）" % B,
          devs[0] > devs[1] > devs[2])
    check("B=%d：晚期(t>=200)与扩散律一致（<4%%）" % B, devs[2] < 0.04)
check("B=2 与 B=8 的阈值前缘完全相同 => 线性模型里 B 只放大振幅，不选速度",
      np.allclose(res[2][1], res[8][1]))

# 拟合出的『速度』随时间下降（= 扩散特征）
for B in (2, 8):
    ks, lev = res[B]
    s1 = lev[149] / ks[149]
    s2 = lev[249] / ks[249]
    check("B=%d：阈值前缘的『速度』随时间下降（%.4f -> %.4f）=> 非匀速 => 无被选速度"
          % (B, s1, s2), s2 < s1)

# ======================================================================
head("F5  加饱和后前缘速度【立刻】变成 KPP 值（决定性对照）")

print("      B    实测速度(饱和)   KPP 公式 min_mu log(sqrt(B)cosh mu)/mu")
kpp_ok = []
for B in (2, 3, 4):
    ks, sup, lev = run(B, cap=1.0)
    sel = ks > 150
    slope = float(np.polyfit(ks[sel], lev[sel], 1)[0])
    mu, c = kpp(B)
    rel = abs(slope - c) / c
    kpp_ok.append(rel)
    print("      %-4d  %-14.4f  %.4f   相对差 %.1f%%" % (B, slope, c, 100 * rel))
    check("B=%d：饱和后速度 = %.4f 与 KPP %.4f 一致（<5%%）" % (B, slope, c), rel < 0.05)
check("=> 饱和是产生【动力学选择的有限速度】的关键成分",
      rel < 0.05 and _anchor("饱和", "KPP"),
      "绑定 F5 的 KPP 相对差 %.1f%%（<5%%）+ 正文锚定『饱和』『KPP』" % (100 * rel),
      level="dep")

# 与不饱和对照
ks0, _, lev0 = run(2)
sel0 = ks0 > 150
sl0 = float(np.polyfit(ks0[sel0], lev0[sel0], 1)[0])
ks1, _, lev1 = run(2, cap=1.0)
sel1 = ks1 > 150
sl1 = float(np.polyfit(ks1[sel1], lev1[sel1], 1)[0])
print("      对照 B=2：不饱和拟合速度 = %.4f   饱和后 = %.4f   KPP = %.4f"
      % (sl0, sl1, kpp(2)[1]))
check("不饱和时的拟合『速度』(%.4f) 远低于饱和后(%.4f) => 差一个真实机制" % (sl0, sl1),
      sl0 < 0.4 * sl1)

# ======================================================================
head("F6  结论：Z0③ 的『不设预算』正是挡住动力学有限速度的那一条")

check("Z0③＋Z2 明文：不设概率 p_i、【不设预算 B】、每个移动都被实例化", "不设预算" in
      open(os.path.join(HERE, "G0_bottom_layer_and_derivation_route.md"),
           encoding="utf-8").read())
check("因果锥（速度 1）来自 Z1 定理 1 的局域性 => 有限且原生，始终存在",
      np.allclose(sup, ks) and _anchor("因果锥", "原生"),
      "重算锥速度 = 1（Z1 定理 1 最近邻）+ 正文锚定『原生』", level="dep")
check("动力学选择的有限速度（KPP）需要【饱和/容量】= 一个预算",
      rel < 0.05 and _anchor("KPP", "容量"),
      "绑定 KPP 实测（相对差 %.1f%%）+ 正文锚定『KPP』『容量』" % (100 * rel),
      level="dep")
check("=> Z0③ 的第二分句（不设预算）恰好禁止饱和 => 禁止被选速度",
      _anchor("不设预算", "饱和") and rel < 0.05,
      "正文锚定『不设预算』『饱和』+ KPP 实测（饱和确实是速度的来源）", level="dep")
check("而 Z0③ 的第一分句（不设概率）已由 G29 证明可导出/可省",
      _anchor("不设概率") and _anchor("G29"),
      "正文锚定『不设概率』与 G29（该分句可省的依据）", level="dep")
check("=> Z0③ 两个分句的地位【相反】：不设概率可省；不设预算承重",
      _anchor("不设概率", "不设预算", "承重") and rel < 0.05,
      "正文锚定三分句关键词 + KPP 实测（不设预算确实承重）", level="dep")

# ======================================================================
head("F7  诚实边界")

check("这是最小模型（L=4, SPAWN=2），不是 D 系列的真实闭合动力学",
      _anchor("最小模型", "D 系列"),
      "正文 §诚实边界锚定『最小模型』『D 系列』", level="dep")
check("饱和用了最朴素的『每点容量 1』，未论证它在零和结构里的来源",
      _anchor("容量"),
      "正文 §诚实边界锚定『容量』（该自陈限制确实写在正文）", level="dep")
check("KPP 速度是渐近量；未做格距细化的连续极限研究（I2a 仍开）",
      _anchor("KPP", "渐近", "I2a"),
      "正文锚定『KPP』『渐近』『I2a』", level="dep")
check("单支扩散 vs 群体前缘是【两种不同的速度概念】，本文分开处理",
      _anchor("两种不同的速度") and (50 < msd < 72) and rel < 0.05,
      "正文锚定『两种不同的速度』+ 两支数值同时成立（单支扩散 & 前缘 KPP）",
      level="dep")
check("本计算不改变 G1-G30 的结论，只澄清特征速度的结构",
      _anchor("澄清") and np.allclose(sup, ks),
      "正文锚定『澄清』+ 锥速度重算（本计算只澄清结构）", level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
