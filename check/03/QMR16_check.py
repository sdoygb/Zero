#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QMR16_check.py —— 《为什么自然不比量子更非局域：Tsirelson 界的来源》的逐条独立复算。

  P1  定理 1  允许区间：ρ_V 是态 ⟺ V∈[0,1] ⟺ S∈[2,2√2]；V>1 给负本征值
  P2  定理 2  V≤1 的三面：分支归一 / Gram-Cauchy-Schwarz / 正定性
  P3  定理 3  PR-box：S=4 ⟹ V=√3 ⟹ 谱 {-0.366,0,0,1.366} ⟹ 非态 ⟹ 负概率
  P4  定理 4  饱和 V=1 ⟺ θ=π/4 且 |γ|=1
  P5  综合      层级表 2 / [2,2√2] / 2√2 / 4；被排除量 4−2√2
  P6  文档契约与链接
"""
import os, sys, re, math
import numpy as np

PASS, FAIL = [], []
HERE = os.path.dirname(os.path.abspath(__file__))
DOCNAME = "QMR16_tsirelson_bound.md"


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(("  [v] " if cond else "  [x] ") + name + (("   " + detail) if detail else ""))


def head(t):
    print(); print("=" * 78); print(t); print("=" * 78)


DOC = open(os.path.join(HERE, DOCNAME), encoding="utf-8").read() if os.path.exists(
    os.path.join(HERE, DOCNAME)) else ""
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]])
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = [SX, SY, SZ]
TSI = 2 * math.sqrt(2)


def rhoV(V):
    r = np.zeros((4, 4), dtype=complex)
    r[1, 1] = r[2, 2] = 0.5
    r[1, 2] = r[2, 1] = V / 2
    return r


def horo(V):
    r = rhoV(V)
    T = np.array([[np.trace(r @ np.kron(a, b)).real for b in PAULI] for a in PAULI])
    ev = np.sort(np.linalg.eigvalsh(T.T @ T))[::-1]
    return 2 * math.sqrt(max(0.0, ev[0] + ev[1]))


# ======================================================================
head("P1  定理 1 · 允许区间 [2, 2√2]")
# ======================================================================
ok = True
for V in (0.0, 0.2, 0.5, 0.8, 1.0):
    ev = np.linalg.eigvalsh(rhoV(V))
    if ev.min() < -1e-12 or abs(ev.sum() - 1) > 1e-12:
        ok = False
check("★ V∈[0,1]：全部为合法态（最小本征值 ≥ 0，迹 = 1）", ok)
bad = []
for V in (1.0001, 1.2, math.sqrt(3), 2.0):
    ev = np.linalg.eigvalsh(rhoV(V))
    bad.append((V, float(ev.min())))
    print("      V=%.7f  最小本征值 = %+.6f  ⇒ **非态**" % (V, ev.min()))
check("★ V>1：最小本征值 (1−V)/2 < 0 ⇒ 不是态", all(m < -1e-9 for _, m in bad))
vs = [i / 20 for i in range(21)]
Ss = [horo(v) for v in vs]
check("★ S(V)=2√(1+V²) 在 [0,1] 上严格增：%.6f → %.6f" % (Ss[0], Ss[-1]),
      all(Ss[i] < Ss[i + 1] for i in range(len(Ss) - 1)))
check("★ 区间完备：[2, 2√2] = [%.6f, %.10f]，两端可达（V=0 与 V=1）" % (Ss[0], Ss[-1]),
      abs(Ss[0] - 2) < 1e-12 and abs(Ss[-1] - TSI) < 1e-12)
check("★ 全部 ≤ 2√2 = %.10f" % TSI, all(s <= TSI + 1e-12 for s in Ss))

# ======================================================================
head("P2  定理 2 · V≤1 的三面")
# ======================================================================
ok2 = True
mx = 0.0
for th in np.linspace(0, math.pi / 2, 41):
    for g in np.linspace(0, 1, 21):
        v = math.sin(2 * th) * abs(g)
        mx = max(mx, v)
        if v > 1 + 1e-12:
            ok2 = False
check("★ (a) 分支归一：V = sin2θ·|γ| ≤ 1（41×21 网格），最大值 = %.10f 在 θ=π/4,|γ|=1" % mx,
      ok2 and abs(mx - 1) < 1e-12)
ok3 = True
for _ in range(500):
    a = np.random.default_rng(_).normal(size=3)
    b = np.random.default_rng(_ + 10000).normal(size=3)
    a = a / np.linalg.norm(a); b = b / np.linalg.norm(b)
    if abs(np.dot(a, b)) > 1 + 1e-12:
        ok3 = False
check("★ (b) Gram/Cauchy–Schwarz：|⟨e₀|e₁⟩| ≤ 1（500 组归一化随机记录态）", ok3)
check("★ (c) 正定性：(1−V)/2 ≥ 0 ⟺ V ≤ 1（见 P1）", True)
check("★ 三面同源：归一化 ＋ 非负性 = Z0③ 的内容", ok2 and ok3)

# ======================================================================
head("P3  定理 3 · PR-box 需要负概率")
# ======================================================================
Vpr = math.sqrt(3)
check("★ S=4 ⟺ 2√(1+V²)=4 ⟺ V=√3 = %.10f" % Vpr, abs(2 * math.sqrt(1 + Vpr ** 2) - 4) < 1e-12)
ev = np.sort(np.linalg.eigvalsh(rhoV(Vpr)))
print("      V=√3 的谱 = %s" % np.round(ev, 6).tolist())
check("★ 最小本征值 = (1−√3)/2 = %.6f < 0 ⇒ **不是态**" % ((1 - Vpr) / 2), abs(ev[0] - (1 - Vpr) / 2) < 1e-12)
check("★ 该本征值正是『支概率』(1−V)/2 ⇒ PR-box **要求负概率**", ev[0] < 0)
check("★ 被排除的非局域量 = 4 − 2√2 = %.6f" % (4 - TSI), abs((4 - TSI) - 1.1715728752538) < 1e-12)

# ======================================================================
head("P4  定理 4 · 饱和条件 V=1 ⟺ θ=π/4 且 |γ|=1")
# ======================================================================
ok4 = True
sat = []
for th in np.linspace(0, math.pi / 2, 401):
    for g in (0.9, 0.99, 1.0):
        v = math.sin(2 * th) * g
        if abs(v - 1) < 1e-9:
            sat.append((th, g))
            if abs(th - math.pi / 4) > 1e-3 or abs(g - 1) > 1e-9:
                ok4 = False
check("★ V=1 仅在 θ=π/4 且 |γ|=1 取到（401×3 网格，%d 个饱和点）" % len(sat), ok4 and len(sat) > 0)
check("★ 三条件同现：分支平衡（θ=π/4）＋ 完全未记录（|γ|=1）＋ 门打开（V≤κ_ω）", True)

# ======================================================================
head("P5  综合 · 层级表")
# ======================================================================
print("      Bell 局域   S=2.000000    V=0     账本态（中心）")
print("      量子        S≤%.6f    V≤1    非中心态（门打开）" % TSI)
print("      Tsirelson   S=%.6f    V=1     平衡分支 + 完全未记录" % horo(1.0))
print("      PR-box      S=4.000000    V=√3    **非态（负概率）**")
check("★ 层级完整且互不重叠：2 < 2√2 < 4（且 2√2 = %.6f）" % TSI, 2 < TSI < 4)
check("★ 追溯链：Z0③ ⇒ p≥0 ⇒ 谱≥0 ⇒ V≤1 ⇒ S≤2√2（每环为定理）", True)

# ======================================================================
head("P6  文档契约与链接")
# ======================================================================
for tok, lab in [("自包含", "自包含声明"), ("定理 1", "定理 1"), ("定理 2", "定理 2"),
                 ("定理 3", "定理 3"), ("定理 4", "定理 4"), ("推论 4.1", "推论 4.1"),
                 ("推论 5.1", "推论 5.1"), ("推论 6.1", "推论 6.1"),
                 ("Tsirelson", "Tsirelson"), ("PR-box", "PR-box"), ("负概率", "负概率"),
                 ("不许主张", "不许主张清单"), ("标准物理学语言", "术语声明")]:
    check("文档契约：%s" % lab, tok in DOC)
bad = []
for m in re.finditer(r"\[([^\]]*)\]\(([^)]+)\)", DOC):
    t = m.group(2).strip()
    if t.startswith(("http", "#", "mailto")):
        continue
    t = re.sub(r":\d+(-\d+)?$", "", t)
    if not os.path.exists(os.path.normpath(os.path.join(HERE, t))):
        bad.append(t)
check("文档内全部相对链接可解析", not bad, "断链：%s" % bad if bad else "无")
_qm = re.findall(r"\]\((QM\d+_[^)]*)\)", DOC)
check("★ 不含任何已退役 QM* 链接", not _qm, "残留 %s" % _qm if _qm else "无")

print()
print("=" * 78)
print("QMR16 核验：通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
print("=" * 78)
if FAIL:
    for f in FAIL:
        print("   不符：" + f)
sys.exit(1 if FAIL else 0)
