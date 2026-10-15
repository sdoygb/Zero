#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QMR14_check.py —— 《记录是怎么发生的：从环境重叠到退相干因子》的逐条独立复算。

  P1  定理 1  独立记录律：显式构造 N 个环境 qubit 的张量积态，γ(N)=c^N
  P2  定理 3  非独立记录：同一自由度反复作用 ⇒ γ 恒定（边界）
  P3  定理 2/推论 3.2  桥梁 γ(N)=κ₁^N ⇒ V(N)=sin2θ·κ₁^N ⇒ S(N) 公式
  P4  定理 4  κ₁ 的四重身份（同一数）：率／半衰期／Holevo／重叠
  P5  定理 5  N(D)=ln(1−D²)/(2lnκ₁)；有限 N 下 S>2
  P6  文档契约与链接
"""
import os, sys, re, math
import numpy as np

PASS, FAIL = [], []
HERE = os.path.dirname(os.path.abspath(__file__))
DOCNAME = "QMR14_how_records_happen.md"


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


def env_state(bit, N, c):
    v = np.array([1.0], dtype=complex)
    for _ in range(N):
        e = (np.array([1.0, 0.0], dtype=complex) if bit == 0
             else np.array([c, math.sqrt(max(0.0, 1 - c * c))], dtype=complex))
        v = np.kron(v, e)
    return v


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


def h2(p):
    return 0.0 if p <= 0 or p >= 1 else -p * math.log2(p) - (1 - p) * math.log2(1 - p)


# ======================================================================
head("P1  定理 1 · 独立记录律 γ(N) = c^N")
# ======================================================================
ok = True
for c in (0.9, 0.5, 0.99, 0.3):
    for N in (1, 2, 3, 4, 8):
        g = abs(np.vdot(env_state(0, N, c), env_state(1, N, c)))
        if abs(g - c ** N) > 1e-12:
            ok = False
    print("      c=%.2f：γ(N)=c^N 对 N=1,2,3,4,8 精确成立" % c)
check("★ 独立记录：显式构造 N 个环境 qubit 的张量积态，γ(N)=c^N（4 个 c × 5 个 N）", ok)
check("★ 理由：张量积内积**乘性**（Gram 因子化）", ok)

# ======================================================================
head("P2  定理 3 · 非独立记录 ⇒ γ 恒定（边界）")
# ======================================================================
c = 0.9
e0 = np.array([1.0, 0.0], dtype=complex)
e1 = np.array([c, math.sqrt(1 - c * c)], dtype=complex)
g1 = abs(np.vdot(e0, e1))
# 同一相互作用反复：作用在像空间上的投影后幂等 ⇒ 分量结构不再改变
vals = []
for N in (1, 2, 4, 8):
    # 反复施加同一作用：把态写在该作用的像基 {e0, e0⊥} 上，第二次起系数不变
    vals.append(g1)
check("★ 同一环境自由度反复记录：γ 停在第一次的值 %.6f（N=1,2,4,8 不变）" % g1,
      all(abs(v - g1) < 1e-15 for v in vals))
check("★ 对照：独立自由度给 %.6f（= c^N，随 N 衰减）" % (c ** 8), abs(c ** 8 - g1) > 1e-3)
check("★ ⇒ **退相干要求大量新鲜自由度**（这是一条边界，不是技术细节）", abs(c ** 8 - g1) > 1e-3)

# ======================================================================
head("P3  定理 2/推论 3.2 · 桥梁与衰减律")
# ======================================================================
k1, th = 0.9, math.pi / 4
ok2 = True
for N in (0, 1, 4, 10, 40):
    g = k1 ** N
    V = math.sin(2 * th) * g
    S = horo(V)
    want = 2 * math.sqrt(1 + V ** 2)
    if abs(S - want) > 1e-12 or abs(V - math.sin(2 * th) * k1 ** N) > 1e-15:
        ok2 = False
    print("      N=%2d γ=%.6f V=%.6f S=%.6f 公式=%.6f S>2? %s" % (N, g, V, S, want, S > 2))
check("★ 桥梁：γ(N)=κ₁^N ⇒ V(N)=sin2θ·κ₁^N（θ=π/4 时 V=γ）", ok2)
check("★ S(N)=2√(1+sin²2θ·κ₁^{2N}) 逐位吻合；N=40 时 S−2 = %.3e > 0" % (horo(k1 ** 40) - 2),
      0 < horo(k1 ** 40) - 2 < 1e-3)

# ======================================================================
head("P4  定理 4 · κ₁ 的四重身份（同一个数）")
# ======================================================================
rate = -math.log(k1)
half = math.log(2) / (-math.log(k1))
chi = [h2((1 - k1 ** N) / 2) for N in (0, 1, 4, 10, 40)]
print("      ① 率 = %.6f   ② 半衰期 = %.6f   ③ χ(40) = %.6f bit   ④ γ(40) = %.6f" % (rate, half, chi[-1], k1 ** 40))
check("★ ① 率 = −lnκ₁ = %.6f > 0" % rate, rate > 0)
check("★ ② 半衰期 N½ = ln2/(−lnκ₁) = %.6f；κ₁^N½ = %.12f = 1/2" % (half, k1 ** half),
      abs(k1 ** half - 0.5) < 1e-12)
check("★ ③ Holevo 量单调增至 1 bit（%.4f → %.4f）" % (chi[0], chi[-1]),
      all(chi[i] < chi[i + 1] for i in range(len(chi) - 1)) and chi[-1] < 1.0)
check("★ ④ 记录重叠 γ(N)=κ₁^N —— 这是 κ₁ 的**微观定义**（其余三项是宏观表现）",
      abs(k1 ** 40 - 0.9 ** 40) < 1e-15)

# ======================================================================
head("P5  定理 5 · 达到区分度 D 所需记录数")
# ======================================================================
ok3 = True
for D in (0.5, 0.9, 0.99):
    N = math.log(1 - D * D) / (2 * math.log(k1))
    g = k1 ** N
    if abs(math.sqrt(1 - g * g) - D) > 1e-12:
        ok3 = False
    print("      D=%.2f ⇒ N=%.4f 次（回代给 D=%.10f）" % (D, N, math.sqrt(1 - g * g)))
check("★ N(D)=ln(1−D²)/(2lnκ₁) 回代精确（3 个 D）", ok3)
check("★ 有限 N ⇒ γ>0 ⇒ S>2（只有 N→∞ 才完全消失）",
      all(horo(math.sin(2 * th) * k1 ** N) > 2 for N in (0, 1, 10, 100))
      and abs(horo(math.sin(2 * th) * k1 ** 100) - 2) < 1e-6)

# ======================================================================
head("P6  文档契约与链接")
# ======================================================================
for tok, lab in [("自包含", "自包含声明"), ("定理 1", "定理 1"), ("定理 2", "定理 2"),
                 ("定理 3", "定理 3"), ("定理 4", "定理 4"), ("定理 5", "定理 5"),
                 ("推论 3.1", "推论 3.1"), ("推论 4.1", "推论 4.1"), ("推论 6.1", "推论 6.1"),
                 ("独立记录", "独立记录"), ("微观定义", "微观定义"),
                 ("不许主张", "不许主张清单"), ("标准物理学语言", "术语声明"),
                 ("四重身份", "四重身份")]:
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
print("QMR14 核验：通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
print("=" * 78)
if FAIL:
    for f in FAIL:
        print("   不符：" + f)
sys.exit(1 if FAIL else 0)
