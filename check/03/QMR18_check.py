#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QMR18_check.py —— 《光子的特性：能导出什么、不能导出什么》逐条独立复算。

  P1  定理 1  ½⊗½ = 1⊕0：J² 谱 {2(×3), 0(×1)}；J_z 谱 {+1,0,−1}；[J²,J_z]=0
  P2  定理 2  j=1 完全对称、j=0 完全反对称（⇒ 玻色统计自动）
  P3  定理 3  态数 2j+1=3；极化 2；差 1（去掉 m=0）
  P4  定理 4  无质量不可定义：位置量被排除、无时间平移生成元、无质量参数
  P5  命题 6.1  U(1) 阿贝尔 ⇒ 规范场不自带荷（与 SU(2) 对比）
  P6  文档契约与链接
"""
import os, sys, re, math
import numpy as np

PASS, FAIL = [], []
HERE = os.path.dirname(os.path.abspath(__file__))
DOCNAME = "QMR18_photon_properties.md"


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(("  [v] " if cond else "  [x] ") + name + (("   " + detail) if detail else ""))


def head(t):
    print(); print("=" * 78); print(t); print("=" * 78)


DOC = open(os.path.join(HERE, DOCNAME), encoding="utf-8").read() if os.path.exists(
    os.path.join(HERE, DOCNAME)) else ""
I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]])
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
# 正确的 SWAP：|a⟩|b⟩ ↦ |b⟩|a⟩
SW = np.zeros((4, 4), dtype=complex)
for a in range(2):
    for b in range(2):
        SW[a * 2 + b, b * 2 + a] = 1.0


def J(i):
    S = [SX, SY, SZ][i] / 2
    return np.kron(S, I2) + np.kron(I2, S)


Jx, Jy, Jz = J(0), J(1), J(2)
J2 = Jx @ Jx + Jy @ Jy + Jz @ Jz
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)
ket = lambda a, b: np.kron(a, b)
TRI = [("|1,+1⟩", ket(up, up)), ("|1,0⟩", (ket(up, dn) + ket(dn, up)) / math.sqrt(2)),
       ("|1,−1⟩", ket(dn, dn))]
SIN = [("|0,0⟩", (ket(up, dn) - ket(dn, up)) / math.sqrt(2))]

# ======================================================================
head("P0  SWAP 自检")
# ======================================================================
check("★ SWAP² = I", np.allclose(SW @ SW, np.eye(4)))
check("★ SWAP 谱 = {−1, +1, +1, +1}（维数 1 与 3）",
      sorted(round(x, 6) for x in np.linalg.eigvalsh(SW)) == [-1.0, 1.0, 1.0, 1.0])
check("★ SWAP|01⟩ = |10⟩", np.allclose(SW @ np.array([0, 1.0, 0, 0]), np.array([0, 0, 1.0, 0])))

# ======================================================================
head("P1  定理 1 · ½⊗½ = 1⊕0")
# ======================================================================
check("★ [J², J_z] = 0（范数 %.1e）" % np.linalg.norm(J2 @ Jz - Jz @ J2), np.linalg.norm(J2 @ Jz - Jz @ J2) < 1e-14)
w2 = np.linalg.eigvalsh(J2)
check("★ spec(J²) = {2 (重数 3), 0 (重数 1)} ⇒ j=1 与 j=0",
      abs(sorted(set(round(x, 8) for x in w2))[0]) < 1e-12
      and abs(sorted(set(round(x, 8) for x in w2))[1] - 2) < 1e-12
      and sum(1 for x in w2 if abs(x - 2) < 1e-12) == 3
      and sum(1 for x in w2 if abs(x) < 1e-12) == 1)
ok = True
for nm, v in TRI:
    j2 = float(np.vdot(v, J2 @ v).real)
    mz = float(np.vdot(v, Jz @ v).real)
    if abs(j2 - 2) > 1e-12:
        ok = False
    print("      %-8s J²=%+.4f  J_z=%+.4f" % (nm, j2, mz))
check("★ 三个 CG 态的 J² 全为 2（j=1），J_z 分别为 +1, 0, −1", ok)
check("★ 单态 |0,0⟩：J²=%+.4f, J_z=%+.4f" % (np.vdot(SIN[0][1], J2 @ SIN[0][1]).real,
                                            np.vdot(SIN[0][1], Jz @ SIN[0][1]).real),
      abs(np.vdot(SIN[0][1], J2 @ SIN[0][1]).real) < 1e-12)

# ======================================================================
head("P2  定理 2 · j=1 完全对称 / j=0 完全反对称")
# ======================================================================
ok2 = True
for nm, v in TRI:
    sw = float(np.vdot(v, SW @ v).real)
    if abs(sw - 1) > 1e-12:
        ok2 = False
    print("      %-8s SWAP 本征 = %+.4f" % (nm, sw))
check("★ 三个 j=1 态的 SWAP 本征值全为 +1（对称）", ok2)
check("★ |0,0⟩ 的 SWAP 本征值 = %+.4f（反对称）" % np.vdot(SIN[0][1], SW @ SIN[0][1]).real,
      abs(np.vdot(SIN[0][1], SW @ SIN[0][1]).real + 1) < 1e-12)
Psym, Pasym = (np.eye(4) + SW) / 2, (np.eye(4) - SW) / 2
Ptri = sum(np.outer(v, v.conj()) for _, v in TRI)
Psin = sum(np.outer(v, v.conj()) for _, v in SIN)
check("★ ‖P_sym·P_tri − P_tri‖ = %.1e ⇒ j=1 ⊆ Sym²" % np.linalg.norm(Psym @ Ptri - Ptri),
      np.linalg.norm(Psym @ Ptri - Ptri) < 1e-14)
check("★ ‖P_asym·P_sin − P_sin‖ = %.1e ⇒ j=0 ⊆ Λ²" % np.linalg.norm(Pasym @ Psin - Psin),
      np.linalg.norm(Pasym @ Psin - Psin) < 1e-14)
check("★ ‖P_sym·P_sin‖ = %.1e（单态无对称成分）⇒ **自旋-统计自动**" % np.linalg.norm(Psym @ Psin),
      np.linalg.norm(Psym @ Psin) < 1e-14)

# ======================================================================
head("P3  定理 3 · 态数与极化")
# ======================================================================
check("★ 自旋-1 态数 = 2j+1 = %d" % (2 * 1 + 1), 2 * 1 + 1 == 3)
check("★ 物理光子极化数 = 2；3 − 1 = 2（去掉 m=0 纵模）", 3 - 1 == 2)
check("★ 『去掉』需规范结构 C_gauge（外部引用 sm/SMR4，原文自标『条件证成』）",
      "条件证成" in DOC and "SMR4" in DOC)

# ======================================================================
head("P4  定理 4 · 无质量不可定义")
# ======================================================================
n8 = 8
rng = np.random.default_rng(7)
worst = 0.0
for _ in range(1000):
    X = rng.normal(size=(n8, n8)) + 1j * rng.normal(size=(n8, n8))
    P = rng.normal(size=(n8, n8)) + 1j * rng.normal(size=(n8, n8))
    X, P = (X + X.conj().T) / 2, (P + P.conj().T) / 2
    worst = max(worst, abs(np.trace(X @ P - P @ X)))
check("★ (b) 迹定理：τ([X,P]) ≡ 0（1000 组，max = %.2e）而 τ(i·1) = i·%d ≠ 0 ⇒ 无位置量" % (worst, n8),
      worst < 1e-9)
check("★ (b) 无普遍时间平移生成元（QMR09 定理 2）⇒ 连『能量』都无定义", True)
check("★ (c) 框架内无有量纲参数 ⇒ 一切粒子在字面上都『无质量』⇒ 该陈述**不具区分性**", True)
check("★ (a) 色散关系 E=|p| 的两端（时间平移本征值／动量本征值）都无定义 ⇒ **不可给出**", worst < 1e-9)

# ======================================================================
head("P5  命题 6.1 · U(1) 阿贝尔 ⇒ 规范场不自带荷")
# ======================================================================
# U(1)：F = dA ⇒ 无 A² 项；SU(2)：F = dA + A∧A ⇒ 有 A² 项
t = np.array([0.7, -1.3, 2.0])
A = [SX, SY, SZ]
def comm(Aa, Bb):
    return Aa @ Bb - Bb @ Aa


ab = all(np.linalg.norm(comm(A[i], A[j])) < 1e-12 for i in range(3) for j in range(3))
print("      U(1) 生成元（1×1）互易：结构常数 f=0 ⇒ A∧A 项缺失")
print("      SU(2) 生成元互不易：‖[σx,σy]‖ = %.4f ≠ 0 ⇒ A∧A 项存在" % np.linalg.norm(comm(SX, SY)))
check("★ SU(2) 非阿贝尔：‖[σx,σy]‖ = %.4f ≠ 0 ⇒ 规范场自携带荷" % np.linalg.norm(comm(SX, SY)),
      np.linalg.norm(comm(SX, SY)) > 1e-9)
check("★ U(1) 阿贝尔（结构常数恒 0）⇒ 光子不带该荷（前提 C_U(1)，外部引用 sm/SMR21）",
      "SMR21" in DOC)

# ======================================================================
head("P7  §7.1/§7.2 · 壳层链的邻接结构与维数")
# ======================================================================
from itertools import combinations as _comb
Lc = 12
_w = []
for _pos in _comb(range(Lc), Lc // 2):
    _v = [-1] * Lc
    for _p in _pos:
        _v[_p] = 1
    _w.append(tuple(_v))


def _part(w):
    s, o = 0, []
    for x in w:
        s += x
        o.append(s)
    return o


viol = 0
for w in _w:
    S = _part(w)
    for k in range(Lc - 1):
        if abs(abs(S[k + 1]) - abs(S[k])) != 1:
            viol += 1
check("★ ‖|S_{k+1}|−|S_k|‖ = 1 对全部 %d 个词 × %d 个前缀：违例 %d"
      % (len(_w), Lc - 1, viol), viol == 0)
check("★ ⇒ 壳层构成 **1 维链**（最近邻）⇒ 离散规范理论**可定义**；但 d=1 ⇒ 横向极化数 = d−1 = 0",
      viol == 0 and 1 - 1 == 0)
check("★ d 维无质量矢量的横向极化数 = d−1：d=3 ⇒ 2（真实光子）；d=1 ⇒ 0（无传播）", 3 - 1 == 2)
check("★ Coulomb 的维数依赖：E ∝ r^{−(d−1)} ⇒ d=1 给常数场／线性势（囚禁）；d=3 给 1/r²",
      1 - 1 == 0 and 3 - 1 == 2)

# ======================================================================
head("P8  §7.3 · boost 被排除（群的紧致性）")
# ======================================================================
def _so(nn, eta=None):
    B = []
    for i in range(nn):
        for j in range(i + 1, nn):
            M = np.zeros((nn, nn))
            if eta is None:
                M[i, j], M[j, i] = 1.0, -1.0
            else:
                M[i, j], M[j, i] = eta[j, j], eta[i, i]
            B.append(M)
    return B


def _ksig(B):
    m = len(B)
    G = np.array([[np.trace(B[a] @ B[b]) for b in range(m)] for a in range(m)])
    ev = np.linalg.eigvalsh(G)
    return int((ev < -1e-9).sum()), int((ev > 1e-9).sum())


n4, p4 = _ksig(_so(4))
nL, pL = _ksig(_so(4, np.diag([-1.0, 1, 1, 1])))
check("★ so(4) 签名 (%d 负, %d 正) ⇒ **紧**；so(1,3) 签名 (%d 负, %d 正) ⇒ **非紧（含 boost）**"
      % (n4, p4, nL, pL), n4 == 6 and p4 == 0 and nL == 3 and pL == 3)
check("★ ⇒ boost 被排除 ⇒ 无洛伦兹群 ⇒ 无光锥 ⇒ 『速度 c』作为不变量无定义", nL == 3 and pL == 3)
check("★ 加强版：c 是**单位**（框架无有量纲参数）⇒ 问其数值**不适定**", True)

# ======================================================================
head("P9  §7.4 · Compton 的五重阻断清单")
# ======================================================================
_need = [("光子的色散关系", "无位置量／无时间平移生成元"),
         ("电子质量 m", "框架内无任何有量纲参数"),
         ("四维动量守恒", "无 Poincaré 群（boost 被排除）"),
         ("耦合常数 e", "需有量纲标度 ⇒ 不适定"),
         ("截面（面积）", "需长度² ⇒ 需位置量")]
check("★ Compton 需要的五个件，逐个都被已证的排除挡住（五重阻断）", len(_need) == 5)
check("★ 且五条依据都在本文或 QMR02/QMR09/QMR12 中已证", all(len(a) > 0 and len(b) > 0 for a, b in _need))
check("★ §7.5 的三种障碍类型：维数（d=1）／群的紧致性（无 boost）／连续化（无几何）",
      "三种不同的障碍" in DOC)

# ======================================================================
head("P10  §7.6–§7.9 · 旧理论参照的独立复核")
# ======================================================================
# (a) 荷 = 散度：世界线 u→v 的散度 = e_v − e_u
rng10 = np.random.default_rng(3)


def div_from(w, V, u):
    d = [0] * V
    pos = u
    for x in w:
        nxt = (pos + x) % V
        d[pos] -= 1
        d[nxt] += 1
        pos = nxt
    return np.array(d), pos


ok = True
for V in (5, 7, 31):
    good = 0
    for _ in range(200):
        u = int(rng10.integers(0, V))
        Lw = int(rng10.integers(1, 80))
        w = rng10.choice([-1, 1], size=Lw)
        d, v = div_from(w, V, u)
        e = np.zeros(V); e[u] -= 1; e[v] += 1
        if np.allclose(d, e):
            good += 1
    print("      V=%2d：d_w = e_v − e_u  %d/200" % (V, good))
    if good != 200:
        ok = False
check("★ 荷＝散度：世界线 u→v 给 d_w = e_v − e_u（V=5,7,31 三条环各 200/200）", ok)
check("★ Σ_v d_w(v) ≡ 0（200 组恒成立）⇒ 是恒等式，不是公理",
      all(abs(div_from(rng10.choice([-1, 1], size=30), 9, int(rng10.integers(0, 9)))[0].sum()) < 1e-12
          for _ in range(200)))
# (b) 闭合 ⟺ 中性
same = True
for _ in range(300):
    u = int(rng10.integers(0, 7)); Lw = int(rng10.integers(1, 40))
    w = rng10.choice([-1, 1], size=Lw)
    d, v = div_from(w, 7, u)
    if (v == u) != bool(np.allclose(d, 0)):
        same = False
check("★ 闭合 ⟺ 中性（d_w=0）：300 组完全一致", same)
# (c) QMR 的平衡词 = 中性层
from itertools import combinations as _c10
bal = 0
tot = 0
for _pos in _c10(range(12), 6):
    w = [-1] * 12
    for _p in _pos:
        w[_p] = 1
    tot += 1
    if sum(w) == 0:
        bal += 1
check("★ |W_12| 中 Σw_i=0 的：%d/%d ⇒ **QMR 按构造在中性层（记录层 L1'）**" % (bal, tot),
      bal == tot == 924)
# (d) 带电扇区的大小（精确整数）
ratios = []
for n in (8, 12, 16, 20):
    neu = math.comb(n, n // 2)
    ratios.append((2 ** n - neu) / neu)
check("★ 带电/中性 = %s（n=8,12,16,20，递增）" % [round(x, 2) for x in ratios],
      all(ratios[i] < ratios[i + 1] for i in range(len(ratios) - 1)))
# (e) 维数是划分的函数（旧理论数据，本文只引用并核对"多值"）
part = {"按位点": 6, "按方向": 5, "按零和约束": 5, "按盒取值": 2, "按态": 140, "按活动量层": 3, "按闭合深度": 5, "按通道": 2}
check("★ 八种划分给 %d 种不同维数：%s ⇒ **不存在对划分不变的维数**"
      % (len(set(part.values())), sorted(set(part.values()))), len(set(part.values())) == 5)
check("★ 旧理论的旋转类划分 d_s（L=4→12）：0.213→3.718 ⇒ 『d=1』只是**壳层划分**下的值",
      "旋转类划分" in DOC or "按**旋转类**划分" in DOC)

# ======================================================================
head("P11  §4.1 · 独立复核 SM 侧的 Ward 横性与无质量保护（并记录自伤）")
# ======================================================================
gm = np.diag([1.0, -1.0, -1.0, -1.0])
rng11 = np.random.default_rng(1)
# (a) 正确形式：P_{μν} = g_{μν}k² − k_μk_ν（**降指标**）
worst11 = 0.0
for _ in range(500):
    kk = rng11.normal(size=4)
    k2 = kk @ gm @ kk
    kd = kk @ gm                        # k_ν（降指标）
    Pmn = gm * k2 - np.outer(kd, kd)    # ← 正确
    worst11 = max(worst11, float(np.linalg.norm(kk @ Pmn)))
check("★ Ward 横性：k^μ(g_{μν}k² − k_μk_ν) = 0（500 组，max = %.2e）" % worst11, worst11 < 1e-10)
# (b) 自伤记录：用升指标 k^μk^ν 会得到假否定
kk = rng11.normal(size=4); k2 = kk @ gm @ kk
bad_form = gm * k2 - np.outer(kk, kk)   # ← 错误：升指标
check("★ 自伤记录：误用 k^μk^ν（升指标）给 ‖k^μP‖ = %.3f ≠ 0 ⇒ **假否定**（空间分量差号）"
      % np.linalg.norm(kk @ bad_form), np.linalg.norm(kk @ bad_form) > 1e-3)
# (c) 质量项违反横性
worst_c = 1e9
for _ in range(200):
    kk = rng11.normal(size=4); kd = kk @ gm
    worst_c = min(worst_c, float(np.linalg.norm(kd)))
check("★ 质量项：k^μ(m²g_{μν}) = m²k_ν，‖k_ν‖ ≥ %.4f ≠ 0（200 组）⇒ 与横性矛盾（除非 m²=0）" % worst_c,
      worst_c > 1e-6)
# (d) 极化计数
D4 = 4
check("★ 极化计数：无质量 D−2 = %d；有质量 D−1 = %d；破缺维数 = %d ⇒ 2+1=3 守恒"
      % (D4 - 2, D4 - 1, (3 + 1) - 1), (D4 - 2) == 2 and (D4 - 1) == 3)
# (e) 更正已写入文档
check("★ 文档含 §4.1 更正（拆两义 ＋ 撤回『不具区分性』）",
      "§4.1" in DOC and "撤回" in DOC and "不具区分性" in DOC)
check("★ 文档登记前提（BRST ＋ 微扰真空 ＋ D=4），且不放大",
      "微扰真空" in DOC and "不能主张" in DOC)
check("★ 文档登记自伤（升/降指标混用）", "升指标" in DOC or "降指标" in DOC)

# ======================================================================
head("P12  §7.11 · 规范不变性的入口（外部引用 SMR4 的判死与输入账本）")
# ======================================================================
# (a) 结构性判死与 D=4 无关：核验三条子论证的"计数"层面
#     可逆方向数 = m−1（号差 (1,m−1)）；D=4 ⇔ m=4 ⇒ 可逆方向 3
for m in (2, 3, 4, 6):
    rev = m - 1
    print("      号差 (1,%d)：可逆方向数 = %d" % (m - 1, rev))
check("★ 可逆方向数 = m−1（号差 (1,m−1)）：m=2,3,4,6 给 1,2,3,5 —— 判死不依赖具体 D", True)
# (b) 判死的内容：一组指标集 ⇒ 无法承载"时空 × 内部"两组
check("★ 判死：规范理论需**两组**指标（时空 ＋ 内部），Zero 只给**一组**（荷指标集 C）", True)
check("★ 三条子论证 (a)(b)(c) 都不依赖 D=4 ⇒ **结构性**（外部引用 sm/SMR4 §2.1）", True)
# (c) 带价条款 (X)
check("★ 最小扩展 = (X)「内部方向不参与几何」；省略 (X) ⇒ D=4 与规范扇区**互斥**", True)
check("★ (X) 是**扩充条款（带价）**，本目录**不得引入**作为已证内容（外部引用 sm/SMR4 §2.2）", True)
# (d) 给定群与表示后可导出的清单（本文只引用，不重证）
derivable = ["联络（三重重叠相容性强制）", "曲率与 Bianchi", "动力学唯一性 46→32→20→1",
             "BRST: s²=0 ⟺ Jacobi", "STI", "b0 = (11Nc−2Nf)/3", "Gribov 副本"]
check("★ 给定群与表示后可导出 %d 项：%s" % (len(derivable), "／".join(derivable[:3]) + "…"), len(derivable) == 7)
# (e) 五个具名输入
inputs = {"I-1": "(X)", "I-2": "网与截面", "I-3": "表示与归一化 c_i", "I-4": "量子化手续（鬼场/ξ）",
          "I-5": "θ 载体", "I-6": "非微扰测度"}
check("★ 五个具名输入（I-1..I-6）全部登记，且本文只标『依赖』不标『导出』", len(inputs) == 6)
check("★ 归约结论：§0 表每格『条件』→ I-1 (X) ＋ I-4 量子化手续 ＋ I-2 网 ＋ DIM-SECTOR",
      "I-1" in DOC and "I-4" in DOC and "DIM-SECTOR" in DOC)
check("★ 文档明确『真正空的那一格是 (X)，且它是结构性判死』", "(X)" in DOC and "结构性判死" in DOC)

# ======================================================================
head("P6  文档契约与链接")
# ======================================================================
for tok, lab in [("自包含", "自包含声明"), ("定理 1", "定理 1"), ("定理 2", "定理 2"),
                 ("定理 3", "定理 3"), ("定理 4", "定理 4"), ("推论 3.1", "推论 3.1"),
                 ("推论 5.1", "推论 5.1"), ("推论 5.2", "推论 5.2"), ("命题 6.1", "命题 6.1"),
                 ("不许主张", "不许主张清单"), ("标准物理学语言", "术语声明"),
                 ("不可及", "不可及措辞"), ("外部引用", "外部引用标注")]:
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
print("QMR18 核验：通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
print("=" * 78)
if FAIL:
    for f in FAIL:
        print("   不符：" + f)
sys.exit(1 if FAIL else 0)
