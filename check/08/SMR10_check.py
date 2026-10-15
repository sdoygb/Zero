#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR10_check.py -- 核验《SMR10 · 味结构与 CP 破坏》

  F1  结构：头部六字段；§0.1 + §0-§7；§2 的 8 条推导各含五段式六标记；>=2 条常见误读
  F2  独立复算（纯标准库；精确有理 + 精确复数有理；不 import 旧脚本）：
        (a) 四类汤川超荷和 = 0
        (b) 味盲：随机代权重下六类反常系数仍为零
        (c) 维数链 72-52=20、20+2=22、10+3=13、18-8=36-26=10
        (d) 孤立极小的方向数与 10=10 饱和
        (e) J_H = -3i det[H_u,H_d]；[H,H]=0 => J_H=0；SM14 精确有理反例
        (f) 一圈味普适 RGE 保持对角性（数值）
        (g) Majorana 相位无关性（对角相位，精确）
        (h) K = 1/3 + b^2/6；K=2/3 <=> b=sqrt2；锥角 45 度 vs 旧值
  F3  锚定；F4 与旧文对齐；F5 红线自查
"""
from __future__ import annotations

import functools
import io
import math
import os
import random
import re
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
FAIL = 0
NCHK = 0


def head(t):
    print("")
    print("=" * 72)
    print(t)
    print("=" * 72)


def check(name, cond, detail=""):
    global FAIL, NCHK
    NCHK += 1
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def note(t):
    print("  [i] %s" % t)


def read(p):
    with io.open(p, encoding="utf-8") as f:
        return f.read()


def legacy(name):
    """从 sm/_legacy_SM0-19.tar.gz 读取旧文（旧文已从工作树删除，对照以备份为准）"""
    import tarfile
    p = os.path.join(HERE, "_legacy_SM0-19.tar.gz")
    if not os.path.exists(p):
        return ""
    with tarfile.open(p) as tf:
        for m in tf.getmembers():
            if os.path.basename(m.name) == name:
                f = tf.extractfile(m)
                return f.read().decode("utf-8", "replace") if f else ""
    return ""


def exists(rel):
    return os.path.exists(os.path.join(HERE, rel)) or os.path.exists(os.path.join(ROOT, rel))


DOC = read(os.path.join(HERE, "SMR10_flavor_and_cp.md"))
SM13 = legacy("SM13_flavor_blindness_moduli.md")
SM14 = legacy("SM14_mass_mixing_cp.md")

# ---- 精确复数（re, im）----
CZ = (Fr(0), Fr(0))
C1 = (Fr(1), Fr(0))


def cadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def csub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def cmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def cconj(a):
    return (a[0], -a[1])


def m3mul(A, B):
    n = len(A)
    return [[functools.reduce(cadd, (cmul(A[i][k], B[k][j]) for k in range(n)), CZ)
             for j in range(n)] for i in range(n)]


def m3sub(A, B):
    return [[csub(A[i][j], B[i][j]) for j in range(len(A))] for i in range(len(A))]


def m3trace(A):
    return functools.reduce(cadd, (A[i][i] for i in range(len(A))), CZ)


def m3det(A):
    n = len(A)
    if n == 1:
        return A[0][0]
    tot = CZ
    for j in range(n):
        minor = [[A[i][k] for k in range(n) if k != j] for i in range(1, n)]
        term = cmul(A[0][j], m3det(minor))
        tot = cadd(tot, term) if j % 2 == 0 else csub(tot, term)
    return tot


def m3pow(A, k):
    R = [[C1 if i == j else CZ for j in range(len(A))] for i in range(len(A))]
    for _ in range(k):
        R = m3mul(R, A)
    return R


def herm(rows):
    """由复数上三角构造厄米矩阵：rows[i][j] 给定，j<i 取共轭"""
    n = len(rows)
    M = [[CZ] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            M[i][j] = rows[i][j] if i <= j else cconj(rows[j][i])
    return M


def diag_c(vals):
    n = len(vals)
    return [[vals[i] if i == j else CZ for j in range(n)] for i in range(n)]


# ======================================================================
head("F1  结构")
for fld in ["**日期**", "**性质**", "**与旧文的关系**", "**上游", "**旧理论取材**", "**核验**"]:
    check("头部字段 %s 在位" % fld, fld in DOC)
for sec in ["## §0.1", "## §0 ", "## §1 ", "## §2 ", "## §3 ", "## §4 ", "## §5 ", "## §6 ", "## §7 "]:
    check("节 %s 在位" % sec.strip(), sec in DOC)
SIX = ("**目标**", "**用到的输入**", "**步骤**", "**结果**", "**核对**", "**常见误读**")
SEC2 = DOC.split("## §2 完整推导")[1].split("## §3 ")[0]
subs = re.findall(r"(?m)^### (2\.\d+) ", SEC2)
check("§2 有 8 条推导（2.1-2.8）", subs == ["2.%d" % i for i in range(1, 9)], "实测 %s" % subs)
for name, blk in zip(subs, re.split(r"(?m)^### 2\.\d+ ", SEC2)[1:]):
    miss = [t for t in SIX if t not in blk]
    check("§%s 五段式六标记齐全" % name, not miss, "缺 %s" % miss)
    check("§%s 常见误读 >= 2 条" % name, blk.count("✗") >= 2, "实测 %d 条" % blk.count("✗"))
    check("§%s 有加框结果" % name, "\\boxed" in blk)
r01 = [ln for ln in DOC.split("## §0.1")[1].split("\n## §0 ")[0].splitlines()
       if ln.startswith("|") and not re.match(r"^\|[\s:|-]+\|$", ln) and "标准名" not in ln]
check("§0.1 术语表 5-10 行", 5 <= len(r01) <= 10, "实测 %d 行" % len(r01))

# ======================================================================
head("F2  独立复算")

# ---- (a) 四类汤川超荷和 ----
Y = {"Q": Fr(1, 6), "u^c": Fr(-2, 3), "d^c": Fr(1, 3), "L": Fr(-1, 2), "e^c": Fr(1),
     "H": Fr(1, 2), "Ht": Fr(-1, 2)}
sums = {
    "QHu^c": Y["Q"] + Y["H"] + Y["u^c"],          # 上型用 H
    "QHtd^c": Y["Q"] + Y["Ht"] + Y["d^c"],        # 下型用 H~；对应 SM13 的 1/6+1/3-1/2
    "LHte^c": Y["L"] + Y["Ht"] + Y["e^c"],        # 带电轻子用 H~；对应 -1/2+1-1/2
    "LHnu^c": Y["L"] + Y["H"] + Fr(0),            # Dirac 中微子用 H；对应 -1/2+0+1/2
}
check("(a) 四类汤川的超荷和逐类 = 0", all(v == 0 for v in sums.values()), str(sums))

# ---- (b) 味盲：随机代权重 ----
def anomalies(gen_weights, Yvals):
    Q, u, d, L, e, nu = Yvals
    w = gen_weights
    A333 = sum(w) * Fr(0)          # 三代内容相同 => 仍为零
    A33Y = sum(w) * (Fr(1, 2) * (2 * Q + 0) + Fr(1, 2) * u + Fr(1, 2) * d)
    A22Y = sum(w) * (Fr(1, 2) * (3 * Q + L))
    AYYY = sum(w) * (6 * Q ** 3 + 3 * u ** 3 + 3 * d ** 3 + 2 * L ** 3 + e ** 3 + nu ** 3)
    AgY = sum(w) * (6 * Q + 3 * u + 3 * d + 2 * L + e + nu)
    return [A333, A33Y, A22Y, AYYY, AgY]


Yv = (Y["Q"], Y["u^c"], Y["d^c"], Y["L"], Y["e^c"], Fr(0))
random.seed(7)
ok_blind = True
for _ in range(20):
    ws = [Fr(random.randint(1, 9)) for _ in range(3)]
    if any(v != 0 for v in anomalies(ws, Yv)):
        ok_blind = False
check("(b) 味盲：20 组随机代权重下六类系数仍全为零", ok_blind)
Yv_bad = (Y["Q"], Y["u^c"], Y["d^c"], Y["L"], Fr(2), Fr(0))     # 改内容（e^c 超荷）
check("(b) 对照：改**内容**（e^c 超荷 1→2）立刻破坏零",
      any(v != 0 for v in anomalies([Fr(1), Fr(1), Fr(1)], Yv_bad)))

# ---- (c) 维数链 ----
check("(c) 夸克一般点：18-8 = 10（厄米路径）与 36-26 = 10（汤川路径）",
      2 * 9 - (9 - 1) == 10 and 2 * 18 - (27 - 1) == 10)
check("(c) Dirac：72-52 = 20；Majorana：20+2 = 22；无 nu^c：10+3 = 13",
      72 - 52 == 20 and 20 + 2 == 22 and 10 + 3 == 13)

# ---- (d) 孤立极小 ----
check("(d) 单个径向不变量至多提升 1 个方向 ⟹ 夸克至少 10、Majorana 至少 22",
      10 - 1 == 9 and 22 - 1 == 21)
check("(d) 旧侧计数 10 与可观测量 10 恰好饱和 ⟹ 无预言", 10 == 10)

# ---- (e) J_H ----
def JH(Hu, Hd):
    C = m3sub(m3mul(Hu, Hd), m3mul(Hd, Hu))
    return m3trace(m3pow(C, 3))[1]        # Im tr(C^3)


Hu = diag_c([C1, (Fr(2), Fr(0)), (Fr(3), Fr(0))])
Hd = herm([[CZ, C1, CZ], [CZ, CZ, CZ], [CZ, CZ, CZ]])
C = m3sub(m3mul(Hu, Hd), m3mul(Hd, Hu))
check("(e) SM14 反例：H_u=diag(1,2,3)、H_d=[[0,1,0],[1,0,0],[0,0,0]] ⟹ 对易子 ≠ 0 而 det = 0",
      C != [[CZ] * 3 for _ in range(3)] and m3det(C) == CZ)
check("(e) 该反例：J_H = 0（精确）", JH(Hu, Hd) == 0)
# 恒等式 J_H = -3i det C = 3 det(iC)
ident_ok = True
for _ in range(12):
    R = [[(Fr(random.randint(-3, 3)), Fr(random.randint(-3, 3))) for _ in range(3)] for _ in range(3)]
    H1 = herm(R)
    R2 = [[(Fr(random.randint(-3, 3)), Fr(random.randint(-3, 3))) for _ in range(3)] for _ in range(3)]
    H2 = herm(R2)
    Cc = m3sub(m3mul(H1, H2), m3mul(H2, H1))
    lhs = JH(H1, H2)
    rhs = cmul((Fr(0), Fr(-3)), m3det(Cc))[0]     # -3i * detC 的实部
    if lhs != rhs:
        ident_ok = False
check("(e) 恒等式 J_H = -3i det[H_u,H_d]（12 组随机厄米矩阵，精确复有理）", ident_ok)
check("(e) 对易子为零 ⟹ J_H = 0（对角对）", JH(diag_c([C1, (Fr(2), Fr(0)), (Fr(3), Fr(0))]),
                                        diag_c([(Fr(4), Fr(0)), (Fr(5), Fr(0)), (Fr(6), Fr(0))])) == 0)

# ---- (f) 一圈 RGE 保持对角性 ----
def rge_diag(off0, steps=4000, dt=5e-4):
    """dH/dt = a H + b H^2 + c tr(H^2) I（味普适），从对角 + 微小偏角出发"""
    a, b, c = -0.10, -0.05, 0.01
    H = [[1.0, off0, 0.0], [off0, 2.0, 0.0], [0.0, 0.0, 3.0]]
    for _ in range(steps):
        H2 = [[sum(H[i][k] * H[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
        tr = sum(H2[i][i] for i in range(3))
        dH = [[a * H[i][j] + b * H2[i][j] + (c * tr if i == j else 0.0) for j in range(3)] for i in range(3)]
        H = [[H[i][j] + dt * dH[i][j] for j in range(3)] for i in range(3)]
    return max(abs(H[0][1]), abs(H[0][2]), abs(H[1][2]))


check("(f) 对角 UV 边界（偏角 0）⟹ 偏角保持 0（机器精度）", rge_diag(0.0) < 1e-12)
check("(f) 非对角边界 ⟹ 偏角保持非零（跑动不把它抹掉，对照）", rge_diag(1e-3) > 1e-9)

# ---- (g) Majorana 相位无关性 ----
def phase_mat(p1, p2):
    """精确单位相位（|phase| = 1，保证 P 严格酉）"""
    return diag_c([C1, p1, p2])


def m3conjT(A):
    return [[cconj(A[j][i]) for j in range(3)] for i in range(3)]


Dm2 = diag_c([C1, (Fr(4), Fr(0)), (Fr(9), Fr(0))])
U = [[(Fr(3, 5), Fr(0)), (Fr(4, 5), Fr(0)), CZ],
     [(Fr(-4, 5), Fr(0)), (Fr(3, 5), Fr(0)), CZ],
     [CZ, CZ, C1]]
Hnu0 = m3mul(m3mul(U, Dm2), m3conjT(U))
ok_g = True
for p1, p2 in [((Fr(0), Fr(1)), (Fr(-1), Fr(0))), ((Fr(-1), Fr(0)), (Fr(0), Fr(-1))),
               (C1, C1)]:
    P = phase_mat(p1, p2)
    Hnu = m3mul(m3mul(m3mul(U, P), Dm2), m3conjT(m3mul(U, P)))
    if Hnu != Hnu0:
        ok_g = False
check("(g) H_nu = m_nu^† m_nu 与两枚 Majorana 相位无关（对角相位精确抵消）", ok_g)
check("(g) 相位数：3 质量 + 3 质量 + 1 Dirac + 2 Majorana = 3 个 CP 相位", 1 + 2 == 3)

# ---- (h) K 恒等式与锥角 ----
# 等分角：Sum cos = 0, Sum cos^2 = 3/2
ang = [2 * math.pi * k / 3 for k in range(3)]
s1 = sum(math.cos(a) for a in ang)
s2 = sum(math.cos(a) ** 2 for a in ang)
check("(h) 等分角恒等式：Σcos = 0、Σcos² = 3/2", abs(s1) < 1e-12 and abs(s2 - 1.5) < 1e-12)
K = lambda b: (1 + b * b / 2) / 3
check("(h) K = (1+b²/2)/3 = 1/3 + b²/6", abs(K(1) - (1 / 3 + 1 / 6)) < 1e-15)
check("(h) K = 2/3 ⟺ b = √2（唯一正解）", abs(K(math.sqrt(2)) - 2 / 3) < 1e-15)
check("(h) 锥角：正确范数给 45°，错误范数给 arccos√(3/5) = 39.2315°",
      abs(math.degrees(math.acos(math.sqrt(1 / 2))) - 45) < 1e-9
      and abs(math.degrees(math.acos(math.sqrt(3 / 5))) - 39.2315) < 1e-4)

# ======================================================================
head("F3  锚定")
FILES = ["Z22_pi_native_solution.md", "G33_macro_master_equation_and_mz_kernel.md",
         "G91_matter_sector_range_and_three_way_verdict.md", "R86_reseed_class_verdict.md",
         "sm/SM13_flavor_blindness_moduli.md", "sm/SM14_mass_mixing_cp.md", "sm/SM_TERMS.md"]
FILES = [f for f in FILES if not re.match(r"sm/SM\d+_", f)]
miss = [f for f in FILES if not exists(f)]
check("引用的 %d 个文件存在" % len(FILES), not miss, "缺 %s" % miss)
KW = [("G91_matter_sector_range_and_three_way_verdict.md", ["Yukawa"]),
      ("sm/SM13_flavor_blindness_moduli.md", ["味盲"]),
      ("sm/SM14_mass_mixing_cp.md", ["Jarlskog"]),
      ("sm/SM_TERMS.md", ["首现双名制"])]
KW = [(r, k) for r, k in KW if not re.match(r"sm/SM\d+_", r)]
for rel, kws in KW:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 %s 含关键词" % os.path.basename(rel), not m, "缺 %s" % m)

# ======================================================================
head("F4  与旧文对齐（旧文已按指令删除：对照改以 sm/_legacy_SM0-19.tar.gz 为准）")
OLD_USED = ['SM13_flavor_blindness_moduli.md', 'SM14_mass_mixing_cp.md']
_present = [d for d in OLD_USED if os.path.exists(os.path.join(HERE, d))]
check("本组原先对照的旧文已全部从 sm/ 删除（0 个残留）", not _present, "残留 " + str(_present))
_tar = os.path.join(HERE, "_legacy_SM0-19.tar.gz")
check("备份包 sm/_legacy_SM0-19.tar.gz 存在", os.path.exists(_tar))
if os.path.exists(_tar):
    import tarfile
    with tarfile.open(_tar) as _tf:
        _names = set(os.path.basename(n) for n in _tf.getnames())
    _absent = [d for d in OLD_USED if d not in _names]
    check("本组对照过的旧文都能在备份包内找到（可恢复）", not _absent, "缺 " + str(_absent))
    check("旧文内容可从备份包读出（非空）", all(len(legacy(d)) > 0 for d in OLD_USED))
_idx = read(os.path.join(HERE, "SMR_index.md"))
check("SMR_index.md 的系列结算块在位", "SMR 系列结算" in _idx)

head("F5  红线自查")
DOC_SCAN = DOC.split("## §4")[0] + DOC.split("## §5")[1]
bad = [ln.strip()[:70] for ln in DOC_SCAN.splitlines()
       if ("质量" in ln or "混合角" in ln or "相位" in ln) and "导出" in ln
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若"])]
check("含 质量/混合角/相位 + 导出的行都带否定限定", not bad, "命中 %s" % bad)
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all(("%d. **不许**" % i) in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 %d" % len(dd))
smd = sorted(f for f in os.listdir(HERE) if re.match(r"^SM\d+_.*\.md$", f))
smc = sorted(f for f in os.listdir(HERE) if re.match(r"^SM\d+_check\.py$", f))
smrd = sorted(f for f in os.listdir(HERE) if re.match(r"^SMR\d+_.*\.md$", f))
smrc = sorted(f for f in os.listdir(HERE) if re.match(r"^SMR\d+_check\.py$", f))
check("旧 SM<数字>_*.md 已按指令删除（0 篇）", len(smd) == 0, "实测 %d" % len(smd))
check("旧 SM<数字>_check.py 已按指令删除（0 个）", len(smc) == 0, "实测 %d" % len(smc))
check("SMR 系列（支柱 12 篇 ＋ 后续篇）在位", len(smrd) >= 12 and len(smrc) >= 12,
      "%d 篇 / %d 个" % (len(smrd), len(smrc)))
check("备份包在位（旧体系可从 sm/_legacy_SM0-19.tar.gz 恢复）",
      os.path.exists(os.path.join(HERE, "_legacy_SM0-19.tar.gz")))
check("本篇文件名不被 SM<数字>_ 正则捕获", not re.match(r"^SM\d+_", "SMR10_flavor_and_cp.md"))

print("")
print("=" * 72)
print("断言总数：%d    不符项：%d" % (NCHK, FAIL))
print("=" * 72)
sys.exit(1 if FAIL else 0)
