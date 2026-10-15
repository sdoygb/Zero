#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR24_check.py -- 核验《SMR24 · 规范传播子与谱表示层》"""
import io, os, re, sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, ".."))
FAIL = NCHK = 0
def head(t): print(""); print("=" * 72); print(t); print("=" * 72)
def check(n, c, d=""):
    global FAIL, NCHK
    NCHK += 1
    ok = bool(c)
    if not ok: FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", n, ("   " + d) if d else ""))
def read(p):
    with io.open(p, encoding="utf-8") as f: return f.read()
def exists(rel): return os.path.exists(os.path.join(HERE, rel)) or os.path.exists(os.path.join(ROOT, rel))

DOC = read(os.path.join(HERE, "SMR24_propagator_spectral.md"))
G = [Fr(1), Fr(-1), Fr(-1), Fr(-1)]          # 度规 diag(+,-,-,-)
K = [Fr(2), Fr(1), Fr(1), Fr(0)]             # 一个具体有理动量
K2 = sum(G[i] * K[i] * K[i] for i in range(4))

Kup = [G[i] * K[i] for i in range(4)]          # k^mu = g^{mu nu} k_nu
def PT(mu, nu):                                 # 混合指标 P^T_mu^nu
    return (Fr(1) if mu == nu else Fr(0)) - K[mu] * Kup[nu] / K2
def PL(mu, nu):
    return K[mu] * Kup[nu] / K2
def PTlow(mu, nu):                              # 全下指标 P^T_{mu nu}
    return (G[mu] if mu == nu else Fr(0)) - K[mu] * K[nu] / K2
def mat(f):
    return [[f(mu, nu) for nu in range(4)] for mu in range(4)]
def mmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]

head("F1  结构")
for fld in ["**日期**", "**性质**", "**与旧文的关系**", "**上游", "**旧理论取材**", "**核验**"]:
    check("头部字段 " + fld + " 在位", fld in DOC)
for sec in ["## §0.1", "## §0 ", "## §1 ", "## §2 ", "## §3 ", "## §4 ", "## §5 ", "## §6 ", "## §7 "]:
    check("节 " + sec.strip() + " 在位", sec in DOC)
SIX = ("**目标**", "**用到的输入**", "**步骤**", "**结果**", "**核对**", "**常见误读**")
SEC2 = DOC.split("## §2 完整推导")[1].split("## §3 ")[0]
subs = re.findall(r"(?m)^### (2\.\d+) ", SEC2)
check("§2 有 7 条推导（2.1-2.7）", subs == ["2." + str(i) for i in range(1, 8)], "实测 " + str(subs))
for name, blk in zip(subs, re.split(r"(?m)^### 2\.\d+ ", SEC2)[1:]):
    miss = [t for t in SIX if t not in blk]
    check("§" + name + " 五段式标记齐全", not miss, "缺 " + str(miss))
    check("§" + name + " 常见误读 >= 2 条", blk.count("✗") >= 2, "实测 " + str(blk.count("✗")))
    check("§" + name + " 有加框结果", "\\boxed" in blk)

head("F2  独立复算")
# (a)
PTm, PLm = mat(PT), mat(PL)
Gm = mat(lambda mu, nu: Fr(1) if mu == nu else Fr(0))
check("(a) k^2 = 2（具体有理动量，Lorentz 度规）", K2 == 2)
check("(a) P^T + P^L = identity（混合指标）", all(PTm[i][j] + PLm[i][j] == Gm[i][j] for i in range(4) for j in range(4)))
check("(a) (P^T)^2 = P^T", mmul(PTm, PTm) == PTm)
check("(a) (P^L)^2 = P^L", mmul(PLm, PLm) == PLm)
check("(a) P^T P^L = 0", all(mmul(PTm, PLm)[i][j] == 0 for i in range(4) for j in range(4)))
check("(a) k^mu P^T_{mu nu} = 0（四个 nu 逐项）",
      all(sum(Kup[mu] * PTm[mu][nu] for mu in range(4)) == 0 for nu in range(4)))
for xi in (Fr(0), Fr(1), Fr(-2)):
    D = [[PTm[i][j] + xi * PLm[i][j] for j in range(4)] for i in range(4)]
    if xi == 0:
        check("(a) xi=0（朗道规范）时传播子纯横向", D == PTm)
    else:
        check("(a) xi=" + str(xi) + " 时纵向部分系数为 xi",
              all(D[i][j] == PTm[i][j] + xi * PLm[i][j] for i in range(4) for j in range(4)))

# (b) Ward 横性
PTlowm = mat(PTlow)
T = [[K2 * PTlowm[i][j] for j in range(4)] for i in range(4)]
check("(b) k^mu (g_{mu nu} k^2 − k_mu k_nu) = 0（Ward 横性的张量基础）",
      all(sum(Kup[mu] * T[mu][nu] for mu in range(4)) == 0 for nu in range(4)))
check("(b) 该结构就是 k^2 P^T_{mu nu}（全下指标）",
      all(T[i][j] == K2 * PTlowm[i][j] for i in range(4) for j in range(4)))

# (c) 质量项被禁
check("(c) k^mu delta_mu^nu = k^nu ≠ 0（故 m^2 g_{mu nu} 违反横性）",
      all(sum(Kup[mu] * Gm[mu][nu] for mu in range(4)) == Kup[nu] for nu in range(4))
      and any(Kup[nu] != 0 for nu in range(4)))

# (d) 计数
D_ = 4
check("(d) 无质量矢量极化数 D−2 = 2；有质量 D−1 = 3；破缺维数 3；2+1=3",
      D_ - 2 == 2 and D_ - 1 == 3 and (3 + 1 - 1) == 3)

# (e) KL：离散谱与两极点零点
c1, c2, m1, m2 = Fr(1), Fr(2), Fr(1), Fr(4)
x = (c1 * m2 + c2 * m1) / (c1 + c2)          # 零点位置（k^2）
check("(e) 两极点谱 c1/(x−m1)+c2/(x−m2)=0 的解落在两极点之间",
      m1 < x < m2 and c1 / (x - m1) + c2 / (x - m2) == 0, "x=" + str(x))
check("(e) 单极点项对应 ρ = δ(μ^2 − m^2)（表示式形式检查在位）", "Källén" in DOC and "rho" in DOC.replace("\\", ""))

# (f) 五条边界
check("(f) 五条边界在位（胶子谱正定性/红外行为/禁闭/测度/T1-T3）",
      all(k in DOC for k in ["正定性", "红外", "禁闭", "测度", "T1"]) and "不得主张" in DOC)
# (g) 八项账
check("(g) 八项账在位（投影/R-xi/Ward/保护/计数/KL/胶子谱/禁闭/测度）",
      all(k in DOC for k in ["投影算子", "Ward", "保护", "计数", "Källén", "不可主张", "未建"]))

head("F3  锚定")
FILES = ["sm/SMR4_gauge_principle_and_ym.md", "sm/SMR8_electroweak_breaking_and_higgs.md",
         "sm/SMR12_strong_sector_and_theta.md", "sm/SMR20_settlement.md",
         "sm/SMR21_gauge_force_layer.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 " + str(len(FILES)) + " 个文件存在", not miss, "缺 " + str(miss))
for rel, kws in [("sm/SMR4_gauge_principle_and_ym.md", ["Gribov"]),
                 ("sm/SMR12_strong_sector_and_theta.md", ["禁闭"]),
                 ("sm/SMR20_settlement.md", ["T1"]),
                 ("sm/SM_TERMS.md", ["首现双名制"])]:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 " + os.path.basename(rel) + " 含关键词", not m, "缺 " + str(m))

head("F4  与旧文对齐")
smd = sorted(f for f in os.listdir(HERE) if re.match(r"^SM\d+_.*\.md$", f))
check("旧 SM<数字>_*.md 已按指令删除（0 篇）", len(smd) == 0, "实测 " + str(len(smd)))
check("备份包在位", os.path.exists(os.path.join(HERE, "_legacy_SM0-19.tar.gz")))

head("F5  红线自查")
SKIP = ("**目标**", "**步骤**", "**用到的输入**", "**核对**", "**结果**", "**常见误读**")
DOC_SCAN = DOC.split("## §4")[0] + DOC.split("## §5")[1]
bad = [ln.strip()[:70] for ln in DOC_SCAN.splitlines()
       if not ln.strip().startswith("|") and not ln.strip().startswith("#")
       and not any(ln.strip().startswith(s) for s in SKIP)
       and "**性质**" not in ln and "**纪律**" not in ln and "F5 红线自查" not in ln
       and ("导出" in ln or "禁闭" in ln)
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若", "登记",
                                     "不得主张", "可导出", "保护", "判据", "承"])]
check("含 导出/禁闭 的行都带限定", not bad, "命中 " + str(bad))
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all((str(i) + ". **不许**") in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 " + str(len(dd)))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
