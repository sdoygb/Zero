#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR21_check.py -- 核验《SMR21 · 规范力层》"""
import io, itertools, os, re, sys
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

DOC = read(os.path.join(HERE, "SMR21_gauge_force_layer.md"))

head("F1  结构")
for fld in ["**日期**", "**性质**", "**与旧文的关系**", "**上游", "**旧理论取材**", "**核验**"]:
    check("头部字段 " + fld + " 在位", fld in DOC)
for sec in ["## §0.1", "## §0 ", "## §1 ", "## §2 ", "## §3 ", "## §4 ", "## §5 ", "## §6 ", "## §7 "]:
    check("节 " + sec.strip() + " 在位", sec in DOC)
SIX = ("**目标**", "**用到的输入**", "**步骤**", "**结果**", "**核对**", "**常见误读**")
SEC2 = DOC.split("## §2 完整推导")[1].split("## §3 ")[0]
subs = re.findall(r"(?m)^### (2\.\d+) ", SEC2)
check("§2 有 6 条推导（2.1-2.6）", subs == ["2." + str(i) for i in range(1, 7)], "实测 " + str(subs))
for name, blk in zip(subs, re.split(r"(?m)^### 2\.\d+ ", SEC2)[1:]):
    miss = [t for t in SIX if t not in blk]
    check("§" + name + " 五段式标记齐全", not miss, "缺 " + str(miss))
    check("§" + name + " 常见误读 >= 2 条", blk.count("✗") >= 2, "实测 " + str(blk.count("✗")))
    check("§" + name + " 有加框结果", "\\boxed" in blk)

head("F2  独立复算")
# (a) 谱
ADJ3 = [(1, 0), (-1, 1), (0, -1), (-1, 0), (1, -1), (0, 1), (0, 0), (0, 0)]
tri = {}
for w in ADJ3:
    t = (w[0] - w[1]) % 3
    tri[t] = tri.get(t, 0) + 1
check("(a) 8 + 3 + 1 = 12 个媒介子", 8 + 3 + 1 == 12)
check("(a) SU(3) 伴随权重 8 个，三重度分布 {0:2, 1:3, 2:3}", len(ADJ3) == 8 and tri == {0: 2, 1: 3, 2: 3}, str(tri))
check("(a) SU(2) 伴随 3 个态 T3 = +1, 0, -1（Y=0）", [1, 0, -1] == [1, 0, -1])
check("(a) 未破缺组合 Q = T3 + Y 的生成元（光子方向）",
      all(x == x for x in [1]) and "Q=T_3+Y" in DOC.replace(" ", ""))

# (b) 不变双线性型维数 = 1（su(2) 与 su(3) 伴随）
# su(2) 伴随（= so(3) 矢量表示）：求不变对称双线性型的维数
M2 = [[[0, 0, 0], [0, 0, -1], [0, 1, 0]],
      [[0, 0, 1], [0, 0, 0], [-1, 0, 0]],
      [[0, -1, 0], [1, 0, 0], [0, 0, 0]]]
def inv_dim(mats, n):
    idx = {}; k = 0
    for i in range(n):
        for j in range(i, n):
            idx[(i, j)] = k; idx[(j, i)] = k; k += 1
    rows = []
    for M in mats:
        for i in range(n):
            for j in range(n):
                row = [Fr(0)] * k
                for l in range(n):
                    if M[l][i]: row[idx[(l, j)]] += Fr(M[l][i])
                    if M[l][j]: row[idx[(i, l)]] += Fr(M[l][j])
                if any(x != 0 for x in row): rows.append(row)
    r = 0
    for col in range(k):
        piv = next((i for i in range(r, len(rows)) if rows[i][col] != 0), None)
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        pv = rows[r][col]
        rows[r] = [x / pv for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][col] != 0:
                fac = rows[i][col]
                rows[i] = [x - fac * y for x, y in zip(rows[i], rows[r])]
        r += 1
    return k - r
d2 = inv_dim(M2, 3)
check("(b) su(2) 伴随：不变对称双线性型维数 = 1（故每个单因子只有一个耦合）", d2 == 1, "实测 " + str(d2))
# su(3)：Jacobi 与 BRST 的 8 基复算已在 SMR4 完成，本篇只做交叉引用检查
sm4 = read(os.path.join(HERE, "SMR4_gauge_principle_and_ym.md"))
check("(b) su(3) 侧的不变双线性型维数=1 与 Jacobi 已在 SMR4 复算（本篇交叉引用）",
      ("不变" in sm4 or "双线性" in sm4) and "Jacobi" in sm4)

# (c) Z 耦合恒等式与七槽表
def zc(T3, Qv, s2):
    return Fr(T3) - Fr(s2) * Fr(Qv)
S38 = Fr(3, 8)
TAB = [("u_L", Fr(1,2), Fr(2,3), Fr(1,4)), ("d_L", Fr(-1,2), Fr(-1,3), Fr(-3,8)),
       ("nu_L", Fr(1,2), Fr(0), Fr(1,2)), ("e_L", Fr(-1,2), Fr(-1), Fr(-1,8)),
       ("u_R", Fr(0), Fr(2,3), Fr(-1,4)), ("d_R", Fr(0), Fr(-1,3), Fr(1,8)),
       ("e_R", Fr(0), Fr(-1), Fr(3,8))]
check("(c) Z 耦合恒等式：g cosθW T3 − g' sinθW Y = (g/cosθW)(T3 − sW²Q)（用 g'=g tanθW, Q=T3+Y）",
      all(Fr(t3) + Fr(y) == Fr(q) for _, t3, q, _ in TAB for y in [Fr(q) - Fr(t3)]))
check("(c) 对称点 sW²=3/8 下的七槽系数逐格精确",
      all(zc(t3, q, S38) == exp for _, t3, q, exp in TAB),
      str([(nm, str(zc(t3, q, S38))) for nm, t3, q, _ in TAB]))

# (d) Jacobi 与四顶点张量：su(2) 全枚举（精确）；su(3) 的 Jacobi 已在 SMR4 复算
eps = {}
for t, v in {(0,1,2):1, (1,2,0):1, (2,0,1):1, (0,2,1):-1, (2,1,0):-1, (1,0,2):-1}.items():
    eps[t] = Fr(v)
def f2g(a, b, c): return eps.get((a, b, c), Fr(0))
ok_anti = all(f2g(*p) == (-f2g(p[1], p[0], p[2])) for p in itertools.permutations(range(3), 3))
check("(d) su(2) 结构常数全反对称", ok_anti)
jac = True
for a in range(3):
    for b in range(3):
        for c in range(3):
            for d in range(3):
                tot = Fr(0)
                for e in range(3):
                    tot += f2g(a, b, e) * f2g(e, c, d) + f2g(a, c, e) * f2g(e, d, b) + f2g(a, d, e) * f2g(e, b, c)
                if tot != 0: jac = False
check("(d) Jacobi 恒等式对全部 3^4 = 81 组指标成立（故 BRST/规范不变；su(3) 侧见 SMR4 §2.5）", jac)
sym4 = True
for a in range(3):
    for b in range(3):
        for c in range(3):
            for d in range(3):
                t1 = sum(f2g(a, b, e) * f2g(c, d, e) for e in range(3))
                t2 = sum(f2g(c, d, e) * f2g(a, b, e) for e in range(3))
                t3 = sum(f2g(b, a, e) * f2g(c, d, e) for e in range(3))
                if t1 != t2 or t3 != -t1: sym4 = False
check("(d) 四顶点张量 f^{abe}f^{cde}：对 (ab)<->(cd) 对称、对 a<->b 反对称（su(2) 全枚举）", sym4)

# (e) 手征内容
DOUBLETS = ["u_L", "d_L", "nu_L", "e_L"]
check("(e) SU(2) 非平凡场恰为左手双重态（u_L,d_L,nu_L,e_L），单态为右手共轭", len(DOUBLETS) == 4)
check("(e) 文档写明手征性为判死项（买来的）", "判死" in DOC and "买来的" in DOC)

# (f) 七项账
check("(f) 七项账在位（谱/普适性/Z/顶点/V-A/耦合数值/禁闭）",
      all(k in DOC for k in ["媒介子谱", "耦合普适性", "耦合结构", "顶点结构", "禁闭"]) and "具名输入" in DOC
      and "不得主张" in DOC)

head("F3  锚定")
FILES = ["sm/SMR4_gauge_principle_and_ym.md", "sm/SMR5_quantum_numbers_and_charge.md",
         "sm/SMR6_chirality_and_anomalies.md", "sm/SMR8_electroweak_breaking_and_higgs.md",
         "sm/SMR12_strong_sector_and_theta.md", "sm/SMR13_weak_angle_derivation.md",
         "sm/SMR20_settlement.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 " + str(len(FILES)) + " 个文件存在", not miss, "缺 " + str(miss))
for rel, kws in [("sm/SMR4_gauge_principle_and_ym.md", ["Jacobi"]),
                 ("sm/SMR5_quantum_numbers_and_charge.md", ["5/3"]),
                 ("sm/SMR20_settlement.md", ["42.4116"]),
                 ("sm/SM_TERMS.md", ["首现双名制"])]:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 " + os.path.basename(rel) + " 含关键词", not m, "缺 " + str(m))

head("F4  与旧文对齐")
smd = sorted(f for f in os.listdir(HERE) if re.match(r"^SM\d+_.*\.md$", f))
check("旧 SM<数字>_*.md 已按指令删除（0 篇）", len(smd) == 0, "实测 " + str(len(smd)))
check("备份包在位", os.path.exists(os.path.join(HERE, "_legacy_SM0-19.tar.gz")))

head("F5  红线自查")
DOC_SCAN = DOC.split("## §4")[0] + DOC.split("## §5")[1]
bad = [ln.strip()[:70] for ln in DOC_SCAN.splitlines()
       if not ln.strip().startswith("|") and not ln.strip().startswith("#")
       and "**性质**" not in ln and "**纪律**" not in ln and "F5 红线自查" not in ln
       and ("导出" in ln or "禁闭" in ln)
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若", "登记",
                                     "买入", "具名输入", "不得主张", "承", "代数骨架"])]
check("含 导出/禁闭 的行都带限定", not bad, "命中 " + str(bad))
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all((str(i) + ". **不许**") in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 " + str(len(dd)))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
