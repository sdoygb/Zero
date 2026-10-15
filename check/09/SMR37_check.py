#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR37_check.py -- 核验《SMR37 · 味结构能否解释扇区选择性》"""
import cmath, io, math, os, re, sys

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

DOC = read(os.path.join(HERE, "SMR37_flavor_structure_selectivity.md"))
W = cmath.exp(2j * math.pi / 3)

def mm(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)] for i in range(n)]
def dag(A): return [[A[j][i].conjugate() for j in range(len(A))] for i in range(len(A[0]))]
def K_of(ms):
    s1 = sum(math.sqrt(max(v, 0.0)) for v in ms)
    return sum(ms) / (s1 * s1)
def tr(A): return sum(A[i][i] for i in range(len(A)))
def tr2(A): return tr(mm(A, A))
def det3(A):
    return (A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
            - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
            + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]))

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
# (a) 基无关：用相似不变量（trace, tr M^2, det）判定
M0 = [[2.0, 0.3 + 0.2j, 0.1 - 0.4j], [0.3 - 0.2j, 5.0, 0.6j], [0.1 + 0.4j, -0.6j, 9.0]]
th = 0.7
U = [[math.cos(th), -math.sin(th) * cmath.exp(1j * 0.3), 0],
     [math.sin(th) * cmath.exp(-1j * 0.3), math.cos(th), 0],
     [0, 0, cmath.exp(1j * 1.1)]]
M1 = mm(mm(U, M0), dag(U))
inv0 = (tr(M0), tr2(M0), det3(M0))
inv1 = (tr(M1), tr2(M1), det3(M1))
check("(a) 相似不变量 (tr M, tr M^2, det M) 在 U M U^dag 下不变",
      all(abs(x - y) < 1e-9 for x, y in zip(inv0, inv1)),
      "%.6f/%.6f" % (inv0[0].real, inv1[0].real))
check("(a) 故特征多项式不变 ⟹ 本征值不变 ⟹ K 基无关（对任意酉 U，含 CKM/PMNS 型）",
      all(abs(x - y) < 1e-9 for x, y in zip(inv0, inv1)))
ml = (0.51099895, 105.6583755, 1776.86)
check("(a) 轻子 K = 0.666660（自算）", abs(K_of(ml) - 0.66666) < 1e-5, "%.8f" % K_of(ml))
# (b) 循环矩阵：解析本征值 + 行列式检验（无需数值本征求解）
W3 = cmath.exp(2j * math.pi / 3)
def circ(a0, a1):
    P = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    P2 = mm(P, P)
    return [[a0 * (1 if i == j else 0) + a1 * P[i][j] + a1.conjugate() * P2[i][j] for j in range(3)] for i in range(3)]
def det_close(M, lam):
    A = [[M[i][j] - (lam if i == j else 0) for j in range(3)] for i in range(3)]
    return abs(det3(A))
ok_b = True
for a0, a1 in [(1.0, 0.4 * cmath.exp(0.3j)), (2.0, 0.9 * cmath.exp(1.1j)), (1.0, (1 / math.sqrt(2)) * cmath.exp(-0.6j))]:
    S = circ(a0, a1)
    for k in range(3):
        lam = a0 + a1 * W3 ** k + a1.conjugate() * W3 ** (-k)
        if det_close(S, lam) > 1e-9: ok_b = False
        if abs(lam - (a0 + 2 * abs(a1) * math.cos(cmath.phase(a1) + 2 * math.pi * k / 3))) > 1e-12: ok_b = False
check("(b) 厄米循环 S 的本征值 = a0 + 2|a1| cos(theta + 2 pi k/3)（行列式检验，三组参数）", ok_b)
def K_circ(a0, a1):
    lam = [a0 + 2 * abs(a1) * math.cos(cmath.phase(a1) + 2 * math.pi * k / 3) for k in range(3)]
    ev = [x * x for x in lam]
    return K_of(sorted(ev))
ok_id = True
for a0, a1 in [(1.0, 0.3), (2.0, 0.5), (1.0, 1 / math.sqrt(2))]:
    if abs(K_circ(a0, a1) - (1 / 3 + (2 / 3) * (abs(a1) ** 2) / (a0 ** 2))) > 1e-12: ok_id = False
check("(b) 恒等式 K = 1/3 + (2/3)(|a1|^2/a0^2)（三组参数）", ok_id)
check("(b) 比值条件 |a1|^2/a0^2 = 1/2 时 K = 2/3（M = S^2，两块功率相等）",
      abs(K_circ(1.0, 1 / math.sqrt(2)) - 2 / 3) < 1e-12, "%.12f" % K_circ(1.0, 1 / math.sqrt(2)))
check("(b) 反例：|a1|/a0 = 1 给 K = 1（退化谱）", abs(K_circ(1.0, 1.0) - 1.0) < 1e-12)
# (c) 两循环矩阵同时对角化 => V = 1
F = [[W ** (j * k) / math.sqrt(3) for k in range(3)] for j in range(3)]
Su = circ(1.0, 0.3 * cmath.exp(0.2j))
Sd = circ(1.4, 0.5 * cmath.exp(-1.0j))
Du = mm(mm(dag(F), Su), F)
Dd = mm(mm(dag(F), Sd), F)
offu = max(abs(Du[i][j]) for i in range(3) for j in range(3) if i != j)
offd = max(abs(Dd[i][j]) for i in range(3) for j in range(3) if i != j)
check("(c) 两个厄米循环矩阵同时被 DFT 矩阵 F 对角化（非对角元 ~ 0）", offu < 1e-12 and offd < 1e-12,
      "%.2e %.2e" % (offu, offd))
check("(c) 故 U_uL = U_dL = F => V_CKM = F^dag F = 1（与观测矛盾）",
      all(abs(mm(dag(F), F)[i][j] - (1 if i == j else 0)) < 1e-12 for i in range(3) for j in range(3)))
check("(c) 结论：至多一个夸克扇区可有循环结构（文档写明只给上限）", "至多一个" in DOC and "上限" in DOC)
# (d) 三行判定
check("(d) 味混合判死三行在位（K 基无关 / 循环性基依赖 / 上限不给观测值）",
      "基无关" in DOC and "基依赖" in DOC and "上限" in DOC)
# (e) 归属表
check("(e) 三层面判死归属在位（群论 / 一圈动力学 / 味混合 => 公设级买入项）",
      all(k in DOC for k in ["群论", "一圈动力学", "味混合", "公设级买入项"]))
# (f) 七项账
check("(f) 七项账在位（基无关/矩阵刻画/上限/判死/归属/相位/振幅）",
      all(k in DOC for k in ["基无关", "矩阵刻画", "上限", "判死", "归属", "相位", "条件证成"]))

head("F3  锚定")
FILES = ["sm/SMR10_flavor_and_cp.md", "sm/SMR22_gauge_flavor_structure.md",
         "sm/SMR30_koide_block_equipartition.md", "sm/SMR34_equal_weight_relaxation.md",
         "sm/SMR35_sector_selectivity.md", "sm/SMR36_dynamics_selectivity.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 " + str(len(FILES)) + " 个文件存在", not miss, "缺 " + str(miss))
for rel, kws in [("sm/SMR22_gauge_flavor_structure.md", ["U_{uL}"]),
                 ("sm/SMR36_dynamics_selectivity.md", ["零次齐次"]),
                 ("sm/SM_TERMS.md", ["首现双名制"])]:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 " + os.path.basename(rel) + " 含关键词", not m, "缺 " + str(m))
found = []
for dp, dn, fn in os.walk(ROOT):
    if ".git" in dp: continue
    for f in fn:
        if re.match(r"^(Q277|Q473|Q628|build_166)", f): found.append(os.path.join(dp, f))
check("防污染：本仓库内无旧理论文件", not found, "命中 " + str(found))

head("F4  与旧文对齐")
smd = sorted(f for f in os.listdir(HERE) if re.match(r"^SM\d+_.*\.md$", f))
check("旧 SM<数字>_*.md 已按指令删除（0 篇）", len(smd) == 0, "实测 " + str(len(smd)))
check("备份包在位", os.path.exists(os.path.join(HERE, "_legacy_SM0-19.tar.gz")))

head("F5  红线自查")
SKIP = ("**目标**", "**步骤**", "**用到的输入**", "**核对**", "**结果**", "**常见误读**", "**性质**")
DOC_SCAN = DOC.split("## §4")[0] + DOC.split("## §5")[1]
bad = [ln.strip()[:70] for ln in DOC_SCAN.splitlines()
       if not ln.strip().startswith("|") and not ln.strip().startswith("#")
       and not any(ln.strip().startswith(s) for s in SKIP)
       and "**纪律**" not in ln and "F5 红线自查" not in ln
       and ("导出" in ln or "已解决" in ln)
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "需", "登记", "开放",
                                     "判死", "可导出", "买入", "公设", "条件"])]
check("含 导出/已解决 的行都带限定", not bad, "命中 " + str(bad))
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all((str(i) + ". **不许**") in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 " + str(len(dd)))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
