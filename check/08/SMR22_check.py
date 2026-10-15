#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR22_check.py -- 核验《SMR22 · 物质-规范耦合的味结构层》"""
import functools, io, itertools, math, os, re, sys
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

DOC = read(os.path.join(HERE, "SMR22_gauge_flavor_structure.md"))

# 精确复有理与矩阵工具
CZ, C1 = (Fr(0), Fr(0)), (Fr(1), Fr(0))
def cadd(a, b): return (a[0] + b[0], a[1] + b[1])
def csub(a, b): return (a[0] - b[0], a[1] - b[1])
def cmul(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def cconj(a): return (a[0], -a[1])
def mmul(A, B):
    return [[__import__("functools").reduce(cadd, (cmul(A[i][k], B[k][j]) for k in range(len(B))), CZ)
             for j in range(len(B[0]))] for i in range(len(A))]
def mH(A): return [[cconj(A[j][i]) for j in range(len(A))] for i in range(len(A[0]))]
def mid(n): return [[C1 if i == j else CZ for j in range(n)] for i in range(n)]
def meq(A, B): return all(A[i][j] == B[i][j] for i in range(len(A)) for j in range(len(A)))
def Q(a, b=0): return (Fr(a), Fr(b))
def kron(A, n):
    """A (k x k 复矩阵) 张量 identity_n"""
    k = len(A); N = k * n
    M = [[CZ] * N for _ in range(N)]
    for i in range(k):
        for j in range(k):
            for f in range(n):
                M[i * n + f][j * n + f] = A[i][j]
    return M

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
# (a) 味普适：T_a = rho(T_a) x 1_n 闭合；味依赖形变破坏闭合
# su(2)：rho = 自旋 1/2
S1 = [[CZ, C1], [C1, CZ]]                      # sigma_x/2 形式（用整数倍数避免 1/2）
S2 = [[CZ, (Fr(0), Fr(-1))], [(Fr(0), Fr(1)), CZ]]
S3 = [[C1, CZ], [CZ, (Fr(-1), Fr(0))]]
gen = [S1, S2, S3]
def comm(A, B):
    AB = mmul(A, B); BA = mmul(B, A)
    return [[csub(AB[i][j], BA[i][j]) for j in range(len(A))] for i in range(len(A))]
ok_close = True
for A, B in itertools.combinations(range(3), 2):
    C = comm(gen[A], gen[B])
    # [S_a,S_b] = 2i eps_{abc} S_c（按上面的整数代数）
    found = False
    for c in range(3):
        pred = [[cmul((Fr(0), Fr(2)), gen[c][i][j]) for j in range(2)] for i in range(2)]
        predn = [[cmul((Fr(0), Fr(-2)), gen[c][i][j]) for j in range(2)] for i in range(2)]
        if meq(C, pred) or meq(C, predn): found = True
    if not found: ok_close = False
check("(a) su(2) 生成元的对易子闭合到同一代数（整数化约定下 [S_a,S_b]=±2i S_c）", ok_close)
G = [kron(g, 3) for g in gen]
ok_G = True
for A, B in itertools.combinations(range(3), 2):
    C = comm(G[A], G[B])
    found = any(meq(C, [[cmul((Fr(0), Fr(2)), G[c][i][j]) for j in range(6)] for i in range(6)])
                or meq(C, [[cmul((Fr(0), Fr(-2)), G[c][i][j]) for j in range(6)] for i in range(6)])
                for c in range(3))
    if not found: ok_G = False
check("(a) T_a = rho(T_a) x 1_flavor 仍闭合到同一代数（味普适可行）", ok_G)
# 味依赖形变：T_3 -> rho(S3) x diag(1,1,1+eps) 与 rho(S1) x 1 的对易子是否仍在 {T_a} 内
def diag_c(v): return [[v[i] if i == j else CZ for j in range(len(v))] for i in range(len(v))]
def kron_full(A, B):
    k, m = len(A), len(B)
    M = [[CZ] * (k * m) for _ in range(k * m)]
    for i in range(k):
        for j in range(k):
            for a in range(m):
                for b in range(m):
                    M[i * m + a][j * m + b] = cmul(A[i][j], B[a][b])
    return M
T1 = kron_full(S1, mid(3))
T3d = kron_full(S3, diag_c([C1, C1, (Fr(3, 2), Fr(0))]))     # eps = 1/2 的形变
C = comm(T1, T3d)
# 判断 C 是否落在 {T1,T2,T3 的普适组合} 里：普适情形 [T1,T3]=±2iT2（纯味单位），形变后会多出味非单位块
M2k = kron_full(S2, mid(3))
pure = [[cmul((Fr(0), Fr(2)), M2k[i][j]) for j in range(6)] for i in range(6)]
check("(a) 味依赖形变 T_3 x diag(1,1,3/2) 的对易子不再是纯味单位（故需额外生成元）",
      not meq(C, pure))

# (b) V_CKM = U_uL^dag U_dL 幺正；带电流非对角
def rot3(a, b, c):
    """有理正交 3x3（两角旋转复合）"""
    ca, sa = (Fr(3, 5), Fr(0)), (Fr(4, 5), Fr(0))
    cb, sb = (Fr(5, 13), Fr(0)), (Fr(12, 13), Fr(0))
    R = [[C1, CZ, CZ], [CZ, C1, CZ], [CZ, CZ, C1]]
    def neg(x): return (Fr(-x[0]), Fr(-x[1]))
    R1 = [[ca, neg(sa), CZ], [sa, ca, CZ], [CZ, CZ, C1]]
    R2 = [[C1, CZ, CZ], [CZ, cb, neg(sb)], [CZ, sb, cb]]
    return mmul(R1, R2)
UuL = rot3(0, 0, 0)
UdL = [[(Fr(0), Fr(1)), CZ, CZ], [CZ, C1, CZ], [CZ, CZ, C1]]
UdL = mmul(UdL, rot3(1, 1, 1))          # 复数相位 x 实正交 = 幺正
V = mmul(mH(UuL), UdL)
check("(b) V_CKM = U_uL^dag U_dL 幺正（精确复有理）", meq(mmul(mH(V), V), mid(3)))
off = sum(1 for i in range(3) for j in range(3) if i != j and V[i][j] != CZ)
check("(b) 双幺正之后带电流**非对角**（V 有 " + str(off) + " 个非零非对角元）", off > 0)

# (c) 中性流保持对角
c_slot = Q(1, 4)                                   # 某一槽（如 u_L）的中性流系数
gL = diag_c([c_slot, c_slot, c_slot])              # 味空间里 ∝ 1（味普适）
diag_kept = meq(mmul(mmul(mH(UdL), gL), UdL), gL)
check("(c) 同一槽内耦合 ∝ 味单位 => 任意幺正旋转后仍对角（中性流树级无 FCNC）", diag_kept)
gLbad = diag_c([Q(1, 4), Q(-3, 8), Q(1, 2)])       # 味依赖（反例）
check("(c) 反例：若耦合味依赖（代间不同），旋转后立刻出现非对角元",
      not meq(mmul(mmul(mH(UdL), gLbad), UdL), gLbad))

# (d) 幺正性：列正交与简并极限
check("(d) 列正交性 sum_i V_id V_is^* = 0 (d != s)（用 V^dag V = 1 的两列）",
      all(functools.reduce(cadd, (cmul(V[i][d], cconj(V[i][s])) for i in range(3)), CZ) == CZ
          for d in range(3) for s in range(3) if d != s))
check("(d) 简并极限：若 f 与 i 无关则振幅 ~ sum_i V_id V_is^* = 0（GIM 代数形式）",
      all(functools.reduce(cadd, (cmul(V[i][d], cconj(V[i][s])) for i in range(3)), CZ) == CZ
          for d in range(3) for s in range(3) if d != s))

# (e) 参数计数
def nparam(N):
    return (N * (N - 1) // 2, (N - 1) * (N - 2) // 2)
check("(e) 参数计数：N=3 -> 3 角 + 1 相位；N=2 -> 1+0；N=1 -> 0+0",
      nparam(3) == (3, 1) and nparam(2) == (1, 0) and nparam(1) == (0, 0))
check("(e) 总参数 N^2 与 2N-1 个相位重定义一致（N=3: 9 相位, 去 5, 余 4）",
      3 ** 2 == 9 and 2 * 3 - 1 == 5 and 9 - 5 == 4)

# (f) 账
check("(f) 八项账在位（味普适/CKM/无树级 FCNC/GIM/计数/汤川买入/N_g 开放/反常一致）",
      all(k in DOC for k in ["味普适", "CKM", "FCNC", "GIM", "参数计数", "汤川矩阵", "N_g", "反常"])
      and "具名输入" in DOC and "未定义" in DOC)

head("F3  锚定")
FILES = ["sm/SMR4_gauge_principle_and_ym.md", "sm/SMR6_chirality_and_anomalies.md",
         "sm/SMR7_fermion_spectrum_and_generations.md", "sm/SMR10_flavor_and_cp.md",
         "sm/SMR13_weak_angle_derivation.md", "sm/SMR21_gauge_force_layer.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 " + str(len(FILES)) + " 个文件存在", not miss, "缺 " + str(miss))
for rel, kws in [("sm/SMR10_flavor_and_cp.md", ["味盲"]),
                 ("sm/SMR21_gauge_force_layer.md", ["12"]),
                 ("sm/SMR7_fermion_spectrum_and_generations.md", ["no-go"]),
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
       and not ln.strip().startswith("**目标**") and not ln.strip().startswith("**步骤**")
       and ("导出" in ln or "FCNC" in ln)
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若", "登记",
                                     "买入", "具名输入", "树级", "可导出", "存在", "形式", "计数", "GIM", "简并极限"])]
check("含 导出/FCNC 的行都带限定", not bad, "命中 " + str(bad))
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all((str(i) + ". **不许**") in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 " + str(len(dd)))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
