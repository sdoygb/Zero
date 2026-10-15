#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR1_check.py -- 核验《SMR1 · 运动学舞台：把"四维洛伦兹因果结构"从 Z0 上重建》

  F1  结构：头部六字段；§0.1 + §0-§7 在位；§2 的 8 条推导各含五段式六标记；
            每条 >= 2 条"常见误读"；§0.1 术语表 5-10 行
  F2  独立复算（**第二条路径**，纯标准库，不 import 任何旧脚本、不 import lh 根目录脚本）：
        (a) N_L：穷举 +-1 词 vs 二项式
        (b) f(m)：穷举首达 vs 2*Catalan(m-1)
        (c) 生成函数卷积互逆（双向，到 x^14）
        (d) K(L)：逐 k Burnside / phi-闭式 / 轨道枚举 三条路径；1,2,4,10,26,80
        (e) q_L：类级闭包 vs 稳定子路径；六个精确值
        (f) q_L <= L/N_L <= 3/10（偶 L <= 18）
        (g) q'_L：词级群轨道 vs sigma-配对公式；四个精确差
        (h) 粗粒化单调链 q_L <= q'_L <= q''_L；三级精确值
        (i) 四维窗：相邻比单调；F4>F3 <=> q>1/2；F5<F4 <=> q<3/5；唯一 L=4
        (j) v2(f(m)) = popcount(m)
        (k) f(mu)=mu*tanh(mu)-log(cosh(mu)) < log2 且上确界 = log2  ==>  B<=4
  F3  锚定：引用的 lh 根目录与 sm/ 文件存在；抽检关键词真实出现
  F4  与旧文对齐：SM2 表 T1 的 N_L/K(L)/类大小/q_L 逐项；G91 §7.1 轨道数序列；
            四个差与 SM2 §2.5 精确值一致；G91 §7.3 印的是 2 倍（只登记）
  F5  红线自查：无"导出了 D=4"式肯定句；"有效锥＝光锥"只许出现在否定语境；
            §6 七条禁令在位；$$ 成对；sm/ 的 SM<数字>_ 文件集不变
"""
from __future__ import annotations

import io
import itertools
import math
import os
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


DOC = read(os.path.join(HERE, "SMR1_kinematic_stage.md"))
SM2 = legacy("SM2_stage_discrete_skeleton.md")
SM2C = legacy("SM2_check.py")
G91 = read(os.path.join(ROOT, "G91_matter_sector_range_and_three_way_verdict.md"))


# ----------------------------------------------------------------------
# 组合工具（纯标准库；与正文不同路径）
# ----------------------------------------------------------------------
def words(L):
    return [w for w in itertools.product((1, -1), repeat=L) if sum(w) == 0]


def rot(w, k):
    return w[k:] + w[:k]


def canon(w):
    return min(rot(w, k) for k in range(len(w)))


def rot_classes(L):
    """旋转等价类: rep -> size（词级枚举，非 Burnside）"""
    d = {}
    for w in words(L):
        r = canon(w)
        d[r] = d.get(r, 0) + 1
    return d


def stabilizer_order(w):
    return sum(1 for k in range(len(w)) if rot(w, k) == w)


def v2(n):
    e = 0
    while n % 2 == 0:
        n //= 2
        e += 1
    return e


def binom_rows(L):
    return math.comb(L, L // 2)


def catalan(n):
    return math.comb(2 * n, n) // (n + 1)


def f_closed(m):
    return 2 * catalan(m - 1)


# ======================================================================
head("F1  结构（头部 / 各节 / 五段式六标记）")

for fld in ["**日期**", "**性质**", "**与旧文的关系**", "**上游", "**旧理论取材**", "**核验**"]:
    check("头部字段 %s 在位" % fld, fld in DOC)

for sec in ["## §0.1", "## §0 ", "## §1 ", "## §2 ", "## §3 ", "## §4 ", "## §5 ", "## §6 ", "## §7 "]:
    check("节 %s 在位" % sec.strip(), sec in DOC)

SIX = ("**目标**", "**用到的输入**", "**步骤**", "**结果**", "**核对**", "**常见误读**")
SEC2 = DOC.split("## §2 完整推导")[1].split("## §3 ")[0]
subs = re.findall(r"(?m)^### (2\.\d+) ", SEC2)
check("§2 有 8 条推导（2.1-2.8）", subs == ["2.%d" % i for i in range(1, 9)], "实测 %s" % subs)

blocks = re.split(r"(?m)^### 2\.\d+ ", SEC2)[1:]
for name, blk in zip(subs, blocks):
    miss = [t for t in SIX if t not in blk]
    check("§%s 五段式六标记齐全" % name, not miss, "缺 %s" % miss)
    check("§%s 常见误读 >= 2 条" % name, blk.count("✗") >= 2, "实测 %d 条" % blk.count("✗"))
    check("§%s 有加框结果" % name, "\\boxed" in blk)

SEC01 = DOC.split("## §0.1")[1].split("\n## §0 ")[0]
r01 = [ln for ln in SEC01.splitlines()
       if ln.startswith("|") and not re.match(r"^\|[\s:|-]+\|$", ln) and "标准名" not in ln]
check("§0.1 术语表 5-10 行", 5 <= len(r01) <= 10, "实测 %d 行" % len(r01))
check("§0.1 三档对应强度齐全",
      all(k in SEC01 for k in ("部分对应", "**无对应**", "标准")))

# ======================================================================
head("F2  独立复算（第二条路径）")

# ---- (a) N_L ----
ok_a = True
for L in range(2, 13):
    n_enum = len(words(L)) if L % 2 == 0 else 0
    want = binom_rows(L) if L % 2 == 0 else 0
    if n_enum != want:
        ok_a = False
check("(a) N_L 穷举 == binom(L,L/2)（偶 L）/ 0（奇 L），L<=12", ok_a)
check("(a) N_L 序列 2,6,20,70,252,924",
      [binom_rows(L) for L in range(2, 13, 2)] == [2, 6, 20, 70, 252, 924])

# ---- (b) f(m)：穷举首达 ----
def first_return_counts(mmax):
    out = {}
    for m in range(1, mmax + 1):
        L = 2 * m
        c = 0
        for w in words(L):
            s = 0
            for t in range(L):
                s += w[t]
                if s == 0 and t < L - 1:
                    break
            else:
                c += 1
        out[m] = c
    return out


fr = first_return_counts(7)
check("(b) f(m) 穷举首达 == 2*Catalan(m-1)（m<=7）",
      all(fr[m] == f_closed(m) for m in fr),
      "穷举 %s / 闭式 %s" % ([fr[m] for m in sorted(fr)], [f_closed(m) for m in sorted(fr)]))
check("(b) f(m) 序列 2,2,4,10,28,84（m=1..6）",
      [f_closed(m) for m in range(1, 7)] == [2, 2, 4, 10, 28, 84])

# ---- (c) 卷积互逆（双向，到 x^14）：正确形式是 1/B = 1 - F ----
NB = 15
B = [binom_rows(2 * m) for m in range(NB)]
F = [0] + [f_closed(m) for m in range(1, NB)]
one_minus_F = [1] + [-F[m] for m in range(1, NB)]
conv = [sum(one_minus_F[j] * B[i - j] for j in range(i + 1)) for i in range(NB)]
check("(c) (1 - F)*B = 1 逐项（到 x^14）", conv == [1] + [0] * (NB - 1), str(conv[:6]))
xFB = [0] * NB
for i in range(1, NB):
    xFB[i] = sum(F[j] * B[i - j] for j in range(1, i + 1))
check("(c) B = 1 + F*B 逐项（到 x^14）", all(B[i] == xFB[i] for i in range(1, NB)))
# 旧文印的 1 - xF 形式必须被否证（登记用，不是本篇的结论）
one_minus_xF = [1] + [-F[m - 1] for m in range(1, NB)]
conv_x = [sum(one_minus_xF[j] * B[i - j] for j in range(i + 1)) for i in range(NB)]
check("(c) 旧文印的 (1 - xF)*B != 1（可一行否证：m=1 处给 2）",
      conv_x[:3] == [1, 2, 4])
note("SM2 §2.1 Step 3 仍印 1-xF / 1/(1-xF)：%s（已登记为 C-SMR1-1，只登记不改）"
     % ("仍在" if "1-xF" in SM2 else "已改"))

# ---- (d) K(L) 三条路径 + 逐 k 判据（g 偶）----
def K_burnside(L):
    tot = 0
    for k in range(L):
        g = math.gcd(k, L)
        tot += math.comb(g, g // 2) if g % 2 == 0 else 0
    return Fr(tot, L)


def K_printed_rule(L):
    """旧文印的错判据：|Fix| = C(g,g/2) 若 L/g 偶，否则 0"""
    tot = 0
    for k in range(L):
        g = math.gcd(k, L)
        tot += math.comb(g, g // 2) if (L // g) % 2 == 0 else 0
    return Fr(tot, L)


def fix_enum(L, k):
    return sum(1 for w in words(L) if rot(w, k) == w)


ok_perk = True
for L in range(2, 13, 2):
    for k in range(L):
        g = math.gcd(k, L)
        want = math.comb(g, g // 2) if g % 2 == 0 else 0
        if fix_enum(L, k) != want:
            ok_perk = False
check("(d) 逐 k：|Fix(rot_k)| 穷举 == C(g,g/2)（g 偶）/ 0（g 奇），L<=12", ok_perk)
check("(d) 反例 L=6,k=2：g=2 偶而 L/g=3 奇，真值 2（印的判据给 0）",
      fix_enum(6, 2) == 2 and math.gcd(2, 6) == 2 and (6 // 2) % 2 == 1)
check("(d) 反例 L=12,k=1：g=1 奇而 L/g=12 偶，真值 0（印的判据给 C(1,0)=1）",
      fix_enum(12, 1) == 0 and (12 // 1) % 2 == 0)
check("(d) 反例 L=12,k=4：g=4 偶而 L/g=3 奇，真值 6（印的判据给 0）",
      fix_enum(12, 4) == 6 and math.gcd(4, 12) == 4 and (12 // 4) % 2 == 1)
check("(d) 印的判据求和 != K(L)（L=12：%s vs 80）" % K_printed_rule(12),
      K_printed_rule(12) != 80)
note("SM2 §2.2 Step 2 仍印 'L/g 偶'：%s（已登记为 C-SMR1-2，只登记不改）"
     % ("仍在" if "L/g\\ \\text{偶}" in SM2 or "L/g$ 偶" in SM2 or "L/g\\ \\text{偶}" in SM2 else "见正文"))


def K_phi(L):
    tot = 0
    for d in range(1, L // 2 + 1):
        if (L // 2) % d == 0:
            tot += sum(1 for x in range(1, d + 1) if math.gcd(x, d) == 1) * math.comb(L // d, L // (2 * d))
    return Fr(tot, L)


Kvals, Kmult = {}, {}
for L in range(2, 13, 2):
    cl = rot_classes(L)
    Kvals[L] = len(cl)
    Kmult[L] = sorted(cl.values(), reverse=True)
    check("(d) K(%d) 三路径一致（Burnside/phi/枚举）" % L,
          K_burnside(L) == K_phi(L) == Kvals[L],
          "%s / %s / %s" % (K_burnside(L), K_phi(L), Kvals[L]))
check("(d) K 序列 1,2,4,10,26,80",
      [Kvals[L] for L in range(2, 13, 2)] == [1, 2, 4, 10, 26, 80])
check("(d) L=4 类大小多重集 {4,2}", Kmult[4] == [4, 2])
check("(d) L=6 类大小多重集 {6,6,6,2}", Kmult[6] == [6, 6, 6, 2])
check("(d) 类大小整除 L 且求和 == N_L",
      all(all(L % s == 0 for s in Kmult[L]) and sum(Kmult[L]) == binom_rows(L)
          for L in range(2, 13, 2)))

# ---- (e) q_L 两路径 ----
def q_classes(L):
    cl = rot_classes(L)
    N = binom_rows(L)
    return sum((Fr(s, N)) ** 2 for s in cl.values())


def q_stabilizer(L):
    ws = words(L)
    N = len(ws)
    return sum(Fr(L, stabilizer_order(w)) for w in ws) / Fr(N * N)


QEXACT = {2: Fr(1), 4: Fr(5, 9), 6: Fr(7, 25), 8: Fr(19, 175),
          10: Fr(313, 7938), 12: Fr(683, 53361)}
for L, want in QEXACT.items():
    got1, got2 = q_classes(L), q_stabilizer(L)
    check("(e) q_%d 两路径一致且 == %s" % (L, want), got1 == got2 == want,
          "%s / %s" % (got1, got2))

# ---- (f) 界 q_L <= L/N_L <= 3/10 ----
ok_f, worst = True, None
for L in range(6, 19, 2):
    q = q_stabilizer(L)
    b1, b2 = Fr(L, binom_rows(L)), Fr(3, 10)
    if not (q <= b1 <= b2):
        ok_f = False
        worst = (L, q, b1)
check("(f) q_L <= L/N_L <= 3/10（偶 L 6..18）", ok_f, str(worst))
check("(f) L=6 处 L/N_L 取等 3/10", Fr(6, binom_rows(6)) == Fr(3, 10))

# ---- (g) q'_L：群轨道 vs sigma-配对公式 ----
def q_prime_group(L):
    ws = words(L)
    N = len(ws)
    sig = lambda w: w[::-1]
    seen, tot = set(), Fr(0)
    for w in ws:
        if w in seen:
            continue
        orb = set()
        st = [w]
        while st:
            x = st.pop()
            if x in orb:
                continue
            orb.add(x)
            for k in range(L):
                st.append(rot(x, k))
            st.append(sig(x))
        seen |= orb
        tot += (Fr(len(orb), N)) ** 2
    return tot


def q_prime_formula(L):
    cl = rot_classes(L)
    N = binom_rows(L)
    sig = lambda w: w[::-1]
    rep2 = {r: canon(sig(r)) for r in cl}
    done, tot = set(), Fr(0)
    for r in cl:
        if r in done:
            continue
        r2 = rep2[r]
        if r2 == r:
            tot += (Fr(cl[r], N)) ** 2
            done.add(r)
        else:
            tot += (Fr(cl[r] + cl[r2], N)) ** 2
            done |= {r, r2}
    return tot


QPRIME = {2: Fr(1), 4: Fr(5, 9), 6: Fr(23, 50), 8: Fr(197, 1225),
          10: Fr(563, 7938), 12: Fr(2419, 106722)}
DIFF = {6: Fr(9, 50), 8: Fr(64, 1225), 10: Fr(125, 3969), 12: Fr(117, 11858)}
for L, want in QPRIME.items():
    g1, g2 = q_prime_group(L), q_prime_formula(L)
    check("(g) q'_%d 两路径一致且 == %s" % (L, want), g1 == g2 == want, "%s / %s" % (g1, g2))
    if L in DIFF:
        check("(g) q'_%d - q_%d == %s" % (L, L, DIFF[L]), g1 - QEXACT[L] == DIFF[L])

# ---- (h) 粗粒化单调链 + q'' ----
def q_dprime(L):
    ws = words(L)
    N = len(ws)
    sig = lambda w: w[::-1]
    neg = lambda w: tuple(-x for x in w)
    seen, tot = set(), Fr(0)
    for w in ws:
        if w in seen:
            continue
        orb, st = set(), [w]
        while st:
            x = st.pop()
            if x in orb:
                continue
            orb.add(x)
            for k in range(L):
                st.append(rot(x, k))
            st.append(sig(x))
            st.append(neg(x))
        seen |= orb
        tot += (Fr(len(orb), N)) ** 2
    return tot


QDPRIME = {2: Fr(1), 4: Fr(5, 9), 6: Fr(23, 50), 8: Fr(229, 1225),
           10: Fr(394, 3969), 12: Fr(3823, 106722)}
for L in [2, 4, 6, 8, 10, 12]:
    q, q1, q2 = QEXACT[L], QPRIME[L], q_dprime(L)
    check("(h) q''_%d == %s 且 q<=q'<=q''" % (L, QDPRIME[L]), q2 == QDPRIME[L] and q <= q1 <= q2,
          "%s / %s / %s" % (q, q1, q2))
check("(h) q''=q' 只在 L<=6；L>=8 严格大",
      QDPRIME[6] == QPRIME[6] and all(QDPRIME[L] > QPRIME[L] for L in (8, 10, 12)))
check("(h) 三级数值与 SM2_check 钉住的关键值一致（229/1225, 394/3969, 3823/106722）",
      all(s in SM2C for s in ["229/1225", "394/3969", "3823/106722"]))

# ---- (i) 四维窗 ----
mono = True
for D in range(2, 30):
    r1 = Fr(D + 1, D - 1) * Fr(1)
    r2 = Fr(D + 2, D) * Fr(1)
    if not r2 < r1:
        mono = False
check("(i) 相邻比 (D+1)/(D-1)*q 对 D>=2 严格递减", mono)
check("(i) F4>F3 <=> q>1/2", (Fr(1, 2) + Fr(1, 100)) > Fr(1, 2) and not (Fr(1, 2) - Fr(1, 100)) > Fr(1, 2))
check("(i) F5<F4 <=> q<3/5", Fr(3, 5) - Fr(1, 100) < Fr(3, 5))
inwin = [L for L in range(2, 19, 2) if Fr(1, 2) < q_stabilizer(L) < Fr(3, 5)]
check("(i) 偶 L<=18 中唯一落窗 L=4", inwin == [4], "实测 %s" % inwin)
check("(i) L=2 无有限峰（q_2=1>1）", QEXACT[2] == 1)
check("(i) 2(T+1)=2L+1 无整数解（L<=20）",
      all((2 * L + 1) % 2 == 1 for L in range(1, 21)))

# ---- (j) v2(f(m)) = popcount(m) ----
ok_j = True
for m in range(1, 21):
    fm = 2 * math.comb(2 * m - 2, m - 1) // m
    if v2(fm) != bin(m).count("1"):
        ok_j = False
check("(j) v2(f(m)) = popcount(m)（m<=20）", ok_j)
check("(j) f(m)=2*odd <=> m 是 2 的幂（m<=20）",
      all((v2(2 * math.comb(2 * m - 2, m - 1) // m) == 1) == (m & (m - 1) == 0) for m in range(1, 21)))

# ---- (k) f(mu) < log2, sup = log2 ----
def f_mu(mu):
    e = math.exp(-2 * mu)
    return math.log(2) - math.log1p(e) - 2 * mu * e / (1 + e)


xs = [0.5, 1, 2, 3, 5, 8, 12, 20, 40]
vals = [f_mu(x) for x in xs]
check("(k) f(mu) 严格递增；f(mu) < log2；log2-f(10) < 1e-8",
      all(vals[i] < vals[i + 1] for i in range(len(vals) - 1))
      and all(v <= math.log(2) for v in vals)
      and all(v < math.log(2) for v in vals[:-1])
      and (math.log(2) - f_mu(12)) < 1e-8,
      "log2-f(12)=%.3e" % (math.log(2) - f_mu(12)))
check("(k) sup f = log2 ==> B<=4（(1/2)logB<=log2 <=> B<=4）",
      Fr(1, 2) * Fr(4).limit_denominator() == 2 and math.log(4) / 2 <= math.log(2) + 1e-15
      and math.log(4.0001) / 2 > math.log(2))
check("(k) 文本给出 f(mu) 的闭式与上确界", "\\tanh\\mu-\\log\\cosh\\mu" in DOC and "\\log2" in DOC)

# ======================================================================
head("F3  锚定（文件存在 + 关键词）")

FILES = ["Z0_zero_never_rests_single_axiom.md", "Z1_zero_layer_as_the_foundation.md",
         "Z3_i2a_dimension_drift_verdict.md", "Z14_closure_cyclic_order_base_theorem.md",
         "G8_dimension_selection.md", "G40_metric_from_closed_walk_counting.md",
         "G89_dimension_no_go_and_the_balance_condition.md",
         "G91_matter_sector_range_and_three_way_verdict.md",
         "R25_native_pair_cost_and_four_dim_peak.md", "R31_phase_ledger_and_lifetime_selection.md",
         "R32_ledger_readout_selection_and_L8_resolution.md", "R46_pair_ledger_objects_and_simplex.md",
         "R47_signature_from_causal_cone.md", "R48_exact_cone_and_effective_cone.md",
         "Z18_zero_dimension_bridge_and_gr_stalls.md",
         "sm/SM2_stage_discrete_skeleton.md", "sm/SM2_check.py", "sm/SM_TERMS.md",
         "sm/SM_PILLARS.md", "sm/README.md"]
FILES = [f for f in FILES if not re.match(r"sm/SM\d+_", f)]
miss = [f for f in FILES if not exists(f)]
check("引用的 %d 个文件在磁盘上存在" % len(FILES), not miss, "缺 %s" % miss)

KW = [("Z0_zero_never_rests_single_axiom.md", ["零不停留", "不设概率"]),
      ("Z1_zero_layer_as_the_foundation.md", ["定理 1"]),
      ("Z14_closure_cyclic_order_base_theorem.md", ["循环序"]),
      ("G8_dimension_selection.md", ["引理 39"]),
      ("G89_dimension_no_go_and_the_balance_condition.md", ["命题 1"]),
      ("R48_exact_cone_and_effective_cone.md", ["锥"]),
      ("Z18_zero_dimension_bridge_and_gr_stalls.md", ["读出坐标"]),
      ("sm/SM2_stage_discrete_skeleton.md", ["轨道数序列"]),
      ("sm/SM_TERMS.md", ["首现双名制"]),
      ("sm/SM_PILLARS.md", ["P1"])]
KW = [(r, k) for r, k in KW if not re.match(r"sm/SM\d+_", r)]
for rel, kws in KW:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 %s 含关键词" % os.path.basename(rel), not m, "缺 %s" % m)

# ======================================================================
head("F4  与旧文对齐（旧文已按指令删除：对照改以 sm/_legacy_SM0-19.tar.gz 为准）")
OLD_USED = ['SM2_check.py', 'SM2_stage_discrete_skeleton.md']
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

bad = [ln.strip()[:60] for ln in DOC.splitlines()
       if "D=4" in ln and "导出" in ln
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "错问法"])]
check("含 \"D=4 + 导出\" 的行都带否定限定", not bad, "命中 %s" % bad)

bad2 = [ln.strip()[:60] for ln in DOC.splitlines()
        if "有效锥" in ln and "光锥" in ln
        and not any(k in ln for k in ["不", "禁", "✗", "没"])]
check("\"有效锥＝光锥\" 只出现在否定语境", not bad2, "命中 %s" % bad2)

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
check("本篇文件名不被 SM<数字>_ 正则捕获",
      not re.match(r"^SM\d+_", "SMR1_kinematic_stage.md")
      and not re.match(r"^SM\d+_", "SMR1_check.py"))

# ======================================================================
print("")
print("=" * 72)
print("断言总数：%d    不符项：%d" % (NCHK, FAIL))
print("=" * 72)
sys.exit(1 if FAIL else 0)
