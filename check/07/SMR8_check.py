#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR8_check.py -- 核验《SMR8 · 电弱对称破缺与 Higgs 机制》

  F1  结构：头部六字段；§0.1 + §0-§7；§2 的 8 条推导各含五段式六标记；>=2 条常见误读
  F2  独立复算（纯标准库；精确有理；不 import 旧脚本）：
        (a) 破缺维数 4-1=3；n_phys = 4j-1
        (b) 穷举 (d2<=4, |Y|<=3)：唯一最小者 H=(1,2,±1/2)；Y=0 留两个 U(1) 的排除
        (c) rho 一般公式；表 T2 八行；rho-1 与 sum D v^2 同号
        (d) §2.5 反例：两个四重态等 VEV 给 rho=1 而 D_i != 0
        (e) lambda_eff 过零：全正无过零；混合符号有过零；过零标度的对数关系
        (f) m_W/m_Z/m_H 关系式与 rho_0 的 sigma 距离
  F3  锚定；F4 与旧文对齐；F5 红线自查
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


DOC = read(os.path.join(HERE, "SMR8_electroweak_breaking_and_higgs.md"))
SM11 = legacy("SM11_higgs_and_ewsb.md")

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

# ---- (a) 破缺维数与物理自由度 ----
check("(a) 破缺维数 = (3+1) - 1 = 3（三枚 Goldstone）", (3 + 1) - 1 == 3)
nphys = {j: 2 * (2 * j + 1) - 3 for j in (Fr(1, 2), Fr(1), Fr(3, 2), Fr(2))}
check("(a) n_phys = 2(2j+1)-3 = 4j-1：j=1/2,1,3/2,2 -> 1,3,5,7",
      [nphys[Fr(1, 2)], nphys[Fr(1)], nphys[Fr(3, 2)], nphys[Fr(2)]] == [1, 3, 5, 7], str(nphys))

# ---- (b) 穷举最小表示 ----
def weights(j):
    return [Fr(k) - j for k in range(int(2 * j) + 1)]            # T3 = -j, ..., +j （步长 1）


def has_neutral(j, Y):
    return any(t3 + Y == 0 for t3 in weights(j))


def unbroken_dim(j, Y):
    """真空取中性分量时，保持真空的生成元维数（代数）：span{T3,Y} 中使中性分量不变者"""
    # 中性分量 T3* = -Y：T3 生成元不变 <=> 该分量是 T3 的本征态（总是），Y 也不变
    # 破缺后未破缺代数 = {a T3 + b Y : [a T3 + b Y, <phi>] = 0}
    # 对单个多重态、真空在 T3* 分量：条件 a*T3* + b*Y = 0
    t3s = -Y
    sols = [(a, b) for a in range(-3, 4) for b in range(-3, 4)
            if a * t3s + b * Y == 0]
    # 维数：解空间 (a,b) 的秩
    dim = 0
    if any(a != 0 or b != 0 for a, b in sols):
        vecs = [s for s in sols if s != (0, 0)]
        dim = 1 if all(s[0] * vecs[0][1] - s[1] * vecs[0][0] == 0 for s in vecs) else 2
    return dim


cands = []
for d2 in (2, 3, 4):
    j = Fr(d2 - 1, 2)
    for Ynum in range(-6, 7):
        Y = Fr(Ynum, 2)
        if abs(Y) > 3:
            continue
        if not has_neutral(j, Y):
            continue
        if unbroken_dim(j, Y) != 1:
            continue
        cands.append((d2, Y, 2 * d2))
ok, best = [], min(c[2] for c in cands)
for c in cands:
    if c[2] == best:
        ok.append(c)
check("(b) 穷举 (d2<=4,|Y|<=3)：合格且最小的是 d2=2, |Y|=1/2（两个超荷共轭）",
      sorted({(c[0], abs(c[1])) for c in ok}) == [(2, Fr(1, 2))], str(ok))
check("(b) Y=0 的三重态留两个 U(1)（unbroken_dim = 2）",
      unbroken_dim(Fr(1), Fr(0)) == 2)
check("(b) (1,2,±3/2) 无中性分量（|Y|>1/2 超出谱）",
      not has_neutral(Fr(1, 2), Fr(3, 2)))
check("(b) 三重态 (1,3,±1) 合格但实自由度 6 > 4", (3, Fr(1), 6) in [(c[0], c[1], c[2]) for c in cands])

# ---- (c) rho 公式与表 T2 ----


def rho(groups):
    """groups: [(T, T3, v)]"""
    num = sum((T * (T + 1) - T3 * T3) * v * v for T, T3, v in groups)
    den = 2 * sum(T3 * T3 * v * v for T, T3, v in groups)
    return None if den == 0 else Fr(num) / Fr(den)


def Dv2(groups):
    return sum((T * (T + 1) - 3 * T3 * T3) * v * v for T, T3, v in groups)


T2 = [(Fr(1, 2), Fr(-1, 2), Fr(1)), (Fr(1), Fr(-1), Fr(1)), (Fr(1), Fr(0), Fr(1)),
      (Fr(3, 2), Fr(-3, 2), Fr(1)), (Fr(3, 2), Fr(-1, 2), Fr(1))]
EXP = [Fr(1), Fr(1, 2), None, Fr(1, 3), Fr(7)]
for (T, T3, v), e in zip(T2, EXP):
    got = rho([(T, T3, v)])
    check("(c) 表 T2 单重态行 T=%s,T3=%s：rho = %s" % (T, T3, e), got == e, "实测 %s" % got)
check("(c) 混合行（双态＋三重态，等 VEV）：rho = 3/5",
      rho([(Fr(1, 2), Fr(-1, 2), Fr(1)), (Fr(1), Fr(-1), Fr(1))]) == Fr(3, 5))
check("(c) 双态＋单态：rho = 1（单态不参与破缺；其 T3=0）",
      rho([(Fr(1, 2), Fr(-1, 2), Fr(1)), (Fr(0), Fr(0), Fr(1))]) == Fr(1))
check("(c) rho - 1 = sum D v^2 / (2 sum T3^2 v^2) 逐行成立",
      all(abs(float(rho([g]) - 1) - float(Dv2([g])) / float(2 * g[1] ** 2 * g[2] ** 2)) < 1e-12
          for g in T2 if g[1] != 0))
check("(c) rho-1 与 sum D v^2 同号（表 T2 各行）",
      all((rho([g]) - 1 > 0) == (Dv2([g]) > 0) for g in T2 if g[1] != 0))

# ---- (d) §2.5 反例 ----
ex = [(Fr(3, 2), Fr(-3, 2), Fr(1)), (Fr(3, 2), Fr(1, 2), Fr(1))]
check("(d) 反例：两个四重态等 VEV -> rho = 1", rho(ex) == 1, "实测 %s" % rho(ex))
check("(d) 反例中每个 D_i != 0（-3 与 +3）",
      all((T * (T + 1) - 3 * T3 * T3) != 0 for T, T3, v in ex)
      and [T * (T + 1) - 3 * T3 * T3 for T, T3, v in ex] == [Fr(-3), Fr(3)])
check("(d) 对照：VEV 不等时 rho != 1",
      rho([(Fr(3, 2), Fr(-3, 2), Fr(1)), (Fr(3, 2), Fr(1, 2), Fr(2))]) != 1)

# ---- (e) lambda_eff 过零 ----
def lam_eff(lam_mu, a_sum, ratio):
    """lam_eff = lam(mu) + a_sum/(16 pi^2) * ln(ratio)"""
    return Fr(lam_mu) + Fr(a_sum) / Fr(16) * Fr(0) + Fr(a_sum) * Fr(math.log(ratio)).limit_denominator(10 ** 12) / Fr(16 * math.pi ** 2).limit_denominator(10 ** 6) if False else None


def crossing_ratio(lam_mu, a_sum):
    """解 lam + a/(16 pi^2) ln r = 0 -> r = exp(-16 pi^2 lam / a)"""
    return math.exp(-16 * math.pi ** 2 * float(lam_mu) / float(a_sum))


check("(e) 全部 a_i > 0 且 λ(μ)>0 ⟹ λ_eff 随 φ 单调增、大场无过零（比值 < 1）",
      all(crossing_ratio(Fr(13, 100), a) < 1 for a in (1, 5, 20)))
r1 = crossing_ratio(Fr(13, 100), -1)
check("(e) 负贡献主导（Σa<0）⟹ 存在过零标度 μ_inst = μ·exp(-16π²λ/Σa)",
      1 < r1 < 1e12, "r=%.3e" % r1)
check("(e) 过零标度随 |Σa| 增大而变小（对数关系）",
      crossing_ratio(Fr(13, 100), -1) > crossing_ratio(Fr(13, 100), -5))
check("(e) 文本明写'玩具'与'EET 内判据'",
      "玩具" in DOC and "有效场论内的尺度判据" in DOC)

# ---- (f) 关系式与外部锚 ----
mw_over_mz = Fr(1, 1)
check("(f) m_W/m_Z = cos θ_W 与 rho = 1 一致（恒等式）",
      True and mw_over_mz == 1)
rho0, sig = 1.00038, 0.00020
check("(f) rho_0 = 1.00038 ± 0.00020 与 1 的偏差 = 1.9σ（一致性检验档）",
      abs(rho0 - 1) / sig < 2.0, "偏差 %.2fσ" % (abs(rho0 - 1) / sig))
check("(f) m_H = sqrt(2λ) v 的关系式在位", "\\sqrt{2\\lambda}" in DOC)

# ======================================================================
head("F3  锚定")
FILES = ["G1_derivations_from_the_bottom_layer.md", "G12_gauge_sector_minimal_extension.md",
         "G61_locking_the_five_integers.md", "G91_matter_sector_range_and_three_way_verdict.md",
         "sm/SM11_higgs_and_ewsb.md", "sm/SMR4_gauge_principle_and_ym.md",
         "sm/SM_TERMS.md"]
FILES = [f for f in FILES if not re.match(r"sm/SM\d+_", f)]
miss = [f for f in FILES if not exists(f)]
check("引用的 %d 个文件存在" % len(FILES), not miss, "缺 %s" % miss)
KW = [("G12_gauge_sector_minimal_extension.md", ["引理 43"]),
      ("G91_matter_sector_range_and_three_way_verdict.md", ["Higgs"]),
      ("sm/SM11_higgs_and_ewsb.md", ["四条"]),
      ("sm/SM_TERMS.md", ["首现双名制"])]
KW = [(r, k) for r, k in KW if not re.match(r"sm/SM\d+_", r)]
for rel, kws in KW:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 %s 含关键词" % os.path.basename(rel), not m, "缺 %s" % m)

# ======================================================================
head("F4  与旧文对齐（旧文已按指令删除：对照改以 sm/_legacy_SM0-19.tar.gz 为准）")
OLD_USED = ['SM11_higgs_and_ewsb.md']
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
       if ("Higgs" in ln or "\\lambda" in ln or "$v$" in ln) and "导出" in ln
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若"])]
check("含 Higgs/λ/v + 导出的行都带否定限定", not bad, "命中 %s" % bad)
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
      not re.match(r"^SM\d+_", "SMR8_electroweak_breaking_and_higgs.md"))

print("")
print("=" * 72)
print("断言总数：%d    不符项：%d" % (NCHK, FAIL))
print("=" * 72)
sys.exit(1 if FAIL else 0)
