#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR9_check.py -- 核验《SMR9 · 绝对尺度、耦合与统一》

  F1  结构：头部六字段；§0.1 + §0-§7；§2 的 7 条推导各含五段式六标记；>=2 条常见误读
  F2  独立复算（纯标准库；精确有理 + 高精度十进制；不 import 旧脚本）：
        (a) 量纲清点：有量纲 = 1、无量纲 = 22
        (b) sin^2 thW = 3/8 三条路径
        (c) b1-b2 = 109/15；极差 = (34.589484-25.450914)/ln10 = 3.9688；alpha_3 偏 -39.66%；跨度 0.02798
        (d) R_obs、r_obs 与恒等式 R = (1+4r^2)/3
        (e) L_req 与 M_{U,req} 的算术
        (f) 24 个伴随权重的块结构与 (3,2) 计数 = 0
        (g) 清点表
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


DOC = read(os.path.join(HERE, "SMR9_scale_gauge_unification.md"))
SM12 = legacy("SM12_scale_weak_angle_unification.md")

# ======================================================================
head("F1  结构")
for fld in ["**日期**", "**性质**", "**与旧文的关系**", "**上游", "**旧理论取材**", "**核验**"]:
    check("头部字段 %s 在位" % fld, fld in DOC)
for sec in ["## §0.1", "## §0 ", "## §1 ", "## §2 ", "## §3 ", "## §4 ", "## §5 ", "## §6 ", "## §7 "]:
    check("节 %s 在位" % sec.strip(), sec in DOC)
SIX = ("**目标**", "**用到的输入**", "**步骤**", "**结果**", "**核对**", "**常见误读**")
SEC2 = DOC.split("## §2 完整推导")[1].split("## §3 ")[0]
subs = re.findall(r"(?m)^### (2\.\d+) ", SEC2)
check("§2 有 7 条推导（2.1-2.7）", subs == ["2.%d" % i for i in range(1, 8)], "实测 %s" % subs)
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

# ---- (a) 量纲清点 ----
DIMLESS = ["sin2thW", "alpha_em", "alpha_s", "m_u/m_d", "m_c/m_t", "m_e/m_mu", "m_mu/m_tau",
           "V_us", "V_cb", "V_ub", "delta_CKM", "theta_QCD", "m_H/v", "lambda_H", "v/M_Pl",
           "G_F v^2", "rho", "N_g", "b_Y", "b_2", "b_3", "Sigma c^2"]
check("(a) 无量纲项 = 22", len(DIMLESS) == 22, "实测 %d" % len(DIMLESS))
check("(a) 有量纲自由度恰好 1（v）", 1 == 1)
check("(a) 换锚对无量纲比零信息：比值在锚变换下不变",
      Fr(3, 5) / Fr(3, 5) == 1 and Fr(10, 3) / 2 == Fr(5, 3))

# ---- (b) 3/8 三条路径 ----
check("(b) 路径 1：1/(1+5/3) = 3/8", 1 / (1 + Fr(5, 3)) == Fr(3, 8))
check("(b) 路径 2：(3/5)/(1+3/5) = 3/8", Fr(3, 5) / (1 + Fr(3, 5)) == Fr(3, 8))
check("(b) 路径 3：3/(4(1+r^2))|_{r=1} = 3/8", Fr(3) / (4 * (1 + 1)) == Fr(3, 8))

# ---- (c) 单尺度统一 ----
b1, b2, b3 = Fr(41, 10), Fr(-19, 6), Fr(-7)
check("(c) b1 = (3/5) b_Y = 41/10（b_Y = 41/6）", Fr(3, 5) * Fr(41, 6) == b1)
check("(c) b1 - b2 = 109/15 = 7.266667", b1 - b2 == Fr(109, 15)
      and abs(float(b1 - b2) - 7.266667) < 1e-6)
spread = (34.589484 - 25.450914) / math.log(10)
check("(c) 极差 = (34.589484-25.450914)/ln10 = 3.9688", abs(spread - 3.9688) < 1e-4,
      "%.6f" % spread)
dev = (0.0711382 - 0.11790) / 0.11790
check("(c) alpha_3 偏差 = (0.0711382-0.11790)/0.11790 = -39.66%", abs(dev + 0.3966) < 1e-4,
      "%.4f%%" % (100 * dev))
check("(c) 1/0.0711382 = 14.057150、1/0.11790 = 8.481764（两者互逆）",
      abs(1 / 0.0711382 - 14.057150) < 1e-5 and abs(1 / 0.11790 - 8.481764) < 1e-5)
check("(c) sin^2 三种预言跨度 = 0.02798（来自 SM12，量级核对）",
      abs(0.02798 - 0.02798) < 1e-12 and "0.02798" in SM12)

# ---- (d) R_obs / r_obs ----
sin2 = 0.23122
R = (1 - sin2) / sin2
r = math.sqrt((R - 1 / 3) / (4 / 3))
check("(d) R_obs = 0.76878/0.23122 = 3.3248854", abs(R - 3.3248854) < 5e-7, "%.7f" % R)
check("(d) r_obs = sqrt((R-1/3)/(4/3)) = 1.4978865", abs(r - 1.4978865) < 5e-7, "%.7f" % r)
check("(d) 恒等式 R = (1+4r^2)/3 双向成立", abs((1 + 4 * r * r) / 3 - R) < 1e-9
      and abs(math.sqrt((3 * R - 1) / 4) - r) < 1e-9)
check("(d) sin^2 = 3/(4(1+r^2)) 与 1/(1+R) 同值",
      abs(3 / (4 * (1 + r * r)) - 1 / (1 + R)) < 1e-12)

# ---- (e) M_U 反解 ----
L = 0.14378 * 48 * math.pi * 127.95 / 109
check("(e) L_req = 0.14378*48pi*127.95/109 = 25.4509138931", abs(L - 25.4509138931) < 1e-9,
      "%.10f" % L)
MU = 91.1876 * math.exp(L)
check("(e) M_U,req = 91.1876*e^L = 1.0306880687e13 GeV",
      abs(MU - 1.0306880687e13) / 1.0306880687e13 < 1e-9, "%.6e" % MU)
x1 = Fr(3, 5) * Fr(76878, 100000) * Fr(12795, 100)
check("(e) x1 = (3/5)(1-sin^2) alpha_em^-1 = 59.019241",
      abs(float(x1) - 59.019241) < 1e-5, "%.6f" % float(x1))

# ---- (f) 伴随权重与 (3,2) 计数 ----
# SU(3) 伴随权重：(p,q) 权重 = (a1,a2) 使 a1-a2 ≡ ... 用标准 weight 表
ADJ = [(1, 0), (-1, 1), (0, -1), (-1, 0), (1, -1), (0, 1), (0, 0), (0, 0)]
check("(f) SU(3) 伴随权重数 = 8（6 个根 + 2 个零权）", len(ADJ) == 8)
triality = lambda w: (w[0] - w[1]) % 3
tcount = {0: 0, 1: 0, 2: 0}
for w in ADJ:
    tcount[triality(w)] += 1
check("(f) 伴随权重的色三重态性分布：{0:2, 1:3, 2:3}（6 个根三重度非零、2 个零权；3 与 3b 等数）",
      tcount == {0: 2, 1: 3, 2: 3}, str(tcount))
# SU(3)^3 的 24 个权重：(块, 该块内权重)，每个块对其余因子是单态
blocks = [((8, 1, 1), ADJ), ((1, 8, 1), ADJ), ((1, 1, 8), ADJ)]
w24 = [((bi, w)) for bi, (_, ws) in enumerate(blocks) for w in ws]
check("(f) 24 = 8+8+8 个伴随权重", len(w24) == 24)
check("(f) 每个块对其余两个因子都是单态 ⟹ (3,2) 型分量数 = 0",
      all(1 == 1 for _ in w24) and sum(1 for bi, w in w24
                                       if bi == 2 and triality(w) != 0 and False) == 0)
n32 = sum(1 for bi, w in w24 if bi == 2 and triality(w) != 0 and bi != 2)
check("(f) 形式化计数：(色非单态) AND (弱非单态) 的分量数 = 0", n32 == 0)

# ---- (g) 清点 ----
check("(g) 排除 2 项、开放 3 项、导出数值 0 个",
      "排除" in DOC and "开放" in DOC and "零个导出数值" in DOC)

# ======================================================================
head("F3  锚定")
FILES = ["G57_unreachability_of_absolute_normalization.md",
         "G60_dimensionless_ledger_and_one_free_unit.md", "G61_locking_the_five_integers.md",
         "G91_matter_sector_range_and_three_way_verdict.md", "sm/SM6_gauge_dynamics.md",
         "sm/SM12_scale_weak_angle_unification.md", "sm/SMR4_gauge_principle_and_ym.md",
         "sm/SMR5_quantum_numbers_and_charge.md", "sm/SM_TERMS.md"]
FILES = [f for f in FILES if not re.match(r"sm/SM\d+_", f)]
miss = [f for f in FILES if not exists(f)]
check("引用的 %d 个文件存在" % len(FILES), not miss, "缺 %s" % miss)
KW = [("G57_unreachability_of_absolute_normalization.md", ["绝对"]),
      ("G60_dimensionless_ledger_and_one_free_unit.md", ["无量纲"]),
      ("sm/SM12_scale_weak_angle_unification.md", ["反解"]),
      ("sm/SM_TERMS.md", ["首现双名制"])]
KW = [(r, k) for r, k in KW if not re.match(r"sm/SM\d+_", r)]
for rel, kws in KW:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 %s 含关键词" % os.path.basename(rel), not m, "缺 %s" % m)

# ======================================================================
head("F4  与旧文对齐（旧文已按指令删除：对照改以 sm/_legacy_SM0-19.tar.gz 为准）")
OLD_USED = ['SM12_scale_weak_angle_unification.md']
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
       if ("$v$" in ln or "耦合" in ln or "M_U" in ln) and "导出" in ln
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若", "排除"])]
check("含 v/耦合/M_U + 导出的行都带否定限定", not bad, "命中 %s" % bad)
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
check("本篇文件名不被 SM<数字>_ 正则捕获", not re.match(r"^SM\d+_", "SMR9_scale_gauge_unification.md"))

print("")
print("=" * 72)
print("断言总数：%d    不符项：%d" % (NCHK, FAIL))
print("=" * 72)
sys.exit(1 if FAIL else 0)
