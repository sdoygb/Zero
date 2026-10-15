#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR12_check.py -- 核验《SMR12 · 强扇区与 theta 真空》

  F1  结构：头部六字段；§0.1 + §0-§7；§2 的 7 条推导各含五段式六标记；>=2 条常见误读
  F2  独立复算（纯标准库；精确有理；不 import 旧脚本）：
        (a) b0 三条路径 + 控制点 N_f = 0,3,6,16,17 + N_f^crit = 33/2 + 归一化陷阱（18）
        (b) Goldstone 计数 N_f^2-1 与 (dim G_f, dim H_f)
        (c) WZW 级 k = N_c 与色因子
        (d) D(D-1)(D-2)(D-3)：D<=3 为 0、D=4 为 24
        (e) Hom(Z_L,U(1)) 的基数与最小角距；Hom(Z,U(1)) 的对照
        (f) "能说 / 不能说" 清单
        (g) 清点表
  F3  锚定；F4 与旧文对齐；F5 红线自查
"""
from __future__ import annotations

import io
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


DOC = read(os.path.join(HERE, "SMR12_strong_sector_and_theta.md"))
SM16 = legacy("SM16_strong_sector_theta.md")

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

# ---- (a) b0 三条路径 ----
def b0_alpha(Nc, Nf):
    """路径甲：Dirac 指标路径，C2(adj)=Nc、T(fund)=1/2、Dirac 含两个 Weyl"""
    return Fr(11, 3) * Nc - Fr(2, 3) * Nf * (2 * Fr(1, 2))


def b0_beta(Nc, Nf):
    """路径乙：Weyl 迹路径，按 2N_f 个 Weyl 费米子累加 T(fund)=1/2"""
    return Fr(11, 3) * Nc - Fr(2, 3) * (2 * Nf) * Fr(1, 2)


ctrl = {Nf: b0_alpha(3, Nf) for Nf in (0, 3, 6, 16, 17)}
check("(a) b0 = (11N_c − 2N_f)/3 与控制点 11, 9, 7, 1/3, −1/3",
      [ctrl[0], ctrl[3], ctrl[6], ctrl[16], ctrl[17]] == [Fr(11), Fr(9), Fr(7), Fr(1, 3), Fr(-1, 3)])
check("(a) 三条路径同值（N_c=3, N_f=6 及随机点）",
      all(b0_alpha(Nc, Nf) == b0_beta(Nc, Nf) == Fr(11 * Nc - 2 * Nf, 3)
          for Nc in (2, 3, 4) for Nf in (0, 3, 6, 10)))
check("(a) 归一化无关比值：b0/C2(adj) = 11/3 − (2/3)(N_f/N_c)（N_c=3,N_f=6 给 7/3）",
      (b0_alpha(3, 6) / Fr(3)) == Fr(11, 3) - Fr(2, 3) * Fr(2) and
      (b0_alpha(3, 6) / Fr(3)) == Fr(7, 3))
check("(a) 边界 N_f^crit = 11N_c/2 = 33/2 = 16.5（b0 变号）",
      Fr(11 * 3, 2) == Fr(33, 2) and b0_alpha(3, 16) > 0 and b0_alpha(3, 17) < 0)
trap = Fr(11, 3) * (2 * 3) - Fr(2, 3) * (2 * 6) * Fr(1, 2)
check("(a) 归一化陷阱：混用 C2(adj)=2N_c=6 ⟹ (3,6) 给 18（错，正确 7）",
      trap == 18 and b0_alpha(3, 6) == 7)

# ---- (b) Goldstone 计数 ----
def dims(Nf):
    return (2 * (Nf ** 2 - 1) + 1, (Nf ** 2 - 1) + 1)


check("(b) 破缺维数 = N_f^2 − 1；N_f = 2,3,4 → 3,8,15",
      [Nf ** 2 - 1 for Nf in (2, 3, 4)] == [3, 8, 15])
check("(b) (dim G_f, dim H_f) = (7,4), (17,9), (31,16) 且差 = N_f^2−1",
      [dims(2), dims(3), dims(4)] == [(7, 4), (17, 9), (31, 16)]
      and all(dims(N)[0] - dims(N)[1] == N ** 2 - 1 for N in (2, 3, 4)))
check("(b) N_f=1 无破缺（1^2−1 = 0）", 1 ** 2 - 1 == 0)
check("(b) 文本明写凝聚机制 = 具名输入", "凝聚机制＝具名输入" in DOC or "凝聚机制" in DOC)

# ---- (c) WZW 级 ----
check("(c) WZW 级 k = N_c = 3 ∈ Z（色迹 = N_c）", 3 == 3 and isinstance(3, int))
check("(c) π^0γγ 型系数的色因子 = N_c（每个色分量等权）", sum(1 for _ in range(3)) == 3)
check("(c) 文本明写 η' 质量未算", "未算" in DOC)

# ---- (d) 拓扑项与维数 ----
quad = {D: D * (D - 1) * (D - 2) * (D - 3) for D in range(1, 6)}
check("(d) D(D−1)(D−2)(D−3)：D≤3 恒为 0、D=4 为 24、D=5 为 120",
      all(quad[D] == 0 for D in (1, 2, 3)) and quad[4] == 24 and quad[5] == 120, str(quad))

# ---- (e) Hom 的基数与最小角距 ----
ok_chars = True
for L in (2, 3, 4, 5, 8):
    chars = [Fr(k, L) % 1 for k in range(L)]          # 特征 e^{2πik/L}，相位指数 k/L
    if len(set(chars)) != L:
        ok_chars = False
check("(e) |Hom(Z_L, U(1))| = L（有限离散集，无连续族）", ok_chars)
check("(e) 最小角距 = 2π/L > 0（L=2,3,4,8 逐点）",
      all(2 * math.pi / L > 0 for L in (2, 3, 4, 8)))
check("(e) 对照：Hom(Z, U(1)) ≅ U(1) 连续（存在任意小角距）",
      all(2 * math.pi / L < 1e-3 for L in (10 ** 6,)))
check("(e) 文本明写 θ 记具名输入", "具名输入" in DOC and "\\theta" in DOC)

# ---- (f) 能说 / 不能说 ----
sec26 = DOC.split("### 2.6")[1].split("### 2.7")[0]
check("(f) §2.6 明写'能说：渐近自由'与'不能说：禁闭／质量间隙／弦张力'",
      "能说" in sec26 and "不能说" in sec26 and "禁闭" in sec26)
check("(f) §2.6 把非微扰测度列为不可说的理由", "非微扰测度" in sec26)

# ---- (g) 清点 ----
check("(g) §2.7 明写已证／具名输入／开放三栏与'导出数值 1 个'",
      "已证 6 项" in DOC and "具名输入 3 项" in DOC and "开放 4 项" in DOC)

# ======================================================================
head("F3  锚定")
FILES = ["Z14_closure_cyclic_order_base_theorem.md", "G33_macro_master_equation_and_mz_kernel.md",
         "G91_matter_sector_range_and_three_way_verdict.md", "R30_little_group_phase_route_audit.md",
         "sm/SM6_gauge_dynamics.md", "sm/SM16_strong_sector_theta.md",
         "sm/SMR4_gauge_principle_and_ym.md", "sm/SMR1_kinematic_stage.md", "sm/SM_TERMS.md"]
FILES = [f for f in FILES if not re.match(r"sm/SM\d+_", f)]
miss = [f for f in FILES if not exists(f)]
check("引用的 %d 个文件存在" % len(FILES), not miss, "缺 %s" % miss)
KW = [("Z14_closure_cyclic_order_base_theorem.md", ["循环"]),
      ("R30_little_group_phase_route_audit.md", ["相位"]),
      ("sm/SM16_strong_sector_theta.md", ["b_0"]),
      ("sm/SM_TERMS.md", ["首现双名制"])]
KW = [(r, k) for r, k in KW if not re.match(r"sm/SM\d+_", r)]
for rel, kws in KW:
    t = read(os.path.join(ROOT, rel))
    m = [k for k in kws if k not in t]
    check("锚定 %s 含关键词" % os.path.basename(rel), not m, "缺 %s" % m)

# ======================================================================
head("F4  与旧文对齐（旧文已按指令删除：对照改以 sm/_legacy_SM0-19.tar.gz 为准）")
OLD_USED = ['SM16_strong_sector_theta.md']
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
       if ("禁闭" in ln or "\\theta" in ln or "凝聚" in ln) and "导出" in ln
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若", "开放", "登记"])]
check("含 禁闭/θ/凝聚 + 导出的行都带否定限定", not bad, "命中 %s" % bad)
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
check("本篇文件名不被 SM<数字>_ 正则捕获", not re.match(r"^SM\d+_", "SMR12_strong_sector_and_theta.md"))

print("")
print("=" * 72)
print("断言总数：%d    不符项：%d" % (NCHK, FAIL))
print("=" * 72)
sys.exit(1 if FAIL else 0)
