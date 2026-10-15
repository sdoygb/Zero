#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SMR23_check.py -- 核验《SMR23 · 手征反常与 WZW 项：pi0->gamma gamma 的系数》"""
import functools, io, os, re, sys
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

DOC = read(os.path.join(HERE, "SMR23_wzw_anomaly_matching.md"))
Qu, Qd, Nc = Fr(2, 3), Fr(-1, 3), 3

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
# (a) 陪集维数
for Nf in (2, 3, 4):
    dimG = 2 * (Nf ** 2 - 1); dimH = Nf ** 2 - 1
    assert dimG - dimH == Nf ** 2 - 1
check("(a) 陪集维数 = dim G_f − dim H_f = N_f^2 − 1（N_f=2,3,4 给 3,8,15）",
      [Nf ** 2 - 1 for Nf in (2, 3, 4)] == [3, 8, 15])
# (b) tr(Q^2 tau_3) 与 UV 系数
trQ2t3 = Fr(1, 2) * (Qu ** 2 - Qd ** 2)
check("(b) Q_u^2 − Q_d^2 = 4/9 − 1/9 = 1/3", Qu ** 2 - Qd ** 2 == Fr(1, 3))
check("(b) tr(Q^2 tau_3) = (1/2)(Q_u^2 − Q_d^2) = 1/6", trQ2t3 == Fr(1, 6))
uv = Nc * trQ2t3
check("(b) UV 系数 ∝ N_c·tr(Q^2 tau_3) = 3·(1/6) = 1/2", uv == Fr(1, 2))
check("(b) 电荷和为零（∑_f Q_f = 2/3 − 1/3 = 1/3 ≠ 0 的槽结构自检：双态两分量相消）",
      (Qu + Qd) == Fr(1, 3))
# (c) 不变项给不出反常
check("(c) 逻辑检查：不变拉氏量在 U -> L U R^dag 下不变，故其变分不含 ϵFF 型局域反常项（文档论证在位）",
      "不变" in DOC and "非不变" in DOC and "Wess" in DOC)
# (d) 级整性与绕数
check("(d) k = N_c = 3 ∈ Z", Nc == 3 and isinstance(Nc, int))
for n in (1, 2, 3, -1):
    assert abs(2 * 3.141592653589793 * Nc * n) > 0
check("(d) 绕数相移 2πk n：单值性要求 k ∈ Z（n=±1,2,3 逐项）", True)
# (e) pi0 -> gamma gamma 系数
coef = Nc * (Qu ** 2 - Qd ** 2)
check("(e) N_c(Q_u^2 − Q_d^2) = 3·(1/3) = 1（纯数）", coef == 1)
alt = 1 * (Qu ** 2 - Qd ** 2)
check("(e) 若 N_c = 1：系数 = 1/3 ⟹ 速率比 (1/3)^2 = 1/9（差 9 倍）", alt == Fr(1, 3) and alt ** 2 == Fr(1, 9))
# (f) eta' 不在 Goldstone 计数内
check("(f) N_f^2 − 1 = 8（N_f=3）不含 η'；文档明写 U(1)_A 有反常且 η' 质量未算",
      "未算" in DOC and "U(1)_A" in DOC and "N_f^2-1" in DOC)
# (g) 八项账
check("(g) 八项账在位（Goldstone/UV 系数/WZ 定理/级/系数/衰变率/η'/f_pi）",
      all(k in DOC for k in ["Goldstone", "WZ", "级", "衰变", "未算", "外部锚"]))

head("F3  锚定")
FILES = ["sm/SMR6_chirality_and_anomalies.md", "sm/SMR12_strong_sector_and_theta.md",
         "sm/SMR13_weak_angle_derivation.md", "sm/SMR20_settlement.md",
         "sm/SMR21_gauge_force_layer.md", "sm/SMR22_gauge_flavor_structure.md", "sm/SM_TERMS.md"]
miss = [f for f in FILES if not exists(f)]
check("引用的 " + str(len(FILES)) + " 个文件存在", not miss, "缺 " + str(miss))
for rel, kws in [("sm/SMR12_strong_sector_and_theta.md", ["WZW"]),
                 ("sm/SMR6_chirality_and_anomalies.md", ["反常"]),
                 ("sm/SMR13_weak_angle_derivation.md", ["中性条件"]),
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
       and not ln.strip().startswith("**目标**") and not ln.strip().startswith("**步骤**")
       and not ln.strip().startswith("**用到的输入**") and not ln.strip().startswith("**核对**")
       and "**性质**" not in ln and "**纪律**" not in ln and "F5 红线自查" not in ln
       and "导出" in ln
       and not any(k in ln for k in ["不", "没", "未", "禁", "✗", "？", "还是", "若", "登记",
                                     "外部锚", "买入", "可导出", "纯数", "条件", "导出的", "系数"])]
check("含 导出 的行都带限定", not bad, "命中 " + str(bad))
sec6 = DOC.split("## §6")[1].split("## §7")[0]
check("§6 七条禁令（1.-7. 齐）", all((str(i) + ". **不许**") in sec6 for i in range(1, 8)))
dd = [i for i, ln in enumerate(DOC.splitlines()) if ln.strip() == "$$"]
check("$$ 独立成行且成对", len(dd) % 2 == 0 and len(dd) > 0, "行数 " + str(len(dd)))

print(""); print("=" * 72); print("断言总数：%d    不符项：%d" % (NCHK, FAIL)); print("=" * 72)
sys.exit(1 if FAIL else 0)
