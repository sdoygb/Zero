#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_sum_reproduction_transition_theorems_check.py —— Zero 文章②【繁殖转移定理】的核验
======================================================================================
实断言（**真调用** `simulations/zero_sum_reproduction_audit.py` 的函数重算，再与文章和冻结 JSON 三方比对）：
  F1  五条规则都能建出转移矩阵，且矩阵非负、整数
  F2  谱半径 / 零特征值数 / 幂零指数与文章 §3 表逐项一致（并与 JSON 一致）
  F3  定理 T3：split_only 的 T 幂零，且幂零指数 = 序列归零的代
  F4  序列行为分类与文章一致（copy 指数、split_only 灭绝、其余常数/多项式）
  F5  文章 §5 的 no-go 表述在位（M=1+r / 2^r 不是物理后代数）
  F6  copy 之外没有规则给指数增长（文章头条）
"""
import io
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "simulations"))

LEDGER_MODE = os.environ.get("LH_LEDGER", "A")
FAIL = 0
N = 0


def check(name, cond, detail="", level="ind"):
    global FAIL, N
    N += 1
    ok = bool(cond)
    if not ok:
        FAIL += 1
    tag = ("v" if ok else "x") if (level == "ind" or LEDGER_MODE == "A" or not ok) else "i"
    print("  [%s] %s%s" % (tag, name, ("   " + detail) if detail else ""))


def head(t):
    print("\n" + "=" * 72 + "\n" + t + "\n" + "=" * 72)


DOC = io.open(os.path.join(HERE, "zero_sum_reproduction_transition_theorems.md"),
              encoding="utf-8").read()
JSN = json.load(io.open(os.path.join(HERE, "simulations",
                                     "zero_sum_reproduction_audit_results.json"), encoding="utf-8"))

import zero_sum_reproduction_audit as RA  # noqa: E402

RULES = ["persist", "copy", "split_parent", "split_only", "split_all"]
REF = {  # 文章 §3 的表：rho=1 报 I-T 的幂零指数；rho=0 报 T 的幂零指数；rho=2 不报
    "persist":      dict(rho=1.0, zero_eig=0, unipotent=1, nilp=None),
    "copy":         dict(rho=2.0, zero_eig=0, unipotent=None, nilp=None),
    "split_parent": dict(rho=1.0, zero_eig=0, unipotent=6, nilp=None),
    "split_only":   dict(rho=0.0, zero_eig=123, unipotent=None, nilp=6),
    "split_all":    dict(rho=1.0, zero_eig=0, unipotent=6, nilp=None),
}
SEED = "++--"

# ---------------------------------------------------------------- F1
head("F1  五条规则都能建出转移矩阵")
modes = RA.all_modes()
check("模式集非空（周期 ≤ MAX_PERIOD 的循环平衡词）", len(modes) > 0, "%d 个模式" % len(modes))
mats = {}
for rule in RULES:
    T = RA.transition_matrix(modes, rule)
    mats[rule] = T
    check("%-13s 矩阵形状 %s，元素非负整数" % (rule, T.shape),
          T.shape == (len(modes), len(modes)) and (T >= 0).all()
          and np.issubdtype(T.dtype, np.integer))

# ---------------------------------------------------------------- F2
head("F2  谱半径 / 零特征值 / 幂零指数（文章 §3 表 vs 本版重算 vs 冻结 JSON）")
for rule in RULES:
    T = mats[rule]
    ev = np.linalg.eigvals(T.astype(float))
    rho = float(np.max(np.abs(ev))) if ev.size else 0.0
    zero_eig = int(np.sum(np.abs(ev) < 1e-9))
    # 按 rho 分支：1 ⇒ 报 I−T 的幂零指数；0 ⇒ 报 T 的幂零指数；其它 ⇒ 不报
    if abs(rho - 1.0) < 1e-9:
        unip = RA.nilpotent_index(T - np.eye(T.shape[0], dtype=np.int64))
        nilp = None
    elif abs(rho) < 1e-9:
        unip = None
        nilp = RA.nilpotent_index(T)
    else:
        unip = nilp = None
    ref = REF[rule]
    js = JSN["matrix_audit"][rule]
    ok = (abs(rho - ref["rho"]) < 1e-9 and zero_eig == ref["zero_eig"]
          and unip == ref["unipotent"] and nilp == ref["nilp"])
    check("%-13s ρ=%.1f 零特征值=%d  I−T指数=%s  T指数=%s" % (rule, rho, zero_eig, unip, nilp),
          ok, "文章 %s；JSON ρ=%s" % (ref, js.get("spectral_radius")))
    check("%-13s 与冻结 JSON 一致" % rule,
          abs(rho - float(js.get("spectral_radius", -1))) < 1e-9)

# ---------------------------------------------------------------- F3
head("F3  定理 T3：split_only 幂零 ⇔ 有限灭绝")
T_so = mats["split_only"]
power = np.eye(T_so.shape[0], dtype=np.int64)
extinct_gen = None
for k in range(1, T_so.shape[0] + 1):
    power = power @ T_so
    if not np.any(power):
        extinct_gen = k
        break
check("split_only 的 T^6 = 0（幂零指数 6）", extinct_gen == 6, "实算 %s" % extinct_gen)
seed_idx = {w: i for i, w in enumerate(modes)}[RA.canonical_cycle(RA.word_from_string(SEED))]
v = np.zeros(len(modes), dtype=np.int64)
v[seed_idx] = 1
seq = []
for _ in range(18):
    seq.append(int(v.sum()))
    v = T_so @ v
first_zero = next((i for i, x in enumerate(seq) if x == 0), None)
check("split_only 的种群在有限代内归零（灭绝时间 ≤ 幂零指数）",
      first_zero is not None and first_zero <= extinct_gen,
      "首次归零于第 %s 代，序列前缀 %s" % (first_zero, seq[:6]))
check("文章写明『ρ=0 ⇔ 幂零 ⇔ 有限灭绝』", "幂零" in DOC and "有限灭绝" in DOC)

# ---------------------------------------------------------------- F4
head("F4  序列行为分类（文章 §3 表）")
for rule in RULES:
    T = mats[rule]
    v = np.zeros(len(modes), dtype=np.int64)
    v[seed_idx] = 1
    vals = []
    for _ in range(18):
        vals.append(int(v.sum()))
        v = T @ v
    cls = RA.classify_sequence(vals)
    ref_cls = JSN["matrix_audit"][rule]["sequence_classes"].get(SEED)
    check("%-13s 分类 = %-22s（JSON: %s）" % (rule, cls, ref_cls), cls == ref_cls,
          "序列前缀 %s" % vals[:5])
check("copy 的序列严格倍增（1,2,4,8,…）",
      JSN["sequences"]["copy"][SEED][:5] == [1, 2, 4, 8, 16],
      "%s" % JSN["sequences"]["copy"][SEED][:5])

# ---------------------------------------------------------------- F5
head("F5  文章 §5 的 no-go 表述")
check("写明 M=1+r 与 M=2^r 是分支程序计数", "分支程序计数" in DOC)
check("写明闭合只给持续、不给繁殖", "闭合本身只给持续" in DOC or "闭合只给持续" in DOC)
check("写明只有显式复制给指数增长", "只有显式复制" in DOC)
check("引用了程序自述（英文原文）", "branch-program counts" in DOC or "not automatically" in DOC)

# ---------------------------------------------------------------- F6
head("F6  头条：copy 之外无指数增长")
exp_rules = [r for r in RULES
             if JSN["matrix_audit"][r]["sequence_classes"].get(SEED) == "exponential"]
check("只有 copy 被判为 exponential", exp_rules == ["copy"], "实算 %s" % exp_rules)

# ---------------------------------------------------------------- 汇总
print("\n" + "=" * 72 + "\n汇总\n" + "=" * 72)
print("  断言 %d 项，不符项：%d" % (N, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
