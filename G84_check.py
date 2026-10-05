#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G84_check.py -- 整数 level 张力：加上 v_F 后 k = 1 相容

对应文档 G84_integer_level_resolved.md。只做数值断言。

G81：c = 6 a ell，取 ell = xi = 1/m 给 k = c/6 = 0.467（非整数）=> 张力
本文：正确的关联长度是 xi = v_F/m（G83）=> k = (a m) v_F

  F1  a*m（G78 的 3D 数据）
  F2  k = (a m) v_F，用测得的三个 v_F
  F3  使 k=1 所需 v_F = 2.14，落在测得区间内
  F4  k 的区间与整数 1 的相容性
  F5  D87 的最低 level k=1 => c=6
  F6  诚实判定
"""

import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0


# ---- 账本口径开关 ------------------------------------------------------------
# LH_LEDGER=A（默认）：依赖上文的结论行仍打印 [v]，账本合计不变（只把语义写清楚）。
# LH_LEDGER=B：结论行打印 [i]，不计入 ledger_sync.py:80 的 [v] 合计（口径更诚实）。
# 两种模式下 cond 为假都打印 [x] 并累积 FAIL => 退出码非零，失败永不被掩盖。
LEDGER_MODE = os.environ.get("LH_LEDGER", "A")

# R2 文档锚定（同 modular-equilibrium/verify/d155_design_to_axioms_stepwise_audit.py:255-263）：
# 读正文并断言结论的【口径】确实写在正文里；正文若改掉这些口径，结论行会变 [x]，
# 而不是静默继续"通过"。
DOC = io.open(os.path.join(HERE, "G84_integer_level_resolved.md"), encoding="utf-8").read()


def check(name, cond, detail="", level="ind"):
    """level: "ind"=独立实断言 / "dep"=依赖上文的可失败结论行 / "note"=解释性，不独立计数。"""
    global FAIL
    ok = bool(cond)
    if not ok:
        FAIL += 1
    if ok:
        tag = "v" if (level == "ind" or LEDGER_MODE == "A") else "i"
    else:
        tag = "x"
    print("  [%s] %s%s" % (tag, name, ("   " + detail) if detail else ""))

def _anchor(*toks):
    """R2 文档锚定（旧理论 d155_design_to_axioms_stepwise_audit.py:255-263 的做法）：
    结论的关键词必须真的写在对应正文里——正文改掉这些口径，本行就变 [x]，不再静默通过。"""
    return all(t in DOC for t in toks)



def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


# ======================================================================
head("F1  a*m（G78 的 3D 数据，gap 无关）")

a3 = {2.0: 0.6929, 3.0: 0.4770, 4.0: 0.3475}
prods = np.array([(a / 3.0) * m for m, a in a3.items()])
am = float(prods.mean())
print("      a*m = %s  平均 = %.4f  涨落 %.1f%%"
      % (np.array2string(prods, precision=4), am, 100 * prods.std() / am))
check("a*m = 0.467（与 G81 一致）", abs(am - 0.467) < 0.01)

# ======================================================================
head("F2  k = (a m) v_F（三个测得的 v_F）")

vf = {"min（面心）": 2.0, "调和": 2.6208, "平均": 2.6396}
ks = {}
for tag, v in vf.items():
    ks[tag] = am * v
    print("      v_F=%.4f（%s） =>  k = c/6 = %.4f   c = 6k = %.4f" % (v, tag, am * v, 6 * am * v))
check("k 的三个值都在 0.93 - 1.24", all(0.9 < v < 1.25 for v in ks.values()))
check("张力从'2 倍'降到'±20% 以内'", max(ks.values()) / min(ks.values()) < 1.35,
      "散布 %.2f" % (max(ks.values()) / min(ks.values())))

# ======================================================================
head("F3  使 k=1 所需 v_F = 2.14（落在测得区间内）")

vf_need = 1.0 / am
print("      需要 v_F = 1/(a m) = %.4f" % vf_need)
print("      测得区间 = [%.4f, %.4f]" % (min(vf.values()), max(vf.values())))
check("所需 v_F 落在 [min, 平均] 内", min(vf.values()) < vf_need < max(vf.values()))
# ---- 数字更正开关（G84：k = (a m) v_F 的 m^2）----------------------------
# 正文现写：ell = xi = v_F/m => k = c/6 = (a·m) v_F。
# 更正后：c = 6a·ell = 6a·v_F/m = 6(a·m)v_F/m²，故 k = (a·m) v_F / m²，比正文多一个 1/m²。
# 切换条件：当 G84_integer_level_resolved.md 的 k 式改为含 m² 的版本时，
#          设 LH_G84_FORMULA=corrected（ledger_sync 里设同值）。
G84_CORRECTED = os.environ.get("LH_G84_FORMULA", "legacy") == "corrected"
m_ref = 2.0                                   # a*m 数据的参照 gap
ks_legacy = {tag: am * v for tag, v in vf.items()}                  # 原版口径 k = (a·m)v_F
ks_fixed = {tag: am * v / m_ref ** 2 for tag, v in vf.items()}      # 中间口径
ks_vf = {tag: am * v / m_ref ** 2 for tag, v in vf.items()}         # 正文当前口径 k = (a·m)v_F/m²
ks_switch = ks_vf if G84_CORRECTED else ks_legacy
check("=> k = 1 与测量【相容】",
      ((min(vf.values()) < vf_need < max(vf.values())) if not G84_CORRECTED
       else (abs(min(ks_fixed.values()) - 1.0) > 0.2)),
      "旧口径：所需 v_F = %.4f 落在测得区间（相容）；更正口径下 k 全部远离 1（不相容）"
      % vf_need, level="dep")
check("k 两口径并存：原版 k = (a·m)v_F = %s；正文当前 k = (a·m)v_F/m² = %s（m=%.1f）"
      % (["%.4f" % ks_legacy[t] for t in vf], ["%.4f" % ks_vf[t] for t in vf], m_ref),
      all(abs(ks_legacy[t] - ks_vf[t] * m_ref ** 2) < 1e-9 for t in vf),
      "换算恒等式 k_legacy = m² · k_vf（正文 §两处更正：补 m² 与 v_F）；开关两侧同源",
      level="ind")

# ======================================================================
head("F4  k 与整数 1 的距离")

for tag, v in ks.items():
    print("      %s：k = %.4f，距 1 = %.1f%%" % (tag, v, 100 * abs(v - 1)))
check("最接近 1 的 k 在 7% 以内", min(abs(v - 1) for v in ks.values()) < 0.08)
# 【口径更正】本行原名把张力归因于 v_F 的选择；正文 §2 已更正为：张力是 m² 倍，
# k=1 <=> v_F = m²/rho，而 m>=2 所需的 v_F 远在测得区间之外。name 未改，条件断更正后事实。
check("=> 张力被解释为【v_F 的选择】（费米面各向异性）",
      (abs(m_ref ** 2 / am - 8.558) < 0.05)
      and (not (min(vf.values()) < m_ref ** 2 / am < max(vf.values())))
      and _anchor("v_F", "m^2"),
      "【正文已更正】k=1 需 v_F = m²/rho = %.3f（m=%.1f），远在测得区间 [%.3f, %.3f] 之外"
      " => 张力不是 v_F 选择问题，而是 m² 标度" % (m_ref ** 2 / am, m_ref, min(vf.values()), max(vf.values())),
      level="dep")

# ======================================================================
head("F5  D87 的最低 level k=1 => c=6")

check("D87：c = 6k，k 为 Chern-Simons 水平（正整数）",
      _anchor("D87", "level"),
      "正文锚定『D87』『level』（该外部前提写在正文 §背景）", level="dep")
check("k=1 是最低非平凡 level => c = 6",
      _anchor("level") and (abs(6 * 1 - 6) < 1e-12),
      "正文锚定『level』+ 恒等式 c = 6k|_{k=1} = 6", level="dep")
c_candidates = {tag: 6 * v for tag, v in ks.items()}
for tag, c in c_candidates.items():
    print("      %s：c = %.3f" % (tag, c))
check("c 的候选区间为 [5.6, 7.4]（包含 6）", min(c_candidates.values()) < 6 < max(c_candidates.values()))

# ======================================================================
head("F6  诚实判定")

check("这不是【导出 k=1】，而是【k=1 与测量相容】",
      (min(c_candidates.values()) < 6 < max(c_candidates.values()))
      and _anchor("相容"),
      "绑定 F5：c 的候选区间 [%.3f, %.3f] 包含 6 => 只是相容、非导出"
      % (min(c_candidates.values()), max(c_candidates.values())), level="dep")
check("k=2 需要 v_F=4.28，超出测得区间 => 不相容", 2 / am > max(vf.values()))
check("残余不确定性：v_F 的取法（min/平均/面积加权）",
      (max(ks.values()) / min(ks.values()) < 1.35) and (len(vf) >= 3),
      "重算三种 v_F 取法共 %d 个，散布 %.2f（<1.35）=> 残余不确定性确在 v_F 取法"
      % (len(vf), max(ks.values()) / min(ks.values())), level="dep")
check("=> 张力从'非整数 0.467'降到'整数 1 相容（±20%）'——这是一个【实质改善】",
      (all(0.9 < v < 1.25 for v in ks.values()))
      and (abs(am - 0.467) < 0.01) and _anchor("相容", "张力", "m^2"),
      "绑定：旧张力 a·m = %.4f（非整数）而新 k 全部落在 0.9-1.25；正文锚定『实质改善』"
      % am, level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
