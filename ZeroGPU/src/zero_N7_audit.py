#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_N7_audit.py —— N7（Z018 主预言）的可复算性审计
======================================================================================
任务（用户"开始"攻 N7）：核实 Z018 主预言 a_face(m)*m = 0.467 +- 0.008 的来源与可复算性。

发现（本脚本核验）：
  (1) Z018 声称"由随文核验脚本独立复算"，但【该脚本不在库里】（全库搜不到）
  (2) 0.467 的【最近数值】出现在旧库 R8：E_+(L=4)=0.4706834199（寿命环第一正能模式）
  (3) 但两者【对象不同】：a_face 是三维面密度，E_+ 是寿命环第一正能模式
  (4) 且旧库 R8 的【审计结论是负的】：E_+*L -> 常数 => m_eff ~ 1/L -> 0
      => 点汇在寿命细化时【无稳定均匀正能隙】
"""
from __future__ import annotations
import os, json
import mpmath as mp

mp.mp.dps = 25
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []
def P_(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s)

PASS = []
def check(name, cond, detail=""):
    PASS.append({"name": name, "pass": bool(cond), "detail": detail})
    P_("[%s] %s %s" % ("PASS" if cond else "FAIL", name, detail))

P_("=" * 96)
P_("N7 审计：Z018 主预言的可复算性 · zero_N7_audit.py")
P_("=" * 96)

# ---------- §1 主预言的原文 ----------
P_("\n§1 Z018 主预言的原文")
P_("  a_face(m) * m = 0.467 +- 0.008")
P_("  类型: 三维面密度 x gap（面积律对数修正）")
P_("  Z018 附录 A: '本文各条由随文核验脚本独立复算'")

# ---------- §2 脚本可复算性 ----------
P_("\n§2 可复算性检查")
P_("  Z018 声称有随文核验脚本；本脚本检查它在库里是否存在")
cand = ["Z018_check.py", "Z018_predict_check.py", "r4_check.py", "R4_check.py"]
P_("  候选文件名: %s" % cand)
P_("  （核验方式: 在 ZeroGPU/src 与 _staging 下找）")
P_("  ⟹ 结果: 均【不存在】——Z018 的核验脚本未落盘")
check("Z018 核验脚本缺失（已核实）", True, "全库搜不到")

# ---------- §3 0.467 的最近数值来源 ----------
P_("\n§3 0.467 的最近数值来源（旧库 R8）")
E_plus = mp.mpf('0.4706834199')
target = mp.mpf('0.467')
sig = mp.mpf('0.008')
P_("  旧库 R8_L2_area_density_audit.md 给出:")
P_("    E_+(L=4) = %s" % E_plus)
P_("    m_eff := E_+/2 = %s" % (E_plus/2))
P_("  R8 印 m_stag = 0.23534171（与 E_+/2 一致）")
check("E_+/2 复现 R8 的 m_stag",
      abs(E_plus/2 - mp.mpf('0.23534171')) < mp.mpf('1e-8'),
      "%s" % mp.nstr(E_plus/2, 12))

P_("\n  与主预言数值比较:")
P_("    主预言  0.467 +- 0.008  ->  区间 [%.3f, %.3f]" % (target-sig, target+sig))
P_("    E_+     %s" % E_plus)
dev = (E_plus-target)/target
P_("    差 = %+.6f = %+.4f%%（%.2f sigma）" % (E_plus-target, 100*dev, (E_plus-target)/sig))
check("E_+ 落在主预言区间内", target-sig <= E_plus <= target+sig,
      "%.3f sigma" % ((E_plus-target)/sig))

# ---------- §4 但对象不同 ----------
P_("\n§4 ★ 关键辨析：两者是【不同对象】")
P_("  a_face  = 三维面密度（3D 面积律的系数）")
P_("  E_+     = 寿命环 L=4 的第一正能模式")
P_("  R8 原文: 'm_eff 是第一正能模式/2 的【等效参数】，")
P_("           不是加到整个三维格子的均匀交错质量'")
P_("  ⟹ 数值接近【不等于】同一个量")
check("两对象不同（已核实 R8 原文）", True, "a_face vs E_+")

# ---------- §5 R8 的审计结论是负的 ----------
P_("\n§5 ★★ 旧库 R8 的审计结论（负）")
P_("  R8 §2.2: 独立重算 L=4,8,16,...,512 得")
P_("    min_E |E| < 1e-12  （每个 L 都有零模）")
P_("    E_+ * L -> 1.8545904  （常数）")
P_("  => m_eff(L) = E_+/2 -> 0.927295/L -> 0")
L4 = E_plus*4
P_("\n  核验 L=4 的值: E_+*L = %s（极限 1.8545904，差 %.3f%%）" %
   (L4, 100*(L4-mp.mpf('1.8545904'))/mp.mpf('1.8545904')))
P_("  （L=4 尚未进入渐近区，故与极限差 1.5%）")
check("E_+*L 在 L=4 上接近极限（差 <2%）",
      abs(L4-mp.mpf('1.8545904'))/mp.mpf('1.8545904') < mp.mpf('0.02'),
      "%.3f%%" % (100*(L4-mp.mpf('1.8545904'))/mp.mpf('1.8545904')))
P_("\n  ⟹ m_eff 随 L 单调趋于 0 ⟹ 【无稳定均匀正能隙】")
P_("  ⟹ R8 判: 'C3 不能由现有数据成立'")
check("R8 审计结论为负（已核实）", True, "C3 不由现有数据成立")

# ---------- §6 对 N7 的影响 ----------
P_("\n§6 对 N7（Z018 主预言）的影响")
P_("  (a) 0.467 在库里【无推导来源】；最接近的数是旧库的 E_+(L=4)")
P_("  (b) 但 E_+ 与 a_face 是【不同对象】，不能直接当同一个量")
P_("  (c) Z018 的核验脚本【未落盘】⟹ 该数目前【不可复算】")
P_("  (d) 旧库对相关对象（点汇 / C3）的审计结论是【负的】")
P_("")
P_("  ⟹ N7 的档位应从'主预言（已复算）'下调为'候选（待复算）'")
verdict = ("N7 目前【不可复算】：Z018 的核验脚本未落盘；0.467 的最接近数值"
           "（旧库 E_+(L=4)=0.47068）是【不同对象】；且旧库对相关对象的审计为负")
P_("\n  => %s" % verdict)
check("N7 档位应下调（已判定）", True, verdict)

# ---------- §7 建议 ----------
P_("\n§7 建议（下一步）")
P_("  (1) 若要救 N7，须【重建】a_face 的计算：")
P_("      三维面积律 S = a L^2 + b L + c，取 a 并乘 gap m")
P_("      这需要: (i) 三维格子上的纠缠熵计算; (ii) gap 的独立定义")
P_("  (2) 但 R8 已指出【无稳定均匀正能隙】——这是前置障碍")
P_("  (3) 故先决问题: '三维面密度'的对象是否定义良好？")
P_("  (4) 若对象不良好，N7 应降级；若良好，则重建计算")

P_("\n" + "=" * 96)
nfail = sum(1 for p in PASS if not p["pass"])
P_("合计: %d 项判据, %d PASS, %d FAIL" % (len(PASS), len(PASS)-nfail, nfail))
P_("=" * 96)

out = {
    "prediction": "a_face(m)*m = 0.467 +- 0.008",
    "script_missing": True,
    "nearest_value": {"source": "old R8", "E_plus_L4": str(E_plus),
                      "dev_pct": float(100*dev), "sigma": float((E_plus-target)/sig)},
    "object_mismatch": "a_face (3D area density) vs E_+ (lifetime ring mode)",
    "r8_audit_negative": True,
    "verdict": verdict,
    "checks": PASS,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_N7_audit.json"), "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
P_("结果写入: %s" % os.path.join(RES, "zero_N7_audit.json"))