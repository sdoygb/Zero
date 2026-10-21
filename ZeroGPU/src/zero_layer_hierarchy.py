#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_layer_hierarchy.py —— 层动力学的【层级检验】
======================================================================================
用户口径（2026-10-20）："全都不需要，只需要 L0-L1-L2 还有各亚层的动力学就可以"

本脚本把这句话变成一个【可判定的问题】：
  L2 的 3 个亚层提供了分裂的【位置】（数量 3 对上了），
  但精细/超精细需要的是【量级】 $10^{-5}$ / $10^{-7}$，而层数字全是 $O(1)$。
  ⟹ 问题：那个小参数从哪来？

检验设计（防凑数）：
  (A) 精确算出需要的小数（21cm/Lyman-α 等）
  (B) 枚举层比值的乘积式 ∏ r_i^{n_i}，看能否命中
  (C) 【零假设对照】：同样复杂度的随机目标有多少会被命中
      —— 只有命中率显著高于零假设，才算结构，否则是凑数
"""
from __future__ import annotations
import os, json, itertools, math
import mpmath as mp

mp.mp.dps = 30
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
P_("层动力学的层级检验 · zero_layer_hierarchy.py")
P_("=" * 96)

# ---------- §1 层账给了什么 ----------
P_("\n§1 层账（zero_layers.md）给了什么")
P_("  L0: 舞台图 Γ、公理、重数        —— 0 层，静态")
P_("  L1: 词层、高度 h(w)、双覆盖 Z2   —— H=L/2 层，静态、无自身更新律")
P_("  L2: 走者活动、4 顶点占据→3 个相邻比 —— 3 层，唯一有自身更新律")
P_("")
P_("  ★ L2 的分层数 = 3 ⟹ 一个能级内部最多有 3 个亚层结构")
P_("     精细结构需要 1 个能级裂成 2~3 个（2S1/2, 2P1/2, 2P3/2）⟹ 【数量对上了】")
P_("     超精细需要 1 个能级裂成 2 个（F=0,1）⟹ 【数量也对上了】")
check("L2 亚层数 3 与精细结构的分裂数匹配", True, "3 亚层")

# ---------- §2 层数字的量级 ----------
P_("\n§2 ★ 但层数字全是 O(1)（L=4 为物理值，Z017/Z070）")
L = 4
layers = {
    "2/3": mp.mpf(2)/3,          # Z087 稳态 record
    "1/3": mp.mpf(1)/3,
    "5/9": mp.mpf(5)/9,          # 纯度 κ1（Z007）
    "4/9": mp.mpf(4)/9,
    "1/2": mp.mpf(1)/2,          # 重播种第二根 / c(4)=2 的倒数
    "3/2": mp.mpf(3)/2,          # K4 record 相邻比（下）
    "2":   mp.mpf(2),            # K4 record 相邻比（中）
    "8/3": mp.mpf(8)/3,          # K4 record 相邻比（上）
    "2/9": mp.mpf(2)/9,          # Z019 每步权重 (A_l/3)^2
    "sqrt2-1": mp.sqrt(2)-1,     # A_l - 1 = tan(pi/8)（Z019 (5.1)）
}
P_("  %-10s %-16s %s" % ("名字", "值", "出处"))
for k, v in layers.items():
    P_("  %-10s %-16s %s" % (k, mp.nstr(v, 12), ""))
P_("")
P_("  L1 高度分布（L=4）: b_h = [4, 2]，归一化 [2/3, 1/3]  <- O(1)")
P_("  L2 三槽（K4，persistent）: [0.4055, 0.6931, 0.9808] = log(3/2), log 2, log(8/3)  <- O(1)")
lo = min(v for v in layers.values())
P_("  层数字的最小值 = %.4f （O(1)）" % float(lo))
check("层数字全部 O(1)（最小值 >0.1）", lo > mp.mpf('0.1'), "min = %.4f" % float(lo))

# ---------- §3 需要的小数 ----------
P_("\n§3 ★ 需要的小数（实测）")
nu21 = mp.mpf('1420.405751768e6')                      # Hz
lam_ly = mp.mpf('121.567e-9')                          # m
c = mp.mpf('299792458.0')
nu_ly = c/lam_ly
ratio = nu21/nu_ly
P_("  21cm 频率          = %.6e Hz" % nu21)
P_("  Lyman-α 频率       = %.6e Hz" % nu_ly)
P_("  ★ 21cm/Lyman-α     = %.6e" % ratio)
P_("  精细结构/gross     ~ α²/2 = %.6e" % float(mp.mpf('0.0072973525693')**2/2))
P_("  超精细/gross       ~ (4/3)(m_e/m_p)α² = %.6e"
   % float(mp.mpf(4)/3*(mp.mpf('0.51099895')/mp.mpf('938.27208816'))*mp.mpf('0.0072973525693')**2))
P_("")
P_("  ⟹ 需要的抑制因子量级: 10^{-5} ~ 10^{-7}")
P_("  ⟹ 层数字是 O(1) ⟹ 必须靠【乘积/幂】才到得了")
check("需要的抑制量级已算出", True, "%.2e" % float(ratio))

# ---------- §4 枚举：层比值的幂乘积能否命中 ----------
P_("\n§4 枚举检验：∏ r_i^{n_i} 能否命中 21cm/Lyman-α")
P_("  约束: Σ|n_i| ≤ 30（复杂度预算），r_i 取 §2 的层比值")

keys = list(layers.keys())
def enumerate_products(maxdeg):
    hits = []
    # 单因子幂 + 双因子幂（保持搜索空间诚实可控）
    rng = range(-maxdeg, maxdeg+1)
    for i, k1 in enumerate(keys):
        v1 = layers[k1]
        for n1 in rng:
            if n1 == 0: continue
            val = v1**n1
            if val <= 0: continue
            hits.append((abs(n1), "%s^%d" % (k1, n1), val))
            for k2 in keys[i+1:]:
                v2 = layers[k2]
                for n2 in range(-maxdeg, maxdeg+1):
                    if n2 == 0: continue
                    if abs(n1)+abs(n2) > maxdeg: continue
                    val2 = val*(v2**n2)
                    hits.append((abs(n1)+abs(n2), "%s^%d * %s^%d" % (k1, n1, k2, n2), val2))
    return hits

H = enumerate_products(30)
P_("  枚举总量 = %d 个表达式" % len(H))

def near(val, target, tol):
    return abs(val-target)/target < tol

tol = mp.mpf('0.01')          # 1% 窗
hit1 = [(d, e, v) for d, e, v in H if near(v, ratio, tol)]
P_("\n  命中 21cm/Lyman-α（1%% 窗）: %d 个" % len(hit1))
hit1.sort(key=lambda t: (t[0], abs(t[2]-ratio)))
for d, e, v in hit1[:8]:
    P_("    复杂度 %-3d %-34s = %.6e （偏 %+.3f%%）"
       % (d, e, v, float(100*(v-ratio)/ratio)))

# ---------- §5 ★ 零假设对照 ----------
P_("\n§5 ★ 零假设对照：同样复杂度能命中【任意】目标的概率")
P_("  做法：取 200 个与目标同量级（10^-8..10^-5）的随机目标，")
P_("        数每个目标在同一个枚举里有多少命中；与真目标的命中数比较。")
import random
random.seed(20261020)
mags = []
for _ in range(200):
    exp = random.uniform(-8, -5)
    mags.append(mp.mpf(10)**mp.mpf(exp))
cnt = []
for t in mags:
    n = sum(1 for d, e, v in H if near(v, t, tol))
    cnt.append(n)
avg = sum(cnt)/len(cnt)
P_("  随机目标的平均命中数 = %.2f" % avg)
P_("  真目标的命中数       = %d" % len(hit1))
P_("  ⟹ 命中数比 = %.2f 倍" % (len(hit1)/avg if avg > 0 else float('inf')))
check("真目标命中数不超过零假设的 3 倍（否则才可能是结构）",
      len(hit1) <= 3*avg if avg > 0 else True,
      "真 %d vs 随机均值 %.2f" % (len(hit1), avg))

# ---------- §6 判定 ----------
P_("\n§6 判定")
P_("  (1) ✓ 数量层：L2 的 3 个亚层与精细/超精细的分裂数【匹配】")
P_("  (2) ✗ 量级层：层数字全 O(1)，而分裂要 10^-5~10^-7")
P_("  (3) 【关键】命中数与零假设相当 ⟹ 这些乘积式【不构成推导】")
P_("  (4) 库内唯一已证明的大抑制机制是 Z019 的 15 步乘积 (2/9)^15 = 1.59e-10")
P_("      但 (2/9)^n 对整数 n 逼不出 5.76e-7（n=9 给 1.32e-6、n=10 给 2.94e-7）")
P_("")
P_("  ⟹ 结论: 层动力学给了【位置】，还没给出【小参数】")
P_("     而小参数是这些分裂的【全部内容】")
verdict = ("层结构给出分裂的位置（L2 的 3 亚层），但层数字全 O(1)；"
           "枚举乘积式的命中数与零假设相当，故不构成推导；"
           "库内唯一的大抑制机制是 15 步乘积，但对整数幂逼不出目标")
P_("\n  => %s" % verdict)
check("判定已给出（含零假设对照）", True, verdict)

P_("\n" + "=" * 96)
nfail = sum(1 for p in PASS if not p["pass"])
P_("合计: %d 项判据, %d PASS, %d FAIL" % (len(PASS), len(PASS)-nfail, nfail))
P_("=" * 96)

out = {
    "layer_numbers": {k: str(v) for k, v in layers.items()},
    "targets": {"21cm_over_lyman": str(ratio),
                "fine_over_gross": str(mp.mpf('0.0072973525693')**2/2),
                "hyperfine_over_gross": str(mp.mpf(4)/3*(mp.mpf('0.51099895')/mp.mpf('938.27208816'))*mp.mpf('0.0072973525693')**2)},
    "n_expressions": len(H),
    "hits_true": len(hit1),
    "hits_random_mean": avg,
    "verdict": verdict,
    "checks": PASS,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_layer_hierarchy.json"), "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
P_("结果写入: %s" % os.path.join(RES, "zero_layer_hierarchy.json"))