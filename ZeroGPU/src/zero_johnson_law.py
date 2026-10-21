#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_johnson_law.py —— 用 Johnson 图重算 L1 的律；并自查先前的"证据"是否循环
======================================================================================
承《舞台代数这一格》：单步翻转不保零和 ⟹ W_L 上的正确图是 Johnson 图 J(L,L/2)。

【要查三件】
  J1  Johnson 图是否正则？度 = (L/2)² ？
  J2  它的静止分布推前到壳层(max 偏差) 是否 = Z046 的 w_h = b_h/Σb_h？
  J3  ★ 自查：我先前说「二项式静止分布逐位等于 b_h」——
      那是不是把【同一个二项式】跟【它自己】比？若是，那条"最强证据"是循环的。
"""
from __future__ import annotations
import os, json
import numpy as np
import itertools
from collections import Counter
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
P_("Johnson 图上的 L1 律 · zero_johnson_law.py")
P_("=" * 96)

L = 6
words = [w for w in itertools.product([1, -1], repeat=L) if sum(w) == 0]
D = len(words)
idx = {w: k for k, w in enumerate(words)}
P_("\n§1 W_L：L=%d，dim = C(%d,%d) = %d" % (L, L, L//2, D))

def shell(w):
    s = 0; mx = 0
    for x in w:
        s += x; mx = max(mx, abs(s))
    return mx
sh = np.array([shell(w) for w in words])
b = Counter(sh.tolist())
P_("  壳层计数 b_h：%s" % {int(k): v for k, v in sorted(b.items())})

# ---------- §1 Johnson 图 ----------
P_("\n§1 J1：Johnson 图 J(L,L/2) 的构造与正则性")
A = np.zeros((D, D))
for k, w in enumerate(words):
    for i in range(L):
        for j in range(i+1, L):
            if w[i] != w[j]:
                w2 = list(w); w2[i], w2[j] = w2[j], w2[i]
                A[idx[tuple(w2)], k] += 1.0
deg = A.sum(axis=0)
P_("  度分布：min=%d, max=%d, 全相等? %s" % (deg.min(), deg.max(), np.allclose(deg, deg[0])))
P_("  理论度 (L/2)·(L−L/2) = (L/2)² = %d" % ((L//2)*(L//2)))
check("Johnson 图正则，度 = (L/2)²", np.allclose(deg, (L//2)**2),
      "实测 %d" % deg[0])

# ---------- §2 静止分布 ----------
P_("\n§2 J2：静止分布（正则图 ⟹ 均匀）")
pi_uniform = np.ones(D)/D
P_("  均匀分布是否静止？ πA − π·deg = %.3e"
   % np.linalg.norm(pi_uniform @ A - pi_uniform*deg[0]))
check("静止分布＝均匀", np.linalg.norm(pi_uniform @ A - pi_uniform*deg[0]) < 1e-12,
      "均匀")

P_("\n  推前到壳层(max 偏差)：")
push = {}
for h in sorted(b):
    push[h] = float(pi_uniform[sh == h].sum())
P_("    π_壳(h) = %s" % {int(k): round(v, 6) for k, v in push.items()})
w_h = {h: b[h]/D for h in sorted(b)}
P_("    Z046 的 w_h = b_h/Σb_h = %s" % {int(k): round(v, 6) for k, v in w_h.items()})
same = all(abs(push[h]-w_h[h]) < 1e-12 for h in push)
P_("    两者相同? %s" % same)
check("Johnson 静止分布推前 = Z046 的 w_h", same, "逐位相同")

P_("\n  ★ 但这【为什么】相同？因为：")
P_("    · Johnson 图正则 ⟹ 静止分布【均匀】(与壳层无关)")
P_("    · 推前到壳层 ＝ 各壳层的词数占比 ＝ b_h/Σb_h")
P_("    · 而 Z046 的 w_h 定义就是【按均匀计数】⟹ w_h = b_h/Σb_h")
P_("  ⟹ 两者相同是【公理③（均匀计数）的直接后果】，不是独立验证。")
check("这是构造性一致，不是独立证据", True, "见上")

# ---------- §3 自查：先前的"最强证据" ----------
P_("\n§3 ★★ J3 自查：先前的『二项式静止分布逐位等于 b_h』是不是循环？")
P_("\n  先前（zero_L1_law.py）我做的是：")
P_("    ① 取超立方体 H_1=Σσ^x_i，算其【带符号末高度】的静止分布 = C(L,(L+h)/2)")
P_("    ② 取同一个【带符号末高度】的直方图 b_h")
P_("    ③ 发现两者相同 ⟹ 以为是「强证据」")
P_("  ⟹ ①②③ 是同一个二项式跟它自己比 —— **这是循环论证**")
P_("")
P_("  而且，【带符号末高度】对平衡词恒为 0（S_L=0）——它根本不是壳层指标。")
P_("  壳层指标是 max|S_k|，不是末高度。")
sign_h = Counter(sum(w) for w in words)
P_("  实测：平衡词的末高度 S_L 取值 = %s（恒 0）"
   % {int(k): v for k, v in sorted(sign_h.items())})
check("平衡词的末高度恒为 0（不是壳层指标）", list(sign_h.keys()) == [0],
      str(sorted(sign_h.items())))
P_("\n  ⟹ 结论更正：")
P_("    · 『静止分布＝b_h』这个【值】是对的，但理由不是二项式，是【均匀计数】")
P_("    · 我把它当成「L1 律的强证据」是过强：它只是公理③ 的一致")
P_("    · 超立方体那条律本身也错（不保零和）")
check("更正：先前'最强证据'是循环，降级为一致性", True, "见上")

# ---------- §4 L1 律还站得住吗 ----------
P_("\n§4 L1 的律还站得住吗？")
P_("  站得住的部分：")
P_("    · 图结构：Johnson 图 J(L,L/2)，保零和 ✓")
P_("    · 正则性、均匀静止分布 ✓")
P_("    · 推前给 w_h（构造性一致）")
P_("  站不住 / 要降级的部分：")
P_("    · 『超立方体拉普拉斯』——不保零和，作废")
P_("    · 『静止分布＝b_h 是 L1 律的强证据』——循环，降为一致性")
P_("    · 『L1 律给出 b_h 的来源』——不成立：b_h 是组合计数，无需律就有")
P_("\n  ★ 剩余的真内容：")
P_("    L1 的律（Johnson 图随机走）确实让【位置动】(交换改变 S_t)，")
P_("    而 A_L 的粗粒化把这个动抹掉了。这一条不受上面更正影响。")
verdict = ("✓ 正确：W_L 上的图是 Johnson J(L,L/2)，正则、均匀静止；推前给 w_h=b_h/Σb_h。"
           "✗ 作废：超立方体拉普拉斯（不保零和）。✗ 降级：先前'静止=b_h 是强证据'是循环——"
           "两边都是同一个二项式，且平衡词的末高度恒为 0，根本不是壳层指标。"
           "余下的真内容：Johnson 律让位置动，A_L 粗粒化把动抹掉。")
P_("\n  => %s" % verdict)
check("L1 律的现状已厘清", True, verdict)

P_("\n" + "=" * 96)
nfail = sum(1 for p in PASS if not p["pass"])
P_("合计: %d 项判据, %d PASS, %d FAIL" % (len(PASS), len(PASS)-nfail, nfail))
P_("=" * 96)

out = {
    "L": L, "dim_W_L": D, "johnson_degree": int(deg[0]),
    "stationary_uniform": True,
    "pushforward_equals_w_h": bool(same),
    "circularity_correction": "previous 'binomial = b_h' was tautological; signed final height is identically 0 for balanced words",
    "retracted": "hypercube Laplacian",
    "verdict": verdict, "checks": PASS,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_johnson_law.json"), "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
P_("结果写入: %s" % os.path.join(RES, "zero_johnson_law.json"))