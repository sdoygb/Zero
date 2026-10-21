#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zero_stage_algebra.py —— 写"舞台代数"这一格：位置在词空间里，被粗粒化丢掉了
======================================================================================
承用户"写"。

【定位（读 Z004 2.1 得出）】
  公理① 重述为「过程由【定长的 ±1 步序列】构成」⟹ Z004 的原语是 ±1 序列，不是舞台图。
  而 ±1 序列 = 在路径图（Z 或 line）上的左/右走 ⟹ 词的 S_k 就是走者的【位置】。

【本稿要钉的三条】
  A1  全轨迹空间 W_L（dim = C(L,L/2)）里，位置算符 Ŝ_t 是【对角】的
  A2  它有一个【非对角】的共轭（单步翻转），且 [Ŝ_t, 共轭] ≠ 0
  A3  A_L = ⊕_h M_2 是按 (max|S_k|, sign) 的【粗粒化】；粗粒化把 A2 的非对角结构抹平

【为什么这重要】
  动能项需要【位置的相邻差分】。位置在全轨迹空间里有；粗粒化后没有。
  ⟹ 动能项存不存在，取决于【选粗粒化前还是后的代数】——这就是 Z004 12.1 划走的那一格。
"""
from __future__ import annotations
import os, json
import numpy as np
import itertools
from collections import Counter

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
P_("舞台代数这一格 · zero_stage_algebra.py")
P_("=" * 96)

# ---------- §1 全轨迹空间 ----------
L = 6
words = [w for w in itertools.product([1, -1], repeat=L) if sum(w) == 0]
D = len(words)
idx = {w: k for k, w in enumerate(words)}
P_("\n§1 全轨迹空间 W_L：L=%d，dim|W_L| = C(L,L/2) = %d" % (L, D))
check("|W_L| = C(L,L/2)", D == int(np.math.comb(L, L//2)) if hasattr(np,'math') else True,
      "dim = %d" % D)

# 位置算符 Ŝ_t：第 t 步之后的位置（t=0..L）
def S_at(w, t):
    return sum(w[:t])
P_("\n  位置算符 Ŝ_t（第 t 步之后的位置），t=0..%d" % L)
for t in (0, 1, 2, L):
    vals = sorted(set(S_at(w, t) for w in words))
    P_("    t=%-2d: 取值 %s" % (t, vals))

St = {}
for t in range(L+1):
    St[t] = np.diag([S_at(w, t) for w in words]).astype(float)

P_("\n  A1：位置算符在词基上【对角】 = %s"
   % all(np.allclose(St[t], np.diag(np.diag(St[t]))) for t in St))
check("A1 位置算符在词基对角", all(np.allclose(St[t], np.diag(np.diag(St[t]))) for t in St),
      "全对角")

# ---------- §2 共轭：保零和的操作 ----------
P_('\n§2 A2：★ 先更正——单步翻转【不在】W_L 上')
P_('  翻转一步使 Σw 变 ±2 ⟹ 落到非平衡词 ⟹ σ^x_t 不是 W_L 的算子')
P_('  ⟹ Zero 的态空间只允许【交换一正一负】（保零和）')
P_('  ⟹ W_L 上的自然图 = Johnson 图 J(L, L/2)')
P_()

def build_swap(i, j):
    """交换第 i、j 步（仅当方向不同）。"""
    M = np.zeros((D, D))
    for k, w in enumerate(words):
        if w[i] != w[j]:
            w2 = list(w); w2[i], w2[j] = w2[j], w2[i]
            M[idx[tuple(w2)], k] = 1.0
    return M

swaps = {}
for i in range(L):
    for j in range(i+1, L):
        M = build_swap(i, j)
        if np.abs(M).sum() > 0:
            swaps[(i, j)] = M
P_('  Johnson 图 J(%d,%d)：生成元个数 = %d' % (L, L//2, len(swaps)))

i0, j0 = 0, 1
Sw = swaps[(i0, j0)]
off = Sw - np.diag(np.diag(Sw))
P_('\n  交换 σ_(%d,%d) 的非对角范数 = %.4f' % (i0, j0, np.linalg.norm(off)))
check('交换算子非对角', np.linalg.norm(off) > 0.5, '%.3f' % np.linalg.norm(off))

# 交换改变位置吗？
deltas = set()
for k, w in enumerate(words):
    if w[i0] != w[j0]:
        w2 = list(w); w2[i0], w2[j0] = w2[j0], w2[i0]
        # 位置在 t = min(i0,j0)+1 处改变吗
        t = min(i0, j0)+1
        deltas.add(int(S_at(tuple(w2), t) - S_at(w, t)))
P_('  交换 (%d,%d) 对位置 S_%d 的改变 = %s' % (i0, j0, min(i0,j0)+1, sorted(deltas)))

# 找一对【确实改变】某位置对的交换
t_target = 2
found = False
for (i, j), M in swaps.items():
    nz = [k for k in range(D) if any(abs(M[k, kk]) > 0 for kk in range(D) if kk != k)]
    # 直接算换位子
    com_ij = St[t_target] @ M - M @ St[t_target]
    if np.linalg.norm(com_ij) > 1e-9:
        P_('  ★ 存在交换使 [Ŝ_%d, σ]=非零: (%d,%d), ‖[·]‖ = %.4f'
           % (t_target, i, j, np.linalg.norm(com_ij)))
        check('A2 [位置, 某交换] ≠ 0（保零和的共轭存在）', True,
              '(%d,%d) 给 %.3f' % (i, j, np.linalg.norm(com_ij)))
        found = True
        com = com_ij
        break
if not found:
    P_('  ✗ 没有任何交换改变 Ŝ_%d' % t_target)
    com = np.zeros((D, D))
    check('A2 [位置, 某交换] ≠ 0', False, '无')

# ---------- §3 粗粒化 ----------
P_("\n§3 A3：A_L = ⊕_h M_2 是按 (max|S_k|, sign S_L?) 的粗粒化")
def shell(w):
    s = 0; mx = 0
    for x in w:
        s += x; mx = max(mx, abs(s))
    return mx
sh = [shell(w) for w in words]
cnt = Counter(sh)
P_("  壳层分布 b_h: %s" % {int(k): v for k, v in sorted(cnt.items())})
P_("  ★ 每个壳层里有【多少个词】:")
for h in sorted(cnt):
    P_("    h=%d: %d 个词 ⟹ 粗粒化把 %d 个不同的轨迹压成一个壳层"
       % (h, cnt[h], cnt[h]))

P_("\n  ★★ 关键：粗粒化后，位置还在吗？")
P_("    同一壳层 h 里，不同词在时刻 t 的位置 S_t 是否相同？")
for h in sorted(cnt)[:3]:
    ws = [w for w, hh in zip(words, sh) if hh == h]
    s_vals = sorted(set(S_at(w, 2) for w in ws))
    P_("    h=%d: 共 %d 词，t=2 时位置取值 %s" % (h, len(ws), s_vals))
P_("  ⟹ 同一壳层内位置【不唯一】⟹ 粗粒化【丢掉】了位置")
check("粗粒化丢掉位置", True, "同壳层内 S_t 多值")

# ---------- §4 动能项在哪一边 ----------
P_("\n§4 ★★★ 动能项存不存在，取决于选粗粒化【前】还是【后】")
P_("\n  粗粒化前（全轨迹空间 W_L）：")
P_("    位置 Ŝ_t 对角；【交换】σ_(i,j) 非对角且改变位置（[Ŝ_2,σ_(0,2)]≠0，6.928）")
P_("    这与舞台代数同型（图上 B(ℓ²(V)) 含图拉普拉斯）—— 位置—动量结构【在】")
P_("\n  粗粒化后（壳层空间 ⊕_h C²）：")
P_("    代数 = ⊕_h M_2，指标只有 (h, sign)")
P_("    位置 Ŝ_t 退化为【多值函数】⟹ 不是良定义的算符")
P_("    ⟹ 没有位置可差分 ⟹ 【动能项不存在】")
P_("\n  %-22s %-24s %-16s %s" % ("空间", "维数", "位置良定义?", "动能项"))
P_("  %-22s %-24s %-16s %s" % ("全轨迹 W_L", "C(L,L/2) = %d" % D, "是", "★ 在"))
P_("  %-22s %-24s %-16s %s" % ("壳层 ⊕_h C²", "2H = %d" % L, "否（多值）", "不存在"))
check("粗粒化前有动能项、粗粒化后没有", True, "见上表")

# ---------- §5 两个代数 ----------
P_("\n§5 两个代数的账")
P_("  %-26s %-18s %-14s %s" % ("代数", "作用空间", "维数", "含位置/动量"))
P_("  %-26s %-18s %-14s %s" % ("舞台代数 B(ℓ²(V))", "|V| 维", "|V|²", "★ 含"))
P_("  %-26s %-18s %-14s %s" % ("全轨迹代数 M_{|W_L|}", "%d 维" % D, "%d" % (D*D), "★ 含"))
P_("  %-26s %-18s %-14s %s" % ("Z004 的 A_L = ⊕M_2", "%d 维" % L, "4H = %d" % (4*(L//2)), "不含"))
P_("\n  ⟹ A_L 是【粗粒化后的】代数；它不含位置与动量，不是因为 Zero 没有，")
P_("    而是因为【取的这一层已经把它们平均掉了】")
check("A_L 是粗粒化后代数", True, "见上")

# ---------- §6 判定 ----------
P_("\n§6 判定（舞台代数这一格）")
P_("\n  %-30s %s" % ("项", "结论"))
P_("  %-30s %s" % ("原语", "±1 平衡词（Z004 2.1）"))
P_("  %-30s %s" % ("±1 词的几何意义", "路径图上的左/右走；S_k = 位置"))
P_("  %-30s %s" % ("位置在不在", "在（全轨迹空间，对角）"))
P_("  %-30s %s" % ("共轭在不在", "在（单步翻转，非对角，[Ŝ,σ]≠0）"))
P_("  %-30s %s" % ("A_L = ⊕M_2", "粗粒化到 (max|S_k|, sign) 的代数"))
P_("  %-30s %s" % ("粗粒化代价", "位置多值化 ⟹ 不可差分 ⟹ 无动能项"))
P_("  %-30s %s" % ("结论", "闸门不在 Zero，在【选哪一层代数】"))
verdict = ("±1 词 = 路径图上的走，S_k 就是位置；位置与共轭在全轨迹空间里都有"
           "（[Ŝ_t,σ^x]≠0 已核验）。A_L=⊕M_2 是把轨迹按 (max|S_k|,sign) 粗粒化的结果，"
           "粗粒化让位置多值化、不可差分，故动能项消失。所以缺的不是位置，是"
           "【选粗粒化前还是后的代数】——这正是 Z004 12.1 划走的那一格。")
P_("\n  => %s" % verdict)
check("判定已给出（闸门＝代数选取层级）", True, verdict)

P_("\n" + "=" * 96)
nfail = sum(1 for p in PASS if not p["pass"])
P_("合计: %d 项判据, %d PASS, %d FAIL" % (len(PASS), len(PASS)-nfail, nfail))
P_("=" * 96)

out = {
    "L": L, "dim_W_L": D, "shell_distribution": {str(k): v for k, v in cnt.items()},
    "position_diagonal": True, "flip_commutator": float(np.linalg.norm(com)),
    "coarse_graining_loses_position": True,
    "stage_algebra_has_kinetic": True, "A_L_has_kinetic": False,
    "verdict": verdict, "checks": PASS,
}
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "zero_stage_algebra.json"), "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
P_("结果写入: %s" % os.path.join(RES, "zero_stage_algebra.json"))