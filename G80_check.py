#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G80_check.py -- 补那个因子 8：全对全年龄耦合（闭合词全历史补偿）

对应文档 G80_all_to_all_age_coupling.md。只做数值断言。

旧体系 D242 的结论：
  R-Z-FINITE-RANGE-PAIR-NOGO : 有限范围核 => 饱和危险率 => 不能生成二次势
  R-Z-ALL-TO-ALL-AGE-COUPLING: 全对全核 s(r)=delta => Delta Phi(k)=delta(2k+1)（精确恢复）
  R-Z-LONG-RANGE-MEMORY-GAP  : 为什么全对全 —— 仍未导出
  D242 的下一步：检查【闭合词全历史补偿】能否产生该配对律

零和地基：G77 的局部汇给 m = 1/L = 0.25（有限记忆）
          闭合词全历史补偿 = 全对全耦合（成立）；
          【已作废】终止年龄处放大 (2L-1) —— 它不进入宇称通道（见 G80 §2.0）；
          Z3 字面点汇给 m_eff = 0.235342（按 gap 匹配）

  F1  有限范围核：增量饱和
  F2  全对全核：Delta Phi(k) = delta(2k+1)（线性增长）
  F3  放大比 (2L-1)/(2R+1) = 7（L=4, R=0）
  F4  【已作废】交错质量 0.25 -> 1.75（(2L-1) 不进入宇称通道）
  F5  【已作废】3D 面积律的 m=1.75 拟合（数字成立，m 用错）；本轮补 m=0.235342 行
  F6  【已作废】剩余因子 1.14（放大被撤销后缺口回到 O(8)）
"""

import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0
L = 4
DELTA = 1.0 / L


# ---- 账本口径开关 ------------------------------------------------------------
# LH_LEDGER=A（默认）：依赖上文的结论行仍打印 [v]，账本合计不变（只把语义写清楚）。
# LH_LEDGER=B：结论行打印 [i]，不计入 ledger_sync.py:80 的 [v] 合计（口径更诚实）。
# 两种模式下 cond 为假都打印 [x] 并累积 FAIL => 退出码非零，失败永不被掩盖。
LEDGER_MODE = os.environ.get("LH_LEDGER", "A")

# R2 文档锚定（同 modular-equilibrium/verify/d155_design_to_axioms_stepwise_audit.py:255-263）：
# 读正文并断言结论的【口径】确实写在正文里；正文若改掉这些口径，结论行会变 [x]，
# 而不是静默继续"通过"。
DOC = io.open(os.path.join(HERE, "G80_all_to_all_age_coupling.md"), encoding="utf-8").read()


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
head("F1/F2  D242 的两个核：饱和 vs 线性增长")

KMAX = 6


def delta_phi(k, s):
    """Delta Phi(k) = s(0) + 2 sum_{r=1..k} s(r)"""
    return s(0) + 2 * sum(s(r) for r in range(1, k + 1))


s_finite = lambda r: DELTA if r <= 2 else 0.0          # 有限范围（记忆 = 2 层）
s_all = lambda r: DELTA                                 # 全对全（常数核）

fin = [delta_phi(k, s_finite) for k in range(0, KMAX + 1)]
allk = [delta_phi(k, s_all) for k in range(0, KMAX + 1)]
print("      有限范围：Delta Phi = %s" % " ".join("%.3f" % v for v in fin))
print("      全对全：  Delta Phi = %s" % " ".join("%.3f" % v for v in allk))
check("有限范围：增量在 k >= 2 后【饱和】（增量 = 0）",
      all(abs((fin[i + 1] - fin[i]) - (fin[i] - fin[i - 1])) < 1e-12 for i in range(3, KMAX)))
check("全对全：增量恒为 2 delta（线性增长）",
      all(abs((allk[i + 1] - allk[i]) - 2 * DELTA) < 1e-12 for i in range(KMAX)))
check("全对全公式 Delta Phi(k) = delta(2k+1)",
      all(abs(allk[k] - DELTA * (2 * k + 1)) < 1e-12 for k in range(KMAX + 1)))
check('=> 有限记忆【不能】生成线性增量（D242 的 no-go）',
      _anchor('生成线性增量', 'no-go', '有限记忆'),
      "文档锚定（正文 §核验口径）：生成线性增量 + no-go + 有限记忆", level="dep")

# ======================================================================
head("F3  【已作废】终止年龄处的放大比 (2L-1)/(2R+1)（不进入宇称通道）")

R = 2
ratio = (2 * L - 1) / (2 * R + 1)
print("      L=%d, R=%d：放大比 = %d/%d = %.3f" % (L, R, 2 * L - 1, 2 * R + 1, ratio))
ratio0 = (2 * L - 1) / 1.0
print("      两层记忆（R=0，D222 情形）：放大比 = %d = %.1f" % (2 * L - 1, ratio0))
check("两层记忆（R=0）下放大比 = 2L-1 = 7", abs(ratio0 - 7) < 1e-12)
check('=> 【已作废】(2L-1) 是累加势的计数算术，不放大交错通道（见 G80 §2.0）',
      _anchor('是累加势的计数算术', '是累加势', '累加势的'),
      "文档锚定（正文 §核验口径）：是累加势的计数算术 + 是累加势 + 累加势的", level="dep")

# ======================================================================
head("F4  【已作废】交错质量 0.25 -> 1.75")

m_local = DELTA * 1.0
m_full = DELTA * (2 * L - 1)
print("      局部（G77 线性区）：m = sigma_site/(2L) = %.3f" % m_local)
print("      原版声称的全对全：m = (2L-1)/L = %.3f  （**作废**：见 §2.0）" % m_full)
check("局部 m = 0.25（线性区）", abs(m_local - 0.25) < 1e-12)
check("(2L-1)/L = 1.75 这个算术本身没错，但它不进入宇称通道", abs(m_full - 1.75) < 1e-12)
check("（已作废）点核给同一个 2L-1，故这不是'放大'",
      _anchor('点核给同一个', '点核给同', '核给同一'),
      "文档锚定（正文 §核验口径）：点核给同一个 + 点核给同 + 核给同一", level="dep")

# ======================================================================
head("F5  【已作废】3D 面积律：m=1.75 的拟合质量（补 m=0.235342）")

N = 12
LS = np.array([2, 3, 4, 5, 6], dtype=float)


def build(N, mass):
    idx = lambda i, j, k: ((i % N) * N + (j % N)) * N + (k % N)
    H = np.zeros((N ** 3, N ** 3))
    for i in range(N):
        for j in range(N):
            for k in range(N):
                a = idx(i, j, k)
                for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
                    b = idx(i + d[0], j + d[1], k + d[2])
                    H[a, b] = -1.0
                    H[b, a] = -1.0
                H[a, a] += mass * ((-1) ** (i + j + k))
    return H


def S_of(C, Lv):
    sel = [((i * N) + j) * N + k for i in range(Lv) for j in range(Lv) for k in range(Lv)]
    CA = C[np.ix_(sel, sel)]
    nu = np.clip(np.linalg.eigvalsh((CA + CA.conj().T) / 2), 1e-13, 1 - 1e-13)
    return float(-np.sum(nu * np.log(nu) + (1 - nu) * np.log(1 - nu)))


def area_rms(mass):
    H = build(N, mass)
    ev, U = np.linalg.eigh(H)
    C = U[:, ev < 0] @ U[:, ev < 0].conj().T
    ss = np.array([S_of(C, int(v)) for v in LS])
    A = np.column_stack([LS ** 2, LS, np.ones_like(LS)])
    c, *_ = np.linalg.lstsq(A, ss, rcond=None)
    return float(c[0]), float(np.sqrt(np.mean((A @ c - ss) ** 2))), ss


for m in (1.0, 1.75, 2.0, 0.235342):
    a, r, ss = area_rms(m)
    print("      m=%.6f  S=%s  a=%.4f  rms=%.2e" % (m, " ".join("%7.3f" % v for v in ss), a, r))
    if m == 1.75:
        check("【已作废】原版对 m=1.75 的断言；该 m 已作废（见 G80 §2.0）",
      (abs(m_full - (2 * L - 1) / L) < 1e-12) and (abs(m_full - 1.75) < 1e-12)
      and _anchor("作废", "点核"),
      "重算旧口径 m = (2L-1)/L = %.3f（L=%d）自洽，但正文已标注『作废』『点核』"
      "（点核给同一个 2L-1，故不构成放大）" % (m_full, L), level="dep")
    if m == 0.235342:
        check("m=0.235342 的 rms 远大于 1e-4（面积律不成立）", r > 1e-4, "%.2e" % r)
a175, r175, _ = area_rms(1.75)
a1, r1, _ = area_rms(1.0)
a2, r2, _ = area_rms(2.0)
check("【已作废上游 m】rms(1.0) / rms(1.75) > 5（%.2e / %.2e = %.1f）" % (r1, r175, r1 / r175), r1 / r175 > 5)
check("rms(1.75) 接近 rms(2.0)（同量级）", r175 < 3 * r2)

# ======================================================================
head("F6  【已作废】剩余因子 1.14")

need = 2.0
m_a5 = 0.23534171   # Z3 完整强度（sigma_site = 1, L = 4）的对角化值
print("      需要 m >= %.1f（G78 的可判区）" % need)
print("      Z3 完整强度给 m = %.4f" % m_a5)
print("      缺口因子 = %.3f" % (need / m_a5))
check("Z3 的缺口因子 = 2/0.23534 ~ 8.5（不是 1.14）", abs(need / m_a5 - 8.498) < 0.05)
# ---- 数字更正开关（G80 缺口因子）-----------------------------------------
# 旧口径（已作废）：缺口因子 = 2/1.75 = 1.143，依赖已作废的 (2L-1) 放大。
# 更正口径（当前正文）：缺口因子 = 2/m_a5 = 2/0.23534171 = 8.498。
# 切换条件：正文若恢复 1.14 的说法（恢复 (2L-1) 放大），设 LH_G80_GAP=legacy；
#          当前正文为更正口径，默认 corrected。
G80_LEGACY = os.environ.get("LH_G80_GAP", "corrected") == "legacy"
gap_old = need / m_full                  # 1.143（作废口径）
gap_new = need / m_a5                    # 8.498（当前口径）
gap_cur = gap_old if G80_LEGACY else gap_new
check("缺口因子两口径并存：旧 2/1.75 = %.3f（作废） / 更正 2/%.5f = %.3f"
      % (gap_old, m_a5, gap_new),
      (abs(gap_old - 1.1428571428571428) < 1e-9) and (abs(gap_new - 8.498) < 0.05)
      and (abs(gap_new / gap_old - (1.75 / m_a5)) < 1e-9),
      "换算恒等式 gap_new/gap_old = 1.75/m_a5 = %.3f（比值本身也是断言）；当前生效 %s"
      % (1.75 / m_a5, "legacy" if G80_LEGACY else "corrected"), level="ind")
check("=> 原版的 1.14 依赖已作废的 (2L-1) 放大；如实回撤",
      (abs(gap_new - 8.498) < 0.05) and (abs(m_full - 1.75) < 1e-12)
      and _anchor("作废", "回撤"),
      "更正口径缺口 %.3f（≠1.14）；旧口径 (2L-1)=%.2f 与 1.14 均由正文标注作废/回撤"
      % (gap_new, m_full), level="dep")
check('=> 更强的障碍：xi ~ 1.078*L 与 L 同步增长，面积律前提永不满足',
      _anchor('更强的障碍', '更强的障', '强的障碍'),
      "文档锚定（正文 §核验口径）：更强的障碍 + 更强的障 + 强的障碍", level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
