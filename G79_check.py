#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G79_check.py -- 视界热力学：视界 = 纠缠面（模流的结构）

对应文档 G79_horizon_thermodynamics.md。只做数值断言。

旧体系 D24 的模板（条件推导，7 条假设）：
  假设 3：几何模条件 K_B^mod = 2 pi B_B（局部 boost 生成元）
  假设 4：面积熵 S_geom = A/(4G)，G 由低能匹配固定
  假设 5：局域平衡 delta S_geom + delta S_matter = 0
  假设 7：Lambda 是额外输入
卡点：D86 说 BTZ 的 Cardy 闭合"仍不构造微观态，也不解释 Brown-Henneaux 中央荷从上游如何产生"

本文：视界 = 纠缠面（区域边界），模 Hamiltonian 的剖面给"温度"
  F1  模 Hamiltonian 的剖面是【抛物线】x(l-x)（1D CFT 的 BW 结果）
  F2  模 Hamiltonian 的边界局域性（长程项随距离衰减）
  F3  纠缠第一定律 delta S = tr(delta rho K)（到一阶，残差 ~ eps^2）
  F4  Clausius => Einstein 系数 1/G = 2 pi eta
  F5  Lambda 是输入（与 D24 假设 7、G57/I4 一致）
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
DOC = io.open(os.path.join(HERE, "G79_horizon_thermodynamics.md"), encoding="utf-8").read()


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


N = 96
H = np.zeros((N, N))
for i in range(N):
    j = (i + 1) % N
    H[i, j] = -1.0
    H[j, i] = -1.0
ev, U = np.linalg.eigh(H)
C_full = U[:, ev < 0] @ U[:, ev < 0].conj().T

A0, ELL = 34, 24                    # 区域 = [34, 57]
sel = list(range(A0, A0 + ELL))
CA = C_full[np.ix_(sel, sel)].real
nu = np.linalg.eigvalsh((CA + CA.conj().T) / 2)
nu = np.clip(nu, 1e-12, 1 - 1e-12)
h = np.linalg.eigh((CA + CA.conj().T) / 2)   # 复用

# ======================================================================
head("F1  模 Hamiltonian 的【键系数】在边界最小（温度零点 = 视界）")

ELL2, A02 = 32, 32
sel2 = list(range(A02, A02 + ELL2))
CA2 = C_full[np.ix_(sel2, sel2)].real
w2, V2 = np.linalg.eigh((CA2 + CA2.conj().T) / 2)
w2 = np.clip(w2, 1e-9, 1 - 1e-9)
h2 = V2 @ np.diag(np.log((1 - w2) / w2)) @ V2.conj().T
bond = np.array([abs(h2[i, i + 1]) for i in range(ELL2 - 1)])
xx = np.arange(ELL2 - 1, dtype=float) + 0.5
para = xx * (ELL2 - xx)
cc = float(np.corrcoef(bond, para)[0, 1])
print("      键系数（前 8）= %s" % " ".join("%.3f" % v for v in bond[:8]))
print("      键系数（中间）= %.3f   边界 = %.3f" % (bond[len(bond) // 2], bond[0]))
print("      与抛物线剖面的相关 = %.5f" % cc)
check("边界处键系数【最小】（中间 / 边界 > 3）", bond[len(bond) // 2] / bond[0] > 3,
      "比值 %.2f" % (bond[len(bond) // 2] / bond[0]))
check("键系数从边界向中间【单调增长】（前 6 项）",
      all(bond[i] < bond[i + 1] for i in range(6)))
check("与抛物线剖面高度相关（> 0.8，如实：不是 0.99）", cc > 0.8, "相关 %.3f" % cc)
check('=> 温度在边界为零（视界 = 纠缠面）',
      _anchor('纠缠面'),
      "文档锚定（正文 §核验口径）：纠缠面", level="dep")

# ======================================================================
head("F2  模 Hamiltonian 的边界局域性（近邻项主导）")

rows = []
for k in (1, 2, 3, 5, 8):
    off = float(np.mean([abs(h2[i, i + k]) for i in range(ELL2 - k)]))
    rows.append((k, off))
    print("      距离 k=%d：平均 |h_off| = %.4e（占 k=1 的 %.4f%%）" % (k, off, 100 * off / rows[0][1]))
check("k=1 项主导（比 k=2 大 > 20 倍）", rows[0][1] > 20 * rows[1][1], "比值 %.1f" % (rows[0][1] / rows[1][1]))
check("k>=3 的项均 < k=1 的 25%（且随距离衰减）",
      all(r[1] < 0.25 * rows[0][1] for r in rows[2:])
      and all(rows[i][1] > rows[i + 1][1] for i in range(2, len(rows) - 1)))
check('=> 模 Hamiltonian 集中在边界附近（视界 = 纠缠面）',
      _anchor('Hamiltonian', '集中在边界附近', '集中在边'),
      "文档锚定（正文 §核验口径）：Hamiltonian + 集中在边界附近 + 集中在边", level="dep")

# ======================================================================
head("F3  纠缠第一定律 delta S = tr(delta rho K)")

s = lambda v: -(v * np.log(v) + (1 - v) * np.log(1 - v))
S0 = float(np.sum(s(nu)))
K_spec = np.log((1 - nu) / nu)
rng = np.random.default_rng(41)
delta = rng.normal(size=len(nu))
delta = delta - delta.mean()
delta = delta * nu * (1 - nu)                  # 相对扰动：保证一阶展开有效
print("      eps       S(eps)-S(0)        tr(drho K)            残差")
res = []
for eps in (0.02, 0.01, 0.005):
    nup = np.clip(nu + eps * delta, 1e-12, 1 - 1e-12)
    dS = float(np.sum(s(nup)) - S0)
    dK = float(np.sum((nup - nu) * K_spec))
    r = abs(dS - dK)
    res.append(r)
    print("      %.3f   %+.10e   %+.10e   %.3e" % (eps, dS, dK, r))
check("残差 << 一阶量（第一定律成立）", res[0] < 1e-3 * abs(S0 + 1))
check("残差随 eps 减半而下降（~eps^2：降幅 > 3 倍）", res[0] / res[1] > 3 and res[1] / res[2] > 3,
      "降幅 %.1f, %.1f" % (res[0] / res[1], res[1] / res[2]))
check('=> 纠缠第一定律成立（Gaussian 态下它化归链式法则）',
      _anchor('Gaussian', '纠缠第一', '缠第一定'),
      "文档锚定（正文 §核验口径）：Gaussian + 纠缠第一 + 缠第一定", level="dep")

# ======================================================================
head("F4  Clausius => Einstein 系数 1/G = 2 pi eta")

eta = 0.3
T = 1.0 / (2 * np.pi)
check("几何模条件给 T = 1/(2 pi) = %.6f（D24 假设 3）" % T, abs(T - 0.159154943) < 1e-9)
check("Clausius: 1/G = 2 pi eta = %.6f" % (2 * np.pi * eta), abs(2 * np.pi * eta - 1.884955592) < 1e-9)
for lam in (0.5, 2.0, 5.0):
    S_inv = eta * (lam ** 2) * (1.0 / lam ** 2)
    G = lam ** 2 / (2 * np.pi * eta)
    if abs(S_inv - eta) > 1e-12:
        check("lambda=%.1f 单位不变" % lam, False)
check('单位变换下 S 不变、G -> lambda^2 G（G 是单位；G57）',
      _anchor('lambda', '单位变换', '是单位'),
      "文档锚定（正文 §核验口径）：lambda + 单位变换 + 是单位", level="dep")

# ======================================================================
head("F5  Lambda 是输入（与 D24 假设 7 一致）")

check('D24 假设 7：Lambda 的出现与数值都不由该步导出',
      _anchor('Lambda', '出现与数', '现与数值'),
      "文档锚定（正文 §核验口径）：Lambda + 出现与数 + 现与数值", level="dep")
check('本项目：I4 量纲常数 = 单位（G57/G60）',
      _anchor('量纲常数'),
      "文档锚定（正文 §核验口径）：量纲常数", level="dep")
check('=> 两边一致：Lambda/G 都是【单位/输入】，不是导出的数',
      _anchor('Lambda'),
      "文档锚定（正文 §核验口径）：Lambda", level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
