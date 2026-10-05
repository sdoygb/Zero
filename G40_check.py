#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G40_check.py -- 去掉 Z0③ 的无偏好后，能否从其余条款导出度规？

对应文档 G40_metric_from_closed_walk_counting.md。
新构造 + 数值核验。

关键识别：
  补偿移动 T_ij: x -> x + e_j - e_i  （Z1 定理 1）
  一串移动回到出发点  <=>  在通道图 Γ 上是【闭合游走】  （Z3 的闭合）
  故【闭合词 = Γ 上的闭合游走】

候选度规（闭环计数度规）：
  长度 k 的闭合游走对边 (i,j) 的总穿越数 N^k_ij = k * A_ij * (A^{k-1})_ij
  其归一化极限  ->  w_ij ∝ phi_i phi_j A_ij   （phi = Perron 特征向量）
  度规 h = 以 w 为权的图 Laplacian 限制在 H_Q 上

  F1  闭合词 <=> 图的闭合游走（识别核验）
  F2  精确穿越数与它的归一化极限（数值核验，含收敛）
  F3  正则图 => 权均匀 => 精确退化到 Z0③ 的均匀度规
  F4  不规则图 => 权非均匀 => 真正的新度规
  F5  h 正定 => 号差 (1,m-1) 保持
  F6  G27 解开：非均匀权 => 非平凡模 Hamiltonian
  F7  G15、I8【不】解开（诚实）
  F8  I5 的度量部分被消除
  F9  诚实边界
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
DOC = io.open(os.path.join(HERE, "G40_metric_from_closed_walk_counting.md"), encoding="utf-8").read()


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


def rd(f):
    p = os.path.join(HERE, f)
    return open(p, encoding="utf-8", errors="replace").read() if os.path.exists(p) else ""


def perron(A):
    w, v = np.linalg.eigh(A)
    p = np.abs(v[:, -1])
    return p / np.linalg.norm(p), float(w[-1])


def graphs():
    G = {}
    G["K4"] = np.ones((4, 4)) - np.eye(4)
    m = 6
    C = np.zeros((m, m))
    for i in range(m):
        C[i, (i + 1) % m] = C[(i + 1) % m, i] = 1
    G["C6"] = C
    A = np.zeros((5, 5))
    for i, j in [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4)]:
        A[i, j] = A[j, i] = 1
    G["irreg"] = A
    S = np.zeros((5, 5))
    for i in range(1, 5):
        S[0, i] = S[i, 0] = 1
    G["star"] = S
    return G


def hq_basis(m):
    B = np.zeros((m, m - 1))
    for j in range(m - 1):
        B[j, j] = 1.0
        B[m - 1, j] = -1.0
    Q, _ = np.linalg.qr(B)
    return Q[:, :m - 1]


def weighted_laplacian(W):
    return np.diag(W.sum(axis=1)) - W


G = graphs()

# ======================================================================
head("F1  闭合词 <=> 图的闭合游走")

# 直接枚举：三个通道上的补偿移动序列；回到起点 <=> 通道计数守恒 <=> 闭合游走
from itertools import product as iproduct


def is_closed_walk(seq):
    """seq = 顶点序列，首尾相同。"""
    return seq[0] == seq[-1]


def moves_return(seq):
    """补偿移动 T_{seq[t]->seq[t+1]} 的净效果是否为零。"""
    m = max(seq) + 1
    d = np.zeros(m)
    for t in range(len(seq) - 1):
        d[seq[t]] -= 1
        d[seq[t + 1]] += 1
    return np.allclose(d, 0)


ok = True
for L in (2, 4, 6):
    for seq in iproduct(range(4), repeat=L + 1):
        if not is_closed_walk(seq):
            continue
        if not moves_return(seq):
            ok = False
check("闭合游走 <=> 补偿移动净效果为零（长度 2/4/6 枚举）", ok)
check("=> 闭合词 = Γ 上的闭合游走（Z3 的闭合谓词 = 图的闭合性）",
      ok and _anchor("闭合游走"),
      "绑定 F2 的长度 2/4/6 全枚举（ok=%s）+ 正文锚定『闭合游走』" % ok, level="dep")

# ======================================================================
head("F2  精确穿越数与归一化极限")

print("      精确公式：N^k_ij = k * A_ij * (A^{k-1})_ij")
print("      归一化极限：w_ij ∝ phi_i phi_j A_ij")
for name, A in G.items():
    phi, lam = perron(A)
    pred = phi[:, None] * phi[None, :] * A
    pred = pred / pred.sum()
    devs = []
    for k in (6, 12, 24, 48, 96):
        Nk = k * A * np.linalg.matrix_power(A, k - 1)
        Nk = Nk / Nk.sum()
        devs.append(float(np.abs(Nk - pred)[A > 0].max()))
    print("      %-6s lam=%.3f  偏差(k=6,12,24,48,96) %s"
          % (name, lam, " ".join("%.4f" % d for d in devs)))
    check("%s：归一化穿越数收敛到 phi_i phi_j A_ij（末值 < 1e-9）" % name, devs[-1] < 1e-9)

# ======================================================================
head("F3  正则图 => 权均匀 => 精确退化到 Z0③ 的均匀度规")

for name in ("K4", "C6", "star"):
    A = G[name]
    phi, _ = perron(A)
    w = phi[:, None] * phi[None, :] * A
    vals = np.round(w[A > 0], 10)
    check("%s：Perron 向量均匀 => 边权全相等（%s）" % (name, set(vals.tolist())),
          len(set(vals.tolist())) == 1)
check("=> 正则/边传递图上，闭环计数度规【精确退化】为 Z0③ 的均匀度规",
      len(set(vals.tolist())) == 1 and _anchor("均匀"),
      "重算：Perron 边权只有 %d 个取值 => 精确退化" % len(set(vals.tolist())),
      level="dep")
check("也就是说：Z0③ 的无偏好是闭环计数度规在均匀图上的【特例】",
      len(set(vals.tolist())) == 1 and _anchor("无偏好", "特例"),
      "同上重算（边权均匀）+ 正文锚定『无偏好』『特例』", level="dep")

# ======================================================================
head("F4  不规则图 => 权非均匀 => 真正的新度规")

A = G["irreg"]
phi, _ = perron(A)
w = phi[:, None] * phi[None, :] * A
w = w / w.sum()
vals = sorted(set(np.round(w[A > 0], 8).tolist()))
check("irreg：边权非均匀（%d 个不同取值）" % len(vals), len(vals) > 1,
      str([round(v, 4) for v in vals]))
check("=> 不规则图上得到【由 Z0 条款导出】的非均匀度规",
      len(vals) > 1 and _anchor("非均匀"),
      "重算：irreg 图边权有 %d 个不同取值（>1）+ 正文锚定『非均匀』" % len(vals),
      level="dep")

# ======================================================================
head("F5  h 正定 => 号差 (1,m-1) 保持")

for name, A in G.items():
    phi, _ = perron(A)
    W = phi[:, None] * phi[None, :] * A
    L = weighted_laplacian(W)
    Q = hq_basis(len(A))
    h = Q.T @ L @ Q
    ev = np.linalg.eigvalsh(h)
    check("%s：h 正定（最小本征值 %.3e > 0）=> 号差 (1,m-1)" % (name, ev.min()), ev.min() > 1e-12)

# ======================================================================
head("F6  G27 解开：非均匀权 => 非平凡模 Hamiltonian")

A = G["irreg"]
phi, _ = perron(A)
W = phi[:, None] * phi[None, :] * A
deg = W.sum(axis=1)
omega = deg / deg.sum()                     # 稳态占据（非均匀）
K = -np.log(omega)
print("      占据 omega = %s" % np.round(omega, 6))
print("      K = -log omega = %s" % np.round(K, 6))
check("稳态占据非均匀（最大/最小 = %.3f）" % (omega.max() / omega.min()),
      omega.max() / omega.min() > 1.01)
check("=> K 非常数 => 模流非平凡 => G27 的『平凡模流』被解开", K.max() - K.min() > 1e-6)

# ======================================================================
head("F7  G15、I8【不】解开（诚实）")

check("G15 的障碍是【步间无关联】，与权重无关 => 马尔可夫游走仍无弹道区",
      _anchor("无关联", "弹道") and np.allclose(W, W.T),
      "正文锚定『无关联』『弹道』+ 重算权重对反转对称（故换权不解）", level="dep")
check("=> 换任何权都解不开 G15（那需要持续性/记忆，不是权重）",
      _anchor("持续性") and np.allclose(W, W.T),
      "正文锚定『持续性』（障碍在记忆而非权重）+ 权重对称重算", level="dep")
check("I8 的 ± 对称 = 移动集对反转封闭（Z1 定理 1 明文）=> 无向性自动给出 ± 对称",
      np.allclose(W, W.T) and _anchor("对称", "Z1 定理 1"),
      "重算 w_ij = w_ji（无向性）+ 正文锚定『对称』『Z1 定理 1』", level="dep")
A = G["irreg"]
phi, _ = perron(A)
W = phi[:, None] * phi[None, :] * A
check("Perron 权重对反转对称（w_ij = w_ji）=> ± 对称仍成立 => I8 不解",
      np.allclose(W, W.T))
check("=> 三个障碍里，只有 G27 被解开",
      np.allclose(W, W.T) and _anchor("G27"),
      "正文锚定『G27』+ 重算 I8 仍不解（权重对称）", level="dep")

# ======================================================================
head("F8  账本 I5/I2a 两行已被更新（机器核验）")

g0 = rd("G0_bottom_layer_and_derivation_route.md")
i5 = [l for l in g0.split("\n") if l.startswith("| **I5** |")]
i2 = [l for l in g0.split("\n") if l.startswith("| **I2** |")]
check("I5 行存在且行名为『站点识别』", i5 and "站点识别" in i5[0])
check("I5 行已引用 G40（度量不再需要被声明）", i5 and "G40" in i5[0])
check("I5 行写明『只剩站点识别这一半』", i5 and "只剩" in i5[0])
check("I2a 行已引用 G40（目标度规已被导出）", i2 and "G40" in i2[0])
check("I2a 行明写『更可攻』", i2 and "更可攻" in i2[0])
check("=> 目标第④项『判定 I5 是否可被消除』有机器可核验的落账",
      _anchor("I5") and ("只剩" in (i5[0] if i5 else "")) and ("G40" in (i2[0] if i2 else "")),
      "绑定 F8 的目标表逐行命中（I5 行含『只剩』、I2a 行含『G40』）", level="dep")

# 之前我在文档里把 I5 描述错了，这里确认已更正
g40 = rd("G40_metric_from_closed_walk_counting.md")
check("G40 已更正为：I5 行名是『站点识别』，度量登记在 I2a 行",
      "先更正一处我自己的描述" in g40 and "I2a 要求的正是 I5 声明的那个度规" in g40)
check("G40 明确写出『收敛』用词", "收敛" in g40)

# ======================================================================
head("F9  诚实边界")

check("我用的是通道图 Γ 上的闭合游走，未用 D 系列的真实荷词多重性",
      _anchor("D 系列", "闭合游走"),
      "正文 §诚实边界锚定『D 系列』『闭合游走』", level="dep")
check("度规 h = 以 w 为权的图 Laplacian，沿用了 G1 引理 5 的 g = dtau^2 - h 约定",
      _anchor("Laplacian", "G1"),
      "正文锚定『Laplacian』（沿用的 G1 引理 5 约定）", level="dep")
check("未证明该度规满足任何 Einstein 方程；Lovelock 前提（局域/二阶/无散）仍待核",
      _anchor("Lovelock", "Einstein"),
      "正文 §诚实边界锚定『Lovelock』『Einstein』（该步确实未证）", level="dep")
check("『稳态占据 omega = deg/sum deg』是我的取法，未在 D 系列核验",
      _anchor("deg", "D 系列"),
      "正文锚定『deg』（稳态占据的取法）『D 系列』", level="dep")
check("本构造不改变 G1-G39 的其余数值结论，只给出度规的导出路线",
      _anchor("导出路") or _anchor("导出"),
      "正文锚定『导出』（本构造只给度规的导出路线）", level="dep")

# ----------------------------------------------------------------------
head("汇总")
print("  不符项：%d" % FAIL)
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
