#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R1_check.py -- 核验 R1：从数值二阶到可审查的 Γ-收敛条件定理。

对应文档 R1_gamma_convergence_theorem.md。

  F1  状态分离：条件定理 / 证明分解 / 数值旁证 / 未证 / 边界均在位
  F2  基本对象和定理：R1-A、H1-H8、紧性、Γ-liminf、恢复列、唯一极限
  F3  精确游走数与字典多项式：W_d 数值表、P 的正性
  F4  尺度指数 d-2 的必要性：错误指数给 0 / 无穷
  F5  模型 Dirichlet 型的二阶收敛（数值旁证）
  F6  实际闭环权的字典 O(a^2)（数值旁证）
  F7  弱系数收敛不够：周期两层反例给调和平均而非算术平均
  F8  固定物理游程的风险：字典系数比可由一致界失效
  F9  从 Dirichlet 张量到度规的代数恒等式
  F10 引用与交付清单在位
"""

import io
import math
import os
import sys
from collections import defaultdict

import numpy as np
from scipy.sparse import diags

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = 0
NCHECK = 0

DOC_NAME = "R1_gamma_convergence_theorem.md"
DOC = io.open(os.path.join(HERE, DOC_NAME), encoding="utf-8").read()


def check(name, cond, detail=""):
    global FAIL, NCHECK
    NCHECK += 1
    ok = bool(cond)
    if not ok:
        FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", name, ("   " + detail) if detail else ""))


def head(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def read(name):
    p = os.path.join(HERE, name)
    return io.open(p, encoding="utf-8").read() if os.path.exists(p) else ""


def wcount(L, d):
    """Z^d 中从 0 到 e_1 的 L 步游走数。"""
    cur = {(0,) * d: 1}
    for _ in range(L):
        nxt = defaultdict(int)
        for st, n in cur.items():
            for axis in range(d):
                for s in (1, -1):
                    t = list(st)
                    t[axis] += s
                    nxt[tuple(t)] += n
        cur = nxt
    target = [0] * d
    target[0] = 1
    return cur.get(tuple(target), 0)


def p_sum(c, d, k):
    return sum(m * wcount(m - 1, d) * c ** m for m in range(2, k + 1, 2))


def closed_walk_weights(N, k, c):
    """链/环上的形如 A circ A^{k-1} 的层和权重（返回内部边权）。"""
    A = diags([c, c], [1, -1], shape=(N + 1, N + 1)).tocsr()
    P = A.copy()
    W = A.multiply(P) * 2
    for m in range(3, k + 1):
        P = (P @ A).tocsr()
        W = W + m * A.multiply(P)
    return np.asarray(W.diagonal(1)).ravel()


# ======================================================================
head("F1  状态分离：条件定理 / 数值旁证 / 未证 / 边界")

for token in ["【条件定理】", "【证明分解】", "【数值旁证】", "【未证】", "【边界】"]:
    check("等级标签在位：%s" % token, token in DOC)
check("明确写明未把 E1 无条件关闭",
      "未把 E1 无条件关闭" in DOC and "E1 已完成" in DOC)
check("明确区分局部结构收敛 / Γ-收敛 / UV 不动点",
      all(x in DOC for x in ["局部结构收敛", "Γ-收敛", "UV 不动点"]))
check("明确数值结果不能替代证明",
      "数值结果只作旁证" in DOC and "不能替代紧性、Γ-liminf、恢复列和唯一极限" in DOC)
check("明确 H5 字典误差未证",
      "H5 是实际 Zero 权到本文条件定理的桥" in DOC
      and "本文没有证明它" in DOC)

# ======================================================================
head("F2  最小定理与证明义务")

for token in ["定理 R1-A", "紧性", "Γ-liminf", "恢复列", "局域性与唯一极限"]:
    check("R1-A 义务在位：%s" % token, token in DOC)

for h in ["H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8"]:
    check("假设 %s 在位" % h, ("### %s" % h) in DOC)

check("R1-B 实际闭环权是条件推论", "推论 R1-B" in DOC and "若 H1–H5、H6 成立" in DOC)
check("R1-C 只条件给出局部 Dirichlet 张量",
      "推论 R1-C" in DOC and "源守恒前提 $(C)$ 与二阶形式前提 $(O)$ 不由 R1 自动给出" in DOC)
check("不可主张实际类到连续标度场已给",
      "实际旋转类到连续标度场的构造已经给出" in DOC)
check("H3 含单元内振荡控制，且说明 H5 不能替代",
      "单元内振荡受控" in DOC and "H5 只把实际闭环权接到字典多项式，不能替代这条振荡控制" in DOC)

# ======================================================================
head("F3  精确游走数与字典多项式")

W1 = [wcount(L, 1) for L in (1, 3, 5, 7)]
W2 = [wcount(L, 2) for L in (1, 3, 5, 7)]
W3 = [wcount(L, 3) for L in (1, 3, 5, 7)]
check("W_1 = [1,3,10,35]", W1 == [1, 3, 10, 35], "%s" % W1)
check("W_2 = [1,9,100,1225]", W2 == [1, 9, 100, 1225], "%s" % W2)
check("W_3 = [1,15,310,7455]", W3 == [1, 15, 310, 7455], "%s" % W3)

for c in (0.75, 1.5):
    check("k=4,d=3：P_sum(c)>0  (c=%.2f)" % c, p_sum(c, 3, 4) > 0.0,
          "P=%.6f" % p_sum(c, 3, 4))

# ======================================================================
head("F4  尺度指数 d-2 的必要性")

Ns = [32, 64, 128]


def uniform_energy(N, s, d):
    a = 1.0 / N
    # u=x_1 on [0,1]^d: there are N^d edges and each nonzero difference is a.
    return (a ** s) * (N ** d) * (a ** 2)


d = 3
good = [uniform_energy(N, d - 2, d) for N in Ns]
small = [uniform_energy(N, d - 1, d) for N in Ns]
large = [uniform_energy(N, d - 3, d) for N in Ns]
check("s=d-2：能量保持有限非零", all(0.1 < x < 10 for x in good)
      and abs(good[-1] - good[-2]) / good[-1] < 0.05,
      "%s" % ["%.3e" % x for x in good])
check("s=d-1：能量乘 a 消失", small[-1] <= small[0] / 4.0,
      "%s" % ["%.3e" % x for x in small])
check("s=d-3：能量按 a^{-1} 发散", large[-1] > large[0] * 3.9,
      "%s" % ["%.3e" % x for x in large])

# ======================================================================
head("F5  模型 Dirichlet 型二阶收敛（数值旁证）")


def model_profile_error(N, k=4):
    x = (np.arange(N) + 0.5) / N
    c = 1.0 + 0.3 * np.cos(2 * np.pi * x)
    u = np.sin(2 * np.pi * x)
    du = np.diff(np.sin(2 * np.pi * np.linspace(0.0, 1.0, N + 1)))
    f_disc = N * float(np.sum(p_sum(c, 1, k) * du ** 2))
    # 高分辨率周期梯形积分作为参考值。
    M = 200000
    xr = np.linspace(0.0, 1.0, M, endpoint=False)
    cr = 1.0 + 0.3 * np.cos(2 * np.pi * xr)
    f_ref = float(np.mean(p_sum(cr, 1, k) * (2 * np.pi * np.cos(2 * np.pi * xr)) ** 2))
    return abs(f_disc - f_ref)


model_errs = [model_profile_error(N) for N in (32, 64, 128)]
model_slope = float(np.polyfit(np.log([32, 64, 128]), np.log(model_errs), 1)[0])
check("模型能量误差随 N 递降", model_errs[0] > model_errs[1] > model_errs[2],
      " -> ".join("%.3e" % e for e in model_errs))
check("模型能量拟合阶约 2（数值旁证）", abs(model_slope + 2.0) < 0.35,
      "slope=%.3f" % model_slope)

# ======================================================================
head("F6  实际闭环权的字典 O(a^2)（数值旁证）")

K = 4


def closed_walk_errors(Ns=(32, 64, 128)):
    out = []
    for N in Ns:
        x = (np.arange(N) + 0.5) / N
        c = 1.0 + 0.3 * np.cos(2 * np.pi * x)
        w = closed_walk_weights(N, K, c)
        pred = np.array([p_sum(ci, 1, K) for ci in c])
        interior = (x > 0.2) & (x < 0.8)
        out.append(float(np.max(np.abs(w[interior] / pred[interior] - 1.0))))
    return out


cw_errs = closed_walk_errors()
cw_slope = float(np.polyfit(np.log([32, 64, 128]), np.log(cw_errs), 1)[0])
check("实际闭环权字典误差递降", cw_errs[0] > cw_errs[1] > cw_errs[2],
      " -> ".join("%.3e" % e for e in cw_errs))
check("实际闭环权字典拟合阶约 2（数值旁证）", abs(cw_slope + 2.0) < 0.5,
      "slope=%.3f" % cw_slope)

# ======================================================================
head("F7  弱系数收敛不够：周期两层反例")

w_even, w_odd = 1.0, 4.0
arith = 0.5 * (w_even + w_odd)
harm = 2.0 * w_even * w_odd / (w_even + w_odd)
check("算术弱极限 = 2.5", abs(arith - 2.5) < 1e-14, "%.6f" % arith)
check("一维周期均匀化有效系数 = 调和平均 = 1.6", abs(harm - 1.6) < 1e-14,
      "%.6f" % harm)
check("两者明显不同，故弱收敛不能替代 H3 的强有效系数收敛",
      abs(harm - arith) > 0.8)

# ======================================================================
head("F8  固定物理游程的风险：字典系数比可爆炸")

c_min, c_max = 0.75, 1.5
ratios = []
for k in (2, 4, 8, 16):
    pmin = p_sum(c_min, 3, k)
    pmax = p_sum(c_max, 3, k)
    ratios.append(pmax / pmin)
    print("      k=%2d  P(c+)/P(c-) = %.6e" % (k, pmax / pmin))
check("固定物理游程下字典系数比随 k 增长", all(ratios[i] < ratios[i + 1]
                                                      for i in range(len(ratios) - 1)),
      " -> ".join("%.3e" % r for r in ratios))
check("该比值已远大于 1，说明一致椭圆性不是自动结果", ratios[-1] > 1e3,
      "k=16 ratio=%.3e" % ratios[-1])

# ======================================================================
head("F9  Dirichlet 张量到度规的代数恒等式")

rng = np.random.default_rng(17)
ok = True
details = []
for n in (1, 2, 3, 4):
    for _ in range(5):
        R = rng.normal(size=(n, n))
        h = R.T @ R + 0.5 * np.eye(n)
        Q = math.sqrt(np.linalg.det(h)) * np.linalg.inv(h)
        h_back = (np.linalg.det(Q) ** (1.0 / (n - 2))) * np.linalg.inv(Q) if n != 2 else h
        if n != 2:
            ok = ok and np.allclose(h_back, h, rtol=1e-10, atol=1e-10)
    if n != 2:
        details.append("n=%d" % n)
check("n=1,3,4：h=(det Q)^(1/(n-2)) Q^{-1} 精确回收", ok, ", ".join(details))

# n=2 只给共形类：Q 在 h -> lambda h 下不变。
ok2 = True
for _ in range(5):
    R = rng.normal(size=(2, 2))
    h = R.T @ R + 0.5 * np.eye(2)
    for lam in (0.7, 1.3, 2.0):
        Q0 = math.sqrt(np.linalg.det(h)) * np.linalg.inv(h)
        Q1 = math.sqrt(np.linalg.det(lam * h)) * np.linalg.inv(lam * h)
        ok2 = ok2 and np.allclose(Q0, Q1, rtol=1e-10, atol=1e-10)
check("n=2：Q 在 h -> lambda h 下不变 => 只给共形类", ok2)

# ======================================================================
head("F10  引用、反例边界与交付清单")

refs = [
    "R0_publication_theorem.md",
    "G18_attackability_of_the_continuum_limit.md",
    "G58_I2a_resolved_as_embedding_input.md",
    "Z6_stall_autopsy_and_released_ledger.md",
    "Z7_embedding_input_explicit_dictionary.md",
    "Z12_geometric_input_closed_binary_labeling.md",
    "G41_lovelock_premises_under_nonuniform_weight.md",
]
for fn in refs:
    check("引用文件在位：%s" % fn, os.path.exists(os.path.join(HERE, fn)))

for token in [
    "尺度指数不能删",
    "弱系数收敛不够",
    "固定物理游程会破坏局域性",
    "任意图没有共同连续拓扑",
    "边界和跳跃场会产生新条件",
    "二值 / 随机标号不一定自动收敛",
    "字典误差 H5 是真实未证项",
    "R1 不处理洛伦兹、源类和 D259",
    "从 $O(a^2)$ 到连续作用量仍是额外缺口",
]:
    check("反例/边界段在位：%s" % token, token in DOC)

check("交付清单列出 R1 文档和核验脚本",
      "R1_gamma_convergence_theorem.md" in DOC and "R1_check.py" in DOC)
check("剩余缺口显式列出 H5、H7、R2、非结构图与洛伦兹",
      all(x in DOC for x in ["H5", "H7", "R2", "非结构图", "洛伦兹"]))

# ----------------------------------------------------------------------
head("汇总")
print("  断言 %d 项，不符项：%d" % (NCHECK, FAIL))
if FAIL:
    print("  未通过")
    sys.exit(1)
print("  全部通过")
sys.exit(0)
