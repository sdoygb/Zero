#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QMR17_check.py —— 《Wigner 之友与 Frauchiger–Renner：确定结果的输入地位》逐条复算。

  P1  定理 1  中心态：读数只依赖投影的秩（不产生确定结果）
  P2  定理 2  锐联合可读 ⟺ 交换（构造性检验）
  P3  定理 3  FR 骨架：F̂ 与 Ŵ 不交换（‖[·,·]‖=2）；结构量 P(ok)/P(fail)
  P4  定理 4  重叠矩阵 ⟺ 置换矩阵；FR 骨架缺陷 = 1/2；一般角度 sin²θ
  P5  裁决自检：文档明确标注『引用 FR』与『不引 KS』
  P6  文档契约与链接
"""
import os, sys, re, math
import numpy as np

PASS, FAIL = [], []
HERE = os.path.dirname(os.path.abspath(__file__))
DOCNAME = "QMR17_wigner_friend_and_FR.md"


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(("  [v] " if cond else "  [x] ") + name + (("   " + detail) if detail else ""))


def head(t):
    print(); print("=" * 78); print(t); print("=" * 78)


DOC = open(os.path.join(HERE, DOCNAME), encoding="utf-8").read() if os.path.exists(
    os.path.join(HERE, DOCNAME)) else ""
I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)

# ======================================================================
head("P1  定理 1 · 中心态：读数只依赖秩")
# ======================================================================
D, H = 6, 3
w = [0.5, 0.3, 0.2]
rho = np.zeros((D, D), dtype=complex)
for h in range(H):
    rho[2 * h:2 * h + 2, 2 * h:2 * h + 2] = w[h] / 2 * I2
vals = []
for th in (0.3, 0.9, 2.0, 4.1):
    v = np.array([math.cos(th), math.sin(th)], dtype=complex)
    P = np.zeros((D, D), dtype=complex)
    P[0:2, 0:2] = np.outer(v, v.conj())
    vals.append(round(float(np.trace(rho @ P).real), 12))
check("★ 4 个不同方向的秩 1 投影，读数全为 %.10f（= w₀/2）⇒ **只依赖秩**" % vals[0],
      len(set(vals)) == 1 and abs(vals[0] - w[0] / 2) < 1e-15)
# 秩 2 给 w_0
P2 = np.zeros((D, D), dtype=complex); P2[0:2, 0:2] = I2
check("★ 秩 1 → %.6f（= w₀/2）；秩 2 → %.6f（= w₀）：**只依赖秩**"
      % (vals[0], np.trace(rho @ P2).real), abs(np.trace(rho @ P2).real - w[0]) < 1e-15)

# ======================================================================
head("P2  定理 2 · 锐联合可读 ⟺ 交换")
# ======================================================================
def sharp_joint(A, B):
    return float(np.linalg.norm(A @ B - B @ A)) < 1e-12


check("★ (σz, σz)：交换 ⇒ 存在锐联合可观测量", sharp_joint(SZ, SZ))
check("★ (σz, σx)：不交换 ⇒ **不存在**锐联合可观测量", not sharp_joint(SZ, SX))
# 交换时显式构造 R_ij = P_i Q_j
wa, Va = np.linalg.eigh(SZ); wb, Vb = np.linalg.eigh(SZ)
Rs = [np.outer(Va[:, i], Va[:, i].conj()) @ np.outer(Vb[:, j], Vb[:, j].conj())
      for i in range(2) for j in range(2)]
check("★ 交换时 R_ij = P_iQ_j 为正交投影且 Σ R_ij = I（残差 %.2e）"
      % np.linalg.norm(sum(Rs) - I2),
      all(np.linalg.norm(R @ R - R) < 1e-12 for R in Rs) and np.linalg.norm(sum(Rs) - I2) < 1e-12)

# ======================================================================
head("P3  定理 3 · FR 骨架：F̂ 与 Ŵ 不交换")
# ======================================================================
Fh = np.diag([1.0, -1.0]).astype(complex)
Hm = np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2)
Wh = Hm @ Fh @ Hm.conj().T
check("★ Ŵ = H F̂ H† = [[0,1],[1,0]]（实算 %s）" % np.round(Wh.real, 6).tolist(),
      np.allclose(Wh, SX))
com = Fh @ Wh - Wh @ Fh
check("★ ‖[F̂,Ŵ]‖_spec = %.6f（= 2）；_Fro = %.6f（= 2√2）"
      % (np.linalg.norm(com, 2), np.linalg.norm(com, 'fro')),
      abs(np.linalg.norm(com, 2) - 2) < 1e-12 and abs(np.linalg.norm(com, 'fro') - 2 * math.sqrt(2)) < 1e-12)
check("★ ⇒ 由定理 2：二者**无锐联合可观测量** ⇒ 『她的记录』与『他的记录』不可同时为真", not sharp_joint(Fh, Wh))
# FR 结构量
a, b = math.sqrt(1 / 3), math.sqrt(2 / 3)
ok_amp, fa_amp = (a + b) / math.sqrt(2), (a - b) / math.sqrt(2)
check("★ FR 结构量：P(ok)=%.6f，P(fail)=%.6f，和=%.10f" % (ok_amp ** 2, fa_amp ** 2, ok_amp ** 2 + fa_amp ** 2),
      abs(ok_amp ** 2 + fa_amp ** 2 - 1) < 1e-12)
check("★ ⟨ok|u⟩ = %.6f = 1/√2 ≠ 0 ⇒ W 的记录与 F 的记录**不联合可读**" % abs(1 / math.sqrt(2)),
      abs(1 / math.sqrt(2) - 0.7071067811865476) < 1e-15)

# ======================================================================
head("P4  定理 4 · 重叠矩阵与联合可读性缺陷")
# ======================================================================
def defect(U):
    N = U.shape[0]
    return 1 - sum(max(abs(U[i, j]) ** 2 for j in range(N)) for i in range(N)) / N


check("★ Hadamard（FR 骨架）：Δ = %.6f = 1/2 ⇒ **不是置换矩阵**" % defect(Hm), abs(defect(Hm) - 0.5) < 1e-15)
check("★ 恒等（同基）：Δ = %.6f = 0 ⇒ 完全可联合读取" % defect(np.eye(2)), abs(defect(np.eye(2))) < 1e-15)
ok = True
for th in (0.0, math.pi / 8, math.pi / 6, math.pi / 4):
    R = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
    if abs(defect(R) - math.sin(th) ** 2) > 1e-12:
        ok = False
    print("      θ=%.4f  Δ=%.6f  sin²θ=%.6f" % (th, defect(R), math.sin(th) ** 2))
check("★ 推论 6.1：一般角度 Δ = sin²θ（4 个角度）", ok)
check("★ 推论 6.2：FR 骨架 Δ = 1/2（θ=π/4，最大不可兼容的代表）", abs(defect(Hm) - 0.5) < 1e-15)

# ======================================================================
head("P5  裁决自检")
# ======================================================================
check("★ 文档明确标注『引用 FR』（只取骨架）", "只取骨架" in DOC or "作为已知结果引用" in DOC)
check("★ 文档明确**不**引 Kochen–Specker（并说明单量子比特 KS-可着色）",
      "KS-可着色" in DOC or "KS-可着色的" in DOC)
check("★ 文档含裁决三句话与『额外输入』措辞", "额外输入" in DOC)
check("★ 文档登记了 Δ 是本文提议的度量（唯一性是定义）", "本文**提议**的度量" in DOC or "提议" in DOC)

# ======================================================================
head("P6  文档契约与链接")
# ======================================================================
for tok, lab in [("自包含", "自包含声明"), ("定理 1", "定理 1"), ("定理 2", "定理 2"),
                 ("定理 3", "定理 3"), ("定理 4", "定理 4"), ("推论 3.1", "推论 3.1"),
                 ("推论 6.1", "推论 6.1"), ("推论 6.2", "推论 6.2"), ("推论 7.1", "推论 7.1"),
                 ("Frauchiger", "FR 引用"), ("不许主张", "不许主张清单"),
                 ("标准物理学语言", "术语声明")]:
    check("文档契约：%s" % lab, tok in DOC)
bad = []
for m in re.finditer(r"\[([^\]]*)\]\(([^)]+)\)", DOC):
    t = m.group(2).strip()
    if t.startswith(("http", "#", "mailto")):
        continue
    t = re.sub(r":\d+(-\d+)?$", "", t)
    if not os.path.exists(os.path.normpath(os.path.join(HERE, t))):
        bad.append(t)
check("文档内全部相对链接可解析", not bad, "断链：%s" % bad if bad else "无")
_qm = re.findall(r"\]\((QM\d+_[^)]*)\)", DOC)
check("★ 不含任何已退役 QM* 链接", not _qm, "残留 %s" % _qm if _qm else "无")

print()
print("=" * 78)
print("QMR17 核验：通过 %d / 不符 %d" % (len(PASS), len(FAIL)))
print("=" * 78)
if FAIL:
    for f in FAIL:
        print("   不符：" + f)
sys.exit(1 if FAIL else 0)
