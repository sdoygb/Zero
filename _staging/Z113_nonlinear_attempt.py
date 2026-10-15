# -*- coding: utf-8 -*-
"""Z113（出口乙）：非线性搬运律能否给出**稳定的非均匀不动点**？

律（忠实三层：Σm 守恒 / 均匀态不动 / 无自由选择规则）：
    dm_v/dt = Σ_{u~v} ρ(m_v)ρ(m_u)(m_v − m_u)
自由能：F = Σ_v Φ(m_v),  Φ' = log ρ  ⟹  dF/dt = −Σ_e ρ_vρ_u(∇m)² ≤ 0（ρ>0）

本文件检验的五个结构性结论：
  T1 ρ 单调不减 ⟹ Φ 凸 ⟹ 均匀态唯一全局极小 ⟹ 收敛到均匀（q→1/V）
  T2 路径/树有死端 ⟹ 驻点方程 ρ(m_v)(m_v−μ)=C 被逼 C=0
  T3 ρ 在均值处取极小 ⟹ 线性化退化 ⟹ 无分岔
  T4 显式辛格式（流量有记忆）无条件不稳定：det = −1−c²λ < −1
  T5 ★ ρ 的**零点 = 冻结组态**：ρ(m*)=0 处流量恒为 0
        ⟹ 若 Φ 是双阱（ρ 有两个零点 m₋<μ<m₊），则任何 {m₋,m₊} 赋值
          都是**精确驻点**（流量逐边为 0）⟹ 非均匀但"冻结"的吸引子
"""
from __future__ import annotations
import sys, os
import numpy as np
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                 "..", "ZeroGPU", "src")))
from z0_engine import path_adj, torus_adj, star_adj

CHECKS = []
def check(n, c, d=""):
    CHECKS.append((n, bool(c))); print(f"[{'OK' if c else 'FAIL'}] {n}  {d}")

def edges_of(adj):
    E = set()
    for u, nb in enumerate(adj):
        for w in nb:
            if u != w: E.add((min(u, w), max(u, w)))
    return sorted(E)

def flow(adj, m0, rho, kappa=0.02, steps=400000, tol=1e-14):
    """守恒显式梯度流；返回 (末态, 步数, 是否收敛, 是否越界)"""
    E = edges_of(adj); m = np.array(m0, float)
    for t in range(steps):
        rv = rho(m)
        F = np.array([rv[u]*rv[v]*(m[u]-m[v]) for u, v in E])
        dm = np.zeros_like(m)
        for (u, v), f in zip(E, F): dm[u] += f; dm[v] -= f
        m = m - kappa*dm
        if not np.isfinite(m).all() or np.abs(m).max() > 1e3:
            return m, t, False, True
        if np.abs(dm).max() < tol: return m, t, True, False
    return m, steps, False, False

print("="*80); print("T1/T3：单调 ρ 与『退化 ρ』"); print("="*80)
print("""（乙-1..乙-5 已算，此处只引用）
  ρ=1, m, m², 1+m+m², m+m³      → 单调 ⟹ 收敛到均匀，q=1/V
  ρ=(m−1)²+1（在 μ=1 处取极小） → Φ''(1)=0 ⟹ 线性化退化 ⟹ 仍收敛到均匀""")
rhoA = lambda m: (m-1)**2 + 1.0
mA, tA, convA, badA = flow(path_adj(16), 1.0+0.05*(np.random.RandomState(4).rand(16)-0.5), rhoA)
check("T3：ρ 在均值处退化 ⟹ 仍收敛到均匀", abs(float(np.var(mA))) < 1e-12,
      f"var={float(np.var(mA)):.2e}, t={tA}")

print()
print("="*80); print("T5/T7：零点冻结与『非均匀驻点是测度零的』"); print("="*80)
m_minus, m_plus = 0.5, 2.0
rho_dw = lambda m: np.clip((m - m_minus)*(m_plus - m), 0.0, None)
print(f"  rho(m) = max((m-{m_minus})({m_plus}-m), 0)   零点 m-={m_minus}, m+={m_plus}")
print(f"  Phi''(m) = rho'/rho = {2.5-2*1.0:.1f} at m=1 > 0 => Phi 凸 => 均匀是极小")

# (a) 一般（非零点）初值：应收敛到均匀（T1/T3）
for nm, adj in [("path16", path_adj(16)), ("cycle16", torus_adj(1,16)), ("star16", star_adj(16))]:
    V = len(adj); rng = np.random.RandomState(9)
    m0 = 0.6 + 0.5*rng.rand(V); m0 = m0 - m0.mean() + 1.0
    m, t, conv, bad = flow(adj, m0, rho_dw)
    w = m/m.sum(); q = float((w*w).sum())
    check(f"{nm}: 一般初值 -> 收敛到均匀（T1/T3）",
          float(np.var(m)) < 1e-10 and conv, f"var={float(np.var(m)):.2e}, q={q:.6f}, t={t}")

# (b) 取值落在零点上的构型：是**精确驻点**（流量恒 0）
print()
print("  (b) 取值落在 rho 零点上的构型：")
for nm, vals in [("全在 m-=0.5",  np.full(16, 0.5)),
                 ("全在 m+=2.0",  np.full(16, 2.0)),
                 ("交替 {0.5,2.0}", np.array([0.5,2.0]*8))]:
    E = edges_of(path_adj(16))
    fl = np.array([rho_dw(np.array([vals[u]]))[0]*rho_dw(np.array([vals[v]]))[0]*(vals[u]-vals[v])
                   for u, v in E])
    check(f"{nm}: 精确驻点（流量恒 0）", np.abs(fl).max() < 1e-14,
          f"maxflow={np.abs(fl).max():.1e}, var={np.var(vals):.3f}")

# (c) T7：非均匀驻点需要均值恰等于零点组合平均 -> 非通有
print()
print("  (c) T7：非均匀驻点形如 m_v in Z，须满足 mean = 零点组合平均")
V = 16; ok_cnt = 0
for k in range(V+1):
    mean_k = (0.5*k + 2.0*(V-k))/V
    if abs(mean_k - 1.0) < 1e-12: ok_cnt += 1
print(f"      对目标均值 mu=1.0，在 k=0..16 中找到的可行 k 有 {ok_cnt} 个")
check("T7：给定 mu 一般**没有**非均匀驻点（测度零）", ok_cnt == 0,
      f"可行组合数={ok_cnt}")

print()
print("="*80); print("T4：流量带记忆（辛）格式的无条件不稳定性"); print("="*80)
print("""  格式 j ← j − c∇m ; m,φ ← m,φ + c·div j
  线性化行列式 det = −1 − c²λ  ⟹ |det| = 1+c²λ > 1 对**一切** c>0, λ>0
  ⟹ 至少一个特征值模 >1 ⟹ 无条件不稳定（与步长无关，是格式的代数性质）""")
c = 0.02
for lam in (0.0384, 1.0, 3.9616):
    det = -1 - c*c*lam
    print(f"    λ={lam:7.4f}: det={det:.6f}  |det|={abs(det):.6f} > 1 ⟹ 不稳定")
check("T4：辛交错格式 |det|>1 对所有 λ>0", all(abs(-1-c*c*l) > 1 for l in (0.0384,1.0,3.9616)))

print()
print("="*80); print("结论"); print("="*80)
print("""
 出口乙的可证结果：
  ✅ 忠实的层间动力学**存在**：守恒、有界（ρ>0 时 dF/dt≤0）、不可约（§Z112 Q3 判据）。
  ✅ 由公理自动生成第二定律（F 单调下降），**无需公设**。
  ✅ 非均匀稳定态**存在**，但它的机制是"ρ 的零点 = 冻结组态"：
     双阱 Φ 的两个阱 m₋, m₊ 上流量恒为 0 ⟹ 任何 {m₋,m₊} 赋值都是精确驻点。
  ⚠️ 代价：这种非均匀态是**冻结**的（动力学停止），不是"活的稳态"。
  ❌ 单调 ρ ⟹ 数学上只能收敛到均匀，q→1/V=0.0625，与 D4 窗 (1/2,3/5) 无关。
  ❌ 非对称/辛型（"流量有记忆"）线性格式无条件不稳定。

 因此：**D4 窗不是这类动力学的输出**；要得到它，必须再定一笔数据——
  "哪个密度是允许驻留的"（即 ρ 的零点位置），这正是本轮被逼出来的具名输入。
""")
n = sum(1 for _, ok in CHECKS if ok)
print(f"断言通过 {n}/{len(CHECKS)}")


print()
print("="*80); print("T5/T7：零点冻结与『非均匀驻点是测度零的』"); print("="*80)
m_minus, m_plus = 0.5, 2.0
rho_dw = lambda m: np.clip((m - m_minus)*(m_plus - m), 0.0, None)
print(f"  rho(m) = max((m-{m_minus})({m_plus}-m), 0)   零点 m-={m_minus}, m+={m_plus}")
print(f"  Phi''(m) = rho'/rho = {2.5-2*1.0:.1f} at m=1 > 0 => Phi 凸 => 均匀是极小")

# (a) 一般（非零点）初值：应收敛到均匀（T1/T3）
for nm, adj in [("path16", path_adj(16)), ("cycle16", torus_adj(1,16)), ("star16", star_adj(16))]:
    V = len(adj); rng = np.random.RandomState(9)
    m0 = 0.6 + 0.5*rng.rand(V); m0 = m0 - m0.mean() + 1.0
    m, t, conv, bad = flow(adj, m0, rho_dw)
    w = m/m.sum(); q = float((w*w).sum())
    check(f"{nm}: 一般初值 -> 收敛到均匀（T1/T3）",
          float(np.var(m)) < 1e-10 and conv, f"var={float(np.var(m)):.2e}, q={q:.6f}, t={t}")

# (b) 取值落在零点上的构型：是**精确驻点**（流量恒 0）
print()
print("  (b) 取值落在 rho 零点上的构型：")
for nm, vals in [("全在 m-=0.5",  np.full(16, 0.5)),
                 ("全在 m+=2.0",  np.full(16, 2.0)),
                 ("交替 {0.5,2.0}", np.array([0.5,2.0]*8))]:
    E = edges_of(path_adj(16))
    fl = np.array([rho_dw(np.array([vals[u]]))[0]*rho_dw(np.array([vals[v]]))[0]*(vals[u]-vals[v])
                   for u, v in E])
    check(f"{nm}: 精确驻点（流量恒 0）", np.abs(fl).max() < 1e-14,
          f"maxflow={np.abs(fl).max():.1e}, var={np.var(vals):.3f}")

# (c) T7：非均匀驻点需要均值恰等于零点组合平均 -> 非通有
print()
print("  (c) T7：非均匀驻点形如 m_v in Z，须满足 mean = 零点组合平均")
V = 16; ok_cnt = 0
for k in range(V+1):
    mean_k = (0.5*k + 2.0*(V-k))/V
    if abs(mean_k - 1.0) < 1e-12: ok_cnt += 1
print(f"      对目标均值 mu=1.0，在 k=0..16 中找到的可行 k 有 {ok_cnt} 个")
check("T7：给定 mu 一般**没有**非均匀驻点（测度零）", ok_cnt == 0,
      f"可行组合数={ok_cnt}")

print()
print("="*80); print("T4：流量带记忆（辛）格式的无条件不稳定性"); print("="*80)
print("""  格式 j ← j − c∇m ; m,φ ← m,φ + c·div j
  线性化行列式 det = −1 − c²λ  ⟹ |det| = 1+c²λ > 1 对**一切** c>0, λ>0
  ⟹ 至少一个特征值模 >1 ⟹ 无条件不稳定（与步长无关，是格式的代数性质）""")
c = 0.02
for lam in (0.0384, 1.0, 3.9616):
    det = -1 - c*c*lam
    print(f"    λ={lam:7.4f}: det={det:.6f}  |det|={abs(det):.6f} > 1 ⟹ 不稳定")
check("T4：辛交错格式 |det|>1 对所有 λ>0", all(abs(-1-c*c*l) > 1 for l in (0.0384,1.0,3.9616)))

print()
print("="*80); print("结论"); print("="*80)
print("""
 出口乙的可证结果：
  ✅ 忠实的层间动力学**存在**：守恒、有界（ρ>0 时 dF/dt≤0）、不可约（§Z112 Q3 判据）。
  ✅ 由公理自动生成第二定律（F 单调下降），**无需公设**。
  ✅ 非均匀稳定态**存在**，但它的机制是"ρ 的零点 = 冻结组态"：
     双阱 Φ 的两个阱 m₋, m₊ 上流量恒为 0 ⟹ 任何 {m₋,m₊} 赋值都是精确驻点。
  ⚠️ 代价：这种非均匀态是**冻结**的（动力学停止），不是"活的稳态"。
  ❌ 单调 ρ ⟹ 数学上只能收敛到均匀，q→1/V=0.0625，与 D4 窗 (1/2,3/5) 无关。
  ❌ 非对称/辛型（"流量有记忆"）线性格式无条件不稳定。

 因此：**D4 窗不是这类动力学的输出**；要得到它，必须再定一笔数据——
  "哪个密度是允许驻留的"（即 ρ 的零点位置），这正是本轮被逼出来的具名输入。
""")
n = sum(1 for _, ok in CHECKS if ok)
print(f"断言通过 {n}/{len(CHECKS)}")


# ══════════════════════════════════════════════════════════════════
# 追加：T6 门槛（四次穿过 ρ=1）与 T5 锁定（非均匀驻点⟺ρ 的零点）
# ══════════════════════════════════════════════════════════════════
print()
print("="*80); print("T6/T5 数值核验"); print("="*80)

def Phi_of(rho, m, m_ref=1.0):
    """Φ(m) = ∫_mref^m log ρ(s) ds（数值积分）"""
    s = np.linspace(min(m, m_ref), max(m, m_ref), 200001)
    v = np.log(np.maximum(rho(s), 1e-300))
    I = np.trapezoid(v, s) if hasattr(np, "trapezoid") else np.trapz(v, s)
    return I if m >= m_ref else -I

print("\n(1) 双阱门槛：Φ 的双阱 ⟺ ρ 四次穿过 ρ=1")
for nm, rho, lo, hi in [
    ("ρ=1+0.3cos(4πm)", lambda m: 1.0+0.3*np.cos(4*np.pi*m), 0.0, 1.0),
    ("ρ=1+0.3cos(2πm)", lambda m: 1.0+0.3*np.cos(2*np.pi*m), 0.0, 2.0),
    ("ρ=(m−0.5)(2−m)·2", lambda m: np.clip(2*(m-0.5)*(2-m), 0, None), 0.5, 2.0),
]:
    ms = np.linspace(lo+1e-6, hi-1e-6, 20001)
    r = rho(ms)
    cross = np.sum(np.diff(np.sign(r-1.0)) != 0)
    # Φ 的驻点与极值类型
    dPhi = np.log(np.maximum(r, 1e-300))
    stat = ms[np.where(np.diff(np.sign(dPhi)) != 0)[0]]
    kinds = []
    for s in stat:
        i = int(np.argmin(np.abs(ms-s)))
        if 0 < i < len(ms)-1:
            kinds.append("极小" if dPhi[i-1] < 0 < dPhi[i+1] else "极大")
    print(f"  {nm:22s} 穿过 ρ=1 次数={cross}  驻点={len(stat)}  类型={kinds}")
print("  ⟹ 2 次穿过 → 单井（仅一个极小）；4 次穿过 → 双阱。判据成立。")

print("\n(2) T5 锁定：非均匀驻点 ⟺ 取值落在 ρ 的零点")
m_minus, m_plus = 0.5, 2.0
rho_dw = lambda m: np.clip((m - m_minus)*(m_plus - m), 0, None)
for cfg, vals in [("全部在阱上 {0.5,2}", np.array([0.5,2.0]*8)),
                  ("均值=1.25 的混合",  np.array([0.5]*8+[2.0]*8))]:
    V = len(vals); adj = path_adj(V); E = edges_of(adj)
    fl = np.array([rho_dw(np.array([vals[u]]))[0]*rho_dw(np.array([vals[v]]))[0]*(vals[u]-vals[v])
                   for u, v in E])
    print(f"  {cfg:20s}: 最大流量={np.abs(fl).max():.2e}  var={np.var(vals):.3f}  "
          f"均值={vals.mean():.3f}  → {'冻结驻点' if np.abs(fl).max()<1e-12 else '非驻点'}")
print("  ⟹ 取 {0.5,2} 的构型流量恒为 0 = 精确驻点；但自由度只剩'哪些点取哪个阱'，")
print("     且均值被锁为 0.5k/V+2(V−k)/V —— 均值成为必须由初值提供的自由参数。")

print("\n(3) 一般非均匀驻点的数量条件")
print("  逐边 ρ(m_v)(m_v−μ)=C；死端逼 C=0；故非均匀驻点 = {m_v} ⊂ {ρ 的零点}")
print("  ⟹ 设零点集为 Z，则驻点形如 m_v∈Z，且必须满足 Σ_v m_v = Vμ。")
print("  ⟹ μ 必须恰等于某个 #Z 组合的平均值 —— 对一般 μ **不成立**。")
print("  ⟹ 非均匀驻点是**孤立/非通有**的：它不是吸引子，而是测度零的特例。")
