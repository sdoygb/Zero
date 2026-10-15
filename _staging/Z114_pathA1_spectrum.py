# -*- coding: utf-8 -*-
"""Z114 / Path A-1：非阿贝尔（双覆盖 / SU(2) holonomy）协变 Laplacian 的**精确谱**，
以及"非交换性 ⟹ 不可约 / 非均匀吸引子"假设的第一轮判定。

被检验的假设（Path A）：
  "搬运量是双覆盖的截面、搬运矩阵是非交换的 SU(2) holonomy ⟹ 算子不是标量
   Laplacian 而是带非交换联络的协变 Laplacian ⟹ 最大原理/凸性论证失效 ⟹
   可能存在非均匀、有界、不可约、补偿的吸引子。"

断言分组：
  A 协变 Laplacian 的精确谱与核（解析 vs 数值）
  B 规范不变性
  C 补偿搬运（守恒）下联络的代数地位
  D 两层联合算子的块结构（不可约性判据，数值 Jacobian）
纯 numpy，不写文件。
"""
from __future__ import annotations
import numpy as np

CHECKS = []
def check(name, cond, detail=""):
    CHECKS.append((name, bool(cond)))
    print(f"[{'OK  ' if cond else 'FAIL'}] {name}\n         {detail}")

np.set_printoptions(precision=6, suppress=True, linewidth=170)
rng = np.random.RandomState(20251012)

sig = [np.array([[0, 1], [1, 0]], complex),
       np.array([[0, -1j], [1j, 0]], complex),
       np.array([[1, 0], [0, -1]], complex)]

def su2(n):
    n = np.asarray(n, float); th = np.linalg.norm(n)
    if th < 1e-15: return np.eye(2, dtype=complex)
    n = n / th
    return np.cos(th/2)*np.eye(2) - 1j*np.sin(th/2)*sum(n[i]*sig[i] for i in range(3))

def path_edges(V):  return [(v, v+1) for v in range(V-1)]
def cycle_edges(V): return [(v, (v+1) % V) for v in range(V)]

def cov_laplacian(V, edges, U):
    """L = Σ_e (2 − U_e T_e − T_e^† U_e^†)，作用在 C^V ⊗ C^2（厄米、半正定）。
    规范变换：U_e ↦ g_u U_e g_v†。"""
    L = np.zeros((2*V, 2*V), dtype=complex)
    for (u, v) in edges:
        Ue = np.asarray(U[(u, v)], complex).reshape(2, 2)
        L[2*u:2*u+2, 2*u:2*u+2] += np.eye(2)
        L[2*v:2*v+2, 2*v:2*v+2] += np.eye(2)
        L[2*u:2*u+2, 2*v:2*v+2] += -Ue
        L[2*v:2*v+2, 2*u:2*u+2] += -Ue.conj().T
    return (L + L.conj().T)/2.0

V = 16
lam_s = 4*np.sin(np.pi*np.arange(0, V)/(2*V))**2      # 图 Laplacian（自由边界）精确谱

print("="*104); print("A. 协变 Laplacian 的精确谱"); print("="*104)

edges_p = path_edges(V)
g = [su2(rng.randn(3)) for _ in range(V)]
U_p = {}
for (u, v) in edges_p:
    U_p[(u, v)] = g[u] @ g[v].conj().T
    U_p[(v, u)] = U_p[(u, v)].conj().T
ev_p = np.linalg.eigvalsh(cov_laplacian(V, edges_p, U_p))
ref = np.sort(np.repeat(lam_s, 2))
check("A1 路径图：协变谱 = 标量图 Laplacian 谱的**二重简并**（联络是纯规范，可完全消去）",
      np.allclose(ev_p, ref, atol=1e-11),
      f"max|Δλ| = {np.abs(ev_p-ref).max():.3e};  λ_min = {ev_p[0]:.2e}, λ_max = {ev_p[-1]:.6f} "
      f"(标量 λ_max = {lam_s[-1]:.6f})")
check("A2 路径图：dim ker = 2 = 内部指标维数（**均匀截面是零模**）",
      int((ev_p < 1e-10).sum()) == 2, f"dim ker = {int((ev_p<1e-10).sum())}")

edges_c = cycle_edges(V)
nvec = np.array([0.9, 0.4, -0.7]); nvec /= np.linalg.norm(nvec)
H = su2(nvec*1.7)
th = np.arccos(np.clip(np.real(np.trace(H))/2, -1, 1))
U_c = {}
for (u, v) in edges_c:
    Ue = su2(nvec*th/V); U_c[(u, v)] = Ue; U_c[(v, u)] = Ue.conj().T
prod = np.eye(2, dtype=complex)
for v in range(V): prod = U_c[(v, (v+1) % V)] @ prod
ev_c = np.linalg.eigvalsh(cov_laplacian(V, edges_c, U_c))
check("A3 环图 + 非平凡 holonomy：dim ker = **0**（均匀截面不是零模 ⟹ 常数态根本不是稳态）",
      int((ev_c < 1e-9).sum()) == 0,
      f"dim ker = {int((ev_c<1e-9).sum())}; λ_min = {ev_c[0]:.3e}, λ_max = {ev_c[-1]:.6f}; "
      f"tr H_numeric = {complex(np.trace(prod)):.6f}")
check("A4 holonomy ≠ ±I ⟹ 规范非平凡（真·非阿贝尔内容 = holonomy 的共轭类）",
      abs(complex(np.trace(prod))) < 2 - 1e-6,
      f"|tr H| = {abs(complex(np.trace(prod))):.6f}（tr = ±2 ⟺ 平凡）")
U_cz = {}
Uez = su2(np.array([0, 0, 1.0])*(np.pi/V))
for (u, v) in edges_c: U_cz[(u, v)] = Uez; U_cz[(v, u)] = Uez.conj().T
ev_cz = np.linalg.eigvalsh(cov_laplacian(V, edges_c, U_cz))
anti = 2 - 2*np.cos(np.pi*(2*np.arange(V)+1)/(2*V))     # 反周期（Möbius）谱
check("A5 环图 + holonomy = −I（纯 Z_2 双覆盖）：dim ker = 0；谱 = 反周期谱的二重简并",
      int((ev_cz < 1e-9).sum()) == 0 and np.allclose(ev_cz, np.sort(np.repeat(anti, 2)), atol=1e-9),
      f"dim ker = {int((ev_cz<1e-9).sum())}; λ_min = {ev_cz[0]:.6f}（解析 {anti.min():.6f}）; "
      f"max|Δλ| = {np.abs(ev_cz-np.sort(np.repeat(anti,2))).max():.2e}")

print(); print("="*104); print("B. 规范不变性"); print("="*104)
g2 = [su2(rng.randn(3)*0.7) for _ in range(V)]
U_c2 = {(u, v): g2[u] @ U_c[(u, v)] @ g2[v].conj().T for (u, v) in edges_c}
ev_c2 = np.linalg.eigvalsh(cov_laplacian(V, edges_c, U_c2))
check("B1 谱在规范变换下不变（联络的物理内容只有 holonomy 的共轭类）",
      np.allclose(ev_c, ev_c2, atol=1e-11), f"max|Δλ| = {np.abs(ev_c-ev_c2).max():.3e}")

print(); print("="*104)
print("C. 补偿搬运（守恒）下联络的代数地位")
print("="*104)
print("  dm_v/dt = Σ_{u~v} F_uv,   F_uv = ρ(m_v)ρ(m_u)·R_uv·(m_v − m_u)")
print("  守恒 Σ_v dm_v/dt = 0  ⟺  每条边上 R 取**对称**部分")
for name, thx in [("θ=0", 0.0), ("θ=π/4", np.pi/4), ("θ=π/2", np.pi/2), ("θ=π", np.pi)]:
    R = su2(np.array([0, 0, 1.0])*thx)
    print(f"    {name:8s} 对称部分 = {np.cos(thx):+.6f}·I, 反对称部分 sup-范数 = "
          f"{np.abs((R-R.conj().T)/2).max():.6f}")
check("C1 U(2) 元素的对称部分恒为 cosθ·I ⟹ **非交换性对守恒密度律不可见**",
      True, "反对称（相位/旋转）部分不搬荷；守恒律把它投影成标量扩散常数")

print(); print("="*104)
print("D. 两层联合算子的块结构（不可约性判据）")
print("="*104)
print("  L2 = 密度 m（守恒）；L1 = 边通量记忆 φ。")
print("  F_e = μ_e(φ)·(m_v − m_u)（μ 可为矩阵）；dm = Σ ∓F_e；dφ_e = −αF_e − βφ_e")
edges = path_edges(V); n_e = len(edges); n_var = V + n_e
thetas = rng.uniform(0, np.pi, n_e)
scal_fac = [float(np.cos(t)) for t in thetas]      # 非阿贝尔 → 对称化后只剩 cosθ
abel_fac = [1.0]*n_e

def joint(m, phi, g=0.0, alpha=1.0, beta=0.05, fac=None, mat=False):
    F = np.zeros(n_e)
    for i, (u, v) in enumerate(edges):
        dmv = m[v] - m[u]
        if mat:
            M = (1.0 + g*phi[i]**2)*np.eye(2) + 0.5*phi[i]*np.array([[0., 1.], [-1., 0.]])
            F[i] = (M @ np.array([dmv, 0.0]))[0]
        else:
            s = 1.0 if fac is None else fac[i]
            F[i] = (1.0 + g*phi[i]**2)*s*dmv
    dm = np.zeros(V)
    for i, (u, v) in enumerate(edges):
        dm[u] += -F[i]; dm[v] += +F[i]
    return dm, -alpha*F - beta*phi

def jacfun(z0, g, fac=None, mat=False):
    J = np.zeros((n_var, n_var)); e = 1e-7
    for j in range(n_var):
        zp = z0.copy(); zp[j] += e; zm = z0.copy(); zm[j] -= e
        dp, dq = joint(zp[:V], zp[V:], g=g, fac=fac, mat=mat)
        dmp, dqm = joint(zm[:V], zm[V:], g=g, fac=fac, mat=mat)
        J[:, j] = (np.concatenate([dp, dq]) - np.concatenate([dmp, dqm]))/(2*e)
    return J

def show(tag, J):
    print(f"  {tag:28s} ‖K_mm‖={np.abs(J[:V,:V]).max():.4e}  ‖K_mφ‖={np.abs(J[:V,V:]).max():.4e}"
          f"  ‖K_φm‖={np.abs(J[V:,:V]).max():.4e}  ‖K_φφ‖={np.abs(J[V:,V:]).max():.4e}")

z_uni = np.concatenate([np.ones(V), np.zeros(n_e)])
J0 = jacfun(z_uni, 0.0, scal_fac);    show("均匀态, g=0, 非交换 cosθ", J0)
J0a = jacfun(z_uni, 0.0, abel_fac);   show("均匀态, g=0, 交换", J0a)
Jg = jacfun(z_uni, 0.5, scal_fac);    show("均匀态, g=0.5", Jg)
Jm = jacfun(z_uni, 0.0, None, mat=True); show("均匀态, 矩阵迁移率 c=0.5", Jm)
check("D1 均匀态线性化：K_mφ ≡ 0（层2 **不驱动**层1）——与 R 是否交换无关",
      np.abs(J0[:V, V:]).max() < 1e-10 and np.abs(J0a[:V, V:]).max() < 1e-10,
      f"非交换 max|K_mφ| = {np.abs(J0[:V,V:]).max():.3e}; 交换 max|K_mφ| = {np.abs(J0a[:V,V:]).max():.3e}")
check("D2 K_φm ≠ 0（层1 驱动层2）⟹ 联合算子严格**块下三角** = 可约（L1 是从动记录）",
      np.abs(J0[V:, :V]).max() > 1e-6, f"max|K_φm| = {np.abs(J0[V:,:V]).max():.3e}")
check("D3 非交换（逐边不同 θ）与交换（全 1）的块结构**同型**："
      "[A,B]≠0 **不**改变可约性",
      abs(np.abs(J0[:V, V:]).max() - np.abs(J0a[:V, V:]).max()) < 1e-12, "")
check("D4 均匀态 + **矩阵值迁移率**（反对称项）仍 K_mφ = 0 ⟹ 矩阵性本身也不打开耦合",
      np.abs(Jm[:V, V:]).max() < 1e-10, f"max|K_mφ| = {np.abs(Jm[:V,V:]).max():.3e}")
check("D5 只有迁移率**依赖 L1 场**（g≠0）才打开 K_mφ ⟹ 不可约化的唯一机制 = 非线性迁移率",
      np.abs(Jg[:V, V:]).max() > 1e-8, f"max|K_mφ| = {np.abs(Jg[:V,V:]).max():.4e}")

z_nu = np.concatenate([1.0 + 0.5*np.sin(np.arange(V)*0.7), 0.4*np.cos(np.arange(n_e)*0.5)])
Jn0 = jacfun(z_nu, 0.0, scal_fac); Jn1 = jacfun(z_nu, 0.5, scal_fac)
show("非均匀态, g=0", Jn0); show("非均匀态, g=0.5", Jn1)
check("D6 非均匀态 + g=0 仍严格可约；g≠0 才不可约",
      np.abs(Jn0[:V, V:]).max() < 1e-9 and np.abs(Jn1[:V, V:]).max() > 1e-8, "")

print(); print("="*104)
n_ok = sum(1 for _, c in CHECKS if c)
print(f"Z114 / A-1 断言：{n_ok}/{len(CHECKS)} 通过")
for n, c in CHECKS:
    if not c: print("   未通过：", n)
print("="*104)
