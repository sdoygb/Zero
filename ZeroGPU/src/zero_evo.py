# -*- coding: utf-8 -*-
"""
zero_evo.py —— Zero 多参数长期演化引擎（统一状态空间）
=====================================================

背景（为什么要有这个文件）
--------------------------
`ZeroGPU/src/` 里已有若干**忠实**引擎，但每个只带一个扇区：

  * `z0_engine.py`      走者计数（整数、闭合、终端、重播种）
  * `loyalty_engine.py` 层间闭环（top2 记忆、T=4 毁灭、稳态两分）
  * `d248_capacity_freeze_check.py`  内容 + 容量（Z139 冻结）
  * `z_integrate.py`    场/力回写（单一标量力）
  * `k4_quantum_*.py`   相位/波（单一频率或单一源）

它们各自把"舞台 + 记忆 + 内容 + 势 + 相位 + 几何 + 容量"切成互不相见的碎片，
所以每问一个数就要重写一次闭环。本文件把这些扇区放进**一个状态对象、一个步进器**，
让它们同时长期演化。

状态向量（全部在一个 `ZeroEvo` 里）
----------------------------------
  W[u,w]       L2 走者：起点 u、当前位置 w 的（浮点）重数
  seed[u]      L1→L2 的播种剖面（top2 记忆）
  record[u]    L1 累积闭合（按顶点）
  terminal[u]  终端层累积（容量事件写入）
  n[v]         内容密度（L1 宏观；均值 = 1 附近，logistic 有界）
  phi[v]       势：L phi = -(n - nbar)           （Z033 定理 16）
  psi[v]       复相位/波幅，i dpsi/dt = H psi    （Z152 色散；Z131 模流频率）
  h[v], hv[v]  几何标量扇区（度规扰动迹）与速度  （Z106/Z108/Z109 的标量代理）
  以及 t, gen, events, history

演化通道（每个通道都有台账出处；"探索"者显式标注）
------------------------------------------------
  C1 内容线性增长      dn/dt = mu (n - nbar)                 Z139 §1（x = q - qbar）
  C2 logistic 饱和     dn/dt += gamma n (1 - n/K)            探索（硬容量的光滑版）
  C3 势驱动交换        dn/dt += kappa_phi (n - nbar)         Z144 定理 144.2（phi = -L^+ x）
  C4 压强/扩散         dn/dt += -kappa_p L n                 Z139 §1
  C5 几何回写          dn/dt += -kappa_h L h                 Z093/Z109 的回写形式
  C6 泊松势            L phi = -(n - nbar)                   Z033 定理 16
  C7 几何波            d2h/dt2 = c_h^2 L h - Gamma hv + g_src (n-nbar)
                                                              Z106/Z108/Z109（标量扇区）
  C8 相位/波           i dpsi/dt = (c_w^2 L + omega0^2) psi + beta diag(n-nbar) psi
                       + g_wave (n-nbar)                     Z152/Z158/Z131（势耦合为探索）
  C9 层间闭环          n <- (1-alpha) n + alpha nbar (c/Σc)  Z129 边③
  C10 回写播种         seed <- seed * exp(eta_n x + eta_h h + eta_psi |psi|^2)
                                                              Z093/Z088 边②（乘法回写为探索）
  C11 容量             |n - nbar|_max > K_hard ⟹ 写终端/清空 L2/重播种
                                                              Z022 定理 4；Z139 §1(1.1)
  C12 周期毁灭         gen % T == 0 ⟹ W = 0                  Z017 定理 7；Z022

"探索"标记的含义：该通道的**形式**不是台账里的原式，而是为了让多个扇区能同时
演化所选的最小正则化（例如 logistic 是硬容量的光滑版）。它们全部有开关，默认 0；
任何结论必须注明开了哪些通道。

用法
----
    cd ZeroGPU/src && /usr/bin/python3 zero_evo.py --selftest
    from zero_evo import ZeroEvo, Params
    e = ZeroEvo(zs.cycle(16), Params(mu=1.0, kappa_p=0.5, capacity="hard"))
    e.run(20000, every=10)
    print(e.summary())

自检 20 条（13 组判据，见 `selftest()`）：与 `loyalty_engine` 逐位对照、Z139 冻结、Z148 比例继承、
泊松残差、交换守恒与两支谱、波保范与色散、几何能量守恒、容量记账。
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zero_sim as zs  # noqa: E402

RES = os.path.normpath(os.path.join(HERE, "..", "results"))


# ══════════════════════════════════════════════════════════════════
#  参数
# ══════════════════════════════════════════════════════════════════
@dataclass
class Params:
    """引擎参数。每个字段末尾注明台账出处或"探索"。"""

    # ── 舞台 / 记忆 ─────────────────────────────────────────────
    L: int = 4                  # 词长（Z021/Z022）
    T: int = 4                  # 毁灭周期，0 = 不毁灭（Z022 定理 4）
    retention: str = "top2"     # last | top2 | all（Z017 定理 7）
    micro: bool = True          # 是否跑 L2 走者扇区
    micro_semantics: str = "persistent"  # persistent（Z020 口径）| fresh（词层口径）

    # ── 内容（L1 宏观密度）──────────────────────────────────────
    mu: float = 0.0             # 线性增长（Z139 的 +x）
    gamma: float = 0.0          # logistic 增长（探索：硬容量的光滑版）
    K: float = 1.0              # logistic 承载（仅 gamma≠0 时用）
    kappa_phi: float = 0.0      # 势驱动交换（Z144；>0 聚拢）
    kappa_p: float = 0.0        # 压强/扩散（Z139 的 -κLx）
    kappa_h: float = 0.0        # 几何回写（Z093/Z109 形式）
    kappa_psi: float = 0.0      # 波强回写内容（探索：辐射压形式）

    # ── 几何标量扇区 ───────────────────────────────────────────
    g_src: float = 0.0          # 内容 -> 几何源
    c_h: float = 1.0            # 几何波速
    Gamma_h: float = 0.0        # 几何阻尼

    # ── 相位 / 波扇区 ──────────────────────────────────────────
    beta: float = 0.0           # 势 -> 相位（探索：局部红移形式）
    g_wave: float = 0.0         # 内容 -> 波源（Z158 的标量源）
    omega0: float = 0.0         # 波的质量/偏移
    c_w: float = 1.0            # 波速（Z152 的 c）
    Gamma_w: float = 0.0        # 波阻尼（探索：Z152 吸收边界的体积版）

    # ── 层间回写 ───────────────────────────────────────────────
    alpha: float = 0.0          # 闭合 -> 宏观密度（Z129 边③；1 = 完全替换）
    eta_n: float = 0.0          # 密度 -> 播种（探索：乘法回写）
    eta_h: float = 0.0          # 几何 -> 播种（探索）
    eta_psi: float = 0.0        # 波强 -> 播种（探索）

    # ── 容量 ───────────────────────────────────────────────────
    capacity: str = "none"      # none | soft | hard（Z022/Z139）
    K_hard: float = 1.0         # 硬阈值：max|n-nbar| > K_hard 触发
    reset_amp: float = 0.05     # 重播种幅度（Z139 的 amp0）

    # ── 数值 ───────────────────────────────────────────────────
    integrator: str = "euler"   # euler（台账口径，Z139/Z144）| rk4
    dt: float = 0.05            # 每个"代"推进的宏观时间
    substeps: int = 1           # 每代子步数（>=0；0 = 自动）
    norm_hi: float = 1e200      # 走者扇区的公共归一化阈值
    clip_lo: float = 1e-12      # 内容密度下限（保持正性）
    wave_on: bool = True        # 是否推进相位扇区（关掉可省时间）
    poisson_on: bool = True     # 是否每步解泊松势（关掉则 phi 不更新）


# ══════════════════════════════════════════════════════════════════
#  工具
# ══════════════════════════════════════════════════════════════════
def lap_matrix(adj) -> np.ndarray:
    """图 Laplacian L = D - A（半正定，行和为零）。"""
    V = len(adj)
    L = np.zeros((V, V))
    for i, nb in enumerate(adj):
        L[i, i] = len(nb)
        for j in nb:
            L[i, j] -= 1.0
    return L


def _gini(x: np.ndarray) -> float:
    x = np.sort(np.asarray(x, dtype=float))
    n = x.size
    s = x.sum()
    if n == 0 or s == 0:
        return 0.0
    idx = np.arange(1, n + 1)
    return float((2 * idx - n - 1) @ x / (n * s))


def _entropy(p: np.ndarray) -> float:
    p = np.asarray(p, dtype=float)
    p = p[p > 0]
    if p.size == 0:
        return 0.0
    p = p / p.sum()
    return float(-(p * np.log(p)).sum())


# ══════════════════════════════════════════════════════════════════
#  引擎
# ══════════════════════════════════════════════════════════════════
class ZeroEvo:
    """统一多参数演化引擎。

    一个 `step()` = 一个"代"（Z022 的周期单位）：
      ① 宏观扇区推进 dt（内容/几何/相位，按需自动子步）
      ② 容量检查（硬阈值）
      ③ L2 走者推进 L 步（播种 -> 走 -> 闭合 -> 记录 -> 存活）
      ④ 闭合回写宏观密度（边③）
      ⑤ t += dt
    """

    def __init__(
        self,
        adj,
        params: Params | None = None,
        seed=None,
        n0=None,
        psi0=None,
        h0=None,
        hv0=None,
        rng=None,
    ):
        self.adj = [tuple(sorted(a)) for a in adj]
        self.V = len(self.adj)
        self.p = params or Params()
        self.A = zs.adj_matrix(self.adj)
        self.Lap = lap_matrix(self.adj)
        self.A_L = None
        self.diag_AL = None
        if self.p.micro_semantics == "fresh":
            self.A_L = np.linalg.matrix_power(self.A, self.p.L)
            self.diag_AL = np.diag(self.A_L)
        self.Lp = np.linalg.pinv(self.Lap)
        self.rng = np.random.default_rng(20261018 if rng is None else rng)

        # 谱（用于 Fiedler 投影、CFL、波色散）
        w, Q = np.linalg.eigh(self.Lap)
        self.lap_eigs = w
        self.lap_vecs = Q
        nz = np.where(w > 1e-9)[0]
        self.fiedler = Q[:, nz[0]].copy() if nz.size else np.zeros(self.V)
        if nz.size:
            lmin = float(w[nz[0]])
            mask = (w > 1e-9) & (np.abs(w - lmin) < 1e-8)
            self.fiedler_sub = Q[:, mask].copy()
        else:
            self.fiedler_sub = np.zeros((self.V, 0))
        self.lam_max = float(w.max()) if w.size else 1.0
        self.lam_min_nz = float(w[nz[0]]) if nz.size else 1.0

        # ── 走者扇区 ──
        self.W = np.zeros((self.V, self.V), dtype=float)
        self.fresh_surv = np.zeros(self.V)
        if seed is None:
            seed = np.ones(self.V)
        elif isinstance(seed, str) and seed == "geom":
            seed = np.array([2.0 ** (self.V - 1 - i) for i in range(self.V)])
        self.seed = np.asarray(seed, dtype=float).copy()
        if self.seed.sum() <= 0:
            self.seed = np.ones(self.V)
        self.seed0 = self.seed.copy()
        self.record = np.zeros(self.V)
        self.terminal = np.zeros(self.V)
        self.layers: list[np.ndarray] = []
        self.gen_closures = np.zeros(self.V)
        self.log_scale = 0.0

        # ── 内容扇区 ──
        if n0 is None:
            x = self.rng.normal(size=self.V)
            x -= x.mean()
            x *= 0.05 / max(1e-12, np.abs(x).max())
            n0 = 1.0 + x
        self.n = np.asarray(n0, dtype=float).copy()
        self.phi = self.poisson(self.n)

        # ── 相位扇区 ──
        if psi0 is None:
            psi0 = self.rng.normal(size=self.V) + 1j * self.rng.normal(size=self.V)
        self.psi = np.asarray(psi0, dtype=complex).copy()
        sn = np.linalg.norm(self.psi)
        if sn > 0:
            self.psi /= sn
        H0 = self.p.c_w ** 2 * self.Lap + self.p.omega0 ** 2 * np.eye(self.V)
        self.wave_eigs, self.wave_vecs = np.linalg.eigh(H0)

        # ── 几何扇区 ──
        self.h = np.zeros(self.V) if h0 is None else np.asarray(h0, dtype=float).copy()
        self.hv = np.zeros(self.V) if hv0 is None else np.asarray(hv0, dtype=float).copy()

        # ── 记账 ──
        self.t = 0.0
        self.gen = 0
        self.events: list[dict[str, Any]] = []
        self.reset_patterns: list[np.ndarray] = []
        self.mode_power = np.zeros(self.V)
        self.mode_power_n = 0
        self.terminal_ok = True
        self.hist: dict[str, list] = {
            "gen": [], "t": [], "n_mean": [], "n_std": [], "n_max": [], "n_min": [],
            "q": [], "entropy": [], "gini": [], "pr": [],
            "lock_cn": [], "lock_cs": [], "fiedler": [], "low": [], "high": [], "lam_peak": [],
            "phi_energy": [], "h_energy": [], "h_max": [],
            "psi_norm": [], "psi_coh": [], "top1": [], "top2": [],
            "events": [],
        }
        self._record()

    # ───────────────────────── 基本量 ─────────────────────────
    def poisson(self, x: np.ndarray) -> np.ndarray:
        """解 L phi = -(x - xbar)；返回零和 phi（Z033 定理 16）。"""
        x = np.asarray(x, dtype=float)
        z = x - x.mean()
        phi = -(self.Lp @ z)
        phi -= phi.mean()
        return phi

    def _mean_n(self) -> float:
        return float(self.n.mean())

    def _shape_n(self) -> np.ndarray:
        m = self._mean_n()
        return self.n / m if m != 0 else self.n.copy()

    def _shape_c(self) -> np.ndarray:
        s = self.gen_closures.sum()
        return self.gen_closures / s if s > 0 else np.ones(self.V) / self.V

    def _seed_used(self) -> np.ndarray:
        """播种剖面（含回写通道）。eta 全 0 时就是 top2 记忆。"""
        s = self.seed.copy()
        if s.sum() <= 0:
            s = np.ones(self.V)
        if self.p.eta_n == 0 and self.p.eta_h == 0 and self.p.eta_psi == 0:
            return s
        mn = max(self._mean_n(), 1e-30)
        x = self.n / mn - 1.0
        h = self.h - self.h.mean()
        a = np.abs(self.psi) ** 2
        ma = max(a.mean(), 1e-30)
        mod = self.p.eta_n * x + self.p.eta_h * h + self.p.eta_psi * (a / ma - 1.0)
        mod = np.clip(mod, -5.0, 5.0)
        return s * np.exp(mod)

    # ───────────────────────── 微观：L2 走者 ─────────────────────────
    def _micro_gen(self):
        """一代 L2 走者：播种、走 L 步、闭合、记录、存活、top2 重播种、周期毁灭。"""
        if self.p.micro_semantics == "fresh":
            self._micro_gen_fresh()
            return
        V = self.V
        diag = np.arange(V)
        self.W[diag, diag] += self._seed_used()
        for _ in range(self.p.L):
            self.W = self.W @ self.A
        c = np.diag(self.W).copy()
        self.gen_closures = c
        self.record = self.record + c
        self.W[diag, diag] = 0.0

        self.layers.append(c.copy())
        if len(self.layers) > 2:
            self.layers.pop(0)
        if self.p.retention == "last":
            ret = c.copy()
        elif self.p.retention == "top2":
            ret = c + (self.layers[-2] if len(self.layers) >= 2 else 0.0)
        elif self.p.retention == "all":
            ret = self.record.copy()
        else:
            raise ValueError(f"未知 retention {self.p.retention}")
        self.seed = ret
        self.gen += 1
        if self.p.T and self.gen % self.p.T == 0:
            self.W[:] = 0.0

        # 走者扇区一起归一化（保持 seed 与 W 的相对注入幅度）；
        # L1 历史（record/terminal）另归一化（否则 record 会被 W 的指数增长压到 0）。
        m_w = max(self.seed.sum(), self.W.sum())
        if m_w > self.p.norm_hi:
            self.seed /= m_w
            self.W /= m_w
            self.layers = [lay / m_w for lay in self.layers]
            self.log_scale += math.log10(m_w)
        m_r = max(self.record.sum(), self.terminal.sum())
        if m_r > self.p.norm_hi:
            self.record /= m_r
            self.terminal /= m_r

    def _micro_gen_fresh(self):
        """词层口径：每代走者从**当前位置**起一条新的长度 L 词。

        与 persistent 的差别：闭合词的长度严格为 L；上一代未闭合的走者
        在下一代按**当前顶点**重新起词（不再沿用原起点），故 L1 高度层
        b_h 与「长度 L 平衡词」逐位对应。
        """
        V = self.V
        s_start = self.fresh_surv + self._seed_used()
        F = s_start[:, None] * self.A_L          # 起点 u -> 终点 w 的流量
        clo = s_start * self.diag_AL             # 长度 L 的闭合（回到起点）
        self.gen_closures = clo
        self.record = self.record + clo
        self.fresh_surv = F.sum(axis=0) - clo    # 未闭合者按终点重标号
        self.W = F.copy()
        self.layers.append(clo.copy())
        if len(self.layers) > 2:
            self.layers.pop(0)
        if self.p.retention == "last":
            ret = clo.copy()
        elif self.p.retention == "top2":
            ret = clo + (self.layers[-2] if len(self.layers) >= 2 else 0.0)
        elif self.p.retention == "all":
            ret = self.record.copy()
        else:
            raise ValueError(f"未知 retention {self.p.retention}")
        self.seed = ret
        self.gen += 1
        if self.p.T and self.gen % self.p.T == 0:
            self.fresh_surv[:] = 0.0
        m_w = max(self.seed.sum(), self.fresh_surv.sum(), self.W.sum())
        if m_w > self.p.norm_hi:
            self.seed /= m_w
            self.fresh_surv /= m_w
            self.W /= m_w
            self.layers = [lay / m_w for lay in self.layers]
            self.log_scale += math.log10(m_w)
        m_r = max(self.record.sum(), self.terminal.sum())
        if m_r > self.p.norm_hi:
            self.record /= m_r
            self.terminal /= m_r

    # ───────────────────────── 宏观：内容 ─────────────────────────
    def _content_rhs(self, n: np.ndarray, h: np.ndarray) -> np.ndarray:
        x = n - n.mean()
        rhs = np.zeros_like(n)
        if self.p.mu:
            rhs += self.p.mu * x                    # Z139 的 +x
        if self.p.gamma:
            rhs += self.p.gamma * n * (1.0 - n / self.p.K)   # 探索：光滑容量
        if self.p.kappa_phi:
            rhs += self.p.kappa_phi * x             # Z144：-κ L phi = +κ x
        if self.p.kappa_p:
            rhs += -self.p.kappa_p * (self.Lap @ n)  # Z139 压强
        if self.p.kappa_h:
            rhs += -self.p.kappa_h * (self.Lap @ h)  # 几何回写
        if self.p.kappa_psi:
            a = np.abs(self.psi) ** 2
            rhs += self.p.kappa_psi * (a - a.mean())  # 探索：辐射压/波强回写
        return rhs

    def _content_step(self, dt: float):
        n, h = self.n, self.h
        if self.p.integrator == "euler":
            # 台账口径：Z139 的 x <- x + dt(x - kappa L x)、Z144 的 q <- q + dt(-L phi)
            self.n = np.maximum(n + dt * self._content_rhs(n, h), self.p.clip_lo)
            return
        k1 = self._content_rhs(n, h)
        k2 = self._content_rhs(np.maximum(n + 0.5 * dt * k1, self.p.clip_lo), h)
        k3 = self._content_rhs(np.maximum(n + 0.5 * dt * k2, self.p.clip_lo), h)
        k4 = self._content_rhs(np.maximum(n + dt * k3, self.p.clip_lo), h)
        self.n = np.maximum(n + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4), self.p.clip_lo)

    # ───────────────────────── 宏观：几何 ─────────────────────────
    def _geom_step(self, dt: float):
        if self.p.g_src == 0 and self.p.Gamma_h == 0 and self.p.c_h == 0:
            return
        if self.p.g_src == 0 and not np.any(self.h) and not np.any(self.hv):
            return
        h, hv = self.h, self.hv
        x = self.n - self.n.mean()
        c2 = self.p.c_h ** 2
        # 图上波方程：h'' = -c^2 L h + ...（L 半正定；Z152 的 omega^2 = c^2 lambda）
        acc = -c2 * (self.Lap @ h) - self.p.Gamma_h * hv + self.p.g_src * x
        hv = hv + 0.5 * dt * acc
        h = h + dt * hv
        acc = -c2 * (self.Lap @ h) - self.p.Gamma_h * hv + self.p.g_src * x
        hv = hv + 0.5 * dt * acc
        self.h, self.hv = h, hv

    # ───────────────────────── 宏观：相位/波 ─────────────────────────
    def _wave_step(self, dt: float):
        # 自由演化始终执行（H0 = c_w^2 L + omega0^2 为正，酉流稳定）
        if not self.p.wave_on:
            return
        if self.p.beta:
            v = self.p.beta * (self.n - self.n.mean())
            self.psi = self.psi * np.exp(-0.5j * dt * v)
        self.psi = self.wave_vecs @ (np.exp(-1j * self.wave_eigs * dt) *
                                     (self.wave_vecs.T @ self.psi))
        if self.p.beta:
            v = self.p.beta * (self.n - self.n.mean())
            self.psi = self.psi * np.exp(-0.5j * dt * v)
        if self.p.g_wave:
            self.psi = self.psi + (-1j * dt * self.p.g_wave) * (self.n - self.n.mean())
        if self.p.Gamma_w:
            self.psi = self.psi * math.exp(-self.p.Gamma_w * dt)

    # ───────────────────────── 子步数 ─────────────────────────
    def _nsub(self) -> int:
        if self.p.substeps > 0:
            return int(self.p.substeps)
        rate = 0.0
        rate = max(rate, abs(self.p.mu) + abs(self.p.kappa_phi) + 2 * abs(self.p.gamma))
        rate = max(rate, abs(self.p.kappa_p) * self.lam_max)
        rate = max(rate, abs(self.p.kappa_h) * self.lam_max)
        rate = max(rate, abs(self.p.c_h) * math.sqrt(self.lam_max))
        emax = abs(self.p.c_w) ** 2 * self.lam_max + abs(self.p.omega0) ** 2
        emax += abs(self.p.beta) * max(1.0, float(np.abs(self.n - self.n.mean()).max()))
        rate = max(rate, emax)
        if self.p.kappa_psi:
            rate = max(rate, abs(self.p.kappa_psi) *
                       max(1.0, float((np.abs(self.psi) ** 2).max())))
        return int(max(1, min(400, math.ceil(self.p.dt * rate / 0.4))))

    # ───────────────────────── 宏观推进 ─────────────────────────
    def _macro_advance(self):
        nsub = self._nsub()
        dt = self.p.dt / nsub
        for _ in range(nsub):
            self._content_step(dt)
            self._geom_step(dt)
            self._wave_step(dt)
            if self.p.poisson_on:
                self.phi = self.poisson(self.n)

    # ───────────────────────── 容量 ─────────────────────────
    def _reset_pattern(self) -> np.ndarray:
        if self.terminal.sum() > 0:
            p = self.terminal / self.terminal.sum()
        elif self.record.sum() > 0:
            p = self.record / self.record.sum()
        else:
            p = np.ones(self.V) / self.V
        return p

    def _capacity(self):
        if self.p.capacity == "none":
            return
        x = self.n - self.n.mean()
        xmax = float(np.abs(x).max())
        if self.p.capacity == "soft":
            # 光滑版：只是把越界部分压回（不写终端、不清空）
            if xmax > self.p.K_hard:
                self.n = self.n.mean() + x * (self.p.K_hard / xmax)
            return
        if xmax <= self.p.K_hard:
            return
        # 硬容量（Z022 §277 / Z139 §1(1.1)）：写终端、清空 L2、重播种
        p = self.n / self.n.sum()
        self.terminal = self.terminal + p
        size = float(np.clip(np.abs(x) - self.p.K_hard, 0.0, None).sum())
        self.events.append({"gen": self.gen, "t": self.t, "size": size,
                            "max_abs_x": xmax})
        # 重播种图案 = 最近两个事件图案之和（Z139 的 last2 口径）
        pat = p - p.mean()
        npn = np.linalg.norm(pat)
        if npn > 0:
            pat = pat / npn
        self.reset_patterns.append(pat.copy())
        if len(self.reset_patterns) > 2:
            self.reset_patterns.pop(0)
        z = np.sum(self.reset_patterns, axis=0)
        nz = np.linalg.norm(z)
        if nz > 0:
            z = z / nz
        self.n = self.n.mean() * (1.0 + self.p.reset_amp * z)
        self.n = np.maximum(self.n, self.p.clip_lo)
        self.W[:] = 0.0
        self.fresh_surv[:] = 0.0
        self.layers = []
        self.seed = np.ones(self.V) * (1.0 + self.p.reset_amp * z)
        self.seed = np.maximum(self.seed, 1e-12)

    # ───────────────────────── 边③ ─────────────────────────
    def _closure_update(self):
        a = self.p.alpha
        if a <= 0:
            return
        c = self.gen_closures
        s = c.sum()
        if s <= 0:
            return
        pc = c / s                      # 概率向量（均值 1/V）
        self.n = (1.0 - a) * self.n + a * (self.n.mean() * self.V) * pc
        self.n = np.maximum(self.n, self.p.clip_lo)

    # ───────────────────────── 一步 ─────────────────────────
    def step(self):
        self._macro_advance()
        self._capacity()
        if self.p.micro:
            self._micro_gen()
            self._closure_update()
        else:
            self.gen += 1
        self.t += self.p.dt
        return self

    # ───────────────────────── 观测 ─────────────────────────
    def snapshot(self) -> dict[str, float]:
        n = self.n
        mn = float(n.mean())
        p = n / n.sum() if n.sum() > 0 else np.ones(self.V) / self.V
        x = n - mn
        xnorm = float(np.linalg.norm(x))
        pc = self._shape_c()
        ps = self.seed / self.seed.sum() if self.seed.sum() > 0 else np.ones(self.V) / self.V
        pr = float(n.sum() ** 2 / max(1e-30, (n * n).sum()))
        h = self.h
        h_energy = 0.5 * float(self.hv @ self.hv)
        h_energy += 0.25 * self.p.c_h ** 2 * float(np.sum(self.A * (h[:, None] - h[None, :]) ** 2))
        psi_abs = np.abs(self.psi)
        psi_norm = float(np.linalg.norm(self.psi))
        psi_coh = float(abs(self.psi.sum()) / max(1e-30, psi_abs.sum()))
        if xnorm > 0 and self.fiedler_sub.size:
            fid = float(np.linalg.norm(self.fiedler_sub.T @ x) / xnorm)
        else:
            fid = 0.0
        if self.mode_power_n > 0:
            P = self.mode_power / self.mode_power.sum() if self.mode_power.sum() > 0 else self.mode_power
        else:
            P = np.zeros(self.V)
        low = float(P[self.lap_eigs < 1.0].sum())
        high = float(P[self.lap_eigs > 2.5].sum())
        order = np.argsort(-p)
        return {
            "n_mean": mn,
            "n_std": float(n.std()),
            "n_max": float(n.max()),
            "n_min": float(n.min()),
            "q": float(p @ p),
            "entropy": _entropy(p),
            "gini": _gini(n),
            "pr": pr,
            "lock_cn": float(np.linalg.norm(pc - p)),
            "lock_cs": float(np.linalg.norm(pc - ps)),
            "fiedler": fid,
            "low": low,
            "high": high,
            "lam_peak": float(P @ self.lap_eigs) if P.sum() > 0 else 0.0,
            "phi_energy": float(-0.5 * x @ self.phi),
            "h_energy": h_energy,
            "h_max": float(np.abs(h).max()),
            "psi_norm": psi_norm,
            "psi_coh": psi_coh,
            "top1": int(order[0]),
            "top2": int(order[1]) if self.V > 1 else int(order[0]),
        }

    def _record(self):
        x = self.n - self.n.mean()
        c = self.lap_vecs.T @ x
        self.mode_power += c * c
        self.mode_power_n += 1
        s = self.snapshot()
        for k in ("gen", "t"):
            self.hist[k].append(self.gen if k == "gen" else self.t)
        for k, v in s.items():
            self.hist[k].append(v)
        self.hist["events"].append(len(self.events))

    def run(self, gens: int, every: int = 1, verbose: bool = False):
        every = max(1, int(every))
        t0 = time.time()
        for i in range(gens):
            self.step()
            if (i + 1) % every == 0:
                self._record()
            if verbose and (i + 1) % max(1, gens // 10) == 0:
                print(f"  gen {self.gen} t={self.t:.3f} max={self.n.max():.4f} "
                      f"q={self.hist['q'][-1]:.5f} ev={len(self.events)} "
                      f"[{time.time()-t0:.1f}s]", flush=True)
        return self

    # ───────────────────────── 谱/诊断 ─────────────────────────
    def dominant_frequency(self, key: str = "n_std", frac: float = 0.5) -> dict[str, float]:
        """对历史标量做 FFT，返回主频（每单位宏观时间）与前三个峰。"""
        y = np.asarray(self.hist[key], dtype=float)
        if y.size < 8:
            return {"freq": 0.0, "amp": 0.0, "ratio": 0.0}
        y = y - y.mean()
        # 线性去趋势
        tt = np.arange(y.size)
        y = y - np.polyval(np.polyfit(tt, y, 1), tt)
        n = y.size
        Y = np.abs(np.fft.rfft(y * np.hanning(n)))
        # 采样间隔：t 的差分中位数
        tg = np.asarray(self.hist["t"], dtype=float)
        if tg.size >= 2:
            dts = np.diff(tg)
            dt_sample = float(np.median(dts))
        else:
            dt_sample = self.p.dt
        freqs = np.fft.rfftfreq(n, d=dt_sample)
        k = int(np.argmax(Y[1:]) + 1) if n > 2 else 0
        med = float(np.median(Y[1:])) if n > 2 else 1.0
        return {"freq": float(freqs[k]), "amp": float(Y[k]),
                "ratio": float(Y[k] / max(1e-30, med)),
                "period": float(1.0 / freqs[k]) if freqs[k] > 0 else 0.0}

    def lyapunov(self, gens: int = 2000, eps: float = 1e-7, renorm_every: int = 20) -> float:
        """双轨迹最大李雅普诺夫指数（只扰动内容扇区）。"""
        other = self.clone()
        other.n = other.n + eps * np.cos(np.arange(self.V) * 0.7)
        other.n = np.maximum(other.n, self.p.clip_lo)
        d0 = float(np.linalg.norm(other.n - self.n))
        logs = []
        for i in range(gens):
            self.step()
            other.step()
            d = float(np.linalg.norm(other.n - self.n))
            if (i + 1) % renorm_every == 0:
                logs.append(math.log(max(d, 1e-300) / max(d0, 1e-300)))
                # 重新对齐
                other.n = self.n + (other.n - self.n) * (eps / max(d, 1e-300))
                d0 = eps
        if not logs:
            return 0.0
        return float(np.mean(logs) / (renorm_every * self.p.dt))

    def clone(self) -> "ZeroEvo":
        e = ZeroEvo.__new__(ZeroEvo)
        e.__dict__ = {}
        for k, v in self.__dict__.items():
            if isinstance(v, np.ndarray):
                e.__dict__[k] = v.copy()
            elif isinstance(v, list):
                e.__dict__[k] = [x.copy() if isinstance(x, np.ndarray) else x for x in v]
            elif isinstance(v, dict):
                e.__dict__[k] = {kk: (list(vv) if isinstance(vv, list) else vv)
                                 for kk, vv in v.items()}
            elif isinstance(v, Params):
                e.__dict__[k] = Params(**asdict(v))
            else:
                e.__dict__[k] = v
        return e

    def summary(self) -> dict[str, Any]:
        s = self.snapshot()
        s.update({"V": self.V, "params": asdict(self.p), "events": len(self.events),
                  "log_scale": self.log_scale})
        return s

    # ───────────────────────── 存/取 ─────────────────────────
    def save(self, path: str):
        arrays = {k: v for k, v in self.__dict__.items() if isinstance(v, np.ndarray)}
        arrays["layers"] = np.array(self.layers) if self.layers else np.zeros((0, self.V))
        arrays["events_json"] = np.array(json.dumps(self.events, ensure_ascii=False))
        arrays["params_json"] = np.array(json.dumps(asdict(self.p), ensure_ascii=False))
        arrays["hist_json"] = np.array(json.dumps(self.hist, ensure_ascii=False))
        arrays["meta_json"] = np.array(json.dumps(
            {"t": self.t, "gen": self.gen, "log_scale": self.log_scale,
             "adj": [list(a) for a in self.adj]}, ensure_ascii=False))
        np.savez_compressed(path, **arrays)

    @classmethod
    def load(cls, path: str):
        z = np.load(path, allow_pickle=False)
        meta = json.loads(str(z["meta_json"]))
        p = Params(**json.loads(str(z["params_json"])))
        e = cls(meta["adj"], p)
        for k in e.__dict__:
            if k in z and isinstance(e.__dict__[k], np.ndarray):
                e.__dict__[k] = z[k].copy()
        e.t = meta["t"]
        e.gen = meta["gen"]
        e.log_scale = meta["log_scale"]
        e.events = json.loads(str(z["events_json"]))
        e.hist = json.loads(str(z["hist_json"]))
        e.layers = [lay.copy() for lay in z["layers"]] if z["layers"].size else []
        e.seed0 = e.seed.copy()
        return e


# ══════════════════════════════════════════════════════════════════
#  自检
# ══════════════════════════════════════════════════════════════════
CK: list[dict[str, Any]] = []


def ck(name: str, cond: bool, detail: str = "") -> bool:
    CK.append({"name": name, "pass": bool(cond), "detail": detail})
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ""))
    return bool(cond)


def _loyalty_compare():
    """与 loyalty_engine 逐位对照（相同 seed/参数/代数）。"""
    from loyalty_engine import LoyaltyEngine
    from z0_engine import path_adj
    V = 16
    adj = path_adj(V)
    seed = np.array([2.0 ** (V - 1 - i) for i in range(V)])
    le = LoyaltyEngine(adj, L=4, seed=seed, T=4, retention="top2")
    ev = ZeroEvo(adj, Params(L=4, T=4, retention="top2", micro=True), seed=seed)
    for _ in range(120):
        le.run_generation()
        ev.step()
    # 相对比（两个引擎的归一化口径不同，比形状）
    r_le = le.record / le.record.sum()
    r_ev = ev.record / ev.record.sum()
    c_le = le.gen_closures / le.gen_closures.sum()
    c_ev = ev.gen_closures / ev.gen_closures.sum()
    s_le = le.seed / le.seed.sum()
    s_ev = ev.seed / ev.seed.sum()
    d1 = float(np.abs(r_le - r_ev).max())
    d2 = float(np.abs(c_le - c_ev).max())
    d3 = float(np.abs(s_le - s_ev).max())
    ck("V1 与 loyalty_engine 逐位对照（record/closure/seed 形状）",
       max(d1, d2, d3) < 1e-12, f"max diff = {max(d1,d2,d3):.3e}")


def _z0_compare():
    """与 z0_engine（整数、无毁灭）对照。"""
    from z0_engine import ZeroCore, path_adj
    V = 8
    adj = path_adj(V)
    core = ZeroCore(adj, L=4, retention="top2")
    ev = ZeroEvo(adj, Params(L=4, T=0, retention="top2", micro=True),
                 seed=np.ones(V))
    # 只比 5 代：z0_engine 的 seed 做 gcd 约简（整数性能口径），会改变
    # fresh 与幸存走者的相对幅度；比值口径前 5 代尚未触发，可逐位比。
    for _ in range(5):
        core.run_generation()
        ev.step()
    r0 = np.array(core.record, dtype=float)
    r0 /= r0.sum()
    r1 = ev.record / ev.record.sum()
    d = float(np.abs(r0 - r1).max())
    ck("V2 与 z0_engine（整数核心）对照 5 代 record 形状", d < 1e-12,
       f"max diff = {d:.3e}")


def _poisson_check():
    rng = np.random.default_rng(7)
    for adj, name in [(zs.path(16), "path16"), (zs.cycle(16), "cycle16"),
                      (zs.k4(), "K4"), (zs.heawood(), "Heawood")]:
        e = ZeroEvo(adj, Params(micro=False))
        x = rng.normal(size=e.V)
        x -= x.mean()
        phi = e.poisson(x)
        res = float(np.abs(e.Lap @ phi + x).max())
        ck(f"V3 泊松残差 {name}", res < 1e-11, f"max|L phi + x| = {res:.3e}")


def _exchange_check():
    """Z144 定理 144.2/144.3：交换精确守恒；吸引支 λ=1+dt。"""
    adj = zs.cycle(8)
    e = ZeroEvo(adj, Params(mu=0.0, kappa_phi=1.0, kappa_p=0.0, micro=False,
                           dt=0.01, substeps=1),
                n0=1.0 + 0.05 * np.cos(2 * np.pi * np.arange(8) / 8))
    x0 = e.n - e.n.mean()
    e.run(50, every=1)
    x1 = e.n - e.n.mean()
    total_drift = abs(e.n.sum() - 8.0)
    growth = float(np.linalg.norm(x1) / np.linalg.norm(x0))
    pred = (1.0 + 1.0 * 0.01) ** 50
    ck("V4 交换型更新精确守恒（Z144）", total_drift < 1e-12,
       f"|Σn-8| = {total_drift:.3e}")
    ck("V5 吸引支增长因子 λ=1+dt（Z144 λ+）", abs(growth / pred - 1) < 1e-10,
       f"实测 {growth:.10f} vs 预测 {pred:.10f}")


def _capacity_freeze():
    """复现 Z139：容量把指数增长截断成周期，图案落 Fiedler。"""
    from z0_engine import torus_adj
    for adj, name in [(zs.cycle(16), "cycle16"), (torus_adj(2, 4), "torus4x4")]:
        lam = np.linalg.eigvalsh(lap_matrix(adj))
        nz = lam[lam > 1e-9]
        lmin, lmax = nz.min(), nz.max()
        kappa = 0.5 * (1.0 / lmax + 1.0 / lmin)
        e = ZeroEvo(adj, Params(mu=1.0, kappa_p=kappa, capacity="hard",
                               K_hard=1.0, reset_amp=0.05, micro=False,
                               dt=0.02, substeps=1),
                    n0=1.0 + 0.05 * np.cos(2 * np.pi * np.arange(len(adj)) / len(adj)))
        e.run(15000, every=10)
        # 冻结：末 20 个事件图案的相邻差
        pats = []
        # 用事件间末态近似：直接看末 200 步的 n 形状变化
        snapshots = []
        for _ in range(3000):
            e.step()
            if len(e.events) and (len(snapshots) == 0 or
                                  e.events[-1]["gen"] != getattr(e, "_last_ev", -1)):
                snapshots.append(e.n.copy())
                e._last_ev = e.events[-1]["gen"]
        if len(snapshots) >= 3:
            diffs = [float(np.linalg.norm(snapshots[i + 1] - snapshots[i]))
                     for i in range(len(snapshots) - 1)]
            last = diffs[-1]
        else:
            last = float("nan")
        x = e.n - e.n.mean()
        fid = float(np.linalg.norm(e.fiedler_sub.T @ x) / max(1e-30, np.linalg.norm(x)))
        ck(f"V6 容量冻结复现 {name}", (len(snapshots) >= 3 and last < 1e-5),
           f"末相邻图案差 = {last:.3e}（事件数 {len(e.events)}）")
        ck(f"V7 冻结图案落 Fiedler {name}", fid > 0.9999, f"投影 = {fid:.6f}")


def _ratio_inheritance():
    """复现 Z148：比例被继承（path16，5 个 seed）。"""
    from z0_engine import path_adj
    V = 16
    adj = path_adj(V)
    seeds = {
        "geom2": np.array([2.0 ** (V - 1 - i) for i in range(V)]),
        "geom1.5": np.array([1.5 ** (V - 1 - i) for i in range(V)]),
        "linear": np.array([float(V - i) for i in range(V)]),
        "quad": np.array([float((V - i) ** 2) for i in range(V)]),
        "uniform": np.ones(V),
    }
    ok = True
    detail = []
    for name, s in seeds.items():
        e = ZeroEvo(adj, Params(L=4, T=4, retention="top2", micro=True, alpha=0.0),
                    seed=s)
        e.run(4000, every=100)
        w = e.record / e.record.sum()
        ratio = w[7] / w[8]
        target = s[7] / s[8]
        rel = abs(ratio / target - 1)
        detail.append(f"{name}:{ratio:.5f}/{target:.5f}")
        ok = ok and rel < 1e-4
    ck("V8 比例继承复现（Z148）", ok, "; ".join(detail))


def _wave_check():
    """波扇区：保范、色散（单模相位频率 = 本征值）。"""
    adj = zs.cycle(16)
    e = ZeroEvo(adj, Params(micro=False, c_w=0.5, omega0=0.3, dt=0.01),
                psi0=np.ones(16, dtype=complex) / 4.0)
    # 取一个本征模
    k = 3
    v = np.exp(2j * np.pi * k * np.arange(16) / 16)
    v /= np.linalg.norm(v)
    e.psi = v.copy()
    # 直接算 H0 v 的 Rayleigh 商
    H0 = e.p.c_w ** 2 * e.Lap + e.p.omega0 ** 2 * np.eye(16)
    E = float(np.real(np.vdot(v, H0 @ v)))
    n0 = np.linalg.norm(e.psi)
    e.run(200, every=50)
    n1 = np.linalg.norm(e.psi)
    ck("V9 波扇区保范（酉演化）", abs(n1 / n0 - 1) < 1e-11,
       f"|d|psi||/|psi| = {abs(n1/n0-1):.3e}")
    # 相位频率：psi(t) ~ e^{-iEt} v；比较相位
    phase = float(np.angle(np.vdot(v, e.psi)))
    pred = (-E * 200 * 0.01) % (2 * np.pi)
    dphi = abs((phase - pred + np.pi) % (2 * np.pi) - np.pi)
    ck("V10 波单模相位频率 = 本征值", dphi < 1e-8,
       f"实测 {phase:.10f} vs 预测 {pred:.10f}（Δ={dphi:.2e}）")


def _geometry_check():
    """几何扇区：无源无阻尼时能量守恒；单模频率 = c_h sqrt(λ)。"""
    adj = zs.cycle(16)
    e = ZeroEvo(adj, Params(micro=False, c_h=0.5, Gamma_h=0.0, g_src=0.0, dt=0.01),
                h0=np.cos(2 * np.pi * 2 * np.arange(16) / 16))
    e.hv = np.zeros(16)
    E0 = e.snapshot()["h_energy"]
    for _ in range(2000):
        e.step()
    E1 = e.snapshot()["h_energy"]
    ck("V11 几何扇区能量守恒（无源无阻尼）", abs(E1 / E0 - 1) < 1e-5,
       f"E0={E0:.8f} E1={E1:.8f} rel={abs(E1/E0-1):.3e}")
    # 单模频率：用 FFT 测（k=2 的余弦模，λ = 2 - 2cos(2π k/V)）
    xs = []
    for _ in range(20000):
        e.step()
        xs.append(float(e.h @ np.cos(2 * np.pi * 2 * np.arange(16) / 16)))
    xs = np.asarray(xs)
    Y = np.abs(np.fft.rfft((xs - xs.mean()) * np.hanning(xs.size)))
    f = np.fft.rfftfreq(xs.size, d=0.01)
    k = int(np.argmax(Y[1:]) + 1)
    lam2 = 2.0 - 2.0 * math.cos(2 * math.pi * 2 / 16)
    pred = 0.5 * math.sqrt(lam2) / (2 * math.pi)   # FFT 给的是 cycles/time
    ck("V12 几何单模色散 ω = c_h sqrt(λ)", abs(f[k] / pred - 1) < 0.03,
       f"实测 {f[k]:.5f} vs 预测 {pred:.5f}（分辨率 {f[1]:.4f}）")


def _lock_check():
    """基线：无耦合时 L1/L2 自指锁死（Z149）。"""
    adj = zs.cycle(8)
    e = ZeroEvo(adj, Params(micro=True, alpha=1.0, retention="last"),
                seed=np.arange(1, 9, dtype=float))
    e.run(200, every=50)
    s = e.snapshot()
    ck("V13 无耦合时 L1/L2 锁死（α=1, last）",
       s["lock_cn"] < 1e-12 and s["lock_cs"] < 1e-12,
       f"lock_cn={s['lock_cn']:.3e} lock_cs={s['lock_cs']:.3e}")


def _capacity_accounting():
    """容量记账：事件数×写入的图案归一化；n 始终为正。"""
    adj = zs.cycle(12)
    e = ZeroEvo(adj, Params(mu=1.0, kappa_p=0.2, capacity="hard", K_hard=1.0,
                           micro=False, dt=0.02, substeps=1),
                n0=1.0 + 0.1 * np.cos(2 * np.pi * np.arange(12) / 12))
    e.run(8000, every=10)
    ck("V14 容量事件后 n > 0 且有限", bool(np.all(e.n > 0) and np.all(np.isfinite(e.n))),
       f"n_min={e.n.min():.3e} events={len(e.events)}")
    ck("V15 容量事件至少发生一次", len(e.events) > 0, f"events={len(e.events)}")


def selftest(write_json: bool = True) -> int:
    print("=" * 78)
    print("zero_evo 自检：统一引擎与旧引擎/台账逐项对照")
    print("=" * 78)
    t0 = time.time()
    _loyalty_compare()
    _z0_compare()
    _poisson_check()
    _exchange_check()
    _capacity_freeze()
    _ratio_inheritance()
    _wave_check()
    _geometry_check()
    _lock_check()
    _capacity_accounting()
    n_fail = sum(1 for c in CK if not c["pass"])
    print("-" * 78)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-n_fail} / FAIL {n_fail}；"
          f"用时 {time.time()-t0:.1f}s")
    if write_json:
        os.makedirs(RES, exist_ok=True)
        fp = os.path.join(RES, "zero_evo_selftest.json")
        with open(fp, "w", encoding="utf-8") as f:
            json.dump({"checks": CK, "n_pass": len(CK) - n_fail, "n_fail": n_fail,
                       "wall_s": round(time.time() - t0, 2)}, f,
                      ensure_ascii=False, indent=1)
        print(f"结果写入 {fp}")
    return n_fail


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        raise SystemExit(1 if selftest() else 0)
    print(__doc__)
