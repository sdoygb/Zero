#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""modular_probe.py — 「模流测试粒子」第一版：三件可算的事。

背景（库内）：
  · Z005 定理 1–3：模流 σ_t(A)=ρ^{it}Aρ^{-it}；账本态是中心态 ⟹ 流平凡。
  · Z033 定理 17：位置算符/时间平移/两质量解集为空 ⟹ v(r) 不是现成对象。
  · 本脚本不引入新公理，只把"模流 + 局部探针"这一读法算到底，看它给什么。

三部分：
  A. 二能级态：模流频率 = log(λ₊/λ₋)（与 Hamiltonian 无关）——速度的"量子"。
  B. Gibbs 态 + 线性势：模流就是 H 演化，加速度 a=F/m* 严格线性，低加速度端无偏离。
  C. 中心态：σ_t=id，探针不动（Z005 定理 3 的数值对照）。

用法：python3 ZeroGPU/src/modular_probe.py [--json]
退出码：0 = 全部判据通过；1 = 有不符（脚本/公式错，不是物理判定）。
"""
import json
import math
import sys

import numpy as np

J = 1.0          # 跃迁幅
HBAR = 1.0


# ───────────────────────── A. 二能级：模流频率 ─────────────────────────

def modular_frequency(p, theta):
    """ρ = U diag(p,1-p) U†，U = exp(-i θ σ_y/2)。返回模流在 σ_x 上的频率。"""
    c, s = math.cos(theta / 2), math.sin(theta / 2)
    U = np.array([[c, -s], [s, c]], dtype=complex)
    rho = U @ np.diag([p, 1 - p]) @ U.conj().T
    # K = -log ρ（2x2 用特征分解）
    w, V = np.linalg.eigh(rho)
    K = V @ np.diag(-np.log(w)) @ V.conj().T
    # 数值：e^{iKt} σ_x e^{-iKt}，取 ⟨σ_x⟩(t) 的振荡频率
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    ts = np.linspace(0, 200, 20001)
    vals = []
    for t in ts:
        Ut = V @ np.diag(np.exp(1j * (-np.log(w)) * t)) @ V.conj().T
        A = Ut @ sx @ Ut.conj().T
        vals.append(float(np.real(np.trace(A @ sx))) / 2.0 + 0.0)
    vals = np.array(vals)
    # 用 FFT 峰值求频率（避开端点）
    dt = ts[1] - ts[0]
    nfft = len(vals) * 2000          # 零填充：频率分辨率 ~1.6e-4
    spec = np.abs(np.fft.rfft(vals - vals.mean(), n=nfft))
    freqs = np.fft.rfftfreq(nfft, dt)
    peak = freqs[np.argmax(spec[1:]) + 1]
    omega = 2.0 * math.pi * peak          # rfftfreq 给的是普通频率 f=ω/2π
    return float(omega), (float(np.log(p / (1 - p))), float(vals.std()))


# ───────────────────────── B. 一维链 + 线性势 ─────────────────────────

def tight_binding(N, F):
    H = np.zeros((N, N), dtype=complex)
    for x in range(N):
        H[x, x] = F * (x - N / 2)
        if x + 1 < N:
            H[x, x + 1] = -J
            H[x + 1, x] = -J
    return H


def gaussian_packet(N, x0, sigma):
    x = np.arange(N)
    psi = np.exp(-((x - x0) ** 2) / (4 * sigma ** 2)).astype(complex)
    return psi / np.linalg.norm(psi)


def accel_from_packet(F, N=256, x0=None, sigma=6.0, T=8.0):
    """Gibbs 态：模流 = e^{-iβHt}（β=1）；测 ⟨x⟩(t) 的二次项。"""
    x0 = x0 if x0 is not None else N // 3
    H = tight_binding(N, F)
    E, U = np.linalg.eigh(H)
    psi0 = gaussian_packet(N, x0, sigma)
    coeff = U.conj().T @ psi0
    ts = np.linspace(0, T, 401)
    xs = np.arange(N) - N / 2
    mean_x = []
    for t in ts:
        psi = U @ (np.exp(-1j * E * t) * coeff)
        mean_x.append(float(np.sum(xs * np.abs(psi) ** 2)))
    mean_x = np.array(mean_x)
    # 二次拟合 x(t)=x0+v0 t+½ a t²
    A = np.vstack([np.ones_like(ts), ts, 0.5 * ts ** 2]).T
    coef, *_ = np.linalg.lstsq(A, mean_x, rcond=None)
    return float(coef[2])


# ───────────────────────── 主程序 ─────────────────────────


# ───────────────── D. 空间非中心态：K 的程函（局域还是长程？）─────────────────

def packet_state_K(N, x1, x2, sigma, p, eps):
    """ρ = (1-ε)(p|ψ1><ψ1|+(1-p)|ψ2><ψ2|) + ε I/N；返回去迹的 K=-log ρ。"""
    x = np.arange(N)
    def g(c):
        v = np.exp(-((x - c) ** 2) / (4 * sigma ** 2)).astype(complex)
        return v / np.linalg.norm(v)
    a, b = g(x1), g(x2)
    rho = (1 - eps) * (p * np.outer(a, a.conj()) + (1 - p) * np.outer(b, b.conj())) \
        + eps * np.eye(N) / N
    w, V = np.linalg.eigh(rho)
    K = V @ np.diag(-np.log(w)) @ V.conj().T
    return K - (np.trace(K) / N) * np.eye(N)


def offdiag_decay(K, ds):
    N = K.shape[0]
    out = {}
    for d in ds:
        vals = [abs(K[i, i + d]) for i in range(N - d)]
        out[d] = float(np.max(vals)) if vals else 0.0
    return out

def main():
    out = {"checks": [], "A": [], "B": []}

    def chk(name, ok, detail):
        out["checks"].append({"name": name, "ok": bool(ok), "detail": detail})
        print(f"{'PASS' if ok else 'FAIL'}  {name} — {detail}")

    print("=" * 74)
    print("A. 二能级态：模流频率 = log(λ₊/λ₋)，与 Hamiltonian / 基无关")
    print("=" * 74)
    p = 2.0 / 3.0
    target = math.log(p / (1 - p))
    for theta in (0.0, 0.3, 0.8, 1.4):
        w_hat, (w_exact, amp) = modular_frequency(p, theta)
        out["A"].append({"p": p, "theta": theta, "omega_fit": w_hat,
                         "omega_exact": w_exact, "amp": amp})
        print(f"  p=2/3 θ={theta:.1f}  ω_fit={w_hat:.6f}  log(p/(1-p))={w_exact:.6f}  "
              f"|Δ|={abs(w_hat - w_exact):.2e}")
    devs = [abs(r["omega_fit"] - r["omega_exact"]) for r in out["A"]]
    chk("A1 频率 = log(p/(1-p))", max(devs) < 2e-2, f"max|Δ|={max(devs):.2e}")
    chk("A2 频率与基的取向 θ 无关",
        max(r["omega_fit"] for r in out["A"]) - min(r["omega_fit"] for r in out["A"]) < 2e-2,
        f"spread={max(r['omega_fit'] for r in out['A']) - min(r['omega_fit'] for r in out['A']):.2e}")

    print("\n" + "=" * 74)
    print("B. Gibbs 态 + 线性势：模流 = H 演化；a/F 应与 F 无关（m*=1/(2J) ⟹ a/F=-2）")
    print("=" * 74)
    for F in (1e-1, 1e-2, 1e-3, 1e-4):
        a = accel_from_packet(F, T=8.0 if F >= 1e-2 else 20.0)
        ratio = a / F
        out["B"].append({"F": F, "a": a, "ratio": ratio})
        print(f"  F={F:.0e}  a={a:+.6f}  a/F={ratio:.6f}  （偏离 -2：{ratio + 2:+.2e}）")
    weak = [r["ratio"] for r in out["B"] if r["F"] <= 1e-2]
    spread = max(weak) - min(weak)
    chk("B1 弱场 a/F 与 F 无关（无低加速度阈值）", spread < 5e-3, f"spread={spread:.2e}")
    chk("B2 a/F ≈ -2（由 m*=1/(2J) 标定；宽度效应 <1%）", all(abs(r + 2) < 0.02 for r in weak),
        f"max|a/F+2|={max(abs(r + 2) for r in weak):.2e}")
    chk("B3 强场（F=0.1）偏离仅来自能带/强场（不是低加速度）",
        abs(out["B"][0]["ratio"] + 2) > 0.02, f"a/F={out['B'][0]['ratio']:.3f}")

    print("\n" + "=" * 74)
    print("C. 中心态：K ∝ I ⟹ σ_t = id，探针不动（Z005 定理 3）")
    print("=" * 74)
    N = 64
    K = np.log(N) * np.eye(N, dtype=complex)      # ρ = I/N ⟹ K = log(N)·I
    K = K - (np.trace(K) / N) * np.eye(N)          # 去掉常数（模流只认 K 的差）
    x = np.diag(np.arange(N) - N / 2).astype(complex)
    comm = float(np.max(np.abs(K @ x - x @ K)))
    out["C"] = {"commutator": comm}
    print(f"  ‖[K,x]‖max = {comm:.2e}（中心态：0 ⟹ 无漂移）")
    chk("C1 中心态 ⟹ [K,x]=0（无运动）", comm < 1e-10, f"‖[K,x]‖={comm:.1e}")


    print("\n" + "=" * 74)
    print("D. 空间非中心态（两个分离的 Gaussian，权重 2/3–1/3）：K 的程函")
    print("=" * 74)
    N = 192
    Kd = packet_state_K(N, 64, 128, sigma=8.0, p=2 / 3, eps=1e-3)
    dec = offdiag_decay(Kd, [1, 2, 4, 8, 16, 32, 64])
    print("  |K_ij| 的最大值（按距离 d）：")
    for d, v in dec.items():
        print(f"    d={d:3d}  {v:.3e}")
    out["D"] = {"offdiag": dec}
    r8 = dec[16] / dec[8] if dec[8] else 0.0
    r16 = dec[32] / dec[16] if dec[16] else 0.0
    r32 = dec[64] / dec[32] if dec[32] else 0.0
    print(f"  比值 r(d)=K(2d)/K(d)：r8={r8:.3f} r16={r16:.3f} r32={r32:.3f}"
          f"（幂律 d^-n 应恒定；Gaussian 应递减）")
    chk("D1 K 的程函由波包宽度定、Gaussian 衰减（不是 1/r 幂律尾）",
        r32 < r8 / 10.0 and dec[64] / dec[1] < 1e-2,
        f"r8={r8:.2f} r16={r16:.2f} r32={r32:.3f}; K64/K1={dec[64] / dec[1]:.1e}")

    n_ok = sum(1 for c in out["checks"] if c["ok"])
    print(f"\n判据 {n_ok}/{len(out['checks'])} 通过")
    if "--json" in sys.argv:
        print(json.dumps(out, ensure_ascii=False, indent=1))
    return 0 if n_ok == len(out["checks"]) else 1


if __name__ == "__main__":
    sys.exit(main())
