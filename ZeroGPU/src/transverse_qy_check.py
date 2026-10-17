#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
transverse_qy_check.py —— 横向分量 q_y：二维色散与对称轴上的 q_y=0
====================================================================
对象：Z152 §5.2 的开放项——三维/方格舞台的横向分量 q_y。

二维色散（分离）：λ(qx,qy)=4(sin²(qx/2)+sin²(qy/2))，ω²=c²λ。
源在 y=y0 对称 ⇒ 对称行上 ∂_yψ=0 ⇒ qy=0 ⇒ 纵向 qx=2arcsin(f/(2c)) 精确。

核验：
  A. 二维散射（源 (y0,x0)，二次吸收层）：对称行 qx 与 2arcsin(f/2c) 一致；
  B. 对称列 qy ≈ 0；
  C. 偏离对称行满足二维色散。
"""
from __future__ import annotations
import os
import json
import math
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "results"))
LOG = []


def P(*a):
    s = " ".join(str(x) for x in a)
    LOG.append(s)
    print(s)


def build_2d(ny, nx, c2=1.0):
    N = ny * nx
    rows, cols, vals = [], [], []
    def idx(y, x):
        return y * nx + x
    for y in range(ny):
        for x in range(nx):
            i = idx(y, x)
            deg = 0
            for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                yy, xx = y + dy, x + dx
                if 0 <= yy < ny and 0 <= xx < nx:
                    deg += 1
                    rows.append(i); cols.append(idx(yy, xx)); vals.append(-c2)
            rows.append(i); cols.append(i); vals.append(c2 * deg)
    return sp.csr_matrix((vals, (rows, cols)), shape=(N, N))


def absorber(ny, nx, m=24, sig0=1.5):
    yy, xx = np.mgrid[0:ny, 0:nx]
    de = np.minimum(np.minimum(yy, ny - 1 - yy), np.minimum(xx, nx - 1 - xx))
    return np.where(de < m, sig0 * ((m - de) / m) ** 2, 0.0)


def solve(ny, nx, q, c2=1.0, src=None, m=24, sig0=1.5):
    H = build_2d(ny, nx, c2)
    N = ny * nx
    E = 2 * c2 * (1 - math.cos(q))
    y0, x0 = (ny // 2, nx // 5) if src is None else src
    b = np.zeros(N, complex)
    b[y0 * nx + x0] = 1.0
    A = (E * sp.eye(N, format="csr") - H + 1j * sp.diags(absorber(ny, nx, m, sig0).ravel())).tocsc()
    return spla.splu(A).solve(b).reshape(ny, nx)


def phase_grad(v, i0, i1):
    seg = v[i0:i1]
    ph = np.unwrap(np.angle(seg))
    return float(np.polyfit(np.arange(len(seg)), ph, 1)[0])


def main():
    P("=" * 80)
    P("横向分量 q_y：二维色散与对称轴 (transverse_qy_check)")
    P("=" * 80)

    ny, nx = 120, 160
    c2 = 1.0
    f = math.log(2)
    qth = 2 * math.asin(f / 2)
    psi = solve(ny, nx, qth)
    y0, x0 = ny // 2, nx // 5
    qx = phase_grad(psi[y0, :], x0 + 15, nx - 30)
    qy = phase_grad(psi[:, x0 + 40], 30, ny - 30)
    P(f"q_th = 2arcsin(f/2c) = {qth:.6f}")
    P(f"对称行 qx(实测) = {qx:.6f}  相对差 = {abs(qx-qth)/qth*100:.4f}%")
    P(f"对称列 qy(实测) = {qy:.6f}")

    checks = []

    def chk(name, cond, detail=""):
        checks.append({"name": name, "pass": bool(cond), "detail": detail})
        P(("  [v] " if cond else "  [x] ") + name + (f"   {detail}" if detail else ""))

    chk("A  对称行 qx = 2arcsin(f/2c)（<2%）", abs(qx - qth) / qth < 2e-2,
        f"{abs(qx-qth)/qth*100:.3f}%")
    chk("B  对称列 qy ≈ 0（<0.1）", abs(qy) < 0.1, f"{qy:.4f}")

    npass_off = 0
    ntot = 0
    for dy in (8, 12, 16):
        y = y0 + dy
        qxo = phase_grad(psi[y, :], x0 + 15, nx - 30)
        qyo = phase_grad(psi[:, x0 + 50], y0 + 3, min(ny - 3, y0 + 3 + 2 * dy))
        lam = 4 * (math.sin(qxo / 2) ** 2 + math.sin(qyo / 2) ** 2)
        om2 = c2 * lam
        ntot += 1
        if abs(om2 - f * f) / (f * f) < 0.10:
            npass_off += 1
        P(f"  dy={dy}: qx={qxo:.4f} qy={qyo:.4f} c²λ={om2:.5f} f²={f*f:.5f}")
    chk("C  偏离行满足二维色散（<10%，≥2/3）", npass_off >= 2, f"{npass_off}/{ntot}")

    npass = sum(1 for c in checks if c["pass"])
    P("\n" + "=" * 80)
    P(f"核验：通过 {npass} / {len(checks)}")
    P("=" * 80)

    out = {"q_th": qth, "qx_sym": qx, "qy_sym": qy,
           "off": npass_off, "checks": checks, "pass": npass, "total": len(checks)}
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "transverse_qy_check.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(RES, "transverse_qy_check.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    P("\n已存 transverse_qy_check.json / .txt")


if __name__ == "__main__":
    main()
