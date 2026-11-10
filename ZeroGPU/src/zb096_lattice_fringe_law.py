#!/usr/bin/env python3
"""判据 D3：格上条纹间距的闭式 Δy = 2πaL/(d sin q) 的数值核验

只用 Zero 的对象：二维格点 + Zb014 定理 1.2 的受迫稳态 + Zb018 的读出计数。
对每个 (q, d) 解两个受迫问题，取屏上一行，量两路相位差的斜率得条纹周期，
与三候选公式比较：
  naive_lattice = 2πL/(d q)          （用波数 q）
  naive_cont    = 2πcL/(d ω)         （用连续波长）
  closed        = 2πaL/(d sin q)     （Zb096 定理 2.4）
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
import math, json
from pathlib import Path

OUT = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb096_lattice_fringe_law.json")

def period(q, d, L_screen=410, Ny=420, edge=34, xsrc=70, g0=3e-4, gmax=0.25, c=1.0):
    om = 2*c*math.sin(q/2)
    Nx = xsrc + L_screen + edge + 6
    n = Nx*Ny; idx = lambda x, y: x*Ny + y
    Lx = sp.kron(sp.diags([1.,-2.,1.],[-1,0,1],shape=(Nx,Nx)), sp.identity(Ny))
    Ly = sp.kron(sp.identity(Nx), sp.diags([1.,-2.,1.],[-1,0,1],shape=(Ny,Ny)))
    Lap = (-(Lx+Ly)).tocsr()
    gx = np.arange(Nx)[:,None]; gy = np.arange(Ny)[None,:]
    dd = np.minimum(np.minimum(gx, Nx-1-gx), np.minimum(gy, Ny-1-gy))
    gam = g0 + np.where(dd < edge, gmax*((edge-np.minimum(dd,edge))/edge)**2, 0.0)
    M = Lap - om**2*sp.identity(n,format='csr',dtype=complex) + 2j*sp.diags(gam.ravel().astype(complex))
    ys = Ny//2; xs = Nx-edge-6; L = xs-xsrc
    U = []
    for sy in (ys-d//2, ys+d//2):
        rhs = np.zeros(n, dtype=complex); rhs[idx(xsrc,sy)] = 1.0
        U.append(spla.spsolve(M.tocsc(), rhs).reshape(Nx,Ny)[xs,:])
    y = np.arange(Ny)-ys
    m = np.abs(y) < min(70, Ny//2-edge-14)
    dphi = np.unwrap(np.angle(U[0][m]*np.conj(U[1][m])))
    slope = abs(np.polyfit(y[m], dphi, 1)[0])
    return 2*math.pi/slope, L

def main():
    rows = []
    print(f"{'q':>7} {'sin q':>7} {'d':>3} {'L':>4} {'measured':>9} {'2piL/(d q)':>11} {'2piL/(d sin q)':>15} {'ratio':>7}")
    for q in (0.35, 0.5, 0.6, 0.707832, 0.8, 0.9):
        for d in (30, 50):
            T, L = period(q, d)
            closed = 2*math.pi*L/(d*math.sin(q))
            naive = 2*math.pi*L/(d*q)
            rows.append(dict(q=q, sinq=math.sin(q), d=d, L=L, measured=T, closed=closed,
                             naive_lattice=naive, ratio=T/closed))
            print(f"{q:7.4f} {math.sin(q):7.4f} {d:3d} {L:4d} {T:9.2f} {naive:11.2f} {closed:15.2f} {T/closed:7.4f}")
    ratios = [r["ratio"] for r in rows]
    naive_ratios = [r["measured"]/r["naive_lattice"] for r in rows]
    print(f"\n闭式：比值 {min(ratios):.4f}–{max(ratios):.4f}（均值 {np.mean(ratios):.4f}，标准差 {np.std(ratios):.4f}）")
    print(f"用 q   ：比值 {min(naive_ratios):.4f}–{max(naive_ratios):.4f}（随 q 单调掉）")
    json.dump(rows, open(OUT,"w"), indent=2)
    print("落盘", OUT)

if __name__ == "__main__":
    main()
