#!/usr/bin/env python3
"""闭合那 8%：格点上用 eikonal（驻相点）算相位，与格点解逐点对账

eikonal：点源到屏点的相位 Φ(Δ) = q·Δ，其中 q=(qx,qy) 满足
   sin qx = λ·Δx/r,  sin qy = λ·Δy/r,  4[sin²(qx/2)+sin²(qy/2)] = ω²
（即群速度 ∇_qω ∝ (sin qx, sin qy) 平行于位移 Δ）。
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
from scipy.optimize import brentq
import math, json
from pathlib import Path

OUT = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb096_eikonal_check.json")

def eikonal_phase(Lx, Ly, q0, c=1.0):
    """从原点到 (Lx,Ly) 的相位（驻相点）"""
    om2 = 4*c*c*math.sin(q0/2)**2
    r = math.hypot(Lx, Ly)
    def f(lam):
        ax = min(0.999999, lam*Lx/r); ay = min(0.999999, lam*Ly/r)
        qx = math.asin(ax); qy = math.asin(ay)
        return 4*c*c*(math.sin(qx/2)**2 + math.sin(qy/2)**2) - om2
    lam = brentq(f, 1e-9, 0.999999)
    qx = math.asin(min(0.999999, lam*Lx/r)); qy = math.asin(min(0.999999, lam*Ly/r))
    return qx*Lx + qy*Ly, qx, qy

def solve_field(q, d, slit_w=7, Nx=460, Ny=760, edge=34, xsrc=70, g0=3e-4, gmax=0.25):
    c=1.0; om=2*c*math.sin(q/2)
    n=Nx*Ny; idx=lambda x,y:x*Ny+y
    Lx=sp.kron(sp.diags([1.,-2.,1.],[-1,0,1],shape=(Nx,Nx)), sp.identity(Ny))
    Ly=sp.kron(sp.identity(Nx), sp.diags([1.,-2.,1.],[-1,0,1],shape=(Ny,Ny)))
    Lap=(-(Lx+Ly)).tocsr()
    gx=np.arange(Nx)[:,None]; gy=np.arange(Ny)[None,:]
    dd=np.minimum(np.minimum(gx,Nx-1-gx), np.minimum(gy,Ny-1-gy))
    gam=g0+np.where(dd<edge, gmax*((edge-np.minimum(dd,edge))/edge)**2, 0.0)
    M=Lap-om**2*sp.identity(n,format='csr',dtype=complex)+2j*sp.diags(gam.ravel().astype(complex))
    ys=Ny//2; xs=Nx-edge-6; L=xs-xsrc
    U=[]
    for sy in (ys-d//2, ys+d//2):
        rhs=np.zeros(n,dtype=complex)
        for j in range(-(slit_w//2), slit_w//2+1): rhs[idx(xsrc,sy+j)]=1.0
        U.append(spla.spsolve(M.tocsc(),rhs).reshape(Nx,Ny)[xs,:])
    return U[0],U[1],L,ys

def main():
    q=0.707832; d=150
    U1,U2,L,ys=solve_field(q,d)
    y=np.arange(U1.size)-ys
    m=np.abs(y)<=110
    ph_num=np.unwrap(np.angle(U1[m]*np.conj(U2[m])))
    ph_eik=np.array([eikonal_phase(L, yy+d/2, q)[0]-eikonal_phase(L, yy-d/2, q)[0] for yy in y[m]])
    # 比较斜率（去掉常数）
    s_num=np.polyfit(y[m],ph_num,1)[0]
    s_eik=np.polyfit(y[m],ph_eik,1)[0]
    resid=ph_num-(ph_eik-ph_eik.mean()+ph_num.mean())
    print(f"q={q:.6f} λ=2π/q={2*math.pi/q:.4f} λ_eff=2π/sin q={2*math.pi/math.sin(q):.4f} d={d} L={L}")
    print(f"  数值相位差斜率 = {s_num:+.6f}")
    print(f"  eikonal 斜率   = {s_eik:+.6f}")
    print(f"  比 = {s_num/s_eik:.4f}   ⟹ 对应周期 数值 {2*math.pi/abs(s_num):.3f} / eikonal {2*math.pi/abs(s_eik):.3f}")
    print(f"  逐点残差 RMS = {np.std(resid):.4f} rad（量程 {np.ptp(ph_num):.1f} rad）")
    # 小角闭式
    closed=2*math.pi*L/(d*math.sin(q)); naive=2*math.pi*L/(d*q)
    print(f"  小角闭式 2πL/(d sin q) = {closed:.3f}；教科书 λL/d = {naive:.3f}；差 {closed/naive-1:+.2%}")
    # 亮纹位置：ΔΦ=2πn
    fringes=[]
    for n in range(-4,5):
        f=lambda yy: (eikonal_phase(L, yy+d/2, q)[0]-eikonal_phase(L, yy-d/2, q)[0]) - 2*math.pi*n
        try:
            fringes.append((n, brentq(f, -140, 140)))
        except Exception: pass
    print("  eikonal 亮纹位置:", [(n, round(v,2)) for n,v in fringes])
    json.dump(dict(q=q, lam=2*math.pi/q, lam_eff=2*math.pi/math.sin(q), d=d, L=L,
                   slope_num=s_num, slope_eikonal=s_eik, ratio=s_num/s_eik,
                   resid_rms=float(np.std(resid)), closed=closed, naive=naive,
                   fringes=[(n,float(v)) for n,v in fringes]), open(OUT,"w"), indent=2)
    print("落盘", OUT)

if __name__=="__main__":
    main()
