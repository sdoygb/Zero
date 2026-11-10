#!/usr/bin/env python3
"""单光子干涉条纹：纯格点计算（Zb014 受迫稳态 + Zb018 读出计数 + Zb096 定理 4.1）

舞台：二维方格子（a=1），c=1（每步位移，Zb001）。两个"缝"= 两个受迫顶点。
方程（Zb014 定理 1.2 的稳态式，本征阻尼 γ）：
    (c² L − ω² + 2 i γ ω) U = g e_s ,   ω = 2c sin(q/2)  (Zb005 定理 2.2)
舞台边界是耗散端口（边缘 γ 增大），屏取在吸收层之前。
屏读数 N(y) = |U₁+U₂|²（Zb018 定理 2.2 的读出计数 + Zb096 定理 4.1）。
单光子：每次事件写一份整数记录，位置按 p(y) ∝ N(y)（Zb096 推论 4.3）。
判读条纹用两路相位差的斜率 dΔφ/dy（不找峰），并与 λL/d 对照。
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
import json, math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTD = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results")
OUTD.mkdir(parents=True, exist_ok=True)

def build(Nx, Ny, c=1.0, omega=0.693148, g0=4e-4, gmax=0.2, edge=30):
    n=Nx*Ny; idx=lambda x,y:x*Ny+y
    rows,cols,vals=[],[],[]
    for x in range(Nx):
        for y in range(Ny):
            i=idx(x,y); deg=0
            for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
                xx,yy=x+dx,y+dy
                if 0<=xx<Nx and 0<=yy<Ny:
                    rows.append(i);cols.append(idx(xx,yy));vals.append(-c*c);deg+=1
            rows.append(i);cols.append(i);vals.append(deg*c*c)
    L=sp.csr_matrix((vals,(rows,cols)),shape=(n,n))
    gam=np.full(n,g0)
    for x in range(Nx):
        for y in range(Ny):
            d=min(x,Nx-1-x,y,Ny-1-y)
            if d<edge: gam[idx(x,y)]+=gmax*((edge-d)/edge)**2
    M=L-omega**2*sp.identity(n,format='csr',dtype=complex)+2j*sp.diags(gam.astype(complex))
    return M, idx

def run(q, d, Nx=420, Ny=400, edge=32, xsrc=60):
    c=1.0; omega=2*c*math.sin(q/2)
    M, idx = build(Nx, Ny, omega=omega, edge=edge)
    ys=Ny//2; xs=Nx-edge-6; Ld=xs-xsrc
    Us=[]
    for sy in (ys-d//2, ys+d//2):
        rhs=np.zeros(Nx*Ny,dtype=complex); rhs[idx(xsrc,sy)]=1.0
        Us.append(spla.spsolve(M,rhs).reshape(Nx,Ny))
    return Us, xs, ys, Ld, omega

def main():
    c=1.0; q=0.707832
    lam=2*math.pi/q; lam_cont=2*math.pi*c/(2*c*math.sin(q/2))
    d=40
    Us, xs, ys, Ld, omega = run(q, d)
    y=np.arange(Ny_dim:=Us[0].shape[1])-ys
    U1,U2=Us[0][xs,:], Us[1][xs,:]
    I=np.abs(U1+U2)**2; I/=I.max()
    m=np.abs(y)<=110
    dphi=np.unwrap(np.angle(U1[m]*np.conj(U2[m])))
    slope=abs(np.polyfit(y[m],dphi,1)[0])
    period=2*math.pi/slope
    pred_lat=lam*Ld/d; pred_cont=lam_cont*Ld/d
    print(f"q={q:.6f} λ格={lam:.4f} λ连续={lam_cont:.4f} d={d} L={Ld}")
    print(f"相位差斜率={slope:.5f} ⟹ 条纹周期={period:.2f} 格")
    print(f"λ格·L/d={pred_lat:.2f}  λ连续·L/d={pred_cont:.2f}"
          f"  （数值比格上公式大 {period/pred_lat-1:+.1%}，比连续公式大 {period/pred_cont-1:+.1%}）")
    # 单光子：按 p ∝ I 抽整数记录
    p=I/I.sum()
    rng=np.random.default_rng(20261110)
    res={"q":q,"lambda_lattice":lam,"lambda_continuum":lam_cont,"d":d,"L":Ld,
         "slope":slope,"period":period,"pred_lattice":pred_lat,"pred_continuum":pred_cont,
         "omega":omega,"single":{}}
    fig,axes=plt.subplots(4,1,figsize=(9,9),sharex=True)
    axes[0].plot(y[m], I[m], lw=0.8); axes[0].set_ylabel("N(y)  continuous solution")
    axes[0].set_title(f"Zb096 single-photon fringes, pure lattice solution  q={q:.4f}  lambda={lam:.2f} grid  d={d}  L={Ld}")
    for k,N in enumerate((200, 2000, 20000)):
        counts=rng.multinomial(N,p)
        res["single"][N]={"counts":counts[m].tolist()}
        axes[k+1].bar(y[m], counts[m], width=0.9)
        axes[k+1].set_ylabel(f"records, N={N}")
    axes[-1].set_xlabel("screen site y (grid units)")
    plt.tight_layout(); plt.savefig(OUTD/"zb096_single_photon_fringe.png", dpi=130)
    json.dump(res, open(OUTD/"zb096_single_photon_fringe.json","w"), indent=2)
    print("图与数据已落盘")

if __name__=="__main__":
    main()
