#!/usr/bin/env python3
"""把条纹画出来：亮度全部由 Zero 的格点解算出来，不拟合任何照片。

链：Zb014 定理 1.2 的受迫稳态 (c²L − ω² + 2iγω)U = g e_s（每个缝是若干相邻受迫顶点）
   → 屏读数 N(y) = |U₁+U₂|²（Zb018 定理 2.2、Zb096 定理 4.1）
   → 单光子：按 p(y) ∝ N(y) 抽整数份记录（Zb096 推论 4.3）
亮度换算：显示用 I^γ（γ=0.5，模拟相机的非线性响应），γ 是显示选择，不改物理。
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
import math, json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['PingFang HK','Hiragino Sans GB','Arial Unicode MS','STHeiti']
plt.rcParams['axes.unicode_minus'] = False

OUTD = Path("/Users/oygb/Downloads/Zero/ZeroGPU/results")
OUTD.mkdir(parents=True, exist_ok=True)

def solve(q, d, slit_w=5, Nx=460, Ny=440, edge=34, xsrc=70, g0=3e-4, gmax=0.25):
    c = 1.0; om = 2*c*math.sin(q/2)
    n = Nx*Ny; idx = lambda x, y: x*Ny + y
    Lx = sp.kron(sp.diags([1.,-2.,1.],[-1,0,1],shape=(Nx,Nx)), sp.identity(Ny))
    Ly = sp.kron(sp.identity(Nx), sp.diags([1.,-2.,1.],[-1,0,1],shape=(Ny,Ny)))
    Lap = (-(Lx+Ly)).tocsr()
    gx = np.arange(Nx)[:,None]; gy = np.arange(Ny)[None,:]
    dd = np.minimum(np.minimum(gx, Nx-1-gx), np.minimum(gy, Ny-1-gy))
    gam = g0 + np.where(dd < edge, gmax*((edge-np.minimum(dd,edge))/edge)**2, 0.0)
    M = Lap - om**2*sp.identity(n,format='csr',dtype=complex) + 2j*sp.diags(gam.ravel().astype(complex))
    ys = Ny//2; xs = Nx-edge-6; L = xs-xsrc
    U = np.zeros((Nx,Ny), dtype=complex)
    for sy in (ys-d//2, ys+d//2):
        rhs = np.zeros(n, dtype=complex)
        for j in range(-(slit_w//2), slit_w//2+1):
            rhs[idx(xsrc, sy+j)] = 1.0
        U += spla.spsolve(M.tocsc(), rhs).reshape(Nx,Ny)
    return U, xs, ys, L, om

def main():
    q=0.707832; d=150; slit_w=7
    U, xs, ys, L, om = solve(q, d, slit_w, Nx=460, Ny=760)
    Nx, Ny = U.shape; edge = 34
    lam = 2*math.pi/q; lam_eff = 2*math.pi/math.sin(q)
    I = np.abs(U)**2
    # 屏：取一行
    Is = I[xs,:]; y = np.arange(I.shape[1]) - ys
    # 可见度（直接由算出的场量）
    m = np.abs(y) <= 110
    Ipk = Is[m].max()
    # 找相邻极值
    pk = [i for i in range(5, m.sum()-5) if Is[m][i]==Is[m][i-5:i+6].max()]
    vn = [i for i in range(5, m.sum()-5) if Is[m][i]==Is[m][i-5:i+6].min()]
    V = (Is[m][pk].mean()-Is[m][vn].mean())/(Is[m][pk].mean()+Is[m][vn].mean()) if pk and vn else float('nan')
    print(f"q={q:.6f} λ={lam:.3f} λ_eff=2π/sin q={lam_eff:.3f} d={d} 缝宽={slit_w} L={L}")
    print(f"屏上峰数={len(pk)}  暗谷数={len(vn)}  可见度 V={V:.3f}")
    print(f"定理2.4 预测条纹周期 = {2*math.pi*L/(d*math.sin(q)):.2f} 格")

    # ---- 图 1：全场的明暗（涟漪盆视角）+ 屏上的条纹（实验室视角）
    fig = plt.figure(figsize=(12.5, 5.2))
    ax1 = fig.add_subplot(1,2,1)
    x0 = 150
    ext = I[x0:Nx-edge, edge:Ny-edge]
    ext = ext/np.percentile(ext, 99.5)
    ax1.imshow(np.clip(ext.T,0,1)**0.5, origin='lower', cmap='inferno',
               extent=[x0, Nx-edge, edge-Ny//2, Ny-edge-Ny//2], aspect='auto')
    ax1.set_title("(a) 舞台上算出的光强 |U|^2（线性色标，γ=0.5 显示）", fontsize=10)
    ax1.set_xlabel("x（格）"); ax1.set_ylabel("y（格）")
    # 屏的实验室视角：竖条 = 屏，亮度 = 算出的强度
    ax2 = fig.add_subplot(1,2,2)
    strip = np.tile((Is[m]/Ipk)[:,None]**0.7, (1,60))
    ax2.imshow(strip, origin='lower', cmap='gray', aspect='auto',
               extent=[0, 1, y[m][0], y[m][-1]])
    ax2.set_title("(b) 屏上看到的条纹（亮度 = 算出的 |U1+U2|^2）", fontsize=10)
    ax2.set_xlabel("屏的一小段（沿缝方向，均匀）"); ax2.set_ylabel("y（格）")
    fig.suptitle(f"Zb096 双缝条纹：纯格点解  q={q:.4f}  λ={lam:.2f} 格  d={d}  缝宽 {slit_w}  L={L}", fontsize=11)
    plt.tight_layout()
    plt.savefig(OUTD/"zb096_fringe_lab_view.png", dpi=140)
    plt.close()

    # ---- 图 2：亮度的定量剖面
    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(y[m], Is[m]/Ipk, lw=1.0, color='k')
    ax.set_xlabel("y（格）"); ax.set_ylabel("归一化光强")
    ax.set_title(f"屏上光强剖面（算出来的，不是拟合）：可见度 V={V:.3f}")
    ax.grid(alpha=.3)
    plt.tight_layout(); plt.savefig(OUTD/"zb096_fringe_profile.png", dpi=140)
    plt.close()

    # ---- 图 3：单光子的整数记录（同一个屏）
    p = Is.copy(); p[~m]=0; p = p/p.sum()
    rng = np.random.default_rng(20261110)
    fig, axes = plt.subplots(1,3, figsize=(12,4.6), sharey=True)
    for k,(N,title) in enumerate(((200,"单光子 200 次"),(2000,"2000 次"),(20000,"20000 次"))):
        cnt = rng.multinomial(N, p)
        strip = np.tile(cnt[m][:,None].astype(float), (1,60))
        axes[k].imshow(strip, origin='lower', cmap='gray', aspect='auto',
                       extent=[0,1,y[m][0],y[m][-1]], vmax=max(1,np.percentile(cnt[m],99.5)))
        axes[k].set_title(title, fontsize=10); axes[k].set_xlabel("屏")
    axes[0].set_ylabel("y（格）")
    fig.suptitle("同一个屏，按算出的 p(y) ∝ |U1+U2|^2 逐份记录（Zb096 推论 4.3）", fontsize=11)
    plt.tight_layout(); plt.savefig(OUTD/"zb096_fringe_single_photon.png", dpi=140)
    plt.close()

    json.dump(dict(q=q, lam=lam, lam_eff=lam_eff, d=d, slit_w=slit_w, L=L,
                   period_pred=2*math.pi*L/(d*math.sin(q)), V=V, n_peaks=len(pk), n_troughs=len(vn),
                   y=y[m].tolist(), I=(Is[m]/Ipk).tolist()),
              open(OUTD/"zb096_fringe_image.json","w"), indent=2)
    print("图与数据已落盘")

if __name__ == "__main__":
    main()
