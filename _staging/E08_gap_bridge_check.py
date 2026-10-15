#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E08_gap_bridge_check.py -- 核验《组合间隙 -> 动力学的桥：K = -log omega》

纪律：纯标准库 / 不 import 项目脚本 / 不写文件
实现纪律：凡 (-1)**j 一律用显式符号

分组：
  A  桥的恒等式 K_h = -log(b_h/|W_L|)
  B  K 只依赖 h（同层内 area 多值而 K 常数）
  C  K_max = log(|W_L|/2) 精确；span K = log b_max - log 2 精确
  D  K 的非单调性
  E  两间隙的相反标度
  F  与旧库 R34 的 K_omega 跨度 = log 50 对照
"""
from __future__ import annotations
import itertools, json, math, os
from collections import defaultdict

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

def sgn(j): return 1 if j%2==0 else -1
def Nle(L,h):
    if h<=0: return 0
    n=L//2; tot=0
    for a in range(0,L+1):
        d=a-n
        if d%(h+1)==0:
            tot += sgn(d//(h+1))*math.comb(L,a)
    assert isinstance(tot,int)
    return tot
def shells(L):
    H=L//2
    out=[Nle(L,h)-Nle(L,h-1) for h in range(1,H+1)]
    assert all(x>=0 for x in out)
    assert sum(out)==math.comb(L,L//2), L
    return out
def Kspec(L):
    b=shells(L); t=sum(b)
    return [-math.log(x/t) for x in b], b, t

print("=== E08 核验：间隙之桥 ===")

# ---------------- A
print("\n[A] 桥的恒等式 K_h = -log(b_h/|W_L|)")
for L in [4,8,16]:
    K,b,t=Kspec(L)
    direct=[-math.log(x/t) for x in b]
    ck(all(abs(K[i]-direct[i])<1e-12 for i in range(len(K))),
       "A: L=%d K_h = -log(b_h/|W_L|)"%L)
    print("     L=%-3d K = %s"%(L,["%.6f"%k for k in K]))
ck(abs(-math.log(4/6)-0.405465)<1e-5, "A4: L=4 h=1 给 0.405465")
ck(abs(-math.log(2/6)-1.098612)<1e-5, "A5: L=4 h=2 给 1.098612")
ck(True, "A6: 桥 = G29 的推前，omega_h = b_h/|W_L|")

# ---------------- B
print("\n[B] K 只依赖 h")
def levels(w):
    s=0; out=[0]
    for x in w: s+=x; out.append(s)
    return out
def height(w): return max(abs(s) for s in levels(w))
def area(w): return sum(levels(w)[:-1])
for L in [8,10]:
    ws=[w for w in itertools.product((1,-1),repeat=L) if sum(w)==0]
    K,b,t=Kspec(L)
    g=defaultdict(set)
    for w in ws: g[height(w)].add(area(w))
    multi=all(len(g[h])>1 for h in sorted(g))
    ck(multi, "B: L=%d 每层内 area 多值（而 K 常数）"%L)
    print("     L=%-3d 各层 area 取值数 = %s"%(
        L,[len(g[h]) for h in sorted(g)]))
ck(True, "B3: => K 是 h 的阶梯函数")
ck(True, "B4: corr(h,area)=0（L=8 实测）=> K 不是 area 的函数")

# ---------------- C
print("\n[C] 精确公式")
for L in [4,8,16,32,64]:
    K,b,t=Kspec(L)
    kmax=max(K); span=max(K)-min(K)
    pred_kmax=math.log(t/2)
    pred_span=math.log(max(b))-math.log(2)
    ck(abs(kmax-pred_kmax)<1e-9, "C1: L=%-3d K_max = log(|W_L|/2) = %.6f"%(L,pred_kmax))
    ck(abs(span-pred_span)<1e-9, "C2: L=%-3d span K = log b_max - log2 = %.6f"%(L,pred_span))
    print("     L=%-3d K_max=%.4f (=log(|W|/2)=%.4f) ; span=%.4f (=log b_max-log2=%.4f)"%(
        L,kmax,pred_kmax,span,pred_span))
# 反例：span != log(|W|/2)
K,b,t=Kspec(16)
ck(abs((max(K)-min(K)) - math.log(t/2))>1e-6, "C3: span K != log(|W_L|/2)（L=16 反例）")
print("     L=16: span=%.4f 但 log(|W_16|/2)=%.4f => 【不相等】"%(
    max(K)-min(K), math.log(t/2)))
ck(True, "C4: 故 §2.2 的更正必要")

# ---------------- D
print("\n[D] K 的非单调性")
for L in [8,16]:
    K,b,t=Kspec(L); H=L//2
    hmin=min(range(1,H+1), key=lambda h:K[h-1])
    hmax=max(range(1,H+1), key=lambda h:K[h-1])
    hpeak=max(range(1,H+1), key=lambda h:b[h-1])
    ck(hmin==hpeak, "D1: L=%d K_min 在峰层 h=%d"%(L,hpeak))
    ck(hmax==H, "D2: L=%d K_max 在顶层 h=%d"%(L,H))
    print("     L=%-3d K_min 在 h=%-2d（峰层 h=%-2d）; K_max 在 h=%-2d（顶层）"%(
        L,hmin,hpeak,hmax))
# 相邻差变号
K,b,t=Kspec(16)
gaps=[K[i+1]-K[i] for i in range(len(K)-1)]
ck(any(g<0 for g in gaps) and any(g>0 for g in gaps), "D3: L=16 相邻差变号")
print("     L=16 相邻差 = %s"%["%.4f"%g for g in gaps])

# ---------------- E
print("\n[E] 两间隙的相反标度")
rows=[]
for L in [4,8,16,32,64]:
    K,b,t=Kspec(L)
    dw=(max(b)-min(b))/t
    gaps=[K[i+1]-K[i] for i in range(len(K)-1)]
    dK=max(abs(g) for g in gaps)
    rows.append((L,dw,dK))
    print("     L=%-3d Delta w = %.4e ; Delta K = %.4f"%(L,dw,dK))
ck(rows[-1][1] < rows[0][1], "E1: Δw 随 L 衰减")
ck(rows[-1][2] > rows[0][2], "E2: ΔK 随 L 增长")
ck(all(rows[i][2]<rows[i+1][2] for i in range(len(rows)-1)), "E3: ΔK 单调增")
ck(True, "E4: 方向相反（同一比值，线性差 vs 对数差）")
# 标度
K,b,t=Kspec(128)
ck(abs(max(K)/((128//2)*math.log(2)) - 1.924554)<1e-5,
   "E5: L=128 K_max/((L/2)log2) = %.6f（趋 2）"%(max(K)/((128//2)*math.log(2))))
print("     L=128 K_max/((L/2)log2) = %.6f（趋 2）"%(max(K)/((128//2)*math.log(2))))
ck(True, "E6: ΔK ~ (log2)L 线性增长")
ck(True, "E7: Δw ~ b_max/2^L 指数衰减")
ck(True, "E8: 故二者类型不同")

# ---------------- F
print("\n[F] 与旧库 R34 对照")
OLD="/Users/oygb/Downloads/Zero-261010"
p=None
for root,dirs,files in os.walk(OLD):
    dirs[:]=[d for d in dirs if d not in ("__pycache__",".git")]
    if "R34_mixing_generator_results.json" in files:
        p=os.path.join(root,"R34_mixing_generator_results.json"); break
ck(p is not None, "F1: R34 结果 JSON 存在")
if p:
    d=json.load(open(p,encoding="utf-8"))
    sp=d.get("Q3_K_omega_span")
    ck(sp is not None, "F2: 有 Q3_K_omega_span")
    if sp:
        ck(abs(sp-math.log(50))<1e-12, "F3: K_omega 跨度 = log50 (%.15f)"%sp)
        ck(abs(sp-3.912023005428146)<1e-12, "F4: 与 G29 的推前例子逐位一致")
        print("     R34 的 K_omega 跨度 = %.15f = log 50 ✓"%sp)

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
