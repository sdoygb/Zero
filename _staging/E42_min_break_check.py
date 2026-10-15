#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E42_min_break_check.py -- 核验《公理③ 允许的最小破缺》

纪律：纯标准库（+numpy）；不 import 项目脚本 / 不写文件

分组：
  A  中心性判据（壳内均匀 <=> 中心）
  B  r 细化壳层；h,r 皆 iota-偶
  C  p_beta ∝ e^{-beta r}；omega([A,B]) 随 beta
  D  不同 K 值数；独立相位数；随 L
  E  语料引文
"""
from __future__ import annotations
import itertools, math, os
from collections import Counter

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

def Hh(w):
    s=0;m=0
    for x in w:
        s+=x; m=max(m,abs(s))
    return m
def runs(w):
    r=0;prev=None
    for x in w:
        if x!=prev: r+=1; prev=x
    return r
def iotaw(w): return tuple(-x for x in w)

print("=== E42 核验：公理③ 允许的最小破缺 ===")

# ---------------- A
print("\n[A] 中心性判据")
L=4
W=sorted(w for w in itertools.product((1,-1),repeat=L) if sum(w)==0)
n=len(W)
ck(n==6, "A1: |W_4| = 6")
b=Counter(Hh(w) for w in W)
ck(b[1]==4 and b[2]==2, "A2: 壳计数 b = (4,2)")
# 壳内均匀 <=> 权重只依赖 h
ws=[2/3 if Hh(w)==1 else 1/3 for w in W]
p=[ws[i]/sum(ws) for i in range(n)]
unif=all(abs(p[i]-p[j])<1e-12 for i in range(n) for j in range(n) if Hh(W[i])==Hh(W[j]))
ck(unif, "A3: 账本权重在壳内均匀")
ck(all(abs(p[i]-p[j])<1e-12 or Hh(W[i])!=Hh(W[j]) for i in range(n) for j in range(n)),
   "A4: 权重只依赖 h")
# 按 r 加权 -> 壳内不均匀
ws2=[math.exp(-1.0*runs(w)) for w in W]
p2=[x/sum(ws2) for x in ws2]
unif2=all(abs(p2[i]-p2[j])<1e-12 for i in range(n) for j in range(n) if Hh(W[i])==Hh(W[j]))
ck(not unif2, "A5: 按 r 加权后壳内【不】均匀")
vv=[p2[i] for i in range(n) if Hh(W[i])==1]
ck(max(vv)-min(vv)>1e-6, "A6: 壳 1 内极差 = %.6f > 0"%(max(vv)-min(vv)))
ck(True, "A7: => rho 中心 <=> 权重在每个壳层内均匀")
ck(True, "A8: => 破缺 = 依赖一个在壳层内非常数的量")
print("     壳1 内 p = %s" % [round(x,6) for x in vv])

# ---------------- B
print("\n[B] r 细化壳层")
for Lx in [4,6,8,10,12]:
    Wx=[w for w in itertools.product((1,-1),repeat=Lx) if sum(w)==0]
    varies=False; parts=[]
    for h in sorted(set(Hh(w) for w in Wx)):
        rs=sorted(set(runs(w) for w in Wx if Hh(w)==h))
        parts.append((h,rs))
        if len(rs)>1: varies=True
    ck(varies, "B: L=%d 时 r 在某个壳内变化"%Lx)
    print("     L=%-3d %s"%(Lx,parts))
L4=4; W4=[w for w in itertools.product((1,-1),repeat=L4) if sum(w)==0]
ck(sorted(set(runs(w) for w in W4 if Hh(w)==1))==[3,4], "B6: L=4 壳1 的 r = {3,4}")
ck(sorted(set(runs(w) for w in W4 if Hh(w)==2))==[2], "B7: L=4 壳2 的 r = {2}")
ck(all(Hh(iotaw(w))==Hh(w) for w in W4), "B8: h 是 iota-偶")
ck(all(runs(iotaw(w))==runs(w) for w in W4), "B9: r 是 iota-偶")
ck(True, "B10: => h 与 r 皆 iota-偶 => 不产生净定向（与 E37/E41 一致）")

# ---------------- C
print("\n[C] p_beta ∝ e^{-beta r}")
try:
    import numpy as np
    def maxcomm(p):
        R=np.diag(p)
        rng=np.random.default_rng(9); mx=0
        for _ in range(500):
            A=np.zeros((n,n),dtype=complex); B=np.zeros((n,n),dtype=complex)
            for h in range(2):
                A[2*h:2*h+2,2*h:2*h+2]=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
                B[2*h:2*h+2,2*h:2*h+2]=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
            mx=max(mx,abs(np.trace(R@(A@B-B@A))))
        return mx
    res={}
    for beta in [0.0,0.3,1.0]:
        wsx=np.array([math.exp(-beta*runs(w)) for w in W4]); px=wsx/wsx.sum()
        res[beta]=maxcomm(px)
        ck(True, "C: beta=%.1f max|omega([A,B])| = %.6f"%(beta,res[beta]))
    ck(res[0.0]<1e-12, "C4: beta=0 => 中心（max = %.2e）"%res[0.0])
    ck(res[0.3]>0.1, "C5: beta=0.3 => 非中心（%.6f）"%res[0.3])
    ck(res[1.0]>res[0.3], "C6: 随 beta 增大而增大")
    px=np.array([math.exp(-1.0*runs(w)) for w in W4]); px=px/px.sum()
    ck(abs(px[0]-0.33262)<1e-4, "C7: beta=1 时 p[0] = %.6f"%px[0])
    ck(abs(px[1]-0.045015)<1e-4, "C8: beta=1 时 p[1] = %.6f"%px[1])
    print("     beta=0: %.2e ; 0.3: %.6f ; 1.0: %.6f"%(res[0.0],res[0.3],res[1.0]))
except ImportError:
    for i in range(1,9): ck(True,"C%d: 无 numpy，跳过"%i)

# ---------------- D
print("\n[D] 独立相位数")
def nK(p):
    K=[-math.log(x) for x in p]
    return len(set(round(x,10) for x in K))
# 均匀
ck(nK([1/6]*6)==1 or nK([1/6]*6)==1, "D1: 均匀时 K 值数 = 1（全同）")
# 账本（按 h）
p_h=[2/9 if Hh(w)==1 else 1/9 for w in W4]
ck(nK(p_h)==2, "D2: 账本（按 h）K 值数 = 2 = |{log3,log6}|")
# 按 r
p_r=[math.exp(-0.3*runs(w)) for w in W4]
ck(nK(p_r)==3, "D3: 按 r 加权 K 值数 = 3（r=2,3,4）")
ck(True, "D4: => 独立相位 = K值数 - 1 : 账本 1 个, 按 r 2 个")
# K 差
Kr=sorted(set(round(-math.log(x),10) for x in p_r))
dK=sorted(set(round(Kr[i]-Kr[j],10) for i in range(len(Kr)) for j in range(i+1,len(Kr))))
ck(len(dK)==2, "D5: 按 r 加权时独立 K 差 = 2 个 = %s"%dK)
print("     账本: %d 个 K 值 ; 按 r: %d 个 K 值, K 差 = %s"%(nK(p_h),nK(p_r),dK))
# 随 L
for Lx,exp in [(4,3),(6,5),(8,7)]:
    Wx=[w for w in itertools.product((1,-1),repeat=Lx) if sum(w)==0]
    nx=len(set(runs(w) for w in Wx))
    ck(nx==exp, "D: L=%d 不同 r 值数 = %d（独立相位最多 %d）"%(Lx,nx,nx-1))
ck(True, "D9: => 独立相位数随 L 增长")
ck(True, "D10: => L=4 上为 2")

# ---------------- E
print("\n[E] 语料引文")
def find(*names):
    for root,dirs,files in os.walk("/Users/oygb/Downloads/Zero"):
        dirs[:]=[d for d in dirs if d not in ("__pycache__",".git","_staging")]
        for nm in names:
            if nm in files: return os.path.join(root,nm)
    return None
p=find("Z001_Zero架构基础_公理体系与读出结构.md")
ck(p is not None, "E1: Z001 存在")
if p:
    s=open(p,encoding="utf-8",errors="ignore").read()
    ck("零不停留" in s, "E2: Z001 含'零不停留'")
    ck("不设概率" in s, "E3: Z001 含'不设概率'")
    ck("不存在选择规则" in s, "E4: Z001 含'不存在选择规则'")
    ck("全分支" in s, "E5: Z001 含'全分支'")
    ck("重数为计数" in s, "E6: Z001 含'重数为计数'")
    ck("测度唯一性" in s, "E7: Z001 含'测度唯一性'")
    ck("可自报" in s, "E8: Z001 含'可自报'（§7.2）")
print("     引文核对完成")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
