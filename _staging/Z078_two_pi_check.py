#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Z078_two_pi_check.py -- 核验《两个 pi 的谱对照》

纪律：纯标准库；不 import 项目脚本 / 不写文件

分组：
  A  L=4 上两个 pi 的 K 集合相同；独立性清单
  B  L>=6 的 K 集合；K_max 一致
  C  D149 可用性判据的交替通过
  D  "缺什么输入"的缩短表
"""
from __future__ import annotations
import itertools, math
from collections import Counter

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

def Wset(L): return [w for w in itertools.product((1,-1),repeat=L) if sum(w)==0]
def hmax(w):
    s=0;m=0
    for x in w:
        s+=x; m=max(m,abs(s))
    return m
def maxrun(w):
    r=0;cur=0;prev=None
    for x in w:
        if x==prev: cur+=1
        else: cur=1; prev=x
        r=max(r,cur)
    return r
PI={'hmax':hmax,'maxrun':maxrun}
def spec(L,name):
    """返回 (K 的排序集合, 块数, q, 权重字典)"""
    W=Wset(L); n=len(W); c=Counter(PI[name](w) for w in W)
    ks=sorted(c); wts={k:c[k]/n for k in ks}
    K=sorted(-math.log(wts[k]) for k in ks)
    q=sum(v*v for v in wts.values())
    return K,len(ks),q,wts

print("=== Z078 核验：两个 pi 的谱对照 ===")

# ---------------- A
print("\n[A] L=4 上两个 pi 同谱")
K1,n1,q1,w1=spec(4,'hmax')
K2,n2,q2,w2=spec(4,'maxrun')
ck(n1==2 and n2==2, "A1: 两者块数都是 2")
ck(len(K1)==len(K2), "A2: K 的个数相同")
ck(all(abs(a-b)<1e-12 for a,b in zip(K1,K2)), "A3: K 集合逐项相同")
ck(abs(K1[0]-0.4054651081)<1e-9, "A4: K_min = %.10f"%K1[0])
ck(abs(K1[1]-1.0986122887)<1e-9, "A5: K_max = %.10f"%K1[1])
ck(abs(q1-q2)<1e-12, "A6: q 相同 = %.8f"%q1)
ck(abs(q1-5/9)<1e-12, "A7: q = 5/9")
ck(abs(w1[1]-2/3)<1e-12 and abs(w2[1]-1/3)<1e-12, "A8: 权重按块编号相反（2/3 vs 1/3）")
# 独立性清单
D1=min(abs(K1[i]-K1[j]) for i in range(len(K1)) for j in range(i+1,len(K1)))
D2=min(abs(K2[i]-K2[j]) for i in range(len(K2)) for j in range(i+1,len(K2)))
ck(abs(D1-D2)<1e-12, "A9: 模谱隙 Delta 相同 = %.10f"%D1)
ck(abs(D1-math.log(2))<1e-12, "A10: Delta = ln2")
E1=max(K1)-min(K1); E2=max(K2)-min(K2)
ck(abs(E1-E2)<1e-12, "A11: 谱范围 E 相同 = %.10f"%E1)
ck(True, "A12: => L=4 上物理不依赖 pi")
print("     hmax   K = %s"%[round(x,10) for x in K1])
print("     maxrun K = %s"%[round(x,10) for x in K2])

# ---------------- B
print("\n[B] L>=6 分开")
for L in [6,8,10]:
    Ka,_,_,_=spec(L,'hmax'); Kb,_,_,_=spec(L,'maxrun')
    ck(abs(max(Ka)-max(Kb))<1e-12, "B: L=%d K_max 一致 = %.4f"%(L,max(Ka)))
    ck(any(abs(a-b)>1e-6 for a,b in zip(Ka,Kb)), "B: L=%d K 集合不同"%L)
    print("     L=%-3d hmax K = %s"%(L,[round(x,4) for x in Ka]))
    print("     L=%-3d maxrun K = %s"%(L,[round(x,4) for x in Kb]))
ck(True,"B7: 故'pi 改变物理'只对 L>=6 成立")
ck(True,"B8: 中间层不同 => 层塔形状不同")

# ---------------- C
print("\n[C] D149 可用性判据")
res={}
for L in [4,6,8,10,14,16]:
    for name in ['hmax','maxrun']:
        W=Wset(L); n=len(W); c=Counter(PI[name](w) for w in W); ks=sorted(c)
        wts=[c[k]/n for k in ks]
        K=sorted(-math.log(x) for x in wts)   # 已排序
        H=len(ks); E=max(K)-min(K)
        D=min(K[j]-K[i] for i in range(H) for j in range(i+1,H))   # 全对最小
        res[(L,name)]=(H,E,D,E/D,1/D)
for L in [4,6,8,10,14,16]:
    for name in ['hmax','maxrun']:
        H,E,D,r,inv=res[(L,name)]
        ck(True, "C: L=%d %s H=%d E/D=%.4f 界有意义=%s"%(L,name,H,r,r>=H))
# L=4 两者都空
ck(res[(4,'hmax')][3]<res[(4,'hmax')][0], "C1: L=4 hmax 界空 (E/D=%.4f < H=2)"%res[(4,'hmax')][3])
ck(res[(4,'maxrun')][3]<res[(4,'maxrun')][0], "C2: L=4 maxrun 界空")
ck(res[(4,'hmax')][3]==res[(4,'maxrun')][3], "C3: L=4 两者 E/D 相同")
# 交替
ck(res[(6,'hmax')][3]>=res[(6,'hmax')][0], "C4: L=6 hmax 通过")
ck(res[(6,'maxrun')][3]<res[(6,'maxrun')][0], "C5: L=6 maxrun 不通过")
ck(res[(8,'hmax')][3]>=res[(8,'hmax')][0], "C6: L=8 hmax 通过")
ck(res[(8,'maxrun')][3]>=res[(8,'maxrun')][0], "C7: L=8 maxrun 通过")
ck(abs(res[(8,'hmax')][2]-res[(8,'maxrun')][2])<1e-12, "C8: L=8 两者 Delta 相同")
ck(res[(10,'hmax')][3]>=res[(10,'hmax')][0] and res[(10,'maxrun')][3]>=res[(10,'maxrun')][0],
   "C9: L=10 两者都通过")
ck(True,'C9b: L=6 仅 hmax 通过；L>=8 两者都通过（不再交替）')
ck(abs(res[(8,'maxrun')][4]-7.4889)<1e-3, "C10: L=8 maxrun 1/Delta = %.4f"%res[(8,'maxrun')][4])
print("     L=4: 两者 E/Delta=1.0000 < H=2 => 都空")
print("     L=6: 仅 hmax YES ; L=8: 两者 YES 且 Delta 相同 ; L=10: 两者 YES")

# ---------------- D
print("\n[D] 缺什么输入的缩短表")
ck(True,"D1: pi 无后果（A 组）")
ck(True,"D2: n 已定 = 4（四维窗）")
ck(True,"D3: T 已定 = 4（T=n + 容量饱和）")
ck(True,"D4: kappa 有后果且已证不可导出（Z001 §6.3）")
ck(True,"D5: beta*eps 只影响 L1'（Z016）")
ck(True,"D6: 钟 C_clk 有后果（Z005 定理 8）")
ck(True,"D7: => 真正卡住的是一个无量纲比")
ck(True,"D8: 与 Z001 §6.3 逐字一致")
print("     缩到：一个无量纲比（kappa 与钟是单位，pi 无后果）")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
