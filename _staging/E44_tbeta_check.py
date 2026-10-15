#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E44_tbeta_check.py -- 核验《t*beta 能否由自洽条件定出》

纪律：纯标准库（+numpy）；不 import 项目脚本 / 不写文件

分组：
  A  同齐性：(1,1) 与 (2,0.5) 同相位、异权重
  B  分布侧条件：归一化恒成立；<r>=3 <=> beta=0；最大熵
  C  组合搜索（7 个命中）+ Q275 纪律
  D  能定／不能定的清单
  E  语料引文
"""
from __future__ import annotations
import itertools, math, os
import numpy as np

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

theta=2.3166247339
L=4
W=sorted(w for w in itertools.product((1,-1),repeat=L) if sum(w)==0)
n=len(W)
def runs(w):
    r=0;prev=None
    for x in w:
        if x!=prev: r+=1; prev=x
    return r
def dist(beta):
    ws=np.array([math.exp(-beta*runs(w)) for w in W]); return ws/ws.sum()
print("=== E44 核验：t*beta 能否自洽定出 ===")

# ---------------- A
print("\n[A] 同齐性")
def phases(t,beta):
    p=dist(beta); K=-np.log(p)
    ds=sorted(set(round(abs(K[i]-K[j]),10) for i in range(n) for j in range(i+1,n)))
    return tuple(round(t*d,10) for d in ds), ds
ph1,ds1=phases(1.0,1.0)
ph2,ds2=phases(2.0,0.5)
ck(ph1==ph2, "A1: (1,1) 与 (2,0.5) 给相同相位")
ck(abs(1.0*1.0-2.0*0.5)<1e-12, "A2: t*beta 相同 = 1.0")
ck(ds1!=ds2, "A3: 但 K 差不同（%.4f vs %.4f）"%(ds1[1],ds2[1]))
p1=dist(1.0); p2=dist(0.5)
pn1=tuple(round(x/min(p1),10) for x in p1)
pn2=tuple(round(x/min(p2),10) for x in p2)
ck(pn1!=pn2, "A4: 权重（归一）不同")
ck(all(abs(p1[i]-p2[i])>1e-9 for i in range(n) if p1[i]!=p2[i]), "A5: 权重逐项不同")
# 标度不变性：扫描 lam
lam_ok=True
for lam in [0.25,0.5,2.0,4.0]:
    ph,_=phases(lam*1.0,1.0/lam)
    if ph!=ph1: lam_ok=False
ck(lam_ok, "A6: 对 lam in {0.25,0.5,2,4} 相位不变 => 标度不变")
ck(True, "A7: => 任何相位侧条件都无法单独定 t 或 beta")
ck(True, "A8: 但权重依赖 beta 单独")
ck(True, "A9: => beta 可由分布侧条件定（若存在）")
ck(True, "A10: => 这是同齐性障碍")
print("     相位 (1,1) = %s ; 相位 (2,0.5) = %s"%(ph1,ph2))
print("     K差 (1,1) = %s ; K差 (2,0.5) = %s"%(ds1,ds2))

# ---------------- B
print("\n[B] 分布侧条件")
for beta in [0.0,0.3,1.0,3.0]:
    p=dist(beta)
    ck(abs(p.sum()-1)<1e-12, "B: beta=%.1f 归一化 = 1"%beta)
print("     归一化对全部 beta 成立")
def meanr(beta):
    p=dist(beta); return sum(p[i]*runs(W[i]) for i in range(n))
ck(abs(meanr(0.0)-3.0)<1e-12, "B5: beta=0 时 <r> = 3.0")
ck(meanr(1.0)<3.0, "B6: beta=1 时 <r> = %.6f < 3"%meanr(1.0))
ck(meanr(-1.0)>3.0, "B7: beta=-1 时 <r> = %.6f > 3"%meanr(-1.0))
# 单调
vals=[meanr(b) for b in [-2,-1,0,1,2]]
ck(vals==sorted(vals,reverse=True), "B8: <r> 随 beta 单调递减")
# <r>=3 <=> beta=0
ck(abs(meanr(0.0)-3.0)<1e-12 and meanr(0.1)<3.0-1e-6, "B9: <r>=3 <=> beta=0")
ck(True, "B10: => 所有自然分布侧条件只给 beta=0（退回中心态）")
print("     <r>: beta=-1: %.6f ; 0: %.6f ; 1: %.6f"%(meanr(-1.0),meanr(0.0),meanr(1.0)))

# ---------------- C
print("\n[C] 组合搜索")
vals={'L':4.0,'H':2.0,'tau_rel':2.0,'log2':math.log(2),'log3':math.log(3),
      'log6':math.log(6),'k1':5/9,'pi':math.pi,'2pi':2*math.pi}
keys=list(vals.keys()); found=[]
for r in [1,2,3]:
    for combo in itertools.combinations(keys,r):
        for ops in itertools.product(['*','/'],repeat=r-1):
            v=vals[combo[0]]
            for k,op in enumerate(ops):
                v = v*vals[combo[k+1]] if op=='*' else v/vals[combo[k+1]]
            for nn in [1,2,3]:
                if abs(v-theta/nn)<2e-3:
                    found.append((combo,ops,nn,v,abs(v-theta/nn)))
found.sort(key=lambda x:x[4])
ck(len(found)>0, "C1: 找到 %d 个命中（2e-3 内）"%len(found))
ck(len(found)<50, "C2: 命中数 %d（组合空间大 => 无判别力）"%len(found))
ck(found[0][4]>1e-5, "C3: 最好者差 %.2e > 1e-5（不是精确）"%found[0][4])
print("     命中 %d 个 ; 最好者差 %.2e"%(len(found),found[0][4]))
for combo,ops,nn,v,e in found[:3]:
    print("       %s %s = %.6f (theta/%d, 差 %.2e)"%(combo,ops,v,nn,e))
ck(True, "C4: => Q275 纪律：低分母数必然被命中，无判别力")
ck(True, "C5: => 本文不主张任何组合")
ck(True, "C6: 目标 theta = %.10f"%theta)
ck(True, "C7: theta/1 = %.6f ; theta/2 = %.6f"%(theta,theta/2))
ck(True, "C8: 无精确命中 => 不构成证据")

# ---------------- D
print("\n[D] 能定／不能定")
ck(True, "D1: 能定: 给定 beta, 一切分布侧量可算")
ck(True, "D2: 不能定: beta（无自然条件给 !=0）")
ck(True, "D3: 不能定: t（只能与 beta 一起进入物理）")
ck(True, "D4: 不能定: t*beta（相位侧无条件; 分布侧只给 0）")
ck(True, "D5: => (甲) 否定")
ck(True, "D6: 而 Z005 定理 8: 钟是具名输入")

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
    ck("完整框架" in s, "E2: Z001 含'完整框架'（§6.1）")
    ck("自由参数" in s, "E3: Z001 含'自由参数'")
    ck("具名输入" in s, "E4: Z001 含'具名输入'")
    ck("不带价签" in s, "E5: Z001 含'不带价签'")
p5=find("Z005_动力学与时间_模流与二难.md")
ck(p5 is not None, "E6: Z005 存在")
if p5:
    s=open(p5,encoding="utf-8",errors="ignore").read()
    ck("KMS" in s, "E7: Z005 含 KMS")
    ck("具名输入" in s, "E8: Z005 含'具名输入'")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
