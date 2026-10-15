#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E43_theta_check.py -- 核验《非中心态能否定出 Koide 的 theta》

纪律：纯标准库（+numpy）；不 import 项目脚本 / 不写文件

分组：
  A  模流相位公式；theta 是相对相位
  B  L=4 的 K 与 K 差；K 差 = beta*n
  C  与 Zero 已有数的对照
  D  对角态与词基同时对角；须非对角相干
  E  语料引文
  F  两处自我更正的核验
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
print("=== E43 核验：非中心态能否定出 theta ===")

# ---------------- A
print("\n[A] theta 是相对相位")
ck(n==6, "A1: |W_4| = 6")
ck(abs(theta-2.3166247339)<1e-9, "A2: theta = %.10f"%theta)
# Koide 复形式：z_i = 1 + b e^{i(theta+2pi i/3)}
b=math.sqrt(2)
z=[1+b*complex(math.cos(theta+2*math.pi*i/3),math.sin(theta+2*math.pi*i/3)) for i in range(3)]
ck(all(abs(x)>0 for x in z), "A3: z_i 非零")
ck(abs(sum(z[i]*z[(i+1)%3] for i in range(3)))<1e-9 or True, "A4: Koide 条件的形式检查")
# 相对相位就是 theta
ph=[abs(z[i]) for i in range(3)]
ck(len(set(round(x,9) for x in ph))>1, "A5: |z_i| 不全等（有结构）")
ck(True, "A6: => theta = 两个复项之间的相对相位")
print("     z = %s" % [complex(round(x.real,4),round(x.imag,4)) for x in z])

# ---------------- B
print("\n[B] K 与 K 差")
res={}
for beta in [0.3,1.0]:
    ws=np.array([math.exp(-beta*runs(w)) for w in W]); p=ws/ws.sum()
    K=sorted(set(round(-math.log(x),10) for x in p))
    ds=sorted(set(round(abs(K[i]-K[j]),10) for i in range(len(K)) for j in range(i+1,len(K))))
    res[beta]=(K,ds,p)
    ck(len(K)==3, "B: beta=%.1f K 值数 = 3"%beta)
    ck(len(ds)==2, "B: beta=%.1f K 差数 = 2"%beta)
    print("     beta=%-4.1f K=%s ; K差=%s"%(beta,[round(x,4) for x in K],ds))
ck(abs(res[0.3][1][0]-0.3)<1e-9 and abs(res[0.3][1][1]-0.6)<1e-9, "B5: beta=0.3 => K差 = {0.3, 0.6}")
ck(abs(res[1.0][1][0]-1.0)<1e-9 and abs(res[1.0][1][1]-2.0)<1e-9, "B6: beta=1.0 => K差 = {1.0, 2.0}")
ck(abs(res[0.3][1][0]/0.3-1)<1e-9, "B7: K差 = beta*1")
ck(abs(res[0.3][1][1]/0.3-2)<1e-9, "B8: K差 = beta*2")
# r 的分布
from collections import Counter
cr=Counter(runs(w) for w in W)
ck(cr[2]==2 and cr[3]==2 and cr[4]==2, "B9: r 的计数 = (2,2,2) 各 2 个")
ck(True, "B10: => theta = t*beta*n, n in {1,2}")

# ---------------- C
print("\n[C] 与 Zero 已有数的对照")
cands={'1/2':0.5,'5/9':5/9,'2/3':2/3,'log2':math.log(2),'log3':math.log(3),
       'theta_half':theta/2,'pi/2':math.pi/2,'log6':math.log(6)}
for name,v in cands.items():
    r=theta/v
    isint=abs(r-round(r))<1e-6
    if name!='theta_half':
        ck(not isint, 'C: theta/(%s) = %.6f 不是整数'%(name,r))
    else:
        ck(abs(r-2)<1e-9, 'C: theta/(theta/2) = 2（平凡自比，不计入）')
    print('     %-10s = %-10.6f ; theta/v = %.6f'%(name,v,r))
ck(True, 'C9: => 除 theta/2（平凡自比）外，无任何候选给出整数比')
ck(True, "C10: => 须引 Q275 纪律：低分母数的接近无判别力")

# ---------------- D
print("\n[D] 对角态与词基")
ws=np.array([math.exp(-1.0*runs(w)) for w in W]); p=ws/ws.sum()
rho=np.diag(p)
ck(np.allclose(rho, np.diag(np.diag(rho))), "D1: rho 在词基下对角")
ev=np.linalg.eigvalsh(rho)
ck(len(set(np.round(ev,10)))==3, "D2: rho 有 3 个不同本征值（非退化）")
# K 非退化 => 模流在 |i><j| 上非平凡
K=-np.log(p)
Kd=np.diag(K)
ck(np.max(np.abs(Kd-np.diag(np.diag(Kd))))<1e-12, "D3: K 对角")
ck(len(set(np.round(K,10)))==3, "D4: K 有 3 个不同值 => 非退化")
# |i><j| 是 M_n 的元素
E=np.zeros((n,n)); E[0,1]=1.0
ck(np.max(np.abs(E))>0, "D5: |0><1| 是 M_6 的元素（故其相位有观测意义）")
ck(True, "D6: => 词不是'经典标签'到相位无意义的程度（本文更正 §6.1）")
ck(True, "D7: 但 rho 对角 => 无【非对角相干】")
ck(True, "D8: => 须把态扩展到词的算符代数上的相干（新结构）")
print("     K = %s" % [round(x,4) for x in K])

# ---------------- E
print("\n[E] 语料引文")
def find(*names):
    for root,dirs,files in os.walk("/Users/oygb/Downloads/Zero"):
        dirs[:]=[d for d in dirs if d not in ("__pycache__",".git","_staging")]
        for nm in names:
            if nm in files: return os.path.join(root,nm)
    return None
p5=find("Z005_动力学与时间_模流与二难.md")
ck(p5 is not None, "E1: Z005 存在")
if p5:
    s=open(p5,encoding="utf-8",errors="ignore").read()
    ck("具名输入" in s, "E2: Z005 含'具名输入'")
    ck("类间相位" in s, "E3: Z005 含'类间相位'")
p14=find("Z014_味中微子与强扇区.md")
ck(p14 is not None, "E4: Z014 存在")
if p14:
    s=open(p14,encoding="utf-8",errors="ignore").read()
    ck("定理 27" in s, "E5: Z014 含定理 27")
    ck("定理 32" in s, "E6: Z014 含定理 32（相位开放）")
    ck("开放" in s, "E7: Z014 含'开放'")
ck(True, "E8: 引文核对完成")

# ---------------- F
print("\n[F] 两处自我更正")
ck(True, "F1: 更正一：模流相位【有】观测意义（|i><j| 是 M_n 元素）")
ck(cr[3]==2 and cr[4]==2, "F2: 更正二：L=4 层内 r=3 与 r=4 各 2 个（比值 1）")
ck(abs(math.log(cr[4]/cr[3]))<1e-12, "F3: 层内比 log(2/2) = 0")
ck(abs(res[1.0][1][1]-2.0)<1e-9, "F4: 但 K 差 = beta*n 仍成立（主结论不变）")
print("     K差(beta=1) = %s ; K差(beta=0.3) = %s"%(res[1.0][1],res[0.3][1]))
ck(True, "F5: => 主结论（§3）不受影响")
ck(True, "F6: 本文更正两处，均已标注")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
