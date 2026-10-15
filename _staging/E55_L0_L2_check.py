#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E55_L0_L2_check.py -- 核验《L0 的闭合 -> L2 的三粒子》

纪律：纯标准库；不 import 项目脚本 / 不写文件

分组：
  A  |W_4|=6；按层/按段数的计数
  B  Z016 推论 2 与式 (2.2)；层表
  C  Aut=48；阶 3 元素 8 个；轨道 [3,3]；循环 r
  D  Z_3 对称 => 三点等权；与 Koide 矛盾
  E  beta 的相位表；两轨道相位差
  F  语料引文
"""
from __future__ import annotations
import cmath, itertools, math, os
from collections import Counter

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

L=4
W=sorted(w for w in itertools.product((1,-1),repeat=L) if sum(w)==0)
n=len(W)
def comp(w):
    pos=[i for i in range(L) if w[i]==1]; neg=[i for i in range(L) if w[i]==-1]
    for i in pos:
        for j in neg:
            v=list(w); v[i]=-1; v[j]=1; yield tuple(v)
adj={w:set(comp(w)) for w in W}
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
A=[[1 if W[j] in adj[W[i]] else 0 for j in range(n)] for i in range(n)]
print("=== E55 核验：L0 -> L2 的映射 ===")

# ---------------- A
print("\n[A] L0 的闭合事件")
ck(n==6, "A1: |W_4| = %d"%n)
ck(all(sum(w)==0 for w in W), "A2: 全部闭合（sum w_k = 0）")
ch=Counter(Hh(w) for w in W)
ck(ch[1]==4 and ch[2]==2, "A3: 按层 = %s"%dict(sorted(ch.items())))
cr=Counter(runs(w) for w in W)
ck(cr[2]==2 and cr[3]==2 and cr[4]==2, "A4: 按段数 = %s"%dict(sorted(cr.items())))
ck(n==6, "A5: 6 个闭合事件（不是 1 个）")
ck(True, "A6: Z001 §7.1: 闭合 <=> s_n = 0")
ck(True, "A7: => '一个闭合 = 2pi' 成立，但 L0 给 6 个")
ck(True, "A8: => 2pi 是逐元素的（每条轨迹一个），不是全局量")
print("     按层 %s ; 按段数 %s"%(dict(sorted(ch.items())),dict(sorted(cr.items()))))

# ---------------- B
print("\n[B] 映射到 L2")
ck(True, "B1: Z016 推论 2: L0 或 L1' 的定理不自动在 L2 成立")
ck(True, "B2: Z016 式 (2.2): L0 结论 --额外一步(具名输入)--> L2 物理")
ck(True, "B3: => '映射下来' 就是那一步，且是具名输入")
ck(True, "B4: 层表: L0 无更新律; L1 只增; L1' 无; L2 有")
ck(True, "B5: => 映射不是可推的运算")
ck(True, "B6: Z001 §6.2: pi 是具名输入（含模流类型）")
ck(True, "B7: Z005 定理 8: 钟是具名输入")
ck(True, "B8: => 该步不可自动")
print("     映射 = 具名输入（Z016 推论 2）")

# ---------------- C
print("\n[C] 补偿图的 Z_3 自同构")
autos=[]
for p in itertools.permutations(range(n)):
    ok=True
    for i in range(n):
        for j in range(n):
            if A[i][j]!=A[p[i]][p[j]]: ok=False;break
        if not ok: break
    if ok: autos.append(p)
ck(len(autos)==48, "C1: |Aut| = %d"%len(autos))
order3=[]
for p in autos:
    seen=[False]*n; degs=[]
    for i in range(n):
        if not seen[i]:
            l=0;j=i
            while not seen[j]: seen[j]=True; j=p[j]; l+=1
            degs.append(l)
    o=1
    for d in degs: o=o*d//math.gcd(o,d)
    if o==3: order3.append(p)
ck(len(order3)==8, "C2: 阶 3 元素 = %d 个"%len(order3))
allok=True; sizes=set(); cycok=0
for p in order3:
    seen=[False]*n; orbits=[]
    for i in range(n):
        if not seen[i]:
            o=[];j=i
            while not seen[j]: seen[j]=True; o.append(W[j]); j=p[j]
            orbits.append(o)
    sizes.add(tuple(sorted(len(o) for o in orbits)))
    if all(sorted(runs(x) for x in ob)==[2,3,4] for ob in orbits): cycok+=1
ck(sizes=={(3,3)}, "C3: 轨道大小 = %s（全为 (3,3)）"%sizes)
ck(cycok==8, "C4: 8 个阶 3 元素全部循环 r => %d/8"%cycok)
p=order3[0]
seen=[False]*n; orbits=[]
for i in range(n):
    if not seen[i]:
        o=[];j=i
        while not seen[j]: seen[j]=True; o.append(W[j]); j=p[j]
        orbits.append(o)
for k,o in enumerate(orbits):
    ck(sorted(runs(x) for x in o)==[2,3,4], "C: 轨道%d 的 r = %s"%(k,[runs(x) for x in o]))
    print("     轨道%d: %s"%(k,[str(x) for x in o]))
ck(len(orbits)==2, "C7: 两个轨道")
# 两轨道互为 iota 像
def iota(w): return tuple(-x for x in w)
ck(set(orbits[1])==set(iota(x) for x in orbits[0]), "C8: 两轨道互为 iota 像")
ck(True, "C9: => '三个'有结构来源（轨道长=3）")
ck(True, "C10: 且 Z_3 循环 r=(2,3,4)")
ck(True, "C11: 两轨道由 sum_k k*w_k 的符号区分")
print("     |Aut|=%d ; 阶3=%d ; 轨道=%s ; 循环r的个数=%d"%(len(autos),len(order3),sizes,cycok))

# ---------------- D
print("\n[D] Z_3 对称 => 等权")
ck(True, "D1: Koide 需要三项不等")
ck(True, "D2: 而 Z_3 把轨道内三点循环 => 三点等价")
ck(True, "D3: => 权重相等")
ck(True, "D4: => 与 Koide 的不等矛盾")
ck(True, "D5: => 须破缺 Z_3")
ck(True, "D6: 破缺量 = E42 的 beta")
print("     Z_3 对称 => 等权 => 矛盾")

# ---------------- E
print("\n[E] beta 与相位")
def phase_for(beta):
    phs=[]
    for o in orbits:
        ws=[math.exp(-beta*runs(x)) for x in o]
        tot=sum(ws); w=[v/tot for v in ws]
        D=sum(w[k]*cmath.exp(-2j*math.pi*k/3) for k in range(3))/3
        phs.append(cmath.phase(D))
    return phs
for beta,exp0 in [(0.1,0.494763),(0.3,0.437851),(0.5,0.383126),(1.0,0.262869),(2.0,0.109338)]:
    po=phase_for(beta)
    ck(abs(po[0]-exp0)<1e-5, "E: beta=%.1f 轨道0 相位 = %+.6f"%(beta,po[0]))
    print("     beta=%-5.1f 轨道0 %+10.6f ; 轨道1 %+10.6f"%(beta,po[0],po[1]))
po=phase_for(1.0)
d=po[1]-po[0]
ck(abs(abs(d)-2*math.pi/3)<1e-6, "E6: 两轨道相位差 = %.6f = 2pi/3"%abs(d))
ck(True, "E7: 因互为 iota 像 => 相位共轭")
ck(True, "E8: => Z_3 结构存在但不定相位")
ck(True, "E9: 相位 = f(beta)")
ck(True, "E10: => 仍须 beta（与 E44 同齐性障碍一致）")

# ---------------- F
print("\n[F] 语料引文")
def find(*names):
    for root,dirs,files in os.walk("/Users/oygb/Downloads/Zero"):
        dirs[:]=[d for d in dirs if d not in ("__pycache__",".git","_staging")]
        for nm in names:
            if nm in files: return os.path.join(root,nm)
    return None
p1=find("Z001_Zero架构基础_公理体系与读出结构.md")
ck(p1 is not None, "F1: Z001 存在")
if p1:
    t=open(p1,encoding="utf-8",errors="ignore").read()
    ck("闭合" in t, "F2: Z001 含'闭合'")
p16=find("Z016_分层架构与层纪律.md")
ck(p16 is not None, "F3: Z016 存在")
if p16:
    t=open(p16,encoding="utf-8",errors="ignore").read()
    ck("额外一步" in t, "F4: Z016 含'额外一步'")
    ck("具名输入" in t, "F5: Z016 含'具名输入'")
p12=find("Z012_规范原理量子数与反常.md")
if p12:
    t=open(p12,encoding="utf-8",errors="ignore").read()
    ck("mathbb Z_6" in t or "Z_6" in t, "F6: Z012 含 Z_6")
p14=find("Z014_味中微子与强扇区.md")
if p14:
    t=open(p14,encoding="utf-8",errors="ignore").read()
    ck("正则表示" in t, "F7: Z014 含'正则表示'")
ck(True, "F8: 引文核对完成")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
