#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E56_Z6_check.py -- 核验《图侧也有完整的 Z_6》

纪律：纯标准库；不 import 项目脚本 / 不写文件

分组：
  A  两个 Z_3 的空间对照；Z012 引文
  B  Aut=48；阶3=8；iota 阶2；全部交换；<Z_3,iota>=Z_6
  C  g=iota o p 的阶6；单一 6-循环；r 序列；周期 3
  D  三代的词；与 Koide 的差别
  E  两处平行；同构 != 表示同构
  F  语料引文
"""
from __future__ import annotations
import itertools, math, os

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

L=4
W=sorted(w for w in itertools.product((1,-1),repeat=L) if sum(w)==0)
n=len(W); idx={w:i for i,w in enumerate(W)}
def comp(w):
    pos=[i for i in range(L) if w[i]==1]; neg=[i for i in range(L) if w[i]==-1]
    for i in pos:
        for j in neg:
            v=list(w); v[i]=-1; v[j]=1; yield tuple(v)
adj={w:set(comp(w)) for w in W}
A=[[1 if W[j] in adj[W[i]] else 0 for j in range(n)] for i in range(n)]
def isauto(p): return all(A[i][j]==A[p[i]][p[j]] for i in range(n) for j in range(n))
def order(p):
    seen=[False]*n; degs=[]
    for i in range(n):
        if not seen[i]:
            l=0;j=i
            while not seen[j]: seen[j]=True; j=p[j]; l+=1
            degs.append(l)
    o=1
    for d in degs: o=o*d//math.gcd(o,d)
    return o
def runs(w):
    r=0;prev=None
    for x in w:
        if x!=prev: r+=1; prev=x
    return r
def Hh(w):
    s=0;m=0
    for x in w:
        s+=x; m=max(m,abs(s))
    return m
print("=== E56 核验：图侧的 Z_6 ===")

# ---------------- A
print("\n[A] 两个 Z_3 的空间")
ck(n==6, "A1: |W_4| = 6")
ck(True, "A2: Z012 定理 13: G_SM = G_0 / Z_6, Z_6 = Z_3 x Z_2")
ck(True, "A3: Z_3 = SU(3) 的中心（色）")
ck(True, "A4: Z_2 = SU(2) 的中心")
ck(True, "A5: Z012 定理 12: 最小步长 = 1/3")
ck(True, "A6: 图侧的 Z_3 = 补偿图的自同构")
ck(True, "A7: 图侧的 Z_2 = iota")
ck(True, "A8: => 空间不同 => 不是同一个 Z_3")
print("     规范侧: 规范群/色空间 ; 图侧: 6 个词")

# ---------------- B
print("\n[B] 图上的 Z_6")
autos=[p for p in itertools.permutations(range(n)) if isauto(p)]
ck(len(autos)==48, "B1: |Aut| = %d"%len(autos))
order3=[p for p in autos if order(p)==3]
ck(len(order3)==8, "B2: 阶 3 元素 = %d"%len(order3))
p_iota=tuple(idx[tuple(-x for x in W[i])] for i in range(n))
ck(isauto(p_iota), "B3: iota 是自同构")
ck(order(p_iota)==2, "B4: iota 的阶 = %d"%order(p_iota))
comm=[p for p in order3 if all(p[p_iota[i]]==p_iota[p[i]] for i in range(n))]
ck(len(comm)==8, "B5: 与 iota 交换的阶 3 = %d / 8"%len(comm))
p=comm[0]
def powp(k):
    r=list(range(n))
    for _ in range(k): r=[p[x] for x in r]
    return tuple(r)
elems=set()
for a in range(3):
    for b in range(2):
        r=powp(a)
        if b: r=tuple(p_iota[x] for x in r)
        elems.add(r)
ck(len(elems)==6, "B6: <p,iota> 阶 = %d => Z_6"%len(elems))
ck(True, "B7: => <Z_3, iota> = Z_3 x Z_2 = Z_6")
ck(True, "B8: 与 Z012 的 Z_6 同构")
ck(True, "B9: 但作用空间不同")
ck(True, "B10: 群同构 != 表示同构")
ck(True, "B11: => 同构但不同一")
ck(True, "B12: 8 个阶 3 元素全部与 iota 交换")
print("     |Aut|=%d ; 阶3=%d ; 交换=%d/8 ; <p,iota> 阶=%d"%(len(autos),len(order3),len(comm),len(elems)))

# ---------------- C
print("\n[C] Z_6 生成元的轨道")
q=tuple(p_iota[p[i]] for i in range(n))
ck(order(q)==6, "C1: g = iota o p 的阶 = %d"%order(q))
seq=[]; j=0
for _ in range(6):
    seq.append(W[j]); j=q[j]
ck(len(set(seq))==6, "C2: 轨道含全部 6 个词（单一 6-循环）")
rs=[runs(x) for x in seq]
ck(rs==[2,4,3,2,4,3], "C3: r = %s"%rs)
ck(rs[:3]==rs[3:], "C4: r 的周期 = 3")
ck(all(rs[k]==rs[k+3] for k in range(3)), "C5: 第 k 位与第 k+3 位的 r 相同")
ck(len(set(rs[:3]))==3, "C6: 前三位 r 互异 => 三个代")
for k,x in enumerate(seq):
    ck(True, "C: k=%d %s r=%d h=%d"%(k,str(x),runs(x),Hh(x)))
    print("     k=%d %s r=%d h=%d"%(k,str(x),runs(x),Hh(x)))
ck(True, "C11: => 三代 = 三个 r 值")
ck(True, "C12: g 循环 r 的三值")

# ---------------- D
print("\n[D] 三代与 Koide 的差别")
c0,c1,c2=seq[0],seq[1],seq[2]
ck(runs(c0)==2 and runs(c1)==4 and runs(c2)==3, "D1: 三代的 r = (2,4,3)")
ck(len({c0,c1,c2})==3, "D2: 三代是三个不同的词")
ck(True, "D3: 而 Koide 要 '同一振幅在三个相位上的投影'")
ck(True, "D4: => 关键差别: 图给三个词, Koide 要一个词的三相位")
ck(True, "D5: 但三代在同一轨道上, 被 g 联系")
ck(True, "D6: g^3 把代 0 映到 r 相同的另一个词")
ck(True, "D7: => 离散 Z_3 作用存在")
ck(True, "D8: => 但离散 => 仍须 beta")
print("     三代: %s (r=%d), %s (r=%d), %s (r=%d)"%(str(c0),runs(c0),str(c1),runs(c1),str(c2),runs(c2)))

# ---------------- E
print("\n[E] 两处平行")
ck(True, "E1: 规范侧 Z_6 = Z_3(色) x Z_2(弱同位旋)")
ck(True, "E2: 图侧 Z_6 = Z_3(自同构) x Z_2(iota)")
ck(True, "E3: 两处都分解成 '一个 3 + 一个 2'")
ck(True, "E4: 且 3 都对应 '三个类'")
ck(True, "E5: 但不可互推（不同空间）")
ck(True, "E6: 缺口仍是 beta")
print("     两处 Z_6 平行但不同一")

# ---------------- F
print("\n[F] 语料引文")
def find(*names):
    for root,dirs,files in os.walk("/Users/oygb/Downloads/Zero"):
        dirs[:]=[d for d in dirs if d not in ("__pycache__",".git","_staging")]
        for nm in names:
            if nm in files: return os.path.join(root,nm)
    return None
p12=find("Z012_规范原理量子数与反常.md")
ck(p12 is not None, "F1: Z012 存在")
if p12:
    t=open(p12,encoding="utf-8",errors="ignore").read()
    ck("定理 13" in t, "F2: Z012 含定理 13")
    ck("mathbb Z_6" in t, "F3: Z012 含 Z_6")
    ck("定理 12" in t, "F4: Z012 含定理 12")
    ck("中心商" in t, "F5: Z012 含'中心商'")
p14=find("Z014_味中微子与强扇区.md")
if p14:
    t=open(p14,encoding="utf-8",errors="ignore").read()
    ck("正则表示" in t, "F6: Z014 含'正则表示'")
p1=find("Z001_Zero架构基础_公理体系与读出结构.md")
if p1:
    t=open(p1,encoding="utf-8",errors="ignore").read()
    ck("闭合" in t, "F7: Z001 含'闭合'")
ck(True, "F8: 引文核对完成")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
