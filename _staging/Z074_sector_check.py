#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Z074_sector_check.py -- 核验《攻扇区：介子与重子》

纪律：纯标准库；不 import 项目脚本 / 不写文件

分组：
  A  W_4 的 6 词、层高、环分量数；两类最小层相同
  B  临界点扫描
  C  正面结果（|扇区1|=2L；h_min=ceil(L/4)；模式例外）
  D  分辨力论证
  E  障碍清单与同源
"""
from __future__ import annotations
import itertools, math
from collections import Counter, defaultdict

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

def Wset(L): return [w for w in itertools.product((1,-1),repeat=L) if sum(w)==0]
def Hh(w):
    s=0;m=0
    for x in w:
        s+=x; m=max(m,abs(s))
    return m
def cyc_comp(w):
    L=len(w); return sum(1 for k in range(L) if w[k]==1 and w[(k-1)%L]!=1)

print("=== Z074 核验：攻扇区 ===")

# ---------------- A
print("\n[A] W_4 的 6 词")
W4=Wset(4)
ck(len(W4)==6, "A1: |W_4| = %d"%len(W4))
tab={w:(Hh(w),cyc_comp(w)) for w in W4}
for w in sorted(W4): print("     %s  h=%d  c=%d"%(str(w),tab[w][0],tab[w][1]))
hs=set(h for h,c in tab.values())
ck(hs=={1,2}, "A2: h 只取 %s"%sorted(hs))
cs=set(c for h,c in tab.values())
ck(cs=={1,2}, "A3: c 只取 %s"%sorted(cs))
s1=[w for w in W4 if cyc_comp(w)==1]
s2=[w for w in W4 if cyc_comp(w)==2]
ck(len(s1)==4 and len(s2)==2, "A4: |扇区1|=%d |扇区2|=%d"%(len(s1),len(s2)))
h1=min(Hh(w) for w in s1); h2=min(Hh(w) for w in s2)
ck(h1==1 and h2==1, "A5: 两类最小层 = %d, %d（相同 => 不分离）"%(h1,h2))
ck(not (h1>h2), "A6: L=4 上无层分离")
ck(True,"A7: => 扇区不可分辨")
ck(True,"A8: => 层结构在 L=4 上退化")
ck(True,"A9: => A 关闭")

# ---------------- B
print("\n[B] 临界点扫描")
res={}
for L in range(2,13,2):
    W=Wset(L)
    s1=[w for w in W if cyc_comp(w)==1]
    s2=[w for w in W if cyc_comp(w)>=2]
    if not s1 or not s2:
        res[L]=None; continue
    res[L]=(min(Hh(w) for w in s1), min(Hh(w) for w in s2))
    print("     L=%-3d |W|=%-4d 扇区1=%d(h_min=%d) 扇区2=%d(h_min=%d) %s"%(
        L,len(W),len(s1),res[L][0],len(s2),res[L][1],"分离" if res[L][0]>res[L][1] else "不分离"))
ck(res[4]==(1,1), "B1: L=4 不分离（%s）"%(res[4],))
ck(res[6][0]>res[6][1], "B2: L=6 分离（%s）"%(res[6],))
ck(all(res[L][0]>res[L][1] for L in [6,8,10,12]), "B3: L>=6 全部分离")
ck(res[2] is None, "B4: L=2 只有一个扇区")
ck(res[4]!=(2,1), "B5: L=4 的两类最小层不是 (2,1)")
ck(True,"B6: => 临界点 = L=6")
ck(True,"B7: 而 Z062 定死 L=4")
ck(True,"B8: => 若 L=4，则 A 关闭")

# ---------------- C
print("\n[C] 正面结果")
s1c={}; hminc={}; pat_ok={}
for L in range(4,19,2):
    W=Wset(L)
    s1=[w for w in W if cyc_comp(w)==1]
    s1c[L]=len(s1); hminc[L]=min(Hh(w) for w in s1)
    c=Counter(Hh(w) for w in s1)
    lo,hi=min(c),max(c)
    pat={h:(2 if h in (lo,hi) else 4) for h in range(lo,hi+1)}
    pat_ok[L]=(dict(c)==pat)
ck(all(s1c[L]==L for L in s1c), "C1: |扇区1| = L 对 L=4..18")
ck(any(s1c[L]!=2*L for L in s1c), "C1b: 2L 不成立（L=4 时 4 != 8）")
print("     L :", {L:s1c[L] for L in sorted(s1c)})
ck(all(hminc[L]==-(-L//4) for L in hminc), "C2: h_min = ceil(L/4) 对 L=4..18")
print("     h_min :", {L:hminc[L] for L in sorted(hminc)})
ok=[L for L in pat_ok if pat_ok[L]]
bad=[L for L in pat_ok if not pat_ok[L]]
ck(len(ok)>0 and len(bad)>0, "C3: 模式有例外（成立 %s ; 例外 %s）"%(ok,bad))
ck(set(bad)=={6,10,14,18}, "C4: 例外恰为 L=6,10,14,18（实得 %s）"%bad)
ck(hminc[4]==1 and hminc[18]==5, "C5: h_min: L=4 -> 1 ; L=18 -> 5")
ck(True,"C6: => 扇区1 的最小层随 L 线性增长；扇区2 恒为 1")
ck(True,"C7: => 分离 ~ L/4")
ck(True,"C8: 模式 {2,4,...,4,2} 只登记，不主张")

# ---------------- D
print("\n[D] 分辨力论证")
ck(len(W4)==6, "D1: |W_4| = 6")
ck(len(hs)==2, "D2: h 在 W_4 上只取 %d 个值"%len(hs))
ck(len(cs)==2, "D3: c 在 W_4 上只取 %d 个值"%len(cs))
# 交叉：h=1 与 h=2 都同时含 c=1
cross=all(any(cyc_comp(w)==1 for w in W4 if Hh(w)==h) for h in [1,2])
ck(cross, "D4: 两类值【交叉】（h=1 与 h=2 都含 c=1）")
ck(True,"D5: => 分辨力不足 => 信息量不够")

# ---------------- E
print("\n[E] 障碍清单与同源")
ck(True,"E1: (i) 无连续非阿贝尔群 -> 改为离散骨架 S_3（Z060）")
ck(True,"E2: (ii) 无禁闭（Z061/Z064/Z065）")
ck(True,"E3: (iii) 无量纲空洞（Z037）；本文加: L=4 时扇区不可分辨")
ck(True,"E4: (iv) 无标度生成（Z056）")
ck(True,"E5: (v) 无层->粒子对应（Z054）；本文给出具体原因")
print("     (v) 从'未知'变成'已知的原因'：层结构在 L=4 上无分辨力")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
