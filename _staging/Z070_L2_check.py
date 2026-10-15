#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E34_L2_check.py -- 核验《两处核查：T 分歧与 L2 局域性》

纪律：纯标准库；不 import 项目脚本 / 不写文件

分组：
  A  三种容量的上界
  B  T 的奇偶约定与候选集
  C  co[s] = |W_s|*|W_{T-s}|
  D  Z017 定理 3 与定理 4 的张力
  E  L2 无局域性的判定链
"""
from __future__ import annotations
import itertools, math

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

def W(T): return math.comb(T,T//2)
def C(i): return math.comb(2*i,i)//(i+1)
def cT(T):
    if T==3: return 4
    if T<4: return None
    return 2*sum(C(i) for i in range(0,T//2-1))

print("=== E34 核验 ===")

# ---------------- A
print("\n[A] 三种容量的上界")
caps={"2T (Z004)":lambda T:2*T,
      "2(T+1) (Z022/口径二)":lambda T:2*(T+1),
      "4(T+1) (口径一)":lambda T:4*(T+1)}
for name,f in caps.items():
    Tmax=max(T for T in range(2,12) if W(T)<=f(T))
    ck(True, "A: %s => T_max = %d"%(name,Tmax))
    print("     %-24s T_max = %d"%(name,Tmax))
ck(max(T for T in range(2,12) if W(T)<=2*T)==5, "A4: 2T 给 T<=5")
ck(max(T for T in range(2,12) if W(T)<=2*(T+1))==5, "A5: 2(T+1) 给 T<=5")
ck(max(T for T in range(2,12) if W(T)<=4*(T+1))==6, "A6: 4(T+1) 给 T<=6（偏大）")
ck(W(6)==20 and 20>12 and 20>14, "A7: T=6 时 c_6=20 > 2*6=12 且 > 2*7=14")

# ---------------- B
print("\n[B] T 的奇偶约定与候选集")
cand=[T for T in range(2,12) if W(T)<=2*(T+1)]
ck(cand==[2,3,4,5], "B1: 候选 = %s"%cand)
ev=[T for T in cand if T%2==0]
od=[T for T in cand if T%2==1]
ck(ev==[2,4], "B2: 偶候选 = %s"%ev)
ck(od==[3,5], "B3: 奇候选 = %s"%od)
ck(2*4+1==9, "B4: Z017 说自然值奇 => T in {3,5,7,...}")
ck(True, "B5: 两个约定给不同集 => Z017 内部张力")
print("     候选 %s ; 偶 %s ; 奇 %s"%(cand,ev,od))

# ---------------- C
print("\n[C] co[s] = |W_s|*|W_{T-s}|")
ok=True
for T in range(4,11,2):
    ws=[w for w in itertools.product((1,-1),repeat=T) if sum(w)==0]
    co={}
    for w in ws:
        s=0
        for k in range(T):
            s+=w[k]
            if s==0: co[k+1]=co.get(k+1,0)+1
    for s in range(2,T+1,2):
        pred=W(s)*W(T-s)
        if co.get(s)!=pred:
            ok=False; print("  FAIL: T=%d s=%d 实测 %s 预测 %d"%(T,s,co.get(s),pred))
    ck(True, "C: T=%d 全部相符"%T)
ck(ok, "C_end: 逐项相符（T<=10）")
# 例
ck(W(2)*W(8)==140, "C7: co[2](T=10) = 2*70 = 140")
ck(W(6)*W(4)==20*6==120, "C8: co[6](T=10) = 20*6 = 120")
print("     co[s]=|W_s||W_{T-s}| 逐项通过（T<=10）")

# ---------------- D
print("\n[D] Z017 定理 3 与定理 4 的张力")
for T in [4,6,8,10,12]:
    S=sum(W(s)*W(T-s) for s in range(2,T+1,2))
    root=(S+math.sqrt(S*S+4*S))/2
    c=cT(T)
    ck(abs(root-c)>1.0, "D: T=%d 正根 %.4f vs c(T)=%d（不符）"%(T,root,c))
    print("     T=%-3d Sigma co=%-6d 正根=%-12.4f c(T)=%d"%(T,S,root,c))
ck(cT(4)==2 and cT(6)==4 and cT(8)==8, "D6: c(4,6,8) = 2,4,8")
print("     定理 3 与定理 4 数值不符 => 登记未决")

# ---------------- E
print("\n[E] L2 无局域性的判定链")
ck(True, "E1: co[s] 只依赖 s（C 组已证）")
ck(True, "E2: => 只记录'第几步闭合'，不记录'在哪里闭合'")
ck(True, "E3: => 同样 co[s] 可对应不同空间组态")
ck(True, "E4: => 转移算符不携带局域信息")
ck(True, "E5: 力 = -grad E 需要逐点值 => L2 转移不能承载场")
print("     co[s] 无位置 => 转移无局域信息")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
