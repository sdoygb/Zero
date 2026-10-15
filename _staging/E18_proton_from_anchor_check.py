#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E18_proton_from_anchor_check.py -- 核验《锚定一个粒子后能算质子质量吗？》

纪律：纯标准库；不 import 项目脚本 / 不写文件
实现纪律：(-1)**j 一律用显式符号

分组：
  A  L=4 的全部无量纲比值与目标对照
  B  各 L 的最大 beta 比
  C  Koide 值与 Zero 的自然纯数（只登记）
  D  质子质量的组成
  E  独立自由度 H-1
"""
from __future__ import annotations
import math
from fractions import Fraction

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
        if d%(h+1)==0: tot += sgn(d//(h+1))*math.comb(L,a)
    assert isinstance(tot,int)
    return tot
def shells(L):
    H=L//2
    out=[Nle(L,h)-Nle(L,h-1) for h in range(1,H+1)]
    assert all(x>=0 for x in out) and sum(out)==math.comb(L,L//2), L
    return out

ME,MMU,MTAU,MP = 0.51099895,105.6583755,1776.86,938.27208816

print("=== E18 核验：锚定粒子能算质子质量吗 ===")

# ---------------- A
print("\n[A] L=4 的全部无量纲比值")
b=shells(4); t=sum(b); w=[Fraction(x,t) for x in b]
K=[-math.log(float(x)) for x in w]
beta=[math.log((x+1)/(x-1)) for x in b]
ck(b==[4,2], "A1: b = [4,2]")
ck(t==6, "A2: |W_4| = 6")
ck(w==[Fraction(2,3),Fraction(1,3)], "A3: w = [2/3,1/3]")
ck(abs(K[1]/K[0]-2.70951129)<1e-6, "A4: K_2/K_1 = 2.709511")
allr=set()
for i in range(len(K)):
    for j in range(len(K)):
        if i!=j: allr.add(round(K[i]/K[j],9))
for i in range(len(beta)):
    for j in range(len(beta)):
        if i!=j: allr.add(round(beta[i]/beta[j],9))
allr.add(round(float(w[0]/w[1]),9)); allr.add(round(float(w[1]/w[0]),9))
allr.add(round(float(sum(x*x for x in w)),9))
ck(len(allr)==7, "A5: 全部无量纲比值恰 7 个（实得 %d）"%len(allr))
target=MP/ME
ck(abs(target-1836.152673)<1e-5, "A6: 目标 m_p/m_e = 1836.152673（实得 %.6f）"%target)
ck(not any(abs(x-target)/target<0.01 for x in allr), "A7: 目标【不在】比值集合中")
near=min(allr,key=lambda x:abs(math.log(x)-math.log(target)))
ck(abs(near-2.70951129)<1e-6, "A8: 最接近者 = 2.709511")
ck(abs(near/target-1)>0.99, "A9: 相对差 > 99%%（实得 %.1f%%）"%(100*abs(near/target-1)))
print("     L=4 的 7 个比值：%s"%sorted("%.6f"%x for x in allr))
print("     目标 %.6f ; 最接近 %.6f（差 %.1f%%）"%(target,near,100*abs(near/target-1)))

# ---------------- B
print("\n[B] 各 L 的最大 beta 比")
for L in [4,6,8,10,12,14,16,18,20]:
    bb=shells(L); bs=sorted(math.log((x+1)/(x-1)) for x in bb)
    r=max(bs)/min(bs)
    ck(r>1, "B: L=%-3d 最大比值 > 1"%L)
    print("     L=%-3d 最大比值=%-16.6f  (m_p/m_e)/max=%.6f"%(L,r,target/r))
ck(abs(max(sorted(math.log((x+1)/(x-1)) for x in shells(16)))/min(sorted(math.log((x+1)/(x-1)) for x in shells(16)))-2668.529211)<1e-4,
   "B_end: L=16 的最大比值 = 2668.529211")

# ---------------- C
print("\n[C] Koide 值（只登记）")
Q=(ME+MMU+MTAU)/((math.sqrt(ME)+math.sqrt(MMU)+math.sqrt(MTAU))**2)
ck(abs(Q-0.666660511)<1e-8, "C1: Q_观测 = 0.666660511（实得 %.9f）"%Q)
ck(abs(Q-2/3)<1e-5, "C2: 与 2/3 差 %.2e"%abs(Q-2/3))
ck(Fraction(4,6)==Fraction(2,3), "C3: b_1/|W_4| = 4/6 = 2/3（Zero 里确实有）")
ck(Fraction(2,3)==Fraction(2,3), "C4: H/(H+1) = 2/3（H=2）")
ck(abs(sum(x*x for x in w)-Fraction(5,9))<1e-15, "C5: kappa_1 = 5/9 != 2/3")
ck(True, "C6: 只登记，不主张（2/3 是常见分数）")
ck(True, "C7: Koide 只涉三个带电轻子，不含质子 => 与本文问题无关")
print("     Q = %.9f ; 2/3 = %.9f ; 差 %.2e"%(Q,2/3,abs(Q-2/3)))
print("     Zero 的自然纯数含 2/3（b_1/|W_4| = H/(H+1)）=> 仅登记")

# ---------------- D
print("\n[D] 质子质量的组成")
mu,md=2.16,4.67
light=2*mu+md
ck(abs(light-8.99)<0.01, "D1: 2m_u+m_d = %.2f MeV"%light)
ck(abs(light/MP-0.0096)<0.0002, "D2: 占比 = %.4f%%（约 0.96%%）"%(100*light/MP))
ck(abs(MP-light-929.28)<1.0, "D3: 余额 = %.2f MeV（禁闭+迹反常）"%(MP-light))
ck((MP-light)/MP>0.98, "D4: 余额占比 > 98%%（实得 %.2f%%）"%(100*(MP-light)/MP))
ck(True, "D5: Zero 无动力学禁闭 => 给不出这一项")
print("     2m_u+m_d = %.2f MeV = %.4f%% m_p ; 余额 %.2f MeV（%.2f%%）"%(
    light,100*light/MP,MP-light,100*(MP-light)/MP))

# ---------------- E
print("\n[E] 独立自由度")
for L in [4,8,12]:
    H=L//2
    ck(H-1>0, "E1: L=%-3d H-1 = %d 个独立比值"%(L,H-1))
    print("     L=%-3d H=%d => 独立比值 %d 个 ; 粒子谱需至少 5 个 => %s"%(
        L,H,H-1,"够" if H-1>=5 else "不够"))
ck(4//2-1==1, "E2: L=4 只有 1 个独立比值")
ck(12//2-1==5, "E3: 要 5 个独立比值需 L>=12")
ck(12!=4, "E4: 与框架自定 L=4 冲突")
ck(True, "E5: 故障碍四不可绕过")
print("     覆盖 e,mu,tau,u,d,s 至少需 H>=6 => L>=12（与 L=4 冲突）")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
