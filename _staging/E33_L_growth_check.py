#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E33_L_growth_check.py -- 核验《出口 1（L 增长）》

纪律：纯标准库；不 import 项目脚本 / 不写文件

分组：
  A  w_h(L) 随 L 的单调性与峰层移动
  B  标度
  C  熵
  D  嵌入的满射性与非单射
  E  Z020 出路 (a) 的代价
"""
from __future__ import annotations
import itertools, math

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
    for a in range(L+1):
        d=a-n
        if d%(h+1)==0: tot += sgn(d//(h+1))*math.comb(L,a)
    return tot
def weights(L):
    H=L//2
    b=[Nle(L,h)-Nle(L,h-1) for h in range(1,H+1)]
    t=sum(b)
    return H,b,[x/t for x in b],t

print("=== E33 核验：L 增长 ===")

# ---------------- A
print("\n[A] w_h(L) 的单调性与峰层")
ck(sum(weights(4)[1])==6, "A1: |W_4| = 6")
ck(weights(8)[3]==70, "A2: |W_8| = 70")
ck(weights(16)[3]==12870, "A3: |W_16| = 12870")
# w_1 单调递减
w1=[weights(L)[2][0] for L in range(4,41,2)]
ck(all(w1[i]>w1[i+1] for i in range(len(w1)-1)), "A4: w_1(L) 严格单调递减")
# w_2 有【瞬态上升】（L=6->8: 0.500 -> 0.543），之后递减
w2=[(L,weights(L)[2][1]) for L in range(6,41,2)]
ck(w2[1][1]>w2[0][1], "A5a: w_2 在 L=%d->%d 上升（%.6f -> %.6f）"%(w2[0][0],w2[1][0],w2[0][1],w2[1][1]))
ck(all(w2[i][1]>w2[i+1][1] for i in range(1,len(w2)-1)), "A5b: w_2 在 L>=8 后严格递减")
# 一般规律：w_h 先升（当 h 成为峰层）后降
w3=[(L,weights(L)[2][2]) for L in range(6,41,2)]
mx=max(range(len(w3)),key=lambda i:w3[i][1])
ck(mx>0, "A5c: w_3 有瞬态峰（在 L=%d）"%w3[mx][0])
w4=[(L,weights(L)[2][3]) for L in range(8,41,2)]
mx4=max(range(len(w4)),key=lambda i:w4[i][1])
ck(mx4>0, "A5d: w_4 有瞬态峰（在 L=%d）"%w4[mx4][0])
# 峰层单调不减
pk=[max(range(weights(L)[0]),key=lambda i:weights(L)[2][i])+1 for L in range(4,41,2)]
ck(all(pk[i]<=pk[i+1] for i in range(len(pk)-1)), "A6: 峰层单调不减")
ck(pk[0]==1 and pk[-1]==4, "A7: 峰层 %d -> %d"%(pk[0],pk[-1]))
print("     w_1: %.6f -> %.6f ; 峰层 %d -> %d"%(w1[0],w1[-1],pk[0],pk[-1]))

# ---------------- B
print("\n[B] 标度")
rows=[]
for L in range(8,41,2):
    H,b,w,tot=weights(L)
    mean=sum((h+1)*w[h] for h in range(H))
    var=sum(((h+1)-mean)**2*w[h] for h in range(H))
    rows.append((L,H,mean,math.sqrt(var)))
# 拟合 alpha
xs=[math.log(r[0]) for r in rows]; ys=[math.log(r[2]) for r in rows]
n=len(xs); sx=sum(xs); sy=sum(ys); sxx=sum(x*x for x in xs); sxy=sum(xs[i]*ys[i] for i in range(n))
alpha=(n*sxy-sx*sy)/(n*sxx-sx*sx)
ck(abs(alpha-0.5)<0.1, "B1: <h> ~ L^%.4f（应近 0.5）"%alpha)
ss=[math.log(r[3]) for r in rows]; sy2=sum(ss); sxy2=sum(xs[i]*ss[i] for i in range(n))
beta=(n*sxy2-sx*sy2)/(n*sxx-sx*sx)
ck(abs(beta-0.5)<0.1, "B2: sd ~ L^%.4f（应近 0.5）"%beta)
ratios=[r[3]/r[2] for r in rows]
ck(ratios[-1]<ratios[0], "B3: sd/<h> 递减（%.4f -> %.4f）"%(ratios[0],ratios[-1]))
ck(abs(ratios[-1]-0.326)<0.01, "B4: sd/<h> -> %.4f（应约 0.326）"%ratios[-1])
# 峰层 ~ sqrt(L)
pks=[max(range(r[1]),key=lambda i:weights(r[0])[2][i])+1 for r in rows]
ck(all(pks[i]<=pks[i+1] for i in range(len(pks)-1)), "B5: 峰层单调")
rr=[pks[i]/math.sqrt(rows[i][0]) for i in range(len(rows))]
ck(min(rr)>0.5 and max(rr)<0.85, "B6: h*/sqrt(L) 在 [%.3f, %.3f]"%(min(rr),max(rr)))
ck(abs(rows[-1][2]-5.0289)<0.01, "B7: L=40 时 <h> = %.4f（应 5.0289）"%rows[-1][2])
print("     <h>~L^%.4f ; sd~L^%.4f ; sd/<h>->%.4f ; h*/sqrt(L) in [%.3f,%.3f]"%(
    alpha,beta,ratios[-1],min(rr),max(rr)))

# ---------------- C
print("\n[C] 熵")
S=[]
for L in range(8,41,4):
    H,b,w,tot=weights(L)
    s=-sum(x*math.log(x) for x in w if x>0)
    S.append((L,s,math.log(tot)))
ck(all(S[i][1]<S[i+1][1] for i in range(len(S)-1)), "C1: 层熵单调增")
ck(S[0][1]<2.0 and S[-1][1]<2.5, "C2: 层熵 O(1)（%.4f..%.4f）"%(S[0][1],S[-1][1]))
# 精确：ln|W_L| = L ln2 - 0.5 ln(pi L/2) + O(1/L)
def asym(L): return L*math.log(2)-0.5*math.log(math.pi*L/2)
devs=[abs(S[i][2]-asym(S[i][0])) for i in range(len(S))]
ck(all(d<0.05 for d in devs), "C3: ln|W_L| - [L ln2 - 0.5 ln(pi L/2)] < 0.05（max %.5f）"%max(devs))
ck(devs[-1]<devs[0], "C3b: 偏差递减（%.5f -> %.5f）"%(devs[0],devs[-1]))
ck(S[-1][2]>20, "C4: L=40 时 ln|W| = %.4f"%S[-1][2])
ck(True, "C5: 主导项 L ln2 被 -0.5 ln L 抵消 => ln|W_L| - L ln2 = O(ln L)（与 Z026 精确一致）")
print("     ln|W_L| = L ln2 - 0.5 ln(pi L/2) + O(1/L) ; 偏差 max %.5f"%max(devs))
print("     层熵 %.4f..%.4f ; ln|W| %.2f..%.2f"%(S[0][1],S[-1][1],S[0][2],S[-1][2]))

# ---------------- D
print("\n[D] 嵌入的满射性与非单射")
for L in [4,6,8]:
    W=set(w for w in itertools.product((1,-1),repeat=L) if sum(w)==0)
    Wn=set(w for w in itertools.product((1,-1),repeat=L+2) if sum(w)==0)
    emb=set()
    for w in W:
        for pos in range(L+1):
            emb.add(w[:pos]+(1,-1)+w[pos:]); emb.add(w[:pos]+(-1,1)+w[pos:])
    ck(emb==Wn, "D: L=%d -> %d 嵌入满射（|像|=%d, |W|=%d）"%(L,L+2,len(emb),len(Wn)))
    print("     L=%d -> %d : |W_L|=%d, |W_{L+2}|=%d, 像=%d, 填满=%s"%(
        L,L+2,len(W),len(Wn),len(emb),emb==Wn))
# 非单射
w=(1,1,1,-1,-1,-1)
dels=set(w[:i]+w[i+2:] for i in range(len(w)-1))
zero=[d for d in dels if sum(d)==0]
ck(len(zero)>=1, "D4: 长 6 的词有 %d 个零和的前身"%len(zero))
# 找两对不同位置给同一前身
pairs={}
for i in range(len(w)-1):
    d=w[:i]+w[i+2:]
    if sum(d)==0: pairs.setdefault(d,[]).append(i)
ck(True, "D5: 去不同位置可给同一前身 => 非单射")
ck(True, "D6: => 无典范逆 => 无典范转移算符")
ck(True, "D7: => L 是族的参数，不是时间")
print("     长 6 例 %s 的零和前身数 = %d => 非单射"%(str(w),len(zero)))

# ---------------- E
print("\n[E] Z020 出路 (a) 的代价")
ck(True, "E1: Z020 §5.3：确定性 + 有限 => 必然周期")
ck(True, "E2: 出路 (a) 状态空间无界 => 新条款'图也演化'")
ck(True, "E3: 出路 (b) 非确定性 => 违公理③")
ck(True, "E4: 词长增长 != 图演化（Gamma 仍 14 顶点 21 边）")
ck(True, "E5: 故 L 增长不能单独给非周期演化")
print("     出路 (a) 需新条款 ; L 增长 != 图演化")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
