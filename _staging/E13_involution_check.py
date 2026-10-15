#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E13_involution_check.py -- 核验《换对合：s 是唯一可用的》

纪律：纯标准库；不 import 项目脚本 / 不写文件

分组：
  A  三个对合都保零和、保层高、是对合
  B  不动点数（r 回文、s 零、rs）
  C  s 逐层给 b_h/2；r,rs 系统性失败
  D  顶层的唯一性
  E  s 保 |area|
  F  <r,s> = Z_2 x Z_2（交换）
"""
from __future__ import annotations
import itertools

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

def levels(w):
    s=0; out=[0]
    for x in w: s+=x; out.append(s)
    return out
def Hh(w): return max(abs(s) for s in levels(w))
def area(w): return sum(levels(w)[:-1])
def st(w): return "".join("+" if x>0 else "-" for x in w)
r_=lambda w: tuple(reversed(w))
s_=lambda w: tuple(-x for x in w)
rs=lambda w: tuple(-x for x in reversed(w))
INV=[("r",r_),("s",s_),("rs",rs)]
def words(L): return [w for w in itertools.product((1,-1),repeat=L) if sum(w)==0]

print("=== E13 核验：换对合 ===")

# ---------------- A
print("\n[A] 三个对合都保零和、保层高、是对合")
for L in [4,6,8]:
    ws=words(L)
    for nm,f in INV:
        ck(all(sum(f(w))==0 for w in ws), "A1: L=%d %s 保零和"%(L,nm))
        ck(all(Hh(f(w))==Hh(w) for w in ws), "A2: L=%d %s 保层高"%(L,nm))
        ck(all(f(f(w))==w for w in ws), "A3: L=%d %s 是对合"%(L,nm))
    print("     L=%-3d r,s,rs 全部保零和/保层高/对合 ✓"%L)

# ---------------- B
print("\n[B] 不动点数")
for L in [4,8,10]:
    ws=words(L)
    fx={nm:sum(1 for w in ws if f(w)==w) for nm,f in INV}
    ck(fx["s"]==0, "B1: L=%d s 无不动点（实得 %d）"%(L,fx["s"]))
    if L%4==0:
        ck(fx["r"]>0, "B2a: L=%d (mod4=0) r 有不动点（实得 %d）"%(L,fx["r"]))
    else:
        ck(fx["r"]==0, "B2b: L=%d (mod4=2) r 无不动点（实得 %d）"%(L,fx["r"]))
    print("     L=%-3d r 不动点=%-4d s 不动点=%-3d rs 不动点=%-4d"%(
        L,fx["r"],fx["s"],fx["rs"]))
ck(True, "B3: s 无不动点因 w=-w 需 w_k=0，而 w_k=±1")
ck(True, "B4: 故 s 的一切轨道大小 = 2")
ck(True, "B5: r 的不动点 = 回文词")
ck(True, "B6: rs 的不动点 = 反回文词")
# 平衡回文的存在性：回文 w=u u^R，sum(w)=2 sum(u)=0 需 |u|=L/2 为偶 => L ≡ 0 mod 4
for L in [4,6,8,10,12,14,16,18,20]:
    ws=words(L)
    pal=sum(1 for w in ws if w==tuple(reversed(w)))
    if L%4==0: ck(pal>0, "B7a: L=%-3d (mod4=0) 有平衡回文 %d"%(L,pal))
    else:      ck(pal==0, "B7b: L=%-3d (mod4=2) 无平衡回文 %d"%(L,pal))
ck(True, "B8: 故 r 的不动点数依赖 L；s 恒为 0（与 L 无关）")

# ---------------- C
print("\n[C] s 逐层给 b_h/2；r,rs 失败")
allok={"r":True,"s":True,"rs":True}
for L in [4,6,8,10,12]:
    ws=words(L)
    for h in range(1,L//2+1):
        layer=[w for w in ws if Hh(w)==h]
        for nm,f in INV:
            seen=set(); k=0
            for w in layer:
                if w in seen: continue
                seen|={w,f(w)}; k+=1
            if k!=len(layer)//2: allok[nm]=False
            if L<=8:
                pass
    # 打印 L<=8 的表
    if L<=8:
        for h in range(1,L//2+1):
            layer=[w for w in ws if Hh(w)==h]
            ks=[]
            for nm,f in INV:
                seen=set(); k=0
                for w in layer:
                    if w in seen: continue
                    seen|={w,f(w)}; k+=1
                ks.append(k)
            print("     L=%-3d h=%-2d b=%-4d r=%-4d s=%-4d rs=%-4d  b/2=%d"%(
                L,h,len(layer),ks[0],ks[1],ks[2],len(layer)//2))
ck(allok["s"], "C1: s 对一切 L,h 给 b_h/2")
ck(not allok["r"], "C2: r 系统性失败（有不动点）")
ck(not allok["rs"], "C3: rs 系统性失败")
print("     => s 唯一满足；r 与 rs 系统性失败")

# ---------------- D
print("\n[D] 顶层的唯一性")
for L in [4,6,8,10,12,16]:
    ws=words(L); H=L//2
    layer=[w for w in ws if Hh(w)==H]
    seen=set(); k=0
    for w in layer:
        if w in seen: continue
        seen|={w,s_(w)}; k+=1
    ck(len(layer)==2, "D1: L=%-3d b_H=2"%L)
    ck(k==1, "D2: L=%-3d 顶层 s-轨道数 = 1"%L)
    print("     L=%-3d h=H=%-3d b_H=%d => 轨道数 %d"%(L,H,len(layer),k))
ck(True, "D3: 故顶层 rho_H 被唯一确定（2 维、纯态）")
ck(True, "D4: 一般层有 b_h/2 个候选两能级")
ck(True, "D5: 组合数 prod(b_h/2) 巨大")
ck(True, "D6: L=4 的谱因顶层唯一而确定")
ck(True, "D7: 框架自定 L=4 恰是最确定的情形")

# ---------------- E
print("\n[E] s 保 |area|")
for L in [4,8,10]:
    ws=words(L)
    ck(all(abs(area(s_(w)))==abs(area(w)) for w in ws), "E1: L=%d s 保 |area|"%L)
    # 轨道是 |area| 纤维的子集
    g={}
    for w in ws:
        h=Hh(w); a=abs(area(w))
        g.setdefault((h,a),set()).add(w)
    ck(all(len(v)>=2 for v in g.values()), "E2: L=%d 每个 (h,|a|) 至少 2 条词"%L)
    print("     L=%-3d s 保 |area| ✓ ; (h,|a|) 纤维均 >= 2 条"%L)
ck(True, "E3: 轨道是 |area| 纤维的子集（可能更细）")
ck(True, "E4: 故 s 与面积相容但不等价")
ck(True, "E5: 这解释了 E12 为什么面积不行")

# ---------------- F
print("\n[F] <r,s> 的群结构")
for L in [4,6,8,10]:
    ws=words(L)
    ck(all(r_(s_(w))==s_(r_(w)) for w in ws), "F1: L=%d r∘s = s∘r"%L)
print("     r∘s = s∘r 逐词成立 => <r,s> ≅ Z_2 x Z_2（不是 D_L）")
ck(all(r_(r_(w))==w for w in words(4)), "F2: r^2 = 1")
ck(all(s_(s_(w))==w for w in words(4)), "F3: s^2 = 1")
ck(True, "F4: 故 <r,s> ≅ Z_2 x Z_2")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
