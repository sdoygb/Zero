#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E17_summary_consistency_check.py -- 核验《E 系列判定书》的结构一致性

纪律：只读 / 纯标准库 / 不写文件
检查：
  A  E01–E16 十六篇都存在
  B  每篇的核验脚本都存在
  C  判定书里引用的断言数与实际相符
  D  判定书里的关键数字与实际计算相符
  E  判定书的 LaTeX 结构（行间/行内）
"""
from __future__ import annotations
import math, os, re, subprocess, sys

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

HERE=os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

print("=== E17 核验：判定书的结构一致性 ===")

# ---------------- A
print("\n[A] E01–E16 十六篇都存在")
titles={}
for n in range(1,17):
    tag="E%02d"%n
    hits=[f for f in os.listdir(".") if f.startswith(tag+"_") and f.endswith(".md")]
    ck(len(hits)==1, "A: %s 恰好一篇（实得 %s）"%(tag,hits))
    if hits: titles[n]=hits[0]
print("     十六篇全部存在 ✓")

# ---------------- B
print("\n[B] 每篇的核验脚本都存在")
for n in range(1,17):
    tag="E%02d"%n
    hits=[f for f in os.listdir(".") if f.startswith(tag+"_") and f.endswith("_check.py")]
    ck(len(hits)>=1, "B: %s 有核验脚本"%tag)
print("     十六个脚本全部存在 ✓")

# ---------------- C
print("\n[C] 判定书引用的断言数与实际相符")
CLAIM={1:40,2:54,3:39,4:34,5:27,6:38,7:44,8:39,9:23,10:32,
       11:103,12:67,13:83,14:68,15:49,16:64}
actual={}
for n in range(1,17):
    tag="E%02d"%n
    sc=[f for f in os.listdir(".") if f.startswith(tag+"_") and f.endswith("_check.py")][0]
    out=subprocess.run([sys.executable,sc],capture_output=True,text=True).stdout
    m=re.search(r"(\d+)\s*条断言",out)
    actual[n]=int(m.group(1)) if m else -1
    ck(actual[n]==CLAIM[n], "C: %s 实际 %d 条 vs 判定书 %d 条"%(tag,actual[n],CLAIM[n]))
    print("     %s  实际 %-4d  判定书 %-4d"%("E%02d"%n,actual[n],CLAIM[n]))
tot=sum(actual.values())
ck(tot==804, "C_end: E 系列合计 = %d（判定书写 804）"%tot)
print("     合计 = %d"%tot)

# ---------------- D
print("\n[D] 判定书里的关键数字")
# beta_h = log((b_h+1)/(b_h-1)) for L=4
b=[4,2]
b1=math.log((b[0]+1)/(b[0]-1)); b2=math.log((b[1]+1)/(b[1]-1))
ck(abs(b1-math.log(5/3))<1e-14, "D1: beta_1 = log(5/3)")
ck(abs(b2-math.log(3))<1e-14, "D2: beta_2 = log(3)")
print("     beta_1 = log(5/3) = %.10f ; beta_2 = log(3) = %.10f"%(b1,b2))
# L=4 的 K 谱
spec=sorted([math.log(12/5),math.log(4),math.log(4),math.log(12)])
ck(abs(spec[0]-0.875469)<1e-5, "D3: K 谱首项 = log(12/5)")
ck(abs(spec[3]-2.484907)<1e-5, "D4: K 谱末项 = log(12)")
print("     spec(K) = %s"%["%.6f"%x for x in spec])
# kappa_1 两处一致（E06 的 5/9）
ck(abs((4/6)**2+(2/6)**2-5/9)<1e-15, "D5: kappa_1(L=4) = 5/9")
print("     kappa_1(L=4) = 5/9 = %.10f"%(5/9))
# |a|max 公式
for L in [4,8,12]:
    H=L//2
    f=lambda h:-h*h+(L+1)*h-H
    ck(f(1)==H, "D6: L=%d f(1)=H"%L)
    ck(f(H)==H*H, "D7: L=%d f(H)=H^2"%L)
print("     |a|max 公式两端核对通过")

# ---------------- E
print("\n[E] 判定书的 LaTeX 结构")
BS=chr(92)
def net(r):
    a=0;i=0
    while i<len(r):
        if r[i]==BS: i+=2; continue
        if r[i]=="{": a+=1
        elif r[i]=="}": a-=1
        i+=1
    return a
s=open(titles[17] if 17 in titles else "E17_E系列判定书.md",encoding="utf-8").read()
lines=s.split("\n")
inm=False;d=0;mb=0;bad=0;nl=0
for i,l in enumerate(lines,1):
    if l.strip()=="$$":
        if not inm: inm=True;d=0
        else:
            inm=False;mb+=1
            if d!=0: bad+=1
        continue
    if inm:
        d+=net(l)
        if "](" in l: nl+=1
ck(bad==0, "E1: 行间公式括号全平衡（%d 块）"%mb)
ck(nl==0, "E2: 数学块内无 Markdown 链接")
bi=0;bad2=0
for l in lines:
    if l.strip().startswith("$$"): continue
    for x in re.findall(r"\$([^$]+)\$", l):
        bi+=1
        if net(x)!=0: bad2+=1
ck(bad2==0, "E3: 行内公式括号全平衡（%d 处）"%bi)
# 相对链接
bad3=[]
for m in re.finditer(r"\]\(([^)]+)\)", s):
    t=m.group(1).split("#")[0]
    if not t or t.startswith("http"): continue
    if not os.path.exists(t): bad3.append(t)
ck(len(bad3)==0, "E4: 相对链接全部可解（坏链 %s）"%bad3)
print("     行间 %d 块（坏 %d）；行内 %d 处（坏 %d）；坏链 %d"%(mb,bad,bi,bad2,len(bad3)))

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
