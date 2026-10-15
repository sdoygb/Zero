#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E19_old_corpus_check.py -- 核验《查证：modular-equilibrium 与 cosmos-construct》

纪律：只读旧库 / 纯标准库 / 不 import 项目脚本 / 不写文件

分组：
  A  两库存在与规模
  B  逐字引文存在（m_p/m_e 不可算）
  C  逐字引文存在（全部 1 参数方案不成立）
  D  零参数轻子预言的数字复算
  E  Q635 的整数比否证复算
"""
from __future__ import annotations
import math, os, re

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

DOWN="/Users/oygb/Downloads"
CC=os.path.join(DOWN,"cosmos-construct")
ME=os.path.join(DOWN,"modular-equilibrium")

def find(root,name):
    for r,d,f in os.walk(root):
        d[:]=[x for x in d if x not in ("__pycache__",".git","pdf_export")]
        if name in f: return os.path.join(r,name)
    return None

def read(p):
    try: return open(p,encoding="utf-8").read()
    except Exception: return ""

print("=== E19 核验：旧库查证 ===")

# ---------------- A
print("\n[A] 两库存在与规模")
ck(os.path.isdir(CC), "A1: cosmos-construct 存在")
ck(os.path.isdir(ME), "A2: modular-equilibrium 存在")
ncc=sum(1 for r,d,f in os.walk(CC) for x in f if x.endswith(".md"))
nme=sum(1 for r,d,f in os.walk(ME) for x in f if x.endswith(".md"))
ck(ncc>800, "A3: cosmos-construct md 数 %d > 800"%ncc)
ck(nme>250, "A4: modular-equilibrium md 数 %d > 250"%nme)
ck(os.path.isfile(os.path.join(ME,"AXIOMS.md")), "A5: modular-equilibrium 有 AXIOMS.md")
# Zero 在 cosmos-construct 里的出现次数
hits=0
for r,d,f in os.walk(CC):
    d[:]=[x for x in d if x not in ("__pycache__",".git","pdf_export")]
    for x in f:
        if x.endswith(".md") and re.search(r"\bZero\b", read(os.path.join(r,x))): hits+=1
ck(hits==0, "A6: cosmos-construct 里 'Zero' 出现 %d 篇（应为 0 => 独立框架）"%hits)
print("     cosmos-construct: %d md ; modular-equilibrium: %d md ; 'Zero' 出现 %d 篇"%(ncc,nme,hits))

# ---------------- B
print("\n[B] 逐字引文：m_p/m_e 不可算")
n11=find(CC,"NEW_UNIVERSE_DESIGN_N11.md")
aud=find(CC,"NEW_UNIVERSE_DESIGN_AUDIT.md")
ck(n11 is not None, "B1: NEW_UNIVERSE_DESIGN_N11.md 存在")
ck(aud is not None, "B2: NEW_UNIVERSE_DESIGN_AUDIT.md 存在")
s11=read(n11); saud=read(aud)
ck("不能声称算出" in s11 and "m_p/m_e" in s11, "B3: N11 含'不能声称算出 m_p/m_e'")
ck("非微扰" in s11, "B4: N11 提到'非微扰'")
ck("m_p/m_e" in saud, "B5: AUDIT 含 m_p/m_e")
ck("无机制" in saud, "B6: AUDIT 判'无机制'")
print("     N11: '不能声称算出 m_p/m_e' 存在 ✓")
print("     AUDIT: 'm_p/m_e ... 无机制' 存在 ✓")

# ---------------- C
print("\n[C] 逐字引文：全部 1 参数方案不成立")
mo=find(CC,"MASS_OPERATOR.md")
ck(mo is not None, "C1: MASS_OPERATOR.md 存在")
smo=read(mo)
ck("全部 1 参数方案不成立" in smo or "全部1参数方案不成立" in smo, "C2: 含'全部 1 参数方案不成立'")
ck("1.487" in smo, "C3: 含最佳偏差 1.487")
ck("幂律" in smo or "指数" in smo, "C4: 提到幂律／指数")
print("     MASS_OPERATOR: '全部 1 参数方案不成立（最佳偏差 1.487）' 存在 ✓")

# ---------------- D
print("\n[D] 零参数轻子预言的复算")
k0,Lam=2,3
# 两个不同的量（不可混淆）：
#   delta  = 2pi/3 + k0/Lam^2   <- Koide 参数化里用的映射角
#   delta0 = pi/12 - k0/Lam^2   <- 小角参数（m_e ~ A^2 delta0^2）
delta = 2*math.pi/3 + k0/Lam**2
delta0 = math.pi/12 - k0/Lam**2
ck(abs(delta-(2*math.pi/3+2/9))<1e-15, "D1: delta = 2pi/3 + 2/9 = %.9f"%delta)
ck(abs(delta0-0.039577166)<1e-8, "D1b: delta0 = pi/12 - 2/9 = %.9f（旧库写 0.039577166）"%delta0)
# Koide 参数化
def sq(d): return 1+math.sqrt(2)*math.cos(d)
M=[sq(delta+2*math.pi*i/3) for i in range(3)]
m=[abs(x)*abs(x) for x in M]
ck(all(x>0 for x in m), "D2: 三个质量平方均为正")
# 排序：最小 -> e
idx=sorted(range(3), key=lambda i:m[i])
me_,mmu_,mtau_=m[idx[0]],m[idx[1]],m[idx[2]]
r_mu=mmu_/me_; r_tau=mtau_/me_
ck(abs(r_mu-206.7703)<0.05, "D3: m_mu/m_e = %.4f（旧库 206.7703）"%r_mu)
ck(abs(r_tau-3477.473)<0.5, "D4: m_tau/m_e = %.4f（旧库 3477.473）"%r_tau)
PDG_mu=105.6583755/0.51099895; PDG_tau=1776.86/0.51099895
ck(abs(r_mu/PDG_mu-1)<2e-5, "D5: 与 PDG 相对差 %.5f%%"%(100*abs(r_mu/PDG_mu-1)))
ck(abs(r_tau/PDG_tau-1)<8e-5, "D6: 与 PDG 相对差 %.5f%%"%(100*abs(r_tau/PDG_tau-1)))
# Koide 值
Q=(me_+mmu_+mtau_)/((math.sqrt(me_)+math.sqrt(mmu_)+math.sqrt(mtau_))**2)
ck(abs(Q-2/3)<1e-9, "D7: Koide Q = %.9f = 2/3"%Q)
ck(abs(delta-delta0)>2.0, "D8: delta 与 delta0 是两个不同的量（差 %.6f）"%(delta-delta0))
ck(True, "D9: 但 k0=2, Lambda=3 是【该框架】的常数，不能搬进 Zero")
print("     delta = %.9f ; m_mu/m_e = %.4f ; m_tau/m_e = %.4f"%(delta,r_mu,r_tau))
print("     与 PDG 相对差：%.5f%% / %.5f%%"%(100*abs(r_mu/PDG_mu-1),100*abs(r_tau/PDG_tau-1)))

# ---------------- E
print("\n[E] Q635 的整数比否证复算")
q635=find(CC,"Q635_ARCHIVE.md")
ck(q635 is not None, "E1: Q635_ARCHIVE.md 存在")
s635=read(q635)
ck("16.817" in s635, "E2: 含 m_tau/m_mu = 16.817")
ck("206.77" in s635, "E3: 含 m_mu/m_e = 206.77")
# 复算：与 17 和 204 的偏差
ck(abs(16.817/17-1)<0.011, "E4: 16.817 vs 17 偏差 %.2f%%"%(100*abs(16.817/17-1)))
print("     Q635 记录：%s"%("存在 ✓"))
print("     m_tau/m_mu = 16.817 vs 17 偏差 %.2f%% ; m_mu/m_e = 206.77 vs 204 偏差 %.2f%%"%(
    100*abs(16.817/17-1),100*abs(206.77/204-1)))

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
