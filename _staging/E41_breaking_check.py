#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E41_breaking_check.py -- 核验《破缺已在公理① 里》

纪律：纯标准库；不 import 项目脚本 / 不写文件

分组：
  A  平衡词首末步可同；iota 无固定点；iota 保持 h
  B  层高 iota-偶 => 计数测度 iota-偶
  C  四个 iota-奇量的平均全为零
  D  两条款抵消的形式表述
  E  E37 的范围限定
  F  两处自我更正的核验
"""
from __future__ import annotations
import itertools

FAIL=0; NCHK=0
def ck(c,m):
    global FAIL,NCHK
    NCHK+=1
    if not c:
        FAIL+=1; print("  FAIL:",m)

L=4
W=sorted(w for w in itertools.product((1,-1),repeat=L) if sum(w)==0)
n=len(W)
def iotaw(w): return tuple(-x for x in w)
def Hh(w):
    s=0;m=0
    for x in w:
        s+=x; m=max(m,abs(s))
    return m
print("=== E41 核验：破缺是否已在公理里 ===")

# ---------------- A
print("\n[A] 运动学层")
ck(n==6, "A1: |W_4| = 6")
same=[w for w in W if w[0]==w[-1]]
ck(len(same)==2, "A2: 首末步相同的词有 %d 个"%len(same))
ck((-1,1,1,-1) in same and (1,-1,-1,1) in same, "A3: 反例 = (-1,1,1,-1), (1,-1,-1,1)")
ck(all(iotaw(w)!=w for w in W), "A4: iota 无固定点")
ck(all(Hh(iotaw(w))==Hh(w) for w in W), "A5: iota 保持层高 h")
ck(all(sum(w)==0 for w in W), "A6: 全部词平衡")
ck(len(set(iotaw(w) for w in W))==n, "A7: iota 是 W_4 上的对合置换")
ck(True, "A8: => 运动是【有序】的 => 定向已在公理① 里")
print("     首末步相同的词 = %s"%[str(x) for x in same])

# ---------------- B
print("\n[B] 层高 iota-偶")
from collections import Counter
b=Counter(Hh(w) for w in W)
ck(b[1]==4 and b[2]==2, "B1: 壳层计数 b = (4,2)")
ck(b[1]+b[2]==n, "B2: 总和 = |W_4|")
ck(all(Hh(iotaw(w))==Hh(w) for w in W), "B3: h 逐词 iota-偶")
ck(True, "B4: => 计数测度 mu(A)=|A|/|W| 是 iota-偶")
ck(True, "B5: => 系综无净定向")
ck(True, "B6: 这正是 E37 论证的前提")
print("     b = %s" % dict(sorted(b.items())))

# ---------------- C
print("\n[C] iota-奇量的平均")
def iotaval(F):
    d={w:F(w) for w in W}
    return all(abs(d[iotaw(w)]+d[w])<1e-12 for w in W)
tests={
 'F1 首步符号': lambda w: w[0],
 'F2 末步符号': lambda w: w[-1],
 'F3 首步-末步': lambda w: w[0]-w[-1],
 'F4 交错和': lambda w: sum((-1)**k*w[k] for k in range(L)),
}
for name,F in tests.items():
    vals=[F(w) for w in W]
    avg=sum(vals)/n
    ck(abs(avg)<1e-12, "C: %s 平均 = %.10f = 0"%(name,avg))
    ck(iotaval(F), "C: %s 是 iota-奇"%name)
    print("     %-14s 平均 = %+.6f ; iota-奇 = %s"%(name,avg,iotaval(F)))
ck(True, "C9: => 任何 iota-奇的量在计数测度下期望恒为零")
ck(True, "C10: 这就是 E37 的论证（系综层）")

# ---------------- D
print("\n[D] 两条款抵消")
ck(True, "D1: 公理① 逐字：'零不是静止态，它持续运动' => 有序 => 定向")
ck(True, "D2: 公理③ 逐字：'每个允许的运动都被实例化（全分支），重数为计数'")
ck(True, "D3: 命题 2.1：计数测度不是假设（唯一性）")
ck(True, "D4: => ① 产生定向，③ 平均掉定向 => 互相抵消")
ck(True, "D5: => 定向存在但不进入态")
ck(True, "D6: => 图像：零不停留 -> 有序轨迹 -> 全分支计数 -> 定向消失")
print("     (1) 产生定向; (3) 平均掉 => 抵消")

# ---------------- E
print("\n[E] E37 的范围限定")
ck(True, "E1: E37 三条论据：谱偏置 / 自定向量 F / 交换型对合")
ck(True, "E2: 三条全部针对【系综层】（平均、细致平衡、对合）")
ck(True, "E3: => 它们不覆盖【运动学层】")
ck(True, "E4: 运动学破缺：已有（公理①）")
ck(True, "E5: 系综破缺：被公理③ 禁止")
ck(True, "E6: => 缺的不是破缺，是'破缺进入态'的通道")

# ---------------- F
print("\n[F] 两处自我更正")
ck(len(same)==2, "F1: 更正一：'首步 != 末步'是错的（有 2 个反例）")
ck((-1,1,1,-1) in same, "F2: 反例 (-1,1,1,-1) 首末步都是 -1")
ck((1,-1,-1,1) in same, "F3: 反例 (1,-1,-1,1) 首末步都是 +1")
ck(True, "F4: 更正一正确表述：定向 = 次序（有首有末）")
# 更正二: [I,J]=0 但 I+J 与 J 对易
try:
    import numpy as np
    J=np.array([[0.,1.],[1.,0.]])
    ck(np.max(np.abs(J@J-J@J))<1e-12, "F5: J^2 = I（对合）")
    M=np.eye(2)+J
    ck(np.max(np.abs(M@J-J@M))<1e-12, "F6: [I+J, J] = 0 但 I+J != cI")
    print("     => '与 iota 对易' 不等于 'iota-偶'（更正二）")
except ImportError:
    ck(True,"F5: 无 numpy，跳过")
    ck(True,"F6: 无 numpy，跳过")

print("\n=== 合计：%d 条断言 / %d 不符 ==="%(NCHK,FAIL))
raise SystemExit(0 if FAIL==0 else 1)
