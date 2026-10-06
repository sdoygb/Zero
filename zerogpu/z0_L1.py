"""
z0_L1.py --- **真正把 L1 跑起来**：闭合写入精确词与闭合类，然后让 L1 繁殖
依据（语料）：写入律 lambda: w -> [w]（R89）；循环序/旋转类 = Z-E*（Z14）；五条子代规则（繁殖定理）
L1 状态 = {(精确词 w, 闭合类 [w], 重数)}；L2 = 环上忠实行走，闭合事件驱动 L1 写入。
"""
import numpy as np, json, os
from itertools import combinations
from collections import Counter, defaultdict
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),"results")

def mo(s):
    T=len(s); return min(tuple(s[i:]+s[:i]) for i in range(T))
def cuts(w):
    s=0;o=[]
    for i,x in enumerate(w):
        s+=x
        if s==0 and 0<i+1<len(w): o.append(i+1)
    return o
def ch(w,r):
    c=cuts(w)
    if r=="persist": return [w]
    if r=="copy": return [w,w]
    if r=="split_parent":
        o=[w]
        for k in c: o+=[mo(w[:k]),mo(w[k:])]
        return o
    if r=="split_only":
        o=[]
        for k in c: o+=[mo(w[:k]),mo(w[k:])]
        return o
    o=[]
    for mask in range(1<<len(c)):
        ks=[c[i] for i in range(len(c)) if mask>>i&1]; prev=0
        for k in ks+[len(w)]: o.append(mo(w[prev:k])); prev=k
    return o

# ---------- L2：环上的忠实行走（全分支，归一化），闭合事件写 L1 ----------
def run_L2(L=12, steps=400):
    N=np.zeros(L); N[0]=1.0          # 从 0 出发的分支分布
    start_pos=np.zeros(0)             # 简化：只跟踪"从 0 出发"的闭合
    L1=Counter(); L1cls=Counter(); log=[]
    for t in range(1,steps+1):
        N=np.array([N[(i-1)%L]/2 + N[(i+1)%L]/2 for i in range(L)])   # 每步向两侧分支
        N/=N.max()
        c=N[0]                        # 回到起点的量 = 闭合发生率
        if c>1e-12:
            # 闭合词：长度 t 的平衡词集合（全分支 ⟹ 全部实现），按路径均匀 ⟹ 权重 ∝ 轨道大小
            for ones in combinations(range(t), t//2):
                w=[1]*t
                for i in ones: w[i]=-1
                w=tuple(w)
                if w[0]!=-1 and mo(w)!=w: continue      # 只取规范代表（旋转类）
                L1[w]+=1; L1cls[mo(w)]+=1
            log.append((t, len(L1), len(L1cls)))
    return L1,L1cls,log

L1,L1cls,log=run_L2(L=12, steps=12)
print("### L1 被 L2 的闭合事件写入了什么（环 L=12，长度<=12 的闭合）")
print(f"{'长度t':>6} {'L1 精确词数':>12} {'L1 闭合类数':>12}")
for t,a,b in log: print(f"{t:>6} {a:>12} {b:>12}")
print(f"\nL1 总计：精确词 {len(L1)} 个（重数总和 {sum(L1.values())}） 闭合类 {len(L1cls)} 个")
rs=Counter(len(cuts(w)) for w in L1cls)
print(f"闭合类的内部归零数分布 r: {dict(sorted(rs.items()))}")
print(f"其中可切分的（r>=1）= {sum(v for k,v in rs.items() if k>=1)} / {len(L1cls)}")

# ---------- L1 繁殖：从**可切分**的模式起跑 ----------
print("\n### L1 繁殖（起点选**可切分**的模式，避免上次 ++-- 无切点的错误）")
for seedword in ((1,-1,1,-1), (1,1,-1,-1), (1,-1,-1,1,1,-1)):
    w=mo(seedword); r0=len(cuts(w))
    print(f"\n  起点 {''.join('+' if x>0 else '-' for x in w)}  r={r0}")
    for rule in ("split_parent","split_only","split_all","copy"):
        cur=Counter({w:1}); seen={w}
        print(f"    {rule:<13}", end="")
        for g in range(6):
            nxt=Counter()
            for p,mult in cur.items():
                for c in ch(p,rule): nxt[c]+=mult
            cur=nxt; seen|=set(cur)
            print(f" g{g+1}:{len(cur)}模式/{sum(cur.values())}条", end="")
        print(f"   → 6 代内见过 {len(seen)} 个模式")
json.dump({"L1_words":len(L1),"L1_classes":len(L1cls),"r_dist":dict(rs)},
          open(os.path.join(OUT,"z0_L1.json"),"w"),ensure_ascii=False,indent=1)
