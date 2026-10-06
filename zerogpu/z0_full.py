"""
z0_full.py --- 零和宇宙完整引擎（L0/L1/L2 + 层塔 + 模流）
链条: q_L -> 层权重 -> omega -> 模频率 -> 引擎(词/账本/活动) + 模相位(g) -> 测周期
状态: 跑通; L=4(被R31唯一选出) 下实测周期约 3-4; 瓶颈 = 状态空间只有 4 个闭合词
"""
import math, cmath, numpy as np
from fractions import Fraction
from itertools import combinations
from collections import Counter, defaultdict

L = 4                                   # R31 唯一选择定理选出的寿命

def qL_exact(L):
    seen = {}
    for ones in combinations(range(L), L//2):
        w = [1]*L
        for i in ones: w[i] = -1
        c = min(tuple(w[i:]+w[:i]) for i in range(L)); seen[c] = seen.get(c, 0) + 1
    return sum(Fraction(o, math.comb(L, L//2))**2 * n for o, n in Counter(seen.values()).items())

def prefix_weights(L):
    out = {}
    for k in range(L+1):
        for s in range(-k, k+1, 2):
            rem = L - k
            if (rem - s) % 2: continue
            a = (rem - s)//2
            if 0 <= a <= rem: out[(k, s)] = math.comb(rem, a)
    return out

def build_tower(L):
    q = qL_exact(L); qf = float(q); blk = Counter()
    for (k, s), W in prefix_weights(L).items(): blk[abs(s)] += W * (qf ** abs(s))
    hs = sorted(blk); ws = np.array([blk[h] for h in hs], float); ws /= ws.sum()
    freq = np.array([-math.log(w/2) for w in ws])       # 模频率 -log(w_a/2)
    return q, hs, ws, freq

def bal(w): return sum(w)

def run(T=4000, verbose=True):
    q, hs, ws, freq = build_tower(L)
    L0 = Counter(); L1 = defaultdict(list); L2 = Counter(); L2[()] = 1
    P = []; g = 0
    for t in range(1, T+1):
        newL2 = Counter(); amp = np.zeros(len(hs), dtype=complex); nclo = 0
        for w, mult in L2.items():
            for d in (+1, -1):
                nw = w + (d,)
                if bal(nw) == 0 and len(nw) >= 2:            # 闭合 ⟺ 平衡归零
                    nclo += 1; L0[(nw, d)] += mult           # L0 只增
                    h = max(abs(bal(nw[:k+1])) for k in range(len(nw)))   # 远足高度
                    hh = hs.index(h) if h in hs else len(hs)-1
                    L1[nw] = sorted(L1[nw]+[mult], reverse=True)[:2]      # L1 留两层
                    amp[hh] += mult * ws[hh] * cmath.exp(-1j*freq[hh]*g)  # 模相位(全局 g)
                elif len(nw) < L: newL2[nw] += mult
        if nclo > 0:
            newL2[()] += nclo; g += 1; P.append(abs(amp.sum())**2)
        L2 = newL2
    P = np.array(P)
    if verbose:
        print(f"L={L} q_4={q} 层高={hs}")
        print(f"层权重={np.round(ws,5).tolist()} 模频率={np.round(freq,4).tolist()}")
        print(f"闭合代数 g={g} L0累计={sum(L0.values())} L1词数={len(L1)} 末活动量={sum(L2.values())}")
        print(f"P: 均值={P.mean():.4f} 起伏={P.std()/P.mean():.4f}")
        ac = [float(((P[:-k]-P.mean())*(P[k:]-P.mean())).mean()/P.var()) for k in (1,3,10,100,500)]
        print(f"自相关(1,3,10,100,500)={[round(a,3) for a in ac]}")
    return P, L0, L1

if __name__ == "__main__":
    run()
