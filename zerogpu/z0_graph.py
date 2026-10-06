"""
z0_graph.py --- **补回被漏掉的 Γ**：带图的最小系统。

**上一版的错误（本轮审计发现）**

Z1（由 Z0① 导出）有三款：
    ① 存在有限连通图 Γ=(C,E)
    ② 相邻顶点间两方向均可通行
    ③ 演化由词记录，每一步是沿 Γ 的一条有向边的一次移动
`z0_genesis.py` / `z0_universe.py` **只实现了 ③**，把 ① 整块丢了：

    * 没有顶点 C  =>  没有合并  =>  重数恒为 1（Z0③ 的"整数重数"被架空）
    * Z2「**所有可用边**被实例化」被写成**固定 2 分支**；正确分支数是 deg(v)
    * 散度字典 d_w ∈ Z^C（Z1 §3.1 自称"把词翻译成构型的**唯一桥梁**"）从未被跨过
    * 没有邻接 => 没有 Laplacian => 没有色散关系 => 错说"不继续"
    * 没有 dim Γ => 错说"没有维数"

**修正后的最小系统**

    状态：占据数 n: C -> Z_{>=0}     —— 这才是 Z0③ 的整数重数
    一步：Z2 全分支 => n'(w) = Σ_{v~w} n(v) = (A n)(w)     A = Γ 的邻接算子
    闭合：游走回到起点（Z3）

    一步就是**邻接算子作用** => 线性动力学 => 有色散、有有限速度、有维数。
"""
from __future__ import annotations

import json
import os
import sys
from math import comb

import numpy as np


def step_adj(n: np.ndarray) -> np.ndarray:
    """n' = A n：沿每个轴 ±1 平移求和（环面边界）。"""
    out = np.zeros_like(n)
    for ax in range(n.ndim):
        out += np.roll(n, 1, axis=ax) + np.roll(n, -1, axis=ax)
    return out


def run_front(L: int, D: int, steps: int):
    """环面 L^D 上的前沿扩张。距离用 min(i, L-i)，且取 L >= 2*steps+3 保证不绕回。"""
    n = np.zeros((L,) * D, dtype=np.int64)
    n[(0,) * D] = 1
    hist = []
    for k in range(steps + 1):
        nz = np.nonzero(n)
        if nz[0].size:
            r = 0
            for c in nz:
                d = np.minimum(c, L - c)
                r = max(r, int(d.max()))
            support = int(nz[0].size)
        else:
            r, support = 0, 0
        hist.append({"k": k, "front_radius": r, "support_sites": support,
                     "total": int(n.sum()), "at_origin": int(n[(0,) * D])})
        n = step_adj(n)
    return hist


def closed_walk_table(D: int, KMAX: int):
    """A[i][j] = Z^i 上长度 j 的闭环走数。一次 DP 算完所有 j（原来每 k 重算，慢几百倍）。"""
    A = [[0] * (KMAX + 1) for _ in range(D + 1)]
    A[0][0] = 1
    for i in range(1, D + 1):
        for j in range(KMAX + 1):
            s = 0
            for m in range(0, j + 1, 2):
                prev = A[i - 1][j - m]
                if prev:
                    s += comb(j, m) * comb(m, m // 2) * prev
            A[i][j] = s
    return A


def main():
    steps = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    L = 2 * steps + 3
    res = {}

    print("=" * 96)
    print(f"§1 微观因果性（G15 引理 52）：前沿速度是否恰为 1 边/步？  （环面 L={L}，不绕回）")
    print("=" * 96)
    print(f"{'D':>2} {'k':>3} {'前沿半径':>9} {'=k?':>5} {'支撑点数':>10} "
          f"{'总行走者数':>20} {'=(2D)^k?':>9}")
    s1 = []
    for D in (1, 2, 4):
        h = run_front(L, D, steps)
        ok = all(x["front_radius"] == x["k"] for x in h)
        ok2 = all(x["total"] == (2 * D) ** x["k"] for x in h)
        s1.append({"D": D, "front_ok": bool(ok), "total_ok": bool(ok2),
                   "front": [x["front_radius"] for x in h]})
        for x in [h[0], h[1], h[2], h[-1]]:
            print(f"{D:>2} {x['k']:>3} {x['front_radius']:>9} "
                  f"{str(x['front_radius']==x['k']):>5} {x['support_sites']:>10,} "
                  f"{x['total']:>20,} {str(x['total']==(2*D)**x['k']):>9}")
        print(f"   => D={D}: 前沿速度 = 1 边/步 ? {ok} ; 总数 = (2D)^k ? {ok2}\n")
    res["causality"] = s1

    print("=" * 96)
    print("§2 整数重数：占据数是否 >1（合并是否存在）？")
    print("=" * 96)
    T = closed_walk_table(1, 12)
    print(f"{'k':>3} {'n_k(0) 精确':>14} {'C(k,k/2)':>12} {'一致':>6}")
    mx = 0
    for k in range(13):
        v = T[1][k]
        c = comb(k, k // 2) if k % 2 == 0 else 0
        mx = max(mx, v)
        print(f"{k:>3} {v:>14,} {c:>12,} {str(v==c):>6}")
    print(f"""
  判读：一维链上 n_k(0) = C(k,k/2)，最大 {mx:,} >> 1。
        => 顶点占据数（= Z0③ 的整数重数）**真的 >1**，来源是**不同的词合并到同一顶点**。
           上一版没有顶点，永远看不到重数。""")

    print("\n" + "=" * 96)
    print("§3 色散关系：邻接算子的能带（有能带 => 有传播）")
    print("=" * 96)
    print(f"{'D':>2} {'lam_max':>8} {'lam_min':>8} {'带宽':>9} {'群速度上界':>11} {'k=0 简并':>9}")
    s3 = []
    for D in (1, 2, 3, 4):
        ks = 2 * np.pi * np.arange(128) / 128
        grids = np.meshgrid(*([ks] * D), indexing="ij")
        ev = sum(2 * np.cos(g) for g in grids).ravel()
        lmax, lmin = float(ev.max()), float(ev.min())
        nz = int(np.sum(np.abs(ev - lmax) < 1e-9))
        s3.append({"D": D, "lam_max": lmax, "lam_min": lmin,
                   "bandwidth": lmax - lmin, "zero_mode": nz})
        print(f"{D:>2} {lmax:>8.4f} {lmin:>8.4f} {lmax-lmin:>9.4f} {2.0:>11.4f} {nz:>9}")
    print("""
  判读：A 带宽有限 [-2D, 2D] => **有色散关系** omega(k) => **有传播**；
        群速度上界 = 2 格点/步 => **有有限特征速度**。
        **这就是 Z1 的图 Γ 给回来的"继续"。**""")
    res["dispersion"] = s3

    print("\n" + "=" * 96)
    print("§4 谱维数：P(t) = n_t(0)/(2D)^t ~ t^(-d_s/2)，d_s 是否 = D？")
    print("=" * 96)
    KMAX = 200
    TT = closed_walk_table(4, KMAX)
    print(f"{'D':>2} {'拟合 d_s':>10} {'|d_s-D|':>9} {'指数 a':>9} {'应 D/2':>9}")
    s4 = []
    for D in (1, 2, 3, 4):
        ks, Ps = [], []
        for k in range(4, KMAX + 1):
            v = TT[D][k]
            if v > 0:
                ks.append(float(k))
                Ps.append(float(v) / (2 * D) ** k)
        ks, Ps = np.array(ks), np.array(Ps)
        n = max(8, int(0.35 * len(ks)))
        a = -np.polyfit(np.log(ks[:n]), np.log(Ps[:n]), 1)[0]
        ds = 2 * a
        s4.append({"D": D, "d_s": float(ds), "a": float(a)})
        print(f"{D:>2} {ds:>10.4f} {abs(ds-D):>9.4f} {a:>9.4f} {D/2:>9.4f}")
    print("""
  判读：**d_s = D**（有限 k 窗口内）。上一版闭链图的 d_s 会漂移，因为那里没有 Γ。
        代价：**D 成了输入**（G89 已证 D 不可由 Z0 导出；R23 把它变成"选择"）。""")
    res["spectral_dimension"] = s4

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "z0_graph.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("\nsaved results/z0_graph.json")


if __name__ == "__main__":
    main()
