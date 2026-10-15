#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
force_from_field.py —— 场→力的第二跳
=====================================
E31 §8 的精确边界：力 = 势的梯度(Zero 有 E=fᵀLf) + 时间演化律(Zero 无 Z005)，
所以静态场分布"给了方向，没给力"。第二跳尝试补缺口的两条路：

路A【可算·标量约束力】：把场的偏置读成"约束回复力"标量场
      D(h,k) = p↓ - p↑ = 1 - 2·p↑(h,k)
  含义：闭合(零和)约束把游走"拉回"原点(低层)的回复强度。无向计数，
       绕开定向墙(E31 §8.6)。检查：壁处 D=1(全向下)、对称点 D=0、
       内部 0<D<1；并检查 D 是否只依赖 h(可作势梯度)还是依赖 k。
路B【演化律缺口】：力缺的那一半是"时间演化律"。用户主线的 L0-L1-L2
       双向动力学(Z088 六边闭环)正是演化律的候选——力 = 势的梯度 + 该演化律。

诚实边界：D 是标量强度分布，不是矢量力/规范力；"力=曲率/和乐"仍被定向墙挡。
"""
import sys, json
import numpy as np
from field_sim import field_at

OUT = "/Users/oygb/Downloads/Zero/ZeroGPU/results/force_from_field.json"

def main(L=16):
    up, tot = field_at(L)
    keys = sorted(tot.keys())
    # 计算 D(h,k)=1-2p↑，以及每层/格点数据
    data = {}
    maxD = 0; maxD_key = None
    sym = None            # D=0 的内部对称点
    for (h, k) in keys:
        u = up[(h,k)]; t = tot[(h,k)]; p = u/t
        D = 1.0 - 2.0*p
        is_wall = (k == h) or (abs(p) < 1e-9)
        data[(h,k)] = dict(p=p, D=D, u=int(u), t=int(t), wall=bool(is_wall))
        if not is_wall:
            if abs(D) < 1e-9: sym = (h,k)
            if abs(D) > maxD: maxD, maxD_key = abs(D), (h,k)

    print("="*72)
    print(f"场→力：约束回复力标量场 D(h,k)=1-2·p↑，L={L}")
    print("="*72)
    print(f"{'层h':<5}{'高度k':<6}{'p↑':<9}{'D=1-2p↑':<10}{'属性'}")
    print("-"*56)
    for (h,k) in sorted(keys):
        d = data[(h,k)]
        tag = "壁D=1" if d['wall'] else ("对称D=0" if abs(d['D'])<1e-9 else "")
        print(f"{h:<5}{k:<6}{d['p']:<9.4f}{d['D']:<10.4f}{tag}")

    print("-"*56)
    print(f"载体点数={len(keys)}，壁点={sum(1 for d in data.values() if d['wall'])}")
    print(f"最大 |D|={maxD:.4f} @{maxD_key}（壁处 D=1）")
    print(f"内部对称点 D=0: {sym}")

    # D 是否只依赖 h？(可作势梯度)
    h_dep_k = False
    for h in sorted(set(k[0] for k in keys)):
        dd = {k: data[(h,k)]['D'] for (hh,k) in keys if hh==h}
        vals = set(round(v,6) for v in dd.values())
        if len(vals) > 1: h_dep_k = True
    print(f"D 依赖 k(同层不同值)? {h_dep_k}")
    print("  ⇒ 若依赖 k：D 不是纯 h 的势的梯度(需第二分量才构成 2D 保守力)")

    # 每层 D 的均值(回复力强度随层变化)
    print("\n每层平均回复力 <D>_h：")
    for h in sorted(set(k[0] for k in keys)):
        dd = [data[(h,k)]['D'] for (hh,k) in keys if hh==h and not data[(h,k)]['wall']]
        print(f"  层{h}: <D>={np.mean(dd):.4f}  (内部 {len(dd)} 点)")

    result = {
        "L": L, "n_carrier": len(keys), "wall": sum(1 for d in data.values() if d['wall']),
        "max_D": maxD, "max_D_key": list(maxD_key), "sym_D0": list(sym) if sym else None,
        "D_depends_on_k": h_dep_k,
        "grid": {f"({h},{k})": {"p↑": round(data[(h,k)]['p'],6),
                                "D": round(data[(h,k)]['D'],6),
                                "壁": data[(h,k)]['wall']}
                 for (h,k) in sorted(keys)},
        "conclusion": "D 是标量约束回复力场(闭合约束拉回原点); 非矢量/规范力",
        "boundary": "力=曲率/和乐仍被定向墙挡; 力缺的演化律指向 L0-L1-L2 双向动力学(Z088)"
    }
    with open(OUT,"w") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)
    print(f"\n结果已写入 {OUT}")

if __name__ == "__main__":
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    main(L)
