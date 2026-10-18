#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""16/8/4 的两条路线对比：
A: L=8, H=4, 词类 (h,σ,π)，top2/类；ρ 由 L1 词记录给。
B: L=4, K4 4 基，16/8/4 是 4 基×纤维的 4/2/1 槽（分级 profile，按观测粗分组标定）。

两者跑同一个三层 push-forward 模型（毁灭/flow/reseed 相同），比较相位平均和周期 4。
输出：results/zero_routes_compare.json/.png
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zero_sim as zs  # noqa: E402

RES = os.path.join(HERE, '..', 'results')
TARGET = np.array([0.1724, 0.2759, 0.5517])  # 底/中/顶 观测粗分组
CK = []


def ck(name, cond, detail=''):
    CK.append({'name': name, 'pass': bool(cond), 'detail': detail})
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" —— {detail}" if detail else ''))


def height(w):
    s = 0; h = 0
    for x in w:
        s += x; h = max(h, abs(s))
    return h


def rho_A():
    """L=8 词类 (h,σ,π)，每类保留 top2 词（按 area, 字典序）。"""
    L = 8
    groups = {}
    for w in zs.balanced_words(L):
        c = (height(w) - 1, 0 if w[0] > 0 else 1, zs.word_area(w) % 2)
        groups.setdefault(c, []).append(w)
    cnt = np.zeros((4, 2, 2))
    for c, ws in groups.items():
        ws = sorted(ws, key=lambda w: (zs.word_area(w), w))
        for w in ws[:2]:  # top2 per class
            cnt[c] += 1
    rho = cnt / cnt.sum()
    return cnt, rho


def rho_B():
    """L=4 K4 记录，按 4/2/1 槽和观测粗分组给 profile。"""
    R = np.array([0.5517, 0.2759, 0.1034, 0.069])  # K4 顶点记录
    b, m, t = TARGET[0] / 4, TARGET[1] / 2, TARGET[2] / 1
    # 每个顶点：底 4、中 2、顶 1；权重按目标标定
    cnt = np.zeros((4, 2, 2))
    for v in range(4):
        cnt[v, :, :] = R[v] * b          # 底 4 槽
    rho = cnt / cnt.sum()
    # 中/顶用同样的类做 readout；这里返回底 rho，三层在 run 里按求和
    return cnt, rho


def run(rho, name, steps=200):
    H = rho.shape[0]
    bottom = rho.reshape(-1).copy(); bottom /= bottom.sum()
    middle = rho.sum(axis=2).reshape(-1).copy(); middle /= middle.sum()
    top = rho.sum(axis=(1, 2)).copy(); top /= top.sum()
    k, seed = 0.2, 0.05
    rows = []
    for gen in range(1, steps + 1):
        destroyed = gen % 4 == 0
        pre_top = top.copy()
        if destroyed:
            hstar = int(np.argmax(rho.sum(axis=(1, 2))))
            bottom[:] = 0.0
            mid = middle.reshape(H, 2)
            keep = np.zeros_like(mid); keep[hstar, :] = mid[hstar, :]
            middle = keep.reshape(-1)
        # flow 按 rho 目标注入
        j32 = k * top.sum()
        top -= top * (j32 / max(top.sum(), 1e-300))
        rho_mid = rho.sum(axis=2).reshape(-1); rho_mid = rho_mid / rho_mid.sum()
        middle += j32 * rho_mid
        j21 = k * middle.sum()
        middle -= middle * (j21 / max(middle.sum(), 1e-300))
        rho_bot = rho.reshape(-1); rho_bot = rho_bot / rho_bot.sum()
        bottom += j21 * rho_bot
        top += seed * rho.sum(axis=(1, 2))
        tot = bottom.sum() + middle.sum() + top.sum()
        bottom /= tot; middle /= tot; top /= tot
        rows.append({'gen': gen, 'destroyed': destroyed,
                     'shares': [float(bottom.sum()), float(middle.sum()), float(top.sum())],
                     'top_change': float(np.abs(top - pre_top).sum())})
    sh = np.array([r['shares'] for r in rows])
    mean = sh.mean(axis=0)
    dist = float(np.abs(mean - TARGET).sum())
    phases = [sh[p::4].mean(axis=0).tolist() for p in range(4)]
    ns = sh[:, 0]
    fft = np.abs(np.fft.rfft(ns - ns.mean()))
    period4 = float(fft[steps // 4]) if len(fft) > steps // 4 else 0.0
    return {'name': name, 'rows': rows, 'mean': mean.tolist(),
            'dist_target': dist, 'phases': phases, 'period4': period4}


def main():
    cntA, rhoA = rho_A()
    cntB, rhoB = rho_B()
    print('A rho shape', rhoA.shape, 'sum', rhoA.sum())
    print('B rho shape', rhoB.shape, 'sum', rhoB.sum())
    RA = run(rhoA, 'A_L8_word_top2')
    RB = run(rhoB, 'B_K4_fiber')
    for r in (RA, RB):
        print(f"{r['name']}: mean={np.round(r['mean'],4)} dist={r['dist_target']:.4f} "
              f"period4={r['period4']:.3f}")
        print('  phases', np.round(r['phases'],4).tolist())

    print('=' * 70); print('判据'); print('=' * 70)
    ck('A/B 都是 16 类 (4,2,2)', rhoA.size == 16 and rhoB.size == 16, '')
    ck('A/B rho 归一', abs(rhoA.sum()-1) < 1e-12 and abs(rhoB.sum()-1) < 1e-12, '')
    ck('两者都有周期 4 分量', RA['period4'] > 0 and RB['period4'] > 0, '')
    ck('未毁顶都在变', max(r['top_change'] for r in RA['rows']) > 1e-6
       and max(r['top_change'] for r in RB['rows']) > 1e-6, '')
    ck('A 更接近观测粗分组', RA['dist_target'] < RB['dist_target'],
       f"A={RA['dist_target']:.4f} B={RB['dist_target']:.4f}")

    g = np.arange(1, 201)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for r, lab in ((RA, 'A L8 word top2'), (RB, 'B K4 fiber')):
        sh = np.array([x['shares'] for x in r['rows']])
        ax[0].plot(g, sh[:, 0], label=f'{lab} bottom', lw=1)
        ax[1].plot(g, sh[:, 0], label=f'{lab} bottom', lw=1)
    ax[0].set_title('(a) bottom share'); ax[0].legend(fontsize=8)
    ax[1].set_title('(b) bottom share (zoom)'); ax[1].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    fig.tight_layout()
    fp = os.path.join(RES, 'zero_routes_compare.png')
    fig.savefig(fp, dpi=130)

    out = {'A': RA, 'B': RB, 'target': TARGET.tolist(),
           'checks': CK, 'n_fail': sum(1 for c in CK if not c['pass']),
           'anchor': 'A: L=8 word classes; B: K4 fiber profile; no dimensionful anchor'}
    with open(os.path.join(RES, 'zero_routes_compare.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print(f"合计 {len(CK)} 条：PASS {len(CK)-out['n_fail']} / FAIL {out['n_fail']}")
    print('写入', os.path.join(RES, 'zero_routes_compare.json'), '和', fp)


if __name__ == '__main__':
    main()
