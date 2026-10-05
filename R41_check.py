#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R41_check.py -- 核验维数无关性命题、数值、以及对 R38/R40 注脚的更正。"""
from __future__ import annotations
import io, json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); FAIL = 0
def head(t): print(""); print("="*72); print(t); print("="*72)
def check(n, c, d=""):
    global FAIL
    ok = bool(c)
    if not ok: FAIL += 1
    print("  [%s] %s%s" % ("v" if ok else "x", n, ("   "+d) if d else ""))
def read(n):
    with io.open(os.path.join(HERE, n), encoding="utf-8") as f: return f.read()

DOC = read("R41_dimension_independence.md")
R38 = read("R38_entanglement_from_shared_closure_origin.md")
R40 = read("R40_seeding_does_not_inject_site.md")
STATUS = read("STATUS.md"); INDEX = read("INDEX.md")
RES = json.load(io.open(os.path.join(HERE, "R41_dimension_independence_results.json"), encoding="utf-8"))

I2 = np.eye(2, dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex); sp = np.array([[0,0],[1,0]], dtype=complex)
sx = np.array([[0,1],[1,0]], dtype=complex); sy = np.array([[0,-1j],[1j,0]], dtype=complex)
P = [sx, sy, sz]

def snake(shape):
    idx = []
    if len(shape) == 1:
        idx = [(k,) for k in range(shape[0])]
    elif len(shape) == 2:
        for i in range(shape[0]):
            rng = range(shape[1]) if i % 2 == 0 else range(shape[1]-1, -1, -1)
            for j in rng: idx.append((i, j))
    else:
        for k in range(shape[0]):
            sub = snake(shape[1:]); sub = sub if k % 2 == 0 else sub[::-1]
            for s in sub: idx.append((k,)+s)
    return idx

def build(shape):
    path = snake(shape); N = len(path)
    def cd(j):
        out = np.array([[1.0+0j]])
        for k in range(N):
            out = np.kron(out, sz if k < j else (sp if k == j else I2))
        return out
    vac = np.zeros(2**N, dtype=complex); vac[0] = 1.0
    return path, N, cd, vac

def s_max(rho):
    T = np.array([[np.trace(rho @ np.kron(P[i], P[j])).real for j in range(3)] for i in range(3)])
    u = np.linalg.svd(T, compute_uv=False)
    return 2.0*float(np.sqrt(u[0]**2 + u[1]**2))

def reduce_two(psi, N, keep):
    T = psi.reshape([2]*N); others = [k for k in range(N) if k not in keep]
    rho = np.tensordot(T, T.conj(), axes=(others, others))
    perm = [keep.index(k) for k in sorted(keep)]
    rho = np.transpose(rho, perm + [p+2 for p in perm])
    return rho.reshape(4, 4)

head("F1  文档结构：命题、数值、更正、边界")
check("标题与性质为更正＋命题＋数值检验", "维数无关性" in DOC and "更正" in DOC and "命题" in DOC)
check("(R41-1)(R41-2)(R41-3) 在位",
      "(R41-1)" in DOC and "(R41-2)" in DOC and "(R41-3)" in DOC)
check("写明选维是另一件事", "另一件事" in DOC and "D=4" in DOC)
check("没有声称推出空间维数", "没有推出空间维数" in DOC)
check("R38 注脚已更正为维数无关", "维数无关（更正 2026-10-03" in R38 and "不是机制本身" in R38)
check("R40 注脚已更正为维数无关", "维数无关（更正 2026-10-03" in R40)

head("F2  独立复算：1D/2D/3D 一律 2√2")
for name, shape in (("1D_L4", (4,)), ("2D_2x2", (2,2)), ("2D_2x4", (2,4)), ("3D_2x2x2", (2,2,2))):
    path, N, cd, vac = build(shape)
    s0, s1 = np.array(path[0]), np.array(path[1])
    adj = bool(np.sum(np.abs(s0 - s1)) == 1)
    sup = (cd(0) + cd(1)) @ vac / math.sqrt(2); sup = sup/np.linalg.norm(sup)
    S = s_max(reduce_two(sup, N, [0, 1]))
    check("%s：位点格上相邻且 S=2√2" % name,
          adj and abs(S - 2*math.sqrt(2)) < 1e-9, "S=%.6f" % S)
check("JSON 判定维数无关",
      RES["verdict"]["dimension_independent"] is True
      and RES["verdict"]["all_dimensions_give_2sqrt2"] is True)

head("F3  账本一致")
check("STATUS 已登记 R41", "### 2.41 维数无关性（`R41`）当前状态" in STATUS
      and "R41_dimension_independence.md" in STATUS and "R41_check.py" in STATUS)
check("INDEX 已收录 R41", "R41_dimension_independence.md" in INDEX and "R41_check.py" in INDEX)
print(""); print("  不符项：%d" % FAIL)
if FAIL: print("  未通过"); sys.exit(1)
print("  全部通过"); sys.exit(0)
