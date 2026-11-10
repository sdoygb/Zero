#!/usr/bin/env python3
"""Zero 系综指纹器：PDB（多模型 NMR 系综）→ σ（接触三档）、|R|（跨模型接触开关数）。

用法:
  python3 zerofp_pdb.py 1D3Z.pdb 2MJB.pdb ...      # 打印表格
  python3 zerofp_pdb.py --selftest                  # 自检（断言全部通过才返回 0）

判据参考 _staging/推导/R历史判别_NMR系综检验.md：
  - 接触：CA ≤ 8 Å，序列间隔 ≥ 2
  - |R| = 相邻模型间接触开关的总次数
  - sigma = 全部模型接触按序差分三档（短 2-5 / 中 6-12 / 长 ≥13）的占比
"""
import sys, os, math
from collections import defaultdict

CUT = 8.0

def parse_models(path):
    """返回 [ [(x,y,z)...] ... ]，每个模型一个残基序列（按 (chain, resseq) 排序）。"""
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    models, cur = [], {}
    for line in open(path, errors="ignore"):
        if line.startswith("MODEL"):
            cur = {}
        elif line.startswith("ATOM") and line[12:16].strip() == "CA" and line[16] in (" ", "A"):
            cur[(line[21], line[22:27])] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        elif line.startswith("ENDMDL"):
            if cur: models.append([cur[k] for k in sorted(cur)])
    if not models and cur:            # 单模型 PDB（无 MODEL 记录）
        models.append([cur[k] for k in sorted(cur)])
    return models

def contacts(pts):
    n = len(pts)
    if n < 3: return set()
    D = [[0.0]*n for _ in range(n)]
    for i in range(n):
        xi, yi, zi = pts[i]
        for j in range(i+2, n):
            xj, yj, zj = pts[j]
            dx, dy, dz = xi-xj, yi-yj, zi-zj
            if dx*dx + dy*dy + dz*dz <= CUT*CUT:
                D[i][j] = 1
    return {(i, j) for i in range(n) for j in range(i+2, n) if D[i][j]}

def fingerprint(path):
    ms = parse_models(path)
    if len(ms) < 2:
        return {"path": path, "models": len(ms), "error": "insufficient models (need >=2)"}
    sets = [contacts(p) for p in ms]
    allc = set().union(*sets)
    trans = sum(sum(1 for k in range(1, len(sets)) if ((c in sets[k]) != (c in sets[k-1]))) for c in allc)
    tot = defaultdict(int)
    for s in sets:
        for (i, j) in s:
            d = j - i
            tot["short" if d <= 5 else ("med" if d <= 12 else "long")] += 1
    N = sum(tot.values()) or 1
    return {"path": path, "models": len(ms), "residues": len(ms[0]), "union": len(allc),
            "R": trans, "R_norm": trans/(len(ms)-1),
            "sigma": (tot["short"]/N, tot["med"]/N, tot["long"]/N)}

def selftest(pdbdir=None):
    d = pdbdir or os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "_staging", "推导", "pdb")
    ok = True
    def check(cond, msg):
        nonlocal ok
        print(("PASS " if cond else "FAIL ") + msg)
        ok = ok and cond
    # 1) 多模型：五个泛素系综
    fam = {}
    for pid in ("1D3Z", "2MJB", "1G6J", "2NR2", "2K39"):
        p = os.path.join(d, pid + ".pdb")
        if not os.path.exists(p):
            check(False, f"{pid} 存在"); continue
        r = fingerprint(p); fam[pid] = r
        check("error" not in r and r["models"] >= 2, f"{pid} 多模型解析 (models={r.get('models')})")
        check(abs(sum(r["sigma"]) - 1.0) < 1e-9, f"{pid} sigma 归一 (sum={sum(r['sigma']):.6f})")
        check(r["R_norm"] > 0, f"{pid} |R|>0 (R_norm={r['R_norm']:.2f})")
    # 2) 同一蛋白 sigma 稳定、|R| 变化大
    if len(fam) == 5:
        shorts = [r["sigma"][0] for r in fam.values()]
        longs = [r["sigma"][2] for r in fam.values()]
        check(max(shorts)-min(shorts) < 0.05, f"sigma_short 稳定 (范围 {max(shorts)-min(shorts):.3f})")
        check(max(longs)-min(longs) < 0.05, f"sigma_long 稳定 (范围 {max(longs)-min(longs):.3f})")
        rr = [r["R_norm"] for r in fam.values()]
        check(max(rr)/min(rr) > 5, f"|R| 区分度 (max/min={max(rr)/min(rr):.1f})")
    # 3) 单模型文件 -> insufficient
    p = os.path.join(d, "4AKE.pdb")
    if os.path.exists(p):
        r = fingerprint(p)
        check("error" in r and r["models"] < 2, f"单模型 4AKE 判 insufficient (models={r['models']})")
    # 4) 缺文件 -> 抛 FileNotFoundError
    try:
        fingerprint(os.path.join(d, "__no_such__.pdb")); check(False, "缺文件应抛错")
    except FileNotFoundError:
        check(True, "缺文件抛 FileNotFoundError")
    # 5) 空/坏文件 -> 无模型
    tmp = os.path.join(d, "_empty_test.pdb"); open(tmp, "w").close()
    r = fingerprint(tmp)
    check("error" in r, "空文件判 insufficient")
    os.remove(tmp)
    print("\nSELFTEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1

def main():
    args = sys.argv[1:]
    if not args or args[0] == "--selftest":
        return selftest()
    print(f"{'PDB':10s} {'models':>6s} {'res':>4s} {'union':>6s} {'|R|':>7s} {'|R|/(m-1)':>10s}  sigma(短,中,长)")
    for p in args:
        try:
            r = fingerprint(p)
        except FileNotFoundError:
            print(f"{p:10s} ERROR file not found"); continue
        if "error" in r:
            print(f"{os.path.basename(p):10s} {r['models']:6d}  {r['error']}"); continue
        s = r["sigma"]
        print(f"{os.path.basename(p):10s} {r['models']:6d} {r['residues']:4d} {r['union']:6d} {r['R']:7d} {r['R_norm']:10.2f}  ({s[0]:.3f},{s[1]:.3f},{s[2]:.3f})")
    return 0

if __name__ == "__main__":
    sys.exit(main())
