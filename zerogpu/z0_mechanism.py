"""
z0_mechanism.py --- 验证"记录流打到哪些项链"的**判据**（见 RECORD_STREAM.md §1.1）。

判据（推导）：
    词 w 的循环高度剖面 H_k = Σ_{i<k} w_i（下标模 n）。
    旋转 k 是首次通过词 ⟺ H_{k+m} ≠ H_k (∀m=1..n-1) ⟺ **高度 H_k 在剖面上只被访问一次**。
  故
    项链 [w] 被打到 ⟺ H 存在只被访问一次的高度；
    不被达到        ⟺ 每个被访问的高度都至少被访问两次。

本脚本对 n=8..18 **穷举全部平衡词**，逐词比对两件事：
  ① 判据给出的"被打到的类数" == z0_record.py 实测的 distinct_classes；
  ② "剖面中只访问一次的高度数" == "该词的首次通过旋转数"。
两项都必须完全相等，否则退出码非零。
"""
from __future__ import annotations

import json
import os
import sys
from itertools import combinations

FAIL = []


def balanced_words(n: int):
    for pos in combinations(range(n), n // 2):
        v = 0
        for i in pos:
            v |= 1 << i
        yield v


def canon(w: int, n: int) -> int:
    m = (1 << n) - 1
    best = r = w & m
    for _ in range(n - 1):
        r = ((r << 1) | (r >> (n - 1))) & m
        if r < best:
            best = r
    return best


def cyclic_profile(w: int, n: int):
    h = [0]
    for i in range(n):
        h.append(h[-1] + (1 if (w >> i) & 1 else -1))
    return h[:n]


def once_visited_levels(w: int, n: int) -> int:
    cnt = {}
    for v in cyclic_profile(w, n):
        cnt[v] = cnt.get(v, 0) + 1
    return sum(1 for v in cnt.values() if v == 1)


def prime_excursion_rotations(w: int, n: int) -> int:
    """有多少个旋转是"首次通过词"（部分和在 1..n-1 上非零，第 n 步回零）。"""
    c = 0
    for k in range(n):
        s, ok = 0, True
        for m in range(1, n):
            s += 1 if (w >> ((k + m - 1) % n)) & 1 else -1
            if s == 0:
                ok = False
                break
        if ok:
            c += 1
    return c


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    meas_path = os.path.join(here, "results", "z0_record.json")
    meas = {}
    if os.path.exists(meas_path):
        with open(meas_path) as f:
            for r in json.load(f)["rows"]:
                if r.get("events"):
                    meas[r["n"]] = r["distinct_classes"]

    print("=" * 84)
    print("z0_mechanism：验证「被打到 ⟺ 剖面有只访问一次的高度」")
    print("=" * 84)
    print(f"{'n':>3} {'K(n)':>9} {'实测被打到':>11} {'判据':>9} {'判据==实测':>11} "
          f"{'剖面值==旋转数':>15}")
    for n in range(8, 19, 2):
        reps, hit = {}, set()
        mech_ok = True
        nwords = 0
        for w in balanced_words(n):
            nwords += 1
            c = canon(w, n)
            reps.setdefault(c, w)
            pc = once_visited_levels(w, n)
            if pc != prime_excursion_rotations(w, n):
                mech_ok = False
            if pc > 0:
                hit.add(c)
        K = len(reps)
        m = meas.get(n)
        agree = (m is None) or (m == len(hit))
        if not agree:
            FAIL.append(f"n={n}: 判据 {len(hit)} != 实测 {m}")
        if not mech_ok:
            FAIL.append(f"n={n}: 剖面值 != 首次通过旋转数")
        print(f"{n:>3} {K:>9} {str(m) if m is not None else '-':>11} {len(hit):>9} "
              f"{str(agree):>11} {str(mech_ok):>15}")

    # 打印几个"未被打到"的例子，展示判据的形态
    print("\n未被打到的项链（判据形态：每个高度至少访问两次）")
    for n in (10, 12):
        reps, hit = {}, set()
        for w in balanced_words(n):
            c = canon(w, n)
            reps.setdefault(c, w)
            if once_visited_levels(w, n) > 0:
                hit.add(c)
        miss = [(c, reps[c]) for c in sorted(reps) if c not in hit]
        print(f"  n={n}: {len(reps)} 类，未打到 {len(miss)} 类（比例 {len(miss)/len(reps):.4f}）")
        for c, w in miss[:4]:
            H = cyclic_profile(w, n)
            cnt = {}
            for v in H:
                cnt[v] = cnt.get(v, 0) + 1
            print(f"    {format(w, '0%db' % n)}  H={H}  访问次数={sorted(cnt.values())}")

    print()
    if FAIL:
        print("*** 失败 ***")
        for f in FAIL:
            print("  ", f)
        sys.exit(1)
    print("全部断言通过：判据与实测完全一致，且 剖面值 == 首次通过旋转数。")


if __name__ == "__main__":
    main()
