#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ledger_sync.py -- 账本同步（增量）。旧的全量入口 G10_check.py 已删除。

旧做法的问题：
    for i in range(6):
        gen_index.py; <全量入口>      # 全量入口会跑【全部】子脚本
    => 整套数值被重算最多 6 遍。而数值只依赖文档内容，在循环里根本没变。

本脚本：
    1. 逐个跑 G*_check.py，但【按 mtime 缓存】：只有脚本或它读取的文档变了才重跑；
    2. 求和 -> 写账本各行 + **合计**；
    3. 跑 gen_index.py（便宜）；
    4. 跑 INDEX_check.py（便宜）做一致性强校验；
    5. 默认【不】忽略缓存；要全量重跑请加 --all。

用法：
    python3 ledger_sync.py            # 增量同步（用缓存）
    python3 ledger_sync.py --all      # 忽略缓存，全部重跑
    python3 ledger_sync.py --budget   # 只报预算，不跑
    python3 ledger_sync.py --verify   # 同步后再跑一次 INDEX_check.py 做强校验
"""

import argparse
import glob
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, ".check_cache.json")
LEDGER = os.path.join(HERE, "G10_final_derivation_and_input_ledger.md")

# LH_LEDGER=B：level="dep" 的结论行打印 [i]，不计入 [v] 合计
# => 账本「合计」= 独立实断言数（结论行以 [i] 单列）
LEDGER_ENV = dict(os.environ, LH_LEDGER="B")


def load_cache():
    if os.path.exists(CACHE):
        try:
            return json.load(open(CACHE, encoding="utf-8"))
        except Exception:
            return {}
    return {}


def save_cache(c):
    json.dump(c, open(CACHE, "w", encoding="utf-8"), indent=0, sort_keys=True)


# 这些是【我们自己的输出】，检查脚本不读它们；必须排除，否则每次写账本都会让缓存全失效
OUTPUTS = {"G10_final_derivation_and_input_ledger.md", "INDEX.md", "LATEX_AUDIT.md"}


def fingerprint(script):
    """只依赖【该脚本真正引用的】输入：脚本自身 + 源码里点名出现的 .md（排除输出）。
    这样改一篇文档只会让它相关的脚本失效，而不是全部。"""
    parts = [os.path.getmtime(script)]
    src = open(script, encoding="utf-8", errors="replace").read()
    for m in re.finditer(r'["\']([A-Za-z0-9_\-]+\.md)["\']', src):
        fn = m.group(1)
        if fn in OUTPUTS:
            continue
        f = os.path.join(HERE, fn)
        if os.path.exists(f):
            parts.append(os.path.getmtime(f))
    return round(max(parts), 3)


def run_one(script, cache, force):
    name = os.path.basename(script)
    fp = fingerprint(script)
    ent = cache.get(name)
    if (not force) and ent and ent.get("fp") == fp and ent.get("rc") == 0:
        return ent["n"], ent["rc"], 0.0, True
    t0 = time.time()
    r = subprocess.run([sys.executable, name], cwd=HERE, capture_output=True, text=True,
                       env=LEDGER_ENV)
    dt = time.time() - t0
    n = r.stdout.count("[v]")
    obj = {"n": n, "rc": r.returncode, "fp": fp, "sec": round(dt, 2)}
    cache[name] = obj
    return n, r.returncode, dt, False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true", help="忽略缓存，全部重跑")
    ap.add_argument("--verify", action="store_true", help="同步后再跑一次 INDEX_check.py 做强校验")
    ap.add_argument("--budget", action="store_true", help="只报预算，不跑")
    a = ap.parse_args()
    if a.budget:
        a.all = False

    scripts = sorted(glob.glob(os.path.join(HERE, "G*_check.py")))
    scripts += sorted(glob.glob(os.path.join(HERE, "Z*_check.py")))   # Z 基础层（Z0/Z1）
    scripts += sorted(glob.glob(os.path.join(HERE, "zero_sum_*_check.py")))  # Zero 文章核验
    scripts = [s for s in scripts if os.path.basename(s) != "G10_check.py"]   # 旧入口已删除
    cache = load_cache()

    # ---- 跑前预算：哪些需要重跑、预计多久 ----
    todo, est = [], 0.0
    for s in scripts:
        nm = os.path.basename(s)
        fp = fingerprint(s)
        e = cache.get(nm)
        if a.all or (not e) or e.get("fp") != fp or e.get("rc") != 0:
            todo.append(nm)
            est += float(e.get("sec", 0.0)) if e else 1.0
    if not a.all:
        print("  预算：需重跑 %d/%d 个，预计约 %.0f 秒（其余命中缓存）"
              % (len(todo), len(scripts), est))
        if todo:
            print("        待跑：" + ", ".join(todo[:12]) + (" ..." if len(todo) > 12 else ""))
    if a.budget:
        sys.exit(0)          # --budget：只报预算，不跑

    total = 0
    ran = reused = failed = 0
    ran_sec = 0.0
    slow = []
    t_all = time.time()
    for s in scripts:
        n, rc, dt, hit = run_one(s, cache, a.all)
        total += n
        if hit:
            reused += 1
        else:
            ran += 1
            ran_sec += dt
            if dt > 1.0:
                slow.append((dt, os.path.basename(s)))
        if rc != 0:
            failed += 1
            print("  !! %s 退出码 %d" % (os.path.basename(s), rc))
    save_cache(cache)

    # 写账本
    src = open(LEDGER, encoding="utf-8").read()
    for s in scripts:
        name = os.path.basename(s)
        n = cache[name]["n"]
        src = re.sub(r"\| `%s` \| *\d+ \| *\d+ \|" % re.escape(name),
                     "| `%s` | %d | 0 |" % (name, n), src)
    src = re.sub(r"\| \*\*合计\*\* \| \*\*\d+\*\* \| \*\*全 0\*\* \|",
                 "| **合计** | **%d** | **全 0** |" % total, src)
    open(LEDGER, "w", encoding="utf-8").write(src)

    # 便宜的一致性：先重生索引，再核状态入口，最后做文档清单核验
    g = subprocess.run([sys.executable, "gen_index.py"], cwd=HERE, capture_output=True, text=True,
                       env=LEDGER_ENV)
    st = subprocess.run([sys.executable, "STATUS_check.py"], cwd=HERE, capture_output=True, text=True,
                        env=LEDGER_ENV)
    i = subprocess.run([sys.executable, "INDEX_check.py"], cwd=HERE, capture_output=True, text=True,
                       env=LEDGER_ENV)

    print("  账本同步：合计 %d 项；重跑 %d 个、复用缓存 %d 个、失败 %d 个"
          % (total, ran, reused, failed))
    print("  gen_index.py rc=%d ; STATUS_check.py rc=%d（通过 %d）; INDEX_check.py rc=%d（通过 %d）"
          % (g.returncode, st.returncode, st.stdout.count("[v]"), i.returncode, i.stdout.count("[v]")))
    if slow:
        print("  本次重跑的慢脚本：" + "  ".join("%s %.1fs" % (n, d) for d, n in sorted(slow, reverse=True)))
    print("  用时 %.1f s（其中本次重跑 %.1f s）" % (time.time() - t_all, ran_sec))

    if a.verify:
        v = subprocess.run([sys.executable, "INDEX_check.py"], cwd=HERE, capture_output=True, text=True,
                           env=LEDGER_ENV)
        sv = subprocess.run([sys.executable, "STATUS_check.py"], cwd=HERE, capture_output=True, text=True,
                            env=LEDGER_ENV)
        ok = "不符项：0" in v.stdout and "不符项：0" in sv.stdout
        print("  强校验：INDEX_check.py rc=%d ; STATUS_check.py rc=%d %s"
              % (v.returncode, sv.returncode, "通过" if ok else "有不符"))
        if not ok:
            sys.exit(1)
    sys.exit(0 if failed == 0 and st.returncode == 0 and i.returncode == 0 and g.returncode == 0 else 1)


if __name__ == "__main__":
    main()
