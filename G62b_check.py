#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G62b_check.py -- 独立复核入口（不改动 G62_check.py、不重复其代码）。

本脚本只做一件事：以子进程跑同目录的 G62_check.py，转发其数值核验结论。
G62_check.py 已在 2026-02 的 Gleason 修正中更新：
  F5  【已撤回阈值主张】dim A = 4(T+1)、dim H_T = 2(T+1)、正交极小投影数 = 2(T+1)
      三个对象分别核验；G62 §4 的反例在 L(A_T) 上对每个 T 都成立。
  F6  L(M_2) 的唯一正交对是 {P,1-P} => 加性退化为互补加性；GNS 向量态直接给迹形式。

跑法：
    python3 G62b_check.py            # 转发 G62_check.py 的全部输出
需要：numpy；同目录的 G62_check.py；不读任何 .md。
耗时：约 0.7 秒；通过项：与 G62_check.py 相同（31 项）。
"""

import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
NESTED = os.path.join(HERE, "G62_check.py")

if not os.path.exists(NESTED):
    print("  缺少 %s（本脚本只转发它的核验）" % NESTED)
    sys.exit(2)

r = subprocess.run([sys.executable, NESTED], cwd=HERE, capture_output=True, text=True)
sys.stdout.write(r.stdout)
sys.stderr.write(r.stderr)
n = r.stdout.count("[v]")
print("  G62b_check.py：转发 G62_check.py，通过 %d 项，退出码 %d" % (n, r.returncode))
sys.exit(0 if r.returncode == 0 else 1)
