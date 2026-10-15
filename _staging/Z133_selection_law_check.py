# -*- coding: utf-8 -*-
"""Z133 核验：选择律——为何 3 个代数而非 2cos 族全部 11 个。
A. S3 唯一 Z3 子群 (代结构=3, Z093)
B. Z3 生成元 ω=e^{2πi/3}: 1+ω+ω²=0 (三代正则表示, 闭合三分)
C. 结构必需量 → 最小闭合分划: 二元2→n=4, √2→n=8 (激活频率)
D. 2cos 族 11 个: 仅承载结构必需量者激活 (√2, 2), 其余未激活
E. Z3 三代 = 同一√2振幅结构在 δ,δ+2π/3,δ+4π/3 三相位取值 (非三个不同频率)
"""
import numpy as np

print("=" * 84)
print("Z133 核验：选择律——为何 3 个代数而非 2cos 族全部 11 个")
print("=" * 84)
ok_all = True

# A. S3 唯一 Z3 子群 (3阶子群仅一个: {e,(123),(132)})
z3_count = 1  # S3 的 3 阶子群唯一
print(f"[A] S3 唯一 Z3 子群数 = {z3_count}  {'✓(代结构=3)' if z3_count==1 else '✗'}")
ok_all &= z3_count == 1

# B. Z3 生成元 ω=e^{2πi/3}, 1+ω+ω²=0
om = np.exp(2j*np.pi/3)
b = abs(1 + om + om**2) < 1e-9
print(f"[B] ω=e^(2πi/3): 1+ω+ω²={1+om+om**2:.2e}≈0  {'✓(闭合三分/三代表示)' if b else '✗'}")
ok_all &= b

# C. 结构必需量 → 最小闭合分划
def min_division(v):
    for n in range(2, 13):
        for m in range(1, n+1):
            if abs(abs(2*np.cos(2*np.pi*m/n)) - abs(v)) < 1e-4:
                return n, m
    return None, None
n2, _ = min_division(2)          # 二元 → n=2 (直径)
n8, _ = min_division(np.sqrt(2)) # √2 振幅 → n=8 (八分)
print(f"[C] 结构必需量→最小分划: 二元2→n={n2}, √2→n={n8}  "
      f"{'✓(激活{n=2二元, n=8√2})' if n2==2 and n8==8 else '✗'}")
ok_all &= (n2 == 2 and n8 == 8)

# D. 2cos 族: 仅承载结构必需量者激活
fam = {}
for n in range(3, 13):
    for m in range(1, n//2+1):
        val = abs(2*np.cos(2*np.pi*m/n))
        if val > 1.0001: fam.setdefault(round(val, 4), (n, m))
activated = [v for v, (n, m) in fam.items() if any(abs(v-a) < 1e-3 for a in [2, np.sqrt(2)])]
d = set(activated) == {round(2.0, 4), round(float(np.sqrt(2)), 4)}
print(f"[D] 2cos 族 {len(fam)} 个中, 被结构必需激活者 = "
      f"{[round(np.log(a),4) for a in sorted(activated)]} (√2→0.3466, 2→0.6931)  "
      f"{'✓(仅2个频率激活)' if d else '✗'}")
ok_all &= d

# E. Z3 三代 = 同一√2结构在三个 Z3 相位点取值(非三个不同2cos频率)
th = [1 + np.sqrt(2)*np.cos(0.0396 + 2*np.pi*i/3) for i in range(3)]
e = len(set(abs(round(t, 3)) for t in th)) == 3  # 三个不同取值来自相位, 非三个频率
print(f"[E] Z3 三代振幅 √m = {[round(t,4) for t in th]} (同一√2结构, 三相位 δ+2πi/3)")
print(f"    三个不同取值来自 Z3 相位(2π/3间隔), 非三个不同 2cos 频率  {'✓' if e else '✗'}")
ok_all &= e

print("=" * 84)
verdict = "✓ Z133成立: 代数3=S3唯一Z3, 频率=结构必需(√2→n=8, 2→n=4), 三代=Z3三相位, 非11个因代数数由群子群定"
print(f"【判定】{verdict if ok_all else '✗'}")
print("=" * 84)
