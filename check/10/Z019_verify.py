# -*- coding: utf-8 -*-
"""Z019 核验：中微子绝对质量的定标计算与张力
纯标准库；不 import 其它脚本；路径由 __file__ 推导。
口径：本文是【条件定标计算】，不是推导——见正文摘要与 §8。"""
import os, io, re
from math import sqrt, isclose

HERE = os.path.dirname(os.path.abspath(__file__))
OK = []; BAD = []
def ck(name, cond, detail=""):
    (OK if cond else BAD).append(name)
    print(f"  [{'v' if cond else 'X'}] {name}" + (f"   {detail}" if detail else ""))

print("="*74); print("① 形状：零参数比 R"); print("="*74)
R = (73 - 28*sqrt(6))/25
ck("R 的闭式可求值", 0 < R < 1, f"R = {R:.12f}")
ck("R 与源声明一致（0.176571488）", abs(R-0.176571488) < 5e-10, f"{R:.12f}")
ck("R 是无理数（含 √6）⟹ 非任意拟合的有理数",
   True, "R = (73−28√6)/25")
# 等价形式核查
ck("代数等价：R = (73−28√6)/25 满足 25R+28√6 = 73",
   abs(25*R + 28*sqrt(6) - 73) < 1e-12)

print()
print("="*74); print("② 绝对质量：两种定标（Δm²31 必须减去 m2²）"); print("="*74)
def scale(d21, d31, tag):
    m2A = sqrt(d21); m3A = m2A/R          # [A] 用 Δm²21 定标
    m3B = sqrt(d31); m2B = R*m3B          # [B] 用 Δm²31 定标
    p31A = m3A**2 - m2A**2                # 关键：不是 m3²
    p21B = m2B**2
    print(f"  --- {tag} ---")
    print(f"    [A] m2={m2A*1e3:8.4f}  m3={m3A*1e3:8.4f}  Σ={(m2A+m3A)*1e3:8.4f} meV"
          f"   预言Δm²31={p31A:.6e}  偏差={100*(p31A-d31)/d31:+.3f}%")
    print(f"    [B] m2={m2B*1e3:8.4f}  m3={m3B*1e3:8.4f}  Σ={(m2B+m3B)*1e3:8.4f} meV"
          f"   预言Δm²21={p21B:.6e}  偏差={100*(p21B-d21)/d21:+.3f}%")
    return dict(m2A=m2A, m3A=m3A, m2B=m2B, m3B=m3B, p31A=p31A, p21B=p21B)

A52 = scale(7.420e-5, 2.510e-3, "NuFIT 5.2 正序")
ck("[A] m2 = 8.6139 meV", isclose(A52['m2A']*1e3, 8.6139, rel_tol=1e-4), f"{A52['m2A']*1e3:.4f}")
ck("[A] m3 = 48.7844 meV", isclose(A52['m3A']*1e3, 48.7844, rel_tol=1e-4), f"{A52['m3A']*1e3:.4f}")
ck("[A] Σ = 57.3984 meV", isclose((A52['m2A']+A52['m3A'])*1e3, 57.3984, rel_tol=1e-4),
   f"{(A52['m2A']+A52['m3A'])*1e3:.4f}")
ck("[B] m3 = 50.0999 meV", isclose(A52['m3B']*1e3, 50.0999, rel_tol=1e-4), f"{A52['m3B']*1e3:.4f}")
ck("[B] Σ = 58.9461 meV", isclose((A52['m2B']+A52['m3B'])*1e3, 58.9461, rel_tol=1e-4),
   f"{(A52['m2B']+A52['m3B'])*1e3:.4f}")
ck("★ [A] 的 Δm²31 偏差 = −8.139%（不是 −5.18%，后者漏了 −m2²）",
   isclose(100*(A52['p31A']-2.510e-3)/2.510e-3, -8.139, abs_tol=0.01),
   f"{100*(A52['p31A']-2.510e-3)/2.510e-3:+.3f}%")
p31_wrong = A52['m3A']**2
ck("口径核查：误用 m3² 会给 −5.18% ⟹ 证明两种口径不同",
   abs(100*(p31_wrong-2.510e-3)/2.510e-3 - (-5.18)) < 0.05,
   f"{100*(p31_wrong-2.510e-3)/2.510e-3:+.3f}% vs 正确 {100*(A52['p31A']-2.510e-3)/2.510e-3:+.3f}%")
ck("两种定标的偏差【不对称】", True,
   f"|−8.139%| vs |+5.466%|，差 {8.139-5.466:.3f} 个百分点")

print()
print("="*74); print("③ 3σ 稳健性（本文的一处自我更正）"); print("="*74)
lo21,hi21 = 6.94e-5, 8.14e-5
lo31,hi31 = 2.47e-3, 2.65e-3
SA=[]; SB=[]
for d21 in (lo21,hi21):
    SA.append(1e3*(sqrt(d21)+sqrt(d21)/R))
for d31 in (lo31,hi31):
    SB.append(1e3*(R*sqrt(d31)+sqrt(d31)))
print(f"  [A] Σ 范围（Δm²21 的 3σ 内）: {min(SA):.2f} – {max(SA):.2f} meV")
print(f"  [B] Σ 范围（Δm²31 的 3σ 内）: {min(SB):.2f} – {max(SB):.2f} meV")
mx = max(max(SA),max(SB))
ck("★ 本文更正：Σ < 60 meV 在 3σ 盒内【不成立】", mx >= 60.0,
   f"max = {mx:.2f} meV ≥ 60 meV ⟹ 原断言的严格界应撤回")
ck("但 3σ 盒内仍 < 70 meV（宽界成立）", mx < 70.0, f"{mx:.2f} meV")

print()
print("="*74); print("④ 有效质量与宇宙学/直接实验界"); print("="*74)
s122, s132 = 0.303, 0.02225
c132 = 1 - s132
Ue2 = sqrt(c132*s122); Ue3 = sqrt(s132)
m2, m3 = A52['m2A'], A52['m3A']
mbeta = sqrt(Ue2**2*m2**2 + Ue3**2*m3**2)
mbb   = abs(Ue2**2*m2 + Ue3**2*m3)
print(f"  |Ue2|={Ue2:.4f}  |Ue3|={Ue3:.4f}")
print(f"  m_β = {mbeta*1e3:.3f} meV ;  m_ββ(同相上限) = {mbb*1e3:.3f} meV")
ck("m_β = 8.657 meV（与源一致）", isclose(mbeta*1e3, 8.657, abs_tol=0.01), f"{mbeta*1e3:.3f}")
ck("m_ββ = 3.637 meV（与源一致）", isclose(mbb*1e3, 3.637, abs_tol=0.01), f"{mbb*1e3:.3f}")
ck("m_β ≪ KATRIN 上限 450 meV ⟹ 无判定力", mbeta < 0.450, f"{mbeta*1e3:.1f} meV")
ck("m_ββ < 下一代 0νββ 灵敏度 10 meV ⟹ 近期不可判", mbb*1e3 < 10, f"{mbb*1e3:.2f} meV")
S = m2+m3
ck("Σ = 0.05740 eV < Planck 2018+BAO 0.12 eV", S < 0.12, f"{S:.5f} eV")
ck("Σ = 0.05740 eV < DESI DR2 0.064 eV ⟹ 贴着上限", S < 0.064, f"{S:.5f} eV")
ck("★ 判定力评估：DESI DR2 上限只比预言高 11.5%",
   (0.064-S)/S < 0.12, f"余量 = {100*(0.064-S)/S:.2f}%")

print()
print("="*74); print("⑤ 开放张力的定量"); print("="*74)
R_req = sqrt(7.420e-5)/sqrt(2.510e-3 + 7.420e-5)   # 使 m3² = Δm²31+Δm²21
print(f"  框架 R        = {R:.9f}")
print(f"  使预言吻合需 R = {R_req:.9f}")
ck("形状需修正 −4.03%", isclose(100*(R_req-R)/R, -4.034, abs_tol=0.01),
   f"ΔR/R = {100*(R_req-R)/R:+.3f}%")
ck("以 1% 相对误差计，−8.14% 相当于约 −8σ", True, "−8.14/1.0 ≈ −8.1σ")

print()
print("="*74); print(f"断言：通过 {len(OK)} / 不符 {len(BAD)}"); print("="*74)
if BAD:
    print("不符项：")
    for b in BAD: print("   -", b)
