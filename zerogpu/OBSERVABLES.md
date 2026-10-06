# 按结构分类的可观测量目录（任务 1 产出）

> 生成时间 2026-10-06 09:36:18｜133 条可观测量｜机器 [`observable_sweep.py`](observable_sweep.py) ＋ 分类器 [`structures.py`](structures.py)｜数据 [`results/observable_sweep.json`](results/observable_sweep.json)

> **等级：全部是【数值证据】，不是【导出】。** 本页只贴**结构标签**；物理名字在 §5，**后验**给出且与结构分开。


---

## §1 生成规则（为什么这是「中性」的）

此前我一直在**挑**可观测量去对 QFT（找极点、找谱密度），那不是中性做法。

这里用机械的叉乘枚举：

```
可观测量 = (对象 O) × (算子/构造 K) × (读出 R)
O ∈ {Γ 上的行走占据 n_k(v)；𝒢_T 闭链图；词层 {±1}^n（无 Γ 对照）；缺陷算符；记录层账本}
K ∈ {计数、邻接 A、热核 e^{tA}、谱、商/投影、扰动 A+V|0><0|}
R ∈ {点值、谱、分布、熵、参与比、极值、支撑、矩、能隙、关联、标度}
```

分类器 `structures.classify_curve` 的规则：**所有候选模型都在对数空间算 R²**（衰减曲线跨若干数量级），取最高者，但若更简单的模型 R² 与之相差 < 2e-3 则取更简单的。复杂度次序 `constant < logarithmic < power-law < exponential < stretched`。**所有候选连同 R² 都写进 JSON**，便于复核与改判。


## §2 结构聚类

| 结构标签 | 条数 | 涉及扇区 |
|:--|--:|:--|
| `power-law` | 34 | CLOSURE, DEFECT, WALK, WORD |
| `stretched` | 26 | CLOSURE, DEFECT, LEDGER |
| `gapped` | 14 | CLOSURE |
| `exponential` | 12 | CLOSURE, DEFECT, LEDGER, WORD |
| `constant` | 10 | CLOSURE, DEFECT, LEDGER, WORD |
| `declared` | 10 | DEFECT, WALK |
| `logarithmic` | 10 | LEDGER, WALK, WORD |
| `too-short` | 6 | CLOSURE |
| `continuum/smooth` | 4 | DEFECT |
| `continuum/center-divergent(van-Hove)` | 3 | CLOSURE |
| `continuum/edge-vanishing+center-divergent(van-Hove)` | 3 | CLOSURE |
| `delta` | 1 | CLOSURE |

---

## §3 逐结构目录


### `power-law`（34 条）

| 可观测量 | 扇区 | 拟合 | R² | 备注 |
|:--|:--|:--|--:|:--|
| `CLOSURE:𝒢_T 系综:平均度 vs T` | CLOSURE | b = +1.0649 | 0.99495 |  |
| `DEFECT:Z^1 + V|0><0|:束缚态径向衰减 |ψ(r)| (V=-4.0)` | DEFECT | b = -10.3550 | 0.85828 | D=1；L=201；直线段 → 指数衰减 → 局域长度 ξ |
| `DEFECT:Z^1 + V|0><0|:束缚态径向衰减 |ψ(r)| (V=4.0)` | DEFECT | b = -9.6031 | 0.82432 | D=1；L=201；直线段 → 指数衰减 → 局域长度 ξ |
| `DEFECT:Z^1 + V|0><0|:束缚态能移 ΔE V>0（上带边 +2D）` | DEFECT | b = +1.6102 | 0.99199 | D=1；L=201；ΔE<0 且 |ΔE| 单调增 = 束缚态；标度见结构分类 |
| `DEFECT:Z^1 L=201:极值态 IPR vs V` | DEFECT | b = +0.6857 | 0.93879 | D=1；L=201 |
| `DEFECT:Z^1 L=401:极值态 IPR vs V` | DEFECT | b = +0.7660 | 0.96191 | D=1；L=401 |
| `DEFECT:Z^1 L=801:极值态 IPR vs V` | DEFECT | b = +0.8324 | 0.97648 | D=1；L=801 |
| `DEFECT:Z^2 + V|0><0|:束缚态能移 ΔE V>0（上带边 +2D）` | DEFECT | b = +3.9619 | 0.94884 | D=2；L=101；ΔE<0 且 |ΔE| 单调增 = 束缚态；标度见结构分类 |
| `DEFECT:Z^3 + V|0><0|:束缚态径向衰减 |ψ(r)| (V=1.0)` | DEFECT | b = -0.0325 | 0.85161 | D=3；L=31；直线段 → 指数衰减 → 局域长度 ξ |
| `DEFECT:Z^3 + V|0><0|:束缚态径向衰减 |ψ(r)| (V=4.0)` | DEFECT | b = -1.0175 | 0.99425 | D=3；L=31；直线段 → 指数衰减 → 局域长度 ξ |
| `DEFECT:Z^4 + V|0><0|:束缚态径向衰减 |ψ(r)| (V=1.0)` | DEFECT | b = -0.0207 | 0.88243 | D=4；L=13；直线段 → 指数衰减 → 局域长度 ξ |
| `DEFECT:Z^4 + V|0><0|:束缚态径向衰减 |ψ(r)| (V=4.0)` | DEFECT | b = -0.1626 | 0.89511 | D=4；L=13；直线段 → 指数衰减 → 局域长度 ξ |
| `WALK:Z^1 行走占据 n_k(v):二阶矩 M2=Σr²p` | WALK | b = +1.0000 | 1.00000 | D=1 |
| `WALK:Z^1 行走占据 n_k(v):参与比 IPR=Σn²` | WALK | b = -0.4905 | 0.99984 | D=1 |
| `WALK:Z^1 行走占据 n_k(v):可分辨支撑 #{p>1e-14}` | WALK | b = +0.8497 | 0.99459 | D=1 |
| `WALK:Z^1 行走占据 n_k(v):回返 R(k)=n_k(0)` | WALK | b = -0.4812 | 0.99939 | D=1 |
| `WALK:Z^1 行走占据 n_k(v):峰位概率 max_v p` | WALK | b = -0.4812 | 0.99939 | D=1 |
| `WALK:Z^2 行走占据 n_k(v):二阶矩 M2=Σr²p` | WALK | b = +1.0000 | 1.00000 | D=2 |
| `WALK:Z^2 行走占据 n_k(v):参与比 IPR=Σn²` | WALK | b = -0.9762 | 0.99982 | D=2 |
| `WALK:Z^2 行走占据 n_k(v):可分辨支撑 #{p>1e-14}` | WALK | b = +1.6835 | 0.99614 | D=2 |
| `WALK:Z^2 行走占据 n_k(v):回返 R(k)=n_k(0)` | WALK | b = -0.9531 | 0.99928 | D=2 |
| `WALK:Z^2 行走占据 n_k(v):峰位概率 max_v p` | WALK | b = -0.9531 | 0.99928 | D=2 |
| `WALK:Z^3 行走占据 n_k(v):二阶矩 M2=Σr²p` | WALK | b = +1.0000 | 1.00000 | D=3 |
| `WALK:Z^3 行走占据 n_k(v):参与比 IPR=Σn²` | WALK | b = -1.4604 | 0.99984 | D=3 |
| `WALK:Z^3 行走占据 n_k(v):可分辨支撑 #{p>1e-14}` | WALK | b = +2.5357 | 0.99735 | D=3 |
| `WALK:Z^3 行走占据 n_k(v):回返 R(k)=n_k(0)` | WALK | b = -1.4233 | 0.99935 | D=3 |
| `WALK:Z^3 行走占据 n_k(v):峰位概率 max_v p` | WALK | b = -1.4233 | 0.99935 | D=3 |
| `WALK:Z^4 行走占据 n_k(v):二阶矩 M2=Σr²p` | WALK | b = +1.0000 | 1.00000 | D=4 |
| `WALK:Z^4 行走占据 n_k(v):参与比 IPR=Σn²` | WALK | b = -1.9312 | 0.99986 | D=4 |
| `WALK:Z^4 行走占据 n_k(v):可分辨支撑 #{p>1e-14}` | WALK | b = +3.3884 | 0.99676 | D=4 |
| `WALK:Z^4 行走占据 n_k(v):回返 R(k)=n_k(0)` | WALK | b = -1.8491 | 0.99869 | D=4 |
| `WALK:Z^4 行走占据 n_k(v):峰位概率 max_v p` | WALK | b = -1.8491 | 0.99869 | D=4 |
| `WORD:1D 零和游走:回返 R(2m)=C(2m,m)/4^m` | WORD | b = -0.4825 | 0.99944 |  |
| `WORD:{±1}^n 全空间:平衡词占比 C(n,n/2)/2^n` | WORD | b = -0.4697 | 0.99914 |  |

### `stretched`（26 条）

| 可观测量 | 扇区 | 拟合 | R² | 备注 |
|:--|:--|:--|--:|:--|
| `CLOSURE:𝒢_16:热核 Z(t)=⟨e^{tλ}⟩` | CLOSURE | β = 1.25, c = -5.8104 | 0.99847 | T=16 |
| `CLOSURE:𝒢_18:热核 Z(t)=⟨e^{tλ}⟩` | CLOSURE | β = 1.25, c = -6.6215 | 0.99855 | T=18 |
| `CLOSURE:𝒢_20:热核 Z(t)=⟨e^{tλ}⟩` | CLOSURE | β = 1.25, c = -7.4394 | 0.99844 | T=20 |
| `CLOSURE:𝒢_22:热核 Z(t)=⟨e^{tλ}⟩` | CLOSURE | β = 1.25, c = -8.2638 | 0.99829 | T=22 |
| `CLOSURE:𝒢_24:热核 Z(t)=⟨e^{tλ}⟩` | CLOSURE | β = 1.25, c = -8.8291 | 0.99890 | T=24 |
| `DEFECT:Z^1 + V|0><0|:束缚态径向衰减 |ψ(r)| (V=-1.0)` | DEFECT | β = 0.75, c = +1.3697 | 0.97447 | D=1；L=201；直线段 → 指数衰减 → 局域长度 ξ |
| `DEFECT:Z^1 + V|0><0|:束缚态径向衰减 |ψ(r)| (V=1.0)` | DEFECT | β = 0.75, c = +1.3590 | 0.97360 | D=1；L=201；直线段 → 指数衰减 → 局域长度 ξ |
| `DEFECT:Z^1 L=101:极值态 IPR vs V` | DEFECT | β = 0.25, c = -3.1517 | 0.91001 | D=1；L=101 |
| `DEFECT:Z^2 + V|0><0|:束缚态径向衰减 |ψ(r)| (V=-1.0)` | DEFECT | β = 2.0, c = +0.0013 | 0.95830 | D=2；L=101；直线段 → 指数衰减 → 局域长度 ξ |
| `DEFECT:Z^2 + V|0><0|:束缚态径向衰减 |ψ(r)| (V=-4.0)` | DEFECT | β = 0.75, c = +2.5405 | 0.98193 | D=2；L=101；直线段 → 指数衰减 → 局域长度 ξ |
| `DEFECT:Z^2 + V|0><0|:束缚态径向衰减 |ψ(r)| (V=1.0)` | DEFECT | β = 0.25, c = +0.8060 | 0.99288 | D=2；L=101；直线段 → 指数衰减 → 局域长度 ξ |
| `DEFECT:Z^2 + V|0><0|:束缚态径向衰减 |ψ(r)| (V=4.0)` | DEFECT | β = 0.75, c = +2.5710 | 0.98762 | D=2；L=101；直线段 → 指数衰减 → 局域长度 ξ |
| `DEFECT:Z^2 L=161:极值态 IPR vs V` | DEFECT | β = 0.5, c = -3.1573 | 0.84102 | D=2；L=161 |
| `DEFECT:Z^2 L=21:极值态 IPR vs V` | DEFECT | β = 0.5, c = -1.8401 | 0.85691 | D=2；L=21 |
| `DEFECT:Z^2 L=41:极值态 IPR vs V` | DEFECT | β = 0.5, c = -2.2715 | 0.84884 | D=2；L=41 |
| `DEFECT:Z^2 L=81:极值态 IPR vs V` | DEFECT | β = 0.5, c = -2.7110 | 0.84537 | D=2；L=81 |
| `DEFECT:Z^3 + V|0><0|:束缚态径向衰减 |ψ(r)| (V=-1.0)` | DEFECT | β = 2.0, c = +0.0090 | 0.96792 | D=3；L=31；直线段 → 指数衰减 → 局域长度 ξ |
| `DEFECT:Z^3 + V|0><0|:束缚态径向衰减 |ψ(r)| (V=-4.0)` | DEFECT | β = 0.75, c = +0.6347 | 0.98587 | D=3；L=31；直线段 → 指数衰减 → 局域长度 ξ |
| `DEFECT:Z^3 + V|0><0|:束缚态能移 ΔE V>0（上带边 +2D）` | DEFECT | β = 0.5, c = -5.4161 | 0.92722 | D=3；L=31；ΔE<0 且 |ΔE| 单调增 = 束缚态；标度见结构分类 |
| `DEFECT:Z^3 L=11:极值态 IPR vs V` | DEFECT | β = 0.75, c = -0.9860 | 0.84987 | D=3；L=11 |
| `DEFECT:Z^3 L=17:极值态 IPR vs V` | DEFECT | β = 0.75, c = -1.1735 | 0.84542 | D=3；L=17 |
| `DEFECT:Z^3 L=25:极值态 IPR vs V` | DEFECT | β = 0.75, c = -1.3427 | 0.83959 | D=3；L=25 |
| `DEFECT:Z^3 L=37:极值态 IPR vs V` | DEFECT | β = 0.75, c = -1.5166 | 0.83164 | D=3；L=37 |
| `DEFECT:Z^4 + V|0><0|:束缚态径向衰减 |ψ(r)| (V=-1.0)` | DEFECT | β = 2.0, c = +0.0430 | 0.99272 | D=4；L=13；直线段 → 指数衰减 → 局域长度 ξ |
| `DEFECT:Z^4 + V|0><0|:束缚态径向衰减 |ψ(r)| (V=-4.0)` | DEFECT | β = 1.5, c = +0.1252 | 0.99249 | D=4；L=13；直线段 → 指数衰减 → 局域长度 ξ |
| `LEDGER:类权重 q_L:q_L vs L` | LEDGER | β = 1.25, c = +0.2317 | 0.99966 |  |

### `gapped`（14 条）

| 可观测量 | 扇区 | 拟合 | R² | 备注 |
|:--|:--|:--|--:|:--|
| `CLOSURE:𝒢_12:谱密度 ρ(λ)` | CLOSURE | gap_width=0.657834, gap_rel=0.05, method=dense-hist | — | T=12 |
| `CLOSURE:𝒢_12:谱密度形状` | CLOSURE | gap_width=0.657834, gap_rel=0.05 | — | T=12 |
| `CLOSURE:𝒢_14:度分布 P(deg)` | CLOSURE | gap_width=0.909091, gap_rel=0.0909 | — | T=14 |
| `CLOSURE:𝒢_14:谱密度 ρ(λ)` | CLOSURE | gap_width=1.008278, gap_rel=0.0638, method=dense-hist | — | T=14 |
| `CLOSURE:𝒢_14:谱密度形状` | CLOSURE | gap_width=1.294783, gap_rel=0.082 | — | T=14 |
| `CLOSURE:𝒢_16:度分布 P(deg)` | CLOSURE | gap_width=0.923077, gap_rel=0.0769 | — | T=16 |
| `CLOSURE:𝒢_16:谱密度 ρ(λ)` | CLOSURE | gap_width=2.067672, gap_rel=0.1064, method=dense-hist | — | T=16 |
| `CLOSURE:𝒢_16:谱密度形状` | CLOSURE | gap_width=1.272943, gap_rel=0.065 | — | T=16 |
| `CLOSURE:𝒢_18:度分布 P(deg)` | CLOSURE | gap_width=0.933333, gap_rel=0.0667 | — | T=18 |
| `CLOSURE:𝒢_18:谱密度 ρ(λ)` | CLOSURE | gap_width=1.835379, gap_rel=0.0851, method=dense-hist | — | T=18 |
| `CLOSURE:𝒢_18:谱密度形状` | CLOSURE | gap_width=1.086475, gap_rel=0.05 | — | T=18 |
| `CLOSURE:𝒢_20:度分布 P(deg)` | CLOSURE | gap_width=0.941176, gap_rel=0.0588 | — | T=20 |
| `CLOSURE:𝒢_22:度分布 P(deg)` | CLOSURE | gap_width=0.947368, gap_rel=0.0526 | — | T=22 |
| `CLOSURE:𝒢_24:度分布 P(deg)` | CLOSURE | gap_width=0.952381, gap_rel=0.0476 | — | T=24 |

### `exponential`（12 条）

| 可观测量 | 扇区 | 拟合 | R² | 备注 |
|:--|:--|:--|--:|:--|
| `CLOSURE:𝒢_10:热核 Z(t)=⟨e^{tλ}⟩` | CLOSURE | b = +5.1376 | 0.99711 | T=10 |
| `CLOSURE:𝒢_12:热核 Z(t)=⟨e^{tλ}⟩` | CLOSURE | b = +6.6311 | 0.99681 | T=12 |
| `CLOSURE:𝒢_14:热核 Z(t)=⟨e^{tλ}⟩` | CLOSURE | b = +8.0190 | 0.99651 | T=14 |
| `CLOSURE:𝒢_8:热核 Z(t)=⟨e^{tλ}⟩` | CLOSURE | b = +3.7368 | 0.99732 | T=8 |
| `CLOSURE:𝒢_T 系综:节点数 vs T` | CLOSURE | b = +0.6082 | 0.99906 |  |
| `DEFECT:Z^4 + V|0><0|:束缚态能移 ΔE V>0（上带边 +2D）` | DEFECT | b = +1.1593 | 0.93990 | D=4；L=13；ΔE<0 且 |ΔE| 单调增 = 束缚态；标度见结构分类 |
| `DEFECT:Z^4 L=11:极值态 IPR vs V` | DEFECT | b = +0.5839 | 0.84372 | D=4；L=11 |
| `DEFECT:Z^4 L=13:极值态 IPR vs V` | DEFECT | b = +0.6221 | 0.83584 | D=4；L=13 |
| `DEFECT:Z^4 L=7:极值态 IPR vs V` | DEFECT | b = +0.4792 | 0.86018 | D=4；L=7 |
| `DEFECT:Z^4 L=9:极值态 IPR vs V` | DEFECT | b = +0.5378 | 0.85230 | D=4；L=9 |
| `LEDGER:记录层 L1′:独立闭合类数 vs K_n` | LEDGER | b = +0.6014 | 0.99931 |  |
| `WORD:{±1}^n 全空间:平衡词计数增长` | WORD | b = +0.6628 | 0.99973 |  |

### `constant`（10 条）

| 可观测量 | 扇区 | 拟合 | R² | 备注 |
|:--|:--|:--|--:|:--|
| `CLOSURE:𝒢_T 系综:平均度 - T/2 vs T` | CLOSURE | a = 0.3544 | -0.00000 |  |
| `DEFECT:Z^1 + V|0><0|:束缚态能移 ΔE V<0（下带边 -2D）` | DEFECT | a = -3.243 | 0.00000 | D=1；L=201；ΔE<0 且 |ΔE| 单调增 = 束缚态；标度见结构分类 |
| `DEFECT:Z^2 + V|0><0|:束缚态能移 ΔE V<0（下带边 -2D）` | DEFECT | a = -2.13 | 0.00000 | D=2；L=101；ΔE<0 且 |ΔE| 单调增 = 束缚态；标度见结构分类 |
| `DEFECT:Z^3 + V|0><0|:束缚态能移 ΔE V<0（下带边 -2D）` | DEFECT | a = -1.285 | 0.00000 | D=3；L=31；ΔE<0 且 |ΔE| 单调增 = 束缚态；标度见结构分类 |
| `DEFECT:Z^4 + V|0><0|:束缚态能移 ΔE V<0（下带边 -2D）` | DEFECT | a = -0.5608 | 0.00000 | D=4；L=13；ΔE<0 且 |ΔE| 单调增 = 束缚态；标度见结构分类 |
| `LEDGER:记录层 L1′:对数缺口 ln(1-coverage) vs n` | LEDGER | a = -1.455 | 0.00000 |  |
| `WORD:平衡词 {±1}^10:两点关联 C(r)=⟨w_0 w_r⟩` | WORD | a = -0.1111 | 1.00000 | n=10 |
| `WORD:平衡词 {±1}^12:两点关联 C(r)=⟨w_0 w_r⟩` | WORD | a = -0.09091 | 1.00000 | n=12 |
| `WORD:平衡词 {±1}^14:两点关联 C(r)=⟨w_0 w_r⟩` | WORD | a = -0.07692 | 1.00000 | n=14 |
| `WORD:平衡词 {±1}^8:两点关联 C(r)=⟨w_0 w_r⟩` | WORD | a = -0.1429 | 1.00000 | n=8 |

### `declared`（10 条）

| 可观测量 | 扇区 | 拟合 | R² | 备注 |
|:--|:--|:--|--:|:--|
| `DEFECT:Z^1 缺陷:临界耦合 Vc vs 格点数 N` | DEFECT | — | — | D=1；Vc→0 = 任意弱耦合成键；Vc→常数 = 有阈值 |
| `DEFECT:Z^2 缺陷:临界耦合 Vc vs 格点数 N` | DEFECT | — | — | D=2；Vc→0 = 任意弱耦合成键；Vc→常数 = 有阈值 |
| `DEFECT:Z^3 缺陷:临界耦合 Vc vs 格点数 N` | DEFECT | — | — | D=3；Vc→0 = 任意弱耦合成键；Vc→常数 = 有阈值 |
| `DEFECT:Z^4 缺陷:临界耦合 Vc vs 格点数 N` | DEFECT | — | — | D=4；Vc→0 = 任意弱耦合成键；Vc→常数 = 有阈值 |
| `DEFECT:带边 DOS:带边指数汇总` | DEFECT | — | — |  |
| `DEFECT:缺陷算符:临界耦合汇总` | DEFECT | — | — |  |
| `WALK:Z^1 邻接算子 A 的谱:谱密度形状` | WALK | — | — | D=1；L=1025 |
| `WALK:Z^2 邻接算子 A 的谱:谱密度形状` | WALK | — | — | D=2；L=257 |
| `WALK:Z^3 邻接算子 A 的谱:谱密度形状` | WALK | — | — | D=3；L=65 |
| `WALK:Z^4 邻接算子 A 的谱:谱密度形状` | WALK | — | — | D=4；L=25 |

### `logarithmic`（10 条）

| 可观测量 | 扇区 | 拟合 | R² | 备注 |
|:--|:--|:--|--:|:--|
| `LEDGER:记录层 L1′:闭合类覆盖率 vs n` | LEDGER | b = +0.0317 | 0.33528 |  |
| `WALK:Z^1 行走占据 n_k(v):占据熵 H=-Σp ln p` | WALK | b = +0.5037 | 0.99979 | D=1 |
| `WALK:Z^1 行走占据 n_k(v):超额峰度 M4/M2²` | WALK | b = +0.1027 | 0.79589 | D=1 |
| `WALK:Z^2 行走占据 n_k(v):占据熵 H=-Σp ln p` | WALK | b = +1.0101 | 0.99977 | D=2 |
| `WALK:Z^2 行走占据 n_k(v):超额峰度 M4/M2²` | WALK | b = +0.0644 | 0.84512 | D=2 |
| `WALK:Z^3 行走占据 n_k(v):占据熵 H=-Σp ln p` | WALK | b = +1.5259 | 0.99954 | D=3 |
| `WALK:Z^3 行走占据 n_k(v):超额峰度 M4/M2²` | WALK | b = +0.0502 | 0.87568 | D=3 |
| `WALK:Z^4 行走占据 n_k(v):占据熵 H=-Σp ln p` | WALK | b = +2.0746 | 0.99946 | D=4 |
| `WALK:Z^4 行走占据 n_k(v):超额峰度 M4/M2²` | WALK | b = +0.0541 | 0.93476 | D=4 |
| `WORD:平衡词 {±1}^14:功率谱 S(q)` | WORD | b = +0.0000 | 0.85714 | n=14 |

### `too-short`（6 条）

| 可观测量 | 扇区 | 拟合 | R² | 备注 |
|:--|:--|:--|--:|:--|
| `CLOSURE:𝒢_10:度分布 P(deg)` | CLOSURE | — | — | T=10 |
| `CLOSURE:𝒢_10:谱密度 ρ(λ)` | CLOSURE | method=dense-hist | — | T=10 |
| `CLOSURE:𝒢_10:谱密度形状` | CLOSURE | — | — | T=10 |
| `CLOSURE:𝒢_12:度分布 P(deg)` | CLOSURE | — | — | T=12 |
| `CLOSURE:𝒢_8:度分布 P(deg)` | CLOSURE | — | — | T=8 |
| `CLOSURE:𝒢_8:谱密度 ρ(λ)` | CLOSURE | method=dense-hist | — | T=8 |

### `continuum/smooth`（4 条）

| 可观测量 | 扇区 | 拟合 | R² | 备注 |
|:--|:--|:--|--:|:--|
| `DEFECT:Z^1 带边 DOS:ρ(2D-δ) vs δ` | DEFECT | center_over_median=1.044, theory_exponent=-0.5 | — | D=1 |
| `DEFECT:Z^2 带边 DOS:ρ(2D-δ) vs δ` | DEFECT | center_over_median=0.99, theory_exponent=0.0 | — | D=2 |
| `DEFECT:Z^3 带边 DOS:ρ(2D-δ) vs δ` | DEFECT | center_over_median=1.032, theory_exponent=0.5 | — | D=3 |
| `DEFECT:Z^4 带边 DOS:ρ(2D-δ) vs δ` | DEFECT | center_over_median=1.127, theory_exponent=1.0 | — | D=4 |

### `continuum/center-divergent(van-Hove)`（3 条）

| 可观测量 | 扇区 | 拟合 | R² | 备注 |
|:--|:--|:--|--:|:--|
| `CLOSURE:𝒢_20:KPM 谱密度 ρ(λ)` | CLOSURE | center_over_median=9.677, method=KPM-Jackson | — | T=20 |
| `CLOSURE:𝒢_22:KPM 谱密度 ρ(λ)` | CLOSURE | center_over_median=8.586, method=KPM-Jackson | — | T=22 |
| `CLOSURE:𝒢_24:KPM 谱密度 ρ(λ)` | CLOSURE | center_over_median=17.182, method=KPM-Jackson | — | T=24 |

### `continuum/edge-vanishing+center-divergent(van-Hove)`（3 条）

| 可观测量 | 扇区 | 拟合 | R² | 备注 |
|:--|:--|:--|--:|:--|
| `CLOSURE:𝒢_20:谱密度形状` | CLOSURE | center_over_median=8.234 | — | T=20 |
| `CLOSURE:𝒢_22:谱密度形状` | CLOSURE | center_over_median=9.204 | — | T=22 |
| `CLOSURE:𝒢_24:谱密度形状` | CLOSURE | center_over_median=12.354 | — | T=24 |

### `delta`（1 条）

| 可观测量 | 扇区 | 拟合 | R² | 备注 |
|:--|:--|:--|--:|:--|
| `CLOSURE:𝒢_8:谱密度形状` | CLOSURE | dominant_fraction=0.6 | — | T=8 |

---

## §4 谱形状（分布类可观测量，逐条）

| 对象 | 标签 | 带 | 带边指数(lo/hi) | 其他 |
|:--|:--|:--|:--|:--|
| `WALK:Z^1 邻接算子 A 的谱:谱密度形状` | continuum/edge-divergent | [-1.987, 1.988] | -0.436 / -0.444 |  |
| `WALK:Z^2 邻接算子 A 的谱:谱密度形状` | continuum/smooth | [-3.975, 3.975] | 0.029 / 0.034 |  |
| `WALK:Z^3 邻接算子 A 的谱:谱密度形状` | continuum/edge-vanishing | [-5.956, 5.963] | 0.514 / 0.507 |  |
| `WALK:Z^4 邻接算子 A 的谱:谱密度形状` | continuum/edge-vanishing | [-7.887, 7.950] | 0.773 / 0.973 |  |
| `CLOSURE:𝒢_8:谱密度形状` | delta | [-3.884, 3.884] | — | 二部=True；简并: 5/10 不同值, 最大重数 6, 1%带宽团簇 5；λmax=3.884 |
| `CLOSURE:𝒢_10:谱密度形状` | too-short | — | — | 二部=False；r_distinct=0.3908；简并: 24/26 不同值, 最大重数 2, 1%带宽团簇 19；λmax=5.349 |
| `CLOSURE:𝒢_12:谱密度形状` | gapped | [-6.578, 6.578] | — | 二部=True；r_distinct=0.4422；简并: 65/80 不同值, 最大重数 6, 1%带宽团簇 31；λmax=6.925 |
| `CLOSURE:𝒢_14:谱密度形状` | gapped | [-7.527, 8.269] | — | 二部=False；r_distinct=0.3915；简并: 185/246 不同值, 最大重数 9, 1%带宽团簇 18；λmax=8.401 |
| `CLOSURE:𝒢_16:谱密度形状` | gapped | [-9.792, 9.792] | — | 二部=True；r_distinct=0.3328；简并: 577/810 不同值, 最大重数 74, 1%带宽团簇 7；λmax=9.841 |
| `CLOSURE:𝒢_18:谱密度形状` | gapped | [-10.535, 11.194] | — | 二部=False；r_distinct=0.3835；简并: 1709/2704 不同值, 最大重数 55, 1%带宽团簇 9；λmax=11.249 |
| `CLOSURE:𝒢_20:谱密度形状` | continuum/edge-vanishing+center-divergent(van-Hove) | [-13.154, 13.154] | 0.956 / 0.604 | 二部=True；λmax=12.598 |
| `CLOSURE:𝒢_22:谱密度形状` | continuum/edge-vanishing+center-divergent(van-Hove) | [-13.985, 14.545] | 1.822 / 1.191 | 二部=False；λmax=13.942 |
| `CLOSURE:𝒢_24:谱密度形状` | continuum/edge-vanishing+center-divergent(van-Hove) | [-15.915, 15.915] | 1.429 / 1.484 | 二部=True；λmax=15.243 |

---

## §5 让结构自己报名字（**后验**，与 §1–§4 分开）

§1–§4 全程只用结构标签，没有任何物理名字参与分类。下表是**事后**把结构对到已知物理概念，供人读；**它不改变上面的分类，也不构成断言**。

| 结构（§2 的标签） | 出现的对象 | 已知物理里的同名结构 | 等级 |
|:--|:--|:--|:--|
| `power-law` 指数 $-D/2$ | Γ 上行走的回返 $n_k(0)$、参与比 IPR、峰位概率 | 连续时间随机游走的回返衰减 ⟹ **谱维数** $d_s=D$（扩散） | 【数值证据】＋闭式核验 |
| `power-law` 指数 $+1$ | 二阶矩 $M_2(k)=\sum r^2p$ | **扩散律** $\langle x^2\rangle\propto t$ | 【数值证据】 |
| `logarithmic` | 占据熵 $H(k)=-\sum p\ln p$ | 熵的对数增长 $\tfrac D2\ln k+\text{const}$（维数读出） | 【数值证据】 |
| `constant` | 超额峰度 $M_4/M_2^2$ | **高斯普适性**（中心极限的不动点） | 【数值证据】 |
| `exponential` | 闭链图热核 $Z(t)=\langle e^{t\lambda}\rangle$ | **Perron/指数增长**，增长率 $=\lambda_{\max}$ | 【数值证据】 |
| `gapped` ＋ $\langle r\rangle_{\rm distinct}\approx0.39$ | 闭链图 $\mathcal G_T$ 的谱 | **Poisson 型能级统计（无能级排斥）** ＋ 团簇/简并谱（非随机矩阵） | 【数值证据】 |
| `continuum/edge-divergent` | $D{=}1$ 邻接谱 | **van Hove 奇点**（弧正弦律 $1/\sqrt{4-\lambda^2}$） | 【数值证据】 |
| `continuum/center-divergent(van-Hove)` | $D{=}2,3,4$ 与 $\mathcal G_T$ 的谱 | **van Hove 奇点**（态密度对数发散） | 【数值证据】 |
| `continuum/edge-vanishing` | $D{=}3,4$ 邻接谱 | 带边**幂律消失**（Weyl 律的带边指数） | 【数值证据】 |
| `delta` ＋ `parity` | 词层关联 $C(r)=\langle w_0w_r\rangle=\delta_{r0}$、回返奇偶零 | **白噪声（无关联）** ／ 二部图奇偶选择定则 | 【数值证据】＋精确 |
| `constant`（饱和） | 记录层覆盖率 $\to0.7657$ | **饱和/稳态**（非幂律、非指数） | 【数值证据】 |
| 带边 DOS 指数 $=D/2-1$ | 缺陷算符 $A+V|0\rangle\langle0|$ 的带边 | **晶格 Green 函数在带边的行为**：指数 $<0$ ⟹ 任意弱耦合成键；$>0$ ⟹ **有阈值** | 【数值证据】＋理论一致 |

**注意**：§5 的最后一行的「成键阈值」是[`L2_REGIONS.md`](L2_REGIONS.md) §3.1 的输入；那里用它把 $D\le2$ 与 $D\ge3$ 定性分开。

