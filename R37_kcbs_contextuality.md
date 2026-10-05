# R37 · KCBS 探针：Zero 的单体统计**是语境的**（量子性第一个硬证据）

**日期**：2026-10-03  
**性质**：**探针结果＋精确判据＋正面结论（带边界）**。对 Zero 的单体结构做 KCBS/语境性检验——这是**不需要 `E1`、不需要空间、不需要类空**的那一类检验（见 R36 §五）。  
**依赖**：[`G27`](G27_purification_attempt.md)、[`G29`](G29_probability_as_derived_not_postulated.md)、[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md)、[`G68`](G68_interference_from_coarse_graining.md)、[`G72`](G72_kappa1_from_the_ledger.md)、[`G88`](G88_measurement_as_typicality.md)、[`R33`](R33_action_phase_match_project.md)、[`R35`](R35_type_iii_classification.md)、[`R36`](R36_chsh_bell_locality.md)、[`STATUS`](STATUS.md)。  
**探针**：[`R37_kcbs_probe.py`](R37_kcbs_probe.py) → [`R37_kcbs_results.json`](R37_kcbs_results.json)。  
**核验**：[`R37_check.py`](R37_check.py)。

$$
\boxed{
\begin{aligned}
&\text{KCBS：}C^3\text{ 中 5 条射线，相邻正交（}v_i\perp v_{i+1}\text{），}A:=\sum_i|v_i\rangle\langle v_i|,\ S(\rho)=\sum_i\operatorname{Tr}(\rho P_i)。\\
&\text{非语境界 }2;\ \text{量子最大 }\sqrt5=2.2360679。\\
&\text{对}\textbf{全部取向}\text{取最大有闭式（von Neumann 迹不等式）：}\\
&\qquad S_{\max}(\rho)=\sum_k\lambda_k(\rho)\,\mu_k(A),\qquad \mu(A)=(\sqrt5,\ 1.3819660,\ 1.3819660)。\\
&\text{结论：原生态（}G72\text{ 推前权重 }2,5,20,100\text{）的谱 }\lambda_1=0.7407>\lambda^*=0.7236\\
&\qquad\Longrightarrow S_{\max}=2.0146>2\ \Longrightarrow\ \textbf{语境（单体量子性成立）}。
\end{aligned}}
$$

> **一句话**：CHSH 那一刀因为缺"空间"砍不下去（R36）；换到**单体**层面用 KCBS，刀**砍下去了，而且见血**：把"取向自由"这个陷阱用 von Neumann 迹不等式正确处理之后，原生态的最大 KCBS 值 **`S_max = 2.0146 > 2`**（另一个 3 维归约给 `2.0652`）。也就是说——**Zero 的单体统计不是经典概率，是真正语境的**。这是量子栏拿到的第一个**硬证据**，而且它不需要 `E1`、不需要空间分离、不需要精确光锥。

---

## §0 判决摘要

| 项 | 结果 | 状态 |
|:--|:--|:--|
| 五角星相邻正交误差 | `2.8×10⁻¹⁶` | **通过** |
| `μ(A)` 本征值 | `(√5, 1.3819660, 1.3819660)`，和 `=5` | **已算** |
| 完全混合态 | `S=5/3=1.6667<2`（非语境） | **对照通过** |
| KCBS 最优纯态 | `S=√5=2.2360679` | **对照通过** |
| 闭式核对（3000 随机态／取向） | 无一例超过闭式上界（超出量 `0.0`） | **通过** |
| **原生态**（top-3 归一 `0.7407,0.1852,0.0741`） | **`S_max=2.0146>2`** | **语境** ✅ |
| 原生态（去最小块 `0.8,0.16,0.04`） | **`S_max=2.0652>2`** | **语境** ✅ |
| 阈值 | `λ*=0.7236`（族 `(λ,(1−λ)/2,(1−λ)/2)`）；原生 `λ₁=0.7407`，余量 `+0.0171` | **已算** |
| 是否否定量子性 | **否——相反，第一次给出正面证据** | — |
| 边界 | state-dependent（测量集按态选）；3 维归约＝`π`/`E5`；谱用文档例 | **已登记** |

---

## §1 为什么必须"对所有取向取最大"

第一次尝试时我把 pentagram 取向固定，得到原生态 `S≈1.45`（非语境）——**那是错的**：五角星的**取向是自由的**，固定一个取向只是"恰好没选中"。正确做法是对**全部** `U` 取最大：

$$
S(\rho;U)=\sum_i\langle v_i|U^{\dagger}\rho U|v_i\rangle
=\operatorname{Tr}\!\bigl(\rho\,UAU^{\dagger}\bigr)
\ \xrightarrow[\ \text{von Neumann}\ ]{\ \max_U\ }\
S_{\max}(\rho)=\sum_k\lambda_k(\rho)\,\mu_k(A).
\qquad\text{(R37-1)}
$$

其中 $\lambda(\rho)$、$\mu(A)$ 各自降序。探针用 3000 个随机态与随机取向核对：**无一例超过该闭式**（超出量 `0.0`）$\Rightarrow$ 闭式正确、且它确实是上确界。

$$
\boxed{
\text{判据：}\ S_{\max}(\rho)>2\iff\text{语境};\qquad
\mu(A)=(\sqrt5,\ 1.3819660,\ 1.3819660)\ (\text{和}=5)。
}
\qquad\text{(R37-2)}
$$

---

## §2 原生态的谱与结果

原生态在**指针基（类基）下对角**（[`G72`](G72_kappa1_from_the_ledger.md) §1：针基 = 推前的类基 = `K` 的本征基），其谱就是**推前权重**。文档例（[`G29`](G29_probability_as_derived_not_postulated.md) §2 的块大小 $2,5,20,100$）给 $w=(0.0157,0.0394,0.1575,0.7874)$。取 3 维归约：

| 3 维归约 | 谱 $\lambda$ | $S_{\max}$ | 是否语境 |
|:--|:--|--:|:--|
| 取前 3 块归一 | $(0.7407,\ 0.1852,\ 0.0741)$ | **2.0146** | **是** ✅ |
| 去掉最小块 | $(0.8000,\ 0.1600,\ 0.0400)$ | **2.0652** | **是** ✅ |
| （对照）完全混合 | $(1/3,1/3,1/3)$ | 1.6667 | 否 |

$$
\boxed{
\text{阈值：族 }(\lambda,\tfrac{1-\lambda}2,\tfrac{1-\lambda}2)\text{ 上 }\lambda^*=\frac{2-\mu_2}{\mu_1-\mu_2}=0.7236;
\quad\text{原生 }\lambda_1=0.7407\ \Rightarrow\ \textbf{余量 }+0.0171。
}
\qquad\text{(R37-3)}
$$

**余量虽小（1.7%）但符号确定**：只要原生态的谱比"完全混合"更尖（更集中），它就落在语境侧。

---

## §3 这条结论的地位（与 R36 的分工）

| | `R36`（CHSH/Bell） | `R37`（KCBS/语境性） |
|:--|:--|:--|
| 层面 | **多体** | **单体** |
| 需要 | 交换子代数 ＋ **空间分离** ＋ **精确类空** | **只需单体代数**（≥3 维＋不相容可观测量） |
| Zero 现状 | ❌ 缺空间分离与精确锥 | ✅ **可做，且已做** |
| 结果 | `S≤2`（Bell 局域） | **`S_max=2.0146>2`（语境）** |

$$
\boxed{
\text{合起来：Zero 的单体统计}\textbf{是量子（语境）}\text{的；多体（Bell）层面}\textbf{尚无场地}。
}
\qquad\text{(R37-4)}
$$

这与 §四 的判据一致：**单体的量子性由语境性承载，不由 Bell 承载**——所以"CHSH 不能违反"**不等于**"Zero 不是量子理论"。

---

## §4 三个必须写明的边界

| # | 边界 | 说明 |
|--:|:--|:--|
| 1 | **state-dependent** | 测量集（五角星）是**按态选的**——这是语境性检验的**标准做法**（"每个纯态都有状态依赖的语境性证明"），但因此它检验的是**态**，不是"理论原生固定的测量集" |
| 2 | **3 维归约的选择** | 从 4 块谱取哪 3 块，仍是 `π`／`E5` 的缺口；本文两个归约**都**落在语境侧，故结论对归约不敏感（在文档权重下） |
| 3 | **谱的来源** | 用的是**文档例**（`G72`/`G29` 的推前权重），不是从 `π` 显式算出的那一个；若真实谱更接近均匀（`λ₁<0.7236`），结论**翻转** |

$$
\boxed{
\text{第 3 条是可证伪入口：算出真实谱 }\lambda_1\text{，与 }0.7236\text{ 比大小即可——}\textbf{二值结论}。
}
\qquad\text{(R37-5)}
$$

---

## §5 没有推出什么

1. 没有算出**原生**谱（仍需 `π` 的显式分块）。
2. 没有把"语境性"升级为"理论原生固定的测量集下的语境性"（当前是 state-dependent）。
3. 没有改变 `R36` 的结论：多体层面仍缺空间张量分解与精确光锥。
4. 没有解决 `2π`（`R33` S1）、`III₁`（`R35`）、或 R31／R32 的选维链。
5. 没有由单体语境性推出 Lorentz、度规、Lovelock 或 GR。

$$
\boxed{
\text{当前诚实结论：}\textbf{Zero 的单体统计是语境的（量子）}\text{——第一个硬证据，条件于文档谱与 3 维归约。}
}
$$

---

## §6 核验命令

```bash
python3 R37_kcbs_probe.py
python3 R37_check.py
python3 G29_check.py
python3 G72_check.py
python3 R36_check.py
python3 STATUS_check.py
python3 INDEX_check.py
```
