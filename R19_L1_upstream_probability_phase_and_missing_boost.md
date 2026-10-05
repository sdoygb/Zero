# R19 · L1 上游重排：概率与相位在位，缺的是洛伦兹 boost

**日期**：2026-10-02  
**性质**：L1 上游依赖图重排 ＋ 一条已证 no-go（旋转双覆盖不能提供 boost）＋ 一组验收条件。不关闭 L1，不新增物理参数。  
**唯一目标**：Jacobson 2016 的最后一跳
$$
K_{B,a}\longrightarrow 2\pi B_B .
$$
**依赖**：[`Z0`](Z0_zero_never_rests_single_axiom.md)、[`Z13`](Z13_zero_foundation_missing_principle.md)、[`Z14`](Z14_closure_cyclic_order_base_theorem.md)、[`G27`](G27_purification_attempt.md)、[`G29`](G29_probability_as_derived_not_postulated.md)、[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md)、[`G77`](G77_staggered_coupling_from_A5.md)、[`R12`](R12_zero_native_gap_filling.md)、[`R13`](R13_L1_strong_resolvent_attempt.md)、[`R15`](R15_zcar_double_cover_and_zstress_scale.md)、[`R17`](R17_L1_critical_path_and_L5_gate.md)、[`R18`](R18_L5_cert_dimension_normalization_gate.md)、[`D44`](/Users/oygb/Downloads/modular-equilibrium/derivations/D44_poincare_closure_obstruction.md)。  
**核验**：[`R19_check.py`](R19_check.py)。

$$
\boxed{
\begin{aligned}
&\text{概率、复相位与模流都已经从 Zero 生出来；}\\
&\text{卡住 L1 的是：这些流都是 Euclidean 旋转或交换平移型，}\\
&\text{而 Bisognano--Wichmann 要的是非阿贝尔洛伦兹 boost。}\\
&\text{缺的正是 Wick 解析延拓，加楔形几何，加 }2\pi\text{ 归一化落点。}
\end{aligned}}
$$

> **一句话**：上游一直以为卡在“Zero 还没长出概率和相位”。重新排一遍依赖图后，概率（`G29`）、复相位与模流（`G62`）、拓扑相位 `-1`（`Z14`／`R15`）都已长出；卡住的是这些流的方向与代数——它们是紧的旋转流与交换的平移流，而 L1 的右端 $B_B$ 是非紧、非阿贝尔的 boost。

---

## §0 判决摘要

| 问题 | 本轮判定 | 依据 |
|:--|:--|:--|
| 概率是否已由 Zero 导出 | **是**；条件于一条层结构 $\pi$ | `G29` |
| 复振幅／Hilbert／模流 | **是**；GNS ＋ KMS | `G62` |
| 拓扑相位（`-1` holonomy） | **是**；闭合序双覆盖 | `Z14`、`R15` |
| 这些是否已经是 $B_B$ | **否** | §3，命题 R19.2 |
| L1 的真正缺口 | Wick 解析延拓 ＋ 楔形几何 ＋ $2\pi$ 落点 | §4 |
| 闭合双覆盖能否补上 boost | **不能**（no-go） | §3，命题 R19.2 与推论 R19.3 |
| 下一步 | `Z-BOOST` 三条验收条件 | §5 |
| L1 当前状态 | **仍未关闭** | 全篇 |

---

## §1 上游重排：两条线，一个合并点

### 1.1 概率线（`G29`）

$$
\texttt{Z0③}
\Longrightarrow
\text{全分支计数}
\Longrightarrow
\mu
\xrightarrow{\ \pi\ }
\omega=\pi_*\mu
\Longrightarrow
K=-\log\omega .
$$

要点：计数测度 $\mu$ 由 Z0 的全分支＋整数重数给出；概率是它在层结构 $\pi$ 下的推前，**不是公设**。非均匀的 $\omega$ 给出非平凡的 $K$。

### 1.2 相位线，有两段，必须分开

**(a) 复相位（量子叠加相位）**：闭合序＋原生 $\pm$ 给出 $D_L$，进而给出非对易代数
$$
\mathcal A_T=M_2(\mathbb C)\otimes\mathbb C^{\,T+1},
$$
配一个忠实态后由 GNS 给出复 $*$-代数、复振幅与模流
$$
\sigma_t(A)=\rho^{it}A\rho^{-it},\qquad K=-\log\rho .
$$
（`G27`、`G62`；KMS 数值残差 $\approx3.6\times10^{-16}$。）

**(b) 拓扑相位（spin holonomy）**：闭合词位给出循环序 $\mathbb Z_L\subset SO(2)$，其双覆盖给出中心 $\mathbb Z_2$ 与
$$
\bigl(U_{2\pi/L}\bigr)^{L}=-\mathbb I .
$$
（`Z14`、`R15`。）

### 1.3 合并点

两条线共用两个上游对象：

1. **同一条层结构 $\pi$**：它既给概率测度 $\omega$，又给 GNS 需要的忠实态；
2. **同一条闭合循环序**：它既给非对易代数 $\mathcal A_T$，又给中心 $\mathbb Z_2$ 的 `-1` holonomy。

因此“零和宇宙上游只有一个具名输入”这句话是准确的：那个输入就是 $\pi$（`Z13` 的 E5a）。

---

## §2 节点状态表

| 节点 | 内容 | 当前逻辑状态 | 依据 |
|:--|:--|:--|:--|
| Z0③ 无偏好 | 唯一公理 | 公理 | `Z0` |
| 全分支＋整数重数 | 由 Z0③ 导出 | **已证** | `Z2` |
| 计数测度 $\mu$ | 每支等权 | **已证** | `Z2`、`G29` |
| 层结构 $\pi$ | 微观态到宏观类的粗粒化 | **具名输入**（E5a） | `Z13`、`G29` |
| 概率 $\omega=\pi_*\mu$ | 推前测度 | **条件证成** | `G29` |
| 非对易代数 $\mathcal A_T$ | $M_2\otimes\mathbb C^{T+1}$ | **已证（有限维）** | `G27`、`Z14` |
| 复振幅／Hilbert | GNS | **条件证成** | `G62` |
| 模流 $\sigma_t$ | $\rho^{it}\cdot\rho^{-it}$ | **条件证成** | `G62` |
| 拓扑相位 `-1` | 双覆盖中心 | **已证** | `Z14`、`R15` |
| 自由费米结构／半满 | 单费米点与交错耦合 | **识别** | `G77`、`R13` |
| 归一化落点 | $h$ 与 $\beta$ 的比例常数 | **开放**（实测 $\approx\pi^{2}/3$） | `R15` |
| 收敛分析 | 公共核心与尾项 | **开放** | `R13`（`Z-CORE`／`Z-TAIL`） |
| 几何 boost $B_B$ | 区域保持的非阿贝尔流 | **开放** | `R12` §2.2、`D44` |

---

## §3 重新评估：概率与相位都不是缺口

### 命题 R19.1（上游量子机器已闭合到“只需一个 $\pi$”）【审计结论】

概率、复振幅、Hilbert 空间、模流与拓扑 `-1` 相位都已由 Zero 给出，各自只条件于同一条具名输入 $\pi$。因此 L1 的困难**不是**“Zero 生不出概率和相位”。

**边界**：这不是说 L1 变得容易。它只否定一个错误的定位——把 L1 当作“底层还缺概率／相位”的问题。缺的东西在更上面一层。

### 命题 R19.2（旋转双覆盖不提供洛伦兹 boost）【已证】

设 $\sigma_z=\operatorname{diag}(1,-1)$。闭合结构给出的旋转**升格**为
$$
U(\theta)=\exp\!\Bigl(i\,\frac{\theta}{2}\,\sigma_z\Bigr),
\qquad
U(2\pi)=-\mathbb I .
\qquad\text{(R19-1)}
$$
洛伦兹 **boost** 流为
$$
\Lambda(\eta)=\exp\!\Bigl(\frac{\eta}{2}\,\sigma_z\Bigr),
\qquad \eta\in\mathbb R .
\qquad\text{(R19-2)}
$$
则
$$
\Lambda(\eta)\ne-\mathbb I\ \ \forall\eta\in\mathbb R,
\qquad
\{\Lambda(\eta)\}\ \text{非紧},
\qquad
\{U(\theta)\}\ \text{紧}.
\qquad\text{(R19-3)}
$$
而 Wick 转动给出
$$
\Lambda(i\theta)=U(\theta),
\qquad\text{故 }-\mathbb I=\Lambda(2\pi i).
\qquad\text{(R19-4)}
$$

**证明**：$\Lambda(\eta)$ 的本征值是 $e^{\eta/2}$ 与 $e^{-\eta/2}$，对实 $\eta$ 都严格为正；$-\mathbb I$ 的本征值为 $-1$，故不可能。$\eta\mapsto e^{\eta/2}$ 无界而 $\theta\mapsto U(\theta)$ 取遍紧群 $U(1)$，给出 (R19-3)。最后 $\Lambda(i\theta)=e^{i\theta\sigma_z/2}=U(\theta)$，取 $\theta=2\pi$ 得 (R19-4)。$\square$

**读法**：`-1` 只出现在**虚 rapidity**（Euclidean 侧）。因此“$2\pi$ 旋转 $=-\mathbb I$”这条 Zero 已经证成的信息，是**旋转侧**的 spin 结构；它不能给出 boost 侧的 rapidity 归一化，也就是不能单独给出 $2\pi B_B$ 里的那个 $2\pi$。

### 推论 R19.3（`R12.5` 的精确读法）【审计结论】

`D44`／命题 R12.5 证明半侧模包含只给出相互对易的正能平移
$$
[P_i,P_j]=0 .
\qquad\text{(R19-5)}
$$
R19.2 补上另一半：闭合循环序给出的是**紧的 Euclidean 旋转**，它的中心 $\mathbb Z_2$ 是 spin holonomy，而不是 boost。

两条合起来给出一个 no-go：

> 现有的 Zero 原生流（交换平移 ＋ 紧旋转）**不能**直接给出 $B_B$。它们不是“还没算出来”，而是代数类型不对：一个交换、一个紧；目标是既非交换又非紧。

这正是 `R17-STOP` 允许进入主线的那一类结果——它证明了一条候选路线（用闭合双覆盖直接供给 boost）不可满足，而不是新增一个待办标签。

---

## §4 L1 真正缺的三件

| 缺口 | 内容 | 类型 | 现有依据 |
|:--|:--|:--|:--|
| `Z-WICK` | 从 Euclidean 旋转流到 Lorentzian boost 流的解析延拓 | 原理 | §3；`G62` 已有 GNS 侧解析延拓 |
| `Z-WEDGE` | 1+1D 区间／Euclidean 环到四维楔形／双锥 | 几何 | `R13` §4.3（原 `Z-CONF`） |
| `Z-STRESS-2π` | $h$ 与 $\beta$ 的比例常数必须落在 $2\pi$ | 归一化 | `R15` §4 |

**三者不是新的独立缺项，而是已有标签的重排**：`Z-WICK` 吸收“几何 boost 缺口”，`Z-WEDGE` 就是原 `Z-CONF`，`Z-STRESS-2π` 是原 `Z-STRESS` 的落点要求。如此重排后，与 boost 直接相关的部分从“若干分散标签”收成一个上游对象。

**与收敛分析的关系**：`Z-CORE`／`Z-TAIL`（`R13`）仍是独立的分析型缺口，原样保留。

---

## §5 补法：`Z-BOOST` 的三条验收条件

候选路线（Osterwalder–Schrader 重构型）：

**`Z-BOOST`**：对计数推前态 $\omega$ 与闭合旋转流 $\alpha_\theta$，证明下列三条。

**A1（几何一致性，KMS-$2\pi$）**：几何旋转流 $\alpha_\theta$（几何周期 $2\pi$）与 $\omega$ 的模流 $\sigma_t$ 满足

$$
\sigma_{\theta/2\pi}=\alpha_\theta ,
$$

即 $\omega$ 关于几何旋转流是周期 $2\pi$ 的 KMS 态。

**A2（反射正性）**：在把闭合序读成虚时间的方向上，$\omega$ 满足 Osterwalder–Schrader 反射正性。

**A3（非阿贝尔与区域保持）**：Wick 延拓后的流是区域保持的单参数流，且其生成元 $K_B$ 与平移 $P$ 不对易，

$$
[K_B,P]\ne0 .
$$

**若 A1–A3 成立**，OS 重构给出洛伦兹结构，boost 流以 $\omega$ 为真空／KMS 态，于是 $K_B$ 的模流等于该 boost 流；L1 只剩 $2\pi$ 落点与 `Z-WEDGE` 提升。

**诚实边界**：A1、A2、A3 目前都**未证**。本条给出的是精确目标，不是结论。三点各有明确的可判定性：

1. **A1 已有现成负信号**：`R15` 实测 $h$ 与 $\beta$ 的比例常数稳定在 $\approx\pi^{2}/3=3.2899$，而 $2\pi=6.2832$。二者之比 $\approx6/\pi\approx1.9099$。这说明当前原生归一化尚未落在 A1 要求的点上；它既是否证点，也是 `Z-STRESS-2π` 的精确目标。
2. **A2 是线性代数／正性检验**：在有限维截断上可判定，不需要四维几何。
3. **A3 与 `D44` 相容**：`D44` 只给交换平移；A3 要求的是一条**不同**的流（区域保持 boost），必须构造，不能从半侧平移里取。

---

## §6 下一步与止损

**第一个可判定子问题**（打 A1）：在自由费米子支线上，计算闭合旋转流 $\alpha_\theta$ 与模流 $\sigma_t$ 的关系，判定 $\sigma_{\theta/2\pi}=\alpha_\theta$ 是否成立。等价形式是：$h$ 与 $\beta$ 的比例常数是否等于 $2\pi$。

1. 若 A1 成立，转入 A2（反射正性），它可在有限维截断上判定；
2. 若 A1 失败且无补救，按 `R17-STOP` 判定“闭合旋转路线不能供给 boost 归一化”，据此换锚或改写平衡条件，而不是继续堆标签。

**不做的事**：

1. 不把 A1–A3 写成已证；
2. 不再新增表示层 `Z-*` 标签；`Z-WICK`／`Z-WEDGE`／`Z-STRESS-2π` 只是既有缺口的分面；
3. 不把 $\pi^{2}/3$ 与 $2\pi$ 的差距当作“差一个常数就跳过”——它是本轮定位出来的判定点。

**L1 状态**：$K_B\to2\pi B_B$ 仍未关闭；本轮只改它的上游定位。

> **后续（[`R20`](R20_A1_verdict_shape_holds_constant_fails.md)）**：A1 已在 1+1D 自由费米子支线上数值判定——**形状成立**（近邻与二次型结构都与 $\xi$ 对上，残差外推趋于约 $1\%$），但常数约 $3.3$–$3.4$，不是 $2\pi$，缺口是一个约 $1.85$ 倍的**纯归一化因子**。因此 A1 归约为 `Z-STRESS-2π` 的归一化子问题（从同一格点模型的 $v_F$ 导出），A2／A3 状态不变。

> **后续校正（[`R21`](R21_vf_normalization_resolves_the_r20_factor.md)）**：R20 的常数额不是独立参数。外部连续极限给出正确因子 $2\pi/v_F$；本模型 $v_F=2$，故连续系数应为 $\pi$。A1 的目标式必须把裸几何算子与物理 boost 生成元分开：$K_B\to(2\pi/v_F)B_{\rm bare}$，或等价地把 $v_F$ 吸收进 $B_B$。旋转双覆盖 no-go 不变，J1 仍未关闭。

> **后续校正（[`R22`](R22_principal_symbol_vs_r20_estimator.md)）**：A1 的测试必须限制到 ETP $T_N$ 的近零谱子空间并读取主符号。R20 的键中点 $B_N$ 与任意光滑向量会混入带边分量，故其 $3.394$／$0.9255$ 不是 $\pi$ 的有效估计。A1 的符号层现在条件证成；`Z-CORE`／`Z-TAIL` 与四维几何仍开放。

---

## §7 对既有状态的影响

| 文档 | 更新 |
|:--|:--|
| [`R12`](R12_zero_native_gap_filling.md) | §2.2 的 boost 缺口细化为“紧旋转／交换平移 vs 非紧非阿贝尔 boost” |
| [`R15`](R15_zcar_double_cover_and_zstress_scale.md) | $\pi^{2}/3$ 的读法补充为 A1 的现成否证点候选 |
| [`R13`](R13_L1_strong_resolvent_attempt.md) | 原 `Z-CONF` 归入 `Z-WEDGE`；`Z-CORE`／`Z-TAIL` 不变 |
| [`R17`](R17_L1_critical_path_and_L5_gate.md) | 本条属于 `R17-STOP` 允许的 no-go 类主攻 |
| [`STATUS.md`](STATUS.md) | 新增 §2.19 |
| [`R0`](R0_publication_theorem.md) | 外部／方向审计范围扩到 `R19` |

**范围边界**：命题 R19.2 只使用 $2\times2$ 矩阵的初等事实与 Lie 代数的同构不变量，不依赖任何四维几何；它给出的是**缺口定位**，不是 L1 的证明或对四维 L1 的否证。
