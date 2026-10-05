# R21 · R20 归一化缺口的判定：因子是 $2\pi/v\_F$，不是新增的 $1.85$

**日期**：2026-10-02
**性质**：对 R20 的常数判定作归一化复核。修正 **R20 判定标号 A1**（与 A0–A5 的历史条款 A1 同名而异指）的比较口径；不关闭 L1，不新增物理参数。
**目标**：判定 R20 的“实测系数约 $3.3$–$3.4$、目标 $2\pi$、缺口约 $1.85$ 倍”到底是新常数，还是费米速度／应力归一化。
**依赖**：[`R8`](R8_jacobson_entanglement_equilibrium_completion.md)、[`R13`](R13_L1_strong_resolvent_attempt.md)、[`R15`](R15_zcar_double_cover_and_zstress_scale.md)、[`R19`](R19_L1_upstream_probability_phase_and_missing_boost.md)、[`R20`](R20_A1_verdict_shape_holds_constant_fails.md)、[`G77`](G77_staggered_coupling_from_A5.md)。
**核验**：[`R21_check.py`](R21_check.py)。
**外部定理**：Eisler–Tonni–Peschel, [arXiv:1902.04474](https://arxiv.org/abs/1902.04474)，第 2–3 节与结论。

> **后续校正（[`R22`](R22_principal_symbol_vs_r20_estimator.md)）**：本轮的 $2\pi/v\_F=\pi$ 结论保留，但 R20 的 $3.394$ 不是该常数的有限尺寸估计。外部定理用的 $T\_N$ 是 (R21-3)，R20 的键中点 $B\_N$（$l=N-1$）不是同一个矩阵；任意光滑向量也不属于 $T\_N$ 的近零模。R22 因此把 $0.9255$ 从“有限尺寸残差”降为估计器／测试子空间偏差。

$$

\begin{aligned}
&\text{R20 的形状结论保留：模 Hamiltonian 的剖面就是抛物线。}\\
&\text{原 A1 的常数目标写错了：在 }v_F=2\text{ 的格点模型中应比较 }\pi\text{，不是 }2\pi。\\
&\text{正确关系是 }h\simeq\frac{2\pi}{v_F}B_{\rm bare}\text{；当 }v_F=2\text{ 时，}\frac{2\pi}{v_F}=\pi。\\
&\text{R20 的 }1.85\text{ 倍是“相对错误目标 }2\pi\text{”的读数；其中 }0.9255\text{ 不是待补常数。}
\end{aligned}
$$

> **一句话**：R20 把格点能量归一化与连续共形 boost 的 $2\pi$ 放在同一坐标里比较，漏掉了 $v\_F$。补上 $v\_F$ 后，正确的连续主符号是 $\pi=2\pi/v\_F$；所谓 $1.85$ 倍不是一个待导出的 Zero 常数。R22 进一步指出，R20 的约 $8\%$ 残差来自估计器与测试子空间，不是要补的 $\pi$ 修正。

---

## §0 判决摘要

| 问题 | 本轮判定 | 依据 |
|:--|:--|:--|
| R20 的形状结论 | **保留**：近邻与二次型形状都是抛物线 | `R20` §2–§4 |
| R20 的常数目标 | **修正**：原 A1 用的裸 $2\pi$ 没有包含 $v\_F$ | §2–§3 |
| 正确连续系数 | **条件证成**：半满自由费米子区间为 $\pi=2\pi/v\_F$，其中 $v\_F=2$ | §2，命题 R21.1 |
| “缺口 $1.85$ 倍” | **解释为口径错误**：$1.851\approx2\times0.9255$，但 $0.9255$ 无定理地位（见 R22） | §3，命题 R21.2 |
| 是否新增 Zero 参数 | **否** | §3 |
| 对 L1 的影响 | **不关闭 L1**；归约归一化由“神秘常数”改为“应力／速度口径” | §4 |
| 是否还要强行把裸格点系数调成 $2\pi$ | **不应**；除非显式重标度 $v\_F=1$ 或把 $v\_F$ 吸收进 $B\_B$ 的定义 | §4 |

---

## §1 两种归一化不能混用

### 1.1 R20 比较的量

R20 取无穷半满自由费米子链的相关矩阵

$$
C_{ij}=\frac{\sin\!\bigl(\pi(i-j)/2\bigr)}{\pi(i-j)},
\qquad C_{ii}=\frac12,
\qquad
h=\log\frac{1-C}{C},
\qquad\text{(R21-1)}
$$

并与键权

$$
\beta_b=x_b\frac{l-x_b}{l},
\qquad x_b=b+\frac12
\qquad\text{(R21-2)}
$$

作比较。R20 的结论是 $h\simeq A\beta$，$A\approx3.3$–$3.4$。

### 1.2 外部连续极限实际比较的量

Eisler–Tonni–Peschel 对无限链上半满区间的纠缠 Hamiltonian 给出（只取主导连续项，差一个整体符号／添加常数）

$$
H_{\rm ent}\simeq \pi N\,T,
\qquad
T_{b,b+1}=T_{b+1,b}=\frac{b+1}{N}\left(1-\frac{b+1}{N}\right),
\qquad\text{(R21-3)}
$$

其中 $T$ 是最近邻几何矩阵，而高次 hopping 在连续极限里不贡献局部速度。若把物理坐标写成 $x=(b+1)a$、$l=Na$，则

$$
N\,T_{b,b+1}\simeq \frac{l}{a^2}\frac{x(l-x)}{l}
=\frac{l}{a^2}\beta(x).
\qquad\text{(R21-4)}
$$

因此 R20 的 $h$ 在连续极限应趋向 $\pi$，而不是 $2\pi$。

---

## §2 费米速度从哪里进入

### 2.1 格点色散

取

$$
\varepsilon(k)=-2t\cos(ka),
\qquad
k_Fa=\frac{\pi}{2}.
\qquad\text{(R21-5)}
$$

则

$$
v_F=\left|a\,\frac{d\varepsilon}{dk}\right|_{k_F}
=2ta\,\sin(k_Fa).
\qquad\text{(R21-6)}
$$

在 R13／R20 的标准归一化 $t=1,\ a=1$ 下，

$$
v_F=2.
\qquad\text{(R21-7)}
$$

所以连续关系中的

$$
\frac{2\pi}{v_F}
=
\pi
\qquad\text{(R21-8)}
$$

正好就是外部定理 (R21-3) 的系数。

### 2.2 任意填充的写法

Eisler–Tonni–Peschel 的更一般结论可写成

$$
H_{\rm ent}\simeq
-\frac{\pi N}{\sin(k_Fa)}\,T.
\qquad\text{(R21-9)}
$$

把 (R21-6) 代入：

$$
\frac{\pi}{\sin(k_Fa)}
=
\frac{2\pi ta}{v_F}.
\qquad\text{(R21-10)}
$$

当 $ta=1$ 时，右边就是 $2\pi/v\_F$。因此：

> **归一化因子不是自由常数**。它是连续应力密度与格点 hopping 的换算因子；在自然格点归一化下已经固定为 $2\pi/v\_F$。

### 命题 R21.1（连续极限的归一化校正）【条件证成】

在下列具名输入下：

1. 半满自由费米子链的识别 `G77`；
2. 外部定理 Eisler–Tonni–Peschel 的连续极限 (R21-3)；
3. 格点色散 (R21-5)；

有

$$
h\simeq \frac{2\pi}{v_F}\,B_{\rm bare},
\qquad
B_{\rm bare}\simeq \frac{l}{a^2}\beta,
\qquad
v_F=2ta\sin(k_Fa).
\qquad\text{(R21-11)}
$$

特别地，$t=a=1,\ k\_Fa=\pi/2$ 给出

$$
h\simeq \pi\,\beta.
\qquad\text{(R21-12)}
$$

**边界**：这不是说 R20 的有限 $N$ 数据已经精确等于 $\pi$。它说明正确的连续主符号常数是 $\pi$；R22 进一步证明 R20 的 $3.3$–$3.4$ 不是 $\pi$ 的有效数值估计器。

---

## §3 重新读 R20 的“1.85”

R20 的外推范围为

$$
A_\infty^{\rm LS}\approx3.3048,
\qquad
\rho_\infty\approx3.3802,
\qquad
|\lambda_\infty|\approx3.3943 .
\qquad\text{(R21-13)}
$$

R20 把它们与 $2\pi$ 相比，得到

$$
\frac{2\pi}{|\lambda_\infty|}\approx1.851.
\qquad\text{(R21-14)}
$$

但按 R21.1，应与 $\pi$ 相比。于是

$$
\frac{2\pi}{|\lambda_\infty|}
=
2\cdot\frac{\pi}{|\lambda_\infty|}
\approx 1.851,
\qquad
\frac{\pi}{|\lambda_\infty|}\approx0.9255.
\qquad\text{(R21-15)}
$$

所以“$1.85$ 倍”不是新的物理常数，而可作如下算术分解：

$$
\text{因子 }2\text{（来自 }v_F=2\text{）}
\ \times\
\text{商 }0.9255.
\qquad\text{(R21-16)}
$$

**校正**：这里的 $0.9255$ 不能被解释为已证实的有限尺寸残差。R20 的分母用 $B\_N^{\rm R20}$，而外部定理用 (R21-3) 的 $T\_N$；两者的近零谱子空间不同。R22 的独立复核显示，沿 $T\_N$ 的近零模，$N=24$ 的比值已经是 $-3.14033$。因此 $0.9255$ 只是 R20 估计器的比值商，不是物理修正因子。

### 命题 R21.2（R20 常数缺口的重述）【条件判定】

在 R20 的自由费米子支线上：

1. **原 A1 的裸 $2\pi$ 目标不成立**，因为它没有包含 $v\_F$；
2. **正确连续目标是 $2\pi/v\_F=\pi$**；
3. **R20 的 $1.85$ 倍重述为 $2\times0.9255$**，其中 $2$ 是速度因子，$0.9255$ 是 R20 估计器／测试子空间的比值商，不是已证残差；
4. 因此 `Z-STRESS-2π` 不再是一个“需要从 Zero 生出 $1.85$”的开放常数问题，而是“把 $v\_F$ 与应力归一化纳入 $B\_B$ 的定义”的口径问题。

**边界**：R21 只解释归一化；它没有证明 R20 的有限 $N$ 外推已经精确等于 $\pi$，也没有把 `Z-STRESS-2π` 抬高成四维 L1 定理。R22 修正了 R20 的估计器解释：低能／主符号层支持 $\pi$，任意光滑二次型比值不支持。

---

## §4 对 L1 的净影响

### 4.1 应改写的目标

原目标写作

$$
K_B\longrightarrow 2\pi B_B ,
\qquad\text{(R21-17)}
$$

这只在 $v\_F=1$ 的共形归一化下可直接成立。对 R20 的格点模型，应把裸几何算子与物理 boost 生成元分开：

$$
K_B\longrightarrow \frac{2\pi}{v_F}B_{\rm bare}
\quad\text{或}\quad
K_B\longrightarrow 2\pi B_B^{\rm phys},
\qquad
B_B^{\rm phys}:=\frac{1}{v_F}B_{\rm bare}.
\qquad\text{(R21-18)}
$$

这条改写不新增参数：$v\_F$ 已由 `G77`／`R8-L2` 登记，半满模型中为 $2$。

### 4.2 什么都没关闭

R21 只处理常数口径。L1 仍有以下承重项：

| 缺口 | 当前状态 | 本轮影响 |
|:--|:--|:--|
| `Z-WICK` | 开放 | 不变 |
| `Z-WEDGE` | 开放 | 不变 |
| `Z-CORE`／`Z-TAIL` | 开放 | 不变 |
| `Z-CRIT-DER`／自由费米识别 | 具名输入／识别 | 不变 |
| `Z-STRESS-2π` | **符号层口径已校正**；算子层与四维应力接口仍开放 | 不再把 $1.85$ 或 $0.9255$ 当独立常数 |

因此：

$$

\text{R21 解决了 R20 的常数读数，不解决 L1 的几何／解析延拓问题。}

\qquad\text{(R21-19)}
$$

---

## §5 与 R19／R20 的状态关系

| 文档 | 更新 |
|:--|:--|
| [`R20`](R20_A1_verdict_shape_holds_constant_fails.md) | 形状结论保留；“常数差 $1.85$ 倍”改读为“相对错误目标 $2\pi$ 的读数”；正确连续目标是 $\pi=2\pi/v\_F$ |
| [`R19`](R19_L1_upstream_probability_phase_and_missing_boost.md) | A1 的比较式需将 $B\_B$ 的物理应力／速度归一化写明；旋转双覆盖 no-go 不变 |
| [`R15`](R15_zcar_double_cover_and_zstress_scale.md) | `Z-STRESS` 的常数从估计器依赖读数收敛到 $2\pi/v\_F$ 口径，而不是独立拟合值 |
| [`R22`](R22_principal_symbol_vs_r20_estimator.md) | 保留本轮的 $\pi$；排除 R20 的 $8\%$ 有限尺寸解释，并把符号层与算子层分开 |
| [`R0`](R0_publication_theorem.md) | 外部／方向审计范围扩到 `R22` |
| [`STATUS.md`](STATUS.md) | 新增 §2.21 |

**范围边界**：R21 的比较对象是 1+1D 半满自由费米子区间。它不构造四维楔形，也不替代共形网、反射正性或低维算符控制。

---

## §6 可证伪点与下一步

### 可证伪点

若使用外部定理的 $T\_N$ 定义、限制到其近零谱子空间，并按一阶主符号读取系数仍不给出 $\pi$，则 R21.1 的口径判定应撤回。R22 已给出该检验的首个数值支持；任意光滑二次型不收敛到 $\pi$ 本身不构成本命题的反例。

### 下一步

不再寻找“补上 $1.85$”的机制。主线应记为：

1. **1+1D 归一化**：R22 已把符号层与算子层分开；后续只需在近零谱子空间上证明带端点和尾项误差的条件定理；
2. **四维接口**：在 `Z-WEDGE` 中说明 $v\_F$ 重标度与局部 boost 法向如何提升；
3. **其余门槛**：`Z-WICK`、`Z-CORE`／`Z-TAIL`、`Z-CRIT-DER` 不变。

**L1 状态**：$K\_B\to 2\pi B\_B$ 仍未关闭；R20 的 $1.85$ 倍常数缺口已降为 $2\pi/v\_F$ 的归一化口径，且其 $0.9255$ 读数由 R22 排除为估计器偏差。
