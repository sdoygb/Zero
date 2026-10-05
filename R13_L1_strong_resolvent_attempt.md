# R13 · L1 正路由：从 Zero 临界链到强预解收敛的诚实状态

**日期**：2026-10-02  
**性质**：把 R12 §2.3 的“临界自由费米子正路由”从一条愿景写成可审计的命题，并给出本轮的实际判决。  
**政策**：区分【已证】／【条件证成】／【识别】／【开放】／【排除】；不得把“数值像 BW”“重标度后可能收敛”改写成“已从 Zero 证明强预解收敛”。  
**依赖**：[`Z0`](Z0_zero_never_rests_single_axiom.md)、[`G27`](G27_purification_attempt.md)、[`G40`](G40_metric_from_closed_walk_counting.md)、[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md)、[`G75`](G75_quantum_geometry_modular_readout.md)、[`G77`](G77_staggered_coupling_from_A5.md)、[`G79`](G79_horizon_thermodynamics.md)、[`R8_L1_refutation_attempt.md`](R8_L1_refutation_attempt.md)、[`R12`](R12_zero_native_gap_filling.md)。  
**配套**：[`R13_external_limit_lemmas.md`](R13_external_limit_lemmas.md)（外部严格结果）、[`R13_refutation_attempt.md`](R13_refutation_attempt.md)（对抗审计）、[`R13_numeric_probe.py`](R13_numeric_probe.py)（数值探针）、[`R13_check.py`](R13_check.py)（独立核验）。

$$
\boxed{
\begin{aligned}
&\text{本轮没有证明 L1。}\\
&\text{R12 §2.3 的原样命题（未重标度的 }h_N\text{ 强预解收敛到有限 }2\pi B_B\text{）被排除；}\\
&\text{存活下来的是一条具名的条件 }1+1D\text{ 路线，其四个承重输入均为【识别／开放】。}
\end{aligned}}
$$

> **一句话**：Zero 底层能给出一条临界的闭环计数结构，但“它是一条半满临界自由费米链”是**识别**而不是导出（`G77` 已自认）。在此识别之上，有限维模 Hamiltonian 有精确公式 $h_N=\log((1-C_N)/C_N)$，数值上确实给出棋盘结构、线性增长的近邻权重与接近抛物线的主导剖面；但原样 $h_N$ 的谱半径随 $N$ 线性发散，固定谱参数的预解式趋于零而不是 $(z-2\pi B_B)^{-1}$。要把路线变成定理，必须先补 `Z-CRIT-DER`、`Z-SCALE`、`Z-HILB/CORE`、`Z-STRESS`、`Z-CONF`。

---

## §0 判决摘要

| 命题 | 判定 | 依据 |
|:--|:--|:--|
| 有限维准自由态有 $h=\log((1-C)/C)$ | **已证** | §2.1，引理 R13.1；外部条目 R13.1／R13.3 |
| 区间的 $h_N$ 有棋盘（粒子–空穴）结构 | **已证（有限维）** | §2.2，引理 R13.2；数值探针 F1 |
| 近邻权重中心值随 $N$ 线性增长，斜率固定 | **数值证成（有限 $N$）** | §2.3；探针 $0.8691\,N-0.7814$，$R^2=0.99996$ |
| 未重标度 $h_N$ 强预解收敛到有限 $2\pi B_B$ | **排除** | §3，命题 R13.4、反例 R13.5；谱半径 $\sim1.72N$ |
| 逐矩阵元意义下收敛到局域 $2\pi\int\beta T_{00}$ | **排除** | §3.3；长程尾项与端点层 |
| 二次型／分布意义下在条件输入下收敛 | **条件证成（1+1D）** | §4，定理 R13.6；依赖外部的 Eisler–Tonni–Peschel 与 Kato |
| 由 Zero 原生对象导出自由费米结构与半满 | **排除（作为导出）／识别** | §1.3；`G77` §4 自认 |
| $2\pi$、$T_{00}$ 格点归一化由 Zero 导出 | **排除／识别** | §4.3；`G57`／`G75`／`G79` |
| 一维区间经薄饼极限关闭四维 L1 | **排除** | §4.3；横向 boost 与曲率需 `Z-CONF` |
| `G79` 的“抛物线相关 $0.865$”是稳定误差估计 | **排除** | §5.1；改为截断依赖量，实测 $0.865\to0.964$ |

**净判决**：R12 §2.3 的“唯一可定理化正路由”应降级为**候选条件路线**；它的原样强预解版本已被谱发散排除。L1 保持**开放**。

---

## §1 Z_crit 的构造链：哪几步是导出，哪几步是识别

### 1.1 从 Z0 到闭环计数结构【导出／已有引理】

1. Z0（零不断乱动）导出词、图与全分支计数（[`Z0`](Z0_zero_never_rests_single_axiom.md) §0–§2）。
2. 闭合词自带循环次序，给出环图与站点 $i\in\{0,\dots,L-1\}$（[`G27`](G27_purification_attempt.md)）。
3. 闭环游程计数给出边权 $W_e$ 与图 Laplacian $L_W$（[`G40`](G40_metric_from_closed_walk_counting.md)、[`G46`](G46_k_is_the_lifetime.md)）。
4. 均匀权下该结构是**临界**的：一维纠缠熵按 $\tfrac{c}{3}\ln\ell$ 标度，$c=1$（[`G75`](G75_quantum_geometry_modular_readout.md) §2 的标定与实算）。

到这一步为止用到的都是已有导出，问题只在“临界”二字：G75 的检验是**有限环数值**，不是极限定理。

### 1.2 从闭环计数结构到半满临界链【识别】

把上述环图识别为**一维紧束缚链**：

$$
H=-t\sum_j\bigl(c_j^\dagger c_{j+1}+\mathrm{h.c.}\bigr),
\qquad t=1,
\qquad \text{半满},
$$

其中 $c_j$ 是 CAR 费米算符。这个识别买回了：有限费米速度、单费米点结构、可直接计算的关联矩阵与模 Hamiltonian。

### 1.3 这一步为什么不是导出【排除作为导出】

[`G77`](G77_staggered_coupling_from_A5.md) §4 已明写：

> “把 Z3 的汇写成自由费米模型里的在格能量；这是标准的对应，但**未**从 Z3 严格推出自由费米形式。”

[`G75`](G75_quantum_geometry_modular_readout.md) §7 与 [`G76`](G76_area_law_in_2d.md) §5 同样把交错质量／耦合形式登记为“识别”。因此：

$$
\boxed{
\text{Z0 条款（A0–A5 历史命名）没有 CAR、反对易关系、hopping、半满占据或单费米点选择器；}
\text{“临界自由费米链”是模型识别，不是 Zero 定理。}
}
$$

保持 Z3 字面数据不变，至少还有经典两态链、临界横场 Ising 链（$c=1/2$）、硬核玻色链等互不相容的完成方式。要把它变成定理，需要新增 `Z-CRIT-DER`（见 §4.3）。

---

## §2 已证的有限维层

### 2.1 引理 R13.1（准自由模 Hamiltonian）【已证】

设有限维 CAR 系统处于规范不变准自由态，单粒子关联矩阵

$$
C_{ij}=\omega(c_i^\dagger c_j),
\qquad 0<C<1 .
$$

则存在单粒子厄米算子

$$
h=\log\!\left(\frac{1-C}{C}\right),
\qquad
C=\frac{1}{e^{h}+1},
$$

使得模 Hamiltonian 是二次量子化 $\mathrm d\Gamma(h)$。**证明**：Bogoliubov 变换同时对角化 $C$ 与二次 Hamiltonian；每个模态的占据数必为 Fermi–Dirac 分布 $\nu_i=(e^{\varepsilon_i}+1)^{-1}$，反解即得 $\varepsilon_i=\log((1-\nu_i)/\nu_i)$。$\square$

**意义**：L1 的有限维侧是精确算子对象，不是近似。**不能**由此得到 $N\to\infty$ 的极限：$C_N$ 的谱在临界情形充满 $(0,1)$，$\log$ 把指数小的本征值放大为对数发散。

### 2.2 引理 R13.2（棋盘与粒子–空穴结构）【已证（有限维）】

半满链基态满足粒子–空穴对称，$\operatorname{spec}(C_N)=1-\operatorname{spec}(C_N)$，故 $\operatorname{spec}(h_N)=-\operatorname{spec}(h_N)$。$C_N$ 的符号在每个费米点跳跃，使 $C_N$ 只在奇数距离上有矩阵元；任意矩阵函数保持该分级，所以

$$
(h_N)_{ij}=0\quad\text{当 }i-j\text{ 为偶数},
\qquad
(h_N)_{ij}\ne0\quad\text{当 }i-j\text{ 为奇数}.
$$

高精度数值（探针 F1，`mpmath` 80 位）给出偶数距离最大伪矩阵元 $5.2\times10^{-69}$，奇数距离最小非零元 $1.6\times10^{-10}$。

### 2.3 命题 R13.3（有限尺寸的定量结构）【数值证成】

对无限链上区间 $0,\dots,N-1$（格距 $a$、$l=Na$）：

1. 中心最近邻权重线性增长：$|h_{c,c+1}|=0.8691\,N-0.7814$，$R^2=0.99996$（探针 F2）；
2. 谱半径线性增长：$\lVert h_N\rVert_{\rm spec}=1.7195\,N-2.754$，相关系数 $0.99999$（探针 F6）；
3. 近邻权重沿位置的归一化剖面与抛物线 $\beta(x)=x(l-x)/l$ 高度相关（探针 F3；注意端点层）；
4. 奇数距离尾项按固定格距快速衰减：$|h_{3}|/|h_{1}|\to$ 约 $0.025$，$|h_5|/|h_1|$ 更小（探针 F7）。

这四条只描述**有限 $N$ 的数值结构**；它们不构成极限定理。

---

## §3 原样强预解命题及其失败

### 3.1 命题 R13.4（R12 §2.3 的原文）

> 存在公共 Hilbert 空间与嵌入，使未重标度的 $h_N=\log((1-C_N)/C_N)$ 在 $N\to\infty$ 时强预解收敛到 $2\pi B_B$，从而在 $1+1D$ 关闭 L1。

### 3.2 反例 R13.5（谱发散使固定谱参数预解式失效）【已证】

由 §2.3，$h_N$ 的谱半径以正斜率线性发散。取最小模型

$$
H_N=N H_0,\qquad H_0=I,
$$

则对任意固定非实 $z$，

$$
(z-H_N)^{-1}=\frac{1}{z-N}I\longrightarrow 0,
$$

而零算子不是任何有限自伴算子的预解式。即使 $H_N=N K_N+o(N)$ 且 $K_N\to K_\infty\ne0$，固定 $z$ 的预解式也不等于 $(z-K_\infty)^{-1}$。**结论**：命题 R13.4 按原文为假／未定义；必须先做仿射重整化

$$
h_N^{\rm ren}=a_Nh_N+b_NI,
$$

并把谱参数或低能子空间同时重整化。$\square$

### 3.3 为什么“逐矩阵元”是错的拓扑【排除】

$h_N$ 含长程尾项；逐矩阵元收敛不能控制远距离项的总二次型贡献。最小反例：$H_N=H_0+N^{-1}(J_N-I_N)$，每个固定距离矩阵元 $\to0$，但对归一化常向量 $\langle u_N,N^{-1}(J_N-I_N)u_N\rangle\to1$。

正确的判据只能是**二次型／分布意义**：

$$
\left\langle f,\left(h_N^{\rm ren}-B_B\right)g\right\rangle\to0,
\qquad f,g\in C_c^\infty(0,l),
$$

外加端点和尾项的 uniform 控制。逐矩阵元比较（以及用单个数“相关系数”作证据）都应排除。

---

## §4 修复后的条件定理与缺口清单

### 4.1 定理 R13.6（条件，1+1D）　【条件证成】

**假设**（四条具名输入）：

1. `Z-CRIT-DER`：CAR、hopping、半满、单费米点已由 Zero 或明确输入给出；
2. `Z-SCALE`：格距 $a\to0$、$l=Na$ 固定，且已选定重标度 $a_N,b_N$ 与低能窗口；
3. `Z-HILB/CORE`：公共 Hilbert 空间、细化嵌入 $E_a$、公共核心 $\mathcal D$；
4. `Z-STRESS`：连续 $T_{00}$ 的正规化与 $\beta(x)=x(l-x)/l$ 的几何角标定。

**结论**：存在由 `Z-STRESS` 固定的常数 $\kappa$，使对 $f,g\in C_c^\infty(0,l)$，

$$
\left\langle E_af,\ h_N^{\rm ren}\,E_ag\right\rangle
\longrightarrow
\kappa\int_0^l dx\,\beta(x)\,\tau_{fg}(x),
$$

其中 $\tau_{fg}$ 是自由费米场的局部能量密度二次型；进而关联自伴算子强预解收敛，由 Trotter–Kato 给出模流在紧时间区间上的强收敛。**这是 $1+1D$ 的条件定理，不是 Zero 定理，也不关闭四维 L1。**

### 4.2 定理所依赖的三件外部工具【引用】

1. **Eisler–Tonni–Peschel**（arXiv:1902.04474）：无限自由费米链的连续极限把模 Hamiltonian 写成
   $$
   \mathcal H=2\pi l\int_0^l dx\,\frac{x}{l}\left(1-\frac{x}{l}\right)T_{00}(x),
   $$
   前提是**把全部长程跳跃纳入连续极限**；只比较最近邻项会有约 $8\%$ 的中心偏差（Eisler–Peschel, arXiv:1703.08126）。这解释了本项目实测的近邻剖面偏差约 $3\%$ 与尾项约 $2.5\%$：它们是同一非局域结构的两个投影，必须合并处理。
2. **Kato 二次型收敛定理**：闭、稠定、下半有界二次型的 Kato／Mosco 收敛直接给出强预解收敛（Kato, *Perturbation Theory*, Ch. VI；Reed–Simon I, VIII.7）。这是把“二次型收敛”升级为“模流收敛”的机器。

> **后续校正（[`R22`](R22_principal_symbol_vs_r20_estimator.md) §6）**：一阶目标算子谱为全实轴，**不是**下半有界二次型。后续应使用图二次型的 Mosco 收敛，或直接证明预解一致性；不能把这里的 Kato 工具无条件套到 $A_N\to A_\infty$。ETP 级数只给出统一界与余项因子分解，真正未关闭的是 `R22-STAG`／`R22-UV`。
3. **Brunetti–Guido–Longo**：共形协变因果可加局部代数网上，双锥模流等于区域保持共形流（arXiv:funct-an/9302008）。它是 `Z-CONF` 的模板，但前提是已有连续共形局部网。

### 4.3 缺口清单（原样继承并改写 R12 的 `Z-*` 标签）

| 标签 | 还缺什么 | 关掉哪一条 |
|:--|:--|:--|
| `Z-CRIT-DER` | 从 Z0 条款（A0–A5 历史命名）导出 CAR、hopping、半满与单费米点 | §1.3 的识别／循环 |
| `Z-SCALE` | 固定 $a_N$、$v_F$、重标度与低能窗口，给谱宽度界 | §3.2 的谱发散 |
| `Z-HILB` | 公共 Hilbert 空间、细化嵌入与连续局部代数网 | §3.3 的拓扑 |
| `Z-CORE` | 公共核心与二次型／强预解收敛，含高模消失 | §4.1 的结论 |
| `Z-TAIL` | 长程尾项的 uniform 二次型界与端点层控制 | §3.3 |
| `Z-STRESS` | $T_{00}$ 的正规化与 $2\pi$ 几何角识别 | §0 的归一化行 |
| `Z-CONF` | 一维区间到四维双锥的横向模式、曲率与 boost 交换子 | 四维 L1 |

> **后续归一化校正（[`R21`](R21_vf_normalization_resolves_the_r20_factor.md)）**：§4.2 的连续公式 $2\pi l\int\beta T_{00}$ 中的 $T_{00}$ 已带连续应力归一化，不能把该 $2\pi$ 直接乘到裸格点 hopping 上再与同一量比较。外部定理给出连续系数 $2\pi/v_F$；本模型 $v_F=2$，故裸格点比较的连续目标是 $\pi$。因此 `Z-STRESS` 不应再寻找独立的 $1.85$ 常数，而应固定应力／速度口径。

---

## §5 数值证据与新发现

### 5.1 `G79` 的“抛物线相关 $0.865$”是截断依赖量【更正】

[`G79`](G79_horizon_thermodynamics.md) §1.1 报出近邻键系数与抛物线 $x(\ell-x)$ 的相关 $0.865$。该数字来自 `G79_check.py` 中对关联矩阵本征值的 $\mathrm{clip}(10^{-9},1-10^{-9})$ 正则化。对同一配置（$N=96$ 周期环、区间 $[32,64)$）改变截断：

| 截断 $c$ | 中心键 $\lvert h_{c,c+1}\rvert$ | 相关 |
|--:|--:|--:|
| $10^{-9}$ | $12.888$ | $0.8654$ |
| $10^{-12}$ | $16.876$ | $0.9273$ |
| $10^{-14}$ | $19.313$ | $0.9544$ |
| $10^{-15}$ | $20.571$ | $0.9642$ |

真值（$\lambda_{\min}(C_A)\approx8.4\times10^{-18}$，在双精度下无法直接取）比 $0.964$ 更接近 $1$。**结论**：$0.865$ 不是误差估计，而是正则化与端点效应共同造成的、偏保守的数字；把它当精确核证据或当反证据都不成立。这一条同时更正 `R8_L1`／`R12` 中对 $0.865$ 的引用。

### 5.2 二次型收敛探针

[`R13_numeric_probe.py`](R13_numeric_probe.py) 用最近邻 BW 算子 $B_N$ 与中心键归一化 $\alpha_N$ 做探针：四条光滑测试序列的最大相对差从 $0.0710$（$N=8$）单调降到 $0.0211$（$N=20$），$\alpha_N\to$ 约 $3.32$。这支持 §4.1 的二次型收敛**方向**，但不替代 uniformity 与端点控制。

---

## §6 净推进与未推进

**净推进（可计入）**

1. 把 L1 正路由从“愿景”写成可证伪命题，并给出其**原样版本的反例**（谱发散）。
2. 确认有限维层的两条精确引理（准自由公式、棋盘结构）并给出高精度数值结构。
3. 给出修复版的条件 $1+1D$ 定理与七个具名缺口，接上外部三条定理。
4. 更正确认 `G79` 的 $0.865$ 为截断依赖量（实测 $0.865\to0.964$），把它从承重证据降为一致性检查。

**未推进（继续开放）**

1. 没有从 Zero 导出自由费米结构（`Z-CRIT-DER`）。
2. 没有证明重标度后二次型的 uniform 收敛与端点／尾项控制（`Z-SCALE`／`Z-TAIL`／`Z-CORE`）。
3. 没有固定 $2\pi$ 与 $T_{00}$ 的归一化（`Z-STRESS`）。
4. 没有四维横向 boost、曲率与薄饼极限（`Z-CONF`）。
5. 因此 **J1 仍未关闭，条件恢复判定不变**。

$$
\boxed{
\begin{aligned}
&\text{Current status: L1 open; conditional 1+1D program only.}\\
&\text{The original strong-resolvent claim is excluded, not proved.}
\end{aligned}}
$$

---

## §7 核验

```text
python3 R13_check.py
python3 R13_numeric_probe.py
```
