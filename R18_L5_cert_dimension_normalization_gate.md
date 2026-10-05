# R18 · L5-CERT 第一关：低维通道的目录与归一化判据

**日期**：2026-10-02  
**性质**：L5 的第一次定向攻坚。新增两条已证判据与一条审计判决；不关闭 L1，不新增物理参数，不新增表示层标签。  
**唯一目标**：Jacobson 2016 固定体积首阶平衡的上游门槛 L5，即
$$
D(\sigma_R\|\rho)=o(R^d).
$$
**依赖**：[`R8`](R8_jacobson_entanglement_equilibrium_completion.md)、[`R9`](R9_external_GR_derivations_landscape.md)、[`R12`](R12_zero_native_gap_filling.md)、[`R17`](R17_L1_critical_path_and_L5_gate.md)、[`G76`](G76_area_law_in_2d.md)、[`G77`](G77_staggered_coupling_from_A5.md)、[`Z16`](Z16_zunif_balanced_regular_module.md)。  
**核验**：[`R18_check.py`](R18_check.py)。

$$
\boxed{
\begin{aligned}
&\text{L5 的危险不是一个系数的取值问题，而是两个因子：}\\
&\text{（目录）低维通道子空间 }V_{\rm light}\text{ 的维数；}\\
&\text{（归一化）物理约束映射 }\mathcal N\text{ 在 }V_{\rm light}\text{ 上的单射性。}\\
&\text{固定体积与固定 }T_{00}\text{ 至多给 }m\le2\text{ 条约束，}\\
&\text{故 }\dim V_{\rm light}\ge3\text{ 时 Z-STRESS 必败；}\\
&\text{L5 只剩目录空、对称零化、接触项合并三条路。}
\end{aligned}}
$$

> **一句话**：上一轮把 L5 立成 go/no-go，但没说清“关掉一个系数”还是“清空一个目录”才算关闭。本轮把它算清：危险等于低维通道维数超过物理约束数。因此 Z-STRESS 的归一化不是 L5 的解；L5 只能由空目录、对称零化或接触项合并来关。

---

## §0 判决摘要

| 问题 | 本轮判定 | 依据 |
|:--|:--|:--|
| L5 与 L1 是否互相推出 | **否；两者逻辑独立** | §1 |
| 危险性的精确来源 | 低维通道 $X\ne0$ 与自由振幅 $\varepsilon_R=R^{d/2}$ 同时存在 | §2，(R18-2) |
| L5 的关闭判据 | $\mathcal N$ 在 $V_{\rm light}$ 上单射 | §3，命题 R18.1 |
| Z-STRESS 是否足够 | **否**；$\dim V_{\rm light}\ge3$ 时必败 | §4，推论 R18.2 |
| 关闭 L5 的可用路径 | DIM-CAT／SYM-ZERO／CONTACT 三条 | §5，推论 R18.3 |
| 下一步 | 在 1+1D 支线枚举 $V_{\rm light}$ 并算维数 | §6 |
| L1 当前状态 | **仍未关闭** | 全篇 |

---

## §1 必要性重审：L5 不是 L1 的推论

L1 与 L5 的对象不同：

1. **L1**：参考态（真空）的模 Hamiltonian 是否收敛到几何 boost，
$$
K_B\longrightarrow 2\pi B_B .
$$
2. **L5**：一阶邻域中被激发的低维局域标量算符，是否污染固定体积首阶平衡。

`R12` 已在同一有限维条件类内给出两者的独立 no-go：L5 的 no-go 来自“模能一阶项为零而相对熵仍达 $R^d$ 阶”（命题 R12.2）；L1 的 no-go 来自“中央年龄剖面只给常数剖面”与“半侧模平移只生成交换子代数”（命题 R12.4、R12.5）。两条论证互不引用彼此的结论。

因此
$$
\text{L1}\not\Rightarrow\text{L5},\qquad \text{L5}\not\Rightarrow\text{L1}.
$$

**审计结论**：只关 L1 不足以关闭 Jacobson 2016；反过来，L5 失败也不会使 L1 的数学内容失效。L5 是同一外部锚的第二前提，而不是 L1 的子步骤。按 `R17-STOP`，L5 只需被收成一条已证判据或一个具名输入，不需要继续新增表示层标签。

---

## §2 R12.2 反例的精确结构

`R12` 的三权重反例可以拆成两件互相独立的事：

1. **通道**：存在非零扰动方向 $X$，其相对熵二阶系数严格为正；
2. **振幅**：扰动振幅 $\varepsilon_R$ 没有被任何物理约束钉死，因而可以取成 $\varepsilon_R=R^{d/2}$。

在 $\sigma_R=\rho+\varepsilon_R X$、$X$ 零迹、$\varepsilon_R\to0$ 时，有限维精确展开为
$$
D(\sigma_R\|\rho)=\frac{\varepsilon_R^{2}}{2}\chi_K(X)+O(\varepsilon_R^{3}),
\qquad \chi_K(X)>0 .
\qquad\text{(R18-1)}
$$
取 $\varepsilon_R=R^{d/2}$ 即得
$$
D(\sigma_R\|\rho)=\frac12R^{d}\chi_K(X)+O(R^{3d/2})=\Theta(R^{d})\ne o(R^{d}).
\qquad\text{(R18-2)}
$$
`R12.2` 的 $w=(1,2,4)/7$、$x=(1,-2,1)$ 正是 (R18-1) 的一次取值：$\chi_K(x)=\sum_a x_a^{2}/w_a=91/4$，故系数为 $91/8$。

**关键读法**：危险来自 (R18-2) 中振幅的**自由**，而不是某个算符的特定系数。把振幅钉死，就等于关掉这条通道。这给出本轮要计算的两个量：低维通道的维数，以及物理约束在该通道上的秩。

---

## §3 命题 R18.1（目录与归一化双因子判据）【已证】

设候选态类 $\mathcal C$ 在参考态 $\rho$ 处的实切空间为 $V$，并给定候选连续极限的维度赋值 $\dim:V\to[0,\infty)$。定义**低维通道子空间**
$$
V_{\rm light}:=\{X\in V:\dim X\le d/2\}.
\qquad\text{(R18-3)}
$$
设 $\mathcal N=(N_1,\dots,N_m)$ 为一组被物理钉住的可观测量（体积、能量密度、计数等），记
$$
\mathcal N(X):=\bigl(\operatorname{Tr}(XN_1),\dots,\operatorname{Tr}(XN_m)\bigr).
\qquad\text{(R18-4)}
$$

**(i)（破坏方向）** 若存在 $X\in\ker \mathcal N|_{V_{\rm light}}$、$X\ne0$ 且 $\chi_K(X)>0$，则态族 $\sigma_R:=\rho+R^{d/2}X$ 满足
$$
\mathcal N(\sigma_R-\rho)=0,
\qquad
D(\sigma_R\|\rho)=\frac12R^{d}\chi_K(X)+O(R^{3d/2})=\Theta(R^{d}),
\qquad\text{(R18-5)}
$$
故 $D\ne o(R^d)$，L5 失败。

**(ii)（关闭方向）** 若 $\mathcal N|_{V_{\rm light}}$ 单射，则任何满足 $\mathcal N(\sigma_R-\rho)=O(1)$ 的族都有 $\varepsilon_R=O(1)$，从而 $D(\sigma_R\|\rho)=O(1)=o(R^d)$。

**证明**  
(i) 由 $\mathcal N(X)=0$ 得 (R18-5) 第一式。相对熵展开用 (R18-1)；代入 $\varepsilon_R=R^{d/2}$ 得第二式，且 $\chi_K(X)>0$ 使主项严格为正。$\mathcal N$ 线性，不改变主项阶数。  
(ii) $V_{\rm light}$ 有限维且 $\mathcal N$ 在其上单射，故存在 $c>0$ 使 $\|\mathcal N(X)\|\ge c\|X\|$ 对一切 $X\in V_{\rm light}$ 成立。取 $\|X\|=1$；由 $\mathcal N(\sigma_R-\rho)=O(1)$ 得 $\varepsilon_R=O(1)$，代入 (R18-1) 得 $D=O(1)$。$\square$

**边界**：命题 (i) 只给充分条件。即使 $\mathcal N|_{V_{\rm light}}$ 非单射，也可能因低维系数在候选类上恒为零（目录空）而使该通道整体消失；这一点由 §5 处理。

**备注**：有限维 BKM 二次型在零迹厄米扰动空间上正定，故 $\chi_K(X)>0$ 对一切非零零迹 $X$ 自动成立。命题保留该条件，只为把陈述写成不依赖规范的有限维形式；$V$ 若已被额外物理条件（如模零模、协变性）限制，只需在受限的 $V_{\rm light}$ 上重新应用本判据。

---

## §4 推论 R18.2（Z-STRESS 的约束数上限）【已证】

固定体积给出 $1$ 条标量约束（$\delta V=0$）；固定 $T_{00}$ 的期望给出 $1$ 条（$\delta\langle T_{00}\rangle$ 固定）。因此 `Z-STRESS` 提供的独立约束数 $m\le2$。

由命题 R18.1(ii)，$\mathcal N|_{V_{\rm light}}$ 单射要求 $m\ge\dim V_{\rm light}$。于是：

| 情形 | 判定 |
|:--|:--|
| $\dim V_{\rm light}\ge3$ | $\ker\mathcal N|_{V_{\rm light}}\ne\{0\}$ 必然成立，`Z-STRESS` **不能**关闭 L5 |
| $\dim V_{\rm light}=2$ | 需逐例检验 $(\delta V,\delta T_{00})$ 在 $V_{\rm light}$ 上是否单射；任一低维通道对 $T_{00}$ 中性即失败 |
| $\dim V_{\rm light}\le1$ | `Z-STRESS` 单独**可能**够（是充分性检验，不是推论） |

**推论**：把 L5 当作“给 $T_{00}$ 挑一个好归一化”的问题，在一般情形下消除不掉它；只有低维通道维数很小时，`Z-STRESS` 才可能足够。

**与 R12 §1.3 的关系**：该处的第 2 条（“一致 BKM 界”）在本判据里被收紧为一条可计算条件——**归一化的独立约束数必须不少于低维通道维数**。

---

## §5 推论 R18.3（L5 的唯一可关路径）【审计结论】

由 R18.1、R18.2，L5 只能由下列之一关闭：

1. **DIM-CAT（目录空）**：所有 $\dim\le d/2$ 的规范中性局域标量算符，在候选态类上的系数恒为零；
2. **SYM-ZERO（对称零化）**：存在附加对称性，使每个低维通道的系数严格为零；
3. **CONTACT（接触项合并）**：低维贡献与接触项／局部反项统一合并，仍给 $D=o(R^d)$。

“靠物理归一化把振幅钉死”不构成第四条路，除非归一化的独立约束数不少于 $\dim V_{\rm light}$。

这与 [`R17`](R17_L1_critical_path_and_L5_gate.md) §5 的 `L5-CERT` 三条等价，但把其中第 2 条与 `R12` §1.3 的第 2 条合并为一条可计算判据。**因此 L5 不再是“担心某个危险算符”，而是“数出低维通道维数，再选一条零化或合并机制”。**

---

## §6 判决与下一步

**判决**：L5 的 go/no-go 现在是一个良定的二元问题——先量目录，再选对称。它不再需要新的表示层，也不需要新的物理公理。

**下一步（唯一）**：在 1+1D 自由费米子支线上枚举 $V_{\rm light}$：

1. 取 `R13`／`R17` 已固定的嵌入与重标度；
2. 列出该支线上所有规范中性、局域、标量的二次费米型通道及其维度；
3. 计算 $\dim V_{\rm light}$ 与 $(\delta V,\delta T_{00})$ 在该子空间上的秩；
4. 若 $\dim V_{\rm light}\ge3$ 且无 `SYM-ZERO`，按 `R17-STOP` 放弃 Jacobson 2016 首阶平衡锚，或改写平衡条件，而不是继续堆 `Z-*`。

**不做的事**：

1. 不在 $\dim V_{\rm light}$ 未确定前继续细分 `Z-CAR`；
2. 不把 (R18-5) 写成对四维 L5 的判决；它是结构判据，不是四维反例；
3. 不把 `Z-STRESS` 的常数 $\approx\pi^{2}/3$ 当作 L5 的关闭证据。

**L1 状态**：$K_B\to2\pi B_B$ 仍未关闭；本轮只改 L5 的判定方式。

**判定不变**：外部锚仍是**条件恢复**；`R18` 只把 L5 收成结构判据，不抬高任何已证内容，也不把候选写成定理。

---

## §7 对既有状态的影响

| 文档 | 更新 |
|:--|:--|
| [`R12`](R12_zero_native_gap_filling.md) | 其 §1.3 第 2 条的“一致 BKM 界”细化为“归一化约束数 $\ge\dim V_{\rm light}$” |
| [`R17`](R17_L1_critical_path_and_L5_gate.md) | 其 `L5-CERT` 三条保留；新增“Z-STRESS 不是第四条路”的判据 |
| [`STATUS.md`](STATUS.md) | 新增 §2.18；L5 状态由“开放”细化为“已收为双因子判据，待 1+1D 枚举” |
| [`R0`](R0_publication_theorem.md) | 外部／方向审计范围扩到 `R18` |

**范围边界**：本轮的 (R18-1)–(R18-5) 只使用有限维 BKM 正定性与线性代数，不引用任何四维连续几何；因此它们是结构判据，不是四维 L5 的关闭或反例。
