# R6 · H5 闭合：实际闭环权到字典多项式的一致误差界

**日期**：2026-10-02
**性质**：R1 的 H5 交付；把“实际闭环权 `A circ A^{m-1}` 到字典多项式 `P_{k,d}(c)` 的 `O(a²)` 数值旁证”升级为一个带显式假设、显式常数、且速率尖锐的定理。
**依赖**：[`R1`](R1_gamma_convergence_theorem.md) H5、[`Z7`](Z7_embedding_input_explicit_dictionary.md)（闭式字典）、[`G58`](G58_I2a_resolved_as_embedding_input.md)、[`G18`](G18_attackability_of_the_continuum_limit.md)（引理 67）、[`G41`](G41_lovelock_premises_under_nonuniform_weight.md)。
**等级标签**：【定义】/【引理】/【定理】/【尖锐速率】/【反例】/【数值核验】/【边界】/【结论】。
**核验**：[`R6_check.py`](R6_check.py)、独立复核 [`R6_independent_check.py`](R6_independent_check.py)、对抗审计 [`R6_refutation_attempt.md`](R6_refutation_attempt.md)。

$$
\
\text{在周期环面上、}k\text{ 固定、}c\text{ 为 Lipschitz 梯度时：}\quad
\left|\frac{\widehat W_{a,e}}{P_{k,d}(c_e)}-1\right|\le C\,a^{2} .
\
$$

> **一句话**：H5 **可以证明，不再是数值猜测**。O(a) 项由“关于边中点的反射＋游程反转”配对面**精确消去**；剩下的 O(a²) 由 Taylor 余项控制。而且速率是**尖锐**的：只有连续没有梯度时降为 O(a^{α})，间断时直接 O(1) 失败。

> **【后续状态｜[`R7`](R7_h3_h7_regularity.md)】** 本文 §7 的「H3／H7 仍属 R2」是 R6 写作时状态。随后 [`R7`](R7_h3_h7_regularity.md) 给出类到连续场的显式 mollification 映射，并在 R2 耦合＋I5b＋GDL 下条件闭合 H3；H7 修正为“极限 $Q=\mathcal P(c)$ 为 $C^2$、光滑提升在 $C^2$ 收敛”。当前状态以 [`STATUS.md`](STATUS.md) 为准。

---

## §0 结论（五句）

1. **H5 成立**：在周期立方环面 $\mathbb T^d$、每方向 $N$ 格、$N>2k$、$k$ 固定、$c\ge c\_{\min}>0$、且 $c\in C^{1,1}$（梯度 Lipschitz）时，存在只依赖 $d,k,c\_{\min},\|c\|\_{C^{1,1}}$ 的常数 $C$，使

   $$
   \max_{e\in E_a}\left|\frac{\widehat W_{a,e}}{P_{k,d}(c(m_e))}-1\right|\le C\,a^{2},
   \qquad a=1/N .
   $$

2. **证明的机关是一个精确对称**：一阶项 $\sum\_{\text{walks}}\sum\_j (m\_j-x\_e)$ 因反射—反转配对**恒为零**，因此误差从 O(a) 直接降到 O(a²)。
3. **速率尖锐**：$c\in C^{1,\alpha}$ 给 $O(a^{1+\alpha})$；$c\in C^{0,\alpha}$ 给 $O(a^{\alpha})$；$c$ 仅连续则仍 $\to0$；$c$ 间断则 $O(1)$ 失败。$C^{1,1}$ 恰好复现 H5 要的 O(a²)。
4. **代价是三条前提**：$k$ 必须固定（$k\propto N$ 时崩溃）；必须是周期／公共极限（开链边界边的游走计数不同）；$c$ 必须在同一个公共拓扑上连续。
5. **对 R1 的接续**：H5 关闭后，[`R1`](R1_gamma_convergence_theorem.md) 的**推论 R1-B**（实际闭环权满足 R1-A 的紧性、Γ-liminf、恢复列与唯一极限）不再是条件缺口，唯一的剩余桥是 H3/H7 的“旋转类到光滑场 $c$”构造，那属于 R2。

---

## §1 精确设定

取周期环面 $\mathbb T^d=(\mathbb R/\mathbb Z)^d$，细化尺度 $a=1/N$，顶点集

$$
\Lambda_a=a\mathbb Z^d/\mathbb Z^d,
\qquad
|\Lambda_a|=N^d .
$$

对每条沿坐标轴、长为 $a$ 的边 $e=\{x,y\}$，记其中点 $m\_e=(x+y)/2\bmod1$。加权邻接矩阵为

$$
A_{xy}=c(m_e)\quad(e=\{x,y\}),
\qquad
A_{xy}=0\quad\text{否则}.
$$

实际闭环权与字典多项式为

$$
\widehat W_{a,e}=\sum_{m=2}^{k} m\,(A\circ A^{m-1})_e,
\qquad
P_{k,d}(c)=\sum_{m=2}^{k} m\,W_d(m-1)\,c^{\,m},
$$

其中 $(A\circ A^{m-1})\_e=A\_{xy}(A^{m-1})\_{xy}$ 是逐元素积在边上的取值，$W\_d(L)$ 是 $\mathbb Z^d$ 上从 $0$ 到 $e\_1$ 的**精确 $L$ 步游走数**（纯拓扑整数）。

**正则性约定**：$\|c\|\_{C^{0,\alpha}}=\|c\|\_\infty+[c]\_{C^{0,\alpha}}$，$\|c\|\_{C^{1,\alpha}}=\|c\|\_{C^{0,\alpha}}+\|\nabla c\|\_{C^{0,\alpha}}$；$C^{1,1}$ 即 $\alpha=1$ 的 Lipschitz 梯度。

---

## §2 引理 R6.1（无绕环时的游走计数）

**引理**：若 $N>2k$，则对任意 $m\le k$ 与任意边 $e=\{x,y\}$，$(A^{m-1})\_{xy}$ 展开式中的游走集合与 $\mathbb Z^d$ 中从 $0$ 到 $e\_1$ 的 $(m-1)$ 步游走集合一一对应；特别地基数为 $W\_d(m-1)$，与 $x$、$a$ 无关。

**证明**：$(m-1)$ 步游走每一位移量的绝对值不超过 $m-1\le k-1<N/2$。因此绕环等价是在位移 $\lVert\cdot\rVert\_\infty<N/2$ 内唯一的，游走在环面上不会产生新的相交；把 $x$ 平移到 $0$、$y$ 平移到 $e\_1$ 即得双射。$\blacksquare$

**推论**：$k$ 固定时，$W\_d(m-1)$ 是常数；一旦 $k$ 随 $N$ 增长到 $k\gtrsim N/2$，绕环游程开始进入计数，字典多项式不再等于精确闭环计数。这就是 H5 必须固定 $k$ 的原因。

**绕环阈值（对抗审计 [`R6_refutation_attempt.md`](R6_refutation_attempt.md)）**：本文用的 $N>2k$ 是**充分**条件而非最锐条件。对抗审计给出更省的充分条件 $N\ge k$（用环面二分图的奇偶性），并给出“不设绕环条件就失败”的字面反例：$d=1$、$N=4$、$k=4$、$c\equiv1$ 时，绕环游程进入计数，使相对误差达 $2/7$。因此“$N$ 相对 $k$ 足够大”不是技术冗余，而是 H5 的前提之一。

---

## §3 引理 R6.2（反射—反转配对消去一阶项）

**引理**：对任意边 $e$ 与任意 $L\ge1$，设 $\Omega\_L$ 为 $e$ 端点间长 $L$ 的游走集合，$m\_j(w)$ 为游走 $w$ 第 $j$ 条边的中点。则

$$
\sum_{w\in\Omega_L}\sum_{j=1}^{L}\left(m_j(w)-x_e\right)=0 .
$$

**证明**：关于 $x\_e$ 作点反射 $\rho(p)=2x\_e-p$，再作游程反转 $\text{rev}$，得对合

$$
\tau=\text{rev}\circ\rho .
$$

$\rho$ 交换 $e$ 的两个端点，故 $\rho$ 把从 $x$ 到 $y$ 的游走送到从 $y$ 到 $x$ 的游走；反转后又回到从 $x$ 到 $y$ 的游走，故 $\tau$ 是 $\Omega\_L$ 到自身、且 $\tau^2=\text{id}$ 的双射。$\rho$ 把每条边中点送到其关于 $x\_e$ 的镜像 $2x\_e-m\_j$，反转只重排下标，故 $\tau$ 把中点多重集 $\{m\_j(w)\}$ 送到 $\{2x\_e-m\_j(w)\}$。于是

$$
\sum_j\bigl(m_j(\tau w)-x_e\bigr)
=-\sum_j\bigl(m_j(w)-x_e\bigr),
$$

对 $\Omega\_L$ 求和，左边等于 $-\,$左边，故为零（不动点游程自身配对，贡献也为零）。$\blacksquare$

这是 H5 从 O(a) 降到 O(a²) 的**唯一机关**；删掉它，所有结论退回 O(a)。

---

## §4 定理 R6.3（H5 及其尖锐速率）

**定理**：设 $d\ge1$、$k\ge2$ 固定、$N>2k$、$0<c\_{\min}\le c$。则对每个 $m\le k$、每条边 $e$，

$$
\left|(A\circ A^{m-1})_e-c_e^{\,m}W_d(m-1)\right|
\le W_d(m-1)\,c_e^{\,m}\,\eta_m(a),
$$

其中速率因子为

$$
\eta_m(a)\le
\begin{cases}
C_1\,(ma)^{\alpha}, & c\in C^{0,\alpha},\ \alpha\in(0,1],\\[1mm]
C_2\,(ma)^{1+\alpha}, & c\in C^{1,\alpha},\ \alpha\in(0,1],
\end{cases}
\qquad
C_1,C_2=C(d,k)\,\Bigl(1+\tfrac{\|c\|}{c_{\min}}\Bigr)^{k+1}.
$$

因此

$$
\left|\frac{\widehat W_{a,e}}{P_{k,d}(c_e)}-1\right|
\le \frac{\sum_{m=2}^{k} m\,W_d(m-1)\,c_e^{\,m}\,\eta_m(a)}{P_{k,d}(c_e)}
\le C(d,k)\,\Bigl(1+\tfrac{\|c\|_{C^{1,1}}}{c_{\min}}\Bigr)^{k+1} a^{s},
$$

其中指数 $s$ 依正则性取 $s=\alpha$（$C^{0,\alpha}$）或 $s=1+\alpha$（$C^{1,\alpha}$）。特别地 $c\in C^{1,1}$ 或 $C^2$ 时 $s=2$，即 H5 要求的 O(a²)。

**证明（要点）**：固定 $m$，写 $L=m-1$，取 $e=\{x,y\}$、$x\_e=m\_e$。

1. **游程展开**：由引理 R6.1，$(A^{L})\_{xy}=\sum\_{w\in\Omega\_L}\prod\_{j=1}^{L}c(m\_j(w))$，其中 $|\Omega\_L|=W\_d(L)$，且每条中点满足 $|m\_j-x\_e|\le La$。
2. **Taylor**：$c(m\_j)=c\_e+\nabla c(x\_e)\cdot(m\_j-x\_e)+r\_j$，其中 $|r\_j|\le\tfrac12 H\,|m\_j-x\_e|^{2}\le\tfrac12 H L^2a^2$（$H=[\nabla c]\_{C^{0,1}}$）。
3. **乘积展开**：$\prod\_j c(m\_j)=c\_e^{L}\prod\_j(1+\theta\_j)$，$\theta\_j=\bigl[\nabla c(x\_e)\cdot(m\_j-x\_e)+r\_j\bigr]/c\_e$，且 $|\theta\_j|\le\Theta:=\bigl(GLa+\tfrac12HL^2a^2\bigr)/c\_{\min}$，$G=\|\nabla c\|\_\infty$。于是

   $$
   \prod_j(1+\theta_j)=1+\sum_j\theta_j+\epsilon,
   \qquad |\epsilon|\le e^{\Theta}-1-\Theta\le \Theta^2e^{\Theta}.
   $$

4. **一阶项消去**：由引理 R6.2，$\sum\_{w}\sum\_j\nabla c(x\_e)\cdot(m\_j-x\_e)=\nabla c(x\_e)\cdot 0=0$。故

   $$
   \sum_{w}\sum_j\theta_j=\frac1{c_e}\sum_{w}\sum_j r_j,
   \qquad
   \Bigl|\sum_{w}\sum_j r_j\Bigr|\le W_d(L)\,L\cdot\tfrac12HL^2a^2 .
   $$

5. **合并**：$(A^{L})\_{xy}=c\_e^{L}\bigl[W\_d(L)+E\bigr]$，$|E|\le W\_d(L)\bigl(\tfrac12HL^3a^2/c\_{\min}+\Theta^2e^{\Theta}\bigr)$，即 $\eta\_m(a)=O(a^2)$。
6. **乘 $A\_{xy}=c\_e$ 并对 $m$ 求和**：$\widehat W\_{a,e}=\sum\_m m\,c\_e^{m}\bigl[W\_d(m-1)(1+O(a^2))\bigr]=P\_{k,d}(c\_e)\bigl(1+O(a^2)\bigr)$。
7. **除以字典**：$P\_{k,d}(c\_e)\ge c\_{\min}^2\min\_{m}\,mW\_d(m-1)>0$，故相对误差与绝对误差同阶。$\blacksquare$

**常数说明**：$C(d,k)$ 由 $\sum\_m mW\_d(m-1)$、$k^3$ 与 $e^{\Theta}$ 给出；数值上（§6）相对常数随 $k$ 的经验增长约为 $k^3$，与证明中 $L^3$ 的因子一致。

---

## §5 正则性边界与反例（速率是尖锐的）

| 正则性 | 一阶项是否消去 | 速率 | 数值验证（§6） | H5（$\varepsilon\_a\to0$） |
|:--|:--|:--|:--|:--|
| $C^{1,\alpha}$（梯度 Hölder） | 是（引理 R6.2） | $O(a^{1+\alpha})$ | α=1 实测指数 2 | 成立 |
| $C^{1,1}$＝Lipschitz 梯度 | 是 | $O(a^{2})$ | 比值精确 4.00 | **成立（H5 目标）** |
| $C^{0,\alpha}$（仅 Hölder，无梯度） | 否，配对差 $O(a^\alpha)$ | $O(a^{\alpha})$ | α=1,0.5,0.25 实测指数 1.00,0.50,0.25 | 成立 |
| 仅连续（任意模） | 否 | $O(\omega(ka))\to0$ | — | 成立 |
| 有跳跃（间断） | 分解失效 | $O(1)$ | 实测不收敛 | **失败** |

$$
\
\text{H5 的充要边界：}c\text{ 在细化下连续}\Longrightarrow\varepsilon_a\to0;\quad
c\text{ 间断}\Longrightarrow\text{H5 失败} .
\
$$

**固定 $k$ 是本质的**：$k\propto N$ 时 $L=m-1$ 也 $\propto N$，$\Theta=O(Gk a)=O(G k/N)$ 不再趋于零，且引理 R6.1 失效（绕环游程进入计数）。这与 [`G58`](G58_I2a_resolved_as_embedding_input.md) §2.3、[`Z7`](Z7_embedding_input_explicit_dictionary.md) §2.3 的“$k\propto N$ 无黎曼极限”同向。

**规则立方网格只有轴向边**：$d=2,3$ 的最近邻边都沿坐标轴，引理 R6.2 的反射 $\rho$ 逐坐标作用，故对 $d$ 无额外限制；非结构复形与斜边不在本文范围。

---

## §6 数值核验（[`R6_check.py`](R6_check.py)）

全部在**周期环面**上计算（无开链端点），$N>2k$。

### 6.1 O(a²) 与常数稳定性（$c=1+0.3\cos2\pi x$，$k=8$）

| $N$ | 64 | 128 | 256 | 512 | 1024 | 2048 |
|:--|--:|--:|--:|--:|--:|--:|
| 最大相对误差 | `1.47e-2` | `3.67e-3` | `9.17e-4` | `2.29e-4` | `5.73e-5` | `1.43e-5` |
| 误差 $\cdot N^2$ | 60.23 | 60.10 | 60.09 | 60.09 | 60.09 | 60.09 |

比值稳定为 $4.00$，常数 $\to60.09$：**O(a²) 成立，常数收敛**。

### 6.2 二维、三维同阶

| 维数 | 场 | $N$ 序列 | 最大相对误差 | 误差 $\cdot N^2$ |
|:--:|:--|:--|:--|--:|
| 2D | $1+0.2\cos2\pi x+0.15\cos2\pi y$ | 16→32→64→128 | `7.73e-2`→`1.91e-2`→`4.76e-3`→`1.19e-3` | 19.78→19.55→19.49→19.47 |
| 3D | $1+0.15\cos2\pi x$ | 8→12→16→24 | 见脚本 | 常数稳定 |

### 6.3 尖锐速率（一维，$k=8$）

| 场正则性 | $N$ 加倍时的比值 | 实测指数 | 定理预测 |
|:--|--:|--:|--:|
| $C^\infty$ | 4.00 | 2.00 | 2 |
| $C^{0,1}$ 折点 | 2.03 | 1.00 | 1 |
| $C^{0,1/2}$ | 1.41 | 0.50 | 0.50 |
| $C^{0,1/4}$ | 1.19 | 0.25 | 0.25 |
| 间断 | 1.00 | 0.00 | 失败 |

### 6.4 常数对 $k$ 的经验增长

$k=4,8,12$ 给常数约 $8.4,\,60.1,\,212.2$，经验指数约 $2.8\text{–}3.1$，与证明中的 $k^3$ 因子一致。

---

## §7 对 R1 的接续与诚实边界

1. **关闭的是 H5**：H5 从“未证的数值旁证”升级为定理 R6.3。R1-B 现在只需 H1–H4、H6 即可用实际闭环权。
2. **没有关闭的**：H3/H7（旋转类到光滑场 $c$ 的构造与 $C^2$ 正则迁移）仍属 R2；洛伦兹号差、精确锥、源作用量类别不属于本文。
3. **前提必须显式**：$k$ 固定、周期／公共极限、$c$ 连续且 $c\ge c\_{\min}>0$。删掉任一条，H5 失败或降阶（§5）。
4. **本文只处理规则立方网格**：非结构复形、随机图、D259 图册不在 R6.3 的覆盖内。
5. **数值是旁证不是证明**：§6 的数字与 §4 的证明分离；数值只用于确认常数与速率，不替代引理 R6.1–R6.2。

$$
\ \textbf{H5 解决后，主 }\Gamma\textbf{-收敛链的剩余承重点只剩 H3/H7（类到光滑场的构造与正则迁移）。}\
$$

**改动文件**：`R6_h5_dictionary_error_bound.md`、`R6_check.py`（另有两份独立复核：`R6_independent_check.py`、`R6_refutation_attempt.md`）。未修改既有主链文件。

---

## §8 核验

```
python3 R6_check.py
```

核验本文：

- $W\_d(L)$ 精确游走数与 Z7 表一致；
- 周期环面上 1D／2D／3D 的相对误差按 O(a²) 缩小，误差 $\cdot N^2$ 收敛；
- 反射—反转配对的一阶项数值为零；
- 尖锐速率：$C^\infty\to2$、$C^{0,1}\to1$、$C^{0,1/2}\to0.5$、$C^{0,1/4}\to0.25$、间断 $\to0$；
- 常数随 $k$ 的经验指数落在 $[2.5,3.3]$；
- 把 H5 写成定理而非数值猜测；把 $k$ 固定、周期、连续性写成显式前提。
