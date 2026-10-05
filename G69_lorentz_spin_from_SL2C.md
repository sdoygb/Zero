# G69 · 洛伦兹自旋：$SU(2)$ 升级到 $SL(2,\mathbb{C})$

**日期**：本轮 · **性质**：**升级构造**（空间旋转 → 洛伦兹）＋ **签名进入自旋结构**。
**等级标签**：【导出】/【核验】/【结构】/【结论】。
**核验**：[`G69_check.py`](G69_check.py) —— **独立实断言 21 / 结论行 0 / 不符 0**，退出码 `0`（0.23 秒）

$$
\boxed{\ SL(2,\mathbb{C})\ \text{双覆盖}\ SO^+(1,3);\qquad \text{签名由 Clifford 关系}\ \{\gamma^\mu,\gamma^\nu\}=2\eta^{\mu\nu}\ \text{进入自旋结构}。\ }
$$

---

## §0 目标

[`G66`](G66_SU2_double_cover_from_geometry.md)／[`G67`](G67_reflection_generates_spin_Z2.md) 给的是**空间**旋转的双覆盖 $SU(2)$（自旋 1/2）。本文把它升级到**洛伦兹**：$SL(2,\mathbb{C})$ 双覆盖 $SO^+(1,3)$。

---

## §1 显式构造

$$
X=x^\mu\sigma_\mu\ (\text{Hermite}),\qquad X\longmapsto A\,X\,A^{\dagger},\qquad A\in SL(2,\mathbb{C})
$$

| 核验（40 组随机 $A$） | 结果 |
|:--|:--|
| $AXA^{\dagger}$ 仍 Hermite | ✅ |
| $\det(AXA^{\dagger})=\det X=x^\mu x_\mu$ ⟹ **保 Minkowski 范数** | ✅ |
| $\Lambda^T\eta\Lambda=\eta$ | ✅ |
| $\det\Lambda=1$（固有） | ✅ |
| $\Lambda^0{}_0>0$（正交时） | ✅ |
| **同态** $\Lambda(AB)=\Lambda(A)\Lambda(B)$ | ✅ |
| **核** $=\{\pm I\}$：$\Lambda(-A)=\Lambda(A)$ | ✅ |

$$
\Longrightarrow\ SL(2,\mathbb{C})\ \text{是}\ SO^+(1,3)\ \text{的}\ \mathbf{2\ \text{对}\ 1}\ \text{双覆盖}。
$$

---

## §2 Cartan 分解 = 旋转 × boost

$$
A=UH,\qquad H=\sqrt{A^{\dagger}A}\ (\text{Hermite 正定}),\quad U=AH^{-1}\ (\text{酉})
$$

| 核验（40 组） | 结果 |
|:--|:--|
| $\det U=\det H=1$ | ✅ |
| $U$ 诱导**纯旋转**（$\Lambda^0{}_0=1$、$\Lambda^0{}_i=0$） | ✅ |
| $H$ 诱导 **boost**（$\Lambda^0{}_0\ge1$） | ✅ |

$$
\Longrightarrow\ SL(2,\mathbb{C})=\underbrace{SU(2)}_{\text{旋转}}\cdot\underbrace{H}_{\text{boost}};
\qquad \text{这就是"从 }SU(2)\text{ 升级"的}\textbf{结构性内容}。
$$

（**我第一版这里错了**：用 Cholesky 当"平方根"，而 $L^{\dagger}$ 是上三角、**不 Hermite**；正确做法是 $A^\dagger A$ 的谱分解 $\Rightarrow H=V\sqrt\Lambda V^\dagger$。）

---

## §3 双值性被**继承**（与 [`G66`](G66_SU2_double_cover_from_geometry.md)／[`G67`](G67_reflection_generates_spin_Z2.md) 自洽）

$$
2\pi\ \text{旋转}= \mathbf{-I},\qquad \Lambda(-I)=\mathbf{I}\ (\text{几何上什么都不做})
$$

$$
\Longrightarrow\ \text{洛伦兹情形}\textbf{继承}\text{空间旋转的双值性}——\text{不需要新机制}。
$$

---

## §4 **新东西**：签名进入自旋结构

Weyl 表示的 $4\times4$ $\gamma$ 矩阵（$\gamma^0=\begin{pmatrix}0&I\\I&0\end{pmatrix}$，$\gamma^i=\begin{pmatrix}0&\sigma_i\\-\sigma_i&0\end{pmatrix}$）：

| 核验 | 结果 |
|:--|:--|
| $(\gamma^0)^2=+I$（**类时**） | ✅ |
| $(\gamma^i)^2=-I$（**类空**） | ✅ |
| $\{\gamma^\mu,\gamma^\nu\}=2\eta^{\mu\nu}I$ | ✅ |

$$
\boxed{\ \text{Clifford 关系}\textbf{编码签名} \Longrightarrow \text{类时反射与类空反射的升格}\textbf{平方不同}。\ }
$$

**这给 [`G67`](G67_reflection_generates_spin_Z2.md) 的反射 $\mathbb Z_2$ 加了一条**：在洛伦兹情形它**按签名劈开**（类时 $+1$／类空 $-1$）。

$$
\Longrightarrow\ \text{G1 引理 5 导出的号差}\ (1,3)\ \textbf{被读进}\text{自旋结构——这是洛伦兹与欧氏情形的}\textbf{唯一实质差别}。
$$

---

## §5 升级链（自旋 1/2 的完整来历）

$$
\underbrace{D_L}_{\text{G27}}\ \longrightarrow\ \underbrace{SU(2)\ \text{双覆盖}\ SO(3)}_{\text{G66}:\ 2\pi=-1}\ \longrightarrow\ \underbrace{SL(2,\mathbb{C})\ \text{双覆盖}\ SO^+(1,3)}_{\textbf{本文}}
$$

而链条的**起点**是 [`G11`](G11_dimension_as_consistency.md) 的反射 $\mathbb Z_2$（它定 $D=4$，从而定空间 3 维），其**平方**给自旋 $\mathbb Z_2$（[`G67`](G67_reflection_generates_spin_Z2.md)）。

---

## §6 诚实边界

| 项 | 说明 |
|:--|:--|
| **仍是结构，不是动力学** | 本文与 [`G66`](G66_SU2_double_cover_from_geometry.md)／[`G67`](G67_reflection_generates_spin_Z2.md) 一样，构造的是**代数结构**；零和动力学是否**选**它，未证 |
| 空间 3 维是前提 | $\mathrm{Cl}(3,0)\to\mathrm{Cl}(1,3)$ 的升级以空间 3 维为前提，而 $D=4$ 是【条件】 |
| 未做 | **洛伦兹不变性**本身（[`G13`](G13_foliation_and_lorentz_invariance_gap.md) 的 I6）；自旋与**螺旋度 $\pm2$** 在同一代数里的显式关系；Majorana/Weyl 条件的物理读法 |
| 引用 | $\gamma$ 矩阵的 Weyl 表示是标准构造；本文**核验**其 Clifford 关系 |
| 影响 | 把自旋结构从空间旋转升级到洛伦兹；把**号差**读进自旋（与 [`G1`](G1_derivations_from_the_bottom_layer.md) 引理 5 接口）；不改变 G1–G68 的其余数值结论 |

---

## §7 核验

```
python3 G69_check.py     # 通过 21 / 不符 0，退出码 0（0.23 秒）
```

F1 **保 Minkowski 范数** · F2 **同态／核／固有／正交时** · F3 **Cartan 分解（旋转 × boost）** · F4 **双值性继承** · F5 **Clifford 关系编码签名**。
