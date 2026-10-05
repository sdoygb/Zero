# G67 · 反射 $\mathbb Z_2$ **生成** 自旋 $\mathbb Z_2$：$(b)=(a)^2$

**日期**：本轮 · **性质**：**代数关系的建立**（补上 [`G66`](G66_SU2_double_cover_from_geometry.md) §5 标为"待建"的那一处）＋ 三个 $\mathbb Z_2$ 的统一。
**等级标签**：【导出】/【核验】/【结构】/【结论】。
**核验**：[`G67_check.py`](G67_check.py) —— **独立实断言 18 / 结论行 0 / 不符 0**，退出码 `0`（0.17 秒）

$$
\boxed{\ \text{同一个反射的两个升格之积}\ =\ -1\ =\ \text{自旋 }\mathbb Z_2\ \Longrightarrow\ (b)=(a)^2\ }
$$

---

## §0 要补的是哪一处

[`G66`](G66_SU2_double_cover_from_geometry.md) §0 把几何扇区的三个 $\mathbb Z_2$ 分开了（反射 (a)／自旋 (b)／定向 (c)），但在 §5 如实登记：

> "本文只说明它们是**不同的** $\mathbb Z_2$ 并处同一几何扇区；**没有**给出它们之间更深的代数关系（该关系**待建**）。"

**本文把那处补上**——而且答案比"不同"漂亮得多：**(b) 由 (a) 生成。**

---

## §1 框架：$\mathrm{Cl}(3,0)$ 与 $\mathrm{Pin}(3)$

$$
v\in\mathbb R^3\ \longmapsto\ V=v_x\sigma_x+v_y\sigma_y+v_z\sigma_z,\qquad V^2=|v|^2 I
$$

| 对象 | Clifford 阶 | 作用 |
|:--|:--|:--|
| **反射**（奇元，grade 1） | $V$ | $x\mapsto -VXV$ |
| **旋转**（偶元，grade 2） | $E=VW$（两个反射之积） | $x\mapsto EXE^{\dagger}$ |

**核验**：
- $V^2=|v|^2I$（Clifford 关系）✅
- $-VXV$ 与几何反射 $x-2(x\cdot v)v$ **逐位一致** ✅；对合、$\det=-1$ ✅
- **两个升格 $\pm V$ 给同一个反射**（$R_{-v}=R_v$）✅

---

## §2 **核心**：$(b)=(a)^2$

$$
\boxed{\ V\cdot(-V)\ =\ -V^2\ =\ \mathbf{-1}\ }
$$

**核验**（30 组随机单位 $v$）：偏差 $\mathbf{2.24\times10^{-16}}$ ✅

**读法**：

> **同一个反射，取两次，在旋量上是 $-1$。**

而**在几何上它是恒等**：

| 层面 | $R_v\circ R_{-v}$ |
|:--|:--|
| **几何（O(3)）** | $\mathbf{I}$（恒等）✅ |
| **旋量（$\mathrm{Pin}(3)$）** | $\mathbf{-I}$ ✅ |

$$
\Longrightarrow\ \textbf{双值性的来源是【升格】，不是变换本身}。
$$

**并且这直接给出 (a) 与 (b) 的关系**：

$$
\boxed{\ (b)\ \text{是由}\ (a)\ \text{的}\textbf{两个升格相乘}\text{生成的}}——\text{不是"同一个 }\mathbb Z_2\text{"，而是"}\textbf{平方}\text{"}。
$$

---

## §3 一般情形：两个不同反射之积 = 旋转，角 $=2\arccos(v\cdot w)$

| 核验（40 组随机 $v,w$） | 结果 |
|:--|:--|
| $E=VW$ 酉且 $\det=1$（$\in SU(2)$） | ✅ |
| $E$ 实现 $R_v\circ R_w$（**组合恰是两个反射的复合**） | ✅ |
| 诱导旋转角 $=2\arccos(v\cdot w)$ | 最大偏差 $\mathbf{8.40\times10^{-13}}$ |

**这是"反射生成旋转"的标准机制**：两个夹角 $\varphi$ 的反射复合 = 转角 $2\varphi$ 的旋转。

**而 §2 是它的退化情形**：

$$
v\cdot w=-1\iff w=-v\iff \text{夹角 } \pi\ \Longrightarrow\ \text{转角 } 2\pi\ \Longrightarrow\ E=-I\ \checkmark
$$

$$
\boxed{\ \text{"2π 旋转"恰好是两个【同一反射的升格】之积}——\text{旋量双值性与反射结构}\textbf{同源}。\ }
$$

---

## §4 三个 $\mathbb Z_2$ 的统一（[`G66`](G66_SU2_double_cover_from_geometry.md) 的遗留问题收口）

| $\mathbb Z_2$ | 在这个框架里是什么 |
|:--|:--|
| **(a) 反射** | **奇元**（grade 1）：单反射，$\det=-1$，$V^2=+1$ |
| **(b) 自旋** | **两个升格之积** $V(-V)=-1$：$\mathrm{Pin}(3)$ 的中心元 |
| **(c) 定向** | **分级**（偶元 = 旋转，$\det=+1$；奇元 = 反射，$\det=-1$） |

$$
\boxed{\ \text{三者都在 }\mathrm{Pin}(3)\text{ 里}:\ (a)=\text{奇部},\ (c)=\text{分级},\ (b)=\text{中心}\ (\text{由 }(a)\text{ 平方得到})。\ }
$$

**所以 [`G65`](G65_spin_half_needs_3d_rotations.md) §5 我那句"用同一个 $\mathbb Z_2$"是错的、而 [`G66`](G66_SU2_double_cover_from_geometry.md) 的"完全不同"也是不完整的**；正确的说法是：

$$
\boxed{\ \text{不是同一个，也不是无关};\ \text{而是}\textbf{一个生成另一个}。\ }
$$

---

## §5 这对纲领的意义

| 之前 | **现在** |
|:--|:--|
| 自旋 1/2 需要"额外识别"（[`G64`](G64_is_spin_half_native.md)） | ❌ 错：循环扇区（2 维旋转）里根本没有（[`G65`](G65_spin_half_needs_3d_rotations.md)） |
| 自旋 1/2 与反射 $\mathbb Z_2$ 是**两个不相干的结构**（[`G66`](G66_SU2_double_cover_from_geometry.md)） | ❌ 不完整：**它们同一个 Clifford 结构** |
| — | ✅ **正解**：**反射（奇元）平方出中心 $-1$ = 自旋 $\mathbb Z_2$**；而 G11 的反射 $\mathbb Z_2$ 已经**原生**在手上 |

$$
\Longrightarrow\ \text{自旋 1/2 的来源被追到}\textbf{公理侧的同一个反射结构}——\text{与定维数用的那个是同一个}。
$$

**这正是 [`G65`](G65_spin_half_needs_3d_rotations.md) 想要而没做到的那条链**：

$$
\text{G11 的反射}\ \mathbb Z_2\ \xrightarrow[\text{本文}]{\ (a)^2=(b)\ }\ \text{自旋}\ \mathbb Z_2\ =\ 2\pi\ \text{旋转的}\ -1\ \Longrightarrow\ \text{自旋 1/2}
$$

---

## §6 诚实边界

| 项 | 说明 |
|:--|:--|
| **是代数关系，不是动力学** | 本文给出的是 $\mathrm{Cl}(3,0)/\mathrm{Pin}(3)$ 里的**代数关系**；它说明"**结构上可达**"，**不**说明零和动力学一定选它 |
| **空间 3 维仍是前提** | 整个构造用 $\mathrm{Cl}(3,0)$，即**空间 3 维**——它来自 $D=4$（[`G11`](G11_dimension_as_consistency.md)），而后者是**【条件】** |
| $H^2$ 与 $D=5$ | [`G65`](G65_spin_half_needs_3d_rotations.md) 已诚实登记：$D=5$ 也有 2 维旋量 ⟹ **唯一性仍来自 (a) 定维数** |
| 未做 | **洛伦兹**版本（$\mathrm{Cl}(1,3)$、$SL(2,\mathbb C)$）；与**螺旋度 $\pm2$** 在同一代数里的显式关系 |
| 影响 | 补上 [`G66`](G66_SU2_double_cover_from_geometry.md) §5 的"待建"；**统一**三个 $\mathbb Z_2$；把自旋 1/2 的来源追到公理侧的反射结构；不改变 G1–G66 的其余数值结论 |

---

## §7 核验

```
python3 G67_check.py     # 通过 18 / 不符 0，退出码 0（0.17 秒）
```

F1 **$\mathrm{Cl}(3,0)$ 表示** · F2 **反射与两个升格** · F3 **核心 $V(-V)=-I$** · F4 **反射之积 = 旋转，角 $2\arccos$** · F5 **$v\cdot w=-1\Rightarrow$ 几何恒等／旋量 $-1$** · F6 **分级与 $\mathrm{Pin}(3)$**。
