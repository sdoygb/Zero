# G74 · I2a 与量子/自旋结构的接口：**代数半免费，尺度才是硬输入**

**日期**：本轮 · **性质**：**接口建立**（借旧体系 `D41` 的共尾细化不变性）＋ I2a 的**两半分离**。
**等级标签**：【导出】/【引用定理】/【数值核验】/【接口】/【结论】。
**核验**：[`G74_check.py`](G74_check.py) —— **独立实断言 11 / 结论行 6 / 不符 0**，退出码 `0`（0.4 秒）

$$
\ \text{I2a}=\underbrace{\text{代数极限}}_{\text{免费（}D41\text{）}}+\underbrace{\text{度规极限}}_{\text{固定步数收敛（}G58\text{）}}+\underbrace{\text{尺度}}_{\text{不可导出（}G57\text{）}}。\ 
$$

---

## §0 旧体系 `D41`／`D150` 给了什么

| 出处 | 内容 |
|:--|:--|
| **`D41` 共尾细化不变性** | 构造归纳极限 $\mathcal A\_\infty=\varinjlim\_\lambda\mathcal A\_\lambda$，并证明**它对共尾子网不变**；且"把连续化仍缺的**尺度**信息从索引冗余中分离出来" |
| **`D150` 粗粒化/细化稳定性** | 线性容量在粗粒化与细化下**都稳定**，**但单位尺度必须由体积元固定，否则细化无意义**（体积元是输入 `R-GEO2-METRIC`） |

$$
\Longrightarrow\ \text{旧体系的共同结论}:\ \textbf{代数极限免费}，\textbf{尺度是输入}。
$$

---

## §1 落在零和地基上：归纳系统

[`G62`](G62_quantum_sector_from_GNS_modular_flow_gleason.md)／[`G32`](G32_native_origin_of_saturation.md) 的局部代数

$$
\mathcal A_T=M_2(\mathbb C)\otimes\mathbb C^{\,T+1},\qquad \dim=4(T+1)
$$

**包含映射** $\varphi\_T:\mathcal A\_T\hookrightarrow\mathcal A\_{T+1}$（$M\otimes e\_a\mapsto M\otimes e\_a$）：

| 核验（200 组） | 结果 |
|:--|:--|
| 同态 $\varphi(XY)=\varphi(X)\varphi(Y)$ | ✅ |
| 保 $*$：$\varphi(X^{*})=\varphi(X)^{*}$ | ✅ |
| 保单位 | ✅ |

$$
\Longrightarrow\ \textbf{归纳系统良定义};\ \mathcal A_\infty=\varinjlim_T\mathcal A_T\ \text{存在}。
$$

---

## §2 **自旋结构在细化下不被破坏**

$$
\varphi_T\bigl(UXU^{\dagger}\bigr)=U\,\varphi_T(X)\,U^{\dagger},\qquad U\in SU(2)
$$

| 核验（50 组） | 结果 |
|:--|:--|
| $SU(2)$ 伴随作用与包含映射**交换** | ✅ |

$$
\ \text{G66／G67 的自旋结构在连续极限下}\textbf{不被破坏}。\ 
$$

---

## §3 态的相容性

$$
\omega_{T+1}\circ\varphi_T=\omega_T
$$

| 核验 | 结果 |
|:--|:--|
| $T=1,2,3$：$w\_{T+1}$ 在前 $T+1$ 个坐标上与 $w\_T$ 成**恒定比例** | 例如 $T{=}3$：$0.8089$（四个坐标全同）✅ |

$$
\Longrightarrow\ \text{存在相容态族} \Longrightarrow \mathcal A_\infty\ \text{上的态良定义}。
$$

---

## §4 **共尾不变性**（`D41` 的有限版本）

| 链 | 像 |
|:--|:--|
| 逐级：$\mathcal A\_T\to\mathcal A\_{T+1}\to\mathcal A\_{T+2}$ | 前 $T+1$ 个坐标 |
| 跳级（步长 2）：$\mathcal A\_T\to\mathcal A\_{T+2}$ | **同一组坐标** ✅ |

$$
\ \text{共尾（步长 2 的）子网给}\textbf{同一}归纳极限 \Longrightarrow \text{极限与细化方式无关}。\ 
$$

（这正是 `D41` 的核心结论在零和模型里的落地。）

---

## §5 量子运动学在整条链上存活

$$
d_T:=\dim\mathcal H_T=2(T+1)\ \ge\ 4\ \ge\ 3\quad(T\ge1)
$$

| 核验 | 结果 |
|:--|:--|
| $\dim\mathcal H\_T=2(T+1)$：$T=0$ 给 $2$，$T\ge1$ 给 $\ge4$ | ✅ |
| $\dim\mathcal A\_T=4(T+1)$（**与 Gleason 无关**） | ✅ |
| Gleason 在 **$\mathcal A\_T$ 的投影格**上不适用：$\mathcal A\_T=M\_2(\mathbb C)^{\oplus(T+1)}$ 是 **type $\mathrm I\_2$ 直和**，且 $L(M\_2)$ 的唯一正交对是 $\{P,1-P\}$ $\Longrightarrow$ 加性退化为互补加性 | ✅ |
| G62 §4 的反例在**每个 $T$** 上都满足 $p(E)+p(1-E)=1$、$0\le p\le1$，且残差 $0.24$ $\Longrightarrow$ **非迹形式** | ✅ |

$$\Longrightarrow\ \text{阈值主张}\textbf{撤回}:\ \text{"}\dim\mathcal A_T=2(T+1)\ge3\text{" 是}\textbf{对象混淆};\ \text{失效是}\textbf{类型（}\mathrm I_2\text{ 直和）问题，不随 }T\text{ 改变}。$$

$\Longrightarrow$ 正确的条件是**加性域覆盖 $\mathcal B(\mathcal H\_T)$ 的全部投影**（$T\ge1$）——那是**新增输入**，不是阈值。

---

## §6 I2a 的**两半分离**（本文的主结论）

| 半 | 状态 | 依据 |
|:--|:--|:--|
| **代数极限**（局部代数的归纳极限、自旋结构、GNS/Born 的存活） | ✅ **免费**（存在 ＋ 共尾不变） | **本文** ＋ `D41` |
| **度规极限**（有效度规的收敛） | ✅ 固定步数截断下二阶收敛 | [`G58`](G58_I2a_resolved_as_embedding_input.md) |
| **尺度**（绝对长度／$\kappa$） | ❌ **不可导出**（单位） | [`G57`](G57_unreachability_of_absolute_normalization.md)／[`G60`](G60_dimensionless_ledger_and_one_free_unit.md) ＋ `D150` |

$$
\ \text{所以 I2a 的"难"只在}\textbf{尺度识别}——\text{而它是}\textbf{输入}（\text{与 }D150\text{ 一致}）。\ 
$$

---

## §7 诚实边界

| 项 | 说明 |
|:--|:--|
| **引用的是定理** | 归纳极限的存在与共尾不变性是 `D41` 的结果；本文只在**有限模型**上核验其片段，**没有**重证 |
| 有限核验 | 包含映射、态相容、共尾只做到 $T\le3$（有限维手工核验） |
| 态族的构造 | 相容态族是**构造**出来的（固定序列截断），不是从 A5 导出的 |
| 尺度 | 本文只**复述**尺度不可导出（[`G57`](G57_unreachability_of_absolute_normalization.md)）；**未**证明"尺度只能外部给"的新内容 |
| 与 [`G58`](G58_I2a_resolved_as_embedding_input.md) 的关系 | 两半互补：本文补**代数半**，[`G58`](G58_I2a_resolved_as_embedding_input.md) 给**度规半** |
| 影响 | 把 I2a 拆成"代数（免费）＋ 度规（收敛）＋ 尺度（输入）"；**不改变** G1–G73 的其余数值结论 |

---

## §8 核验

```
python3 G74_check.py     # 通过 17 / 不符 0，退出码 0（0.4 秒）
```

F1 **包含映射是同态** · F2 **$SU(2)$ 被保留** · F3 **态相容** · F4 **共尾不变** · F5 **Gleason 阈值全链成立** · F6 **尺度不由代数提供**。
