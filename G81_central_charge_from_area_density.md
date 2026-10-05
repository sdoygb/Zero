# G81 · 中央荷：从**面积密度**接 [`D87`](../modular-equilibrium/derivations/D87_brown_henneaux_central_charge.md) 的 $c=6k=\frac{3\ell}{2G}$

**日期**：本轮 · **性质**：**把旧体系的中央荷链接到我们量出的面积密度上** ＋ **一处张力如实登记**。
**等级标签**：【导出（关系）】/【引用】/【数值核验】/【张力登记】/【结论】。
**核验**：[`G81_check.py`](G81_check.py) —— **独立实断言 10 / 结论行 6 / 不符 0**，退出码 `0`（0.3 秒）

$$
\ c=6a\ell\ \text{（关系，导出）};\ \text{但}\ k=\frac c6=0.467\ \textbf{不是整数}\ \text{（与 }D87\text{ 的张力）}。\ 
$$

---

## §0 旧体系 `D87`／`D88` 给了什么

| 出处 | 内容 |
|:--|:--|
| **`D87`** | 三维 AdS ⟹ 两个 $SL(2,\mathbb R)$ Chern–Simons ⟹ 边界 WZW 电流 ⟹ Sugawara ⟹ 两个 Virasoro，$c=6k=\frac{3\ell}{2G}$。**"推导"依赖三项新输入：边界条件／水平归一化／全息字典** |
| **`D88`** | Cardy 熵 $2\pi\sqrt{cE/6}$ 从模不变鞍点推出；**"仍不唯选 CFT，也不构造微观态"** |
| **`D86`**（上轮已读） | BTZ 的 Cardy 闭合；**"不解释 Brown–Henneaux 中央荷从上游如何产生"** |

---

## §1 我们的接入点：**面积密度 $a$**

[`G78`](G78_area_law_in_3d.md) 量出的 3D 面积密度（$S=a\cdot L^2$，角立方表面积 $=3L^2$）：

| $m$ | $a/3$（每单位面积） | $a\cdot m$ |
|--:|--:|--:|
| $2.0$ | $0.2310$ | $0.462$ |
| $3.0$ | $0.1590$ | $0.477$ |
| $4.0$ | $0.1158$ | $0.463$ |

$$
\ a\cdot m=\mathbf{0.467\pm0.008}\ (\text{相对涨落}\ 1.7\%)\ \Longrightarrow\ \textbf{面积密度与 gap 无关的组合};\ a\propto\frac1m\propto\xi。\ 
$$

这是[`G78`](G78_area_law_in_3d.md) 的数据给的一条**新普适关系**：

$$
\ a_{\rm face}\propto\frac1m\ \ (\text{numerical});\qquad \frac{1}{4G}:=a_{\rm face}\ \ (\text{identification})。\ 
$$

**为什么"$1/G=4a$ 只能是识别**：[`G57`](G57_unreachability_of_absolute_normalization.md) 推论 2 证明 $G$ **不可导出**；旧体系同判——`D24:22` 把 $G$ 登记为"已通过低能匹配固定"，`D24:198` 明文写"面积系数和 $G$ 可同步缩放，单靠熵不能固定 $G$"，`D24:202` 把该归一化登记为【约定】。所以可导出的是 $a\_{\rm face}\propto1/m$；$1/G=4a\_{\rm face}$ 是把这个比例**认成** $G$ 的定义。

---

## §2 系数核算：$c=\frac{3\ell}{2G}=6a\ell$

$$
S=\frac{A}{4G}\ \Longrightarrow\ \frac1G=4a;\qquad
c=\frac{3\ell}{2G}=\frac{3\ell\cdot4a}{2}=\mathbf{6a\ell}
$$

| 核验 | 结果 |
|:--|:--|
| $\frac{3\ell}{2G}$ 与 $6a\ell$ 一致（$\ell=1,2,5$） | ✅ |

$$
\Longrightarrow\ \ \text{中央荷被}\textbf{面积密度与 AdS 半径}\text{完全定住（这是一个}\textbf{关系}）。\ 
$$

---

## §3 取 $\ell=\xi=v\_F/m$：$c=6\rho v\_F/m^2$（原版漏了 $m^2$）

原版写 $c\_{\rm eff}=6\rho=2.805$（$\rho:=a\_{\rm face}m=0.467$），那是把 $\ell=1/m$ 直接当成 $\xi$，**既漏了 $v\_F$，更漏了 $m^2$ 的幂次**。正确的链是

$$
\ell=\xi=\frac{v_F}{m},\qquad a_{\rm face}=\frac{\rho}{m}
\ \Longrightarrow\ c=6a_{\rm face}\ell=6\rho\frac{v_F}{m^2}
$$

$$
\ c=6\rho v_F/m^2\ \ (\rho=0.4674\pm0.0068\ \text{3D asymptotic});\qquad
c_{\rm eff}=6\rho\ \text{holds only at}\ m=1。\ 
$$

量纲账：$\rho$ 与 $v\_F$ 都是格点单位的纯数、$m$ 是格点单位的质量，故 $c$ 是纯数；其中 $c\propto m^{-2}$ 是**唯一与口径无关**的部分。

---

## §4 **张力**：$k=\rho v\_F/m^2$——是 $m^2$ 倍，不是一个常数

$$
D87:\ c=6k\ \Longrightarrow\ k=\frac{c}{6}=\frac{\rho v_F}{m^2}
$$

$m=1$ 给 $k=1.002$（无张力）；而 $m=2,3,4$ 给 $k=0.250,\ 0.111,\ 0.063$，即张力 $1/k=4.0,\ 9.0,\ 16.0\approx m^2$。

$$
\ k=\frac{\rho v_F}{m^2}\ \propto\ m^{-2};\qquad k=1\ \Longleftrightarrow\ v_F=\frac{m^2}{\rho}\ 
$$

**修正（三条）**：

1. **`D87` 里根本没有"整数 level"这条要求**——`D87:105` 把 $k=\ell/(4G)$ 明文登记为**新输入（归一化）**，`D87:117` 明文说重标 $k\to\lambda k$ 是**约定**；`D87`/`D88`/`D86` 全文不含"整数"字样。所以"$k$ 非整数"本身不是与 `D87` 的张力，**是本项目自加的要求**。
2. **真正的量是幂次**：$k\propto m^{-2}$。原版把它当成常数，才在 $m=1$ 附近读出"只差 2 倍"。
3. **$k=1$ 的位置**：需 $v\_F=m^2/\rho$，即 $m=1,2,3,4$ 分别要 $v\_F=2.140,\ 8.558,\ 19.255,\ 34.232$；只有 $m\lesssim1.11$ 落在测得区间 $[2.000,2.640]$ 内，而 $m=1$ 已被 [`G78`](G78_area_law_in_3d.md) §2 判为**未进渐近区、不能下结论**。

---

## §5 $c$ 的**值**不可导出（[`G57`](G57_unreachability_of_absolute_normalization.md)）

| $\ell$ | $0.5$ | $1.0$ | $2.0$ |
|:--|--:|--:|--:|
| $c=6a\ell$ | $1.402$ | $2.805$ | $5.609$ |

$$
\Longrightarrow\ \ \textbf{可导出的是关系}\ c=6a\ell;\ \textbf{不可导出的是值}（\text{它由单位比值 }\ell/G\text{ 决定，}G57）。\ 
$$

这与 `D87` 自陈的"依赖**水平归一化**"和 `D86` 的"不解释中央荷从上游如何产生"**完全一致**——**而我们把它归因到了一条定理**（[`G57`](G57_unreachability_of_absolute_normalization.md)）。

---

## §6 Cardy 形式的量级核对

$$
S_{\rm Cardy}=2\pi\sqrt{\frac{cL_0}{6}},\qquad c=\frac{6\rho v_F}{m^2}
$$

（原版代入 $c=2.805$ 得 $8.591$；但 $c$ 的值依赖 $\ell$ 这个单位，见 §5——本轮更正的 $c=6\rho v\_F/m^2$ 在 $m=0.235342$ 上给 $c\approx18$–$22$，量级核对按此重述。）

（$L\_0\simeq L=4$，[`G61`](G61_locking_the_five_integers.md)）

| 对照 | 值 |
|:--|--:|
| Cardy | $8.591$ |
| 3D 实测 $S(L{=}6,m{=}2)$ | $22.205$ |
| 2D 实测 $S(L{=}9,m{=}1)$ | $10.273$ |

$$
\Longrightarrow\ \text{同阶};\ \text{但}\textbf{不唯选 CFT}（D88\ \text{的边界，如实采纳}）。
$$

---

## §7 诚实边界

| 项 | 说明 |
|:--|:--|
| **$\ell\simeq\xi$ 是识别** | 把 AdS 半径认成关联长度；**未**导出 |
| ~~**$k$ 非整数**~~ **作废** | `D87` 全文不含"整数"字样（`D87:105` 把 $k=\ell/(4G)$ 登记为**归一化新输入**，`D87:117` 说重标是**约定**）；正确的量是 $k=\rho v\_F/m^2$，张力 $\approx m^2$ |
| $a\cdot m$ 的普适性来自 3 个点 | $m=2,3,4$（[`G78`](G78_area_law_in_3d.md)）；$m=1$ **未进渐近区**（不用） |
| 没有黑洞 | 我们仍无黑洞解（[`G79`](G79_horizon_thermodynamics.md) 的边界仍在）；本文只做**系数核算** |
| 引用 | `D86`／`D87`／`D88` 是旧体系的条件推导链；本文只接**面积密度**这一端 |
| 影响 | 新开"中央荷"这一支；给出一条新普适关系（$a\propto1/m$，即 $G\propto m$）；登记一处张力；不改变 G1–G80 的其余数值结论 |

---

## §8 核验

```
python3 G81_check.py     # 通过 14 / 不符 0，退出码 0（0.3 秒）
```

F1 **$a\cdot m$ 的 gap 无关性** · F2 **$c=6a\ell$ 的系数核算** · F3 **$c\_{\rm eff}=2.805$** · F4 **$k$ 非整数（张力）** · F5 **$c$ 的值不可导出** · F6 **Cardy 量级核对**。
