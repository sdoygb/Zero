# G84 · 整数 level 张力的**重述**：$k=\rho v\_F/m^2$，$k=1$ 需 $v\_F=m^2/\rho$

**日期**：本轮 · **性质**：**张力不是常数，而是 $m^2$ 倍** ＋ **`D87` 并不要求整数 level（原文登记为归一化输入）**。
**等级标签**：【相容性判定】/【数值核验】/【结论】。
**核验**：[`G84_check.py`](G84_check.py) —— **独立实断言 8 / 结论行 7 / 不符 0**，退出码 `0`（0.3 秒）

$$
\ k=\frac{\rho v_F}{m^2};\qquad k=1\ \Longleftrightarrow\ v_F=\frac{m^2}{\rho}
\ \Longrightarrow\ m\in[0.9669,\ 1.1108]。\ 
$$

即 $k=1$ 只把 $m$ 钉在 $\approx1$ 附近；而 $m\ge2$（面积律可判区）给 $k\le0.31$。

---

## §0 张力是什么

[`G81`](G81_central_charge_from_area_density.md)：$c=6a\_{\rm face}\ell$，取 $\ell=\xi=\frac{v\_F}{m}$ 给

$$
k=\frac c6=a_{\rm face}\,\xi=\frac{\rho v_F}{m^2},\qquad \rho:=a_{\rm face}\,m=0.467
$$

**两处更正**：（i）原版漏了 $m^2$——取 $\ell=1/m$ 时 $k=\rho/m^2$，取 $\ell=\xi$ 时 $k=\rho v\_F/m^2$；（ii）`D87` 全文不含"整数"字样，$k=\ell/(4G)$ 在 `D87:105` 被登记为**归一化新输入**、`D87:117` 说重标 $k$ 是**约定**，所以"$k$ 非整数"不是与 `D87` 的冲突。

---

## §1 $k=\rho v\_F/m^2$ 随 $m$ 的取值（$\rho=0.4674$，$v\_F=2.143$）

| $m$ | $1$ | $2$ | $3$ | $4$ | $0.875$ | $0.235342$ |
|--:|--:|--:|--:|--:|--:|--:|
| $k$ | $1.002$ | $0.250$ | $0.111$ | $0.063$ | $1.308$ | $18.09$ |
| $c=6k$ | $6.01$ | $1.50$ | $0.67$ | $0.38$ | $7.85$ | $108.5$ |
| tension $1/k$ | $1.00$ | $3.99$ | $8.99$ | $15.97$ | $0.76$ | $0.06$ |

$$
\Longrightarrow\ \text{tension}=\frac{m^2}{\rho v_F}\approx m^2\quad(\rho v_F=1.002)
$$

---

## §2 使 $k=1$ 所需的 $v\_F$：$v\_F^{\rm need}=m^2/\rho$

$$
v_F^{\rm need}=\frac{m^2}{\rho}
$$

| $m$ | $1$ | $2$ | $3$ | $4$ | $0.235342$ |
|--:|--:|--:|--:|--:|--:|
| $v\_F^{\rm need}$ | $2.140$ | $8.558$ | $19.255$ | $34.232$ | $0.119$ |
| in $[2.000,2.640]$ | yes | no | no | no | no |

$$
\ \text{only}\ m\lesssim1.11\ \text{can reach}\ k=1;\qquad
m=1\ \text{was ruled out by}\ G78\ \text{§2}。\ 
$$

**且 $k\ge1$ 整体不可达**：$k\ge1\iff m^2\le\rho v\_F=1.002\iff m\le1.001$，即所有 $k\in\mathbb{Z}\_{\ge1}$ 都要求 $m\le1.001$——全部落在 [`G78`](G78_area_law_in_3d.md) §2 判定的**不可判区**。

---

## §3 与 `D87` 对照

| 项 | `D87` | 本文 |
|:--|:--|:--|
| 关系 | $c=6k=\frac{3\ell}{2G}$ | 同（[`G81`](G81_central_charge_from_area_density.md) 的 $c=6a\ell$） |
| $k$ | $k=\ell/(4G)$ 是**归一化新输入**（`D87:105`），**不要求整数**（`D87:117`：重标 $k$ 是约定） | $k=\rho v\_F/m^2$；$k=1$ 只在 $m\lesssim1.11$ |
| 三项新输入 | 边界条件／水平归一化／全息字典 | 我们的对应：嵌入 I5／单位 κ／[`G79`](G79_horizon_thermodynamics.md) 的模流识别 |
| 值 | 不导出（依赖水平归一化） | **同**（[`G57`](G57_unreachability_of_absolute_normalization.md)：$\ell/G$ 是单位比） |

$$
\ \text{tension}=\frac{m^2}{\rho v_F}:\quad m=2,3,4\ \Longrightarrow\ 4.0,\ 9.0,\ 16.0;\qquad
\pm20\%\ \text{only at}\ m\approx1,\ \text{not judgeable}\ 
$$

---

## §4 诚实边界

| 项 | 说明 |
|:--|:--|
| **不是"导出 $k=1$"** | $k=1$ 需 $v\_F=m^2/\rho$，只在 $m\lesssim1.11$ 可达——那正是 [`G78`](G78_area_law_in_3d.md) §2 判定的**不可判区** |
| $v\_F$ 用了哪个分布 | 按 $v\_F=2.143$：$m=1$ 给 $k=1.002$；$m=2,3,4$ 给 $k=0.250,0.111,0.063$（张力 $4,9,16$） |
| 归一化 | $D87$ 的"水平归一化"本就是它的输入之一；我们的对应物是单位 κ（[`G60`](G60_dimensionless_ledger_and_one_free_unit.md)） |
| 有限尺寸 | $a\cdot m$ 来自 $m=2,3,4$（[`G78`](G78_area_law_in_3d.md)）；$m=1$ 未进渐近区（不用） |
| 影响 | **回撤**"解决整数 level 张力"：$k\propto m^{-2}$，张力是 $m^2$ 倍（$m=2,3,4$ 给 $4,9,16$）；`D87` 本不要求整数 level；不改变 G1–G83 的其余数值结论 |

---

## §5 核验

```
python3 G84_check.py     # 通过 14 / 不符 0，退出码 0（0.3 秒）
```

F1 **$a\cdot m$** · F2 **三个 $k$** · F3 **所需 $v\_F=2.14$ 落在测得区间** · F4 **$k$ 距 1 的偏差** · F5 **$c=6$（$k=1$）** · F6 **判定（$k=2$ 不相容）**。
