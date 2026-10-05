# R47 · 因果结构供出洛伦兹签名：`R46` 的价格**可以付**

**日期**：2026-10-03
**性质**：**判据裁决（正面）＋ 链条收束**。执行 `R46` (R46-5)：检查因果/锥结构能否供出洛伦兹签名。
**依赖**：[`R46`](R46_pair_ledger_objects_and_simplex.md)、[`R45`](R45_ledger_form_scan.md)、[`R44`](R44_survival_vs_contextuality_no_go.md)、[`G59`](G59_I7_settled_native_cone_and_its_residue.md)、[`R33`](R33_action_phase_match_project.md)、[`G57`](G57_dimensional_constants_no_go.md)、[`G60`](G60_dimensionless_ledger_and_one_free_unit.md)、[`STATUS`](STATUS.md)。
**探针**：[`R47_signature_from_causal_cone_probe.py`](R47_signature_from_causal_cone_probe.py) → [`R47_signature_from_causal_cone_results.json`](R47_signature_from_causal_cone_results.json)。
**核验**：[`R47_check.py`](R47_check.py)。

$$
\begin{aligned}
&\text{A smooth cone field} \iff \text{conformal Lorentz structure (the cone is the zero set of a quadratic form,}\\
&\qquad\text{signature } (1,D-1)\text{).}\\
&\text{Numerics: } D=2..6 \text{ all give signature } (1,\,D-1) \ \Longrightarrow\ \text{Lorentzian}.\\
&\text{The effective cone suffices: the signature is a}\ \textbf{discrete invariant},\\
&\qquad\text{so the exponentially small tails of } G59 \text{ only blur the boundary, not the topology.}\\
&\therefore\ \text{the price can be paid: } \binom{D+1}2 \text{ is legitimate} \Longrightarrow D=4 \text{ and contextuality can coexist.}\\
&\text{Residual: the conformal factor / scale is still not derived (} E4/G57\text{).}
\end{aligned}
$$

> **一句话**：上一轮留下的价格——"单纯形字典是欧氏的、签名须外供"——**这一轮付掉了**。理由是一条标准事实加一条不变量论证：**光滑锥场 ⟺ 共形洛伦兹结构**（锥就是某个二次型的零集，其符号差必为 $(1,D-1)$），而**签名是离散不变量**，所以 `G59` 那个"只有有效锥、锥外指数小但非零"的毛病**无害**——指数尾巴只把锥的边界抹糊（数值：锥外占比随阈值从 `0.395` 变到 `0.0515`），**但改不了锥的拓扑**（pointed、convex），签名照旧是洛伦兹。于是 `R44` 的互斥被**完全绕过**：**`D=4` 与单体量子性可以同时到手。**

> **【后续更正｜[`R48`](R48_exact_cone_and_effective_cone.md)（2026-10-03）】** 本文 §3 残留 2「精确光锥仍缺」是**措辞错误**：`R48` 证明**支持锥一直是精确的**（$\text{supp}\rho\_n\subseteq[-n,n]$，$B$ 无关）；`G59` 的"锥外指数小"是**锥边可见性** $\varepsilon(B)=\tfrac12\log(4/B)$ 而非锥的精确性。因此**本文的签名判决不需要"离散不变量"辩护**（它面对的是精确锥），结论本身**加强**且不变。本文其余内容保持原状（按 `STATUS.md` §9 的纪律，原文不原地改写）。

---

## §0 判决摘要

| `D` | 签名 `(正, 负, 零)` | 洛伦兹 |
|--:|:--|:--|
| 2 | `(1, 1, 0)` | ✅ |
| 3 | `(1, 2, 0)` | ✅ |
| 4 | `(1, 3, 0)` | ✅ |
| 5 | `(1, 4, 0)` | ✅ |
| 6 | `(1, 5, 0)` | ✅ |

| 项 | 结论 | 状态 |
|:--|:--|:--|
| 锥 ⟹ 签名 | 符号差 $(1,D-1)$（`D=2..6` 全对） | **已证（数值）** |
| 有效锥够用 | 签名离散 ⇒ 指数尾巴无害 | **已判（不变量论证）** |
| 模糊性的量级 | 锥外占比随阈值 `ε`：`0.395 → 0.302 → 0.222 → 0.154 → 0.052`（`ε: 10⁻¹→10⁻⁶`） | **已量化** |
| **判据 (R46-5)** | **能供出** ⇒ `C(D+1,2)` 合法 | **已裁决** |
| 链条后果 | **`D=4` 与单体语境性可共存** | **已达** |
| 残留 | 共形因子／尺度不导出（与 `G57` 一致） | **已登记** |
| 四维 GR 的**动力学** | — | **未由此推出** |

---

## §1 论证链

**(1) 锥 ⟺ 二次型零集。** 若因果可达结构在每点由一个**尖凸锥**给出，则该锥是某个二次型的零集；其符号差只有两种可能：$(1,D-1)$（洛伦兹）或正定（无锥）。$\Rightarrow$ **有锥就有洛伦兹签名**。

$$
\text{锥的拓扑（pointed，不含直线）}\ \Longleftrightarrow\ \text{恰有一个时间方向}\ \Longleftrightarrow\ \text{符号差 }(1,D-1).
\qquad\text{(R47-1)}
$$

**(2) 有效锥够用（不变量论证）。** `G59` 的锥是**有效**的：锥外影响指数小而非零。但：

$$
\textbf{签名是离散不变量}\ \Longrightarrow\ \text{指数小尾巴只移动边界，不改变拓扑，}\textbf{故签名不变}。
\qquad\text{(R47-2)}
$$

数值佐证：以阈值 $\varepsilon$ 定义"有效锥"，锥外占比随 $\varepsilon$ 变化（`0.395 → 0.052`），**但锥的指向性/凸性不变**，故签名不变。

---

## §2 链条收束：`R44` 的互斥被绕过

$$
\begin{aligned}
&\text{账本多重度} &&\binom{D+1}2\ (\text{单纯形边，}R46)\\
&\text{签名} &&(1,D-1)\ (\text{因果锥，}R47)\\
&\text{峰位} &&D=4\iff q\in(0.6,\ 2/3)\ (R45)\\
&\text{语境性} &&q>3/5\ (R44)\\
&\Longrightarrow\ &\textbf{取 } q\in(0.6,\,2/3)\ \text{即可}\ \textbf{同时}\text{得到 } D=4\ \text{与单体量子性}\ \checkmark
\end{aligned}
\qquad\text{(R47-3)}
$$

**载体复核**（`R45` 已算）：文档例 `(2,5,20,100)`：`q=0.6466` ✅∈`(0.6,2/3)`、`S=2.0328>2` ✅、单纯形字典峰 `D=4` ✅；两层族 `p=0.8`：`q=0.6585` ✅、`S=2.0150` ✅、峰 `D=4` ✅。

$$
\textbf{逃生口成立，且价格已付：}D=4\ \text{与单体量子性共存。}
\qquad\text{(R47-4)}
$$

---

## §3 残留与边界（照旧说清）

| # | 残留 | 说明 |
|--:|:--|:--|
| 1 | **共形因子／尺度** | 锥只给共形结构；尺度仍不导出——**与 `G57` 的"量纲常数不可导出"一致**，不是新缺口 |
| 2 | **锥的精确化** | 签名对有效锥稳健，但**精确光锥**仍缺（`G59` 的开放项未关闭） |
| 3 | **动力学** | 由 (欧氏距离 + 锥) 到 **Einstein 方程**仍需 Lovelock（`G1`）与 L1 的时间识别 |
| 4 | `R32` 唯一性 | `R45-4` 的削弱仍在（39 个形式可选，需另一原则钉死） |
| 5 | `2\pi`、`E1` | 仍未解决 |

$$
\ \text{The signature is derived; the scale and the dynamics are not.}\
$$

---

## §4 没有推出什么

1. **没有**推出 Einstein 方程或任何动力学（`G1` 的 Lovelock 结构仍在原状）。
2. **没有**给出精确光锥（`G59` 仍只有有效锥）。
3. **没有**恢复 `R32` 的账本唯一性。
4. **没有**解决 `2\pi`（`R33` S1）、`E1`、`III_1` 的唯一性。
5. **没有**由四维标签推出度规的具体形式或 GR。

---

## §5 核验命令

```bash
python3 R47_signature_from_causal_cone_probe.py
python3 R47_check.py
python3 R46_check.py
python3 R45_check.py
python3 STATUS_check.py
```
