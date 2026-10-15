# 分布机：$K$ 被消除，但 $T$ 仍未被定

**日期**：本轮 · **状态**：$T$–$K$ 依赖**已消除**（通过取消 $K$）；$T$ 另需来源
**层指标**：演化层（$\mathrm{L2}$）

---

## §0 一句话

$$
\boxed{\;
\textbf{把状态从"走者计数"改成"每点分布"}\;\Longrightarrow\;K\ \textbf{消失}\;\Longrightarrow\;\text{"}T\ \text{依赖}\ K\text{"}\ \textbf{问题消失} .
\;}
$$

$$
\boxed{\;\textbf{但}:\ \text{分布机}\textbf{自动有界}\ \Longrightarrow\ \text{无饱和}\ \Longrightarrow\ T\ \textbf{没有读出量} .\;}
$$

---

## §1 动机：$K$ 从哪来

**旧（计数机）**：

$$
\text{状态}=\text{走者数}\;\xrightarrow{\ \text{每代分支（度 }k\text{）}\ }\;\text{无界}
\;\Longrightarrow\;\text{必须外给容量 }K
\;\Longrightarrow\;T=f(K)\ .
$$

**框架给的形态**（局部代数）：

$$
\mathcal A_i=M_2(\mathbb C)\otimes\mathbb C^{T+1},\qquad
\dim\mathcal A_i=4(T+1),\qquad
\textbf{可区分态}=2(T+1) .
$$

$$
\Longrightarrow\quad \textbf{每点容量}\ \textbf{内生于代数}\;\Longrightarrow\;\text{状态应是}\textbf{分布}，\text{不是计数} .
$$

---

## §2 分布机：实现与核心发现

**状态**：每点上可区分状态的集合（可区分 $\iff$ 起点标签不同）
**推进**：全分支一步；每点超出容量者被粗粒化
**容量**：$\mathrm{cap}=2(T+1)$（**内生**）

### 2.1 ★ 核心发现：自动有界，不需要 $K$

取 $\mathrm{cap}=200$（**完全不绑定**），跑 $50\,000$ 代：

| 量 | 结果 |
|:--|:--|
| `max_site` 稳态分布 | $4$：$49\,994$ 次；$1,2,3$：各 $1$–$3$ 次 |
| `n_states` 范围 | $[19,\ 31]$ |

$$
\boxed{\;\textbf{分布机自动有界};\ \mathrm{cap}=200\ \textbf{完全不绑定};\ \textbf{不需要外部 }K .\;}
$$

$$
\Longrightarrow\quad \boxed{\;\textbf{自然容量}=4\ (\text{与 cap 无关}) .\;}
$$

### 2.2 自然容量由什么定

| $L$ | $2$ | $3$ | $4$ | $5$ | $6$ | $8$ |
|:--|--:|--:|--:|--:|--:|--:|
| 稳态 `max_site` | $4$ | $1$ | $4$ | $1$ | $4$ | $4$ |

$$
\boxed{\;\textbf{只依赖 } L\ \textbf{的奇偶性}:\ \text{偶 }L\Rightarrow4;\ \text{奇 }L\Rightarrow1 .\;}
$$

**奇数 $L$ 退化的原因**：净位移不可能为零 $\Longrightarrow$ **闭合不存在** $\Longrightarrow$ 状态塌缩。

$$
\Longrightarrow\quad \textbf{自然容量是}\textbf{图与动力学}\text{的性质};\ \text{不是自由参数} .
$$

---

## §3 结果：$K$ 消除了，$T$ 仍缺

### 3.1 得到了什么

$$
\boxed{\;\textbf{"}T\ \text{依赖}\ K\text{" 这个问题}\textbf{已消除}\text{——因为 }K\ \textbf{不再存在} .\;}
$$

分布机自动有界，无需任何外部容量参数。

### 3.2 失去了什么

$$
\mathrm{cap}=200\ \text{而自然容量}=4\ \Longrightarrow\ \textbf{没有任何东西饱和}
\;\Longrightarrow\ \textbf{毁灭不触发}\;\Longrightarrow\ T\ \textbf{无读出量} .
$$

$$
\boxed{\;\text{净结果}:\ \text{旧问题消失};\ \text{新问题}:\ T\ \text{需要}\textbf{另一个}\text{来源（目前无）} .\;}
$$

---

## §4 与既有结论的对照

| 项 | 状态 |
|:--|:--|
| 计数机（旧） | $T\propto\log K$，$K$ 无出处 |
| 分布机（新） | **有界、无 $K$**、$T$ 无读出量 |
| 框架的局部代数 | $2(T+1)$ 可区分态——**内生于代数** ✓ |
| 那个 $T$（年龄原子） | **仍未与动力学接上** |

$$
\boxed{\;\textbf{局部代数里的 }T\ \text{与动力学里的周期，仍是}\textbf{两个未接通的量} .\;}
$$

---

## §5 结论

$$
\boxed{\;
\begin{aligned}
&\textbf{实现}:\ \text{状态}=\text{每点分布};\ \text{容量}=\mathrm{cap}=2(T+1)\ \text{内生}.\\[1mm]
&\textbf{发现一}:\ \text{自动有界（}\mathrm{cap}\ \text{不绑定）}\Longrightarrow K\ \textbf{消失}.\\[1mm]
&\textbf{发现二}:\ \text{自然容量}=4;\ \textbf{只依赖 } L\ \text{的奇偶性}.\\[1mm]
&\textbf{结果}:\ \text{"}T\ \text{依赖}\ K\text{" 消失};\ \text{但 } T\ \textbf{无读出量}.\\[1mm]
&\textbf{遗留}:\ \text{代数里的 }T\ \text{与动力学的周期}\textbf{未接通};\ \text{容量公式仍未被触发}.
\end{aligned}
}
$$

**一句话：**

$$
\boxed{\;\textbf{分布机把 K 消掉了（问题消失），代价是毁灭不再触发——T 需要另一条路。}\;}
$$

---

## §6 下一步的三个候选（登记，**不选**）

| # | 路线 | 内容 |
|--:|:--|:--|
| **1** | 让自然容量**成为**容量 | 取 $2(T+1)$ = 自然容量 $=4\Longrightarrow T=1$；与**代数里的 $T$** 无关，需判其地位 |
| **2** | 换触发判据 | 不用"饱和"触发毁灭，改用别的（例如记录层到某比例） |
| **3** | 接受 $T$ 不可由动力学定 | 回到"$T$ 是一笔具名输入" |

$$
\boxed{\;\text{三条都未验证};\ \text{本文}\textbf{不选} .\;}
$$

---

## 附录 A　复现

```bash
cd ZeroGPU/src
/usr/bin/python3 distribution_engine.py    # 7 条断言
```

## 附录 B　本裁决**不**做的事

| # | 不做 |
|--:|:--|
| **1** | **不**声称 $T$ 已有值 |
| **2** | **不**声称 $T=1$（那是候选 1，未验证） |
| **3** | **不**声称容量公式被触发（它未被触发） |
| **4** | **不**声称分布机是唯一忠实实现 |

## 附录 C　符号表

| 符号 | 含义 |
|:--|:--|
| 计数机 | 状态 = 走者数（旧） |
| 分布机 | 状态 = 每点可区分状态集合（新） |
| 自然容量 | 稳态下每点状态数的上界（$\mathrm{cap}$ 不绑定时） |
| $\mathrm{cap}=2(T+1)$ | 内生于局部代数的每点容量 |
