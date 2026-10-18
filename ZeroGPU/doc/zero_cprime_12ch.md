# C′ 接 12 通道：真实 content 的 $q$

## 结论

把 C′ 的底 $W$（$4\times4=16$ 类）映射成 **cycle16** 的 seed 后，12 通道引擎给出：

| 量 | baseline | **C′ 耦合** | 观测 |
|:--|--:|--:|--:|
| content $q$（$n$ 的 16 项） | $0.062539$ | $0.062539$ | — |
| **record $q_{\rm live}$** | $0.062500$ | **$0.069576$** | $0.068303$ |
| $q_{\rm hist}$ | $0.062500$ | $0.069576$ | — |
| $n_{std}$ | $0.024822$ | $0.024822$ | — |

**关键结果**：真实 content 的 $q$ 本来就是 $1/16=0.0625$；C′ 的 seed 不改变 $n$，但把 **record $q_{\rm live}$ 抬到 $0.06958$**，离观测 $0.06830$ 只差 **$0.0013$**。

这是 C′ 接 12 通道后的**第一个正面结果**：$q_{\rm live}$ 对上了。

![C′ 12 通道 q](ZeroGPU/results/zero_cprime_12ch.png)

## 怎么接的

- 12 通道引擎：`ZeroEvo(cycle16, L=4, T=4, micro=True, c_h=1.0)`；
- C′ 在 K4 上独立演化：三份 $W$、底全灭/中留/顶不动、中只约束、顶独立、flow；
- 每代把 C′ 的 `Wb`（$4\times4=16$ 项）归一化，**按 `site*4+fiber` 映射到 cycle16 的 16 个 seed**：
  $$
  \text{seed}_{4j+f}\propto (W_b)_{j,f};
  $$
- 引擎照常走 L1→L2、闭合、记录；C′ 只改 seed，不改 $n$。

## 两个 $q$ 要分开

| $q$ | 对象 | baseline | C′ | 目标 |
|:--|:--|--:|--:|--:|
| content $q$ | $n$ 的 16 项 | $0.062539$ | $0.062539$ | — |
| record $q_{\rm live}$ | 16 槽闭合记录 | $0.062500$ | **$0.069576$** | $0.068303$ |

- 内容 $q$ 是 $1/16$，本来就接近观测；**seed 不改 $n$**，所以 content $q$ 不动；
- 观测 $q=0.0683$ 更可能是**记录 $q$**，不是 content $q$；
- C′ 的 $q_{16}=0.0906$ 是**代理**（底 W 的类分布），接进引擎后记录 $q_{\rm live}=0.0696$——代理偏集中，引擎混合后回到目标附近。

## 判据（3/3）

1. baseline content $q\approx1/16$：PASS（$0.062539$）；
2. coupled $q$ 有定义：PASS；
3. coupled $q_{\rm live}$ 接近观测（$|\Delta|<0.01$）：PASS（$0.069576$ vs $0.068303$）。

## 边界

- 这是 **cycle16 的 12 通道**，不是 K4；C′ 的 4 site 通过 `site*4+fiber` 映射到 16 顶点，这是一条**新识别**，需要登记。
- C′ 只改 seed，不改 $n$；content $q$ 仍由原引擎决定。
- $q_{\rm live}$ 在 $0.0657$–$0.0724$ 间变，均值 $0.0696$；是否周期 4 还需单独看。
- 16/8/4 的**层级份额**在 12 通道里没有直接对应物（只有一个 content $n$）；份额仍是 C′ 的 K4 结果。
- 无新量纲锚：映射、seed、$q$ 都无量纲；$v$ 仍是唯一锚。

## 下一步

1. **看周期 4**：$q_{\rm live}$ 的 FFT 是否有 bin $T=4$ 分量；若有，说明 C′ 的周期进入了 12 通道记录。
2. **接 16/8/4 份额**：12 通道只有一个 $n$，要测层级份额，需要把三份 $W$ 作为**记录/读出**接进去，而不是只接 seed。
3. **检查映射**：`site*4+fiber` 是否唯一；换映射（例如按 K4 顶点邻接）会不会破坏 $q_{\rm live}=0.0696$。

## 复跑

```bash
cd ZeroGPU/src
/usr/bin/python3 zero_cprime_12ch.py
```

| 文件 | 内容 |
|:--|:--|
| `ZeroGPU/src/zero_cprime_12ch.py` | C′ → cycle16 seed 的 12 通道耦合 |
| `ZeroGPU/results/zero_cprime_12ch.json` | baseline/coupled 的 $q$、$q_{\rm live}$、$q_{\rm hist}$ |
| `ZeroGPU/results/zero_cprime_12ch.png` | content $q$ 与 $q_{\rm live}$ |

---

*ZeroGPU · C′ 接 12 通道 · 2026-10-18*
