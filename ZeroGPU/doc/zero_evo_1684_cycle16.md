# 16/8/4 固定底层接 cycle16：也能跑

## 结论

把 16/8/4 固定底层接到 **cycle16**（16 顶点）后，6/6 判据全过：

| 量 | $\alpha=0$ | $\alpha=0.1$ | $\alpha=0.3$ |
|:--|--:|--:|--:|
| content $q$ | $0.0625$ | $0.1123$ | $0.1123$ |
| $n_{std}$ | $0.0248$ | $0.8926$ | $0.8926$ |
| $q_{16}$ | $0.1123$ | $0.1123$ | $0.1123$ |
| 三层份额 $w_b/w_m/w_t$ | $[0.1956,0.4116,0.3928]$ | 同左 | 同左 |
| content 周期 4 | $0$ | $0.037$ | $0.131$ |
| $q_{16}$ 周期 4 | $0.59$ | $0.59$ | $0.59$ |

- **无 NaN、无发散**；content $n_{std}=0.893<1$；
- $\alpha>0$ 后 content 出现周期 4；
- $q_{16}$ 周期 4 很强（bin $0.59$）；
- 16/8/4 底层自己的 $q_{16}=0.1123$、份额 $[0.196,0.412,0.393]$ 不随 $\alpha$ 变。

![16/8/4 接 cycle16](ZeroGPU/results/zero_evo_1684_cycle16.png)

## 怎么接的

- **macro/geometry/wave** 在 cycle16（16 顶点）；
- **micro** 仍是三份 $W$（底 16／中 8／顶 4），class 空间是 $4\times4=16$；
- 每代把 class 记录按

$$
x=\frac{1}{7}\Bigl(W_b+A_m^{\mathsf T}W_m+A_t^{\mathsf T}W_t\Bigr)
$$

变成 16 槽记录，作为 `gen_closures`；
- $\alpha>0$ 时，`_closure_update` 用这 16 槽记录回写 content；
- 16 个 class 与 16 个顶点按 `site*4+fiber` 一一对应。

## 判据（6/6）

1. 无 NaN：PASS；
2. content 有界（$n_{std}<1$）：PASS（最大 $0.894$）；
3. $q_{16}>0$：PASS（$0.1123$）；
4. 三层份额归一：PASS；
5. $\alpha>0$ 时 content 有周期 4：PASS（$0.037$/$0.131$）；
6. $q_{16}$ 有周期 4：PASS（$0.59$）。

## 边界

- 16/8/4 是**输入**，不是从观测推出来的；$q_{16}=0.1123$、份额 $[0.196,0.412,0.393]$ 是这套输入的后果，**不对目标 $0.0683$/$[0.172,0.276,0.552]$ 拟合**。
- 16 class 与 16 顶点的映射是 `site*4+fiber`，这是一条新识别；换映射要重跑。
- $\alpha\ge0.1$ 后 content $q$ 与 $n_{std}$ 饱和；$\alpha$ 只改相位幅度，不改均值。
- cycle16 的 geometry/wave 仍在；但这次没有单独审它们的谱。
- 无新量纲锚：16/8/4、映射、$\alpha$ 都无量纲；$v$ 仍是唯一锚。

## 下一步

1. **长跑**：5000–10000 代，看 $q_{16}=0.1123$、份额 $[0.196,0.412,0.393]$ 是否漂移；$\alpha=0.1$ 是否长期有界。
2. **换映射**：`site*4+fiber` 是按序映射；换成 K4 邻接/顶点角色映射，看 $q_{16}$/周期 4 是否变。
3. **接 geometry/wave 的谱**：cycle16 的 $h$/$\psi$ 是否被 16/8/4 记录带出周期 4；这是 12 通道的下一层。

## 复跑

```bash
cd ZeroGPU/src
/usr/bin/python3 zero_evo_1684_cycle16.py
```

| 文件 | 内容 |
|:--|:--|
| `ZeroGPU/src/zero_evo_1684_cycle16.py` | 16/8/4 固定底层 + cycle16 macro |
| `ZeroGPU/results/zero_evo_1684_cycle16.json` | 各 $\alpha$ 的 $q$、$q_{16}$、份额、周期 4 |
| `ZeroGPU/results/zero_evo_1684_cycle16.png` | $q$ 与份额随 $\alpha$ |

---

*ZeroGPU · 16/8/4 接 cycle16 · 2026-10-18*
