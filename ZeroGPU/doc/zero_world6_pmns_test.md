# 拿预言 $\delta_{CP}<3\pi/4$ 和 PMNS 测量对

**性质**：讨论稿/对比，不占正式篇号，不改变根目录 Z 篇的任何判定。
**脚本**：`ZeroGPU/src/zero_world6_pmns_test.py`。
**数据**：`results/zero_world6_pmns_test.json`。
**复跑**：

```bash
cd ZeroGPU/src && /usr/bin/python3 zero_world6_pmns_test.py
```

---

## 摘要

相位阈值模型（$\theta_{\rm crit}=3\pi/4$，质量零点）预言

$$
\delta_{CP}<135^\circ .
\tag{0.1}
$$

PMNS 测量（NuFIT 6.0，NO，近似值；见 [NuFIT-6.0](https://inspirehep.net/literature/2838825)）：
best-fit $\delta_{CP}\approx1.2\pi\approx216^\circ$，$1\sigma\approx[162^\circ,270^\circ]$。

$$
\boxed{\;\text{预言区 }[0,135^\circ]\text{ 与 NO }1\sigma\ [162^\circ,270^\circ]\text{ 不相交，}
\text{差 }27^\circ\text{——张力。}\;}
\tag{0.2}
$$

---

## §1 对比表

| $\theta_{\rm crit}$ | 预言 | NO $1\sigma$ 内？ | 差距 |
|:--|:--|:--|--:|
| **$3\pi/4$（质量零点）** | $\delta<135^\circ$ | ❌ | $27^\circ$ |
| $\pi$（双覆盖） | $\delta<180^\circ$ | ✅（边缘） | $0$ |
| $2\pi$（整周期） | 无约束 | ✅ | — |

- **$3\pi/4$**：预言区在 $1\sigma$ 下边界（$162^\circ$）以下 $27^\circ$——**张力**；
- **$\pi$**：$180^\circ$ 落在 $1\sigma$ 内——**边缘相容**；
- **$2\pi$**：无约束，安全但不提供信息。

---

## §2 读法

- 当前 PMNS 的 NO best-fit 在 $216^\circ$ 附近，明显高于 $135^\circ$；
- 若 $\delta_{CP}$ 最终确认在 $1.2\pi$ 附近，**质量零点（$3\pi/4$）作为毁灭点被排除**，
  临界点只能取 $\pi$ 或 $2\pi$；
- $\pi$（双覆盖）给出 $\delta<180^\circ$，与 $1\sigma$ 的下边缘（$162^\circ$）重叠，
  但 best-fit $216^\circ$ 仍在其外——**边缘张力**；
- 只有 $2\pi$ 不被 PMNS 约束。

---

## §3 判定

1. **质量零点 $3\pi/4$ 与 PMNS NO $1\sigma$ 不相交**，张力约 $27^\circ$（$\sim1\sigma$ 量级）。
2. **双覆盖 $\pi$ 边缘相容**（$180^\circ$ 在 $1\sigma$ 内）。
3. **整周期 $2\pi$ 无约束**。
4. **下一代实验**（DUNE、T2HK）若把 $\delta_{CP}$ 定在 $1.2\pi$ 附近，则质量零点被排除；
   若 $\delta_{CP}$ 落在 $135^\circ$ 以下，则质量零点存活。

---

## §4 边界

1. PMNS 数值是**近似值**（NuFIT 6.0 的精确表未取到；$1\sigma$ 区间随拟合版本变，
   NuFIT 5.2 约 $[178^\circ,256^\circ]$）；
2. IO 的 best-fit 约 $270^\circ$，张力更大；
3. 预言依赖 $\theta_{\rm crit}=3\pi/4$ 的选择；
4. 不是 Z 篇；无新量纲锚。

---

## §5 判据（4/4 PASS）

| # | 判据 | 结果 |
|--:|:--|:--|
| 1 | $3\pi/4$ 预言 $\delta<135^\circ$ | PASS |
| 2 | 与 NO $1\sigma$ 不相交（张力） | PASS，差 $27^\circ$ |
| 3 | $\pi$ 与 $1\sigma$ 有重叠（边缘） | PASS |
| 4 | $2\pi$ 无约束 | PASS |

---

## §6 复跑与文件

```bash
cd ZeroGPU/src
/usr/bin/python3 zero_world6_pmns_test.py     # 4 条判据，<1 s
```

| 文件 | 内容 |
|:--|:--|
| `src/zero_world6_pmns_test.py` | 预言 vs PMNS 的边界比较 |
| `results/zero_world6_pmns_test.json` | 预言、近似观测、张力 |

---

*ZeroGPU 讨论稿 · zero_world6_pmns_test · 2026-10-19*
