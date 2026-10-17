# wip/ — 未接入构建的 Lean 草稿与探针（2026-10-14 审计清理后保留）

这里的文件**都不在 `lake` 构建里**（`Zero.lean` 不 import 它们），所以 CI 从未覆盖。
按「保内容、去伪」的处置：留在这里待重写或留作记录，**不接入构建**——以免「磁盘上有文件」被误读成「验证过」。

## 文件清单（19 个）

| 文件 | 内容 | 情况 / 重写要点 |
|:--|:--|:--|
| `Z022_Saturation.lean` | 饱和判据 $c_T=\binom{T}{\lfloor T/2\rfloor}\le2(T+1)\iff T\le5$ | 命题真；用了不存在的 `Nat.choose_le_choose` 参数名。重写：先 `#check Nat.choose_le_choose` 拿签名，再补 $T\ge6\Rightarrow\binom{T}{T/2}\ge\binom63=20>2(T+1)$ |
| `Z047_Involution.lean` | 有限集上无不动点对合 ⟹ 基数为偶 | 命题真；归纳（挖掉 $\{x,\sigma x\}$ 后递归）的目标管理有类型错。重写：先证 `s.card = s'.card + 2`（`s' = s.erase x |>.erase (σ x)`），再递归 |
| `Z002_Dimension.lean` | Z002 §6.7 的维数账：$\dim\ker B=\dim\operatorname{im}B\iff m=4$（$B$ 取 1-骨架关联映射；$\ker B$＝圈空间，$\operatorname{im}B$＝零和超平面） | 未接入构建；探针见下一个文件 |
| `Z002_Dimension_probe.lean` | 同一命题的探针版（Mathlib 矩阵／线性映射记法试验） | 探针，不是证明 |
| `_orphan_Z097_QNumbers.lean` | Z097 的量子数关系（$Q=T_3+Y/2$、$Y_R=Y_L+2T_3$，纯有理数算术） | 孤儿文件（未被 import、未审计）；若要保留，改名 `Z097_QNumbers.lean` 并接入 `Zero.lean` |
| `probe2.lean` … `probe15.lean`（14 个） | Mathlib API／记法的探索性片段（`#check`、类型试验、`example` 片段） | 探针，非命题、不接入构建 |

**已修好的先例**：`Z052`（$b$ 偶 ⟹ $\gcd(b+1,b-1)=1$）2026-10-14 重写并接入构建（现为 `Zero/Z052_Gcd.lean`，定理 `z052_even_b_gcd`）——走「公因数整除两数之差 2，而 $b-1$ 奇 ⟹ 公因数只能是 1」，手工构造整除避开 `dvd_sub` 的 API 细节。
