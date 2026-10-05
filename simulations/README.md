# `simulations/` —— Zero 系列（第一代数值探索）的**可运行副本**

**来源**：本目录的 `zero_sum_*.py` / `_results.json` / `.png` 逐字节复制自
`../modular-equilibrium/simulations/`。

> **本版说明（笔记不在此处）**：7 篇 Zero 笔记的**唯一副本在仓库根目录**
> （已登记进 [`INDEX.md`](../INDEX.md) §3）。本目录原先随包复制了它们的 .md 副本，
> 造成**一份内容两处存放**（译文、去 U 改动、链接前缀都要同步两次，已两次出错），
> 故本轮**删除这 7 个 .md 副本**，只保留程序与结果。需要笔记时请读根目录版本：
> [`closure_exit`](../zero_sum_closure_exit.md)、[`closure_time_selection`](../zero_sum_closure_time_selection.md)、
> [`global_R_local_P`](../zero_sum_global_R_local_P.md)、[`open_reservoir`](../zero_sum_open_reservoir.md)、
> [`periodic_destruction`](../zero_sum_periodic_destruction.md)、[`reproduction_audit`](../zero_sum_reproduction_audit.md)、
> [`tri_layer_universe`](../zero_sum_tri_layer_universe.md)。

**为什么搬进来**：本语料根目录的 `zero_sum_*.md/.py` 是同一批文件的副本，但脚本用

```python
ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "simulations" / "zero_sum_closure_exit_results.json"
```

在原址 `modular-equilibrium/simulations/` 下，`ROOT` 解析为 `modular-equilibrium/`，写入正确；
被复制到 `lh/` **根目录**后，`ROOT` 解析为 `/Users/oygb/Downloads/`（不存在）⇒
**10 个实验在 `lh/` 里全部 `rc=1`（FileNotFoundError）**。

放进 `lh/simulations/` 后 `ROOT` 恰好解析为 `lh/`，输出仍写回本目录，**脚本恢复可跑**：

```bash
cd /Users/oygb/Downloads/lh
python3 simulations/zero_sum_closure_exit.py     # rc=0，写回 simulations/ 下的 png 与 json
```

**这一层对理论的贡献（已在 G 系列登记）**：

| 实验 | 贡献 | 被引次数 |
|:--|:--|--:|
| `zero_sum_closure_exit.py` | **`G34`** 的依据：其"继续规则"对照组显示闭合路径会再次闭合、产生更多闭合事件 ⇒ 由此证明「闭合即退出」是 **π 的定义**，不是独立动力学规则 | 7 |
| `zero_sum_periodic_destruction.py` | **`I8`（闭合再播种）** 的登记依据（`G20` §3 引理 74） | 5 |
| `zero_sum_open_reservoir.py` | 参与 `I8` 登记 | 4 |
| `zero_sum_global_R_local_P.py` | 参与 | 3 |
| `zero_sum_tri_layer_universe.py` | 参与；并自述 *"The next missing structure is locality"* | 2 |
| `closure_time_selection` / `cycle_evolution` / `geometry_probe` / `living_universe` / `reproduction_audit` | **未被任何 G/D 文档引用**（`living_universe_results.json` 3.4 MB，是全语料最大文件，亦未被消费） | **0** |

**等级**：`G20` 明确规定本系列结论为**【数值证据】，不是【导出】**；且**不进入** A5
（「闭合再播种」按"不增公理"登记为账本 **I8**）。

**开放项**：本系列**没有逐篇分诊**（D 系列有 [`G9`](../G9_d_series_reference_triage.md) 的对应物，
Zero 系列没有）。上表 5 个零引用实验应显式标为【纯历史】或补做分诊。
