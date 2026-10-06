# 旧仿真审计：这些 `zero_sum_*` 到底在演化什么

审计范围与已核实的事实：

- [`simulations/`](simulations) 共 **32 项** = 10 个 `.py` + 10 个 `_results.json` + 10 个 `.png` + [`README.md`](simulations/README.md) + `__pycache__`。
- [`simulations/README.md:3`](simulations/README.md#L3) 自述这些文件"逐字节复制自 `../modular-equilibrium/simulations/`"。**已用 `cmp` 逐文件核实**：10 个 `.py`、10 个 `_results.json`、9 个 `.png` 全部 `SAME`，只有 [`zero_sum_global_R_local_P.png`](simulations/zero_sum_global_R_local_P.png) 不同（该目录的 json/png 时间戳为 Oct 6 09:50，被重跑过一次）。
- [`/Users/oygb/Downloads/modular-equilibrium/simulations/`](../modular-equilibrium/simulations/) 多出 **3 个 .py**（`d186_random_loop_storage.py`、`d187_zero_quotient_differentiation.py`、`zero_generative_selection.py`）与 **7 个 .md**（就是被 [`README.md:9`](simulations/README.md#L9) 说明"已删除、唯一副本在仓库根目录"的那 7 篇）。
- [`verify/`](verify) 53 项全部是 `d210`–`d259` 的 **静态结构/代数核验脚本**；按 `lifetime|persist|survival|autocorrel|aging|metastab|glass` 全文检索只命中 1 个文件 [`d222_stratified_destruction_and_local_memory.py`](verify/d222_stratified_destruction_and_local_memory.py)，**无任何时间演化仿真**。

---

## 一、`simulations/` 下每个 `.py` 的一句话用途

用途全部取自文件头 docstring，逐行给出处。

| 文件 | 一句话用途（docstring 出处） | 有 `for` 推进时间？ |
|:--|:--|:--|
| [zero_sum_closure_exit.py:3-29](simulations/zero_sum_closure_exit.py#L3-L29) | 在最小离散模型里对比"闭合后**继续**留在活动层"与"闭合即**退出**活动层（写入 P_i 与全局 Z）"两条规则 | ✅ `for epoch in range(EPOCHS+1)` [`:139`](simulations/zero_sum_closure_exit.py#L139) |
| [zero_sum_closure_time_selection.py:3-33](simulations/zero_sum_closure_time_selection.py#L3-L33) | 对比 `persist_only`（w→w）与 `closure_daughter`（w→w+w），检验"更快闭合 ⇒ 更多周期 ⇒ 更多后代"的选择 | ✅ `for step in range(STEPS+1)` [`:145`](simulations/zero_sum_closure_time_selection.py#L145) |
| [zero_sum_cycle_evolution.py:3-48](simulations/zero_sum_cycle_evolution.py#L3-L48) | 在三种闭合代数 `sterile M=1`／`single_cut M=1+r`／`all_cuts M=2^r` 下算谱系指数 λ=log(M)/T，并带一个保荷确定性突变核 | ✅ `for step in range(STEPS+1)` [`:270`](simulations/zero_sum_cycle_evolution.py#L270) |
| [zero_sum_geometry_probe.py:3-32](simulations/zero_sum_geometry_probe.py#L3-L32) | 只问一件事：零和闭环图是否已选定一个稳定的低维几何（邻域增长／直径／热核谱维数／4 维嵌入应力） | ❌ 静态图 + 热核，无轨迹 |
| [zero_sum_global_R_local_P.py:3-33](simulations/zero_sum_global_R_local_P.py#L3-L33) | 检验"R 全局、P 局部"：三场景 `local_P`／`global_P`／`global_R_reconstruct` 只差"站点如何重建活动库" | ✅ `for epoch in range(EPOCHS+1)` [`:194`](simulations/zero_sum_global_R_local_P.py#L194) |
| [zero_sum_living_universe.py:3-32](simulations/zero_sum_living_universe.py#L3-L32) | 把分支湮灭过程 + RAF 型闭合 + Eigen 准种/误差阈值 + 三元基因组 (closure, fidelity, speed) 压成一个最小闭环 | ✅ `for step in range(steps+1)`，`dt=0.025` [`:220`](simulations/zero_sum_living_universe.py#L220)、[`:195`](simulations/zero_sum_living_universe.py#L195) |
| [zero_sum_open_reservoir.py:3-19](simulations/zero_sum_open_reservoir.py#L3-L19) | 持续开放储层：状态是 Z^d 上的整数向量，只有回到零向量才退出；**无有限寿命、无重播种** | ✅ `for epoch in range(EPOCHS)` [`:115`](simulations/zero_sum_open_reservoir.py#L115) |
| [zero_sum_periodic_destruction.py:3-29](simulations/zero_sum_periodic_destruction.py#L3-L29) | 检验"E 演化 → 闭合入 P/Z → 整个 E 定期毁灭 → D 记录未闭合者 → E 从 P 重生"的循环，三规则 `continuous`/`wipe_no_reseed`/`wipe_reseed` | ✅ `for time in range(1, DESTRUCTION_PERIOD*CYCLES+1)` [`:193`](simulations/zero_sum_periodic_destruction.py#L193) |
| [zero_sum_reproduction_audit.py:3-35](simulations/zero_sum_reproduction_audit.py#L3-L35) | 审计：把此前混为一谈的 closure／persistence／splitting／copying **四者分开**，建周期 ≤12 的精确转移矩阵 | ⚠️ 有代循环 `population = matrix @ population` [`:206-208`](simulations/zero_sum_reproduction_audit.py#L206-L208)，但是纯矩阵幂迭代 |
| [zero_sum_tri_layer_universe.py:3-37](simulations/zero_sum_tri_layer_universe.py#L3-L37) | R/P/E/D 四层设计的最小实现：层由一条闭合历史 w 播种，在寿命 L 内生成全部平衡扩展 e，寿命到则剩余路径入 D | ✅ `for time in range(STEPS+1)` [`:198`](simulations/zero_sum_tri_layer_universe.py#L198) |

---

## 二、真正在跑时间演化的那些

**跑时间演化的：8 个**（closure_exit、closure_time_selection、cycle_evolution、global_R_local_P、living_universe、open_reservoir、periodic_destruction、tri_layer_universe），外加 1 个"矩阵幂迭代式"的 reproduction_audit、1 个静态的 geometry_probe。

**全目录共同事实（三条，可直接引用）**

1. **除 `zero_sum_living_universe.py` 外，没有任何一个仿真使用随机数。** 对 `simulations/*.py` 检索 `random|rng|np.random|poisson|binomial|multinomial|choice`，命中全部落在 `living_universe`（[`:198`](simulations/zero_sum_living_universe.py#L198)、[`:251`](simulations/zero_sum_living_universe.py#L251)、[`:258`](simulations/zero_sum_living_universe.py#L258)、[`:268`](simulations/zero_sum_living_universe.py#L268)、[`:275`](simulations/zero_sum_living_universe.py#L275)、[`:299-300`](simulations/zero_sum_living_universe.py#L299-L300)、[`:313`](simulations/zero_sum_living_universe.py#L313)、[`:318`](simulations/zero_sum_living_universe.py#L318)）。`geometry_probe` 唯一的"随机"是固定 seed 的对照随机正则图（[`:136`](simulations/zero_sum_geometry_probe.py#L136)、[`:140`](simulations/zero_sum_geometry_probe.py#L140)），可复现。其余 8 个的 docstring 明确写：`"No probabilities or mutation rates are used. Every allowed next step is instantiated."`（[closure_exit:25-26](simulations/zero_sum_closure_exit.py#L25-L26)、[open_reservoir:16](simulations/zero_sum_open_reservoir.py#L16)、[global_R_local_P:30](simulations/zero_sum_global_R_local_P.py#L30)、[periodic_destruction:26](simulations/zero_sum_periodic_destruction.py#L26)、[reproduction_audit:23](simulations/zero_sum_reproduction_audit.py#L23)、[tri_layer_universe:21-22](simulations/zero_sum_tri_layer_universe.py#L21-L22)、[closure_time_selection:13-14](simulations/zero_sum_closure_time_selection.py#L13-L14)、[cycle_evolution:26](simulations/zero_sum_cycle_evolution.py#L26)）。
2. **"记录权重"在旧仿真里不存在。** 没有任何仿真把权重写进账本。唯一出现的 `weights` 都在 `living_universe` 里作为**抽样概率向量**（[`:136-138`](simulations/zero_sum_living_universe.py#L136-L138) 突变邻居权重、[`:165-174`](simulations/zero_sum_living_universe.py#L165-L174) de-novo 闭合偏置权重、[`:354`](simulations/zero_sum_living_universe.py#L354) 末代性状均值权重）。`periodic_destruction`/`reproduction_audit`/`tri_layer_universe` 显式登记 `"weights": "none"`（[reproduction_audit_results.json:416](simulations/zero_sum_reproduction_audit_results.json#L416)）。
3. **旧仿真没有测量过关联时间、自关联函数、或任何寿命分布。** 对 `simulations/` 检索 `autocorrel|auto_correl|correlation_time|relaxation|tau_` **零命中**。出现的"寿命"只有：
   - `LAYER_LIFETIME = 4` —— 是**硬编码输入常数**，不是测量输出（[zero_sum_closure_exit.py:61](simulations/zero_sum_closure_exit.py#L61) → [json:14](simulations/zero_sum_closure_exit_results.json#L14)；[zero_sum_global_R_local_P.py:67](simulations/zero_sum_global_R_local_P.py#L67) → [json:12](simulations/zero_sum_global_R_local_P_results.json#L12)）；`tri_layer_universe` 则把 L 当扫描参数 `SCENARIOS = (0, 4, 6)`（[zero_sum_tri_layer_universe.py:72](simulations/zero_sum_tri_layer_universe.py#L72)）。
   - `mean_closure_time` —— "平均闭合时长"，是被选择/被统计的对象，不是关联时间（[closure_time_selection.py:154-156](simulations/zero_sum_closure_time_selection.py#L154-L156)；[tri_layer_universe.py:256-260](simulations/zero_sum_tri_layer_universe.py#L256-L260)）。
   - `first_closure_time`/`first_birth_time`/`alive`/`survival_probability` —— 只有 `living_universe` 有，是**灭亡判据**（阈值式），不是寿命分布（[`:387-388`](simulations/zero_sum_living_universe.py#L387-L388)、[`:403-409`](simulations/zero_sum_living_universe.py#L403-L409)、[`:428-429`](simulations/zero_sum_living_universe.py#L428-L429)）。

### 2.1 逐个演化仿真的四问

**① [zero_sum_living_universe.py](simulations/zero_sum_living_universe.py) —— 唯一带随机抽样的**

- **状态变量**：`open_pairs`（未闭合的成对涨落数，整数）+ `counts`（5³=125 个三元基因型的闭环计数，`GENOTYPES = tuple(product(range(1,6), repeat=3))` [`:58`](simulations/zero_sum_living_universe.py#L58)）。
- **更新规则**（每步 `dt=0.025`，共 400 步到 `t_end=14.0`）：**含随机抽样**。源项 `rng.poisson(pair_source*dt)` [`:251`](simulations/zero_sum_living_universe.py#L251)；湮灭 `rng.poisson(annihilation*open_pairs²*dt)` [`:258`](simulations/zero_sum_living_universe.py#L258)；闭合 `rng.poisson(closure_rate*open_pairs²*dt)` 后按 `DE_NOVO` 权重 `rng.choice` [`:268-280`](simulations/zero_sum_living_universe.py#L268-L280)；复制/死亡 `rng.poisson(birth_mean)`/`rng.poisson(death_mean)` [`:296-300`](simulations/zero_sum_living_universe.py#L296-L300)；忠实复制 `rng.binomial`、突变去向 `rng.multinomial` [`:313-320`](simulations/zero_sum_living_universe.py#L313-L320)。**"记录权重"：无**（权重只用于抽样）。
- **寿命/关联测量**：只有灭亡/自持判据（`total_loops ≥ 100 且 major lineages ≥ 2 且 post_cutoff_birth_events > 0`，[`:396-409`](simulations/zero_sum_living_universe.py#L396-L409)），以及首闭合/首出生时间。**无关联函数、无关联时间。**
- **定量结果在哪**：只在同目录 JSON，**没有对应的 .md**（该文件同时被 [`simulations/README.md:42`](simulations/README.md#L42) 列为"未被任何 G/D 文档引用"的 5 个零引用实验之一；其 3.4 MB 结果也是全语料最大文件）。
  - model：`replicates: 20` [json:16](simulations/zero_sum_living_universe_results.json#L16)、`time_end: 10.0` [json:17](simulations/zero_sum_living_universe_results.json#L17)、`source_cutoff: 7.0` [json:18](simulations/zero_sum_living_universe_results.json#L18)、`sample_size_stop: 250000` [json:21](simulations/zero_sum_living_universe_results.json#L21)
  - `annihilation_only`：`survival_probability = 0.0` [json:28](simulations/zero_sum_living_universe_results.json#L28)、`mean_final_loops = 0.0` [json:30](simulations/zero_sum_living_universe_results.json#L30)
  - `sterile_closure`：`survival_probability = 0.0` [json:38](simulations/zero_sum_living_universe_results.json#L38)，但 `mean_final_loops = 98.65` [json:40](simulations/zero_sum_living_universe_results.json#L40)
  - `high_error`：`1.0` / `485.4` / `mean_lineages 83.85` [json:48,50,52](simulations/zero_sum_living_universe_results.json#L48)
  - `living`：`1.0` / `572.15` / `71.5` [json:58,60,62](simulations/zero_sum_living_universe_results.json#L58)
  - `high_fidelity`：`1.0` / `526.7` / `57.2` [json:68,70,72](simulations/zero_sum_living_universe_results.json#L68)

**② [zero_sum_closure_exit.py](simulations/zero_sum_closure_exit.py) —— 确定性全分支**

- **状态**：每站点一个"活动词集合"`active`（±词，8 站点）、历史层 `histories`、全局零层 `zero_layer`、死层 `dead`（[`:130-136`](simulations/zero_sum_closure_exit.py#L130-L136)）。
- **更新规则**：**无随机**。每个词若 `len(word) >= LAYER_LIFETIME(4)` 则入 D；否则同时生成 `word+(+1,)` 与 `word+(-1,)` 两个子词，闭合的写历史+零类，并按场景决定是否留在活动层（[`:172-190`](simulations/zero_sum_closure_exit.py#L172-L190)）。**记录权重：无**。
- **寿命测量**：无关联测量。只有硬编码 `L=4` 和"寿命到后活动层为空"这一布尔判据。
- **结果**：
  - `continue_after_closure` 末态：`final_active_paths = 0` [json:398](simulations/zero_sum_closure_exit_results.json#L398)、`final_history_records = 12` [json:400](simulations/zero_sum_closure_exit_results.json#L400)、`final_zero_events = 12` [json:401](simulations/zero_sum_closure_exit_results.json#L401)、`final_dead_paths = 32` [json:403](simulations/zero_sum_closure_exit_results.json#L403)
  - `exit_after_closure` 末态：`0` / `8` / `8` / `18` [json:774,776,777,779](simulations/zero_sum_closure_exit_results.json#L774)
  - 5 条 verification 全 `true`，含 `"exit_leaves_no_active_paths_after_lifetime": true` [json:787](simulations/zero_sum_closure_exit_results.json#L787)
  - 逐代活动层序列（[zero_sum_closure_exit.md:24-33](zero_sum_closure_exit.md#L24-L33)）：继续 `8,16,16,16,0,0,0`；退出 `8,12,10,6,0,0,0`
  - 结论段 [zero_sum_closure_exit.md:69](zero_sum_closure_exit.md#L69)：*"这里活动层最终清空，是因为本程序另行设定了有限寿命 L=4，并且没有持续开放储层或重新播种。"*；"尚未证明"5 条 [`:83-87`](zero_sum_closure_exit.md#L83-L87)

**③ [zero_sum_periodic_destruction.py](simulations/zero_sum_periodic_destruction.py) —— 毁灭/重生循环**

- **状态**：`active`（活动路径集合）、`histories`（P_i）、`zero_layer`（Z）、`dead`（D_i）（[`:178-179`](simulations/zero_sum_periodic_destruction.py#L178-L179)）。
- **更新规则**：**无随机**。每步 `evolve_one_step`；当 `time % DESTRUCTION_PERIOD(3) == 0` 时把全部活动路径移入 `dead` 并清空 `active`；若 `scenario.reseed` 则 `reseed_from_history` 用每条闭合历史 w 生成 `w+`/`w-` 两条开放种子（[`:193-243`](simulations/zero_sum_periodic_destruction.py#L193-L243)）。**记录权重：无**。
- **寿命测量**：`DESTRUCTION_PERIOD = 3`、`CYCLES = 4` 是输入（[`:60-61`](simulations/zero_sum_periodic_destruction.py#L60-L61)），非测得。
- **结果**：
  - `continuous`：`final_active_paths = 10296` [json:690](simulations/zero_sum_periodic_destruction_results.json#L690)、`history_records = 784` [json:691](simulations/zero_sum_periodic_destruction_results.json#L691)、`zero_modes = 286` [json:692](simulations/zero_sum_periodic_destruction_results.json#L692)、`dead = 0` [json:693](simulations/zero_sum_periodic_destruction_results.json#L693)
  - `wipe_no_reseed`：`0` / `8` / `3` / `32` [json:1396-1399](simulations/zero_sum_periodic_destruction_results.json#L1396)
  - `wipe_reseed`：`1360` / `680` / `210` / `1040` [json:2102-2105](simulations/zero_sum_periodic_destruction_results.json#L2102)
  - 15 条 verification 全 `true` [json:2142-2157](simulations/zero_sum_periodic_destruction_results.json#L2142-L2157)
  - 周期边界表（[zero_sum_periodic_destruction.md:64-71](zero_sum_periodic_destruction.md#L64-L71)）：重播规则 `t=3→16, 6→80, 9→336, 12→1360`

**④ [zero_sum_tri_layer_universe.py](simulations/zero_sum_tri_layer_universe.py) —— R/P/E/D**

- **状态**：`future`（未来各时刻将出生的层）、`layer_deaths`、`path_deaths`、`result`(R)、`history`(P)（[`:190-196`](simulations/zero_sum_tri_layer_universe.py#L190-L196)）。
- **更新规则**：**无随机**。每个出生层在寿命 L 内生成所有 `τ ≤ L` 的平衡扩展（[`:270-286`](simulations/zero_sum_tri_layer_universe.py#L270-L286)）；`layer_deaths[time+L] += layer_count`、`path_deaths[time+L] += layer_count * 2**L`。**记录权重：无。**
- **寿命测量**：`lifetime` 是扫描输入 `SCENARIOS = (0, 4, 6)`（[`:72`](simulations/zero_sum_tri_layer_universe.py#L72)）；测的是"寿命为 L 时的末态计数"，不是寿命分布。
- **结果**：
  - L=4：`final_active_layers = 564.0` [json:329](simulations/zero_sum_tri_layer_universe_results.json#L329)、`final_result_records = 608.0` [json:330](simulations/zero_sum_tri_layer_universe_results.json#L330)、`final_result_modes = 39` [json:331](simulations/zero_sum_tri_layer_universe_results.json#L331)、`final_dead_paths = 720.0` [json:336](simulations/zero_sum_tri_layer_universe_results.json#L336)、`final_mean_closure_time = 2.888` [json:337](simulations/zero_sum_tri_layer_universe_results.json#L337)
  - L=6：`1176.0` / `1188.0` / `43` / `832.0` / `3.532` [json:636-644](simulations/zero_sum_tri_layer_universe_results.json#L636)
  - L=0：`final_active_layers = 1.0`、`result_records = 0.0` [json:168-169](simulations/zero_sum_tri_layer_universe_results.json#L168)
  - 精确计数式 B_L、D_L=2^L（[zero_sum_tri_layer_universe.md:74-86](zero_sum_tri_layer_universe.md#L74-L86)）；数值表 [`:95-97`](zero_sum_tri_layer_universe.md#L95-L97)

**⑤ [zero_sum_cycle_evolution.py](simulations/zero_sum_cycle_evolution.py) —— 谱系增长 + 确定性突变**

- **状态**：按模式索引的 **cohort 数组**（每个模式一个长度为周期的数组，追踪"词内位置"），不是空间态（[`:261-265`](simulations/zero_sum_cycle_evolution.py#L261-L265)）。
- **更新规则**：**无随机**。`next_cohorts[index][1:] += current[:-1]` 推进，`finishers = current[-1]` 乘以 `multiplicity` 后按 `scenario.mutation`（0 或 0.03）**确定性比例拆分**给突变目标（[`:332-352`](simulations/zero_sum_cycle_evolution.py#L332-L352)）。总体超过 `POPULATION_SCALE_LIMIT=1e6` 时**整体重标定**（[`:356-361`](simulations/zero_sum_cycle_evolution.py#L356-L361)）。**记录权重：无**。
- **寿命测量**：无。有的是 `mean_return_density`（=r/T）与 `mean_lambda`。
- **结果**（summary [json:22-83](simulations/zero_sum_cycle_evolution_results.json#L22-L83)）：`sterile_regular` `final_population = 1.0`、`top_lambda = 0.0` [json:28,32](simulations/zero_sum_cycle_evolution_results.json#L28)；`one_cut_clone` `final_population = 1000000.0`、`top_lambda = 0.18310`、`rescale_steps = 38` [json:43,50,51](simulations/zero_sum_cycle_evolution_results.json#L43)；`one_cut_evolution` `final_lineages = 122`、`final_diversity = 0.982`、`top_fraction = 0.6867`、`rescale_steps = 137` [json:59,60,64,66](simulations/zero_sum_cycle_evolution_results.json#L59)；`all_cuts_evolution` `final_lineages = 122`、`final_diversity = 1.2022`、`top_word = "-+-+-+-+-+"`、`top_fraction = 0.5298`、`rescale_steps = 142` [json:74,75,78,79,81](simulations/zero_sum_cycle_evolution_results.json#L74)。**注意 `rescale_steps` 意味着"种群数"是数值重标定的产物。**

**⑥ [zero_sum_closure_time_selection.py](simulations/zero_sum_closure_time_selection.py) —— 闭合时间选择**

- **状态**：全模式 cohort 数组（初始每个模式各 1，[`:136-139`](simulations/zero_sum_closure_time_selection.py#L136-L139)）。
- **更新规则**：**无随机**。闭合后 `next_cohorts[index][0] += finishers`，若 `daughters_per_closure=1` 再加一份（w→w+w）（[`:210-220`](simulations/zero_sum_closure_time_selection.py#L210-L220)）。**记录权重：无。**
- **寿命测量**：测的是 `mean_closure_time`（被选择的闭合时长），不是关联时间。
- **结果**：`persist_only` `final_population = 3016.0`、`final_lineages = 3016`、`mean_closure_time = 17.107` [json:22,23,25](simulations/zero_sum_closure_time_selection_results.json#L22)；`closure_daughter` `final_population = 1e9`、`final_diversity = 2.04e-08`、`mean_closure_time = 2.0000000019`、`top_word = "-+"`、`top_fraction = 0.99999999907`、`rescale_steps = 40` [json:34,36,37,39,40,41](simulations/zero_sum_closure_time_selection_results.json#L34)。公式 λ=ln2/T（[zero_sum_closure_time_selection.md:56](zero_sum_closure_time_selection.md#L56)）。

**⑦ [zero_sum_open_reservoir.py](simulations/zero_sum_open_reservoir.py) —— 无寿命持续**

- **状态**：`Counter` 计数的活动状态 `q ∈ Z^d`（[`:98`](simulations/zero_sum_open_reservoir.py#L98)）。
- **更新规则**：**无随机**。每个状态的全部后继（某个坐标 ±1），回零者计入 `closures` 并退出，其余进入 `next_active`（[`:115-134`](simulations/zero_sum_open_reservoir.py#L115-L134)）。**记录权重：无。**
- **寿命测量**：`"layer_lifetime": null`、`"reseeding": null` [json:6-7](simulations/zero_sum_open_reservoir_results.json#L6-L7) —— 显式无寿命。
- **结果**：d=1 `final_active_paths = 504`、`cumulative_closure_fraction = 0.0419708029` [json:135,137](simulations/zero_sum_open_reservoir_results.json#L135)；d=2 `1214752` / `0.0073788330` [json:258,260](simulations/zero_sum_open_reservoir_results.json#L258)；d=4 `1777511808` / `0.0005702916` [json:381,383](simulations/zero_sum_open_reservoir_results.json#L381)。verification 5 条全 true [json:388-394](simulations/zero_sum_open_reservoir_results.json#L388-L394)。

**⑧ [zero_sum_global_R_local_P.py](simulations/zero_sum_global_R_local_P.py) —— 全局化 vs 局域性**

- **状态**：每站点的 `local_histories`、`active_histories`，加全局 `global_results`(R) 与 `global_histories`（[`:189-191`](simulations/zero_sum_global_R_local_P.py#L189-L191)）。
- **更新规则**：**无随机**。每 epoch 每个活动历史 × 全部扩展 e 生成 `history+e`，写 local P、全局 R；再按场景从 `local_P`／`global_P`／`global_R_reconstruct` 重建活动库（[`:221-247`](simulations/zero_sum_global_R_local_P.py#L221-L247)）。**记录权重：无。**
- **寿命测量**：无；`LAYER_LIFETIME = 4` 是输入（[`:67`](simulations/zero_sum_global_R_local_P.py#L67)）。
- **结果**：`local_P` `final_global_results = 10343`、`final_local_history_overlap = 0.01091`、`final_distinct_top_results = 6` [json:160,162,164](simulations/zero_sum_global_R_local_P_results.json#L160)；`global_P` `10343` / `0.99955` / `1` [json:317,319,321](simulations/zero_sum_global_R_local_P_results.json#L317)；`global_R_reconstruct` `8846` / `0.99901` / `1` [json:474,476,478](simulations/zero_sum_global_R_local_P_results.json#L474)。

**⑨ [zero_sum_reproduction_audit.py](simulations/zero_sum_reproduction_audit.py) —— 唯一给出"灭绝时间"的**

- **状态**：模式集上的整数种群向量（[`:203-204`](simulations/zero_sum_reproduction_audit.py#L203-L204)）。
- **更新规则**：**无随机**。`population = matrix @ population`，18 代（[`:206-208`](simulations/zero_sum_reproduction_audit.py#L206-L208)）。**记录权重：无。**
- **"寿命"测量**：**这是全目录唯一测"灭绝时间"的地方**——`nilpotent_index(T)`（[`:212-220`](simulations/zero_sum_reproduction_audit.py#L212-L220)），即谱系归零所需的代数。
- **结果**（matrix_audit [json:16-71](simulations/zero_sum_reproduction_audit_results.json#L16-L71)）：
  - `persist`：`spectral_radius = 1.0` [json:18](simulations/zero_sum_reproduction_audit_results.json#L18)，序列 `constant` [json:22-25](simulations/zero_sum_reproduction_audit_results.json#L22)
  - `copy`：`spectral_radius = 2.0` [json:29](simulations/zero_sum_reproduction_audit_results.json#L29)，序列 `exponential` [json:32-35](simulations/zero_sum_reproduction_audit_results.json#L32)
  - `split_parent`：`1.0` / `unipotent_nilpotent_index = 6` [json:39,41](simulations/zero_sum_reproduction_audit_results.json#L39)，`polynomial/subexponential` [json:44-46](simulations/zero_sum_reproduction_audit_results.json#L44)
  - **`split_only`：`spectral_radius = 0.0` [json:50](simulations/zero_sum_reproduction_audit_results.json#L50)、`zero_eigenvalues = 123` [json:51](simulations/zero_sum_reproduction_audit_results.json#L51)、`nilpotent_index = 6` [json:52](simulations/zero_sum_reproduction_audit_results.json#L52)，四个种子全部 `"finite extinction"` [json:54-57](simulations/zero_sum_reproduction_audit_results.json#L54)** ← **灭绝时间 = 6 代**
  - `split_all`：`1.0` / `6` / `polynomial` [json:61,63,66](simulations/zero_sum_reproduction_audit_results.json#L61)
  - interpretation [json:737-742](simulations/zero_sum_reproduction_audit_results.json#L737)："closure alone gives persistence, not reproduction" [json:739](simulations/zero_sum_reproduction_audit_results.json#L739)、"if the parent is removed, the process can die in finite generations" [json:740](simulations/zero_sum_reproduction_audit_results.json#L740)

**⑩ [zero_sum_geometry_probe.py](simulations/zero_sum_geometry_probe.py) —— 不演化，但否掉了几何**

- **状态**：静态图（周期 T 的零和循环词旋转类为节点，交换相邻异号步为边）。
- **更新规则**：无。用**热核** `e^{-tL}`（3 个谱点 0.10/0.05/0.02，[`:79`](simulations/zero_sum_geometry_probe.py#L79)）测谱维数——这是"扩散时间"而非状态演化。
- **结果**：`summary` [json:22](simulations/zero_sum_geometry_probe_results.json#L22)：T=8 `nodes = 10`、`edges = 18`、`mean_degree = 3.6`、`diameter = 4`、`mds_4d_stress = 0.2107` [json:26,27,28,29,45](simulations/zero_sum_geometry_probe_results.json#L26)；T=10 `26` / `4.923` / `0.184` [json:50,52,69](simulations/zero_sum_geometry_probe_results.json#L50)。**interpretation [json:270-273](simulations/zero_sum_geometry_probe_results.json#L270)：`"status": "no stable four-dimensional plateau detected"` [json:271](simulations/zero_sum_geometry_probe_results.json#L271)**、`"next_gap": "a zero-sum closure graph alone does not fix the number of effective geometric directions"` [json:273](simulations/zero_sum_geometry_probe_results.json#L273)。

---

## 三、亚稳畴、长寿命、老化、玻璃态 —— 在旧仿真里**未找到**

### 3.1 在 `simulations/` 与根目录 `zero_sum_*` 里

**未找到。** 具体核实：

- 对 [`simulations/`](simulations) 检索 `metastab|glassy|aging|ageing|domain|autocorrel|correlation_time|relaxation|亚稳|玻璃态|老化|畴` —— **零命中**。
- 旧仿真的时间尺度全部是**硬编码输入**（`LAYER_LIFETIME=4`、`EPOCHS=4/6/10`、`STEPS=10/120/300`、`DESTRUCTION_PERIOD=3`、`CYCLES=4`、`dt=0.025`），没有任何一处是"测出来的寿命"。
- 旧仿真的结论词是**"熄灭/持续/增长"**（extinction / persistence / growth），不是"活多久"。

### 3.2 但 `zerogpu/`（更新的一代，Oct 5–6）明确测了这件事，结论是**否定的**

这一代不在你列的三个目录里，但它正是问题 3 的答案所在，必须报告。

**先看"看起来像亚稳畴/老化"的正面数值** —— [zerogpu/EVO_LAYERS.md:53-68](zerogpu/EVO_LAYERS.md#L53-L68)：

| 观测 | 值域 | 出处 |
|:--|:--|:--|
| 区域数（活动层 L2） | **1 – 36** | [EVO_LAYERS.md:57](zerogpu/EVO_LAYERS.md#L57) |
| **空间自相关 ac1** | **0.19 – 0.66** | [EVO_LAYERS.md:60](zerogpu/EVO_LAYERS.md#L60) |
| 出生/死亡（每窗） | 0–4 / 0–4，**constant 8/9** | [EVO_LAYERS.md:61](zerogpu/EVO_LAYERS.md#L61) |
| 记录速率衰减（"会变老的原生时钟"） | REC `2.27→1.38`、TRAF `1.52→1.03`、BOTH `4.71→1.94` | [EVO_LAYERS.md:49](zerogpu/EVO_LAYERS.md#L49) |
| 畴粗化事件 | "区域数 9→3、平均区域尺寸 77→267" | [EVO_LAYERS.md:96](zerogpu/EVO_LAYERS.md#L96) |

原文 [EVO_LAYERS.md:51](zerogpu/EVO_LAYERS.md#L51)：*"**这是"老化"结构，而且是自发的**（没有任何外加时钟）"*；[`:105`](zerogpu/EVO_LAYERS.md#L105) 把"畴式活动层（L2，空间正自相关 0.19–0.66）"列为"长出来的"之一。

**但同一批文档的复核把它全部推翻** —— [zerogpu/EVAL_can_zero_evolve.md](zerogpu/EVAL_can_zero_evolve.md)：

> [`:18`](zerogpu/EVAL_can_zero_evolve.md#L18)：② **结构/演化层** | ❌ **不是演化** | 最慢非平凡弛豫 $\tau_2\le1.18$ 微观步（$L{=}4$）⟹ **无亚稳畴、无长寿命、无老化、无记忆**

> [`:19`](zerogpu/EVAL_can_zero_evolve.md#L19)：③ §40 曾报的"演化" ❌ 是伪谱 × 二次外推 | $2.64\times10^6$ 步 vs 真上界 $1.18$ 步，差 $2.2\times10^6$ 倍

> [`:69-70`](zerogpu/EVAL_can_zero_evolve.md#L69-L70)：系统**几步之内就平衡了**。没有"畴"可以活过两代，因此 §37/§38 的"亚稳畴＝局域慢模"在真谱下**没有对象**。

> [`:110`](zerogpu/EVAL_can_zero_evolve.md#L110)：有限状态 ⟹ 谱是离散的、间隙是 $\Theta(1)$、**没有指数小间隙** ⟹ **没有亚稳性可言**。

> [`:124`](zerogpu/EVAL_can_zero_evolve.md#L124)：**Zero 在电脑里能"跑"，但按现有条款它不会"演化"**；[`:126`](zerogpu/EVAL_can_zero_evolve.md#L126) 这不是算力问题，是零和层约束的**结构性后果**。

**"三方不可兼得"塌成两方**（[EVAL_can_zero_evolve.md:100-110](zerogpu/EVAL_can_zero_evolve.md#L100-L110)）：忠实 ∧ 守恒 ⟹ 自由弛豫（无结构）；要结构必须有偏核（＝具名输入，违反 Z0③）。

**专门的畴寿命测量是直接否定的** —— [zerogpu/z0_life.py:5-6](zerogpu/z0_life.py#L5-L6)：

> ① 细尺度畴**寿命恒为 1 代**（重排，不是成形/存活/解体）
> ② 宏观（块平均）场自相关 $\tau_{1/e}$ 也只有 **1–6 代**（无时间记忆）

判据要求三条同时成立（[z0_life.py:14](zerogpu/z0_life.py#L14)）：变化 `move>0.05` ／ 记忆 `τ_{1/e}≥3` ／ 结构持久 `畴寿命中位 >3`。

实测 [zerogpu/results/z0_life.json](zerogpu/results/z0_life.json)（7 个相点）：

| 相点 | `f` | `move` | `acf_tau` | `life_med` | `life_mean` | `life_max` |
|:--|--:|--:|--:|--:|--:|--:|
| sync | 1.0 [`:6`](zerogpu/results/z0_life.json#L6) | 0.5948 [`:7`](zerogpu/results/z0_life.json#L7) | 6 [`:9`](zerogpu/results/z0_life.json#L9) | **1.0** [`:12`](zerogpu/results/z0_life.json#L12) | 1.0 [`:11`](zerogpu/results/z0_life.json#L11) | 7 [`:14`](zerogpu/results/z0_life.json#L14) |
| async_kac 0.25 | 0.0154 [`:19`](zerogpu/results/z0_life.json#L19) | **0.0** [`:20`](zerogpu/results/z0_life.json#L20) | 1 [`:22`](zerogpu/results/z0_life.json#L22) | **1.0** [`:25`](zerogpu/results/z0_life.json#L25) | 19.56 [`:24`](zerogpu/results/z0_life.json#L24) | 543 [`:27`](zerogpu/results/z0_life.json#L27) |
| async_kac 0.1 | 0.0391 [`:32`](zerogpu/results/z0_life.json#L32) | 0.0274 [`:33`](zerogpu/results/z0_life.json#L33) | 1 [`:35`](zerogpu/results/z0_life.json#L35) | **1.0** [`:38`](zerogpu/results/z0_life.json#L38) | 1.9 [`:37`](zerogpu/results/z0_life.json#L37) | 215 [`:40`](zerogpu/results/z0_life.json#L40) |
| async_kac 0.05 | 0.0762 [`:45`](zerogpu/results/z0_life.json#L45) | 0.0298 [`:46`](zerogpu/results/z0_life.json#L46) | 1 [`:48`](zerogpu/results/z0_life.json#L48) | **1.0** [`:51`](zerogpu/results/z0_life.json#L51) | 1.17 [`:50`](zerogpu/results/z0_life.json#L50) | 111 [`:53`](zerogpu/results/z0_life.json#L53) |
| async_kac 0.02 | 0.2031 [`:58`](zerogpu/results/z0_life.json#L58) | 0.135 [`:59`](zerogpu/results/z0_life.json#L59) | 1 [`:61`](zerogpu/results/z0_life.json#L61) | **1.0** [`:64`](zerogpu/results/z0_life.json#L64) | 1.02 [`:63`](zerogpu/results/z0_life.json#L63) | 34 [`:66`](zerogpu/results/z0_life.json#L66) |
| async_kac 0.01 | 0.3945 [`:71`](zerogpu/results/z0_life.json#L71) | 0.1388 [`:72`](zerogpu/results/z0_life.json#L72) | 1 [`:74`](zerogpu/results/z0_life.json#L74) | **1.0** [`:77`](zerogpu/results/z0_life.json#L77) | 1.02 [`:76`](zerogpu/results/z0_life.json#L76) | 86 [`:79`](zerogpu/results/z0_life.json#L79) |
| async_event | 0.4074 [`:84`](zerogpu/results/z0_life.json#L84) | 0.0392 [`:85`](zerogpu/results/z0_life.json#L85) | 85 [`:87`](zerogpu/results/z0_life.json#L87) | **1.0** [`:90`](zerogpu/results/z0_life.json#L90) | 1.01 [`:89`](zerogpu/results/z0_life.json#L89) | 7 [`:92`](zerogpu/results/z0_life.json#L92) |

**畴寿命中位数在全部 7 个相点都是 1.0**（唯一的长尾出现在 async_kac/0.25，但那个相点 `move = 0.0`，即系统其实不动）。汇总 [json:96-99](zerogpu/results/z0_life.json#L96-L99)：`pass_ = 3, fail = 0`（按 z0_life 的三条判据口径）。

**"长寿命没有出处"的正面 no-go** —— [zerogpu/z0_entropy_barrier.py:22-27](zerogpu/z0_entropy_barrier.py#L22-L27)：

> **单一原生约束无法产生指数寿命：约束只能把状态空间裁成细管／小集，**
> **而沿细管的弛豫至多是多项式的。**
> 要 $\tau_2\sim e^{cV}$ 必须有**能量壁垒**（＝非均匀率 $\kappa$，即偏好，违反 Z0③）或图上的稀有瓶颈（实测 barbell/lollipop 反而更**快**，§47）。
> ⟹ **"零＋约束"这一支里，长寿命没有出处；第二个要素不是原生的。**

配套：[zerogpu/GENERATIVE_AUDIT.md:80](zerogpu/GENERATIVE_AUDIT.md#L80) —— *"**容量越大，寿命越短**（$\tau_2\sim c/V\cdot L^{-1/3}$）"*；[`:90`](zerogpu/GENERATIVE_AUDIT.md#L90) —— *"任何**结构**（寿命层级、质量谱）**不由这一支产生**"*。

### 3.3 一个容易误认成"长寿命"的东西：回零等待时间的重尾

根目录 H 系列测过一个**重尾、期望发散**的等待时间——但它是**随机游走首次回零的等待时间**，不是亚稳陷阱，必须区分开。

- [H1_closure_events_and_last_closure.md:29](H1_closure_events_and_last_closure.md#L29)：*"「闭合有**某个**特征周期」❌ **不存在** | §3：回零次数 $\sim\sqrt\tau$，$E[T]=\infty$"*
- [H2_next_closure_under_anchor.md:44](H2_next_closure_under_anchor.md#L44)：$P(T>t)=\binom{t}{t/2}/2^{t}\sim\sqrt{2/(\pi t)}$
- [H2:51-53](H2_next_closure_under_anchor.md#L51-L53)：中位数 **2 步**；$E[T]$ **发散**；尾部 $\sim\sqrt{2/(\pi t)}$，**重尾、$E[T]=\infty$**
- [H2:74](H2_next_closure_under_anchor.md#L74)：条件化后 $P(T>t\mid T>t_0)=\sqrt{t_0/t}$，中位 $=4t_0$（[`:15`](H2_next_closure_under_anchor.md#L15)）
- [H1:15](H1_closure_events_and_last_closure.md#L15)：$P(\text{上次闭合在最近 }u\text{ 比例内})=\frac2\pi\arcsin\sqrt u$
- 数值核验：[H2_check.py](H2_check.py) 自述"通过 30 / 不符 0（含 30 万条精确反演抽样）"（[H2:8](H2_next_closure_under_anchor.md#L8)）；[H1_check.py](H1_check.py) "通过 39 / 不符 0（含 2^L 全空间精确枚举与 40 万次精确反演抽样）"（[H1:6](H1_closure_events_and_last_closure.md#L6)）

---

## 四、根目录 `zero_sum_*` 里与"演化/寿命/持续性"最相关的 4 组

根目录 `zero_sum_*` 共 **11 个 .md + 14 个 .py**（约 12 组）。挑出最相关的 4 组：

### 4.1 [`zero_sum_persistence_theorems.md`](zero_sum_persistence_theorems.md) —— 熄灭判据（最相关）

**结论段原文**（[§2 定理 T1，`:44-46`](zero_sum_persistence_theorems.md#L44-L46)）：

> $$
> \ \textbf{给定"不断"（}N_k\not\to0\textbf{ 恒不成立）}:\quad \text{必须"无寿命"}\ \vee\ \text{"有重播种"}\ \text{二择一}。\
> $$

[§0，`:13`](zero_sum_persistence_theorems.md#L13)：*"'零不断乱动'里的**'不断'**是一条硬要求：活动层**不可永久为空**。"*

[§6 诚实边界，`:94`](zero_sum_persistence_theorems.md#L94)：*"判据是**对四个被测模型**的枚举，不是'全部可能持续机制'的定理；程序未测别的补充机制（例如自催化环）。"*

**四模型对照表**（[`:38-42`](zero_sum_persistence_theorems.md#L38-L42)）：

| 模型 | 有限寿命 L | 重播种 | 活动层逐代 | 结局 |
|:--|:--:|:--:|:--|:--|
| `closure_exit` / continue | ✅ | ❌ | $8,16,16,16,\mathbf{0},0,0$ | **熄灭** |
| `closure_exit` / exit | ✅ | ❌ | $8,12,10,6,\mathbf{0},0,0$ | **熄灭** |
| `periodic_destruction` / `wipe_no_reseed` | ✅ | ❌ | $10,18,\mathbf{0},0,0,\dots$ | **熄灭** |
| `periodic_destruction` / `wipe_reseed` | ✅ | ✅ | $10,18,16,16,32,80,80,160,336$ | **持续** |
| `open_reservoir`（d=1,2,4） | ❌ | ❌ | 末值 $504$／$1\,214\,752$／$1\,777\,511\,808$ | **持续** |

**T2 爆炸表**（[`:71-73`](zero_sum_persistence_theorems.md#L71-L73)）：$d=1\to504$、$d=2\to1\,214\,752$、$d=4\to1\,777\,511\,808$；[`:76`](zero_sum_persistence_theorems.md#L76) *"'无寿命'这一支虽能持续，却给不出**有界的活动层**"*。

对应 JSON 关键数字：[open_reservoir json:135/258/381](simulations/zero_sum_open_reservoir_results.json#L135)、[periodic_destruction json:690/1396/2102](simulations/zero_sum_periodic_destruction_results.json#L690)、[closure_exit json:398/774](simulations/zero_sum_closure_exit_results.json#L398)。
核验脚本自述 *"独立实断言 22 / 结论行 0 / 不符 0"*（[`:6`](zero_sum_persistence_theorems.md#L6)），实断言清单在 [zero_sum_persistence_theorems_check.py:6-14](zero_sum_persistence_theorems_check.py#L6-L14)。

### 4.2 [`zero_sum_reproduction_transition_theorems.md`](zero_sum_reproduction_transition_theorems.md) —— 唯一的"灭绝时间"定理

**结论段原文**（[§3，`:72-74`](zero_sum_reproduction_transition_theorems.md#L72-L74)）：

> $$
> \ \text{五条规则中，}\textbf{只有显式复制（copy）给指数增长};\ \text{其余至多常数／多项式。}\
> $$

[§4 定理 T3，`:86-92`](zero_sum_reproduction_transition_theorems.md#L86-L92)：

> 故 $N_k=\mathbf 1^\top T^k v_0$ 在**有限代内归零**——程序把这一指数直接算作 `nilpotent_index`：
> $$
> \ \rho(T)=0\ \Longleftrightarrow\ T\ \text{幂零}\ \Longleftrightarrow\ \text{谱系在有限代内灭绝};\quad\text{灭绝时间}=\text{幂零指数}。
> $$
> **这条是**"闭合不等于永续"**的定量形式**：切分而**不复制**，父一撤，谱系就死。

**五规则全表**（[`:64-70`](zero_sum_reproduction_transition_theorems.md#L64-L70)）：

| 规则 | 谱半径 ρ(T) | 零特征值数 | I−T 幂零指数 | T 幂零指数 | 序列行为 |
|:--|--:|--:|--:|--:|:--|
| `persist` | 1.0 | 0 | 1 | — | **constant** |
| `copy` | **2.0** | 0 | — | — | **exponential** |
| `split_parent` | 1.0 | 0 | 6 | — | constant／polynomial |
| `split_only` | **0.0** | **123** | — | **6** | **finite extinction** |
| `split_all` | 1.0 | 0 | 6 | — | constant／polynomial |

对应 JSON：[reproduction_audit json:18,29,39,41,50,51,52,61](simulations/zero_sum_reproduction_audit_results.json#L18)。
开放项（[`:135`](zero_sum_reproduction_transition_theorems.md#L135)）：*"**多类型自催化环**：程序点名但未测——这可能是'繁殖'的真正形态。"*

### 4.3 [`zero_sum_periodic_destruction.md`](zero_sum_periodic_destruction.md) —— 最干净的"续命"对照

**结论段原文**（[§4，`:115`](zero_sum_periodic_destruction.md#L115)）：

> 不重播时，第一次毁灭后系统就**永久失去活动层**。重播时，活动层在每一个周期后都能**重新长出来**。

（[§6，`:183`](zero_sum_periodic_destruction.md#L183)）：

> 因此本次程序验证的是"周期毁灭、保留记录、重新播种"三者可以构成自洽循环，而不是证明该周期由零和约束自动出现。结构登记为 `R-Z-PERIODIC-WASHOUT-RESEED`。

**周期边界数据**（[`:64-71`](zero_sum_periodic_destruction.md#L64-L71)）：

| 规则 | 时间 | 毁灭后活动层 | 重播路径 | 活动路径 | 历史记录 | 零模式 | 死亡路径 |
|:--|--:|--:|--:|--:|--:|--:|--:|
| 不重播 | 3 | 0 | 0 | 0 | 8 | 3 | 32 |
| 重播 | 3 | 0 | 16 | 16 | 8 | 3 | 32 |
| 重播 | 6 | 0 | 80 | 80 | 40 | 13 | 80 |
| 重播 | 9 | 0 | 336 | 336 | 168 | 51 | 272 |
| 重播 | 12 | 0 | 1360 | 1360 | 680 | 210 | 1040 |

**最终状态**（[`:77-79`](zero_sum_periodic_destruction.md#L77-L79)）：不毁灭 10,296／784／286／0；不重播 0／8／3／32；重播 1,360／680／210／1,040。

对应 JSON：[json:690-693](simulations/zero_sum_periodic_destruction_results.json#L690)、[json:1396-1399](simulations/zero_sum_periodic_destruction_results.json#L1396)、[json:2102-2105](simulations/zero_sum_periodic_destruction_results.json#L2102)；15 条 verification 全 true [json:2142-2157](simulations/zero_sum_periodic_destruction_results.json#L2142-L2157)。

### 4.4 [`zero_sum_closure_time_selection.md`](zero_sum_closure_time_selection.md) —— "演化"被压成单一吸引子

**结论段原文**（[§5，`:85-90`](zero_sum_closure_time_selection.md#L85-L90)）：

> 但同一条规则也产生一个强吸引子：**最短的闭合模式几乎完全主导**
> （原文：*"the shortest closure mode dominates almost completely"*）。
>
> 要得到丰富的**分化生态**，下一个缺失的结构**不是**再加一个适应度参数；
> 它必须解释**为什么不同的闭合时间能够共存**——例如通过局域性、相容性，
> 或一个能让多种闭合类型都存活的**自催化网络**。

**数值结果**（[`:65-68`](zero_sum_closure_time_selection.md#L65-L68)）：

| 规则 | 末代平均闭合时间 | 末代多样性 | 主导模式 |
|:--|--:|--:|:--|
| `persist_only` | 17.11 | 8.01 | 无唯一胜者 |
| `closure_daughter` | 2.00 | 近于零 | `+-` |

公式 λ=ln2/T（[`:56`](zero_sum_closure_time_selection.md#L56)）；*"该结果是**种群份额的选择**，而不是字面意义上的灭绝"*（[`:71`](zero_sum_closure_time_selection.md#L71)）。

对应 JSON：[closure_time_selection json:22-29](simulations/zero_sum_closure_time_selection_results.json#L22)（persist_only）、[json:34-41](simulations/zero_sum_closure_time_selection_results.json#L34)（closure_daughter，`top_fraction = 0.99999999907`）。

> **补充（第 5 组，与"寿命"最直接）**：[`zero_sum_closure_exit.md`](zero_sum_closure_exit.md) 是唯一把"寿命"写进结论的。核心 [`:69`](zero_sum_closure_exit.md#L69)：*"这里活动层最终清空，是因为本程序另行设定了**有限寿命 $L=4$**，并且没有持续开放储层或重新播种。它**不表示**所有振动会同时闭合，也不表示一般模型必然清空。"* 尚未证明 5 条 [`:83-87`](zero_sum_closure_exit.md#L83-L87)，其中第 4 条正是 *"离散步序怎样接到连续时间"*。

---

## 总结：旧仿真里，"演化"是被实现了、还是被证明实现不了？

**答案分三层，必须分开说——混在一起会得出错误结论。**

### (a) "演化"作为**可运行的时间推进**：旧仿真里被**实现了**

- `simulations/` 10 个 `.py` 里 **8 个**有真正的时间推进循环（[closure_exit:139](simulations/zero_sum_closure_exit.py#L139)、[closure_time_selection:145](simulations/zero_sum_closure_time_selection.py#L145)、[cycle_evolution:270](simulations/zero_sum_cycle_evolution.py#L270)、[global_R_local_P:194](simulations/zero_sum_global_R_local_P.py#L194)、[living_universe:220](simulations/zero_sum_living_universe.py#L220)、[open_reservoir:115](simulations/zero_sum_open_reservoir.py#L115)、[periodic_destruction:193](simulations/zero_sum_periodic_destruction.py#L193)、[tri_layer_universe:198](simulations/zero_sum_tri_layer_universe.py#L198)），第 9 个有代循环（[reproduction_audit:206](simulations/zero_sum_reproduction_audit.py#L206)）。
- 它们产生逐代状态序列、分支增长、谱系选择、以及**唯一的灭绝时间**（`nilpotent_index = 6`，[reproduction_audit_results.json:52](simulations/zero_sum_reproduction_audit_results.json#L52)）。
- 证据：`zero_sum_periodic_destruction.md:41` 的逐代序列 $10,18,16,16,32,80,80,160,336$ 与 `continuous` 的 $10,18,32,60,110,210$（[L2_period.py:19-20](L2_period.py#L19-L20) 交叉引用）都是真实跑出来的整数序列。

### (b) "演化"作为**有寿命/有记忆的物理过程**：旧仿真里**既未实现，也未测量**

- 旧仿真**没有一处**测量关联时间、自关联函数或寿命分布（对 `simulations/` 检索 `autocorrel|correlation_time|relaxation|亚稳|老化|畴` **零命中**）。
- 唯一的"寿命"`LAYER_LIFETIME = 4` 是**硬编码输入**（[zero_sum_closure_exit.py:61](simulations/zero_sum_closure_exit.py#L61)），不是测量输出。旧仿真测的是**"是否熄灭/是否增长"**，不是**"能活多久"**。
- 旧仿真的"种群增长"大量依赖**数值重标定**（`rescale_steps = 38/137/142`，[cycle_evolution_results.json:51,66,81](simulations/zero_sum_cycle_evolution_results.json#L51)），不是物理后代数——这一点由 [zero_sum_reproduction_audit.md:111](zero_sum_reproduction_audit.md#L111) 自己承认：*"那是一个**隐藏的复制假设**（原文：*"That was a hidden copy assumption."*）"*。
- 旧仿真从来没有"稳定态"概念，因此也谈不上"离开稳定态要多久"。

### (c) "演化"作为**亚稳/长寿命/老化**：旧仿真里**未测**；**"实现不了"是在更新的一代 `zerogpu/` 里才被证明的**

- 这不是旧仿真的结论，而是 [zerogpu/EVAL_can_zero_evolve.md:18](zerogpu/EVAL_can_zero_evolve.md#L18) 的定理级上界：**τ₂ ≤ 1.18 微观步（L=4）⟹ 无亚稳畴、无长寿命、无老化、无记忆**。
- 机制（[`:110`](zerogpu/EVAL_can_zero_evolve.md#L110)）：零和约束 $|n_v|\le L$ 让每点成为**有限的反射随机游走**，有限状态 ⟹ 谱离散、间隙 $\Theta(1)$、**没有指数小间隙** ⟹ 没有亚稳性可言。
- 正面 no-go（[zerogpu/z0_entropy_barrier.py:22-27](zerogpu/z0_entropy_barrier.py#L22-L27)）：**单一原生约束无法产生指数寿命；沿细管的弛豫至多是多项式的。**
- 实验数据（[zerogpu/results/z0_life.json](zerogpu/results/z0_life.json)）：**畴寿命中位数在全部 7 个相点都是 1.0 代**；`acf_tau` 多为 1。
- 处方（[zerogpu/EVAL_can_zero_evolve.md:114-126](zerogpu/EVAL_can_zero_evolve.md#L114-L126)）：要有演化，必须**加有偏核（付价签）**，因为 *"忠实 ∧ 守恒 ⟹ 自由弛豫（无结构）"*（[`:100`](zerogpu/EVAL_can_zero_evolve.md#L100)）；否则 *"诚实的结论是：**Zero 在电脑里能"跑"，但按现有条款它不会"演化"**"*（[`:124`](zerogpu/EVAL_can_zero_evolve.md#L124)）。

### 一句话

> **旧仿真（`simulations/` + 根目录 `zero_sum_*`）把"演化"实现成了"确定性全分支枚举的时间推进"，并能给出熄灭/持续/灭绝时间；但它们从未测量过寿命、关联时间或自关联函数，"有寿命的演化"在旧仿真里既未实现也未测量。而"亚稳畴/长寿命/老化在忠于 Z0③ 的零和动力学下不可能"这一更强的否定结论，是在 `zerogpu/`（Oct 5–6 那一代）才被证明的——τ₂ ≤ 1.18 微观步，畴寿命中位恒为 1 代。**

**边界声明（只报告读到的）**：`zerogpu/` 中的 `EVAL_can_zero_evolve.md` §6 自列 4 条诚实边界，含 *"§2 的上界是**上界**；真值可能更小"*（[`:132`](zerogpu/EVAL_can_zero_evolve.md#L132)）与 *"V=256 的真谱**没有**直接对角化，结论靠变分上界 + 小图验证"*（[`:136`](zerogpu/EVAL_can_zero_evolve.md#L136)）。我**没有**独立复跑任何一个脚本，所有数字均为从 JSON/MD 中读取。
