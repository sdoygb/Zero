# Zero 形式化覆盖矩阵（机器生成，勿手改）

> 由 `ZeroLean/scripts/coverage_matrix.py` 从 kernel 实际状态 + 覆盖账生成，可复跑。
> 判据：三分归属（① 已证 ／ ② 具名输入 ／ ③ 不可形式化），**不存在第四种状态**。
> 依赖清单（每条定理用了哪些公理）见 `zero_lean` 的输出或 `ZeroLean/.lean-report.json`。

| 篇 | kernel 定理 | 覆盖账 F/A/N | 状态 |
|---|---|---|---|
| Z001 | — | 1/1/1 | ⏳ 待形式化（F 类未做） |
| Z002 | z002_energy_nonneg、z002_energy_eq_zero_iff、z002_mulVec_mem_zeroSum | 2/1/0 | ① 已证 |
| Z003 | z003_no_fin_dim_commutator、z003_commutator_ne_scalar | 1/0/2 | ① 已证 |
| Z004 | z004_born_one_forces_parallel | 3/0/0 | ① 已证（部分） |
| Z005 | — | 2/1/0 | ⏳ 待形式化（F 类未做） |
| Z006 | — | 2/1/0 | ⏳ 待形式化（F 类未做） |
| Z007 | z007_povm_probability_sum、z007_no_signal_trace | 2/1/0 | ① 已证 |
| Z008 | — | 3/0/0 | ⏳ 待形式化（F 类未做） |
| Z009 | — | 2/0/0 | ⏳ 待形式化（F 类未做） |
| Z010 | — | 0/2/1 | ② 仅具名输入 |
| Z011 | z011_sum_eq_card_diff、z011_card_true_eq_half、z011_balanced_count | 1/1/0 | ① 已证 |
| Z012 | — | 2/1/0 | ⏳ 待形式化（F 类未做） |
| Z013 | z013_graviton_unique_D | 2/1/0 | ① 已证（部分） |
| Z014 | — | 3/0/0 | ⏳ 待形式化（F 类未做） |
| Z015 | — | 0/0/1 | ③ 不可形式化 |
| Z016 | — | 0/0/2 | ③ 不可形式化 |
| Z017 | — | 3/0/0 | ⏳ 待形式化（F 类未做） |
| Z018 | — | 0/0/2 | ③ 不可形式化 |
| Z019 | — | 1/2/0 | ⏳ 待形式化（F 类未做） |
| Z020 | z020_iterate_eventually_periodic、z020_iterate_periodic_eventually、z020_monotone_on_cycle_constant | 2/1/0 | ① 已证 |
| Z021 | z021_relation_dirs、z021_relation_dirs_forced_three、z021_coupling_two_hypotheses …(+25) | 2/1/0 | ① 已证 |
| Z022 | z022_central_choose_mono、z022_saturation_small、z022_saturation_fails_range …(+3) | 2/1/0 | ① 已证 |
| Z023 | — | 0/1/2 | ② 仅具名输入 |
| Z024 | z024_charge_conservation | 2/0/1 | ① 已证（部分） |
| Z025 | — | 2/1/0 | ⏳ 待形式化（F 类未做） |
| Z026 | — | 2/0/1 | ⏳ 待形式化（F 类未做） |
| Z027 | — | 2/1/0 | ⏳ 待形式化（F 类未做） |
| Z028 | — | 2/0/1 | ⏳ 待形式化（F 类未做） |
| Z029 | — | 2/0/1 | ⏳ 待形式化（F 类未做） |
| Z030 | — | 0/0/2 | ③ 不可形式化 |
| Z031 | — | 2/0/1 | ⏳ 待形式化（F 类未做） |
| Z032 | — | 3/0/0 | ⏳ 待形式化（F 类未做） |
| Z033 | z033_mulVec_mem_zeroSum、z033_sum_zero_of_solvable、z033_diff_constant …(+1) | 1/0/2 | ① 已证 |
| Z034 | — | 2/1/0 | ⏳ 待形式化（F 类未做） |
| Z035 | z035_chsh_classical_bound | 2/1/0 | ① 已证（部分） |
| Z036 | — | 2/0/1 | ⏳ 待形式化（F 类未做） |
| Z037 | — | 0/0/1 | ③ 不可形式化 |
| Z038 | — | 0/0/1 | ③ 不可形式化 |
| Z039 | — | 0/0/1 | ③ 不可形式化 |
| Z040 | — | 2/0/0 | ⏳ 待形式化（F 类未做） |
| Z041 | z041_zero_sum_selection | 1/1/1 | ① 已证 |
| Z042 | z042_step_le_one、z042_neg_one_le_step、z042_step_eq_one_iff …(+20) | 2/0/1 | ① 已证 |
| Z043 | — | 0/0/2 | ③ 不可形式化 |
| Z044 | z044_mod_energy_general、z044_K1、z044_K2 …(+1) | 2/0/0 | ① 已证 |
| Z045 | z045_flow_trivial_of_eq、z045_flow_offdiag、z045_nontrivial | 1/1/0 | ① 已证 |
| Z046 | — | 2/1/0 | ⏳ 待形式化（F 类未做） |
| Z047 | z047_wordStep_neg、z047_sum_neg、z047_prefix_neg …(+3) | 1/1/0 | ① 已证 |
| Z048 | z048_prefix_succ、z048_prefix_sub_even、z048_family_prefix …(+6) | 3/0/0 | ① 已证 |
| Z049 | z049_card_image_pair、z049_wordNeg_ne_self、z049_wordNeg_mem_shell …(+3) | 1/1/0 | ① 已证 |
| Z050 | z050_max_area_upper_bound、z050_max_area_at_top | 1/1/0 | ① 已证 |
| Z051 | z051_partition_normalized、z051_partition_eq_one、z051_entropy_mean | 2/1/0 | ① 已证 |
| Z052 | z052_even_b_gcd | 3/0/0 | ① 已证（部分） |
| Z053 | — | 0/0/1 | ③ 不可形式化 |
| Z054 | — | 1/1/1 | ⏳ 待形式化（F 类未做） |
| Z055 | — | 0/0/1 | ③ 不可形式化 |
| Z056 | — | 2/1/0 | ⏳ 待形式化（F 类未做） |
| Z057 | z057_singlet_integral_charge | 3/0/0 | ① 已证（部分） |
| Z058 | z058_nine_not_dividing_168、z058_168_factor | 2/0/1 | ① 已证 |
| Z059 | z059_skew_hermitian_structure、z059_traceless_diagonal_sum | 1/1/1 | ① 已证 |
| Z060 | — | 1/1/1 | ⏳ 待形式化（F 类未做） |
| Z061 | z061_girth_le_six、z061_six_le_girth、z061_girth_eq_six | 2/0/1 | ① 已证 |
| Z062 | — | 1/1/1 | ⏳ 待形式化（F 类未做） |
| Z063 | z063_center_S3_trivial | 1/0/2 | ① 已证 |
| Z064 | z064_quadform、z064_div_zero_adj_eq、z064_connected_constant …(+1) | 1/0/1 | ① 已证 |
| Z065 | z065_tree_card_edgeFinset、z065_tree_cycle_rank_zero、z065_heawood_regular …(+7) | 1/1/1 | ① 已证 |
| Z066 | z066_card_R1、z066_card_R2、z066_naive_reading_gives_85 …(+4) | 1/0/1 | ① 已证 |
| Z067 | — | 1/0/2 | ⏳ 待形式化（F 类未做） |
| Z068 | — | 2/0/1 | ⏳ 待形式化（F 类未做） |
| Z069 | z069_const_word_sum_ne_zero、z069_has_both_signs、z069_nat_constant_of_no_change …(+4) | 2/0/1 | ① 已证 |
| Z070 | — | 1/0/1 | ⏳ 待形式化（F 类未做） |
| Z071 | — | 1/1/1 | ⏳ 待形式化（F 类未做） |
| Z072 | z072_mem_neighbors、z072_injOn、z072_neighbors_card …(+4) | 1/1/0 | ① 已证 |
| Z073 | — | 2/0/1 | ⏳ 待形式化（F 类未做） |
| Z074 | z074_mod_lt、z074_bool_eq_false_of_ne_true、z074_bool_eq_true_of_ne_false …(+19) | 1/1/0 | ① 已证 |
| Z075 | z075_runs_even_smallL、z075_runs_odd_smallL、z075_runs_L12_even_c3 …(+19) | 2/1/0 | ① 已证 |
| Z076 | z076_kappa_less_than_target | 1/1/0 | ① 已证 |
| Z077 | z077_layer_partition_q | 2/0/0 | ① 已证（部分） |
| Z078 | — | 2/0/0 | ⏳ 待形式化（F 类未做） |
| Z079 | z079_W4_card、z079_layer_table、z079_layer_counts …(+11) | 3/0/0 | ① 已证 |
| Z080 | z080_Mplus_distribution、z080_ballot_formula_L4、z080_mass_sum_L4 …(+6) | 3/0/0 | ① 已证 |
| Z081 | — | 0/1/1 | ② 仅具名输入 |
| Z082 | z082_koide_identity、z082_koide_Q_two_thirds、z082_equal_weights_kill_first_harmonic …(+1) | 1/0/0 | ① 已证 |
| Z083 | z083_telescoping、z083_closed_loop_zero | 1/1/0 | ① 已证 |
| Z084 | — | 0/1/1 | ② 仅具名输入 |
| Z085 | — | 1/1/1 | ⏳ 待形式化（F 类未做） |
| Z086 | z086_mean_deviation_sum_zero、z086_hyperplane_finrank | 1/1/1 | ① 已证 |
| Z087 | — | 2/1/0 | ⏳ 待形式化（F 类未做） |
| Z088 | — | 2/0/1 | ⏳ 待形式化（F 类未做） |
| Z089 | z089_four_vertex_three_regular_is_complete、z089_K4_counts | 2/1/0 | ① 已证 |
| Z090 | z090_omega_cube、z090_omega_ne_one、z090_omega_sum …(+8) | 2/1/1 | ① 已证 |
| Z091 | — | 0/1/0 | ② 仅具名输入 |
| Z092 | — | 0/1/0 | ② 仅具名输入 |
| Z093 | z093_s3_three_conj_classes、z093_s3_order | 1/1/0 | ① 已证 |
| Z094 | — | 0/0/1 | ③ 不可形式化 |
| Z095 | z095_row_sum、z095_two_mul_edges、z095_eigen_product …(+2) | 1/2/0 | ① 已证 |
| Z096 | z096_omega_cubed、z096_omega_eq、z096_omega_re …(+6) | 2/0/0 | ① 已证 |
| Z097 | z097_gmn、z097_right_left_hypercharge | 1/1/0 | ① 已证 |
| Z098 | — | 1/1/0 | ⏳ 待形式化（F 类未做） |
| Z099 | — | 0/1/0 | ② 仅具名输入 |
| Z100 | — | 0/1/0 | ② 仅具名输入 |
| Z101 | — | 0/1/0 | ② 仅具名输入 |
| Z102 | z102_laplacian_on_zero_sum、z102_phi_three | 1/1/0 | ① 已证 |
| Z103 | — | 1/1/0 | ⏳ 待形式化（F 类未做） |
| Z104 | z104_grid_dispersion | 2/0/0 | ① 已证（部分） |
| Z105 | z105_deriv_wave、z105_deriv2_wave、z105_wave_solution …(+1) | 2/0/0 | ① 已证 |
| Z106 | z106_proj_traceless、z106_proj_idempotent、z106_proj_fixes_traceless …(+1) | 3/0/0 | ① 已证 |
| Z107 | z107_symmetric_entries、z107_traceless_diagonal | 2/1/0 | ① 已证 |
| Z108 | — | 0/1/1 | ② 仅具名输入 |
| Z109 | z109_deriv_kinetic、z109_deriv_potential、z109_energy_deriv …(+2) | 2/0/0 | ① 已证 |
| Z110 | — | 1/1/0 | ⏳ 待形式化（F 类未做） |
| Z111 | — | 0/0/2 | ③ 不可形式化 |
| Z112 | — | 0/2/1 | ② 仅具名输入 |
| Z113 | — | 2/1/0 | ⏳ 待形式化（F 类未做） |
| Z114 | z114_neel_eigen、z114_neel_image、z114_gap_alternates …(+3) | 2/1/0 | ① 已证 |
| Z115 | — | 1/2/0 | ⏳ 待形式化（F 类未做） |
| Z116 | z116_exp_conjTranspose_mul_self、z116_normSq_eq_re_dotProduct、z116_normSq_mulVec_of_conjTranspose_mul_self …(+4) | 1/0/1 | ① 已证 |
| Z117 | z117_aliasing、z117_aliasing_funext、z117_beyond_nyquist …(+6) | 2/0/1 | ① 已证 |
| Z118 | z118_dispersion_identity、z118_inverse_of_dispersion、z118_dispersion_of_inverse | 2/1/0 | ① 已证 |
| Z119 | — | 2/1/0 | ⏳ 待形式化（F 类未做） |
| Z120 | — | 3/0/0 | ⏳ 待形式化（F 类未做） |
| Z121 | z121_zero_eq_cos、z121_cos_three_pi_four、z121_cos_five_pi_four | 1/2/1 | ① 已证 |
| Z122 | z122_two_zeros_midpoint、z122_two_zeros_separation | 1/1/0 | ① 已证 |
| Z123 | z123_delta0_pos | 1/1/1 | ① 已证 |
| Z124 | z124_single_color_not_closed、z124_three_colors_closed、z124_singlet_amplitude_zero …(+1) | 3/1/0 | ① 已证 |
| Z125 | — | 0/1/2 | ② 仅具名输入 |
| Z129 | — | 2/1/1 | ⏳ 待形式化（F 类未做） |
| Z130 | — | 0/0/3 | ③ 不可形式化 |
| Z131 | z131_flow_fixes_diagonal、z131_offdiagonal_ratio、z131_frequency_is_log_ratio …(+3) | 4/0/0 | ① 已证 |
| Z132 | — | 2/0/1 | ⏳ 待形式化（F 类未做） |
| Z133 | z133_omega_cubed、z133_omega_ne_one、z133_cos_two_pi_div_three …(+2) | 3/1/0 | ① 已证 |
| Z134 | — | 1/0/3 | ⏳ 待形式化（F 类未做） |
| Z135 | — | 0/0/0 | ❔ 未记账 |

## 汇总

- 库内篇数：**132**
- 有 kernel 定理的篇：**63**
- 手写定理总数：**329**（口径：`.lean` 里的 `theorem/lemma zNNN_*`）
  - ① 已证：55 篇
  - ⏳ 待形式化（F 类未做）：43 篇
  - ③ 不可形式化：13 篇
  - ② 仅具名输入：12 篇
  - ① 已证（部分）：8 篇
  - ❔ 未记账：1 篇

## 缺口（F 类有内容但 kernel 里没有）—— 43 篇

Z001、Z005、Z006、Z008、Z009、Z012、Z014、Z017、Z019、Z025、Z026、Z027、Z028、Z029、Z031、Z032、Z034、Z036、Z040、Z046、Z054、Z056、Z060、Z062、Z067、Z068、Z070、Z071、Z073、Z078、Z085、Z087、Z088、Z098、Z103、Z110、Z113、Z115、Z119、Z120、Z129、Z132、Z134

## 参考
- 覆盖账：`_staging/FORMALIZATION_COVERAGE.md`
- 待办账：`_staging/LEAN_BACKLOG.md`（70 行）
- 审计报告：`_staging/FORMALIZATION_AUDIT.md`
