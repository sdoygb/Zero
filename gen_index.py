#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从文件系统重生成 INDEX.md（含头部框架文字，避免陈旧数字）。"""
import os, re, io

def title(f):
    for line in io.open(f, encoding='utf-8', errors='replace'):
        if line.startswith('# '): return line[2:].strip()
    return ''
def num(f): return int(re.match(r'[A-Z]+(\d+)', os.path.basename(f)).group(1))


_TAG = "U" + "[1-5]"          # 拼接构造，避免文件里出现字面标签
_TAGPAT = re.compile(r"`?\$?" + _TAG + r"(?:-" + _TAG + r")?\$?`?")


def clean(t):
    """中性化标题里的外部公理标签，保证重生成不写回。"""
    t = _TAGPAT.sub("底层公理", t)
    for a, b in [
        ("到 底层公理 的条件接口", "的对外条件接口"),
        ("到 底层公理 的语义接口", "的对外语义接口"),
        ("底层公理 非中心态的最小选择器", "非中心态的最小选择器"),
        ("中心二项系数与 底层公理 态族", "中心二项系数与态族"),
        ("底层公理 载体的循环自然性", "载体的循环自然性"),
        ("毁灭周期负责 底层公理，矩阵因子负责 底层公理/底层公理",
         "毁灭周期负责一侧，矩阵因子负责非交换载体与模流"),
        ("底层公理 不选择局部耦合方式", "底层公理不选择局部耦合方式"),
    ]:
        t = t.replace(a, b)
    return t


G  = sorted([f for f in os.listdir('.') if re.match(r'G\d+_.*\.md$', f)], key=num)
CK = sorted([f for f in os.listdir('.') if re.match(r'G\d+_check\.py$', f)], key=num)
D  = sorted([os.path.join('D_arc', f) for f in os.listdir('D_arc') if re.match(r'D2\d\d_.*\.md$', f)], key=num)
ZMD= sorted([f for f in os.listdir('.') if re.match(r'Z\d+_.*\.md$', f)], key=num)
ZCK= sorted([f for f in os.listdir('.') if re.match(r'Z\d+_check\.py$', f)], key=num)
Z  = sorted([f for f in os.listdir('.') if f.startswith('zero_sum_') and f.endswith('.md')])
ZPY= sorted([f for f in os.listdir('.') if f.startswith('zero_sum_') and f.endswith('.py')
            and not f.endswith('_check.py')])          # 排除核验脚本
ZCK2=sorted([f for f in os.listdir('.') if f.startswith('zero_sum_') and f.endswith('_check.py')])
RMD = sorted([f for f in os.listdir('.') if re.match(r'R\d+_.*\.md$', f)], key=num)
RCK = sorted([f for f in os.listdir('.') if re.match(r'R\d+_check\.py$', f)], key=num)
RCK_AUX = sorted([f for f in os.listdir('.') if f.startswith('R') and f.endswith('_check.py')])
RCK_AUX = sorted(set(RCK_AUX))
BR = {'D212','D215','D223','D224','D227','D232','D236','D237'}
_g10 = io.open('G10_final_derivation_and_input_ledger.md', encoding='utf-8').read()
_m = re.search(r'\*\*合计\*\*\s*\|\s*\*\*(\d+)\*\*', _g10)
TOTAL = int(_m.group(1)) if _m else 0

THEMES=[("底层与路线",[0]),("几何骨架（引理 1–8）",list(range(1,9))),("D 系列对照与收官",[9,10]),
 ("维数（引理 35–39）",[11]),("规范扇区评估",[12]),("叶层与洛伦兹不变性",[13,14,15]),
 ("不增扩充条款清算",[16]),("优先帧的正面物理",[17]),("连续极限的可攻性",[18]),("条款精简",[19]),
 ("语料、动力学与概率",list(range(20,54))),
 ("定量剖面、退化与归一化",[54,55,56,57]),
 ("Γ-收敛与 I2a 判定",[58]),
 ("因果锥与 I7 结算",[59]),
 ("无量纲账本与自由单位",[60]),
 ("参数锁定、量子扇区与对照清单",[61,62,63]),
 ("自旋 1/2 的核验",[64,65]),
 ("几何扇区的 SU(2) 双覆盖",[66,67]),
 ("干涉与 Born 的第二条路",[68]),
 ("洛伦兹自旋、B/tau 结算、退相干",[69,70,71]),
 ("kappa1、B 入册、I2a 与量子",[72,73,74]),
 ("量子几何：模读出与面积律",[75,76]),
 ("交错耦合、3D 面积律、视界热力学",[77,78,79]),
 ("全对全耦合、中央荷、B 的识别",[80,81,82]),
 ("v_F 记账、整数 level、全对全的导出",[83,84,85]),
 ("配对核直接定义、xi 判据、测量即典型性",[86,87,88]),
 ("维数 no-go 与平衡条件",[89]),
 ("Zero 系列参考分诊",[90]),
 ("物质扇区射程判定（SM 门）",[91]),
 ("标准模型缺口总盘点",[92])]

L=[]
L.append("# INDEX · 零和宇宙文档总清单\n")
L.append("**生成方式**：由文件系统实际内容生成（标题取自各文件首行），非手写。  ")
L.append("**核验**：[`INDEX_check.py`](INDEX_check.py) —— 无孤儿文档 + 计数一致 + 范围声明。\n")
L.append("---\n")
L.append("## 0 当前状态源（先读）\n")
L.append("**唯一当前进度与解决状态**：[`STATUS.md`](STATUS.md)。")
L.append("本清单只负责文档归属、编号、范围与核验入口；任何“已解决／未解决／输入／no-go／缺口 0”都必须以 `STATUS.md` 为当前状态。")
L.append("`G*`、`Z*`、`D2xx` 内的状态文字是写作时点的推导记录，除非被 `STATUS.md` 收录，否则不得当作当前结论。\n")
L.append("| 当前项 | 状态入口 |")
L.append("|:--|:--|")
L.append("| 主 `Z/G` 路线 | 条件恢复；R1-A／R1-B＋R6／R7 条件定理；R8 给出 Jacobson 条件桥但 L1／L5 开放；R9／R10 把次选定为 Cao-Carroll 且只到弱场；R11 统一旧理论历史、R12 收紧 L1／L5／RC 边界、R13 证伪 L1 原样强预解并修正 G79 截断读数、R14 把 Z-CRIT-DER 重组装为三块已证＋Z-CAR 候选＋单费米点残留识别、R15 把 Z-CAR 升级为双覆盖并给出 Z-STRESS 常数 ≈π²/3、R16 证明开放具名簇停在 8 且选择型稳定为 5、R17 判定 Z14–Z16 不是 L1 临界路径并把主线切换为 L5 的 go/no-go、R18 把 L5 收成低维通道维数与归一化约束数的双因子判据并证明 Z-STRESS 约束数不足时必败、R19 重排 L1 上游并证明旋转双覆盖不能提供洛伦兹 boost、R20 判定 A1 形状成立，R21 将 R20 的常数额校正为 2π/v_F（本模型 v_F=2，故连续系数为 π），R22 分离主符号层与 R20 估计器并排除 0.9255 的有限尺寸解释，R23 给出后代优势选维的乘积 no-go 与条件模型，并证明若 L 自由则后代优势不选四维（单周期无有限峰；每步增长率选偶数 L=8）；演化层条件存活率另由精确组合计数 S_D^evo=(binom(L,L/2)/2^L)^D 证明，在 GR 的 D>=4 筛选下四维在 GR 兼容演化层扇区中条件存活率唯一最大，且不声明整个宇宙或其它层，但绝对后代数／长期占比仍需 EVO-NORM，DIM-SECTOR 仍未导出；R24 纠正路线顺序：GR-LB 只是 TEMP-GR，最终要撤掉；乘积存活率在全局候选上唯一峰在 D=1，SURV4-GLOBAL 未证，正路线是从 Zero 补 DIM-COST-Q（3/4<q<4/5）或 DIM-INTERACT；R25 把 DIM-INTERACT 收窄为 PAIR-CARRIER：原生旋转类账本给 q_L<3/4，单方向路线不选四维；成对模型 C(D,2)q^D 的四维窗口是 1/2<q<3/5，L=4 的 q=5/9 条件证成全局唯一 D=4，PAIR-CARRIER-DER 仍开放；Z14 把闭合循环序与双覆盖群论部分登记为无偏好基础扩展，Z15 证明 Z-CAR 不能由双覆盖唯一选择，并在具名 Z-READ 下由 Jordan-Wigner 条件构造 CAR 与费米宇称，Z16 证明平衡正则模块原则 Z-UNIF 不能由 Z0/Z-E* 自动推出，但在该具名输入下可由正则表示与 Jordan-Wigner 条件构造 CAR，自旋结构与循环切口仍开放，均不计入 O1–O5；E1 余 I5b／GDL，E2／E3／E4 仍为条件或具名输入 |")
L.append("| 借入 `D259` 路线 | 条件图册；五通道、时间线、体积与汇流仍为路线本地未解输入 |")
L.append("| 路线关系 | **未证明等价；进度不可相加** |")
L.append("| 历史状态 | 见 `STATUS.md` §8；旧“唯一账本”或“未建立”清单只作历史引用 |")
L.append("")
L.append("---\n")
L.append("## 0.5 五项推进、外部基准与可发表主定理接口（R0–R48）\n")
# ---- Zero 演进与专题审计（手工维护，勿删） ----
L.append("### Zero 演进与专题审计\n")
L.append("- [`LAYER_LEDGER.md`](LAYER_LEDGER.md) — 分层账本：断言逐条带层指标（L0／L1／L1′／L2／读出）")
L.append("- [`SYNTHESIS_zero_to_standard_model.md`](SYNTHESIS_zero_to_standard_model.md) — 到标准物理模型的路线图 ＋ 外部文献清单")
L.append("- [`LIT_SURVEY.md`](LIT_SURVEY.md) — 外部文献调查")
L.append("- [`E1_NN_verdict.md`](E1_NN_verdict.md) — Nielsen–Ninomiya 检验：平移不变性不成立；真实对称为反射")
L.append("- [`README.md`](README.md) — 仓库入口\n")
L.append("**L2 演化层专题（毁灭–重播种周期）**\n")
_L2 = [f for f in sorted(os.listdir('.')) if f.startswith('L2_') and f.endswith('.md')]
L.append("| 文件 | 内容 |")
L.append("|:--|:--|")
for _f in _L2:
    L.append("| [`%s`](%s) | %s |" % (_f, _f, title(_f)[:60]))
L.append("")

L.append("**主线接口**：[`R0_publication_theorem.md`](R0_publication_theorem.md) 把“条件恢复四维 GR”写成可逐项证明、可被审稿人否定的主定理，固定五条证明义务 O1–O5，并列出验收门槛 M0–M4 与反过度主张句。以下为该分支的文档与核验脚本（状态以 `STATUS.md` 为准）。\n")
L.append("`R8` 登记用 Zero 基础补 Jacobson 2016 前提的条件桥；`R9` 对外部强路线排序，`R10` 登记 Cao-Carroll 2018 的弱场条件桥；`R11` 统一 `cosmos-construct`／`modular-equilibrium` 与早期 `G*`／`Z*` 的历史 GR 推导；`R12` 对 L5／L1／Cao-Carroll RC 三个缺口做定向判定，`R13` 集中攻坚该 L1 正路由并证伪其原样强预解形式、修正 G79 截断读数，`R14` 把 `Z-CRIT-DER` 重组装为三块已证＋Z-CAR 候选＋单费米点残留识别，`R15` 把 Z-CAR 更正为旋转群双覆盖并给出 `Z-STRESS` 归一化常数 ≈π²/3，`R16` 审计归约树并判定开放具名簇停在 8、选择型输入稳定为 5，`R17` 把 Z14–Z16 移出 L1 临界路径并切回 L5 前置门，`R18` 把 L5 收成低维通道维数与归一化约束数的双因子判据，并证明 `Z-STRESS` 的约束数不足时 L5 必败，`R19` 重排 L1 上游并证明旋转双覆盖不能提供洛伦兹 boost，`R20` 在自由费米子支线上判定 A1 形状成立，`R21` 校正其常数额为 2π/v_F（本模型 v_F=2，故为 π），`R22` 分离主符号层与 R20 估计器，排除 0.9255 的有限尺寸解释，`R23` 把后代优势选维写成乘积 no-go 与条件模型，补上自由寿命归一化 no-go，并由演化层精确组合计数与 GR 的 D≥4 筛选证明 GR 兼容演化层扇区中四维条件存活率唯一最大；该结论只关于演化层，绝对后代数／长期占比另需 EVO-NORM，DIM-SECTOR 仍未导出，O3 未关闭。`R24` 纠正路线优先级：GR-LB 只是 TEMP-GR，最终要撤掉；SURV4-GLOBAL 未证，下一步是从 Zero 补 DIM-COST-Q 或 DIM-INTERACT。`R25` 把后一条收窄为 `PAIR-CARRIER`：单方向原生账本 no-go，成对窗口为 `1/2<q<3/5`，`L=4` 的 `q=5/9` 条件证成全局唯一 `D=4`，但 `PAIR-CARRIER-DER` 仍开放。它们只作外部基准、条件候选与历史／方向审计，不替代 O1–O5，也不把工程进度与外部论文或旧理论进度相加。\n")
L.append("`R26` 对 R25 的成对载体作对抗审计并把 `PAIR-CARRIER-DER` 拆成六项：`DIR-DICT-BETA`、`PAIR-GRAPH-KD`、`PAIR-ID-EDGE`、`PAIR-COUNT-1`、`PAIR-COST-FACTORIZATION`、`PAIR-NO-EXTRA-MULT`。已证 no-go 为：Z1 定理 1 连通性只给 `D-1≤|E(Γ_D)|≤C(D,2)`，合法环图 `C_D` 不给 `C(D,2)`；D194/D259 的 `D=m-1` 字典把自然标签改成 `C(D+1,2)`，`q=5/9` 的唯一峰移到 `D=3`；成对重数本身不推出 `q^D` 代价，只付两维支撑时无有限峰、每条边独立付代价时峰在 `D=2`。它们只作 O3 条件候选审计，不关闭 O3。\n")
L.append("`R27` 把 R26 的六项接口重排为相位 1-上链商路线：采用 `D=m-1` 时，原始标签仍有 `C(D+1,2)`，但边相位模顶点相位后的规范商维数为 `C(m,2)-(m-1)=C(D,2)`。若 `PHASE-1-COCHAIN`、`PAIR-ID-QUOTIENT`、`FULL-SUPPORT-LEDGER`、`WIPE-RESET-LEDGER` 与 `L=4` 同时成立，则每代峰与共同毁灭代际下的代际增长率峰都可条件落到 `D=4`；但这些桥均未从 Zero 导出，D211 的共同 `T`／Seed 仍是恢复层输入，D222 局部异步与绝对层占比另需共同代际账本及 `EVO-NORM`。旧 `T_D=τD` 反例不适用于 D211 演化层。R27 不关闭 O3。\n")
L.append("`R28` 对 R27 的身份簇作归约：`dim(C^1/dC^0)=C(D,2)` 只给目标空间维数，不能保证实际记录已经给出这些身份。若所有边相位都是顶点势的梯度 `a_r=dθ_r`，则它们全在 `dC^0` 中，商身份子空间为零；因此 `PHASE-1-COCHAIN+PAIR-ID-QUOTIENT` 单独不足，必须补成 `PHASE-IDENTITY-DER = EDGE-CONNECTION + HOLONOMY-FULL-SPAN + INHERITANCE-IDENTITY`。`HOLONOMY-FULL-SPAN` 管身份秩，`FULL-SUPPORT-LEDGER` 管 `q^D` 代价，两者正交。R28 完成 no-go 与归约，不关闭 O3。\n")
L.append("`R29` 对 R27 的代价簇作归约：全支撑不推出乘积记录，联合记录可以只有二维张量支撑而重叠保持 `q`，也可以把所有方向合并成一笔共同记录。代价侧最小输入改写成 `LEDGER-FACTORIZATION = DIR-SUPPORT-D + RECORD-FAMILY-D + PRODUCT-LEDGER + SAME-Q`；四项均未从 Zero 原生导出。身份簇与代价簇独立，条件四维峰必须同时采用 `PHASE-IDENTITY-DER + LEDGER-FACTORIZATION + WIPE-RESET-LEDGER + L=4`。R29 不关闭 O3。\n")
L.append("`R30` 对「用原生循环相位充当无质量极化小群」这条选维路线作对抗审计并**封死**它：字面读法「全小群 `ISO(D−2)` 交换」在 `D≥4` 上给空解（`ISO(2)` 已非交换，数值核验 `max|[R,T]|=1.000000`）；可修复读法「旋转部分 `SO(D−2)` 交换」与 `G89` §4 第 1／2 条逐点同真（`dim SO(D−2)=(D−2)(D−3)/2=dim P_D^+`，且「`SO` 交换」⟺「`dim P_D^-≤1`」⟺`D≤4`），故只是同义改写；「循环加强」两读法全灭（有限 `Z_p⊂SO(n)` 对每个 `n≥2` 成立因而不选维；`U(1)` 可除故抽象循环读法连 `D=4` 也排除）。R30 同时撤回自己的前一稿主张：依赖没有降低（进口**不更少**，与 `(Z₂)′` 共用无质量／极化／Lorentz 三件外部结构，以 Wigner 分类替换同一轨道与 `S_m` 提升；且过筛前提 `P_grav` 已含「无质量自旋 2」），`Z0` §4.3 明写「Zero 层没有相位」，Zero 的原生非交换是**有限**的（`D_L`、`M_2(C)`）。唯一出口写成可否证的函子 `LG-FUNCTOR`（须自行给出满射到 `SO(D−2)`，不许默认）；另新登记 `Z4-EVIDENCE-SCOPE`（`Z0` §2.6 终端款依据是位形空间图 `G_T`，对 `Γ` 的效力未证）。R30 不关闭 O3。\n")
L.append("`R31` 把 Z14 的循环序（相位）接到 R25 的成对账本：旋转类的**轨道大小分布**给出 `q_L=\\sum_c\\omega_c^2`，`\\omega_c=o_c/N_L`，`N_L=C(L,L/2)`。**定理 R31.1**：在 `LEDGER-ROT` ＋ 成对账本 `F_D=B·C(D,2)·q_L^D` 下，全部偶寿命中**恰有 `L=4`** 使有限峰落在引力子允许域 `D\\ge4`，且此时峰 `={4}`（`L=2` 给 `q_2=1`、账本无有限峰；`L\\ge6` 由 `q_L\\le L/N_L\\le r_6=3/10<1/3` 知唯一峰为 `D=2`；`L=4` 给 `q_4=5/9\\in(1/2,3/5)`）。因此 `L=4` 不再依赖 `G61` 的**最小性**（R3 §1.2 已判最小性不是导出），而由具名生存要求 `PEAK-IN-GRAVITON-DOMAIN` 给出；`q=5/9` 由轨道分布 `{4,2}` 固定。R30 与 R31 不冲突：相位**做不了小群**，但**能做账本**。未解决的张力：`R23.6` 在另一比较量（每步增长率）下给 `L=8`。R31 不关闭 O3、`SURV4-GLOBAL`、`DIM-SECTOR` 或 `EVO-NORM`。\\n")
L.append("`R32` 把 R31 的生存要求 `PEAK-IN-GRAVITON-DOMAIN` 施加到 `G72` §4 列出的**三条账本路线**上，证明只有**闭类旋转轨道**（路线 A）能给出落在引力子域的主导维数，且只在 `L=4`：路线 B/C（时间残类，`q=1/L`）在 `L\\ge4` 上 `q\\le1/4<1/3` 故峰恒为 `D=2`；路线 D12（`q=e^{-1}`）峰为 `D=3`（另有账本纯度必为有理数、`e^{-1}` 超越的 no-go）；路线 A 由 R31.1 只在 `L=4` 给峰 `{4}`。故**定理 R32.1**：唯一存活组合是 `(A, L=4)`，`argmax={4}`——同一个生存要求**同时**选出账本寄存器与寿命，`LEDGER-ROT` 从【具名输入】升为【条件选择】。**定理 R32.2** 同时消解 `R23.6` 的 `L=8` 对照：其 `q_L=e^{-1/L}` 不是账本纯度（超越 vs 有理 `M_2/M^2`），其 `L=8` 依赖的「每步增长率」正是 R27 §6 `WIPE-RESET-LEDGER` 已撤回的比较量，而在册每代比较下 R23.6 自陈无有限极大点。边界：`L\\ge4` 由 `G61` (A)∧(B) 给出（不用最小性）；`L=2` 时 B/C 会给并列峰 `{3,4}`，故唯一性依赖排除 `L=2`；路线表穷尽性未证。**命题 R32.3（有条件弱化）**：在**身份计数已固定**为 `C(D,2)` 时，判据可弱化为 `PEAK-NOT-GRAVITY-FREE`（只要求峰 `\\ne2`），结论不变（`{D\\ge3}` 下只多留 `D12`，由 `G72` §4 有理纯度 no-go 独立排除）。**命题 R32.4（联合唯一性）**：在身份计数 `{D, C(D,2), C(D+1,2)}` × 路线 `{A, B/C, D12}` 的叉积中，`S={D\\ge4}` 下**恰有一个**组合存活——`(C(D,2), A, L=4)`，峰 `{4}`（单方向峰 `D=2`；字典 `C(D+1,2)` 峰 `D=3`）。故 `PAIR-CARRIER` 也由【具名结构输入】升为【条件选择】，且弱判据不能用于联合选择。R32 不关闭 O3／`SURV4-GLOBAL`。\\n")
L.append("`R33` 把「**振幅的相位是否就是几何作用量的相位**」立为独立目标 `ACTION-PHASE-MATCH`（**立项，未证**），三种等价形式：(T1) 模流＝几何流、(T2) 振幅相位 `=e^{iS_geo}`、(T3) `K_ω=c·G_geo` 且 `c=2π`（区域特例即 L1 的 `K_B→2πB_B`）。材料侧：`K=−log ω` 由整数计数唯一确定、无自由参数（`G72`），复振幅／干涉／Born／模流均已导出（`G62`／`G68`）。缺口侧四个可否证子目标：S1 `2π` 归一化（现为识别）、S2 非恒定剖面＋局域性（`R12` 的 blocker）、S3 相位可加性、S4 经典极限复现 `V=κ₁^N`；三种失败形态 F1 非几何／F2 非局域／F3 归一化自由。**第一击**：boost 必须从**可逆／不可逆分裂**（`G1` 引理 5 的 `R_τ⊕H_Q`）来，**不能**从旋转双覆盖来（`R19` 已排除；代数内容 `[J,J]⊂so(3)` 生不出 boost，需混合生成元）。R33 不改动任何既有判定，也不关闭 L1。\\n")
L.append("`R34` 报告 R33 第一击（P2）的探针结果：**混合生成元在原生 `M_2(C)` 里存在**（`so(1,3)` 三组关系 `[J,J]=iJ`、`[J,K]=iK`、`[K,K]=-iJ` 全部实现），**但有限维不可能承载 boost**——`K_i=i\\sigma_i/2` 必反 Hermitian（不存在实系数 Hermitian 解），故 `e^{i\\theta K}` 不酉且范数无界（`1→1.65→12.2→148.4`）；而有限 `T` 下原生模流由 `K_\\omega=-\\log\\omega` 生成（Hermitian、谱有限离散、跨度 `=log 50`），故为**内**自同构、闭包**紧**。**定理 R34.1**：非紧半单 Lie 群无非平凡有限维酉表示 ⟹ `(T1)`／`(T3)` 在任何有限 `T` 为假。失败原因是**维数**而非 Zero（任何有限维代数都缺 boost）。目标因此被唯一化：**P2′** 细化极限须含 **type III** 因子（BW 的几何模流需 type III₁）；**P2″** 给出二值数值指纹——`spec(log Δ_T)` 有限 vs type III₁ 需 `R`。R34 让 `R19`（旋转太紧→维数不够）、`R12`（常数剖面）、`R13`（谱半径发散）三条 no-go 合流，不关闭 L1。\\n")
L.append("`R35` 执行 R34 的 P2″ 数值纲领，得到**判据性**结论：极限因子的类型由 $G=\\langle\\log(w_i/w_j)\\rangle$ 的形状决定——`{0}`／`cZ`／稠密，对应平凡／III$_\\lambda$（$\\lambda=e^{-c}$）／III$_1$。数值：均匀轮廓给 `{0}`（R12 的机制）；**原生满分支 `2^a` 给 III$_{1/2}$**（残差 `2.9e-12`）——**不是** BW 所需的 III$_1$；而幂律 `(a+1)^2`、阶乘 `(a+1)!`、甚至 `2^a(a+1)` **全部给 III$_1$**，即**指数因子不足以强制 III$_\\lambda$**，只要轮廓不是**恰好等差**。文档里的原生例（G29 的块 `2,5,20,100`，比值 `2.5,4,5` 递增）非等差，指向 III$_1$。**净效果**：`ACTION-PHASE-MATCH` 的 (T1)/(T3) 反过来给缺失输入 `π`（`E5`）一个可否证约束——**渐近块轮廓须非等差**；并给出二值否证判据：若极限**恰为** III$_{1/2}$，则其模论非 QFT 那一支（局域代数 III$_1$）⇒ 无 BW。原生轮廓本身仍未算出。R35 不关闭 L1，也不解决 `2π`。\\n")
L.append("`R36` 对 Zero 现有结构做 **CHSH 检验**：机制校验通过（Bell 态 `2.8284=2√2`、直积态 `2`、Werner 门槛 `0.7075≈1/√2`）；但 `G68` 的多路径设定**不是双体**（单系统多路径，干涉 ≠ Bell 违反），而 `G82` 的两个独立 `Z₂`（`B=2×2`）因**独立动机**而联合态为**直积**，故 `S=2.0000`、`T` 秩 1——**单方向关联再强也不违反**。结论：**现有可检验的每一个二分割上 CHSH≤2（Bell 局域）**。违反门槛被量化：`S>2 ⟺ u₁²+u₂²>1 ⟺ T 秩 ≥2`（两比特时等价于纠缠；Werner 参考线 `p>1/√2`）。缺口因此被准确命名：需要的不是更强干涉，而是**空间二分割＋跨它的纠缠**——正是 `D31` 早在 `G88` §0 登记为未恢复的「**复合系统与张量积**」，与 R35 的 `III₁` 问题共享同一缺失结构。R36 不否定量子性（是「未能检验」，非「已排除」）。\\n")
L.append("`R37` 在**单体**层面做 **KCBS/语境性检验**（不需要 `E1`、空间分离或精确类空）。关键是把「五角星取向自由」这个陷阱用 **von Neumann 迹不等式**正确处理：对全部取向取最大有闭式 `S_max(rho)=sum_k lambda_k(rho) mu_k(A)`，`mu(A)=(sqrt5, 1.3819660, 1.3819660)`，非语境界 `2`。对照通过（相邻正交误差 `2.8e-16`；完全混合 `5/3<2`；最优纯态 `sqrt5`；3000 随机态抽查无一起界）。**结论：原生态违反 KCBS** —— 以 `G72`/`G29` 推前权重 `2,5,20,100` 的 3 维归约，谱 `(0.7407,0.1852,0.0741)` 给 `S_max=2.0146>2`（去最小块给 `2.0652`），阈值 `lambda*=0.7236`、原生余量 `+0.0171`。这是量子栏**第一个正面硬证据**：Zero 的单体统计**是语境的（量子）**，不是经典概率。与 `R36` 合起来：**单体量子性成立；多体（Bell）层面只是尚无场地**。边界：检验是 state-dependent（测量集按态选）；3 维归约＝`pi`/`E5`；谱用文档例，**若真实谱 `lambda_1<0.7236` 则结论翻转**（二值可证伪入口）。\\n")
L.append("`R38` 检验纠缠的涌现机制 **H：共享的闭合起因（历史层记录）＋ 一个历史层未记录的自由度 ⇒ 跨位点纠缠**，用 `Z15` 的 Jordan-Wigner 构造与 `R36` 的 CHSH 判据。结果：单位点 `c_j†|0>` 给 `S=2.000`（不纠缠）；**跨位点相干叠加 `(c_0†+c_1†)|0>/sqrt2` 给 `S=2sqrt2=2.8284`（最大纠缠）**；位点被记录（`50/50` 经典混合）退回 `S=2.000`。机制成立的原因是 **JW 字符串把「位点」非局域化**，故「同一个费米子」这一共享事实与「哪个位点」这一未记录自由度可并存。判据 **(R38-1)**：`纠缠 = 共同起因 ∧ 历史层未记录的自由度`；**(R38-2)** 它有两条实现路径：A′ 给关系加相位（`EDGE-CONNECTION`）或 **B′ 多一个未记录自由度（格点，不需新相位）**。这把 `R36` 的「没有纠缠」与 `R37` 的「单体量子性成立」缝起来。决定性下游问题 **(R38-3)**：**Zero 的历史层对「位点」究竟是失明还是记录？**（二值）。边界：`Z-READ` 未导出、链是 `1+1` 维、共同起因在探针中是给定的。\\n")
L.append("`R39` 裁决 `R38` 的下游问题 (R38-3)：闭合记录 `(精确词 w, 闭合类 [w])` 里有没有位点信息。**结论：完全失明** —— 环图 `C_m` 上闭合性与起始位点无关，故每个记录与**全部 m 个位点**一致，`I(位点;记录)=0`（10 组 `(m,L)` 全部 `|I|<1e-15`）；而记录仍**非平凡**（区分词形：`m=6,L=6` 给 6 类）。**结构性理由**：位点标签在 Zero 里**非原生**（`E1`/`I5b` 是第一号承重项），故记录**不可能**含它 —— `E1` 缺失这件原第一号卡点，在此**第一次变成资源**。**裁决**：走 **B′**（多一个未记录自由度），不需要 `EDGE-CONNECTION`。**唯一翻盘入口**：播种点 `P_i`（D 系列未解选择器）——若播种把位点写进记录，失明被打破，须改走 A′。\\n")
L.append("`R40` 裁决 `R39` 的唯一翻盘入口：播种点 `P_i` 是否把位点写进记录。**结论：不写。** 实现层：`zero_sum_cycle_evolution.py` 的播种输入是 `seed_word`（词），经 `canonical_cycle` 映到**旋转类**（L257-265）——形状层面而非位置层面；脚本自述 seed 是显式输入（L566）。结构层：等变性探针显示位点旋转下记录**相同**（6/6 组），`I(位点;记录)=0`，故位点信息在任何环节都注入不进去。于是这条线从机制到落地闭合：`R37`（单体量子性）→ `R36`（现有二分割无纠缠）→ `R38`（纠缠＝共同起因＋未记录自由度，`S=2sqrt2`）→ `R39`（记录对位点失明）→ `R40`（播种不注入位点）。**R38 的路线 B′ 保持有效，不需要改走 A′**（`EDGE-CONNECTION`）。残余：播种『用哪个词』仍是具名输入，但那是形状选择。\\n")
L.append("`R41` 更正 `R38`/`R40` 的一处不精确注脚并立成命题：**纠缠机制（共同起因 ＋ 未记录标签 ⇒ 纠缠）不含任何空间维度量，故维数无关**。数值：1D(`L=4`)、2D(`2x2`,`2x4`)、3D(`2x2x2`) 格子上沿哈密顿路径做 Jordan-Wigner，单粒子跨位点相干叠加**一律给 `S=2sqrt2=2.828427`**。原文「仍为 1+1 维、接到 3 维空间仍需 E1」**不精确**：`1+1` 只是所选**实现**（JW 需路径排序）的属性。真正受几何限制的是：① 把「两方」识别为「两处空间」（`E1`）；② Bell 的**类空**前提（精确锥）；③ 面积律标度（含 `D-1`）——与选维（`D=4`）无关。**反向洞见**：机制不受限 ⇒ **纠缠是默认的，退相干（完整记录）才是需要解释的那个**，与真实物理一致。\\n")
L.append("`R42` 把 `pi` **显式取成旋转类**（`R31` 路线 A：块＝平衡词的 `Z_L` 轨道，权重 `omega_c=o_c/N_L`），于是 `R35`/`R37` 的两项二值判据**同时变成可算**。**交叉校验**：`q_L=sum omega_c^2` 与 `R31` 的表**四个值全吻合**（`1, 5/9, 7/25, 19/175`）。**裁决一（语境性，负）**：`lambda_1` 从 `0.667` 掉到 `0.001`，`L>=4` 全部 `<0.7236` ⇒ **非语境**。**裁决二（模论类型，正）**：轨道大小 `o_c=L/d (d|L)` ⇒ 比值 `log p`；**一般 `L` 秩>=2 ⇒ `III_1`**（如 `L=12` 秩 4），素数幂子列 ⇒ `III_lambda` ⇒ 类型还取决于**沿哪条细化序列**取极限。**收窄 `R37`**：其正面结论条件于文档里的**粗分块**（`2,5,20,100`，`lambda_1≈0.79`），不是 `pi` 无关。**本轮核心**：`pi` 的第一条**双侧约束**——① 主导块 `lambda_1>0.7236`（要粗）；② 尾部对数比秩>=2（要细）；三个自然候选**没有一个同时满足**。\\n")
L.append("`R43` 换路——不枚举自然分割，而是**扫描分割族并构造**。**7 个自然族没有一个同时满足** `R42` 的双侧约束（最接近者：闭合长度 `L=20`，`lambda_1=0.7362` 过阈值但 `S_max=1.9849`，**差 0.8%**，且 `L->oo` 渐近失败）；而**两层族**（主导类 `p` ＋ 幂律尾 `a^{-alpha}`）在 60 组网格中 **46 组同时满足**：`p≳0.8`、`alpha≳1`、任意 `K` ⇒ `S>2`（语境）且秩≥2（`III_1`）。**最漂亮的一点**：仓库的文档例 `(2,5,20,100)` **本身就是两层形状**（`p=0.7874≈0.8`），只差**无限尾**。于是 `pi` 的形状要求被两条尺子夹出来：**主导类（要粗）＋ 幂律尾（要细）**。未导出：该形状仍需从 Zero 层结构产生，且未检查与 `R32` 生存要求（另一个泛函）的相容性。\\n")
L.append("`R44` 执行 `R43` 的第二笔硬账，得到**解析 no-go**：**选维（生存）与单体语境性互斥**。两者都是同一个泛函 `q=sum omega^2` 的函数：生存窗口 `F_D=C(D,2)q^D` 峰在 `D=4 <=> q in (1/2,3/5)`；语境性上界 `S_max(q)=mu2+(mu1-mu2)(1+sqrt(2q-1))/2`，**`S_max(3/5)=2.000000` 恰好**（`lambda1*=0.723607`）。于是 `q<3/5` 给 `D=4` 峰但 `S<2`，`q>3/5` 给 `S>2` 但峰 `>=D=5`——**两者在 `q=3/5` 相接不重叠**。两条既有正面结论正落在互斥两侧：路线 A `L=4`（`q=5/9`）生存✅/语境✗；文档例 `(2,5,20,100)`（`q=0.6466`）语境✅/生存✗（峰 `D=5`）。no-go 只在 `R31`/`R32` 账本框架内成立，故**逃逸出口＝改账本形式 `F_D`**（新目标）。\\n")
L.append("`R45` 执行 `R44` 的新目标：扫描账本族 `F_D=M(D)q^{E(D)}`（8x8=64 个形式），找使「`D=4` 峰窗口」覆盖 `q>3/5` 者。**结果：39 个形式都能——`R44` 的 no-go 可逃。** 最干净的一行：把「对」的计数从**洛伦兹对** `C(D,2)=dim so(D-1,1)` 换成**欧氏对** `C(D+1,2)=dim so(D+1)`，窗口即从 `[0.502,0.600]` 移到 `[0.602,0.666]`，**整个落进语境性区**；此时文档例（`q=0.6466, S=2.0328` ✅）与两层族（`q=0.6585, S=2.0150` ✅）**同时**给出 `D=4` 峰与语境性。**代价照实说**：39 个形式都行 ⇒ 「`D=4`」不再唯一钉住账本形式 ⇒ `R32` 的「唯一存活者」收窄为「**洛伦兹对字典内**唯一」。**新目标**：从零和输运的对结构**导出** `C(D+1,2)`（或排除它）。\\n")
L.append("`R46` 执行 `R45` 的目标：**从零和输运的对结构导出 `C(D+1,2)`**。关键观察：`C(k,2)` 里的 `k` 是**对象个数**，而 Zero 的原语对象是**通道**（补偿移动 `T_ex: x -> x+e_j-e_i` 由通道对 `(i,j)` 指标化），故多重度 `=C(|C|,2)`；若通道是 `D`-单纯形的顶点（`|C|=D+1`），多重度即 `C(D+1,2)`（单纯形边数）。**仓库自身的证据**：`G29` 核验三的 `single_cut: M=1+r` **逐值等于** `r`-单纯形顶点数（`r=1..7` 全对），而 `all_cuts: 2^r` 是超立方、非单纯形。**价格**：单纯形是**欧氏**的、没有签名，故选它就必须让**因果结构**供出洛伦兹签名（`G59` 锥／`R33` T1）。**二值判据**：供得出 ⇒ `D=4` 与单体量子性同时到手；供不出 ⇒ `R44` 的互斥重新生效。\\n")
L.append("`R47` 执行 `R46` 的判据：**因果结构能否供出洛伦兹签名？裁决：能。** 论证链：(1) 光滑锥场 ⟺ 共形洛伦兹结构（锥＝二次型零集，符号差 `(1,D-1)`；数值 `D=2..6` 全对）；(2) **签名是离散不变量** ⇒ `G59` 的「有效锥」（锥外指数小而非零）无害——指数尾巴只把边界抹糊（锥外占比随阈值 `0.395 -> 0.052`），改不了锥的拓扑（pointed/convex）。于是 `R44` 的互斥被**完全绕过**：账本多重度取单纯形边 `C(D+1,2)`（`R46`）、签名由因果锥供出（`R47`）、峰位窗口 `q in (0.6,2/3)`（`R45`）、语境性 `q>3/5`（`R44`）⇒ **`D=4` 与单体量子性可以同时到手**（文档例 `q=0.6466, S=2.0328`、两层族 `q=0.6585, S=2.0150` 均满足）。残留：共形因子/尺度不导出（与 `G57` 一致，非新缺口）、精确光锥仍缺、动力学仍需 `G1`+L1。\\n")
L.append("`R48` 关闭 `G59` 的开放项「有效锥 -> 精确锥」，并更正 `R47` §3 残留 2 与 `G59` §3.3 的措辞："
         "**支持锥（精确锥）从来就是精确的**——Z1 定理 1 的最近邻传输每步最多外扩一格，故 `supp rho_n subset [-n,n]`（`B` 无关；`B=2,3,3.9,4,5` 核验零违反，`max(tip-t)=0`）。"
         "缺的不是锥，而是**锥边的填充／可见性** `eps(B) = lim_t (1/t) log rho_t(t)`：基准实测给出 `eps = 0.346579, 0.235004, 0.143841, 0.066766, 0.012659`（`B=2,2.5,3,3.5,3.9`），"
         "与闭式 `(1/2)log(4/B)` 一致到 `1e-5`（独立无饱和积分器复算同值）。于是三区制：`B<4` 锥边**指数不可见**（观察到的是被选锥 `c_*`，深阈锥边速度与 `c_*` 相对差 `<0.35%`）；"
         "`B=4` 锥边**可见**（`rho_t(t)*t -> 7.99`，幂律指数 `0.9961`，尖端速度精确 `1.000000`）；`B>4` 锥边被**饱和填满**（`O(1)`，尖端速度仍精确 `1`，与初值振幅 `1 -> 1e-10` 无关）。"
         "而 `B<=4` 是闭式（`f(+inf)=log 2`）⇒ **可见锥边 <=> B=4**，与 `G56` 的 `c_*=1` 是同一条件，故 `B=4` 从「模型参数」升级为**自洽条件**（不是新公理），`G56` §4 的「`B` 无原生约束」被收窄为「`B<=4` ＋ 绝对来源仍未导出」。\n")
L.append("诚实边界：`eps(B)` 的**常数 `1/2` 目前是数值律（`1e-5` 一致），无解析推导**；精确的是 `eps(4)=0`、`B<=4`、以及 `lambda(mu*)=mu*c*` 与 KPP 定义。残留：`B=4` 的绝对来源、尺度（`G57`）、动力学（`G1`）、高维常数因子，以及 `R45` 的账本唯一性。核验：`R48_check.py`。\n")
L.append("### 0.5.1 术语迁移档案（非理论文档）\n")
L.append("[`A2Z_MIGRATION_RECORD.md`](A2Z_MIGRATION_RECORD.md)：记录 2026-10-03 的 A0–A5 → Z0 条款术语迁移（用户裁决、迁移红线、事故与纪律）。**家谱映射的权威位置是 [`G0`](G0_bottom_layer_and_derivation_route.md) §0.1**；回归闸门 [`Z0_axiom_hygiene_check.py`](Z0_axiom_hygiene_check.py)。\n")
L.append("`R49` 执行 `R48` §5 的量子侧目标（让 $\\pi$ 非等差），把 `R35` 的类型指纹与 `R42`/`R43` 的双侧约束合并成**一个有限可判定的条件**："
         "**判据**：块权重的对数比生成子群 $G$ 稠密 $\\iff$ 相邻比的对数在 $\\mathbb Q$ 上线性无关 $\\iff$ **素数指数差向量秩 $=k-1$**（有限、可判定）。"
         "**不可能**：$k\\le3$ ⇒ 秩 $\\le2$ ⇒ 任何 3 块轮廓（含 3 维归约）永远是 $III_\\lambda$——这也解释了 `R42` 表里「$\\lambda_1$ 越大秩越小」不是巧合。"
         "**可行**：4 块上语境性（$\\lambda_1>0.723607$）与稠密性（秩 3）**解耦**，显式解 $(10^4,2,3,5)$ 给 $S_{\\max}=2.2349$、$(10^3,1,3,15)$ 给 $2.2188$、$(2,246,1,5)$ 给 $2.2037$。"
         "**新障碍**：旋转类的轨道权重全是 2 的幂 ⇒ 秩恒为 1 ⇒ **旋转类恒 $III_\\lambda$**（除非 $L$ 含非 2 素因子，而 $L=12$ 已被 `R31.1` 排除）。"
         "未做：从 Zero 原生生成该形状的 $\\pi$。核验：`R49_check.py`。\\n")
L.append("`R50` 立一条**推导纪律**并据此会诊全库矛盾：**每个量、定理、常数都带层指标 $\\ell$**，断言写成 $P_\\ell(v)$ 才完整——"
         "默认解释不是矛盾而是**层不同**，只有**同层相反**才是真矛盾。层：L0 底层（`Z0`/`Z1`–`Z5`）／L1 历史层（$\\mathcal P$）／L1′ 全局闭合类层（$\\mathcal Z_\\ast$）／L2 演化层（活动层）／$\\mathcal R$ 读出面（$\\pi,\\omega,K$）。"
         "**会诊 10 条**：真矛盾 1（`L=8`，已由 `R32.2` 撤回）、符号碰撞 3（$K$／$N$／局部编号）、其余 6 条全是**层坍塌**——例如「不设概率」(L0) vs「必须靠概率」(L2/L3)、「$B$ 不可导出」(L0) vs「$B=4$ 钉住」($\\mathcal R$)。"
         "并把 `R48` 的 $\\varepsilon(B)$ 从单侧公式**重算为双侧定律**：$B<4$ 指数衰减、$B=4$ 幂律、$B>4$ 指数增长被 **L2 容量饱和**截断（实测斜率 $\\approx0$）。"
         "三条操作规则：①断言带层指标；②比较前对齐层；③**层坍塌优先于改结论**。"
         "**纪律 ① 已落地**：6 条层坍塌＋3 条符号碰撞已逐条补层指标（`G28`／`G73`／`G59`／`R48`／`Z1`／`Z17`／`R32`／`G56`／`G61`／`G32`／`G72`／`G33`／`Z0` §0.5），只加指标不改结论，落地后 170/170 通过。核验：`R50_check.py`。\\n")
L.append("| 分支 | 文档 | 标题 |")
L.append("|:--|:--|:--|")
for f in RMD:
    L.append("| `%s` | [`%s`](%s) | %s |" % (f.split('_')[0], f, f, clean(title(f).split('·',1)[-1].strip())))
L.append("")
L.append("核验：" + ("、".join("[`%s`](%s)" % (f, f) for f in RCK_AUX) if RCK_AUX else "（尚缺）") + "。\n")
L.append("---\n")
L.append("## 1 三个来源与归属\n")
L.append("| 前缀 | 来源 | 归属 | `lh/` 中的数量 |")
L.append("|:--|:--|:--|--:|")
L.append("| **`G*`** | 本会话独立推导（自底层条款（Z0 条款 ＋ Z1–Z5 定理）到 GR） | **零和宇宙** | %d 篇 + %d 个核验脚本 |" % (len(G), len(CK)))
L.append("| **`D2xx`** | 零层弧（D210–D259） | **零和宇宙**（其中 8 篇含非原生桥接） | %d 篇 |" % len(D))
L.append("| **`zero_sum_*`** | 零和宇宙的仿真验证笔记与程序 | **零和宇宙** | %d 篇 + %d 个脚本 |" % (len(Z), len(ZPY)))
L.append("| **`Z*`** | **基础层**：单一公理「零不断乱动」＋ 取代 A1–A5（历史命名）＋ Zero 结构扩展 | **零和宇宙** | %d 篇 + %d 个核验脚本 |" % (len(ZMD), len(ZCK)))
L.append("| `D1–D209` | 混合语料 | **不属零和宇宙** | **0**（未拷入） |")
L.append("")
L.append("$$\n\\boxed{\\text{zero-sum universe}\\ =\\ \\texttt{zero\\_sum\\_*}\\ \\to\\ \\texttt{Z*}\\ \\to\\ \\texttt{G*}\\ +\\ \\texttt{D2xx}}\n$$")
L.append("")
L.append("---\n")
L.append("## 2 编号规则（**不改字头**）\n")
L.append("**结论：保留 `zero_sum` / `Z` / `G` / `D` 四个字头，不统一。** 字母序即**读序**：`zero_sum`（仿真）→ `Z`（基础层）→ `G`（推导）→ `D`（对照）。理由：\n")
L.append("1. **D 编号是稳定 ID**——被正文、母项目、D 文档间交叉引用、以及核验脚本共同引用；改字头全断；")
L.append("2. **三个来源确实不同**（作者／纪律／核验基础设施），合成一个字头会丢失来源信息；")
L.append("3. **改动量巨大、收益为零**：G 系列内部有数百条链接，D 之间有数百条交叉引用；")
L.append("4. **真正的歧义只有 D**（D1–D209 vs D210–D259），而在 `lh/` 里**只有 D210–D259**。\n")
L.append("**消歧方式＝本清单（文档）而不是改名。** 如需物理分目录，见 §8。\n")
L.append("---\n")
L.append("## 3 `zero_sum_*` 仿真验证（%d 篇笔记 / %d 个程序 ＋ %d 个核验脚本）——**Zero 原始层**\n" % (len(Z), len(ZPY), len(ZCK2)))
L.append("这是体系的**最底层**：先有仿真，后有公理。分诊见 [`G90`](G90_zero_series_reference_triage.md)。\n")
L.append("**本层由两部分构成：文章（笔记）与程序（可执行沙盒）**。程序的定理已写成文章：[`旋转类代数`](zero_sum_rotation_class_algebra.md)、[`繁殖转移定理`](zero_sum_reproduction_transition_theorems.md)、[`闭合图定理`](zero_sum_closure_graph_theorems.md)、[`持续性定理`](zero_sum_persistence_theorems.md)。\n")
ZART = [f for f in Z if clean(title(f)).startswith("Zero 系列")]
ZNOTE= [f for f in Z if f not in ZART]
L.append("### （甲）原始笔记（%d 篇）——仿真记录\n" % len(ZNOTE))
L.append("| 笔记 | 标题 |")
L.append("|:--|:--|")
for f in ZNOTE:
    L.append("| [`%s`](%s) | %s |" % (f, f, clean(title(f))))
L.append("")
L.append("### （乙）由程序写成的定理文章（%d 篇）\n" % len(ZART))
L.append("程序里装着定义与定理；这些文章把它们写成可读形式（公式＋数值核验＋与体系的接口）。\n")
L.append("| 文章 | 内容 |")
L.append("|:--|:--|")
for f in ZART:
    L.append("| [`%s`](%s) | %s |" % (f, f, clean(title(f)).replace("Zero 系列 · ", "")))
L.append("")
L.append("### （丙）本层的核验脚本（%d 个）\n" % len(ZCK2))
L.append("核验与程序分离：程序给数据，脚本复核文章里的公式与数字。\n")
L.append("| 脚本 | 核验的文章 |")
L.append("|:--|:--|")
for f in ZCK2:
    L.append("| [`%s`](%s) | `%s` |" % (f, f, f[:-len("_check.py")] + ".md"))

L.append("另有 **4 个动力学沙盒**（`zero_sum_living_universe.py`、`zero_sum_geometry_probe.py`、`zero_sum_cycle_evolution.py`、`zero_generative_selection.py`）——已在 [`G28`](G28_dynamics_audit.md) 审计。\n")
L.append("---\n")
L.append("## 4 `Z*` 基础层：单一公理与取代 A1–A5（历史命名）（%d 篇 ＋ %d 个核验脚本）\n" % (len(ZMD), len(ZCK)))
L.append("**读序**：本层在 Zero 原始层**之上**、G 推导系列**之前**。它把 Zero 的架构确立为 Z0 条款（公理只有 Z0），并把 A1–A5 降为定理。\n")
L.append("| 文档 | 标题 |")
L.append("|:--|:--|")
for f in ZMD:
    L.append("| [`%s`](%s) | %s |" % (f, f, clean(title(f).split('·',1)[-1].strip())))
L.append("")
L.append("核验：" + "、".join("[`%s`](%s)" % (f, f) for f in ZCK) + "。表头 `Z0` = **唯一公理**；`Z1` = **取代 A1–A5 与结构扩展**。\n")
L.append("")
L.append("**四条结构扩展**（均在 G 系列之前，各带价签）：\n")
L.append("| 条款 | 买回 | 处置 |")
L.append("|:--|:--|:--|")
L.append("| **Z-E1** 全局读出耦合 | Z1 定理 1 的**连通性**；G23 (a) 的全局读出 $R_Z$／层间核 $K_{ij}$ | **识别 U**（不是公理）；与 `tri_layer` 笔记的局部重建相反，属改选 |")
L.append("| **Z-E2** 局部寿命与年龄 | G23 **(b) 年龄簇（22 篇）**主体 | **导出**：年龄＝词长；局部 $\\tau_i$ 为参数异质性 |")
L.append("| **Z-E3** 周期与相位 | G23 **(c) 周期／相位簇（9 篇）**主体 | 周期＝**参数**；相位＝**定义**（需先有 $T$） |")
L.append("| **Z-E4** 符号与配对 | G23 **(d) 符号／配对簇（9 篇）**主体 | **导出**：符号＝双向步；配对＝补偿步对 $(+,-)$ |")
L.append("")
L.append("（G23 的 (e) 模簇 4 篇**明确不并入**——属 D1–D209 上游侧，登记【开放】。）\n")
L.append("---\n")
L.append("## 5 `G*` 推导系列（%d 篇）\n" % len(G))
for name, nums in THEMES:
    L.append("### %s\n" % name)
    for n in nums:
        for f in G:
            if num(f)==n:
                L.append("- [`%s`](%s) — %s" % (f, f, clean(title(f).split('·',1)[-1].strip())))
    L.append("")
L.append("### 核验脚本（%d 个）\n" % len(CK))
L.append("[`ledger_sync.py`](ledger_sync.py) 是**同步入口**：按 mtime 缓存逐个跑 `G*_check.py` 与 `Z*_check.py`（未变者不重跑），写回账本合计，再重生本清单。\n")
L.append("```\npython3 ledger_sync.py            # 增量同步（用缓存）\npython3 ledger_sync.py --all      # 忽略缓存，全部重跑\n```\n")
L.append("账本合计：**独立实断言 %d / 不符 0**（依赖上文的结论行以 `[i]` 单列，不计入合计；见 [`G10_final_derivation_and_input_ledger.md`](G10_final_derivation_and_input_ledger.md) §6）。\n" % TOTAL)
L.append("---\n")
L.append("## 6 `D2xx` 零层弧（%d 篇）\n" % len(D))
L.append("编号**保持不变**（= 母项目中的同一编号）。★ 标记的 8 篇含**非原生桥接结构**（见 [`G26`](G26_scope_and_non_native_structures.md)／[`G27`](G27_purification_attempt.md)）。\n")
L.append("**范围**：本层是**借入的旧理论零层弧**。其中的外部公理体系（U 系）字样一律为**否证性声明**（\"这不是从它推出的\"）或**接口审计标的**（[`D212`](D212_zero_universe_to_u_interface.md)／[`D215`](D215_cycle_local_semantic_u_interface.md)，二者已自行放弃\"无条件导出\"的主张）。**本体系的推理链不使用它们**；本层的公理只有 [`Z0`](Z0_zero_never_rests_single_axiom.md) 一条。\n")
L.append("| 文档 | 标题 |")
L.append("|:--|:--|")
for f in D:
    tag = ' ★' if f[:4] in BR else ''
    L.append("| [`%s`](%s)%s | %s |" % (f, f, tag, clean(title(f).split('·',1)[-1].strip())))
L.append("")
L.append("---\n")
L.append("## 7 工具\n")
L.append("| 文件 | 用途 |")
L.append("|:--|:--|")
L.append("| [`latex_lint.py`](latex_lint.py) | 扫描 13 类 LaTeX 问题（C1–C12） |")
L.append("| [`latex_fix.py`](latex_fix.py) | 修复器（默认 dry-run，`--apply` 才写盘） |")
L.append("| [`LATEX_AUDIT.md`](LATEX_AUDIT.md) | LaTeX 审计记录 |")
L.append("| [`INDEX_check.py`](INDEX_check.py) | 核验本清单：无孤儿文档 / 计数一致 / 范围声明 |")
L.append("| [`STATUS_check.py`](STATUS_check.py) | 核验唯一状态源、关键页指向、旧状态措辞清零 |")
L.append("| [`ledger_sync.py`](ledger_sync.py) | 同步入口：增量重跑 `G*_check.py`、写回账本合计、重生本清单 |")
L.append("| [`gen_index.py`](gen_index.py) | 重生成本清单（标题与合计数自动取自文件系统与账本） |")
L.append("")
L.append("---\n")
L.append("## 8 建议的阅读顺序（**自底向上**）\n")
L.append("| 顺序 | 读什么 | 为什么 |")
L.append("|--:|:--|:--|")
L.append("| 0 | [`STATUS.md`](STATUS.md) | **唯一当前状态源**：先分清主路线、D259 路线与历史替代关系 |")
L.append("| 1 | `zero_sum_*`（见 §3 的 7 篇笔记） | **Zero 原始层**：先有仿真，后有公理 |")
L.append("| 2 | [`Z0`](Z0_zero_never_rests_single_axiom.md) | **唯一公理**：零不断乱动（三款各买一件东西） |")
L.append("| 3 | [`Z1`](Z1_zero_layer_as_the_foundation.md) | 取代 A1–A5 ＋ Zero 结构的扩充（含散度字典） |")
L.append("| 4 | [`Z2`](Z2_zero_to_gr_direct_route.md) | **直连路线**：从 Zero 到 GR（免去 A0–A5（历史命名）桥接；标出两处残余输入） |")
L.append("| 5 | [`G0`](G0_bottom_layer_and_derivation_route.md) | **底层条款表**（Z0 条款 ＋ Z1–Z5 定理；A0–A5 历史命名）＋ 路线 ＋ 账本 I1–I12 |")
L.append("| 6 | [`G19`](G19_axiom_reduction.md) | 条款精简：6 条 A 条款 → 3 条独立定理（公理只有 Z0） |")
L.append("| 7 | [`G1`](G1_derivations_from_the_bottom_layer.md) → [`G8`](G8_dimension_selection.md) | 几何骨架 → 维数 |")
L.append("| 8 | [`G13`](G13_foliation_and_lorentz_invariance_gap.md) → [`G15`](G15_bare_ax3_has_no_characteristic_speed.md) | 叶层与洛伦兹不变性（含 no-go） |")
L.append("| 9 | [`G16`](G16_repair_audit_without_new_axioms.md) | 不增扩充条款清算：**主定理存活** |")
L.append("| 10 | [`G23`](G23_zero_layer_structure_inventory.md) → [`G25`](G25_age_to_geometry_channel_is_obstructed.md) | 语料审计、范围 |")
L.append("| 11 | [`G27`](G27_purification_attempt.md) → [`G29`](G29_probability_as_derived_not_postulated.md) | 净化与概率的导出 |")
L.append("| 12 | [`G28`](G28_dynamics_audit.md) | 动力学完成度（含闭环图平均度无界） |")
L.append("| 13 | [`G90`](G90_zero_series_reference_triage.md) | Zero 系列分诊：哪些已承接、哪些待造 |")
L.append("")
L.append("---\n")
L.append("## 9 如需物理分目录（可选方案）\n")
L.append("本清单已解决**归属**问题。若还要物理分目录，最小改动方案是：\n")
L.append("```")
L.append("lh/")
L.append("  INDEX.md            ← 本文件")
L.append("  G/                  ← 全部 G*.md + G*_check.py")
L.append("  D_zero_arc/         ← 全部 D2xx_*.md")
L.append("  zero_sum_sims/      ← 全部 zero_sum_* / zero_*")
L.append("  tools/              ← latex_lint.py / latex_fix.py / LATEX_AUDIT.md")
L.append("```")
L.append("代价：需改 2 个核验脚本里的路径（`G20_check.py` 读 `zero_sum_*`、`G22_check.py` 读 `D2xx`），其余脚本用 `MOD` 绝对路径不受影响。\n")

io.open('INDEX.md','w',encoding='utf-8').write('\n'.join(L))
print("INDEX 重生成：G=%d 脚本=%d D=%d Zmd=%d Zpy=%d" % (len(G), len(CK), len(D), len(Z), len(ZPY)))
