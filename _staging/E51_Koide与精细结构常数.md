# E51：Koide 与精细结构常数 —— **同源、同平面，但无代数关系**

**篇性**：**资料查证 ＋ 独立复算篇（暂存讨论稿）**。用户提问："有没有可能，这个小出公式和精细结构常数有联系？查查最旧那个库"。
**结论**：**(一)** 最旧库（`geometric-theory-geometric-theory`，2607 系列）里 $\alpha=1/S_e$，$S_e=$ **六项作用量**在电子构型处的值；**(二)** 而 $S_e^{1/2}$ **直接进入轻子质量公式** $\Longrightarrow$ **$\alpha$ 与轻子质量同源**；**(三)** **但 Koide 是纯比值，$S_e$ 在其中不抵消** $\Longrightarrow$ **与 $\alpha$ 无代数关系**（三重数值核验）；**(四)** 两条**结构性**联系：**同一平面**（$(1,1,1)^\perp$）与**同一工具**（角度与三角函数）；**(五)** **独立复算**：旧库给的 $(57.93^\circ,26.16^\circ,5.91^\circ)$ 给 $S_e=137.024934855$，**差 $-0.0081\%$** $\Longrightarrow$ **不能精确**。
**纪律**：旧库**只读**；**暂存稿**，未落盘。根目录 Z001–Z079 未改动。
**依据**：最旧库 `geometric-theory-geometric-theory/`：`Vol-5_标准模型重建/5.3_精细结构常数_CN_260713.1.md`（六项作用量；$\theta_M+\theta_C+\theta_I=90^\circ$；$P_0$）；`Vol-5_标准模型重建/5.8_轻子质量谱与PMNS矩阵_CN_260713.1.md`（带电轻子质量定理；$S_e^{1/2}$）；`Vol-0_从零开始/0.4_六项代价函数与九素互扼_CN_260713.1.md`（$S_e$ 的诚实标注）；`Vol-2_量纲桥/2.7_多路径验证_CN_260713.1.md`。其他库：`Zero-261010/sm/SMR28–34`；`GeometryAI-Mac-Build/app/articles/12.7`。本系列 E45（纯权重侧）、E47（旧库的相位律）、E50（$2\pi$）。
**核验**：[`E51_alpha_koide_check.py`](E51_alpha_koide_check.py)。

---

## 摘要

$$
\boxed{\;\textbf{（一）}\ \text{最旧库}:\quad \alpha=\frac1{S_e},\qquad S_e=S(\theta_M,\theta_C,\theta_I)\ (\textbf{六项作用量}) .\;}
\tag{0.1}
$$

$$
\boxed{\;\textbf{（二）}\ S_e^{1/2}\ \textbf{进入} \ \text{轻子质量} \Longrightarrow \alpha\ \text{与轻子质量}\ \textbf{同源} .\;}
\tag{0.2}
$$

$$
\boxed{\;\textbf{（三）}\ \text{但 Koide 比值}\ \textbf{不含}\ S_e \Longrightarrow \textbf{无代数关系} .\;}
\tag{0.3}
$$

$$
\boxed{\;\textbf{（四）}\ \text{两条结构联系}:\ \textbf{同一平面}\ ((1,1,1)^\perp)\ +\ \textbf{同一工具} .\;}
\tag{0.4}
$$

$$
\boxed{\;\textbf{（五）}\ \textbf{独立复算}:\ \text{旧库角度给}\ S_e\ \textbf{差}\ -0.0081\%\ \Longrightarrow \textbf{不能精确} .\;}
\tag{0.5}
$$

---

## §1 最旧库的 $\alpha$（**核验 §A**）

$$
\textbf{公理 3（最旧库）逐字}:\quad \theta_M+\theta_C+\theta_I=90^\circ\qquad(\text{Hopf 纤维化的旋转对称性}) .
\tag{1.1}
$$

$$
\textbf{六项作用量逐字}:\quad
S(\theta_M,\theta_C,\theta_I)=\frac{1}{\sin^2\theta_M}+\frac{1}{\sin^2\theta_C}+\frac{1}{\sin^2\theta_I}
+\frac{1}{\sin\theta_M\sin\theta_C}+\frac{1}{\sin\theta_M\sin\theta_I}+\frac{1}{\sin\theta_C\sin\theta_I} .
\tag{1.2}
$$

$$
\textbf{逐字}:\quad \alpha=\frac1{S_e};\qquad S_e=S(P_0);\qquad P_0=(\theta_M^0,\theta_C^0,\theta_I^0)=(57.93^\circ,\ 26.16^\circ,\ 5.91^\circ) .
\tag{1.3}
$$

$$
\textbf{逐字}:\quad \text{"}S_e=137.035999084\ \text{由九素互扼定理（GT-0.4.0.1）从几何第一原理唯一确定"} .
\tag{1.4}
$$

$$
\textbf{本文独立复算（核验 §A）}:\quad S(57.93^\circ,26.16^\circ,5.91^\circ)=\mathbf{137.024934855} .
\tag{1.5}
$$

$$
\textbf{核验 §A}:\quad \text{差}=S-137.035999084=\mathbf{-0.011064}\ (\mathbf{-0.0081\%}) .
\tag{1.6}
$$

$$
\textbf{核验 §A}:\quad \text{六个分项}=1.3926,\ 5.1447,\ 94.3216,\ 2.6767,\ 11.4609,\ 22.0285\ (\text{和}=137.0249) .
\tag{1.7}
$$

$$
\boxed{\;\textbf{核验 §A}:\ \text{旧库给的角度}\ \textbf{不能精确} \text{给出}\ S_e \text{—— 差}\ 0.008\% .\;}
\tag{1.8}
$$

$$
\textbf{核验 §A}:\quad \text{微调可到}\ 9.8\times10^{-5}:\quad (56.44^\circ,\ 27.69^\circ,\ 5.87^\circ)\Longrightarrow S=137.036096974 .
\tag{1.9}
$$

$$
\textbf{诚实标注（最旧库逐字）}:\quad \text{"}S_e\ \text{的精确数值将由第 2 卷（量纲桥）中的 Dixmier 迹重建完成"} .
\tag{1.10}
$$

---

## §2 **关键**：$\alpha$ 与轻子质量**同源**（**核验 §B**）

$$
\textbf{5.8 带电轻子质量定理逐字}:\quad
\begin{cases}
m_e=K\sin^3\theta_M^e=K\sin^3 57.93^\circ=511.0\ \mathrm{keV}\\[1mm]
m_\mu=\dfrac12K\Lambda S_e^{1/2}\sin^3\theta_M^\mu=\dfrac12K\Lambda S_e^{1/2}\sin^3 31.94^\circ=109.17\ \mathrm{MeV}\\[1mm]
m_\tau=8K\Lambda S_e^{1/2}\sin^3\theta_M^\tau=8K\Lambda S_e^{1/2}\sin^3 31.94^\circ=1.747\ \mathrm{GeV}
\end{cases}
\tag{2.1}
$$

$$
\boxed{\;\textbf{核验 §B}:\ m_\mu\ \text{与}\ m_\tau\ \textbf{都含}\ S_e^{1/2}=\alpha^{-1/2} \Longrightarrow \text{最旧库里}\ \alpha\ \text{与轻子质量}\ \textbf{同源} .\;}
\tag{2.2}
$$

$$
\textbf{而}:\quad m_e\ \textbf{不含}\ S_e \Longrightarrow \text{三代}\ \textbf{不同源} \text{（这正是最旧库的弱点）} .
\tag{2.3}
$$

$$
\textbf{最旧库自报偏差逐字}:\quad m_\mu\ +3.3\%;\qquad m_\tau\ -1.7\%\qquad(\textbf{远差于}\ 12.7\ \text{的}\ +0.0010\%/+0.0070\%) .
\tag{2.4}
$$

$$
\textbf{本文独立复算（核验 §B）}:\quad
\begin{cases}
m_e=K\sin^3 57.93^\circ=\mathbf{511.0\ keV}\ (\text{与旧库相符})\\
m_\mu=\dfrac12K\Lambda S_e^{1/2}\sin^3 31.94^\circ=\mathbf{109.1632\ MeV}\ (\text{偏}\ \mathbf{+3.32\%})\\
m_\tau=8K\Lambda S_e^{1/2}\sin^3 31.94^\circ=\mathbf{1746.6\ MeV}\ (\text{偏}\ \mathbf{-1.70\%})
\end{cases}
\tag{2.4a}
$$

$$
\boxed{\;\textbf{核验 §B}:\ \text{复算}\ \textbf{逐条复现} \ \text{旧库自报的}\ +3.3\%/-1.7\% \Longrightarrow \text{旧库在这一步}\ \textbf{自洽} .\;}
\tag{2.4b}
$$

---

## §3 **但 Koide 比值不含 $S_e$** —— 故**无代数关系**（**核验 §C**）

$$
\text{Koide}:\quad Q=\frac{m_e+m_\mu+m_\tau}{(\sqrt{m_e}+\sqrt{m_\mu}+\sqrt{m_\tau})^2}=\frac23\qquad(\textbf{纯比值}) .
\tag{3.1}
$$

$$
\text{而}\ m_\mu\propto S_e^{1/2},\ m_e\propto S_e^{0} \Longrightarrow \textbf{比值里}\ S_e\ \textbf{不抵消} .
\tag{3.2}
$$

$$
\textbf{核验 §C（三重）}:\quad
\alpha\cdot Q=91.3564892402;\qquad
\alpha\cdot\frac29=30.4524442409;\qquad
\alpha\cdot\theta=317.4609849104 .
\tag{3.3}
$$

$$
\boxed{\;\textbf{核验 §C}:\ \text{没有一个是特殊数} \Longrightarrow \textbf{Koide 与}\ \alpha\ \textbf{无代数关系} .\;}
\tag{3.4}
$$

$$
\textbf{核验 §C}:\quad \frac{\alpha}{137}=1.0002627670;\qquad \frac{\theta S_e}{137}=2.3172334665 \Longrightarrow \text{亦无} .
\tag{3.5}
$$

---

## §4 **两条结构性联系**

### 4.1 同一平面

$$
\textbf{12.7 §5.5 逐字}:\quad \theta\ \text{住在垂直于民主方向}\ (1,1,1)\ \text{的}\ 2\text{-平面};\qquad
T\Sigma=(1,1,1)^\perp .
\tag{4.1}
$$

$$
\textbf{最旧库的公理 3}:\quad \theta_M+\theta_C+\theta_I=90^\circ\qquad(\textbf{也是一个"和固定"}) .
\tag{4.2}
$$

$$
\boxed{\;\textbf{两处的"和固定"是同一件事} \Longrightarrow \text{Koide 的}\ \theta\ \text{与扇区角住在}\ \textbf{同一个几何} .\;}
\tag{4.3}
$$

### 4.2 同一工具

$$
\text{最旧库的}\ S_e=\sin\ \text{的函数};\qquad \text{Koide 的}\ \theta=\cos\ \text{的参数} \Longrightarrow \textbf{都用量角与三角函数} .
\tag{4.4}
$$

$$
\Longrightarrow\quad \text{联系是}\ \textbf{工具的} \text{，}\textbf{不是代数的} .
\tag{4.5}
$$

---

## §5 三库对 $\alpha$ 的**三种接法**（对照）

| 库 | 时期 | $\alpha$ 的接法 | 味结构 |
|:--|:--|:--|:--|
| `geometric-theory` | 260713 | $\alpha=1/S_e$（**角度函数**）；$m_\mu\propto S_e^{1/2}$ | $m=K\sin^3\theta_M$ |
| `Zero-261010` / `cosmos-construct` | 2609–2610 | $L=c\,\alpha_{\rm em}^{-1}$（**线性**） | Koide $Q=k_0/\Lambda$ |
| **Zero**（本系列） | 现在 | **未接**（$\pi$ 是具名输入，未用于味） | Koide：形状已证；$\theta$ 开放 |

$$
\boxed{\;\textbf{故三库用}\ \textbf{不同方式} \text{把}\ \alpha\ \text{接进味结构}:\ \text{最旧库走}\ \textbf{"量纲"} \text{，新库走}\ \textbf{"参数"} \text{，Zero}\ \textbf{尚未接} .\;}
\tag{5.1}
$$

$$
\textbf{而 Zero 的缺口（E45–E50）}:\quad \theta=k_0/\Lambda^2=2/9\ \text{（候选）};\quad \text{且须一个}\ \textbf{连续周期结构} .
\tag{5.2}
$$

$$
\textbf{⚠ 本文一处判断的修正（见 E54 §5）}:\quad \text{本文 §3 说}\ \alpha\ \text{是"}\textbf{代数接口} \text{"——}\textbf{只对一半}:
\tag{5.2a}
$$

$$
\textbf{对}:\ S_e\ \text{本身确是}\ \sin\ \text{的代数函数};\qquad
\textbf{错}:\ \text{其自变量（角度）须先被确定，而"和}=\tfrac1D\ \text{转}\ \textbf{要用到}\ 2\pi .
\tag{5.2b}
$$

$$
\Longrightarrow\quad \alpha\ \text{的接口}=\text{"}\textbf{代数}+\text{一个}\ 2\pi\ \text{读法"}\ \text{，}\textbf{不是纯代数} .
\tag{5.2c}
$$

$$
\Longrightarrow\quad \textbf{一个可走的方向}:\quad \text{把}\ \alpha\ \text{接进 Zero 的味结构（}\text{若}\ L=c\,\alpha^{-1}\ \text{在 Zero 里成立}） .
\tag{5.3}
$$

---

## §6 不许主张

1. **不主张** Koide "与 $\alpha$ 有联系"（代数意义）——**已否**（§3.4）；
2. **不主张** $\alpha^{-1}=S_e$ "已被导出"——最旧库自报**待第 2 卷**（§1.10）；
3. **不主张**"同源"就等于"有关系"——本文区分**同源**（同一 $S_e$）与**代数关系**（无）；
4. **不主张**旧库的角度"错了"——本文只说它**不能精确**给 $S_e$（差 $0.008\%$）；
5. **不主张**三库的接法"等价"——它们是**不同**的（§5.1）；
6. **不主张** Zero "应当"接 $\alpha$——本文只指出这是一条**未走的方向**；
7. **不主张**本文改动了任何库或正式篇——**暂存稿**，旧库**只读**。

---

## §7 结论

$$
\boxed{\;
\begin{aligned}
&\textbf{最旧库}:\quad \alpha=\frac1{S_e},\quad S_e=\text{六项作用量};\quad m_\mu,m_\tau\propto S_e^{1/2} .\\[1mm]
&\textbf{故}\ \alpha\ \text{与轻子质量}\ \textbf{同源} \text{—— 但这只在}\ \textbf{绝对值} \text{层面} .\\[1mm]
&\textbf{而 Koide 是纯比值}:\ S_e\ \textbf{不抵消} \Longrightarrow \textbf{无代数关系} \text{（三重核验）} .\\[1mm]
&\textbf{两条结构联系}:\ \textbf{同一平面}\ (1,1,1)^\perp\ +\ \textbf{同一工具} \text{（角度与三角）} .\\[1mm]
&\textbf{独立复算}:\ \text{旧库角度给}\ S_e\ \textbf{差}\ -0.0081\% \Longrightarrow \textbf{不能精确} .
\end{aligned}
}
$$

**一句话：**

$$
\boxed{\;\textbf{$\alpha$ 与 Koide 不是代数上的同一个东西，但它们确实同源 —— 都住在同一个 $(1,1,1)^\perp$ 平面上，且最旧库里 $\alpha^{-1}$ 直接乘进轻子质量。}\;}
$$

---

## 附录 A　核验

`E51_alpha_koide_check.py`：纯标准库；不 import 项目脚本、不写文件（旧库只读）。

| 组 | 内容 | 断言 |
|:--|:--|--:|
| A | 六项作用量的复算；旧库角度给 $S_e$ 的差；微调到 $9.8\times10^{-5}$ | 10 |
| B | 轻子质量定理含 $S_e^{1/2}$；$m_e$ 不含；旧库自报偏差 | 8 |
| C | Koide 与 $\alpha$ 的三重数值核验（无特殊数） | 8 |
| D | 同一平面（$(1,1,1)^\perp$）；同一工具 | 6 |
| E | 三库三种接法的对照 | 6 |
| F | 旧库文件与引文存在性（只读） | 6 |
| | **合计** | **46** |

**核验纪律**：数值只算旁证；旧库**只读**；**本文更正旧库一处数值**（§1.8）。

---

*查证 ＋ 复算稿 · 暂存 · 未改动任何公理 · 旧库只读 · 完*
