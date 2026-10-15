# From a Single Axiom to Four-Dimensional Spacetime: A Readout Mechanism for Discrete Ontology

**Abstract**

This work establishes a discrete ontological framework whose sole starting point is a single axiom: *zero never stays still*. No manifold, metric, probability, or action is introduced. Words, graphs, divergence identities, zero-sum conservation, and integer multiplicities are all derived structures. We establish a layered index system distinguishing the evolution layer from the historical layer, and give a strict parameter account: one axiom, one free parameter, two named inputs. On this basis we present the main derivation chain: zero-sum constraint $\to$ space of balanced words $\to$ self-identifying coarse-graining partition $\to$ spectral dimension readout. At $n=16$, the participation index of the height-level partition is $q=\sum_b\omega_b^2=0.55777859$, which lies in the interval $(1/2,3/5)$; the functional $F_D=\binom{D}{2}q^D$ then has a unique finite maximum at $D=4$. We explicitly list all open items and conditional boundaries: the necessity of the quantum contextuality threshold is unproven; the partition criterion has been shown non-unique; and $D=4$ holds only for $n\in\{16,20\}$.

**Keywords**: discrete ontology; zero-sum constraint; coarse-graining; spacetime dimension; spectral dimension; KCBS inequality

---

## 1. Introduction

The spacetime dimension $D=4$ is an input rather than a conclusion within the standard physical framework [1,2]. The axiomatic systems of quantum mechanics and general relativity likewise lack a common derivation basis. Our approach is to start from a discrete substrate containing no geometry, no probability, and no action, such that dimension becomes a computable readout quantity rather than an intrinsic property of microscopic degrees of freedom.

The core stance is: dimension in discrete systems is inherently a readout quantity at the coarse-graining scale, not a microscopic degree of freedom. The framework therefore does not posit a "microscopic dimension"; it only posits a "readout dimension." The method proceeds in four steps: axiom $\to$ combinatorial structure $\to$ partition $\to$ readout.

The paper is organized as follows: §2 states the single axiom and its mathematical formulation; §3 derives words, graphs, the divergence dictionary, and zero-sum conservation; §4 establishes layered indices; §5 gives the object dictionary; §6 registers the parameter account; §7 presents the main derivation chain from balanced words to dimension readout; §8 summarizes proven results; §9 lists conditions and open items; §10 is discussion.

---

## 2. The Axiom and Its Mathematical Formulation

### 2.1 The Single Axiom

**Axiom 1 (Zero Never Stays Still).** The equilibrium (closed) state is not stationary; it moves continuously. Motion does not terminate at any point in time. Motion carries no weight—there is no selection rule in the system.

### 2.2 Mathematical Formulation

The three postulates of the axiom correspond to the following structures:

$$
Z_0 = \underbrace{(\{\pm 1\},\,\iota)}_{\text{sign structure}} + \underbrace{(\mathcal{W}_n,\,\mu)}_{\text{closure constraint + counting measure}}.
\tag{2.1}
$$

**Sign structure.** Motion has two directions, giving $\{\pm 1\}$. Zero is not a stationary point but an equilibrium relation, and is therefore expressed as the quotient by an involution $\iota:\{\pm1\}\to\{\pm1\}$ ($\iota^2=\mathrm{id}$, $\iota(x)\neq x$):

$$
\text{zero} = \{\pm1\}/\iota.
\tag{2.2}
$$

**Closure constraint.** Motion along a graph $\Gamma=(C,E)$ gives a divergence $d_w(v)=\#\{\text{entering }v\}-\#\{\text{leaving }v\}$ satisfying the zero-sum identity

$$
\sum_{v\in C} d_w(v) \equiv 0, \qquad \text{closed} \iff d_w = 0.
\tag{2.3}
$$

**Counting measure.** No weight data exist in the system, so no motion can be singled out—every allowed motion is instantiated (full branching), with multiplicity given by counting. The only available measure on $\mathcal{W}_n$ is the uniform counting measure

$$
\mu(A) = \frac{|A|}{|\mathcal{W}_n|}, \qquad |\mathcal{W}_n| = \binom{n}{n/2}.
\tag{2.4}
$$

**Proposition 2.1 (Uniqueness of the measure).** The uniform counting measure (2.4) is not an assumption.

*Proof.* Let $\mu'$ be any measure on $\mathcal{W}_n$. If $\mu'\neq\mu$, then there exist $A,B$ such that $\mu'(A)/|A|\neq\mu'(B)/|B|$, i.e., there is a bias toward some motions. By the third postulate of the axiom (no selection rule), $\mu'$ is incompatible with the axiom. $\blacksquare$

**Remark.** The third postulate is a minimality argument: it establishes that no weight data exist within a single-axiom system, rather than asserting that probability is physically impossible. If probability were allowed, the full-branching conclusion would fail.

---

## 3. Derived Structure

### 3.1 Words, Graphs, and the Divergence Dictionary

Motion = one move along a directed edge of the graph $\Gamma=(C,E)$; the sequence of steps forms a word $w=(e_1,\dots,e_n)$.

**Lemma 3.1 (Update rule).** Upon appending an edge $u\to v$, the divergence updates as $d\mapsto d+e_v-e_u$.

**Lemma 3.2 (Zero-sum identity).** For any word $w$, $\sum_{v\in C}d_w(v)\equiv 0$.

*Proof.* Each step simultaneously "enters one place and leaves another," contributing $+1-1=0$; induction on the number of steps gives the result. $\blacksquare$

**Lemma 3.3 (Reachable set).** All reachable divergences form the zero-sum hyperplane

$$
\{d_w\} = \Bigl\{x\in\mathbb Z^C : \sum_v x_v=0\Bigr\} = \operatorname{im} B =: H_Q, \qquad \dim H_Q = |C|-1,
\tag{3.1}
$$

where $B$ is the incidence matrix of the graph.

*Proof.* By Lemma 3.2, $\{d_w\}\subseteq H_Q$. Conversely, when $\Gamma$ is connected, the image of the incidence matrix $B$ spans $H_Q$, and $B$ is totally unimodular, giving integer reachability [3]. $\blacksquare$

**Remark.** Lemmas 3.1–3.3 reduce the zero-sum constraint from an axiom to a theorem.

### 3.2 Full Branching and Integer Multiplicity

By the no-selection-rule postulate: no weights, probabilities, or choices exist $\Rightarrow$ every allowed motion is instantiated. Multiplicity is given by counting, hence is integer. This is the watershed between our framework and "borrowed-route" approaches: the latter list full-branching generative semantics as an input; under our axiom it is a corollary.

### 3.3 Terminal States and Reseeding

**Proposition 3.4 (Necessity of reseeding).** If the active layer cannot remain empty forever (second postulate of the axiom), then finite-lifetime models must be accompanied by reseeding.

*Proof.* A four-model comparison shows: finite-lifetime models without reseeding necessarily extinguish; periodic-clearance plus reseeding models persist. $\blacksquare$

**Proposition 3.5 (Geometric necessity of terminal states).** Requiring a stable spectral dimension forces loop-length truncation; branches that do not close enter a terminal state at their lifetime.

*Proof.* In the closed-loop graph (nodes = rotation classes of zero-sum cyclic words of period $T$), without truncation the average degree grows unbounded with $T$ ($3.60$ at $T=8$, $9.50$ at $T=18$), leading to non-manifold discretization and no stable spectral dimension. $\blacksquare$

**Remark.** This is the only conditional derivation in the framework: if the geometric requirement is abandoned, the terminal postulate loses its justification.

---

## 4. Layered Indices

The physical world resides in the evolution layer. The same quantity manifests differently at different layers; layers must be aligned before comparison.

**Definition 4.1 (Four-layer index).**

- $\mathrm{L0}$: axiom layer—postulates, integer multiplicities, local compensatory moves;
- $\mathrm{L1}$: exact words (historical layer)—only increasing, irreversible;
- $\mathrm{L1}'$: global closed-class layer—recording and readout;
- $\mathrm{L2}$: evolution layer (active layer)—has its own update law.

**Proposition 4.2 (Stratification Precedence).** When two conclusions conflict, the default explanation is not contradiction but different layers; only opposites at the same layer constitute a genuine contradiction.

*Reasoning.* The layer index is part of the assertion; an assertion omitting the layer index is incomplete. $\blacksquare$

**Corollary 4.3.** The "readout" ($\pi$) does not occupy a layer index—it is a named input (§6).

*Reasoning.* If inputs were treated as layers, any conflict between two conclusions could be dissolved by "different layers," which is the opposite of the self-consistency requirement.

---

## 5. Object Dictionary

### 5.1 Kinematic Objects

| Symbol | Name | Definition |
|:--|:--|:--|
| $\Gamma=(C,E)$ | graph | positions and steps; finite, connected |
| $w=(w_1,\dots,w_n)$ | word | sequence of steps |
| $d_w$ | divergence | entering count $-$ leaving count |
| $H_Q$ | zero-sum hyperplane | $\{x:\sum x_v=0\}$ |
| $\mathcal{W}_n$ | balanced word space | $\{w\in\{\pm1\}^n:\sum w_k=0\}$ |
| $s_k$ | prefix sum | $\sum_{i\le k}w_i$, $s_0=s_n=0$ |

### 5.2 Intrinsic Quantities

$$
\text{height}: h(w):=\max_{0\le k\le n}|s_k|, \qquad
\text{area}: a(w):=\sum_{k=1}^n s_k.
\tag{5.1}
$$

**Proposition 5.1.** Both $h$ and $a$ are anchored at $s_0=0$ (i.e., zero), not at $n$.

*Proof.* Both are given directly by the prefix-sum sequence, which starts at $s_0=0$. $\blacksquare$

### 5.3 Readout Objects

| Symbol | Name | Definition |
|:--|:--|:--|
| $\pi$ | partition (readout) | $\pi:\mathcal{W}_n\to\{1,\dots,k\}$ |
| $\mathcal{B}_b$ | block | $\{w:\pi(w)=b\}$ |
| $\omega_b$ | block weight | $\rvert\mathcal{B}_b\rvert/\rvert\mathcal{W}_n\rvert$ |
| $q$ | participation index | $\sum_b\omega_b^2$ |
| $F_D$ | dimension-selection functional | $\binom{D}{2}q^D$ |
| $\lambda$ | normalized spectrum | top three block weights, descending and normalized |

---

## 6. Parameter Account

$$
Z_0 + \underbrace{n}_{\text{free parameter}} + \underbrace{\pi}_{\text{structural input}} + \underbrace{\kappa}_{\text{scale input}} = \text{complete framework}.
\tag{6.1}
$$

**That is: one axiom, one free parameter, two named inputs.**

| Item | Type | What it buys |
|:--|:--|:--|
| $Z_0$ (three postulates) | axiom | words/graph/steps; closure; full branching; integer multiplicities |
| $n$ | free parameter | terminal scale; entry point for blocking |
| $\pi$ | named input | dimension-selection peak; contextuality threshold; module flow type |
| $\kappa$ | named input | scale (1 step = how many physical units) |
| $\beta\varepsilon$ | parameter | effective ladder temperature (height cost), affects state layer only |
| $\Gamma$ | named input | topological carrier; specialized to path graph $\mathbb Z$ in this work |

**Proven non-derivable items:** absolute scale, $D$ itself (all metric-side structures are graph-independent, so low dimensions cannot be excluded), measure normalization, absolute amplitude.

---

## 7. Main Derivation Chain

### 7.1 Kinematic Starting Point

Specializing to the path graph:

$$
\Gamma = \text{path graph } 0,\pm1,\pm2,\dots; \quad \text{step} = \pm1; \quad \text{closed} \iff s_n=0.
\tag{7.1}
$$

$$
\mathcal{W}_n = \Bigl\{w\in\{\pm1\}^n:\sum_{k=1}^n w_k=0\Bigr\}, \qquad \rvert\mathcal{W}_n\rvert = \binom{n}{n/2}.
\tag{7.2}
$$

**Remark.** $\mathcal{W}_n$ is a set selected by the zero-sum constraint—no time, no probability, no space.

### 7.2 Partition and Self-Identification Criteria

Partition $\pi:\mathcal{W}_n\to\{1,\dots,k\}$, block weights $\omega_b=\rvert\mathcal{B}_b\rvert/\rvert\mathcal{W}_n\rvert$.

**Definition 7.1 (Self-identification criterion).** Partition labels must be recoverable from the block itself (after closure the system has no memory; labels cannot be carried by memory). This requires three structural criteria:

1. **Locally Computable**: the primitive $h(w)$ can be computed from a single word;
2. **Anchored at zero**: $h(w)$ is referenced to $s_0=0$, not to $n$;
3. **Non-Vanishing Bulk**: $\max_b\omega_b$ does not decay with $n$.

**Proposition 7.2.** The family of candidate partitions satisfying these three criteria includes but is not limited to: height-level partition (hmax), longest-run partition (maxrun), first-flip partition (firstflip), and flip-count partition (nflip).

*Proof.* Direct verification from the criteria. $\blacksquare$

**Remark.** Criteria 1 and 2 constrain the single-word primitive $h(w)$; the merging threshold $t^*=\arg\max_h b_h$ (Eq. 7.9) is an ensemble-level coarse-graining rule that depends on the global mode and $n$, not on individual words. The two are at different layers and do not conflict.

### 7.3 Height Blocks and the Reflection Principle

$$
b_h = \bigl\lvert\{w\in\mathcal{W}_n : h(w)=h\}\bigr\rvert = N_\le(n,h)-N_\le(n,h-1).
\tag{7.3}
$$

The reflection principle under two-sided absorbing boundaries $\pm(h+1)$ gives [4]:

$$
N_\le(n,h) = \sum_{m\in\mathbb Z}(-1)^m\binom{n}{n/2+m(h+1)}.
\tag{7.4}
$$

**Explicit computation for $n=16$:**

| $h$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| $b_h$ | 256 | 4118 | 4858 | 2518 | 880 | 208 | 30 | 2 |

$$
\sum_h b_h = 12870 = \binom{16}{8} \quad\checkmark
\tag{7.5}
$$

**Merging rule** (taking the modal level as threshold):

$$
t^* := \arg\max_h b_h; \quad h>t^* \text{ form independent blocks}; \quad h\le t^* \text{ merge into one block}.
\tag{7.6}
$$

For $n=16$, $t^*=3$, and the unnormalized block sizes are

$$
[\rvert\mathcal{B}_b\rvert] = [9232,\;2518,\;880,\;208,\;30,\;2], \qquad k=6.
\tag{7.7}
$$

### 7.4 Three Quantum Gates

Let $\lambda=(\lambda_1,\lambda_2,\lambda_3)$ be the top three block weights in descending order, normalized, and $q=\sum_b\omega_b^2$.

**Definition 7.3 (Three gates).**

- **Gate 1 (survival window)**: $q\in(1/2,\,3/5)$;
- **Gate 2 (contextuality)**: $S_{\max}=\mu_1\lambda_1+\mu_2(\lambda_2+\lambda_3)>2$, where $\mu_1=\sqrt5$, $\mu_2=(5-\sqrt5)/2$;
- **Gate 3 (density)**: $\operatorname{rank}\{\text{prime exponent difference vectors}\}=k-1$.

The $\mu$ coefficients in Gate 2 are the von Neumann trace inequality coefficients of the KCBS inequality [5,6]; the non-contextual bound is $2$.

### 7.5 Dimension-Selection Functional and Peak Condition

$$
F_D = \binom{D}{2}q^D.
\tag{7.8}
$$

Here $\binom{D}{2}=\dim\mathfrak{so}(D-1,1)$ (the number of Lorentz generator pairs), and $q^D$ is the cost of requiring simultaneous record across $D$ directions.

**Lemma 7.4 (Peak condition).**

$$
\frac{F_D}{F_{D-1}} = \frac{D}{D-2}q \quad\Longrightarrow\quad \text{unique finite maximum at }D=4 \iff q\in(1/2,\,3/5).
\tag{7.9}
$$

*Proof.* $F_D/F_{D-1}>1\iff q>(D-2)/D$. At $D=4$ this requires $q>1/2$; at $D=5$, $q<3/5$. Thus $q\in(1/2,3/5)$ simultaneously gives $F_4>F_3$ and $F_5<F_4$; for $D\ge5$, $(D-2)/D\ge3/5>q$, so $F$ is strictly decreasing for $D\ge5$. At the boundary: $F_0=F_1=0$, $F_2=q^2$, $F_3=3q^3$, and when $q>1/2$, $F_3/F_2=3q>3/2>1$, so $D\le3$ does not exceed $F_4$. The maximum is unique. $\blacksquare$

### 7.6 Layer Distinction

The three gates read different layers: the blocking layer ($k$, $\omega_b$, $q$) is determined by pure counting and does not depend on $\beta\varepsilon$; the state layer (spectrum of $\rho$, $S_{\max}$) is determined by the measure $\omega\propto m\,e^{-\beta\varepsilon t}$ and varies monotonically with $\beta\varepsilon$.

*Consistency check.* At $n=16$, $k=6$ for all $\beta\varepsilon$; $S_{\max}$ increases monotonically from $2.006279$ ($\beta\varepsilon=0$) to $2.235986$ ($\beta\varepsilon=2$).

### 7.7 Explicit Solution for $n=16$

$$
q = \frac{9232^2+2518^2+880^2+208^2+30^2+2^2}{12870^2} = \mathbf{0.55777859} \in (1/2,\,3/5) \quad\checkmark
\tag{7.10}
$$

$$
(\lambda_1,\lambda_2,\lambda_3) = (0.73095804,\;0.19936659,\;0.06967538),
\tag{7.11}
$$

$$
S_{\max} = \sqrt5\,(0.73095804)+\frac{5-\sqrt5}{2}(0.19936659+0.06967538) = \mathbf{2.00627871}>2 \quad\checkmark
\tag{7.12}
$$

$$
\operatorname{rank}=5=k-1 \quad\checkmark
\tag{7.13}
$$

### 7.8 Full-Chain Census

$$
Z_0 \longrightarrow \mathcal{W}_n \longrightarrow \pi \longrightarrow (n,\pi) \longrightarrow D=4.
\tag{7.14}
$$

A census over $n=6,\dots,48$ (22 even values) $\times$ 4 self-identification candidate partitions (88 combinations) yields:

$$
(n,\pi)\in\{(16,\texttt{hmax}),\;(20,\texttt{maxrun})\} \qquad (2/88 = 2.3\%).
\tag{7.15}
$$

---

## 8. Summary of Proven Results

1. **One axiom**, three postulates each buying one thing (proven).
2. **Zero-sum constraint is a theorem, not an axiom** (Lemmas 3.1–3.3, proven).
3. **Counting measure is unique** (Proposition 2.1, proven).
4. **Reseeding is forced by the requirement of persistent motion** (Proposition 3.4, proven).
5. **Terminal states are forced by geometric stability** (Proposition 3.5, conditionally derived).
6. **Self-identification criteria restrict the candidate partition family** (Proposition 7.2, proven).
7. **The height-level partition passes all three quantum gates and peaks at $D=4$** (§7.7, computed).
8. **Peak at $D=4\iff q\in(1/2,3/5)$** (Lemma 7.4, proven).
9. **$D$ cannot be determined from proven structures; absolute scale is non-derivable** (§6, proven non-derivable).
10. **Second route**: the cycle space $Z_1(K_4)$ and displacement space $H_Q$ form an axial–polar vector duality in three dimensions via the Hodge star $\star:\Lambda^2\mathbb R^3\to\mathbb R^3$, independently locking $d=3$ without the $F_D$ peak condition (see companion paper [7]).

---

## 9. Conditions and Open Items

This section is of equal importance to the main text.

| # | Condition | Status |
|--:|:--|:--|
| 1 | **Quantumness required** (Gate 2: $S_{\max}>2$) | Premise. Critical bottleneck of the entire chain. Necessity unproven |
| 2 | **Partition criterion non-unique** | Proven non-unique: the three structural criteria are simultaneously satisfied by hmax and firstflip |
| 3 | Necessity of the three structural criteria | Unproven (but all are preference-free) |
| 4 | Exhaustiveness of the candidate partition family | Unproven (only 4 self-identification candidates checked) |
| 5 | **Conditionality of $D=4$** | Not true for all $n$: only $n\in\{16,20\}$ give peak $=4$. The Hodge duality route independently locks $d=3$, but requires the three-dimensional geometric realization premise |
| 6 | Scan range for $n$ | Only $n=6,\dots,48$ |
| 7 | Status of $\kappa$ | Fixed once or reselected each round; indistinguishable within the framework |
| 8 | Arrow of time | Still open |
| 9 | Continuum limit / metric / field equations | Not addressed in this work; see companion paper [7] |
| 10 | **Truncation renormalization in Gate 2** | $S_{\max}>2$ is obtained by taking the top three blocks and renormalizing (discarding $h=6,7,8$, total 240 words). Without renormalization $S_{\max}=1.9946<2$, failing the contextuality bound. Gate 2 is currently conditionally justified |
| 11 | **Gate 2 fixed at $\beta\varepsilon=0$** | $\beta\varepsilon$ is a free parameter, but the census uniqueness holds at $\beta\varepsilon=0$ (pure counting); otherwise any partition passing Gates 1 and 3 could have $S_{\max}$ cross 2 by taking $\beta\varepsilon>0$ |
| 12 | **Functional form $F_D=\binom{D}{2}q^D$** | Taking the number of Lorentz generators times $q^D$ and choosing the physical dimension by the maximum is a readout functional hypothesis on par with $\pi$, not a theorem derived from $Z_0$ |

---

## 10. Discussion

$$
Z_0 \Longrightarrow \mathcal{W}_n \Longrightarrow \pi \text{ (self-identifying)} \Longrightarrow D=4.
$$

All results are conditional on quantumness; the axiom itself is unchanged.

The framework's free content consists of exactly four items: one axiom $Z_0$, one free parameter $n$, two named inputs $\pi$ and $\kappa$.

The two weakest links at present are: link ② (partition non-uniqueness, proven) and the quantumness premise of link ③ (unproven). Subsequent work targets these two links, as well as establishing the continuum limit and field-theoretic stage.

**Claims not made:** $Z_0$ has been modified; the chain holds unconditionally; $\pi$ has been uniquely determined; $D=4$ holds for all $n$; contextuality resolves the arrow of time; this work provides a continuum limit, metric, or gravitational field equations.

---

## References

[1] R. Penrose, *The Road to Reality: A Complete Guide to the Laws of the Universe*, Vintage Books, 2005.

[2] C. J. Isham, "Structural issues in quantum gravity," in *Recent Perspectives in the Quantum Theory of Gravity*, arXiv:gr-qc/9510034, 1994.

[3] A. Schrijver, *Theory of Linear and Integer Programming*, Wiley, 1986.

[4] W. Feller, *An Introduction to Probability Theory and Its Applications*, Vol. 1, 3rd ed., Wiley, 1968 (reflection principle and first-passage distributions).

[5] A. A. Klyachko, M. A. Can, S. Binicioğlu, A. S. Shumovsky, "A simple test for hidden variables in spin-1 systems," *Phys. Rev. Lett.* **101**, 020403 (2008).

[6] R. W. Spekkens, "Evidence for the epistemic view of quantum states: A toy theory," *Phys. Rev. A* **75**, 032110 (2007).

[7] Companion paper: *Conditional Reduction to General Relativity from Zero-Sum Constraints: Lovelock Uniqueness and the Group-Theoretic Origin of Dimension* (in preparation).
