# Conditional Reduction to General Relativity from Zero-Sum Constraints: Lovelock Uniqueness and the Group-Theoretic Origin of Dimension

**Abstract**

This work establishes a conditional reduction chain from a discrete ontological substrate to general relativity. The method does not assemble the Einstein field equations piece by piece; instead, it compresses the field equations into three structural properties—locality $(L)$, second-order $(O)$, and conserved source $(C')$—and then invokes the Lovelock uniqueness theorem to force them into Einstein form. The algebraic skeleton: zero-sum constraints fix codimension 1; the incidence matrix of a connected graph spans the zero-sum hyperplane and gives integer reachability; discrete Dirichlet energy is positive semidefinite with kernel the constants; continuum approximation combined with metric-volume consistency recovers the spatial metric; reversible/irreversible splitting directly gives signature $(1,m-1)$. The metric construction: closed-walk counting gives a parameter-free edge weight $w_{ij}\propto\phi_i\phi_jA_{ij}$; finite truncation $k=n$ restores locality. The dimension account is the weakest link: the core postulates do not constrain dimension (proposition: for every $m\ge2$ there exists a model satisfying the postulates); $D=4$ follows from two independent routes—physically, "propagating gravitons exist" requires the unique positive solution of $\frac{D(D-3)}2=2$, and geometrically, the closure–motion correspondence principle $\dim\ker B=\dim\operatorname{im}B$ combined with edge-transitivity uniquely selects $K_4$. The premises of both routes are not consequences of the core postulates; dimension is therefore a conditional input.

**Keywords**: Lovelock theorem; general relativity; zero-sum constraint; discrete geometry; Hodge duality; dimension problem

---

## 1. Introduction

The spacetime dimension $D=4$ is an input rather than a conclusion within the standard physical framework [1,2]. Our approach starts from a discrete substrate containing no manifold, coordinates, metric, probability, or action, such that the shape of the Einstein field equations becomes a consequence derivable from the Lovelock uniqueness theorem.

**Strategy.** We do not assemble GR piece by piece. GR is compressed into three structural properties: $(L)$ locality, $(O)$ second-order, $(C')$ conserved source. The substrate need only provide these three; the shape of the geometric field equations is then locked by the Lovelock theorem [3].

**Scope.** The Lovelock theorem's literal statement is "in 4 dimensions." This lever therefore cannot be used to select $D=4$; dimension must be argued separately (§6).

---

## 2. Core Postulates of the Substrate

The substrate declares only the following; no manifold, coordinates, dimension, distance, curvature, action, or probability weights appear:

- **State**: codimension-1 affine constraint on integer lattice (zero-sum);
- **Transport**: connected, reversible compensatory moves;
- **Energy**: Dirichlet-type quadratic form;
- **Time**: step index $\tau$ and terminal;
- **Record**: multiplicity multiset of closed classes;
- **Weights**: integer multiplicities.

The substrate axiom is $Z_0=(\{\pm1\},\iota)+(\mathcal W_n,\mu)$, where $\iota$ is the direction involution and $\mu$ the uniform counting measure [4].

---

## 3. Algebraic Skeleton

### 3.1 Zero-Sum Codimension and Transport Rank

**Lemma 3.1 (Codimension 1).** The zero-sum constraint $\sum_v x_v=0$ gives only one linear condition:

$$
\dim H_Q=m-1, \qquad \{e_j-e_1\}_{j=2}^m \text{ is a basis}.
$$

**Corollary 3.2.** To obtain "four dimensions," one must first have $m-1=4$ or $m=4$; neither is a consequence of the substrate postulates.

**Lemma 3.3 (Transport rank and integer reachability).**

$$
\operatorname{im}B=H_Q, \qquad \operatorname{rank}B=m-1,
$$

$$
x\in H_Q\cap\mathbb Z^C \implies \exists f\in\mathbb Z^E,\; Bf=x.
$$

*Proof.* $B$ is totally unimodular [5]; integer reachability follows from standard integer linear programming results. $\blacksquare$

### 3.2 Discrete Dirichlet Energy

**Lemma 3.4.** The discrete Dirichlet energy with edge weights $K_e$

$$
E(f)=\sum_e K_e(f_i-f_j)^2 \ge 0, \qquad E(f)=0 \iff f \text{ is constant},
$$

i.e. $E=f^{\mathsf T}Lf$, $\ker L=\operatorname{span}\{\mathbf 1\}$, $L(H_Q)\subseteq H_Q$, $L|_{H_Q}\succ0$.

**Lemma 3.5 (No preference $\Rightarrow$ uniform conductance).** If $\Gamma$ is edge-transitive, then $K_e\equiv\kappa$.

*Remark.* This is conditional: path graphs $P_4$, $P_5$ have 2 edge orbits, and no-preference does not yield uniform conductance.

### 3.3 Local Quadratic Form and Spatial Metric

**Lemma 3.6 (Metric-volume consistency inversion).** Continuum approximation gives $E\to\int Q^{ij}\partial_if\,\partial_jf\,\omega$. Requiring metric-volume consistency

$$
\mu=\operatorname{vol}_h=\sqrt{\det h}\;\omega,
\tag{3.1}
$$

one obtains $Q=\sqrt{\det h}\,h^{-1}$, $\det Q=(\det h)^{d/2-1}$ ($d$ = spatial dimension). For $d\neq2$ the unique inverse is

$$
h=(\det Q)^{\frac{1}{d-2}}\,Q^{-1}.
\tag{3.2}
$$

For $d=2$, $\det Q=1$, and $h\mapsto\Omega^2h$ leaves $Q$ invariant; only the conformal class remains.

*Remark.* (3.1) is the sole additional input in this step; it binds the Dirichlet quadratic form one-to-one to the spatial metric and removes conformal ambiguity for $d\neq2$. The unique inverse holds for all $d\neq2$ and does not select any dimension.

### 3.4 Splitting, Time, and Signature

**Lemma 3.7 (Reversible/irreversible splitting).**

$$
T\cong H_Q\oplus\mathbb R_\tau\cong\mathbb R^m.
$$

The regular bilinear form can only be $g=d\tau^2-h$ or $g=-d\tau^2+h$. Taking the former gives $\operatorname{signature}(g)=(1,m-1)$.

*Mechanism.* The substrate's compensatory moves are reversible (each edge has an inverse move); step progression is irreversible. Spatial directions carry the inverse structure; the step direction does not. The signature originates here—there is no need to single out a time direction among the spatial directions.

*Boundary.* $1\oplus(m-1)$ and $h>0$ are derived; which side carries the minus sign is a sign convention.

### 3.5 Curvature and Two-Derivative Uniqueness

**Lemma 3.8.** Given a smooth manifold $M$ and nondegenerate $g$, the connection $\Gamma^a{}_{bc}$, curvature $R^a{}_{bcd}$, Ricci $R_{ab}$, and Einstein tensor $G_{ab}$ are uniquely determined by standard differential geometry [6].

**Lemma 3.9 (Uniqueness within the two-derivative family).** Let $E_{ab}=\alpha R_{ab}+\beta Rg_{ab}+\gamma g_{ab}$. Using $\nabla^ag_{ab}=0$ and the contracted Bianchi identity $\nabla^aR_{ab}=\tfrac12\nabla_bR$:

$$
\nabla^aE_{ab}=\Bigl(\frac{\alpha}{2}+\beta\Bigr)\nabla_bR.
\tag{3.3}
$$

Thus when $\nabla R\not\equiv0$, divergence-free $\iff\beta=-\alpha/2$, i.e. $E_{ab}\in\operatorname{span}\{g_{ab},G_{ab}\}$.

### 3.6 Lovelock Uniqueness and the GR Criterion

**Lemma 3.10 (Lovelock uniqueness) [3].** In 4 dimensions, a symmetric second-rank tensor constructed from $g$ and its derivatives up to second order, linear in the highest derivative, and satisfying $\nabla^aE_{ab}=0$ must be

$$
E_{ab}=c_1G_{ab}+c_2g_{ab}.
$$

(The Gauss–Bonnet term is a topological term in 4 dimensions and does not enter the field equations.)

**Corollary 3.11 (GR criterion).**

$$
(L)\wedge(O)\wedge(C') \implies G_{ab}+\Lambda g_{ab}=8\pi G\,T_{ab}.
\tag{3.4}
$$

Here $(C')$ requires that $T_{ab}$ be a functional of matter data (the metric only enters index raising/lowering and normalization) and must not contain independent metric degrees of freedom.

*Remark.* The Lovelock theorem's literal premise "in 4 dimensions" means it cannot be used to derive $D=4$; dimension must be argued separately (§6).

### 3.7 Locality Bifurcation and Excluded Routes

**Lemma 3.12.** The discrete Dirichlet energy is local; the effective resistance distance $R_{ij}=(e_i-e_j)^{\mathsf T}L^+(e_i-e_j)$ is nonlocal.

**Lemma 3.13 (Limitations of the algebraic assembly route).** The curvature of an equilateral complex is locked by combinatorial structure: the deficit angles around edge $k$ are $+148.41°,+77.88°,+7.36°,-63.17°$ ($k=3,4,5,6$), with only $k=3,4,5$ giving positive deficit. Algebraic assembly can therefore only land on an equilateral, rigid subset and cannot approximate a general metric.

| Excluded route | Basis |
|:--|:--|
| Algebraic assembly ($K\to h$) | Unique fixed point = uniform $\Rightarrow$ equilateralization $\Rightarrow$ one-parameter rigidity |
| Effective resistance as geometry | Nonlocal; external edges break local isotropy |
| Homogeneity degree for dimension selection | Homogeneity degree is 1 for every $d$; does not select dimension |

---

## 4. Metric Construction and the Lovelock Lever

### 4.1 Common Weight Formula

All three metric routes use the same weight formula

$$
w_{ij}=\phi_i\phi_jA_{ij},
\qquad h=\text{graph Laplacian with weight }w\text{ restricted to }H_Q,
\tag{4.1}
$$

differing only in the source of $\phi$.

**Route 1 (closed-walk counting, $k\to\infty$).** The total traversal count of edge $(i,j)$ by length-$k$ closed walks is

$$
N^{(k)}_{ij}=k\,A_{ij}(A^{k-1})_{ij}.
\tag{4.2}
$$

After normalization, $w_{ij}\propto\phi_i\phi_jA_{ij}$, with $\phi$ the Perron vector. Zero free parameters. But $w$ depends on the Perron eigenvector—a global quantity of the entire graph—so $(L')$ ($g$ is locally determined by data) fails.

**Route 2 (finite truncation $k=n$).** Truncating (4.2) at scale $k$:

$$
w^{(k)}_{ij}=\sum_{m=2}^{k} m\,A_{ij}(A^{m-1})_{ij}.
\tag{4.3}
$$

Locality holds exactly when $k\le$ graph diameter (the effect of distal perturbations on central edge weight ratios is exactly zero). $k>n$ counts unrealizable closed words; $k<n$ misses realizable ones, so $k=n$ is forced.

**Route 3 (locally given $\phi$).** If $\phi$ is given locally, $(L')$ holds exactly, and all three premises are satisfied.

### 4.2 Degrees of Freedom and Parameterization

The number of spatial metric components is

$$
\dim\operatorname{Sym}^2H_Q=\binom m2=\frac{(m-1)m}{2}.
\tag{4.4}
$$

*Key distinction.* (4.4) counts spatial metric $h$ components ($m=4$, 3D space: 6). Our construction $g=d\tau^2-h$ (ADM $N=1,\beta=0$) fixes $g_{00}=1$, $g_{0i}=0$; graph edge weights provide only the 6 components of $h$, never the 10 components of the 4D spacetime metric $g_{ab}$.

*Parameterization is not bijective.* The weight formula (4.1) uses a vertex field $\phi=(\phi_1,\dots,\phi_m)$ to give edge weights $w_{ij}=\phi_i\phi_jA_{ij}$. $\phi$ has only $m$ free parameters, while $\dim\operatorname{Sym}^2H_Q=\binom m2$; on $K_m$, the $\binom m2$ edge weights from $m$ vertex values always satisfy the cross-product constraint $w_{12}w_{34}=w_{13}w_{24}=w_{14}w_{23}=\phi_1\phi_2\phi_3\phi_4$. It therefore spans only a $\le m$-dimensional submanifold of $\operatorname{Sym}^2H_Q$, not a bijection. The injectivity of $w\mapsto h$ (for independent weights) holds, but the $\phi$ parameterization of (4.1) is a local parameterization on a submanifold, not proven to cover the full space.

*Graviton degrees of freedom.* In general 4D GR: 10 metric components $-4$ coordinate gauge $-4$ initial-value constraints $=2$ physical propagating degrees of freedom, i.e. $\frac{D(D-3)}2=2$.

### 4.3 Refinement Convergence and Nonlocality of the Limit

The edge-weight ratio converges with truncation $k$, but extremely slowly, and the limit is the Perron form $P_1=vv^{\mathsf T}$, globally nonlocal. The deviation decay length is

$$
\tau\approx\frac{2N^2}{3\pi^2}.
\tag{4.5}
$$

*Conclusion.* The effective metric has a limit under refinement; the limit is the Perron form, nonlocal $\Rightarrow$ there exists no UV-fixed local metric.

---

## 5. Continuum Limit and Causal Structure

### 5.1 $\Gamma$-Convergence

Lovelock requires $\Gamma$-convergence (local structure convergence), not a UV fixed point (self-similarity) [7]. A standard fact of homogenization theory: non-self-similar coefficient families still $\Gamma$-converge. With fixed step count $k=8$, the profile maximum difference between adjacent $N$ converges at approximately $O(a^2)$.

*Boundary.* $\Gamma$-convergence is not proven (no compactness, no $\Gamma$-liminf); only convergence and order are measured, and only on a 1D chain with a single profile.

### 5.2 Causal Structure: Support Cone and Cone Edge

Finite causal velocity requires saturation (per-point capacity), not current relaxation. Out-of-cone signals are exponentially small but nonzero pulled fronts [8].

**The support cone is always exact**: nearest-neighbor transport expands at most one lattice site per step, $\operatorname{supp}\rho_n\subseteq[-n,n]$, independent of $B$.

**Cone-edge visibility** is derived analytically: on the $e^{-\mu x}$ mode, the one-step mapping characteristic equation is

$$
\Lambda^2=B\cosh^2\mu \implies \lambda(\mu)=\tfrac12\log B+\log\cosh\mu.
\tag{5.1}
$$

The cone-edge density decay rate is

$$
\varepsilon(B)=\tfrac12\log\Bigl(\tfrac4B\Bigr).
\tag{5.2}
$$

*Visible cone edge $\iff B=4$.* $B=4$ is upgraded from a model parameter to a self-consistency condition.

---

## 6. Dimension Account

**This is the weakest link of the entire chain.**

### 6.1 The Core Postulates Do Not Constrain Dimension (No-Go)

**Proposition 6.1 (Classification-level no-go).** Let $\mathcal A$ be the substrate core postulates. Then for every integer $m\ge2$, there exists a model $\mathcal M_m$ of $\mathcal A$ with channel-set cardinality $|C|=m$.

*Proof.* Take $C=\{0,\dots,m-1\}$, $E=\{\{i,i+1\}:i\in\mathbb Z_m\}$ (cycle graph $C_m$). Verify postulate by postulate: $C_m$ is finite and connected; every move has an inverse; all $2m$ directed moves are taken with equal weight; step index progresses irreversibly; closed state $x=0$ exists. Numerical verification ($m=2,\dots,12$): $\operatorname{rank}B=m-1$, $\ker L=\operatorname{span}\{\mathbf1\}$, $L|_{H_Q}\succ0$. $C_m$ is edge-transitive, so the no-preference $\Rightarrow$ uniform conductance premise also holds. $\blacksquare$

$$
\implies \text{any statement depending only on }D\text{ cannot be derived from the substrate core postulates}.
\tag{6.1}
$$

### 6.2 Excluding Low Dimensions

| $D$ | Conclusion | Basis |
|:--:|:--|:--|
| $2$ | No gravitational dynamics | $R_{ab}=\tfrac12Rg_{ab}\implies G_{ab}\equiv0$ |
| $3$ | No graviton | Weyl tensor vanishes identically; curvature fully determined by local matter |
| $4$ | Propagating graviton exists | see §6.3 |

### 6.3 Uniqueness of $D=4$ from Graviton Degrees of Freedom

The $D$-dimensional graviton degrees of freedom are $\frac{D(D-3)}2$:

| $D$ | 2 | 3 | 4 | 5 | 6 | 7 |
|:--|--:|--:|--:|--:|--:|--:|
| $\frac{D(D-3)}2$ | $-1$ | $0$ | $\mathbf2$ | $5$ | $9$ | $14$ |

$$
\frac{D(D-3)}2=2 \iff D^2-3D-4=0 \iff D\in\{-1,4\} \implies D=4.
\tag{6.2}
$$

*Screening chain.* $\{2,3,\dots\}\xrightarrow{\ D(D-3)/2>0\ }\{4,5,\dots\}\xrightarrow{\text{minimality}}4$.

*Boundary.* Minimality, even-dimension reading, and sign–polarization identification are not consequences of the substrate postulates.

### 6.4 Candidate Constraints as Synonymous Reformulations

Of 15 candidate constraints tested item by item, items 1–10 all have survival set $\{4\}$; they are all synonymous reformulations of $(Z_2)'$ and do not constitute independent evidence.

Among them, item 7, $b_1(K_m)=m-1$, is examined in detail below.

### 6.5 Representation-Theoretic Status of the Cycle–Displacement Duality

As an **oriented** cycle space, $Z_1(K_4)$ decomposes as the standard $S_4$-module $(2,1,1)\cong(3,1)\otimes\mathrm{sgn}$, which is the representation of axial vectors (pseudovectors, flux, magnetic field). The displacement space is $H_Q=(3,1)$, the standard polar-vector representation. The two are distinct $S_4$-modules (they differ by the sign representation). **However, in three-dimensional geometric realization** ($K_4$ embedded as a regular tetrahedron in $\mathbb R^3$), the Hodge star operator

$$
\star:\Lambda^2\mathbb R^3\longrightarrow\mathbb R^3
$$

establishes a canonical duality between the axial cycle space and the polar displacement space. This is not a numerical coincidence: it is the group-theoretic counterpart of §6.7(ii) "antisymmetric composition closure"—the identity $\Lambda^2\mathbb R^3\cong\mathbb R^3\otimes\mathrm{sgn}$ is precisely the Hodge duality between closed 2-forms (loops, flux) and 1-forms (displacements) in three dimensions [9].

### 6.6 Dimension-Blind No-Go

All four metric-side structures $u=h/\sqrt n$, $\mathrm H(w)$, $\min(k,n-k)$, $T(M,n)$ do not reference the adjacency matrix $A$ and take the same value on any graph.

$$
\implies \textbf{the metric side cannot exclude low dimensions (e.g., }d=1\text{)}.
\tag{6.3}
$$

### 6.7 The Closure–Motion Correspondence Principle

**The core postulates do not constrain dimension** (Proposition 6.1). What fixes the dimension is this correspondence:

> **Closure–motion correspondence [named input].** The axiom, read at the counting level: $\dim\ker B=\dim\operatorname{im}B$ (degrees of freedom of closed loops = degrees of freedom of motion).

Here $B$ is the incidence matrix of compensatory moves: $\operatorname{im}B=H_Q$ is the zero-sum state, and $\ker B$ is the closed-loop space.

$$
\begin{aligned}
&\text{dictionary}: & d=m-1 \quad\bigl(\text{uniquely fixed by }\binom m2-(m-1)=\binom{m-1}2\bigr);\\
&\text{reading}: & \underbrace{\binom m2-(m-1)}_{\dim\ker B}=\underbrace{m-1}_{\dim\operatorname{im}B=d}\iff m=4\ (\text{unique});\\
&\text{no preference (edge-transitive)}: & \deg\text{ equal}\implies t\,m=4(m-1)\iff(\deg,m)=(3,4)\ (\text{unique nontrivial solution}).
\end{aligned}
\tag{6.4}
$$

**Machine check**: on $K_m$, $\dim\ker B=\binom m2-(m-1)$ and $\dim\operatorname{im}B=m-1$ are equal only when $m=4$ ($m=2$ gives $0\ne1$, $m=3$ gives $1\ne2$, $m\ge5$ always gives $\dim\ker B>\dim\operatorname{im}B$).

Together the three give: **spatial dimension $d=3$, channel count $m=4$, $\deg=3$, spacetime dimension $D_{\rm st}=d+1=4$**.

$$
D_{\rm st}=4\text{ follows from two independent routes: the substrate-side premise is the closure--motion correspondence, the physical-side premise is "propagating gravitons exist".}
\tag{6.5}
$$

*Why the correspondence cannot be a consequence of the postulates.* The substrate gives $\Gamma=(C,E)$ together with "motion"; the postulates only constrain *how* to move—"full branching" means every *allowed* motion is instantiated, and does not constrain which pairs have allowed moves. Any condition of the form "$\Gamma$ must be complete/self-dual" is a choice of $\Gamma$, not a consequence of the postulates. Counterexample: cycle graph $C_m$ ($m\ge4$) satisfies the postulates and edge-transitivity, but $\dim\ker B=1\ne m-1$.

### 6.8 Equivalent Formulations of the Correspondence Principle

The counting content $\binom d2=d$ has three equivalent formulations:

**(i) Pairing uniqueness.** Each edge has exactly one opposite edge $\iff\binom{m-2}2=1\iff m=4$.

**(ii) Antisymmetric composition closure.** $\dim\Lambda^2V=\dim V\iff D=3$: the "relative relation" of two directions is again a direction (cross-product structure). I.e. $\binom D2=D\iff D\in\{0,3\}$, with $D=3$ the unique nontrivial solution.

**(iii) $A_3$ self-duality.** The 6 edge directions $\pm(e_i-e_j)$ of $K_4$ are exactly the 6 positive roots of $A_3$; $\operatorname{Aut}(K_4)=S_4$ is the Weyl group of $A_3$. Among $A_n$, only $A_3$ satisfies "number of relations = dimension": $\binom{n+1}2-n=n\iff n=3$. This is precisely the accidental isomorphism $A_3\cong D_3$ ($\mathfrak{so}(3)\cong\mathbb R^3$ cross-product isomorphism).

$$
\implies \Gamma=K_4:\ 4\text{ vertices},\ \deg=3,\ 6\text{ edges},\ 3\text{ opposite-edge pairs};\ \text{angle}=\arccos(-1/3)=109.4712°;\ \operatorname{Aut}(K_4)=S_4.
\tag{6.6}
$$

---

## 7. Axiomatic Budget

### 7.1 What the Lovelock Lever Presupposes

| # | Presupposition | Nature |
|--:|:--|:--|
| 1 | $D=4$ (literal hard-coding) | Cannot select dimension via this lever |
| 2 | Metric exists and is smooth | Input |
| 3 | Metric-volume consistency $\mu=\operatorname{vol}_h$ | Sole additional input of (3.1) |
| 4 | Second-order $(O)$ | Excludes higher-order terms |
| 5 | Conserved source $(C')$ | $(C)$ as-is is vacuous |
| 6 | Source = relativistic matter | Additional claim |
| 7 | Metric nondegenerate and Lorentzian | Input |
| 8 | 4D diffeomorphism covariance $\nabla^aE_{ab}=0$ | Input: the preferred foliation of $g=d\tau^2-h$ gives only 3D spatial diffeomorphism invariance |

### 7.2 Absence of Intrinsic Scales

$$
\text{The dimensionful degrees of freedom in the geometric sector are exactly }1.
\tag{7.1}
$$

*Evidence.* The layer coefficient (edge $(0,1)$) $\binom{m-1}{m/2}$ is an integer ($m=2,4,6,8$ gives $1,3,10,35$). All raw metric-side data are integers $\Rightarrow$ they cannot carry length.

*Unit-transformation invariance.* Under $W\to sW$ ($s$ spanning 8 orders of magnitude), normalized eigenvalue deviation is $1.89\times10^{-15}$, edge-weight ratio deviation is $0$.

*On $G,\Lambda.* Dimensional constants are units, not parameters [10]. "Deriving the numerical values of $G,\Lambda,\hbar$" is therefore not a well-posed problem.

---

## 8. Main Results and Boundaries

### 8.1 Main Result

$$
\begin{aligned}
&\text{substrate (zero-sum constraint) + (inputs: embedding / source class / scale)}\\
&\quad\Longrightarrow (L)\wedge(O)\wedge(C')\\
&\quad\Longrightarrow G_{ab}+\Lambda g_{ab}=8\pi G\,T_{ab}\quad(\text{in }D=4).
\end{aligned}
\tag{8.1}
$$

$$
\text{the shape of the field equations is uniquely forced; the metric and dimension are inputs; the source is a no-go.}
\tag{8.2}
$$

### 8.2 Layering of "Dimension Is an Input"

"Input" here refers specifically to the Lovelock layer: the shape theorem $E_{ab}=\alpha G_{ab}+\beta g_{ab}$ takes $D$ as an argument, and Lovelock itself does not select dimension. But the value $D=4$ is not a bare free input; it is given by the readout chain: peak condition locks $L=4$ (Gate 1) + closure-equality $\dim\ker B=\dim\operatorname{im}B$ locks $K_4/A_3$. This readout chain is half-locked: the peak condition hard-presses to $L=4$; the step "$L=4$ is exactly the physical spacetime dimension" is only count-level equality—$Z_1(K_4)\cong(2,1,1)\cong(3,1)\otimes\mathrm{sgn}$ and $H_Q=(3,1)$ are distinct $S_4$-modules differing by a parity sign, but in three-dimensional geometry form an axial–polar vector duality via the Hodge star.

### 8.3 Boundaries That Must Be Stated

1. **Division of labor of Lovelock**: it only answers "what must the field equations look like," not "why does geometry have dynamics."
2. **Local $K_4$ fiber vs. macroscopic grid are two structures.** $\dim H_Q=m-1=3$ holds only for the *single-point* local $K_4$ fiber; on the macroscopic $N$-vertex grid, $\dim H_Q=N-1\ne3$. The construction bridge "how local $K_4$ fibers/cells tile a macroscopic 3D manifold grid" is missing; it is an open item.
3. **$L=4$ uses Gate 1 only.** The peak condition locks $L=4$; this $k=2$ partition does not pass the contextuality gate or the density gate, so it serves only as dimension readout, not as a quantumness criterion.

---

## 9. Open Problem

**Continuum limit ($\Gamma$-convergence)** is unproven: no compactness, no $\Gamma$-liminf; only on a 1D chain with a single profile. This is analytic (compactness / variational limit) and requires independent work in homogenization theory.

---

## 10. Conclusion

$$
\begin{aligned}
&\text{substrate (zero-sum constraint)}\\
&\quad\Longrightarrow \dim H_Q=m-1;\ \operatorname{im}B=H_Q;\ \text{discrete Dirichlet energy};\\
&\quad\Longrightarrow\text{(continuum + metric-volume consistency) spatial metric }h;\\
&\quad\Longrightarrow\text{(reversible/irreversible splitting) signature }(1,m-1);\\
&\quad\Longrightarrow\text{(locality}\wedge\text{second-order}\wedge\text{conserved source)}\\
&\quad\Longrightarrow G_{ab}+\Lambda g_{ab}=8\pi G\,T_{ab}\quad(\text{in }D=4).
\end{aligned}
\tag{10.1}
$$

**Shape is unconditional; entry is conditional**: $(L)\wedge(O)\wedge(C')\Rightarrow$ Einstein form is unconditional (Lovelock); dimension, metric, and source are three entries, all non-consequences of the postulates.

Current progress is **conditional recovery**, not **unconditional derivation**. To reach *derivation*, the same ledgers are missing: the premise of dimension, the origin of the metric, and the price tag of the source.

---

## References

[1] R. Penrose, *The Road to Reality*, Vintage Books, 2005.

[2] C. J. Isham, "Structural issues in quantum gravity," arXiv:gr-qc/9510034, 1995.

[3] D. Lovelock, "The Einstein tensor and its derivatives," *J. Math. Phys.* **12**, 498 (1971).

[4] W. Feller, *An Introduction to Probability Theory and Its Applications*, Vol. 1, 3rd ed., Wiley, 1968.

[5] A. Schrijver, *Theory of Linear and Integer Programming*, Wiley, 1986.

[6] M. P. do Carmo, *Riemannian Geometry*, Birkhäuser, 1992.

[7] D. G. Cioranescu, P. Donato, *An Introduction to Homogenization*, Oxford University Press, 1999.

[8] J. D. Murray, *Mathematical Biology*, Vol. II: Spatial Models and Biomedical Applications, 3rd ed., Springer, 2003.

[9] H. Weyl, *Classical Groups: Their Invariants and Representations*, Princeton University Press, 1939.

[10] M. J. Duff, L. B. Okun, G. Veneziano, "Trialogue on the number of fundamental constants," *JHEP* **0203**, 023 (2002), arXiv:physics/0110060.
