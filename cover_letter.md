**To:** The Editors, *Foundations of Physics*

**Date:** October 15, 2026

**Re:** Submission of two interconnected manuscripts — *From Zero-Sum Constraints to Spacetime: A Discrete Substrate Yields Four-Dimensional Geometry and the Einstein Field Equations*

Dear Editors,

We are pleased to submit two interconnected manuscripts for consideration as a combined Research Article in *Foundations of Physics*. The two papers form a single derivation chain from a minimal ontological axiom to four-dimensional spacetime geometry and the Einstein field equations.

**Paper 1** establishes the dimension readout mechanism: starting from the single axiom that the equilibrium state is never static, we derive zero-sum constraints on compensatory moves, the space of balanced walks, and a self-identifying coarse-graining partition. We show that the dimension functional $F_D=\binom{D}{2}q^D$ peaks uniquely at $D=4$ when the participation index $q$ falls in $(1/2, 3/5)$; this condition is satisfied at $n=16$ on the height-level partition ($q=0.55777859$), with the KCBS contextuality gate passing on the top-three-block subspace.

**Paper 2** then invokes the Lovelock uniqueness theorem to compress the gravitational field equations into three structural properties—locality $(L)$, second-order $(O)$, and conserved source $(C')$—forcing them into Einstein form in four dimensions. The spatial metric is recovered via metric-volume consistency from a discrete Dirichlet energy, and the signature $(1,3)$ follows from reversible/irreversible splitting. Dimension is locked by two independent routes: the graviton-degree-of-freedom count, and the closure–motion correspondence that uniquely selects $K_4$ as the local fiber.

We explicitly register all conditional inputs and open items, including the unproven necessity of the quantum contextuality threshold and the postulational status of four-dimensional diffeomorphism covariance. No result is overclaimed.

**About the Lean 4 supplementary material.** The attached `ZeroLean_supplementary.zip` contains 16 Lean 4 proof files that machine-verify the key lemmas and numerical claims in both papers:

| Lean file | What it verifies |
|---|---|
| `Z001_Peak.lean` | Lemma 7.4: $F_D = \binom{D}{2} q^D$ has a unique finite maximum at $D=4$ iff $q \in (1/2, 3/5)$ |
| `Z011_Counting.lean` | Reflection-principle enumeration of balanced walks; block counts $b_h$ at $n=16$ |
| `Z017_Reseeding.lean` | Proposition 3.4: necessity of reseeding for persistent motion |
| `Z024_Conservation.lean` | Zero-sum identity $\sum_v d_w(v) \equiv 0$ for any walk $w$ |
| `Z025_Laplacian.lean` | Graph Laplacian positive semidefiniteness; kernel = constants |
| `Z044_Ledger.lean` / `Z079_Ledger.lean` | Dimension-selection functional $F_D$ and its spectral properties |
| `Z051_Partition.lean` | Self-identification criteria for the blocking partition |
| `Z065_CycleRank.lean` | Cycle-space rank $b_1(\Gamma) = |E| - |C| + 1$ |
| `Z086_Hyperplane.lean` | Zero-sum hyperplane $H_Q = \{x : \sum x_v = 0\}$; $\dim H_Q = |C| - 1$ |
| `Z089_Local.lean` | Locality property $(L)$ of the discrete energy functional |
| `Z090_K4MinusOne.lean` | $K_4$ group-theoretic structure; edge-transitivity; closure–motion correspondence |
| `Z109_Energy.lean` | Discrete Dirichlet energy positive definiteness on $H_Q$ |
| `Z124_Closure.lean` | Closure–motion correspondence principle; $\dim \ker B = \dim \operatorname{im} B$ condition |

These files can be compiled with Lean 4 (toolchain version specified in `lean-toolchain`) and provide independent machine-checkable verification of the central mathematical claims. They are included as supplementary material for reproducibility and transparency.

We confirm that this work is original, has not been published elsewhere, and is not under consideration by any other journal. All authors have approved the submission.

Thank you for your consideration. We look forward to your feedback.

Sincerely,

Guobin Ouyang
Independent Researcher
Foshan, China
