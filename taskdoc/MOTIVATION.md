# EspressoDFT — Motivation

> Frozen clean-room specification `espressodft-v0.3-qe7.5-2026-07-21`.
> Read together with `CONTRACT.md` (the complete V0 public boundary) and
> `TESTS.md` (visible unit and integration anchors). This document defines
> scientific intent and invariants, not implementation technique.

## Why this package exists

Quantum ESPRESSO (QE) is a mature reference implementation for plane-wave
pseudopotential density-functional theory, but its public workflow is centred
on executable-specific namelist files and a Fortran implementation that is not
designed as a differentiable Julia library. EspressoDFT reconstructs one narrow,
scientifically complete vertical slice in Julia under the MIT licence. Its
defining distinction is that the converged Kohn-Sham problem is a differentiable
implicit layer from the first implementation milestone, rather than an SCF
program to which derivatives are attached later.

V0 therefore provides one canonical path that can:

1. read a documented subset of a QE 7.5 `pw.x` input;
2. solve the periodic Kohn-Sham ground state;
3. expose energy, forces, stress, density, eigenvalues, and occupations;
4. differentiate converged energy and density with respect to continuous atomic
   and strain perturbations through backend-neutral Julia AD rules;
5. compute implicit first-order electronic response at a commensurate phonon
   wave vector `q`;
6. construct the dynamical matrix and phonon modes; and
7. compute dielectric and Born effective-charge tensors for the non-analytic
   correction of polar crystals.

QE 7.5 is a black-box numerical oracle and an external compatibility target.
It is not an implementation dependency. The package is implemented from this
specification, published equations, documented formats, and independently
generated observations; QE source and tests are outside the clean-room source
boundary.

EspressoDFT is an independent project and is neither affiliated with nor
endorsed by the Quantum ESPRESSO Foundation. References to Quantum ESPRESSO
identify compatibility semantics and numerical provenance, not project
ownership or official status.

## Scientific model

V0 treats three-dimensional periodic, spin-unpolarized, time-reversal-
symmetric insulating crystals using scalar-relativistic norm-conserving UPF
2.0.1 pseudopotentials and either LDA or PBE exchange-correlation.

For continuous parameters `theta` (including nuclear coordinates `R` and cell
`h`) and electronic variables `z`, the ground state is a constrained stationary
solution of one discretized Kohn-Sham problem,

```text
F(z*(theta), theta) = 0.
```

The primal calculation solves this equation. Its forward derivative solves
`F_z delta_z = -F_theta delta_theta`; its reverse derivative solves the adjoint
linear system involving `F_z'`. These are respectively direct and adjoint forms
of density-functional perturbation theory. Total energy is reported in Hartree.
Forces, stress, response density, force constants, Born charges, and dielectric
response are derivatives of this same stationary problem, not separately fitted
observables or unrelated numerical paths.

An atomic perturbation with reduced reciprocal coordinate `q` couples Bloch
states at `k` and `k+q`. V0 response differentiates the converged stationary
problem implicitly. SCF mixing, iterative diagonalization, and stopping history
are primal algorithms and are not part of the differentiated mathematical map.
The phonon dynamical matrix is the mass-weighted second derivative of the
Born-Oppenheimer energy. Its eigenvalues are signed squared frequencies, and its
eigenvectors are mass-normalized Cartesian displacement patterns.

Differentiability is semantic rather than instruction-by-instruction. FFTs,
local and nonlocal pseudopotential actions, Hartree and exchange-correlation
kernels, orthogonal projectors, and linear solves may use optimized mutable
workspaces, but their public composition has mathematically defined JVP and VJP
rules. Degenerate eigenvectors are never assigned individually meaningful
derivatives; the occupied subspace projector is the gauge-invariant object.

For polar insulators, the electronic dielectric tensor and Born effective
charges provide the direction-dependent long-range non-analytic contribution
near Gamma. The three translational acoustic modes at Gamma must vanish in the
converged and complete-basis limit.

## Why the API has two boundaries

The QE-compatibility boundary preserves documented input semantics so existing
workflows can enter the package without first being rewritten. The native
Julia boundary exposes typed models, bases, converged states, and perturbations
so differentiation and new algorithms do not inherit QE namelist structure.

Both boundaries translate to one canonical model. A calculation created from
a QE input and an equivalent native construction must therefore produce the
same discretized problem and agree within the tolerances in `CONTRACT.md`.
Differentiability is behaviour of this canonical model, not a third input API
and not a second numerical backend.

## V0 invariants

- Electron number: the integrated density equals the model electron count.
- Orthonormality: Bloch orbitals are orthonormal at every electronic k point.
- Variational consistency: energy, forces, stress, and response refer to the
  same converged basis and pseudopotential interpretation.
- Differential consistency: AD gradients of converged energy reproduce public
  forces and stress, while the derivative of density reproduces direct response.
- Adjoint duality: direct and reverse sensitivities obey the real inner-product
  identity to their declared numerical tolerance.
- Solver-history independence: derivatives depend on the converged stationary
  problem, not the number or kind of iterations used to reach it.
- Hermiticity: the Hamiltonian and phonon dynamical matrix are Hermitian up to
  their declared numerical tolerances.
- Periodicity: shifting a reduced atomic coordinate by an integer lattice
  vector does not change observable results.
- Permutation covariance: consistently permuting atoms permutes atom-indexed
  results but leaves scalar observables and phonon spectra unchanged.
- Gauge invariance: unitary rotation within a degenerate occupied subspace does
  not change density, forces, response density, or the dynamical matrix.
- Acoustic sum rule: uniform translation has zero restoring force at Gamma in
  the converged limit.
- q conjugacy: `D(-q)` is the complex conjugate of `D(q)`.
- Native/QE equivalence: equivalent native and parsed-QE inputs produce the
  same canonical model and scientific outputs.

## Explicit non-goals

V0 does not promise ultrasoft or PAW pseudopotentials, non-collinear spin,
spin-orbit coupling, metals or smearing, DFT+U, exact exchange, van der Waals
corrections, molecular dynamics, geometry optimization, NEB, electron-phonon
transport, anharmonic force constants, MPI, GPU execution, binary QE restart
compatibility, byte-for-byte reproduction of QE terminal output, differentiation
with respect to discrete species/k-mesh/band-count choices, or arbitrary third-
and-higher derivatives. V0 guarantees the first derivatives and selected second
or mixed derivatives named in `CONTRACT.md`; the same stationary-response
architecture is intended to support higher-order extensions without making
them part of this frozen task.

These omissions are outside the domain rather than silently approximated.
Unsupported input must fail before numerical work begins with an informative
exception naming the unsupported field or feature.

## Completion criterion

V0 is complete only when every public symbol in `CONTRACT.md` has direct public
tests, every contract clause is represented in the private verification matrix,
the direct/adjoint differential-consistency tests pass, and the Si, NaCl, and
low-symmetry held-out crystal workflows in `TESTS.md` pass against the pinned QE
7.5 oracle at fixed pseudopotential, cutoff, k mesh, q point, and convergence
thresholds.
