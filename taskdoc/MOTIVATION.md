# QuantumDFT — Motivation

> Frozen clean-room specification `quantumdft-v0.2-qe7.5-2026-07-21`.
> Read together with `CONTRACT.md` (the complete V0 public boundary) and
> `TESTS.md` (visible unit and integration anchors). This document defines
> scientific intent and invariants, not implementation technique.

## Why this package exists

Quantum ESPRESSO (QE) is a mature reference implementation for plane-wave
pseudopotential density-functional theory, but its public workflow is centred
on executable-specific namelist files and a Fortran implementation that is not
designed as a differentiable Julia library. QuantumDFT reconstructs one narrow,
scientifically complete vertical slice in Julia under the MIT licence:

1. read a documented subset of a QE 7.5 `pw.x` input;
2. solve the periodic Kohn-Sham ground state;
3. expose energy, forces, stress, density, eigenvalues, and occupations;
4. compute implicit first-order electronic response at a commensurate phonon
   wave vector `q`;
5. construct the dynamical matrix and phonon modes; and
6. compute dielectric and Born effective-charge tensors for the non-analytic
   correction of polar crystals.

QE 7.5 is a black-box numerical oracle and an external compatibility target.
It is not an implementation dependency. The package is implemented from this
specification, published equations, documented formats, and independently
generated observations; QE source and tests are outside the clean-room source
boundary.

## Scientific model

V0 treats three-dimensional periodic, spin-unpolarized, time-reversal-
symmetric insulating crystals using scalar-relativistic norm-conserving UPF
2.0.1 pseudopotentials and either LDA or PBE exchange-correlation.

For nuclear coordinates `R`, cell `h`, and electronic variables `x`, the
ground state is a stationary solution of a discretized Kohn-Sham problem,
written abstractly as `residual(x, R, h) = 0`. Total energy is reported in
Hartree. Forces and stress are derivatives of the same converged discrete
energy, not separately fitted observables.

An atomic perturbation with reduced reciprocal coordinate `q` couples Bloch
states at `k` and `k+q`. V0 response differentiates the converged stationary
problem implicitly; differentiating a finite sequence of SCF iterations is not
part of the contract. The phonon dynamical matrix is the mass-weighted second
derivative of the Born-Oppenheimer energy. Its eigenvalues are signed squared
frequencies, and its eigenvectors are mass-normalized Cartesian displacement
patterns.

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

## V0 invariants

- Electron number: the integrated density equals the model electron count.
- Orthonormality: Bloch orbitals are orthonormal at every electronic k point.
- Variational consistency: energy, forces, stress, and response refer to the
  same converged basis and pseudopotential interpretation.
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
compatibility, or byte-for-byte reproduction of QE terminal output.

These omissions are outside the domain rather than silently approximated.
Unsupported input must fail before numerical work begins with an informative
exception naming the unsupported field or feature.

## Completion criterion

V0 is complete only when every public symbol in `CONTRACT.md` has direct public
tests, every contract clause is represented in the private verification
matrix, and the Si, NaCl, and low-symmetry held-out crystal workflows in `TESTS.md` pass
against the pinned QE 7.5 oracle at fixed pseudopotential, cutoff, k mesh, q
point, and convergence thresholds.
