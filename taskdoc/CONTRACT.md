# EspressoDFT — Contract (frozen differentiable V0)

> Frozen clean-room specification `espressodft-v0.3-qe7.5-2026-07-21`.
> Everything in this file is immutable for the V0 implementation task. A
> change requires a new specification identifier and regeneration of both
> visible and private tests.

## Deliverable and licence

Deliver a standard Julia package named `EspressoDFT`, compatible with Julia
1.12 or later and licensed under MIT. It exports exactly the symbols in the
Public API table below. Internal modules, algorithms, storage, and solver
choices are unrestricted provided the observable contract is met.

The package must provide backend-neutral `ChainRulesCore.frule` and
`ChainRulesCore.rrule` semantics for the differentiable public compositions
defined below. `ChainRulesCore` may be a package dependency; a particular AD
frontend such as Zygote, Enzyme, or ForwardDiff must not be required to use the
ordinary primal API and is not an additional public export.

The package is an independent clean-room implementation, not an official
Quantum ESPRESSO distribution or a project endorsed by its Foundation.

The clean-room implementation boundary consists of these three documents,
published mathematical references, the QE 7.5 user-facing input
documentation, the UPF 2.0.1 format documentation, and black-box numerical
observations. QE source code and QE test implementations are not permitted
implementation inputs.

## V0 domain

The valid domain is:

- three-dimensional periodic crystals;
- CPU execution with `Float64` arithmetic;
- spin-unpolarized, time-reversal-symmetric insulators at zero electronic
  temperature and fixed integer occupation;
- scalar-relativistic norm-conserving UPF 2.0.1 pseudopotentials;
- LDA or PBE exchange-correlation;
- finite plane-wave cutoff and Monkhorst-Pack electronic k mesh;
- phonon q vectors in reduced reciprocal coordinates for which `k+q` is a
  permutation of the full electronic k mesh.

Within this domain, the continuous differential variables in V0 are Cartesian
atomic positions and homogeneous cell strain with fractional positions held
fixed. Atomic species, pseudopotential identity, functional choice, charge,
spin, plane-wave cutoff, reciprocal/FFT index sets, k mesh, band count,
occupation pattern, q-mesh topology, and solver options are nondifferentiable
configuration data.

Inputs outside this domain must be rejected explicitly; silently replacing an
unsupported feature by a supported approximation is a contract failure.

## Units and conventions

- Native lattice vectors are columns of a `3×3` matrix in bohr.
- Native positions are a `3×N` matrix and are fractional by default.
- Native atomic masses are a length-`N` vector in electron masses; the QE
  parser converts `ATOMIC_SPECIES` masses from unified atomic mass units.
- Energy and eigenvalues are in Hartree.
- Forces are Cartesian and in Hartree per bohr.
- Stress is the symmetric tensor `sigma = (1/Omega) dE/depsilon`, in Hartree
  per bohr cubed; positive diagonal entries are tensile.
- Density values are electrons per bohr cubed.
- Electronic k and phonon q coordinates are reduced reciprocal coordinates.
- Atom indices and Cartesian directions are one-based.
- The dynamical matrix uses Cartesian atomic coordinates and
  `D[I,a,J,b] = Phi[I,a,J,b] / sqrt(M[I] M[J])`, with masses in electron masses.
- At Gamma, `response(gs, AtomicDisplacement(I,a,(0,0,0))).delta_density`
  is the derivative of the periodic density grid with respect to a unit
  Cartesian displacement of atom `I` in direction `a`, measured per bohr.
- Phonon frequencies are signed square roots of dynamical-matrix eigenvalues,
  in atomic angular-frequency units; negative values denote imaginary modes.

The QE parser converts QE's documented Ry/Bohr/alat conventions into these
native conventions exactly once at the compatibility boundary.

## Comparison semantics

- `exact` means matching type category, shape, keys, ordering, and exact value.
- `num(atol,rtol)` means elementwise `abs(got-ref) <= atol + rtol*abs(ref)`;
  NaN is never accepted for a valid finite reference.
- `energy` means `5e-7` Hartree per atom absolute plus `5e-8` relative.
- `force` means `5e-6` Hartree per bohr absolute plus `5e-6` relative.
- `stress` means `5e-6` Hartree per bohr cubed absolute plus `5e-6` relative.
- `band` means `5e-6` Hartree absolute plus `5e-6` relative after sorting only
  within numerically degenerate groups.
- `density` means integrated electron number within `2e-7` and selected
  reciprocal-density coefficients within `2e-5` relative/absolute.
- `response-density` means selected complex response-density coefficients
  within `2e-5` relative/absolute.
- `dynamical` means dynamical-matrix entries within `2e-9` absolute plus
  `2e-4` relative. The absolute floor is deliberately below the `10^-6`
  atomic-unit scale of the frozen fixtures, so an all-zero matrix cannot pass.
- `phonon` means matching signed frequencies within `2 cm^-1` after sorting
  within degenerate groups; eigenvectors are compared by subspace projectors.
- `polar` means Born effective charges within `5e-4` and dielectric-tensor
  entries within `5e-3` relative/absolute.
- `ad-gradient` means an AD directional derivative agrees with the corresponding
  public force, stress, or density response within `5e-5` absolute plus `5e-5`
  relative after accounting for units and sign.
- `ad-duality` means direct and adjoint scalar contractions agree within `2e-6`
  absolute plus `2e-5` relative.
- `ad-second` means a derivative of force agrees with the mass-unweighted
  Gamma force-constant entry within `2e-7` absolute plus `5e-4` relative.

Numerical oracle comparisons use the same structure, pseudopotential bytes,
cutoff, k mesh, q, functional, and convergence thresholds. Bitwise equality
with QE is not required.

## Public API

| ID | exported symbol | frozen signature | result/comparison | valid domain |
|---|---|---|---|---|
| `API-001` | `Crystal` | `Crystal(lattice::AbstractMatrix, species::AbstractVector{Symbol}, positions::AbstractMatrix; masses::AbstractVector, positions_are_fractional::Bool=true)` | immutable crystal semantics; `exact` | lattice is finite, nonsingular `3×3`; positions are finite `3×N`; species and positive finite masses both have length `N` |
| `API-002` | `KSModel` | `KSModel(crystal::Crystal; pseudopotentials::AbstractDict, xc::Symbol=:pbe, charge::Real=0, spin::Symbol=:unpolarized)` | opaque model | one NC-UPF path for each species; `xc in (:lda,:pbe)`; `charge==0`; `spin==:unpolarized` |
| `API-003` | `PlaneWaveBasis` | `PlaneWaveBasis(model::KSModel; Ecut::Real, kgrid::NTuple{3,<:Integer}, fft_size::Union{Nothing,NTuple{3,<:Integer}}=nothing)` | opaque basis | `Ecut>0`; positive k-grid dimensions; optional FFT grid represents every density component required by the wavefunction cutoff |
| `API-004` | `SCFOptions` | `SCFOptions(; energy_tolerance::Real=1e-10, density_tolerance::Real=1e-8, maxiter::Integer=100, extra_bands::Integer=4)` | immutable options; `exact` | positive tolerances and iteration count; nonnegative extra bands |
| `API-005` | `QEInput` | type returned by `read_qe_input` | opaque canonical run specification | scoped QE input described below |
| `API-006` | `AtomicDisplacement` | `AtomicDisplacement(atom::Integer, direction::Integer, q::NTuple{3,<:Real})` | immutable perturbation; `exact` | valid atom; direction in `1:3`; finite q commensurate with the basis k mesh |
| `API-007` | `read_qe_input` | `read_qe_input(path::AbstractString)` and `read_qe_input(io::IO)` | `QEInput`; semantic equivalence | scoped QE 7.5 `pw.x` input |
| `API-008` | `run_qe_input` | `run_qe_input(input::QEInput)` and path/IO convenience methods | same opaque ground-state result as `ground_state` | parsed calculation is `scf` and otherwise inside V0 |
| `API-009` | `ground_state` | `ground_state(basis::PlaneWaveBasis; options::SCFOptions=SCFOptions())` | opaque converged result; differentiable implicit layer | valid V0 basis; convergence achieved |
| `API-010` | `energy` | `energy(gs)` | scalar; `energy`, `ad-gradient` | converged ground state |
| `API-011` | `forces` | `forces(gs)` | real `3×N`; `force`, `ad-second` | converged ground state |
| `API-012` | `stress` | `stress(gs)` | symmetric real `3×3`; `stress`, `ad-gradient` | converged ground state |
| `API-013` | `density` | `density(gs)` | named tuple `(values, cell_volume)`; `density`, `ad-duality` | `values` is a real three-dimensional periodic grid |
| `API-014` | `eigenvalues` | `eigenvalues(gs)` | vector per full k point; `band` | same k-point order as canonical full mesh |
| `API-015` | `occupations` | `occupations(gs)` | vector per full k point; `num(1e-12,1e-12)` | occupations include spin degeneracy and integrate to electron number with k weights |
| `API-016` | `response` | `response(gs, perturbation::AtomicDisplacement; tolerance::Real=1e-8, maxiter::Integer=200)` | named tuple `(delta_density, residual_norm, converged)`; `response-density` | positive tolerance and `maxiter`; perturbation in V0 response domain |
| `API-017` | `dynamical_matrix` | `dynamical_matrix(gs, q::NTuple{3,<:Real}; tolerance::Real=1e-8, maxiter::Integer=200)` | complex `3N×3N`; `dynamical` | finite commensurate q; positive `maxiter`; V0 response assumptions |
| `API-018` | `phonon_modes` | `phonon_modes(gs, q::NTuple{3,<:Real}; tolerance::Real=1e-8, maxiter::Integer=200)` | named tuple `(frequencies, eigenvectors)`; `phonon` | same as `dynamical_matrix` |
| `API-019` | `born_effective_charges` | `born_effective_charges(gs; tolerance::Real=1e-8, maxiter::Integer=200)` | real `N×3×3`; `polar` | insulating ground state; positive tolerance and `maxiter` |
| `API-020` | `dielectric_tensor` | `dielectric_tensor(gs; tolerance::Real=1e-8, maxiter::Integer=200)` | symmetric real `3×3`; `polar` | insulating ground state; positive tolerance and `maxiter` |

No other symbol is exported in V0. Ordinary Base and LinearAlgebra methods
needed for use of these objects may be defined but are not additional exports.

### Public read-only object properties

The following `getproperty` names are part of the public boundary. They may be
stored or computed, but must have the stated canonical semantics; all other
properties are implementation details.

| object | public properties |
|---|---|
| `Crystal` | `lattice`, `species`, `positions`, `masses`; canonical native units and shapes from the Units section |
| `KSModel` | `crystal`, `electron_count`, `xc`; electron count is a finite integer-valued real |
| `PlaneWaveBasis` | `model`, `Ecut`, `kpoints::Vector{NTuple{3,Float64}}`, `kweights::Vector{Float64}`, `G_vectors::Vector{Vector{NTuple{3,Int}}}`, `fft_size::NTuple{3,Int}`; k points use reduced reciprocal coordinates and one G-vector list exists per full k point |
| `SCFOptions` | `energy_tolerance`, `density_tolerance`, `maxiter`, `extra_bands` |
| `AtomicDisplacement` | `atom`, `direction`, `q` |
| `QEInput` | `model`, `basis`, `options`; the canonical objects implied by the parsed input |

Arrays returned through these properties are observationally read-only:
mutating a retrieved array must not mutate the originating object.

`density(gs).values[i,j,k]` samples the periodic density at reduced coordinate
`((i-1)/n1,(j-1)/n2,(k-1)/n3)`, where `(n1,n2,n3)=size(values)`. The complex
`response(...).delta_density` uses the same grid and origin. This convention
freezes low-order reciprocal coefficients without freezing an FFT library.

## Differentiability contract

Differentiability is a required behaviour of the existing public API and does
not add exported symbols. For fixed nondifferentiable configuration `c`, define
the public composite

```text
z*(theta) = ground_state(PlaneWaveBasis(KSModel(Crystal(theta; c); c); c); c)
```

where `theta` contains a Cartesian position matrix or a scalar homogeneous-
strain amplitude. A `ChainRulesCore`-compatible AD frontend must be able to
differentiate real scalar functions formed from `energy(z*)`, `forces(z*)`, or
real weighted contractions of `density(z*).values`. The package may implement
custom rules at any granularity, but the following semantics are frozen.

- `DIF-001`: for Cartesian differentiation, `Crystal` is constructed with
  `positions_are_fractional=false`; the gradient of converged energy with
  respect to the `3×N` Cartesian position matrix equals `-forces(gs)`.
- `DIF-002`: for strain differentiation, `h(t)=(I+t*eta)h` with symmetric
  dimensionless `eta` and fixed fractional positions. At `t=0`,
  `dE/dt = -Omega * sum(stress(gs) .* eta)`, following the QE convention in
  which positive stress is compressive.
- `DIF-003`: derivatives use the reciprocal integer G lists and FFT topology
  selected at the primal point. This frozen-topology tangent convention is
  piecewise differentiable and does not define derivatives with respect to
  `Ecut`, `kgrid`, `fft_size`, band count, occupation choices, or q labels.
  Visible finite-difference checks use perturbations that preserve topology.
- `DIF-004`: `ground_state` is differentiated as the converged constrained
  stationary equation `F(z*(theta),theta)=0`. A JVP solves
  `F_z*delta_z=-F_theta*delta_theta`; a VJP solves the corresponding adjoint
  system. Differentiating stored SCF, mixing, or eigensolver iterations is not
  the contracted derivative.
- `DIF-005`: complex-valued orbital, response, FFT, projector, eigenspace, and
  linear-solve rules use the real scalar objective convention
  `real(dot(cotangent, tangent))`. JVP and VJP contractions obey adjoint duality
  under this convention within `ad-duality` tolerance.
- `DIF-006`: eigensolver and orthogonalization derivatives are defined for the
  occupied subspace projector. They are invariant under phase changes and
  arbitrary unitary rotations inside a degenerate occupied subspace; no
  derivative of an individually labelled degenerate eigenvector is promised.
- `DIF-007`: the Gamma density JVP for a Cartesian atomic direction equals
  `response(gs, AtomicDisplacement(...)).delta_density` within `ad-gradient`
  tolerance. Conversely, the VJP of any finite real weighted density-grid
  contraction satisfies `DIF-005` against that direct response.
- `DIF-008`: differentiating `forces(gs)` with respect to Cartesian positions
  yields the negative mass-unweighted Gamma force constants. Equivalently,
  `d forces[I,a] / d R[J,b] = -sqrt(M[I]*M[J]) * D[I,a,J,b]` under the matrix
  indexing convention above, within `ad-second` tolerance.
- `DIF-009`: changing SCF convergence history while reaching the same
  stationary solution does not change derivatives beyond the error implied by
  the primal and adjoint residuals. The reverse pass stores no tape whose size
  grows with the number of completed SCF iterations.
- `DIF-010`: an unconverged primal or adjoint/response solve raises
  `ErrorException` containing `did not converge`; it must not return a zero,
  NaN, truncated-unroll, or falsely converged derivative.

Fractional coordinate wrapping, a change of discrete basis topology, a gap
closure, or a change of occupation is outside a single differentiable chart.
V0 fixtures avoid these boundaries. V0 guarantees the first derivatives and
selected second derivatives above; arbitrary nesting to third and higher order
is not part of this specification.

## Detailed behavioural clauses

### Crystal and model

- `CRY-001`: fractional positions are reduced modulo one without changing
  species order; Cartesian inputs are converted using the supplied lattice.
- `CRY-002`: translating any atom by an integer reduced lattice vector leaves
  all scientific observables invariant.
- `CRY-003`: a singular lattice, non-finite coordinate, non-positive mass, or
  species/position/mass mismatch raises `ArgumentError` before model
  construction.
- `MOD-001`: each element must map to an existing UPF 2.0.1 file whose element,
  valence, functional compatibility, and norm-conserving type are validated.
- `MOD-002`: unsupported `xc`, charge, spin, relativistic mode, USPP, or PAW
  input raises `ArgumentError` naming the unsupported feature.

### Plane-wave discretization and ground state

- `BAS-001`: the reciprocal basis contains exactly the `G` vectors satisfying
  `|k+G|^2/2 <= Ecut`, using one consistent boundary convention.
- `BAS-002`: the FFT grid is sufficient for all represented density Fourier
  components; an explicitly insufficient grid raises `ArgumentError`.
- `BAS-003`: full k-point weights sum to one and are invariant under a
  permutation of the enumeration.
- `SCF-001`: the returned state satisfies both declared convergence tolerances;
  failure within `maxiter` raises `ErrorException` containing `did not converge`.
- `SCF-002`: independent calls with the same basis and options return energy,
  density, eigenvalues, and occupations within their named tolerances.
- `SCF-003`: integrated density equals the pseudopotential valence-electron
  count within `2e-7` electrons.
- `SCF-004`: energy is stationary at the returned electronic ground state;
  its nuclear and strain derivatives agree with the public force and stress
  observables as specified by `OBS-002` and `OBS-003`.
- `SCF-005`: results are invariant, within scientific tolerances, under a
  consistent permutation of atoms.

### Ground-state observables

- `OBS-001`: `energy`, `forces`, and `stress` cannot be evaluated on an
  unconverged or basis-inconsistent state.
- `OBS-002`: forces are derivatives of the converged discrete total energy;
  central finite differences converge to them as displacement shrinks until
  SCF and basis errors dominate.
- `OBS-003`: stress is symmetric and converges to the central finite
  difference of energy with respect to symmetric strain.
- `OBS-004`: `density(gs).values` is real, finite, periodic, and paired with
  the cell volume needed to integrate it.
- `OBS-005`: eigenvalue and occupation outer vectors have one entry per full
  electronic k point and matching band lengths.

### Response and phonons

- `RSP-001`: response differentiates the converged stationary equations; it
  must not depend on how many SCF iterations happened before convergence.
- `RSP-002`: construction rejects `atom<1`, directions outside `1:3`, and
  non-finite q. Evaluation rejects an atom larger than the crystal size and a
  q that does not map electronic states from k to k+q by permuting the full
  k mesh.
- `RSP-003`: response residual norm is no larger than the requested tolerance
  when `converged==true`; exhausting `maxiter` raises `ErrorException`
  containing `did not converge` rather than returning a false result.
- `DYN-001`: `dynamical_matrix(gs,q)` is Hermitian within `5e-8` spectral norm.
- `DYN-002`: `D(-q) = conj(D(q))` within the response tolerance.
- `DYN-003`: at Gamma, the returned Hermitian matrix has the acoustic sum rule
  projected in the mass-weighted translational subspace; its three uniform
  Cartesian translations have restoring spectral norm at most `5e-8`.
- `DYN-004`: consistently permuting atoms applies the corresponding block
  permutation to the dynamical matrix and preserves its spectrum.
- `PHN-001`: `phonon_modes` diagonalizes the returned dynamical matrix using
  mass-normalized Cartesian eigenvectors.
- `PHN-002`: eigenvector phases and bases within a degenerate subspace are not
  fixed; verification compares projectors, not individual vector columns.

### Polar response

- `POL-001`: the acoustic sum over atomic Born effective charges vanishes
  within `5e-5` elementary charges for converged fixtures.
- `POL-002`: the dielectric tensor is symmetric and positive definite for the
  valid insulating fixtures.
- `POL-003`: a non-analytic correction constructed from the reported Born
  charges and dielectric tensor reproduces the direction-dependent Gamma
  limit and LO-TO splitting within `phonon` tolerance.

### QE 7.5 compatibility input

V0 accepts only plain-text `pw.x` SCF inputs with:

- namelists `CONTROL`, `SYSTEM`, and `ELECTRONS`;
- `calculation='scf'`, `ibrav=0`, `occupations='fixed'`, and `nspin=1`;
- `nat`, `ntyp`, `ecutwfc`, optional `ecutrho`, optional `tot_charge=0`, and
  `input_dft` selecting LDA or PBE, plus either `celldm(1)` in bohr or `A` in
  angstrom when a `CELL_PARAMETERS alat` card requires a lattice scale;
- `ATOMIC_SPECIES`, `ATOMIC_POSITIONS` in `crystal`, `bohr`, or `angstrom`,
  `CELL_PARAMETERS` in `bohr`, `angstrom`, or `alat`, and
  `K_POINTS automatic`;
- `prefix`, `pseudo_dir`, `conv_thr`, `electron_maxstep`, `tprnfor`, and
  `tstress` with their documented QE 7.5 semantics.

- `QEI-001`: namelist and card names are case-insensitive; quoted string values
  preserve their documented content semantics.
- `QEI-002`: QE Rydberg quantities are converted to native Hartree quantities
  exactly by the factor two; length units follow the card qualifier; and
  `ATOMIC_SPECIES` masses are converted from unified atomic mass units to
  electron masses.
- `QEI-003`: relative pseudopotential paths are resolved against `pseudo_dir`
  and validated before SCF.
- `QEI-004`: an equivalent native model and parsed QE model give results within
  the corresponding scientific comparison tolerance.
- `QEI-005`: every unsupported card, non-default unsupported keyword, malformed
  count, or inconsistent species declaration raises `ArgumentError` naming the
  field or card. Unknown input is never silently ignored.
- `QEI-006`: terminal prose need not match QE. V0 does not promise QE binary
  restart files or byte-for-byte `prefix.save` compatibility.

## Verification and scope rules

- Every `API-*` symbol and every behavioural clause above maps to at least one
  row in `TESTS.md` and at least one private verification test.
- Private tests may hide structures, numerical values, and parameter
  combinations, but may not introduce behaviour absent from this contract.
- Verification uses QE 7.5.0, identified pseudopotential-family versions and
  file SHA-256 values, and original clean-room inputs.
- Differentiability verification uses a pinned ChainRules-compatible Julia AD
  frontend only as a consumer of the candidate's declared rules. The candidate
  must not depend on that frontend for ordinary primal calculations.
- Passing visible examples is necessary but not sufficient. Completion also
  requires private unit, property, AD duality, differential, and end-to-end
  tests.
- Performance measurements are reported but do not fail V0 correctness unless
  a calculation exceeds its declared resource ceiling or cannot finish.
