# EspressoDFT — Visible unit and integration tests

> Frozen clean-room specification `espressodft-v0.2-qe7.5-2026-07-21`.
> These tests are the visible floor. The private suite uses different crystal
> parameters, pseudopotentials, q points, and numerical observations while
> testing only behaviour declared in `CONTRACT.md`.

## Reproducibility envelope

Visible scientific cases use QE 7.5.0 as the external oracle and scalar-
relativistic norm-conserving UPF files addressed by pseudopotential-family
identifier, element, and SHA-256. Each comparison fixes crystal, UPF bytes,
exchange-correlation functional, cutoff, FFT convention, full k mesh, q,
occupation, and convergence thresholds.

The visible fixtures are:

| fixture | purpose | crystal/model |
|---|---|---|
| `si-visible` | non-polar insulating SCF and phonons | two-atom diamond Si primitive cell; LDA NC-UPF; `Ecut=10 Ha`; `4×4×4` k mesh |
| `nacl-visible` | polar response and LO-TO behaviour | two-atom rocksalt NaCl primitive cell; PBE NC-UPF; converged cutoff and `4×4×4` k mesh |
| `invalid-visible` | rejection before numerical work | malformed cells, missing species, unsupported pseudopotentials, metals, spin, and non-commensurate q |
| `heldout-low-symmetry` | private complete workflow | four-atom non-cubic polar insulator; nonzero forces/stress; anisotropic dielectric tensor; Gamma and commensurate non-Gamma q |

Exact lattice values, pseudopotential identifiers, checksums, and oracle
observations used by the executable suite are recorded in the verification
repository. Visible tests below state the behaviour without disclosing the
private fixtures.

## Unit tests

| test ID | contract coverage | operation | expected behaviour |
|---|---|---|---|
| `UT-001` | `API-001`, `CRY-001` | construct a crystal with explicit masses from fractional and equivalent Cartesian positions | canonical species order, masses, lattice, and fractional coordinates agree |
| `UT-002` | `API-001`, `CRY-002` | translate one fractional position by an integer vector | canonical observable geometry is unchanged |
| `UT-003` | `CRY-003` | use a singular cell, non-finite position, non-positive mass, wrong shape, and species/mass-count mismatch | each case raises `ArgumentError` before model construction |
| `UT-004` | `API-002`, `MOD-001` | construct LDA and PBE models with matching NC-UPF files | construction succeeds and preserves electron count implied by UPF valence |
| `UT-005` | `MOD-002` | use real USPP/PAW artifacts and request spin, nonzero charge, unsupported XC, or mismatched UPF metadata | `ArgumentError` names the rejected feature |
| `UT-006` | `API-003`, `BAS-001` | independently enumerate a finite reciprocal bounding box | the returned set equals every and only `G` satisfying the declared cutoff convention |
| `UT-007` | `BAS-002` | give an insufficient explicit FFT grid | `ArgumentError` is raised; no aliased density grid is constructed |
| `UT-008` | `BAS-003` | construct several Monkhorst-Pack meshes | full k weights sum to one and enumeration permutation changes no weighted sum |
| `UT-009` | `API-004` | construct default and non-default SCF options | values and validation match the frozen signature; invalid values raise `ArgumentError` |
| `UT-010` | `API-006` | construct atomic displacements at Gamma and a commensurate non-Gamma q | one-based atom/direction and reduced q are preserved exactly |
| `UT-011` | `API-006`, `RSP-002` | use invalid atom/direction, atom beyond the crystal, or q that does not permute the k mesh | `ArgumentError` names the invalid component at construction or evaluation |
| `UT-012` | `API-007`, `API-005`, `QEI-001` | parse equivalent mixed-case QE namelists and cards | canonical `QEInput` semantics are identical |
| `UT-013` | `QEI-002` | parse independently written equivalent `bohr`, `angstrom`, `alat+celldm(1)`, Ry, unified-mass, and native atomic-unit quantities | canonical lattice, positions, masses, and cutoff agree within conversion roundoff |
| `UT-014` | `QEI-003` | resolve a relative UPF name through `pseudo_dir` | the validated canonical pseudopotential path identifies the expected bytes |
| `UT-015` | `QEI-005` | parse unsupported calculation, `ibrav`, occupation, spin, card, unknown keyword, malformed count, and inconsistent species cases | every case raises `ArgumentError` naming the field or card |

## Ground-state property tests

| test ID | contract coverage | operation | expected behaviour |
|---|---|---|---|
| `PT-001` | `API-009`, `SCF-001` | solve `si-visible` at frozen tolerances | converged state is returned and reported residuals satisfy both thresholds |
| `PT-002` | `SCF-001`, `OBS-001` | set an intentionally impossible tolerance and one iteration | `ground_state` raises `ErrorException` containing `did not converge`; accessors cannot treat it as converged |
| `PT-003` | `SCF-002` | repeat an independent solve with the same basis and options | energy, density, eigenvalues, and occupations agree within their named tolerances |
| `PT-004` | `API-013`, `SCF-003`, `OBS-004` | integrate `density(gs).values` using its cell volume and grid size | electron count agrees within `2e-7`; values are finite, real, and periodic |
| `PT-005` | `API-010`, `SCF-004` | take a directional derivative of energy along a mixed atomic displacement | derivative agrees with the corresponding projection of public forces within finite-difference tolerance |
| `PT-006` | `API-011`, `OBS-002` | use a non-equilibrium structure and at least three symmetric decreasing steps | finite-difference forces converge toward `forces(gs)` before numerical noise dominates; a constant-zero force fails |
| `PT-007` | `API-012`, `OBS-003` | apply symmetric cell strains and compare energy differences | finite-difference stress converges to `stress(gs)` and the tensor is symmetric |
| `PT-008` | `API-014`, `API-015`, `OBS-005` | inspect eigenvalues and occupations on the full mesh | outer lengths equal number of full k points; band lengths pair; weighted occupations equal electron count |
| `PT-009` | `SCF-005` | permute atoms consistently in the crystal and pseudopotential mapping | energy is unchanged; forces and atom-indexed data undergo only the corresponding permutation |

## Response and phonon property tests

| test ID | contract coverage | operation | expected behaviour |
|---|---|---|---|
| `RT-001` | `API-016`, `RSP-001`, `RSP-003` | solve the same Si displacement response from two independently converged ground states | response density agrees; `converged` is true and residual does not exceed tolerance |
| `RT-002` | `API-016`, `RSP-003` | force response nonconvergence with `maxiter=1` and an impossible tolerance | `ErrorException` containing `did not converge` is raised rather than returning a false result |
| `RT-003` | `API-017`, `DYN-001` | construct `D(q)` for `si-visible` at Gamma and `(1/4,1/4,1/4)` | both matrices have shape `6×6` and are Hermitian within `5e-8` spectral norm |
| `RT-004` | `DYN-002` | compute `D(q)` and `D(-q)` | the latter equals the complex conjugate of the former within response tolerance |
| `RT-005` | `DYN-003` | compute the Gamma dynamical matrix of `si-visible` | the three mass-weighted uniform translations have restoring norm at most `5e-8` |
| `RT-006` | `DYN-004` | permute the two Si atoms | a `3×3` block permutation relates the matrices and spectra are preserved |
| `RT-007` | `API-018`, `PHN-001` | compare `phonon_modes(gs,q)` with an independent Hermitian generalized eigensolve | signed frequencies and mass-normalized eigenspaces agree |
| `RT-008` | `PHN-002` | apply arbitrary phases and rotate a degenerate mode basis | subspace-projector comparison passes although individual columns differ |
| `RT-009` | `API-019`, `POL-001` | compute Born charges for `nacl-visible` | result shape is `2×3×3`; atomic sum vanishes within `5e-5` |
| `RT-010` | `API-020`, `POL-002` | compute the electronic dielectric tensor for `nacl-visible` | tensor is symmetric, positive definite, and matches QE within polar tolerance |
| `RT-011` | `POL-003` | construct direction-dependent non-analytic Gamma limits from reported polar tensors | anisotropic LO and TO limits match pinned QE-derived observations within `2 cm^-1` |

## End-to-end integration tests

| test ID | contract coverage | workflow | acceptance |
|---|---|---|---|
| `IT-001` | `API-001`–`API-004`, `API-009`–`API-015` | native `si-visible`: crystal → model → basis → SCF → all ground-state accessors | energy, forces, stress, selected density coefficients, eigenvalues, and occupations agree with QE under their named tolerances |
| `IT-002` | `API-005`, `API-007`, `API-008`, `QEI-004`, `QEI-006` | parse the equivalent QE input and call `run_qe_input` | native and QE-compatible paths return scientifically equivalent states; no terminal-text equality is required |
| `IT-003` | `API-006`, `API-016`–`API-018` | `si-visible`: converged ground state → atomic responses → non-Gamma dynamical matrix → phonon modes | dynamical-matrix entries and frequencies agree with QE 7.5; finite-difference supercell results provide an independent secondary check |
| `IT-004` | `API-019`, `API-020`, `POL-001`–`POL-003` | `nacl-visible`: ground state → electric/atomic response → Born/dielectric tensors → NAC Gamma limits | polar tensors and LO-TO splitting agree with QE 7.5 |
| `IT-005` | all API IDs and clauses through the mapping below | run the complete native ground-state, density/band, Gamma/non-Gamma phonon, Born, dielectric, and QE-input workflow on `heldout-low-symmetry` | all stored QE observations and structural invariants pass without relying on Si/NaCl golden values |

## Verifier and CI self-tests

| test ID | gate | acceptance |
|---|---|---|
| `VT-001` | coverage structure | every visible ID is defined exactly once as an `@testset` containing a real assertion; comments and manifest-only mentions do not count |
| `VT-002` | oracle comparator mutations | missing keys and above-tolerance mutations in energy, density, bands, dynamical matrices, frequencies, Born charges, and dielectric tensors are rejected |
| `VT-003` | fail-closed candidate | a local package exporting the right names but returning zero/placeholder results is rejected by the candidate runner |
| `VT-004` | cross-platform oracle | pinned QE observations regenerate on Linux and macOS within field-specific oracle reproducibility tolerances |
| `VT-005` | candidate integration | manual, repository-dispatch, or reusable-workflow invocation resolves the requested candidate ref and runs the private suite after contract and oracle gates |

## API and clause coverage matrix

This matrix is the denominator. The private verifier checks it mechanically so
that deleting an API, clause, or test cannot make the reported pass rate look
better.

| contract IDs | visible coverage | required private coverage |
|---|---|---|
| `API-001`, `CRY-001`–`CRY-003` | `UT-001`–`UT-003`, `IT-001` | invalid, boundary, translation, and atom-permutation variants |
| `API-002`, `MOD-001`–`MOD-002` | `UT-004`–`UT-005`, `IT-001` | independent LDA/PBE NC-UPF files and every excluded pseudopotential class |
| `API-003`, `BAS-001`–`BAS-003` | `UT-006`–`UT-008`, `IT-001` | cutoff-boundary, FFT-alias, shifted-mesh, and enumeration cases |
| `API-004` | `UT-009`, `PT-001`–`PT-002` | default, custom, invalid, and nonconvergent cases |
| `API-005`, `API-007`, `API-008`, `QEI-001`–`QEI-006` | `UT-012`–`UT-015`, `IT-002` | held-out units, ordering, malformed input, unsupported fields, and native-equivalence cases |
| `API-006` | `UT-010`–`UT-011`, `IT-003` | boundary atom/direction and multiple commensurate/non-commensurate q values |
| `API-009`, `SCF-001`–`SCF-005` | `PT-001`–`PT-005`, `PT-009`, `IT-001` | held-out structures, convergence histories, orthogonality, density, and permutations |
| `API-010`–`API-015`, `OBS-001`–`OBS-005` | `PT-004`–`PT-009`, `IT-001`–`IT-002` | differential QE goldens, finite differences, shape/type/error contracts |
| `API-016`, `RSP-001`–`RSP-003` | `RT-001`–`RT-002`, `IT-003` | held-out perturbations, solver histories, residual and rejection cases |
| `API-017`, `DYN-001`–`DYN-004` | `RT-003`–`RT-006`, `IT-003` | Gamma/non-Gamma matrices, conjugacy, translation, and permutation properties |
| `API-018`, `PHN-001`–`PHN-002` | `RT-007`–`RT-008`, `IT-003` | signed imaginary modes, degeneracies, phase, and subspace comparisons |
| `API-019`, `API-020`, `POL-001`–`POL-003` | `RT-009`–`RT-011`, `IT-004` | held-out polar crystal, tensor symmetry, sum rule, and direction-dependent NAC |

## Private-suite rules

- Hidden tests use only exported symbols and behaviour in `CONTRACT.md`.
- Hidden numerical values, structures, and parameter combinations are not
  reproduced in this visible document.
- The oracle is regenerated from pinned QE 7.5.0 and pseudopotential hashes
  before candidate tests; oracle drift is a verifier failure, not a candidate
  failure.
- Unit and integration assertions are correctness gates. Performance data is
  reported separately and is not a noisy pass/fail criterion in V0.
