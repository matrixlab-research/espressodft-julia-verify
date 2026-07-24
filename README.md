# EspressoDFT verification

Public black-box verification for frozen specification
`espressodft-v0.3-qe7.5-2026-07-21`.

EspressoDFT is an independent MIT-licensed clean-room project. It is not
affiliated with or endorsed by the Quantum ESPRESSO Foundation; Quantum
ESPRESSO is named only as the external compatibility target and oracle.

This repository has three independent gates:

1. **contract coverage** — every frozen `API-*` and behavioural clause has a
   non-empty mapping to real test IDs;
2. **oracle integrity** — original inputs are rerun with
   `QuantumEspresso_jll` 7.5.0 and versioned NC-UPF artifacts, then compared to
   the checked-in observations; and
3. **candidate verification** — a requested `EspressoDFT.jl` repository/ref is
   injected into a disposable copy of the committed locked environment with
   `PRESERVE_ALL`, then tested only through its documented public boundary.
   The lock includes the candidate's generic numerical dependency closure as
   well as the verifier and oracle dependencies.

The default GitHub candidate job sets `VERIFY_PROFILE=ci`. It runs the frozen
surface, unit, ground-state property, ground-state integration, and a bounded
real-He Gamma-response/density-AD smoke gate, which fit the CI time budget.
The larger response/phonon and differentiability files are retained as the
explicit `full` profile; they are not counted as passing when skipped and are
intended for manually provisioned extended runs.

The full candidate profile also treats the converged ground state as a
differentiable implicit layer. A pinned ChainRules-compatible consumer checks
energy/force and energy/stress gradients, direct/adjoint density-response
duality, selected second derivatives, failure semantics, and bounded pullback
storage. The candidate is not required to depend on the verifier's AD frontend.

The verifier keeps implementation and oracle concerns separate. QE source and
QE tests are not copied here; QE is executed as a pinned black-box oracle.

## Local commands

```bash
python3 ci/check_contract_coverage.py
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. oracle/generate_reference.jl --check
CANDIDATE_REPOSITORY=owner/EspressoDFT.jl CANDIDATE_REF=main \
  VERIFY_PROFILE=ci \
  julia --project=. ci/runcandidate.jl
CANDIDATE_REPOSITORY=owner/EspressoDFT.jl CANDIDATE_REF=main \
  VERIFY_PROFILE=full \
  julia --project=. ci/runcandidate.jl
```

`--update` deliberately requires `ALLOW_ORACLE_UPDATE=1`. An oracle update
must accompany a new specification ID or an explained, reviewed correction.
