# EspressoDFT private verification

Private black-box verification for frozen specification
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
   installed into this project and tested only through its documented public
   boundary.

The candidate gate also treats the converged ground state as a differentiable
implicit layer. A pinned ChainRules-compatible consumer checks energy/force and
energy/stress gradients, direct/adjoint density-response duality, selected
second derivatives, failure semantics, and bounded pullback storage. The
candidate is not required to depend on the verifier's AD frontend.

Private tests hide structures, values, and parameter combinations. They do
not add undocumented semantics. QE source and QE tests are not copied here;
QE is executed as a pinned black-box oracle.

## Local commands

```bash
python3 ci/check_contract_coverage.py
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. oracle/generate_reference.jl --check
CANDIDATE_REPOSITORY=owner/EspressoDFT.jl CANDIDATE_REF=main \
  julia --project=. ci/runcandidate.jl
```

`--update` deliberately requires `ALLOW_ORACLE_UPDATE=1`. An oracle update
must accompany a new specification ID or an explained, reviewed correction.
