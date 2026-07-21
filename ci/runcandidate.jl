#!/usr/bin/env julia

using Pkg

const ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(ROOT)

repository = get(ENV, "CANDIDATE_REPOSITORY", "kunyuan/QuantumDFT.jl")
reference = get(ENV, "CANDIDATE_REF", "main")
url = startswith(repository, "http") ? repository : "https://github.com/$repository.git"

@info "installing candidate" repository reference
Pkg.add(PackageSpec(url=url, rev=reference))
Pkg.instantiate()

include(joinpath(ROOT, "test", "runtests.jl"))
