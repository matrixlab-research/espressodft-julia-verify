#!/usr/bin/env julia

using Pkg

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RUNNER_ENV = mktempdir(prefix="espressodft-candidate-")
Pkg.activate(RUNNER_ENV)
Pkg.develop(PackageSpec(path=ROOT))
# Candidate tests intentionally use fixture and AD-consumer APIs directly.
# Make them runner dependencies instead of relying on transitive loading from
# EspressoDFTVerify; the candidate itself need not depend on Zygote.
Pkg.add([
    PackageSpec(name="PseudoPotentialData",
                uuid=Base.UUID("5751a51d-ac76-4487-a056-413ecf6fbe19")),
    PackageSpec(name="ChainRulesCore",
                uuid=Base.UUID("d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4")),
    PackageSpec(name="Zygote",
                uuid=Base.UUID("e88e6eb3-aa80-5325-afca-941959d7151f")),
])

repository = get(ENV, "CANDIDATE_REPOSITORY", "kunyuan/EspressoDFT.jl")
reference = get(ENV, "CANDIDATE_REF", "main")
candidate_path = get(ENV, "CANDIDATE_PATH", "")

if isempty(candidate_path)
    url = startswith(repository, "http") ? repository : "https://github.com/$repository.git"
    @info "installing candidate" repository reference
    Pkg.add(PackageSpec(url=url, rev=reference))
else
    @info "developing local candidate" candidate_path
    Pkg.develop(PackageSpec(path=abspath(candidate_path)))
end
Pkg.instantiate()

candidate = only(filter(pair -> pair.second.name == "EspressoDFT", Pkg.dependencies()))
@info "resolved candidate" name=candidate.second.name version=candidate.second.version tree_hash=candidate.second.tree_hash

include(joinpath(ROOT, "test", "runtests.jl"))
