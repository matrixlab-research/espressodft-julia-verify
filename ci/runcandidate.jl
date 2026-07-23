#!/usr/bin/env julia

using Pkg

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RUNNER_ENV = mktempdir(prefix="espressodft-candidate-")
const RUNNER_ROOT = joinpath(RUNNER_ENV, "verify")
mkpath(RUNNER_ROOT)
for source in readdir(ROOT; join=true)
    basename(source) == ".git" && continue
    cp(source, joinpath(RUNNER_ROOT, basename(source)); force=true)
end
isfile(joinpath(RUNNER_ROOT, "Manifest.toml")) ||
    error("the committed verification Manifest.toml is required")
Pkg.activate(RUNNER_ROOT)
Pkg.instantiate()

repository = get(ENV, "CANDIDATE_REPOSITORY", "matrixlab-research/EspressoDFT.jl")
reference = get(ENV, "CANDIDATE_REF", "main")
candidate_path = get(ENV, "CANDIDATE_PATH", "")

if isempty(candidate_path)
    url = startswith(repository, "http") ? repository : "https://github.com/$repository.git"
    @info "installing candidate" repository reference
    Pkg.add(PackageSpec(url=url, rev=reference); preserve=Pkg.PRESERVE_ALL)
else
    @info "developing local candidate" candidate_path
    Pkg.develop(PackageSpec(path=abspath(candidate_path)); preserve=Pkg.PRESERVE_ALL)
end
Pkg.instantiate()

candidate = only(filter(pair -> pair.second.name == "EspressoDFT", Pkg.dependencies()))
@info "resolved candidate" name=candidate.second.name version=candidate.second.version tree_hash=candidate.second.tree_hash

include(joinpath(RUNNER_ROOT, "test", "runtests.jl"))
