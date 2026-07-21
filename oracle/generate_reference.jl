#!/usr/bin/env julia

using QuantumDFTVerify
using TOML

const ROOT = normpath(joinpath(@__DIR__, ".."))
const REFERENCE = joinpath(@__DIR__, "reference.toml")

function generate()
    fixtures = Dict{String,Any}()
    for fixture in (si_fixture(), nacl_fixture())
        @info "running pinned QE oracle" fixture=fixture.name
        fixtures[fixture.name] = run_oracle(fixture)
    end
    Dict(
        "spec_id" => SPEC_ID,
        "qe_version" => "7.5.0+0",
        "fixtures" => fixtures,
    )
end

generated = generate()
if "--print" in ARGS
    TOML.print(stdout, generated; sorted=true)
elseif "--update" in ARGS
    get(ENV, "ALLOW_ORACLE_UPDATE", "0") == "1" ||
        error("set ALLOW_ORACLE_UPDATE=1 for an intentional reviewed update")
    open(REFERENCE, "w") do io
        TOML.print(io, generated; sorted=true)
    end
    @info "updated oracle reference" path=REFERENCE
elseif "--check" in ARGS || isempty(ARGS)
    check_reference(load_reference(REFERENCE), generated)
    @info "pinned QE oracle matches checked-in reference" spec=SPEC_ID
else
    error("usage: generate_reference.jl [--check|--print|--update]")
end
