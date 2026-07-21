using Test
using EspressoDFTVerify

const ROOT = normpath(joinpath(@__DIR__, ".."))

@testset "VT-001 structural coverage cannot be satisfied by mentions" begin
    checker = joinpath(ROOT, "ci", "check_contract_coverage.py")
    @test success(`python3 $checker`)
end

@testset "VT-002 oracle comparator rejects scientific mutations" begin
    comparator = EspressoDFTVerify.compare_reference
    @test comparator(1.0, 1.0, "fixtures.x.energy_ha")
    for (path, delta) in (
        ("fixtures.x.energy_ha", 1e-4),
        ("fixtures.x.gamma_eigenvalues_ha[1]", 1e-3),
        ("fixtures.x.density_fourier.1,0,0[1]", 1e-3),
        ("fixtures.x.phonons.q.dynamical_real[1]", 1e-8),
        ("fixtures.x.phonons.q.frequencies_cm1_raw[1]", 2.0),
        ("fixtures.x.phonons.q.born_charges_asr[1]", 1e-3),
        ("fixtures.x.phonons.q.dielectric[1]", 1e-2),
    )
        @test_throws ErrorException comparator(1.0 + delta, 1.0, path)
    end
    @test_throws ErrorException comparator(Dict{String,Any}(), Dict("required" => 1),
                                           "fixtures.x")
end

@testset "VT-003 zero-valued candidate is rejected" begin
    runner = joinpath(ROOT, "ci", "runcandidate.jl")
    sentinel = joinpath(ROOT, "test", "sentinels", "ZeroEspressoDFT")
    command = addenv(`$(Base.julia_cmd()) --project=$ROOT $runner`,
                     "CANDIDATE_PATH" => sentinel)
    mktemp() do path, io
        close(io)
        process = run(pipeline(ignorestatus(command), stdout=path, stderr=path))
        log = read(path, String)
        @test process.exitcode != 0
        @test occursin("Test Failed", log)
        @test occursin(r"(?:UT|PT|RT|IT)-\d{3}", log)
    end
end

@testset "VT-004 oracle CI is cross-platform" begin
    workflow = read(joinpath(ROOT, ".github", "workflows", "verify.yml"), String)
    @test occursin("ubuntu-latest", workflow)
    @test occursin("macos-14", workflow)
    @test occursin("matrix.os", workflow)
end

@testset "VT-005 candidate CI has three invocation paths" begin
    workflow = read(joinpath(ROOT, ".github", "workflows", "verify.yml"), String)
    @test occursin("workflow_dispatch:", workflow)
    @test occursin("repository_dispatch:", workflow)
    @test occursin("workflow_call:", workflow)
    @test occursin("ci/runcandidate.jl", workflow)
end
