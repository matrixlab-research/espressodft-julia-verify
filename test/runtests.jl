using Test
using LinearAlgebra
using PseudoPotentialData
using EspressoDFTVerify
using SHA
using TOML
using EspressoDFT

const EXPECTED_EXPORTS = Set([
    :Crystal, :KSModel, :PlaneWaveBasis, :SCFOptions, :QEInput,
    :AtomicDisplacement, :read_qe_input, :run_qe_input, :ground_state,
    :energy, :forces, :stress, :density, :eigenvalues, :occupations,
    :response, :dynamical_matrix, :phonon_modes, :born_effective_charges,
    :dielectric_tensor,
])

@testset "Frozen public surface" begin
    actual = Set(names(EspressoDFT; all=false, imported=false))
    delete!(actual, :EspressoDFT)
    @test actual == EXPECTED_EXPORTS
end

include("helpers.jl")
include("unit/objects.jl")
include("unit/qe_input.jl")
include("integration/ground_state.jl")
include("integration/response_phonons.jl")
