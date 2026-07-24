@testset "PT-010 bounded SCF history and diagnostics" begin
    gs = candidate_state(si_fixture())
    @test hasproperty(gs, :iterations)
    @test hasproperty(gs, :energy_history)
    @test hasproperty(gs, :density_residual_history)

    iterations = gs.iterations
    energies = gs.energy_history
    residuals = gs.density_residual_history
    @test 2 <= iterations <= 40
    @test length(energies) == length(residuals) == iterations
    @test all(isfinite, energies)
    @test all(x -> isfinite(x) && x >= 0, residuals)
    @test last(energies) ≈ energy(gs) atol=5e-12 rtol=5e-12
    @test last(residuals) <= 1e-8
    @test last(residuals) < first(residuals)

    energies[1] = Inf
    residuals[1] = Inf
    @test isfinite(first(gs.energy_history))
    @test isfinite(first(gs.density_residual_history))
end
