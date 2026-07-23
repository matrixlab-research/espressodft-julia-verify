@testset "RT-012 bounded Gamma response smoke" begin
    fixture = helium_smoke_fixture()
    gs = candidate_state_from_geometry(fixture; options=smoke_options())
    perturbation = AtomicDisplacement(1, 1, (0.0, 0.0, 0.0))
    result = response(gs, perturbation; tolerance=1e-6, maxiter=80)
    @test result.converged
    @test result.residual_norm <= 1e-6
    @test size(result.delta_density) == size(density(gs).values)
    @test all(isfinite, real.(result.delta_density))
    @test all(isfinite, imag.(result.delta_density))
    @test norm(result.delta_density) > 1e-6
    integrated = sum(real.(result.delta_density)) *
                 density(gs).cell_volume / length(result.delta_density)
    @test abs(integrated) <= 2e-7
end

@testset "AD-008 bounded density-response duality smoke" begin
    fixture = helium_smoke_fixture()
    positions = cartesian_positions(fixture)
    gs = candidate_state_from_geometry(
        fixture;
        positions,
        positions_are_fractional=false,
        options=smoke_options(),
    )
    rho = density(gs)
    cotangent = deterministic_density_cotangent(rho.values)
    direct = response(gs, AtomicDisplacement(1, 1, (0.0, 0.0, 0.0)))
    lhs = sum(cotangent .* real.(direct.delta_density)) *
          rho.cell_volume / length(rho.values)
    gradient = only(Zygote.gradient(positions) do x
        candidate_density_contraction(
            fixture, x, cotangent; options=smoke_options())
    end)
    @test all(isfinite, gradient)
    @test abs(lhs) > 1e-8
    @test isapprox(lhs, gradient[1, 1]; atol=2e-6, rtol=2e-5)
end
