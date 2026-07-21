function deterministic_density_cotangent(values)
    seed = reshape(sin.(collect(1:length(values))), size(values))
    seed ./ norm(seed)
end

@testset "AD-001 energy reverse gradient equals negative force" begin
    fixture = distorted_si_fixture()
    positions = cartesian_positions(fixture)
    gradient = only(Zygote.gradient(
        x -> candidate_energy_from_cartesian(fixture, x), positions))
    reference = forces(candidate_state_from_geometry(
        fixture; positions, positions_are_fractional=false))
    @test all(isfinite, gradient)
    @test norm(reference) > 1e-5
    @test gradient ≈ -reference atol=5e-5 rtol=5e-5
end

@testset "AD-002 strain reverse gradient equals stress contraction" begin
    fixture = aln_fixture()
    identity3 = Matrix{Float64}(I, 3, 3)
    strain = [0.7 0.2 0.0; 0.2 -0.4 0.1; 0.0 0.1 0.3]
    derivative = only(Zygote.gradient(0.0) do amplitude
        lattice = (identity3 + amplitude * strain) * fixture.lattice_bohr
        energy(candidate_state_from_geometry(fixture; lattice))
    end)
    gs = candidate_state(fixture)
    expected = abs(det(fixture.lattice_bohr)) * sum(stress(gs) .* strain)
    @test isfinite(derivative)
    @test abs(expected) > 1e-7
    @test isapprox(derivative, expected; atol=5e-5, rtol=5e-5)
end

@testset "AD-003 implicit gradient is independent of SCF history" begin
    fixture = distorted_si_fixture()
    positions = cartesian_positions(fixture)
    loose = SCFOptions(
        energy_tolerance=1e-9,
        density_tolerance=1e-7,
        maxiter=120,
        extra_bands=4,
    )
    tight = SCFOptions(
        energy_tolerance=1e-11,
        density_tolerance=1e-9,
        maxiter=180,
        extra_bands=4,
    )
    loose_gradient = only(Zygote.gradient(
        x -> candidate_energy_from_cartesian(fixture, x; options=loose), positions))
    tight_gradient = only(Zygote.gradient(
        x -> candidate_energy_from_cartesian(fixture, x; options=tight), positions))
    tight_state = candidate_state_from_geometry(
        fixture;
        positions,
        positions_are_fractional=false,
        options=tight,
    )
    @test loose_gradient ≈ tight_gradient atol=5e-5 rtol=5e-5
    @test tight_gradient ≈ -forces(tight_state) atol=5e-5 rtol=5e-5

    basis = candidate_basis_from_geometry(
        fixture; positions, positions_are_fractional=false)
    primal, pullback = ChainRulesCore.rrule(ground_state, basis; options=tight)
    retained = Base.summarysize(pullback)
    state_scale = Base.summarysize(primal) + Base.summarysize(basis)
    @test retained <= 12state_scale
end

@testset "AD-004 Gamma density JVP matches finite difference" begin
    fixture = distorted_si_fixture()
    positions = cartesian_positions(fixture)
    step = 2e-3
    plus_positions = copy(positions)
    minus_positions = copy(positions)
    plus_positions[1, 2] += step
    minus_positions[1, 2] -= step

    basis = candidate_basis_from_geometry(
        fixture; positions, positions_are_fractional=false)
    plus_basis = candidate_basis_from_geometry(
        fixture; positions=plus_positions, positions_are_fractional=false)
    minus_basis = candidate_basis_from_geometry(
        fixture; positions=minus_positions, positions_are_fractional=false)
    @test plus_basis.fft_size == basis.fft_size == minus_basis.fft_size
    @test plus_basis.G_vectors == basis.G_vectors == minus_basis.G_vectors

    gs = ground_state(basis; options=SCFOptions(maxiter=120))
    direct = response(gs, AtomicDisplacement(2, 1, (0.0, 0.0, 0.0)))
    plus = density(ground_state(plus_basis; options=SCFOptions(maxiter=120))).values
    minus = density(ground_state(minus_basis; options=SCFOptions(maxiter=120))).values
    finite_difference = (plus - minus) / (2step)
    @test size(direct.delta_density) == size(finite_difference)
    @test norm(direct.delta_density) > 1e-8
    @test norm(imag.(direct.delta_density)) <= 5e-5 + 5e-5norm(finite_difference)
    @test real.(direct.delta_density) ≈ finite_difference atol=5e-5 rtol=5e-5
end

@testset "AD-005 direct and adjoint density responses are dual" begin
    fixture = distorted_si_fixture()
    positions = cartesian_positions(fixture)
    gs = candidate_state_from_geometry(
        fixture; positions, positions_are_fractional=false)
    rho = density(gs)
    cotangent = deterministic_density_cotangent(rho.values)
    direct = response(gs, AtomicDisplacement(2, 1, (0.0, 0.0, 0.0)))
    lhs = sum(cotangent .* real.(direct.delta_density)) *
          rho.cell_volume / length(rho.values)
    gradient = only(Zygote.gradient(
        x -> candidate_density_contraction(fixture, x, cotangent), positions))
    rhs = gradient[1, 2]
    @test abs(lhs) > 1e-10
    @test isapprox(lhs, rhs; atol=2e-6, rtol=2e-5)
end

@testset "AD-006 force derivative equals Gamma force constant" begin
    fixture = distorted_si_fixture()
    positions = cartesian_positions(fixture)
    force_gradient = only(Zygote.gradient(positions) do x
        gs = candidate_state_from_geometry(
            fixture; positions=x, positions_are_fractional=false)
        forces(gs)[1, 2]
    end)
    gs = candidate_state_from_geometry(
        fixture; positions, positions_are_fractional=false)
    matrix = dynamical_matrix(gs, (0.0, 0.0, 0.0))
    masses = native_masses(fixture)
    expected = -sqrt(masses[2] * masses[1]) * real(matrix[4, 1])
    @test abs(expected) > 1e-8
    @test isapprox(force_gradient[1, 1], expected; atol=2e-7, rtol=5e-4)
end

@testset "AD-007 nonconvergent differentiation fails closed" begin
    fixture = distorted_si_fixture()
    positions = cartesian_positions(fixture)
    impossible = SCFOptions(
        energy_tolerance=eps(),
        density_tolerance=eps(),
        maxiter=1,
        extra_bands=0,
    )
    error = try
        Zygote.gradient(
            x -> candidate_energy_from_cartesian(fixture, x; options=impossible),
            positions,
        )
        nothing
    catch caught
        caught
    end
    @test error isa ErrorException
    @test occursin("did not converge", lowercase(sprint(showerror, error)))
end

@testset "IT-006 held-out differentiable stationary workflow" begin
    fixture = aln_fixture()
    positions = cartesian_positions(fixture)
    gs = candidate_state_from_geometry(
        fixture; positions, positions_are_fractional=false)

    energy_gradient = only(Zygote.gradient(
        x -> candidate_energy_from_cartesian(fixture, x), positions))
    @test energy_gradient ≈ -forces(gs) atol=5e-5 rtol=5e-5

    rho = density(gs)
    cotangent = deterministic_density_cotangent(rho.values)
    direct = response(gs, AtomicDisplacement(2, 1, (0.0, 0.0, 0.0)))
    lhs = sum(cotangent .* real.(direct.delta_density)) *
          rho.cell_volume / length(rho.values)
    density_gradient = only(Zygote.gradient(
        x -> candidate_density_contraction(fixture, x, cotangent), positions))
    @test isapprox(lhs, density_gradient[1, 2]; atol=2e-6, rtol=2e-5)

    force_gradient = only(Zygote.gradient(positions) do x
        state = candidate_state_from_geometry(
            fixture; positions=x, positions_are_fractional=false)
        forces(state)[1, 2]
    end)
    matrix = dynamical_matrix(gs, (0.0, 0.0, 0.0))
    masses = native_masses(fixture)
    expected = -sqrt(masses[2] * masses[1]) * real(matrix[4, 1])
    @test isapprox(force_gradient[1, 1], expected; atol=2e-7, rtol=5e-4)

    basis = candidate_basis_from_geometry(
        fixture; positions, positions_are_fractional=false)
    primal, pullback = ChainRulesCore.rrule(
        ground_state, basis; options=SCFOptions(maxiter=120))
    @test Base.summarysize(pullback) <=
          12(Base.summarysize(primal) + Base.summarysize(basis))
end
