@testset "PT-001 converged Si ground state" begin
    gs = candidate_state(si_fixture())
    @test isfinite(energy(gs))
    @test isapprox(energy(gs), fixture_reference(si_fixture())["energy_ha"];
                   atol=1e-6, rtol=5e-8)
end

@testset "PT-002 SCF nonconvergence is explicit" begin
    options = SCFOptions(energy_tolerance=eps(), density_tolerance=eps(),
                         maxiter=1, extra_bands=0)
    error = try
        ground_state(candidate_basis(si_fixture()); options)
        nothing
    catch caught
        caught
    end
    @test error isa ErrorException
    @test occursin("did not converge", sprint(showerror, error))
end

@testset "PT-003 independently repeated solves agree" begin
    fixture = si_fixture()
    first = candidate_state(fixture)
    second = ground_state(candidate_basis(fixture); options=SCFOptions(maxiter=120))
    @test isapprox(energy(first), energy(second); atol=5e-7, rtol=5e-8)
    @test density(first).values ≈ density(second).values atol=2e-5 rtol=2e-5
    @test all(isapprox(a, b; atol=5e-6, rtol=5e-6) for (a, b) in
              zip(eigenvalues(first), eigenvalues(second)))
    @test all(isapprox(a, b; atol=1e-12, rtol=1e-12) for (a, b) in
              zip(occupations(first), occupations(second)))
end

@testset "PT-004 density electron-number invariant" begin
    fixture = si_fixture()
    gs = candidate_state(fixture)
    rho = density(gs)
    @test ndims(rho.values) == 3
    @test eltype(rho.values) <: Real
    @test all(isfinite, rho.values)
    electron_count = sum(rho.values) * rho.cell_volume / length(rho.values)
    @test isapprox(electron_count, candidate_basis(fixture).model.electron_count; atol=2e-7)
end

@testset "PT-005 directional energy stationarity" begin
    fixture = si_fixture()
    step = 2e-3
    direction = normalize([1.0, -2.0, 0.5])
    plus = candidate_state(shifted_fixture(fixture, 2, step .* direction))
    minus = candidate_state(shifted_fixture(fixture, 2, -step .* direction))
    derivative = (energy(plus) - energy(minus)) / (2step)
    force_projection = dot(forces(candidate_state(fixture))[:, 2], direction)
    @test isapprox(derivative, -force_projection; atol=5e-5, rtol=5e-5)
end

@testset "PT-006 force finite difference" begin
    fixture = distorted_si_fixture()
    analytic = forces(candidate_state(fixture))[1, 2]
    @test abs(analytic) > 1e-5
    errors = Float64[]
    for step in (4e-3, 2e-3, 1e-3)
        plus = candidate_state(shifted_fixture(fixture, 2, (step, 0.0, 0.0)))
        minus = candidate_state(shifted_fixture(fixture, 2, (-step, 0.0, 0.0)))
        finite_difference = -(energy(plus) - energy(minus)) / (2step)
        push!(errors, abs(analytic - finite_difference))
    end
    @test errors[end] <= errors[1] + 2e-6
    @test errors[end] <= 5e-5 + 5e-5abs(analytic)
end

@testset "PT-007 stress finite difference and symmetry" begin
    fixture = aln_fixture()
    # Keep the exact integer G lists selected at the primal point, as required
    # by DIF-003.  At 5e-4 several AlN plane waves cross the hard cutoff and
    # the quotient measures a Pulay topology jump instead of V0 stress.
    step = 1e-4
    strain = zeros(3, 3)
    strain[1, 1] = step
    plus_fixture = strained_fixture(fixture, strain)
    minus_fixture = strained_fixture(fixture, -strain)
    @test candidate_basis(plus_fixture).G_vectors == candidate_basis(fixture).G_vectors
    @test candidate_basis(minus_fixture).G_vectors == candidate_basis(fixture).G_vectors
    plus = candidate_state(plus_fixture)
    minus = candidate_state(minus_fixture)
    volume = abs(det(fixture.lattice_bohr))
    finite_difference = (energy(plus) - energy(minus)) / (2step * volume)
    sigma = stress(candidate_state(fixture))
    @test sigma ≈ sigma' atol=5e-8
    @test abs(sigma[1, 1]) > 1e-6
    @test isapprox(sigma[1, 1], -finite_difference; atol=5e-5, rtol=5e-5)
end

@testset "PT-008 band and occupation structure" begin
    fixture = si_fixture()
    gs = candidate_state(fixture)
    bands = eigenvalues(gs)
    occs = occupations(gs)
    weights = candidate_basis(fixture).kweights
    @test length(bands) == length(weights) == length(occs)
    @test all(length(bands[k]) == length(occs[k]) for k in eachindex(bands))
    electron_count = sum(weights[k] * sum(occs[k]) for k in eachindex(weights))
    @test isapprox(electron_count, candidate_basis(fixture).model.electron_count; atol=1e-10)
end

@testset "PT-009 atom permutation covariance" begin
    fixture = aln_fixture()
    permutation = [3, 1, 4, 2]
    permuted = QEFixture(
        fixture.name * "-permuted", fixture.family, fixture.lattice_bohr,
        fixture.species[permutation], fixture.masses_amu[permutation],
        fixture.positions_fractional[:, permutation], fixture.ecut_ry,
        fixture.kgrid, fixture.input_dft, fixture.q_reduced,
    )
    original_state = candidate_state(fixture)
    permuted_state = candidate_state(permuted)
    @test isapprox(energy(original_state), energy(permuted_state); atol=5e-7, rtol=5e-8)
    @test forces(original_state)[:, permutation] ≈ forces(permuted_state) atol=5e-6 rtol=5e-6
end

@testset "IT-001 native ground-state workflow" begin
    fixture = si_fixture()
    gs = candidate_state(fixture)
    reference = fixture_reference(fixture)
    @test isapprox(energy(gs), reference["energy_ha"]; atol=1e-6, rtol=5e-8)
    @test forces(gs) ≈ rows_to_matrix(reference["forces_ha_per_bohr"])' atol=5e-6 rtol=5e-6
    @test stress(gs) ≈ rows_to_matrix(reference["stress_ha_per_bohr3"]) atol=5e-6 rtol=5e-6
    assert_density_matches(gs, fixture)
    assert_gamma_bands_match(gs, fixture)
    @test !isempty(occupations(gs))
end

@testset "IT-002 QE-compatible and native paths agree" begin
    mktempdir() do directory
        input_path = private_qe_input(si_fixture(), directory)
        parsed = read_qe_input(input_path)
        parsed_state = run_qe_input(parsed)
        path_state = run_qe_input(input_path)
        io_state = run_qe_input(IOBuffer(read(input_path, String)))
        native = candidate_state(si_fixture())
        @test isapprox(energy(parsed_state), energy(native); atol=5e-7, rtol=5e-8)
        @test isapprox(energy(path_state), energy(io_state); atol=5e-7, rtol=5e-8)
        @test forces(parsed_state) ≈ forces(native) atol=5e-6 rtol=5e-6
        @test stress(parsed_state) ≈ stress(native) atol=5e-6 rtol=5e-6
        assert_density_matches(parsed_state, si_fixture())
    end
end
