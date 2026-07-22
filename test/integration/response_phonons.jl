@testset "RT-001 implicit response reproducibility" begin
    fixture = si_fixture()
    perturbation = AtomicDisplacement(2, 1, (0.25, 0.25, 0.25))
    first = response(candidate_state(fixture), perturbation; tolerance=1e-8)
    independent = ground_state(candidate_basis(fixture); options=SCFOptions(maxiter=120))
    second = response(independent, perturbation; tolerance=1e-8)
    @test first.converged && second.converged
    @test first.residual_norm <= 1e-8
    @test second.residual_norm <= 1e-8
    @test ndims(first.delta_density) == 3
    @test all(isfinite, real.(first.delta_density))
    @test all(isfinite, imag.(first.delta_density))
    @test norm(first.delta_density) > 1e-8
    @test first.delta_density ≈ second.delta_density atol=5e-5 rtol=5e-5
end

@testset "RT-002 response nonconvergence is explicit" begin
    perturbation = AtomicDisplacement(2, 1, (0.25, 0.25, 0.25))
    error = try
        response(candidate_state(si_fixture()), perturbation;
                 tolerance=eps(), maxiter=1)
        nothing
    catch caught
        caught
    end
    @test error isa ErrorException
    @test occursin("did not converge", lowercase(sprint(showerror, error)))
end

@testset "RT-003 Hermitian Gamma and non-Gamma matrices" begin
    gs = candidate_state(si_fixture())
    for q in ((0.0, 0.0, 0.0), (0.25, 0.25, 0.25))
        matrix = dynamical_matrix(gs, q; tolerance=1e-8)
        @test size(matrix) == (6, 6)
        @test norm(matrix) > 1e-8
        @test opnorm(matrix - matrix') <= 5e-8
    end
end

@testset "RT-004 q conjugacy" begin
    gs = candidate_state(si_fixture())
    q = (0.25, 0.25, 0.25)
    @test dynamical_matrix(gs, ntuple(i -> -q[i], 3)) ≈
          conj(dynamical_matrix(gs, q)) atol=2e-9 rtol=2e-4
end

@testset "RT-005 projected acoustic translations" begin
    fixture = si_fixture()
    matrix = dynamical_matrix(candidate_state(fixture), (0.0, 0.0, 0.0))
    masses = native_masses(fixture)
    for direction in 1:3
        translation = zeros(6)
        for atom in 1:2
            translation[3(atom - 1) + direction] = sqrt(masses[atom])
        end
        normalize!(translation)
        @test norm(matrix * translation) <= 5e-8
    end
end

@testset "RT-006 atom block permutation" begin
    fixture = si_fixture()
    original = dynamical_matrix(candidate_state(fixture), (0.25, 0.25, 0.25))
    permutation = [2, 1]
    permuted_fixture = QEFixture(
        fixture.name * "-phonon-permuted", fixture.family, fixture.lattice_bohr,
        fixture.species[permutation], fixture.masses_amu[permutation],
        fixture.positions_fractional[:, permutation], fixture.ecut_ry,
        fixture.kgrid, fixture.input_dft, fixture.q_reduced,
    )
    permuted = dynamical_matrix(candidate_state(permuted_fixture), (0.25, 0.25, 0.25))
    indices = vcat(4:6, 1:3)
    @test original[indices, indices] ≈ permuted atol=2e-9 rtol=2e-4
end

@testset "RT-007 phonon diagonalization" begin
    fixture = si_fixture()
    q = (0.25, 0.25, 0.25)
    matrix = dynamical_matrix(candidate_state(fixture), q)
    modes = phonon_modes(candidate_state(fixture), q)
    @test modes.eigenvectors' * modes.eigenvectors ≈ I atol=5e-8
    @test modes.frequencies ≈ sign.(eigvals(Hermitian(matrix))) .*
          sqrt.(abs.(eigvals(Hermitian(matrix)))) atol=5e-8 rtol=5e-5
end

@testset "RT-008 degenerate subspace projector comparison" begin
    modes = phonon_modes(candidate_state(si_fixture()), (0.0, 0.0, 0.0))
    group = modes.eigenvectors[:, 1:3]
    rotated = group * Matrix(qr([1.0 2 3; 2 -1 1; 1 1 -1]).Q)
    @test group * group' ≈ rotated * rotated' atol=5e-8
end

@testset "RT-009 Born effective charges and acoustic sum" begin
    fixture = nacl_fixture()
    charges = born_effective_charges(candidate_state(fixture); tolerance=1e-8)
    expected = fixture_reference(fixture)["phonons"]["0.0,0.0,0.0"]["born_charges_asr"]
    @test size(charges) == (2, 3, 3)
    @test dropdims(sum(charges; dims=1); dims=1) ≈ zeros(3, 3) atol=5e-5
    for atom in 1:2
        assert_polar_known_issue(
            charges[atom, :, :], rows_to_matrix(expected[atom]);
            atol=BORN_ATOL, rtol=BORN_RTOL)
    end
end

@testset "RT-010 dielectric symmetry and positivity" begin
    fixture = nacl_fixture()
    dielectric = dielectric_tensor(candidate_state(fixture); tolerance=1e-8)
    expected = rows_to_matrix(
        fixture_reference(fixture)["phonons"]["0.0,0.0,0.0"]["dielectric"])
    @test dielectric ≈ dielectric' atol=5e-8
    @test isposdef(Hermitian(dielectric))
    @test dielectric ≈ expected atol=DIELECTRIC_ATOL rtol=DIELECTRIC_RTOL
end

@testset "RT-011 direction-dependent non-analytic correction" begin
    fixture = aln_fixture()
    gs = candidate_state(fixture)
    charges = born_effective_charges(gs)
    dielectric = dielectric_tensor(gs)
    masses = native_masses(fixture)
    volume = abs(det(fixture.lattice_bohr))
    analytic = dynamical_matrix(gs, (0.0, 0.0, 0.0))
    spectra = Dict{NTuple{3,Float64},Vector{Float64}}()
    for direction in ((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 1.0, 0.0))
        correction = nonanalytic_correction(charges, dielectric, masses, volume, direction)
        got = sort(signed_frequencies_cm1(analytic + correction))
        expected = sort(reference_nac_frequencies(fixture, direction))
        @test got ≈ expected atol=2.0
        spectra[direction] = got
    end
    @test maximum(abs.(spectra[(1.0, 0.0, 0.0)] .-
                       spectra[(0.0, 0.0, 1.0)])) > 1.0
end

@testset "IT-003 Si response-to-phonon workflow" begin
    fixture = si_fixture()
    gs = candidate_state(fixture)
    q = (0.25, 0.25, 0.25)
    got = dynamical_matrix(gs, q)
    expected = oracle_dynamical(fixture, "0.25,0.25,0.25")
    assert_dynamical_matches(got, expected)
    got_cm1 = phonon_modes(gs, q).frequencies .* CM1_PER_ATOMIC_FREQUENCY
    expected_cm1 = fixture_reference(fixture)["phonons"]["0.25,0.25,0.25"]["frequencies_cm1_raw"]
    @test sort(got_cm1) ≈ sort(expected_cm1) atol=2.0
end

@testset "IT-004 NaCl polar workflow" begin
    fixture = nacl_fixture()
    gs = candidate_state(fixture)
    charges, dielectric = reference_polar_tensors(fixture)
    assert_polar_known_issue(
        born_effective_charges(gs), charges; atol=BORN_ATOL, rtol=BORN_RTOL)
    @test isapprox(dielectric_tensor(gs), dielectric;
                   atol=DIELECTRIC_ATOL, rtol=DIELECTRIC_RTOL)
end

@testset "IT-005 held-out complete API denominator" begin
    fixture = aln_fixture()
    gs = candidate_state(fixture)
    reference = fixture_reference(fixture)
    @test REFERENCE["spec_id"] == SPEC_ID
    @test isapprox(energy(gs), reference["energy_ha"]; atol=2e-6, rtol=5e-8)
    @test forces(gs) ≈ rows_to_matrix(reference["forces_ha_per_bohr"])' atol=5e-6 rtol=5e-6
    @test stress(gs) ≈ rows_to_matrix(reference["stress_ha_per_bohr3"]) atol=5e-6 rtol=5e-6
    assert_density_matches(gs, fixture)
    assert_gamma_bands_match(gs, fixture)

    for (q, key) in (((0.0, 0.0, 0.0), "0.0,0.0,0.0"),
                     ((0.25, 0.0, 0.0), "0.25,0.0,0.0"))
        got = dynamical_matrix(gs, q)
        expected = oracle_dynamical(fixture, key)
        if all(iszero, q)
            expected = project_acoustic_sum_rule(expected, native_masses(fixture))
        end
        assert_dynamical_matches(got, expected)
        got_cm1 = sort(phonon_modes(gs, q).frequencies .* CM1_PER_ATOMIC_FREQUENCY)
        expected_cm1 = sort(signed_frequencies_cm1(expected))
        @test got_cm1 ≈ expected_cm1 atol=2.0
    end

    charges, dielectric = reference_polar_tensors(fixture)
    assert_polar_known_issue(
        born_effective_charges(gs), charges; atol=BORN_ATOL, rtol=BORN_RTOL)
    assert_polar_known_issue(
        dielectric_tensor(gs), dielectric;
        atol=DIELECTRIC_ATOL, rtol=DIELECTRIC_RTOL)

    mktempdir() do directory
        parsed = run_qe_input(private_qe_input(fixture, directory))
        @test isapprox(energy(parsed), energy(gs); atol=5e-7, rtol=5e-8)
    end
end
