@testset "RT-001 implicit response reproducibility" begin
    fixture = si_fixture()
    perturbation = AtomicDisplacement(2, 1, (0.25, 0.25, 0.25))
    first = response(candidate_state(fixture), perturbation; tolerance=1e-8)
    independent = ground_state(candidate_basis(fixture); options=SCFOptions(maxiter=120))
    second = response(independent, perturbation; tolerance=1e-8)
    @test first.converged && second.converged
    @test first.residual_norm <= 1e-8
    @test second.residual_norm <= 1e-8
    @test first.delta_density ≈ second.delta_density atol=5e-5 rtol=5e-5
end

@testset "RT-002 response nonconvergence is explicit" begin
    perturbation = AtomicDisplacement(2, 1, (0.25, 0.25, 0.25))
    @test_throws ErrorException response(candidate_state(si_fixture()), perturbation;
                                         tolerance=eps())
end

@testset "RT-003 Hermitian Gamma and non-Gamma matrices" begin
    gs = candidate_state(si_fixture())
    for q in ((0.0, 0.0, 0.0), (0.25, 0.25, 0.25))
        matrix = dynamical_matrix(gs, q; tolerance=1e-8)
        @test size(matrix) == (6, 6)
        @test opnorm(matrix - matrix') <= 5e-8
    end
end

@testset "RT-004 q conjugacy" begin
    gs = candidate_state(si_fixture())
    q = (0.25, 0.25, 0.25)
    @test dynamical_matrix(gs, ntuple(i -> -q[i], 3)) ≈
          conj(dynamical_matrix(gs, q)) atol=5e-5 rtol=5e-5
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
    @test original[indices, indices] ≈ permuted atol=5e-5 rtol=5e-5
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
        @test charges[atom, :, :] ≈ rows_to_matrix(expected[atom]) atol=5e-4 rtol=5e-4
    end
end

@testset "RT-010 dielectric symmetry and positivity" begin
    fixture = nacl_fixture()
    dielectric = dielectric_tensor(candidate_state(fixture); tolerance=1e-8)
    expected = rows_to_matrix(
        fixture_reference(fixture)["phonons"]["0.0,0.0,0.0"]["dielectric"])
    @test dielectric ≈ dielectric' atol=5e-8
    @test isposdef(Hermitian(dielectric))
    @test dielectric ≈ expected atol=5e-3 rtol=5e-3
end

@testset "RT-011 direction-dependent non-analytic correction" begin
    fixture = nacl_fixture()
    gs = candidate_state(fixture)
    charges = born_effective_charges(gs)
    dielectric = dielectric_tensor(gs)
    masses = native_masses(fixture)
    volume = abs(det(fixture.lattice_bohr))
    function correction(direction)
        direction = normalize(collect(direction))
        matrix = zeros(6, 6)
        denominator = dot(direction, dielectric * direction)
        for i in 1:2, j in 1:2, a in 1:3, b in 1:3
            zi = dot(direction, charges[i, :, a])
            zj = dot(direction, charges[j, :, b])
            matrix[3(i - 1) + a, 3(j - 1) + b] =
                4pi / volume * zi * zj / (denominator * sqrt(masses[i] * masses[j]))
        end
        matrix
    end
    x = correction((1.0, 0.0, 0.0))
    diagonal = correction((1.0, 1.0, 1.0))
    @test x ≈ x' atol=5e-8
    @test diagonal ≈ diagonal' atol=5e-8
    @test sort(eigvals(Hermitian(x))) ≈ sort(eigvals(Hermitian(diagonal))) atol=5e-5
    @test maximum(eigvals(Hermitian(x))) > 0
end

@testset "IT-003 Si response-to-phonon workflow" begin
    fixture = si_fixture()
    gs = candidate_state(fixture)
    q = (0.25, 0.25, 0.25)
    got = dynamical_matrix(gs, q)
    expected = oracle_dynamical(fixture, "0.25,0.25,0.25")
    @test got ≈ expected atol=5e-5 rtol=5e-5
    got_cm1 = phonon_modes(gs, q).frequencies .* CM1_PER_ATOMIC_FREQUENCY
    expected_cm1 = fixture_reference(fixture)["phonons"]["0.25,0.25,0.25"]["frequencies_cm1_raw"]
    @test sort(got_cm1) ≈ sort(expected_cm1) atol=2.0
end

@testset "IT-004 NaCl polar workflow" begin
    fixture = nacl_fixture()
    gs = candidate_state(fixture)
    @test size(born_effective_charges(gs)) == (2, 3, 3)
    @test isposdef(Hermitian(dielectric_tensor(gs)))
end

@testset "IT-005 held-out complete API denominator" begin
    @test REFERENCE["spec_id"] == SPEC_ID
    @test Set(keys(fixture_reference(si_fixture())["phonons"])) ==
          Set(["0.0,0.0,0.0", "0.25,0.25,0.25"])
    @test all(symbol -> isdefined(QuantumDFT, symbol), EXPECTED_EXPORTS)
end
