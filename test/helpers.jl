const REFERENCE = load_reference()
const STATE_CACHE = Dict{String,Any}()
const BASIS_CACHE = Dict{String,Any}()

native_masses(fixture::QEFixture) = fixture.masses_amu .* QuantumDFTVerify.AMU_TO_ELECTRON_MASS

function candidate_crystal(fixture::QEFixture)
    Crystal(
        fixture.lattice_bohr,
        fixture.species,
        fixture.positions_fractional;
        masses=native_masses(fixture),
    )
end

function candidate_basis(fixture::QEFixture)
    get!(BASIS_CACHE, fixture.name) do
        model = KSModel(
            candidate_crystal(fixture);
            pseudopotentials=pseudopotential_paths(fixture),
            xc=fixture.input_dft == "LDA" ? :lda : :pbe,
        )
        PlaneWaveBasis(model; Ecut=fixture.ecut_ry / 2, kgrid=fixture.kgrid)
    end
end

function candidate_state(fixture::QEFixture)
    get!(STATE_CACHE, fixture.name) do
        ground_state(candidate_basis(fixture); options=SCFOptions(
            energy_tolerance=1e-10,
            density_tolerance=1e-8,
            maxiter=120,
            extra_bands=4,
        ))
    end
end

fixture_reference(fixture::QEFixture) = REFERENCE["fixtures"][fixture.name]

function rows_to_matrix(rows)
    reduce(vcat, permutedims.(rows))
end

function oracle_dynamical(fixture::QEFixture, qkey::String)
    entry = fixture_reference(fixture)["phonons"][qkey]
    complex.(rows_to_matrix(entry["dynamical_real"]),
             rows_to_matrix(entry["dynamical_imag"]))
end

function assert_dynamical_matches(got, expected)
    @test size(got) == size(expected)
    @test norm(expected) > 1e-8
    @test isapprox(got, expected; atol=2e-9, rtol=2e-4)
end

candidate_density_coefficients(gs) =
    QuantumDFTVerify.density_coefficients(density(gs).values)

function assert_density_matches(gs, fixture::QEFixture)
    expected = fixture_reference(fixture)["density_fourier"]
    got = candidate_density_coefficients(gs)
    @test Set(keys(got)) == Set(keys(expected))
    for key in keys(expected)
        @test isapprox(got[key], expected[key]; atol=2e-5, rtol=2e-5)
    end
end

function assert_gamma_bands_match(gs, fixture::QEFixture)
    basis = candidate_basis(fixture)
    gamma = findfirst(k -> norm(collect(k)) <= 1e-12, basis.kpoints)
    @test gamma !== nothing
    expected = fixture_reference(fixture)["gamma_eigenvalues_ha"]
    got = eigenvalues(gs)[gamma]
    @test length(got) >= length(expected)
    @test isapprox(got[1:length(expected)], expected; atol=5e-6, rtol=5e-6)
end

function project_acoustic_sum_rule(matrix, masses)
    nat = length(masses)
    translations = zeros(eltype(matrix), 3nat, 3)
    for atom in 1:nat, direction in 1:3
        translations[3(atom - 1) + direction, direction] = sqrt(masses[atom])
    end
    translations = Matrix(qr(translations).Q[:, 1:3])
    projector = I - translations * translations'
    Hermitian(projector * matrix * projector)
end

function signed_frequencies_cm1(matrix)
    values = eigvals(Hermitian(matrix))
    sign.(values) .* sqrt.(abs.(values)) .* CM1_PER_ATOMIC_FREQUENCY
end

function nonanalytic_correction(charges, dielectric, masses, volume, direction)
    nat = length(masses)
    direction = normalize(collect(Float64, direction))
    correction = zeros(Float64, 3nat, 3nat)
    denominator = dot(direction, dielectric * direction)
    for i in 1:nat, j in 1:nat, a in 1:3, b in 1:3
        zi = dot(direction, charges[i, :, a])
        zj = dot(direction, charges[j, :, b])
        correction[3(i - 1) + a, 3(j - 1) + b] =
            4pi / volume * zi * zj / (denominator * sqrt(masses[i] * masses[j]))
    end
    Hermitian(correction)
end

function reference_polar_tensors(fixture::QEFixture)
    entry = fixture_reference(fixture)["phonons"]["0.0,0.0,0.0"]
    dielectric = rows_to_matrix(entry["dielectric"])
    charges = zeros(length(fixture.species), 3, 3)
    for atom in eachindex(fixture.species)
        charges[atom, :, :] = rows_to_matrix(entry["born_charges_asr"][atom])
    end
    charges, dielectric
end

function reference_nac_frequencies(fixture::QEFixture, direction)
    masses = native_masses(fixture)
    analytic = project_acoustic_sum_rule(
        oracle_dynamical(fixture, "0.0,0.0,0.0"), masses)
    charges, dielectric = reference_polar_tensors(fixture)
    correction = nonanalytic_correction(
        charges, dielectric, masses, abs(det(fixture.lattice_bohr)), direction)
    signed_frequencies_cm1(analytic + correction)
end

function shifted_fixture(fixture::QEFixture, atom::Int, cartesian_shift)
    positions = copy(fixture.positions_fractional)
    positions[:, atom] .+= fixture.lattice_bohr \ collect(cartesian_shift)
    QEFixture(
        fixture.name * "-shifted-" * string(hash(cartesian_shift)), fixture.family,
        fixture.lattice_bohr, fixture.species, fixture.masses_amu, positions,
        fixture.ecut_ry, fixture.kgrid, fixture.input_dft, fixture.q_reduced,
    )
end

distorted_si_fixture() = shifted_fixture(si_fixture(), 2, (0.03, -0.02, 0.01))

function strained_fixture(fixture::QEFixture, strain)
    lattice = (I + strain) * fixture.lattice_bohr
    QEFixture(
        fixture.name * "-strained-" * string(hash(strain)), fixture.family,
        lattice, fixture.species, fixture.masses_amu,
        fixture.positions_fractional, fixture.ecut_ry, fixture.kgrid,
        fixture.input_dft, fixture.q_reduced,
    )
end

function private_qe_input(fixture::QEFixture, directory::AbstractString)
    outdir = joinpath(directory, "scratch")
    write_pw_input(joinpath(directory, "pw.in"), fixture,
                   pseudopotential_paths(fixture), outdir)
end


function qe_input_in_cell_units(fixture::QEFixture, directory::AbstractString,
                                unit::Symbol)
    source_path = private_qe_input(fixture, directory)
    text = read(source_path, String)
    unit in (:bohr, :angstrom, :alat) || throw(ArgumentError("unsupported test unit"))
    if unit == :bohr
        return source_path
    end

    scale = unit == :angstrom ? BOHR_TO_ANGSTROM : 1 / norm(fixture.lattice_bohr[:, 1])
    replacement = unit == :angstrom ? "CELL_PARAMETERS angstrom" : "CELL_PARAMETERS alat"
    text = replace(text, "CELL_PARAMETERS bohr" => replacement)
    for vector in eachcol(fixture.lattice_bohr)
        old = join(QuantumDFTVerify.fmt.(vector), " ")
        new = join(QuantumDFTVerify.fmt.(vector .* scale), " ")
        text = replace(text, old => new)
    end
    if unit == :alat
        text = replace(text, "  ibrav = 0\n" =>
                             "  ibrav = 0\n  celldm(1) = $(norm(fixture.lattice_bohr[:, 1]))\n")
    end
    path = joinpath(directory, "pw-$unit.in")
    write(path, text)
    path
end

function caught_argumenterror(f)
    try
        f()
        nothing
    catch error
        error
    end
end
