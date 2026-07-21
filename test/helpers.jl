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

function shifted_fixture(fixture::QEFixture, atom::Int, cartesian_shift)
    positions = copy(fixture.positions_fractional)
    positions[:, atom] .+= fixture.lattice_bohr \ collect(cartesian_shift)
    QEFixture(
        fixture.name * "-shifted-" * string(hash(cartesian_shift)), fixture.family,
        fixture.lattice_bohr, fixture.species, fixture.masses_amu, positions,
        fixture.ecut_ry, fixture.kgrid, fixture.input_dft, fixture.q_reduced,
    )
end

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
