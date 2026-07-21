@testset "UT-001 explicit-mass geometry canonicalization" begin
    fixture = si_fixture()
    fractional = candidate_crystal(fixture)
    cartesian = Crystal(fixture.lattice_bohr, fixture.species,
                        fixture.lattice_bohr * fixture.positions_fractional;
                        masses=native_masses(fixture), positions_are_fractional=false)
    @test fractional.lattice ≈ cartesian.lattice
    @test fractional.positions ≈ cartesian.positions
    @test fractional.species == cartesian.species
    @test fractional.masses == cartesian.masses
    lattice_copy = fractional.lattice
    lattice_copy[1, 1] += 1
    @test fractional.lattice ≈ fixture.lattice_bohr
end

@testset "UT-002 lattice translation canonicalization" begin
    fixture = si_fixture()
    translated = copy(fixture.positions_fractional)
    translated[:, 2] .+= [2, -1, 3]
    crystal = Crystal(fixture.lattice_bohr, fixture.species, translated;
                      masses=native_masses(fixture))
    @test crystal.positions ≈ candidate_crystal(fixture).positions
end

@testset "UT-003 invalid crystal rejection" begin
    fixture = si_fixture()
    masses = native_masses(fixture)
    @test_throws ArgumentError Crystal(zeros(3, 3), fixture.species,
                                       fixture.positions_fractional; masses=masses)
    @test_throws ArgumentError Crystal(fixture.lattice_bohr, fixture.species,
                                       fill(NaN, 3, 2); masses=masses)
    @test_throws ArgumentError Crystal(fixture.lattice_bohr, fixture.species,
                                       fixture.positions_fractional; masses=[masses[1], 0.0])
    @test_throws ArgumentError Crystal(fixture.lattice_bohr, [:Si],
                                       fixture.positions_fractional; masses=masses)
end

@testset "UT-004 valid NC-UPF LDA and PBE models" begin
    for fixture in (si_fixture(), nacl_fixture())
        model = candidate_basis(fixture).model
        @test model.xc == (fixture.input_dft == "LDA" ? :lda : :pbe)
        @test model.electron_count > 0
        @test isinteger(model.electron_count)
    end
end

@testset "UT-005 excluded model features fail closed" begin
    fixture = si_fixture()
    crystal = candidate_crystal(fixture)
    pseudos = pseudopotential_paths(fixture)
    @test_throws ArgumentError KSModel(crystal; pseudopotentials=pseudos, xc=:scan)
    @test_throws ArgumentError KSModel(crystal; pseudopotentials=pseudos, charge=1)
    @test_throws ArgumentError KSModel(crystal; pseudopotentials=pseudos, spin=:polarized)
    @test_throws ArgumentError KSModel(crystal; pseudopotentials=Dict{Symbol,String}(), xc=:lda)
    @test_throws ArgumentError KSModel(crystal; pseudopotentials=pseudos, xc=:pbe)

    mixed = PseudoFamily("sssp.mixed.sr.pbe.v1_3_0.efficiency.upf")
    uspp = pseudofile(mixed, :Si)
    @test occursin("USPP", read(uspp, String))
    @test_throws ArgumentError KSModel(crystal;
                                       pseudopotentials=Dict(:Si => uspp), xc=:pbe)

    oxygen_crystal = Crystal(fixture.lattice_bohr, [:O, :O],
                             fixture.positions_fractional;
                             masses=fill(15.999 * QuantumDFTVerify.AMU_TO_ELECTRON_MASS, 2))
    paw = pseudofile(mixed, :O)
    @test occursin("PAW", read(paw, String))
    @test_throws ArgumentError KSModel(oxygen_crystal;
                                       pseudopotentials=Dict(:O => paw), xc=:pbe)
end

@testset "UT-006 exact plane-wave cutoff enumeration" begin
    basis = candidate_basis(si_fixture())
    @test length(basis.G_vectors) == length(basis.kpoints)
    for (k, vectors) in zip(basis.kpoints, basis.G_vectors)
        reciprocal = 2pi .* inv(basis.model.crystal.lattice)'
        radius = sqrt(2basis.Ecut) / minimum(svdvals(reciprocal))
        ranges = ntuple(axis ->
            floor(Int, -k[axis] - radius) - 1:ceil(Int, -k[axis] + radius) + 1, 3)
        expected = Set{NTuple{3,Int}}()
        for g1 in ranges[1], g2 in ranges[2], g3 in ranges[3]
            G = (g1, g2, g3)
            kinetic = sum(abs2, reciprocal * (collect(k) .+ collect(G))) / 2
            kinetic <= basis.Ecut + 100eps(basis.Ecut) && push!(expected, G)
        end
        @test Set(vectors) == expected
    end
end

@testset "UT-007 insufficient FFT grid rejection" begin
    fixture = si_fixture()
    model = candidate_basis(fixture).model
    @test_throws ArgumentError PlaneWaveBasis(model; Ecut=fixture.ecut_ry / 2,
                                               kgrid=fixture.kgrid, fft_size=(2, 2, 2))
end

@testset "UT-008 full k-mesh weights" begin
    basis = candidate_basis(si_fixture())
    @test length(basis.kpoints) == prod(si_fixture().kgrid)
    @test length(basis.kweights) == length(basis.kpoints)
    @test length(unique(basis.kpoints)) == length(basis.kpoints)
    @test all(weight -> weight > 0, basis.kweights)
    @test sum(basis.kweights) ≈ 1 atol=2e-15
end

@testset "UT-009 SCF option defaults and validation" begin
    defaults = SCFOptions()
    @test defaults.energy_tolerance == 1e-10
    @test defaults.density_tolerance == 1e-8
    @test defaults.maxiter == 100
    @test defaults.extra_bands == 4
    @test_throws ArgumentError SCFOptions(maxiter=0)
    @test_throws ArgumentError SCFOptions(density_tolerance=-1)
end

@testset "UT-010 atomic displacement construction" begin
    displacement = AtomicDisplacement(2, 3, (0.25, 0.25, 0.25))
    @test displacement.atom == 2
    @test displacement.direction == 3
    @test displacement.q == (0.25, 0.25, 0.25)
end

@testset "UT-011 displacement and q rejection" begin
    @test_throws ArgumentError AtomicDisplacement(0, 1, (0.0, 0.0, 0.0))
    @test_throws ArgumentError AtomicDisplacement(1, 4, (0.0, 0.0, 0.0))
    @test_throws ArgumentError AtomicDisplacement(1, 1, (NaN, 0.0, 0.0))
    gs = candidate_state(si_fixture())
    @test_throws ArgumentError response(gs, AtomicDisplacement(3, 1, (0.0, 0.0, 0.0)))
    @test_throws ArgumentError dynamical_matrix(gs, (0.2, 0.0, 0.0))
end
