@testset "UT-012 QE input type and case-insensitive parser" begin
    mktempdir() do directory
        path = private_qe_input(si_fixture(), directory)
        original = read(path, String)
        mixed = replace(original, "&CONTROL" => "&control",
                         "ATOMIC_POSITIONS" => "atomic_positions")
        mixed_path = joinpath(directory, "mixed.in")
        write(mixed_path, mixed)
        a = read_qe_input(path)
        b = read_qe_input(IOBuffer(mixed))
        c = read_qe_input(mixed_path)
        @test a isa QEInput
        @test a.model.crystal.positions ≈ b.model.crystal.positions
        @test b.basis.Ecut == c.basis.Ecut
    end
end

@testset "UT-013 QE unit conversion" begin
    mktempdir() do directory
        parsed = read_qe_input(private_qe_input(si_fixture(), directory))
        fixture = si_fixture()
        @test parsed.model.crystal.lattice ≈ fixture.lattice_bohr
        @test parsed.model.crystal.masses ≈ native_masses(fixture)
        @test parsed.basis.Ecut == fixture.ecut_ry / 2
    end
end

@testset "UT-014 QE pseudo_dir resolution and hash" begin
    mktempdir() do directory
        parsed = read_qe_input(private_qe_input(si_fixture(), directory))
        @test parsed.model.crystal.species == [:Si, :Si]
        @test pseudopotential_hashes(si_fixture())["Si"] ==
              fixture_reference(si_fixture())["pseudopotential_sha256"]["Si"]
    end
end

@testset "UT-015 unsupported QE fields fail closed" begin
    mktempdir() do directory
        input = read(private_qe_input(si_fixture(), directory), String)
        for invalid in (
            replace(input, "calculation = 'scf'" => "calculation = 'relax'"),
            replace(input, "ibrav = 0" => "ibrav = 2"),
            replace(input, "occupations = 'fixed'" => "occupations = 'smearing'"),
            replace(input, "nspin = 1" => "nspin = 2"),
        )
            @test_throws ArgumentError read_qe_input(IOBuffer(invalid))
        end
    end
end
