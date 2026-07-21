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
        fixture = si_fixture()
        parsed = Dict(unit => read_qe_input(qe_input_in_cell_units(
            fixture, joinpath(directory, String(unit)), unit))
            for unit in (:bohr, :angstrom, :alat))
        for input in values(parsed)
            @test input.model.crystal.lattice ≈ fixture.lattice_bohr atol=2e-10
            @test input.model.crystal.positions ≈ fixture.positions_fractional atol=2e-12
            @test input.model.crystal.masses ≈ native_masses(fixture) atol=2e-8
            @test input.basis.Ecut == fixture.ecut_ry / 2
        end
    end
end

@testset "UT-014 QE pseudo_dir resolution and hash" begin
    mktempdir() do directory
        parsed = read_qe_input(private_qe_input(si_fixture(), directory))
        @test parsed.model.crystal.species == [:Si, :Si]
        @test pseudopotential_hashes(si_fixture())["Si"] ==
              fixture_reference(si_fixture())["pseudopotential_sha256"]["Si"]
        input = read(private_qe_input(si_fixture(), directory), String)
        missing = replace(input, r"pseudo_dir = '[^']+'" =>
                                "pseudo_dir = 'definitely-missing-directory'")
        error = caught_argumenterror(() -> read_qe_input(IOBuffer(missing)))
        @test error isa ArgumentError
        @test occursin("pseudo", lowercase(sprint(showerror, error)))
    end
end

@testset "UT-015 unsupported QE fields fail closed" begin
    mktempdir() do directory
        input = read(private_qe_input(si_fixture(), directory), String)
        invalid_cases = (
            (replace(input, "calculation = 'scf'" => "calculation = 'relax'"), "calculation"),
            (replace(input, "ibrav = 0" => "ibrav = 2"), "ibrav"),
            (replace(input, "occupations = 'fixed'" => "occupations = 'smearing'"), "occupations"),
            (replace(input, "nspin = 1" => "nspin = 2"), "nspin"),
            (replace(input, "  nat = 2" => "  nat = 3"), "nat"),
            (replace(input, "  ecutwfc" => "  mystery_keyword = 7\n  ecutwfc"), "mystery_keyword"),
            (input * "\nATOMIC_FORCES\nSi 0 0 0\nSi 0 0 0\n", "atomic_forces"),
            (replace(input, "Si 0.25 0.25 0.25" => "Ge 0.25 0.25 0.25"), "ge"),
        )
        for (invalid, field) in invalid_cases
            error = caught_argumenterror(() -> read_qe_input(IOBuffer(invalid)))
            @test error isa ArgumentError
            @test occursin(field, lowercase(sprint(showerror, error)))
        end
    end
end
