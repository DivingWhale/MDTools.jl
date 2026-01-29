using Test
using MDTools

@testset "MDTools.jl" begin
    @testset "sizeofint" begin
        @test MDTools.sizeofint(1) == 1
        @test MDTools.sizeofint(2) == 2
        @test MDTools.sizeofint(255) == 8
        @test MDTools.sizeofint(256) == 9
    end

    @testset "BitBuffer" begin
        data = UInt8[0b11010110, 0b10101010]
        buffer = MDTools.BitBuffer(data)
        result = MDTools.receivebits(buffer, 4)
        @test result == 0b1101
        result = MDTools.receivebits(buffer, 4)
        @test result == 0b0110
    end

    @testset "Read XTC" begin
        test_file = joinpath(@__DIR__, "test.xtc")
        if isfile(test_file)
            traj = read_xtc(test_file)
            @test traj.natoms == 3726
            @test traj.nframes == 1001
            @test traj.frames[1].step == 0
            @test traj.frames[1].time ≈ 0.0f0
            @test traj.frames[1].natoms == 3726
            @test traj.frames[1].precision ≈ 1000.0f0
            @test size(traj.frames[1].coords) == (3, 3726)
            @test traj.frames[end].step == 5000000
            @test traj.frames[end].time ≈ 10000.0f0
            @test traj.frames[1].box[1, 1] ≈ 7.4124293f0
            @test traj.frames[1].coords[1, 1] ≈ 4.399f0 atol = 0.01
            @test traj.frames[1].coords[2, 1] ≈ 2.44f0 atol = 0.01
            @test traj.frames[1].coords[3, 1] ≈ 5.126f0 atol = 0.01
        else
            @warn "Test file not found: $test_file"
        end
    end

    @testset "eachframe iterator" begin
        test_file = joinpath(@__DIR__, "test.xtc")
        if isfile(test_file)
            count = 0
            for frame in eachframe(test_file)
                count += 1
                if count >= 10
                    break
                end
            end
            @test count == 10
        end
    end

    @testset "Read GRO" begin
        test_file = joinpath(@__DIR__, "test.gro")
        if isfile(test_file)
            top = read_gro(test_file)
            @test top.natoms == 3726
            @test length(top.atoms) == 3726
            @test length(residue_ids(top)) == 258
            @test top.atoms[1].name == "N"
            @test top.atoms[1].resname == "GLY"
            @test top.atoms[1].resid == 35
            @test length(top.box) == 3
            @test top.box[1] ≈ 4.56125f0 atol = 0.01
        else
            @warn "Test file not found: $test_file"
        end
    end

    @testset "Atom Selection" begin
        test_file = joinpath(@__DIR__, "test.gro")
        if isfile(test_file)
            top = read_gro(test_file)
            @test length(select_by_name(top, "CA")) == 258
            @test length(select_by_resname(top, "GLY")) > 0
            @test length(select_by_resid(top, 35)) > 0

            protein = select_protein(top)
            backbone = select_backbone(top)
            sidechain = select_sidechain(top)
            hydrogen = select_hydrogen(top)
            heavy = select_heavy(top)

            @test length(protein) == 3705
            @test length(backbone) > 0
            @test length(sidechain) > 0
            @test length(hydrogen) + length(heavy) == top.natoms

            ca = select_by_name(top, "CA")
            gly = select_by_resname(top, "GLY")
            ca_in_gly = select_and(ca, gly)
            @test length(ca_in_gly) == 24

            n = select_by_name(top, "N")
            c = select_by_name(top, "C")
            n_or_c = select_or(n, c)
            @test length(n_or_c) == length(n) + length(c)

            not_protein = select_not(top, protein)
            @test length(protein) + length(not_protein) == top.natoms
        end
    end

    @testset "Universe" begin
        xtc_file = joinpath(@__DIR__, "test.xtc")
        gro_file = joinpath(@__DIR__, "test.gro")
        if isfile(xtc_file) && isfile(gro_file)
            top = read_gro(gro_file)
            traj = read_xtc(xtc_file)
            u = Universe(top, traj)
            @test u.topology.natoms == u.trajectory.natoms

            ca = select_by_name(top, "CA")
            coords = get_coords(u.trajectory.frames[1], ca)
            @test size(coords) == (3, length(ca))
        end
    end
end
