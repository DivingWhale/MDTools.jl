# Topology data structure and GRO file reading
# Merged from: topology.jl, gro.jl

# =============================================================================
# 数据类型
# =============================================================================

"""
    Atom

Represents topological information for a single atom.

# Fields
- `index::Int`: Atom index (1-based)
- `name::String`: Atom name, e.g., "CA", "N", "O"
- `resname::String`: Residue name, e.g., "GLY", "ALA"
- `resid::Int`: Residue number
- `x::Float32`: Initial x coordinate (nm)
- `y::Float32`: Initial y coordinate (nm)
- `z::Float32`: Initial z coordinate (nm)
"""
struct Atom
    index::Int
    name::String
    resname::String
    resid::Int
    x::Float32
    y::Float32
    z::Float32
end

"""
    Topology

Molecular topology structure containing information for all atoms and fast lookup indices.

# Fields
- `natoms::Int`: Total number of atoms
- `atoms::Vector{Atom}`: List of atoms
- `title::String`: Topology title
- `box::Vector{Float32}`: Box dimensions [x, y, z]
- `atom_names::Dict{String, Vector{Int}}`: Atom name -> list of atom indices
- `residue_names::Dict{String, Vector{Int}}`: Residue name -> list of atom indices
- `residue_ids::Dict{Int, Vector{Int}}`: Residue number -> list of atom indices
"""
struct Topology
    natoms::Int
    atoms::Vector{Atom}
    title::String
    box::Vector{Float32}
    atom_names::Dict{String,Vector{Int}}
    residue_names::Dict{String,Vector{Int}}
    residue_ids::Dict{Int,Vector{Int}}
end

function Topology(atoms::Vector{Atom}, title::AbstractString, box::Vector{Float32})
    natoms = length(atoms)
    atom_names = Dict{String,Vector{Int}}()
    residue_names = Dict{String,Vector{Int}}()
    residue_ids = Dict{Int,Vector{Int}}()

    for atom in atoms
        if !haskey(atom_names, atom.name)
            atom_names[atom.name] = Int[]
        end
        push!(atom_names[atom.name], atom.index)

        if !haskey(residue_names, atom.resname)
            residue_names[atom.resname] = Int[]
        end
        push!(residue_names[atom.resname], atom.index)

        if !haskey(residue_ids, atom.resid)
            residue_ids[atom.resid] = Int[]
        end
        push!(residue_ids[atom.resid], atom.index)
    end

    return Topology(natoms, atoms, String(title), box, atom_names, residue_names, residue_ids)
end

"""
    Universe

Analysis unit combining topology and trajectory.

# Fields
- `topology::Topology`: Topology structure
- `trajectory::XTCTrajectory`: Trajectory data
"""
struct Universe
    topology::Topology
    trajectory::XTCTrajectory
end

# =============================================================================
# 显示方法
# =============================================================================

function Base.show(io::IO, atom::Atom)
    print(io, "Atom($(atom.index), $(atom.name), $(atom.resname)$(atom.resid))")
end

function Base.show(io::IO, top::Topology)
    print(io, "Topology($(top.natoms) atoms, $(length(top.residue_ids)) residues)")
end

function Base.show(io::IO, u::Universe)
    print(io, "Universe($(u.topology.natoms) atoms, $(u.trajectory.nframes) frames)")
end

# =============================================================================
# GRO 文件读取
# =============================================================================

"""
    read_gro(filename::String) -> Topology

Reads a GROMACS GRO format structure file and returns a topology object.

# GRO Format
- Line 1: Title
- Line 2: Number of atoms
- Atom lines: residue number (5), residue name (5), atom name (5), atom number (5), x, y, z
- Last line: Box vector

# Example
```julia
top = read_gro("structure.gro")
println("Atoms: ", top.natoms)
```
"""
function read_gro(filename::String)::Topology
    atoms = Atom[]
    title = ""
    box = Float32[0.0, 0.0, 0.0]
    natoms = 0

    open(filename, "r") do io
        title = strip(readline(io))
        natoms_line = strip(readline(io))
        natoms = parse(Int, natoms_line)
        sizehint!(atoms, natoms)

        for i in 1:natoms
            line = readline(io)
            atom = parse_gro_atom(line, i)
            push!(atoms, atom)
        end

        box_line = strip(readline(io))
        box_values = split(box_line)
        if length(box_values) >= 3
            box[1] = parse(Float32, box_values[1])
            box[2] = parse(Float32, box_values[2])
            box[3] = parse(Float32, box_values[3])
        end
    end

    return Topology(atoms, title, box)
end

function parse_gro_atom(line::String, expected_index::Int)::Atom
    resid_str = strip(line[1:5])
    resid = parse(Int, resid_str)
    resname = strip(line[6:min(10, length(line))])
    atomname = strip(line[11:min(15, length(line))])
    x = parse(Float32, strip(line[21:min(28, length(line))]))
    y = parse(Float32, strip(line[29:min(36, length(line))]))
    z = parse(Float32, strip(line[37:min(44, length(line))]))
    return Atom(expected_index, atomname, resname, resid, x, y, z)
end
