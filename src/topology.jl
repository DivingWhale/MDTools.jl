# 拓扑数据结构和 GRO 文件读取
# 合并自: topology.jl, gro.jl

# =============================================================================
# 数据类型
# =============================================================================

"""
    Atom

表示单个原子的拓扑信息。

# 字段
- `index::Int`: 原子索引 (1-based)
- `name::String`: 原子名称，如 "CA", "N", "O"
- `resname::String`: 残基名称，如 "GLY", "ALA"
- `resid::Int`: 残基编号
- `x::Float32`: 初始 x 坐标 (nm)
- `y::Float32`: 初始 y 坐标 (nm)
- `z::Float32`: 初始 z 坐标 (nm)
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

分子拓扑结构，包含所有原子的信息和快速查找索引。

# 字段
- `natoms::Int`: 原子总数
- `atoms::Vector{Atom}`: 原子列表
- `title::String`: 拓扑标题
- `box::Vector{Float32}`: 盒子尺寸 [x, y, z]
- `atom_names::Dict{String, Vector{Int}}`: 原子名 -> 原子索引列表
- `residue_names::Dict{String, Vector{Int}}`: 残基名 -> 原子索引列表
- `residue_ids::Dict{Int, Vector{Int}}`: 残基号 -> 原子索引列表
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

结合拓扑和轨迹的分析单元。

# 字段
- `topology::Topology`: 拓扑结构
- `trajectory::XTCTrajectory`: 轨迹数据
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

读取 GROMACS GRO 格式的结构文件，返回拓扑对象。

# GRO 格式
- 第1行: 标题
- 第2行: 原子数
- 原子行: 残基号(5), 残基名(5), 原子名(5), 原子号(5), x, y, z
- 最后一行: 盒子矢量

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
