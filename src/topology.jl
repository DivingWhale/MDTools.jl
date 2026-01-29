# 拓扑数据结构

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

"""
    Topology(atoms::Vector{Atom}, title::String, box::Vector{Float32}) -> Topology

从原子列表构建拓扑结构，自动创建查找索引。
"""
function Topology(atoms::Vector{Atom}, title::AbstractString, box::Vector{Float32})
    natoms = length(atoms)

    # 构建查找索引
    atom_names = Dict{String,Vector{Int}}()
    residue_names = Dict{String,Vector{Int}}()
    residue_ids = Dict{Int,Vector{Int}}()

    for atom in atoms
        # 原子名索引
        if !haskey(atom_names, atom.name)
            atom_names[atom.name] = Int[]
        end
        push!(atom_names[atom.name], atom.index)

        # 残基名索引
        if !haskey(residue_names, atom.resname)
            residue_names[atom.resname] = Int[]
        end
        push!(residue_names[atom.resname], atom.index)

        # 残基号索引
        if !haskey(residue_ids, atom.resid)
            residue_ids[atom.resid] = Int[]
        end
        push!(residue_ids[atom.resid], atom.index)
    end

    return Topology(natoms, atoms, String(title), box, atom_names, residue_names, residue_ids)
end

"""
    Universe

结合拓扑和轨迹的分析单元，类似于 MDAnalysis 的 Universe。

# 字段
- `topology::Topology`: 拓扑结构
- `trajectory::XTCTrajectory`: 轨迹数据
"""
struct Universe
    topology::Topology
    trajectory::XTCTrajectory
end

# 显示方法
function Base.show(io::IO, atom::Atom)
    print(io, "Atom($(atom.index), $(atom.name), $(atom.resname)$(atom.resid))")
end

function Base.show(io::IO, top::Topology)
    print(io, "Topology($(top.natoms) atoms, $(length(top.residue_ids)) residues)")
end

function Base.show(io::IO, u::Universe)
    print(io, "Universe($(u.topology.natoms) atoms, $(u.trajectory.nframes) frames)")
end
