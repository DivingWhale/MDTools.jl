# 原子选择器

# 常见氨基酸残基名
const PROTEIN_RESNAMES = Set([
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # 带电/修饰形式
    "HIE", "HID", "HIP", "CYX", "ASH", "GLH", "LYN",
    # N/C 端修饰
    "NALA", "NARG", "NASN", "NASP", "NCYS", "NGLN", "NGLU", "NGLY", "NHIS", "NILE",
    "NLEU", "NLYS", "NMET", "NPHE", "NPRO", "NSER", "NTHR", "NTRP", "NTYR", "NVAL",
    "CALA", "CARG", "CASN", "CASP", "CCYS", "CGLN", "CGLU", "CGLY", "CHIS", "CILE",
    "CLEU", "CLYS", "CMET", "CPHE", "CPRO", "CSER", "CTHR", "CTRP", "CTYR", "CVAL"
])

# 常见水分子残基名
const WATER_RESNAMES = Set(["SOL", "WAT", "HOH", "TIP3", "TIP4", "TIP5", "SPC", "SPCE", "OPC"])

# 常见离子残基名
const ION_RESNAMES = Set(["NA", "CL", "K", "MG", "CA", "ZN", "FE", "NA+", "CL-", "K+", "MG2+", "CA2+"])

# 骨架原子名
const BACKBONE_NAMES = Set(["N", "CA", "C", "O", "H", "HA", "OC1", "OC2", "OXT"])

# 氢原子名（通常以 H 开头）
function is_hydrogen(name::String)::Bool
    return length(name) > 0 && name[1] == 'H'
end

# ============================================================================
# 基础选择函数
# ============================================================================

"""
    select_by_name(top::Topology, names::Vector{String}) -> Vector{Int}
    select_by_name(top::Topology, name::String) -> Vector{Int}

按原子名选择原子，返回原子索引列表。

# Example
```julia
ca_atoms = select_by_name(top, "CA")
backbone = select_by_name(top, ["N", "CA", "C", "O"])
```
"""
function select_by_name(top::Topology, names::Vector{String})::Vector{Int}
    indices = Int[]
    for name in names
        if haskey(top.atom_names, name)
            append!(indices, top.atom_names[name])
        end
    end
    return sort!(unique!(indices))
end

select_by_name(top::Topology, name::String) = select_by_name(top, [name])

"""
    select_by_resname(top::Topology, resnames::Vector{String}) -> Vector{Int}
    select_by_resname(top::Topology, resname::String) -> Vector{Int}

按残基名选择原子，返回原子索引列表。

# Example
```julia
glycines = select_by_resname(top, "GLY")
aromatic = select_by_resname(top, ["PHE", "TYR", "TRP"])
```
"""
function select_by_resname(top::Topology, resnames::Vector{String})::Vector{Int}
    indices = Int[]
    for resname in resnames
        if haskey(top.residue_names, resname)
            append!(indices, top.residue_names[resname])
        end
    end
    return sort!(unique!(indices))
end

select_by_resname(top::Topology, resname::String) = select_by_resname(top, [resname])

"""
    select_by_resid(top::Topology, resids::Vector{Int}) -> Vector{Int}
    select_by_resid(top::Topology, resid::Int) -> Vector{Int}
    select_by_resid(top::Topology, start::Int, stop::Int) -> Vector{Int}

按残基号选择原子，返回原子索引列表。

# Example
```julia
res1 = select_by_resid(top, 1)
res_range = select_by_resid(top, 10, 20)  # 残基 10-20
```
"""
function select_by_resid(top::Topology, resids::Vector{Int})::Vector{Int}
    indices = Int[]
    for resid in resids
        if haskey(top.residue_ids, resid)
            append!(indices, top.residue_ids[resid])
        end
    end
    return sort!(unique!(indices))
end

select_by_resid(top::Topology, resid::Int) = select_by_resid(top, [resid])

function select_by_resid(top::Topology, start::Int, stop::Int)::Vector{Int}
    return select_by_resid(top, collect(start:stop))
end

"""
    select_by_index(top::Topology, start::Int, stop::Int) -> Vector{Int}

按原子索引范围选择原子。

# Example
```julia
first_100 = select_by_index(top, 1, 100)
```
"""
function select_by_index(top::Topology, start::Int, stop::Int)::Vector{Int}
    return collect(max(1, start):min(top.natoms, stop))
end

# ============================================================================
# 逻辑组合选择
# ============================================================================

"""
    select_and(indices1::Vector{Int}, indices2::Vector{Int}) -> Vector{Int}

返回两个选择的交集。

# Example
```julia
ca_in_gly = select_and(select_by_name(top, "CA"), select_by_resname(top, "GLY"))
```
"""
function select_and(indices1::Vector{Int}, indices2::Vector{Int})::Vector{Int}
    return sort!(collect(intersect(Set(indices1), Set(indices2))))
end

"""
    select_or(indices1::Vector{Int}, indices2::Vector{Int}) -> Vector{Int}

返回两个选择的并集。

# Example
```julia
n_or_c = select_or(select_by_name(top, "N"), select_by_name(top, "C"))
```
"""
function select_or(indices1::Vector{Int}, indices2::Vector{Int})::Vector{Int}
    return sort!(unique!(vcat(indices1, indices2)))
end

"""
    select_not(top::Topology, indices::Vector{Int}) -> Vector{Int}

返回选择的补集（不在给定索引中的所有原子）。

# Example
```julia
non_water = select_not(top, select_water(top))
```
"""
function select_not(top::Topology, indices::Vector{Int})::Vector{Int}
    all_indices = Set(1:top.natoms)
    return sort!(collect(setdiff(all_indices, Set(indices))))
end

# ============================================================================
# 预定义选择器
# ============================================================================

"""
    select_all(top::Topology) -> Vector{Int}

选择所有原子。
"""
function select_all(top::Topology)::Vector{Int}
    return collect(1:top.natoms)
end

"""
    select_protein(top::Topology) -> Vector{Int}

选择蛋白质原子（基于常见氨基酸残基名）。
"""
function select_protein(top::Topology)::Vector{Int}
    indices = Int[]
    for (resname, atom_indices) in top.residue_names
        if resname in PROTEIN_RESNAMES
            append!(indices, atom_indices)
        end
    end
    return sort!(indices)
end

"""
    select_backbone(top::Topology) -> Vector{Int}

选择蛋白质骨架原子 (N, CA, C, O)。
"""
function select_backbone(top::Topology)::Vector{Int}
    protein_indices = Set(select_protein(top))
    backbone_indices = Int[]

    for name in ["N", "CA", "C", "O"]
        if haskey(top.atom_names, name)
            for idx in top.atom_names[name]
                if idx in protein_indices
                    push!(backbone_indices, idx)
                end
            end
        end
    end

    return sort!(backbone_indices)
end

"""
    select_sidechain(top::Topology) -> Vector{Int}

选择蛋白质侧链原子（不包含骨架原子）。
"""
function select_sidechain(top::Topology)::Vector{Int}
    protein = select_protein(top)
    backbone = select_backbone(top)
    return select_and(protein, select_not(top, backbone))
end

"""
    select_water(top::Topology) -> Vector{Int}

选择水分子原子。
"""
function select_water(top::Topology)::Vector{Int}
    indices = Int[]
    for (resname, atom_indices) in top.residue_names
        if resname in WATER_RESNAMES
            append!(indices, atom_indices)
        end
    end
    return sort!(indices)
end

"""
    select_ions(top::Topology) -> Vector{Int}

选择离子原子。
"""
function select_ions(top::Topology)::Vector{Int}
    indices = Int[]
    for (resname, atom_indices) in top.residue_names
        if resname in ION_RESNAMES
            append!(indices, atom_indices)
        end
    end
    return sort!(indices)
end

"""
    select_hydrogen(top::Topology) -> Vector{Int}

选择氢原子（原子名以 H 开头）。
"""
function select_hydrogen(top::Topology)::Vector{Int}
    indices = Int[]
    for (name, atom_indices) in top.atom_names
        if is_hydrogen(name)
            append!(indices, atom_indices)
        end
    end
    return sort!(indices)
end

"""
    select_heavy(top::Topology) -> Vector{Int}

选择重原子（非氢原子）。
"""
function select_heavy(top::Topology)::Vector{Int}
    return select_not(top, select_hydrogen(top))
end

# ============================================================================
# 辅助函数
# ============================================================================

"""
    get_atoms(top::Topology, indices::Vector{Int}) -> Vector{Atom}

根据索引获取原子列表。
"""
function get_atoms(top::Topology, indices::Vector{Int})::Vector{Atom}
    return [top.atoms[i] for i in indices if 1 <= i <= top.natoms]
end

"""
    get_coords(frame::XTCFrame, indices::Vector{Int}) -> Matrix{Float32}

从帧中提取选定原子的坐标。

返回 3 x n 矩阵，其中 n 是选定原子数。
"""
function get_coords(frame::XTCFrame, indices::Vector{Int})::Matrix{Float32}
    return frame.coords[:, indices]
end

"""
    get_coords(u::Universe, frame_idx::Int, indices::Vector{Int}) -> Matrix{Float32}

从 Universe 的指定帧中提取选定原子的坐标。
"""
function get_coords(u::Universe, frame_idx::Int, indices::Vector{Int})::Matrix{Float32}
    return u.trajectory.frames[frame_idx].coords[:, indices]
end

"""
    residue_ids(top::Topology) -> Vector{Int}

获取所有残基号的有序列表。
"""
function residue_ids(top::Topology)::Vector{Int}
    return sort!(collect(keys(top.residue_ids)))
end

"""
    residue_names(top::Topology) -> Vector{String}

获取所有残基名的列表。
"""
function residue_names(top::Topology)::Vector{String}
    return collect(keys(top.residue_names))
end

"""
    atom_names(top::Topology) -> Vector{String}

获取所有原子名的列表。
"""
function atom_names(top::Topology)::Vector{String}
    return collect(keys(top.atom_names))
end
