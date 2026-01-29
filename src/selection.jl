# Atom selectors

# Common amino acid residue names
const PROTEIN_RESNAMES = Set([
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # Charged/modified forms
    "HIE", "HID", "HIP", "CYX", "ASH", "GLH", "LYN",
    # N/C terminal modifications
    "NALA", "NARG", "NASN", "NASP", "NCYS", "NGLN", "NGLU", "NGLY", "NHIS", "NILE",
    "NLEU", "NLYS", "NMET", "NPHE", "NPRO", "NSER", "NTHR", "NTRP", "NTYR", "NVAL",
    "CALA", "CARG", "CASN", "CASP", "CCYS", "CGLN", "CGLU", "CGLY", "CHIS", "CILE",
    "CLEU", "CLYS", "CMET", "CPHE", "CPRO", "CSER", "CTHR", "CTRP", "CTYR", "CVAL"
])

# Common water molecule residue names
const WATER_RESNAMES = Set(["SOL", "WAT", "HOH", "TIP3", "TIP4", "TIP5", "SPC", "SPCE", "OPC"])

# Common ion residue names
const ION_RESNAMES = Set(["NA", "CL", "K", "MG", "CA", "ZN", "FE", "NA+", "CL-", "K+", "MG2+", "CA2+"])

# Backbone atom names
const BACKBONE_NAMES = Set(["N", "CA", "C", "O", "H", "HA", "OC1", "OC2", "OXT"])

# Hydrogen atom names (usually start with H)
function is_hydrogen(name::String)::Bool
    return length(name) > 0 && name[1] == 'H'
end

# ============================================================================
# 基础选择函数
# ============================================================================

"""
    select_by_name(top::Topology, names::Vector{String}) -> Vector{Int}
    select_by_name(top::Topology, name::String) -> Vector{Int}

Select atoms by atom name, returns a list of atom indices.

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

Select atoms by residue name, returns a list of atom indices.

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

Select atoms by residue number, returns a list of atom indices.

# Example
```julia
res1 = select_by_resid(top, 1)
res_range = select_by_resid(top, 10, 20)  # residues 10-20
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

Select atoms by atom index range.

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

Returns the intersection of two selections.

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

Returns the union of two selections.

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

Returns the complement of a selection (all atoms not in the given indices).

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

Selects all atoms.
"""
function select_all(top::Topology)::Vector{Int}
    return collect(1:top.natoms)
end

"""
    select_protein(top::Topology) -> Vector{Int}

Selects protein atoms (based on common amino acid residue names).
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

Selects protein backbone atoms (N, CA, C, O).
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

Selects protein sidechain atoms (excluding backbone atoms).
"""
function select_sidechain(top::Topology)::Vector{Int}
    protein = select_protein(top)
    backbone = select_backbone(top)
    return select_and(protein, select_not(top, backbone))
end

"""
    select_water(top::Topology) -> Vector{Int}

Selects water molecule atoms.
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

Selects ion atoms.
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

Selects hydrogen atoms (atom names starting with H).
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

Selects heavy atoms (non-hydrogen atoms).
"""
function select_heavy(top::Topology)::Vector{Int}
    return select_not(top, select_hydrogen(top))
end

# ============================================================================
# 辅助函数
# ============================================================================

"""
    get_atoms(top::Topology, indices::Vector{Int}) -> Vector{Atom}

Gets a list of atoms by indices.
"""
function get_atoms(top::Topology, indices::Vector{Int})::Vector{Atom}
    return [top.atoms[i] for i in indices if 1 <= i <= top.natoms]
end

"""
    get_coords(frame::XTCFrame, indices::Vector{Int}) -> Matrix{Float32}

Extracts coordinates of selected atoms from a frame.

Returns a 3 x n matrix where n is the number of selected atoms.
"""
function get_coords(frame::XTCFrame, indices::Vector{Int})::Matrix{Float32}
    return frame.coords[:, indices]
end

"""
    get_coords(u::Universe, frame_idx::Int, indices::Vector{Int}) -> Matrix{Float32}

Extracts coordinates of selected atoms from a specified frame in a Universe.
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

Gets a list of all residue names.
"""
function residue_names(top::Topology)::Vector{String}
    return collect(keys(top.residue_names))
end

"""
    atom_names(top::Topology) -> Vector{String}

Gets a list of all atom names.
"""
function atom_names(top::Topology)::Vector{String}
    return collect(keys(top.atom_names))
end
