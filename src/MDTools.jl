module MDTools

# XTC 轨迹相关
export read_xtc, XTCFrame, XTCTrajectory, eachframe

# 拓扑相关
export Atom, Topology, Universe
export read_gro

# 选择器
export select_by_name, select_by_resname, select_by_resid, select_by_index
export select_and, select_or, select_not
export select_all, select_protein, select_backbone, select_sidechain
export select_water, select_ions, select_hydrogen, select_heavy
export get_atoms, get_coords
export residue_ids, residue_names, atom_names

# 加载模块
include("xtc.jl")
include("topology.jl")
include("selection.jl")

end # module
