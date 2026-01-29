module XTCReader

# XTC 相关导出
export read_xtc, XTCFrame, XTCTrajectory, eachframe

# 拓扑相关导出
export Atom, Topology, Universe
export read_gro

# 选择器导出
export select_by_name, select_by_resname, select_by_resid, select_by_index
export select_and, select_or, select_not
export select_all, select_protein, select_backbone, select_sidechain
export select_water, select_ions, select_hydrogen, select_heavy
export get_atoms, get_coords
export residue_ids, residue_names, atom_names

# XTC 读取模块
include("types.jl")
include("xdr.jl")
include("bitbuffer.jl")
include("compression.jl")
include("reader.jl")

# 拓扑和选择模块
include("topology.jl")
include("gro.jl")
include("selection.jl")

end # module

