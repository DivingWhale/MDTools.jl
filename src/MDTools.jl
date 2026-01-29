module MDTools

# XTC trajectory related
export read_xtc, XTCFrame, XTCTrajectory, eachframe

# Topology related
export Atom, Topology, Universe
export read_gro

# Selectors
export select_by_name, select_by_resname, select_by_resid, select_by_index
export select_and, select_or, select_not
export select_all, select_protein, select_backbone, select_sidechain
export select_water, select_ions, select_hydrogen, select_heavy
export get_atoms, get_coords
export residue_ids, residue_names, atom_names

# Load modules
include("xtc.jl")
include("topology.jl")
include("selection.jl")

end # module
