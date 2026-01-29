# XTCReader.jl

A Julia package for reading GROMACS XTC trajectory files and GRO structure files, with atom selection capabilities.

## Installation

```julia
using Pkg
Pkg.develop(path="/path/to/XTCReader.jl")
```

## Quick Start

```julia
using XTCReader

# Read topology and trajectory
top = read_gro("structure.gro")
traj = read_xtc("trajectory.xtc")

# Create a Universe (combines topology and trajectory)
u = Universe(top, traj)

# Select atoms
ca_atoms = select_by_name(top, "CA")
backbone = select_backbone(top)
protein = select_protein(top)

# Extract coordinates for selected atoms
coords = get_coords(traj.frames[1], ca_atoms)
```

## Reading XTC Trajectories

### Read entire trajectory

```julia
using XTCReader

traj = read_xtc("trajectory.xtc")

println("Number of atoms: ", traj.natoms)
println("Number of frames: ", traj.nframes)

# Access first frame
frame = traj.frames[1]
println("Step: ", frame.step)
println("Time: ", frame.time, " ps")
println("Box: ", frame.box)
println("Coordinates shape: ", size(frame.coords))  # (3, natoms)
```

### Iterate over frames (memory efficient)

```julia
for frame in eachframe("trajectory.xtc")
    println("Step: ", frame.step, " Time: ", frame.time, " ps")
    # Process frame.coords...
end
```

## Reading GRO Files

```julia
top = read_gro("structure.gro")

println("Number of atoms: ", top.natoms)
println("Number of residues: ", length(residue_ids(top)))
println("First atom: ", top.atoms[1])
```

## Atom Selection

### Basic Selection

```julia
# By atom name
ca = select_by_name(top, "CA")
backbone_atoms = select_by_name(top, ["N", "CA", "C", "O"])

# By residue name
glycines = select_by_resname(top, "GLY")
aromatics = select_by_resname(top, ["PHE", "TYR", "TRP"])

# By residue ID
res1 = select_by_resid(top, 1)
res_range = select_by_resid(top, 10, 50)  # residues 10-50

# By index range
first_100 = select_by_index(top, 1, 100)
```

### Predefined Selectors

```julia
select_all(top)        # All atoms
select_protein(top)    # Protein atoms
select_backbone(top)   # Backbone (N, CA, C, O)
select_sidechain(top)  # Sidechain atoms
select_water(top)      # Water molecules
select_ions(top)       # Ions
select_hydrogen(top)   # Hydrogen atoms
select_heavy(top)      # Heavy (non-hydrogen) atoms
```

### Combining Selections

```julia
# AND (intersection)
ca_in_gly = select_and(select_by_name(top, "CA"), select_by_resname(top, "GLY"))

# OR (union)
n_or_c = select_or(select_by_name(top, "N"), select_by_name(top, "C"))

# NOT (complement)
non_water = select_not(top, select_water(top))
```

### Extracting Coordinates

```julia
# Get coordinates from a frame
ca_atoms = select_by_name(top, "CA")
coords = get_coords(frame, ca_atoms)  # Returns 3×n matrix

# From Universe
coords = get_coords(u, frame_idx, ca_atoms)

# Get atom information
atoms = get_atoms(top, ca_atoms)  # Returns Vector{Atom}
```

## Data Structures

### Atom

```julia
struct Atom
    index::Int       # Atom index (1-based)
    name::String     # Atom name ("CA", "N", etc.)
    resname::String  # Residue name ("GLY", "ALA", etc.)
    resid::Int       # Residue number
    x::Float32       # Initial x coordinate (nm)
    y::Float32       # Initial y coordinate (nm)
    z::Float32       # Initial z coordinate (nm)
end
```

### Topology

```julia
struct Topology
    natoms::Int
    atoms::Vector{Atom}
    title::String
    box::Vector{Float32}  # Box dimensions [x, y, z]
    # ... lookup indices
end
```

### XTCFrame

```julia
struct XTCFrame
    step::Int64           # Frame number (simulation step)
    time::Float32         # Simulation time (ps)
    box::Matrix{Float32}  # 3×3 box matrix (nm)
    natoms::Int32         # Number of atoms
    precision::Float32    # Coordinate precision
    coords::Matrix{Float32}  # 3×natoms coordinate matrix (nm)
end
```

### Universe

```julia
struct Universe
    topology::Topology
    trajectory::XTCTrajectory
end
```

## Utility Functions

```julia
residue_ids(top)    # Get sorted list of residue IDs
residue_names(top)  # Get list of residue names
atom_names(top)     # Get list of atom names
```

## Technical Details

The XTC format uses:
- XDR (External Data Representation) for cross-platform binary serialization
- Lossy coordinate compression with configurable precision
- Run-length encoding for similar consecutive coordinates
- Adaptive bit-width encoding based on coordinate ranges

This implementation is based on the GROMACS source code analysis and has been validated against Chemfiles.jl.

## License

MIT License
