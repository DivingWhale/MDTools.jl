# MDTools.jl

Julia 分子动力学分析工具包，支持 GROMACS XTC 轨迹文件和 GRO 结构文件读取，以及原子选择功能。

## 安装

```julia
using Pkg
Pkg.develop(path="/path/to/MDTools.jl")
```

## 快速开始

```julia
using MDTools

# 读取拓扑和轨迹
top = read_gro("structure.gro")
traj = read_xtc("trajectory.xtc")

# 创建 Universe
u = Universe(top, traj)

# 选择原子
ca_atoms = select_by_name(top, "CA")
backbone = select_backbone(top)
protein = select_protein(top)

# 提取坐标
coords = get_coords(traj.frames[1], ca_atoms)
```

## XTC 轨迹读取

### 读取完整轨迹

```julia
traj = read_xtc("trajectory.xtc")
println("Frames: ", traj.nframes)
println("Atoms: ", traj.natoms)
```

### 迭代帧（内存高效）

```julia
for frame in eachframe("trajectory.xtc")
    println("Step: ", frame.step, " Time: ", frame.time, " ps")
end
```

## GRO 文件读取

```julia
top = read_gro("structure.gro")
println("Atoms: ", top.natoms)
println("Residues: ", length(residue_ids(top)))
```

## 原子选择

### 基本选择

```julia
# 按原子名
ca = select_by_name(top, "CA")
backbone_atoms = select_by_name(top, ["N", "CA", "C", "O"])

# 按残基名
glycines = select_by_resname(top, "GLY")

# 按残基号
res_range = select_by_resid(top, 10, 50)
```

### 预定义选择器

```julia
select_all(top)        # 所有原子
select_protein(top)    # 蛋白质原子
select_backbone(top)   # 骨架原子
select_sidechain(top)  # 侧链原子
select_water(top)      # 水分子
select_ions(top)       # 离子
select_hydrogen(top)   # 氢原子
select_heavy(top)      # 重原子
```

### 组合选择

```julia
# AND (交集)
ca_in_gly = select_and(select_by_name(top, "CA"), select_by_resname(top, "GLY"))

# OR (并集)
n_or_c = select_or(select_by_name(top, "N"), select_by_name(top, "C"))

# NOT (补集)
non_water = select_not(top, select_water(top))
```

## 性能

使用零分配迭代器 `eachframe()`，读取速度约 **16,000+ FPS**（3726 原子），优于 MDAnalysis。

## 代码结构

```
src/
├── MDTools.jl    # 模块入口
├── xtc.jl        # XTC 轨迹读取
├── topology.jl   # 拓扑和 GRO 读取
└── selection.jl  # 原子选择
```

## License

MIT License
