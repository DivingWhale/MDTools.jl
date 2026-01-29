# GRO 文件读取器

"""
    read_gro(filename::String) -> Topology

读取 GROMACS GRO 格式的结构文件，返回拓扑对象。

# GRO 格式
- 第1行: 标题
- 第2行: 原子数
- 第3行到第 n+2 行: 原子数据
  - 列 1-5: 残基号
  - 列 6-10: 残基名 (5字符)
  - 列 11-15: 原子名 (5字符)
  - 列 16-20: 原子号
  - 列 21-28, 29-36, 37-44: x, y, z 坐标 (nm, f8.3 格式)
  - 列 45-52, 53-60, 61-68: 可选速度
- 最后一行: 盒子矢量

# Example
```julia
top = read_gro("structure.gro")
println("Atoms: ", top.natoms)
println("First atom: ", top.atoms[1])
```
"""
function read_gro(filename::String)::Topology
    atoms = Atom[]
    title = ""
    box = Float32[0.0, 0.0, 0.0]
    natoms = 0

    open(filename, "r") do io
        # 读取标题
        title = strip(readline(io))

        # 读取原子数
        natoms_line = strip(readline(io))
        natoms = parse(Int, natoms_line)

        # 预分配数组
        sizehint!(atoms, natoms)

        # 读取原子数据
        for i in 1:natoms
            line = readline(io)
            atom = parse_gro_atom(line, i)
            push!(atoms, atom)
        end

        # 读取盒子向量
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

"""
    parse_gro_atom(line::String, expected_index::Int) -> Atom

解析 GRO 文件中的单行原子数据。
"""
function parse_gro_atom(line::String, expected_index::Int)::Atom
    # GRO 格式使用固定宽度列
    # 列 1-5: 残基号 (右对齐)
    # 列 6-10: 残基名 (左对齐, 但实际可能没有填满)
    # 列 11-15: 原子名 (右对齐)
    # 列 16-20: 原子号 (右对齐)
    # 列 21-28: x (f8.3)
    # 列 29-36: y (f8.3)
    # 列 37-44: z (f8.3)

    # 解析残基号和残基名 (前10个字符)
    resid_resname = line[1:min(10, length(line))]

    # 残基号在前5个字符
    resid_str = strip(resid_resname[1:5])
    resid = parse(Int, resid_str)

    # 残基名在6-10位置
    resname = strip(resid_resname[6:end])

    # 原子名在11-15位置
    atomname = strip(line[11:min(15, length(line))])

    # 原子号在16-20位置
    atomnum_str = strip(line[16:min(20, length(line))])
    atomnum = parse(Int, atomnum_str)

    # 坐标在21-44位置
    x = parse(Float32, strip(line[21:min(28, length(line))]))
    y = parse(Float32, strip(line[29:min(36, length(line))]))
    z = parse(Float32, strip(line[37:min(44, length(line))]))

    return Atom(expected_index, atomname, resname, resid, x, y, z)
end
