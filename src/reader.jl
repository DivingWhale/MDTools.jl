# XTC 文件读取接口

"""
    read_xtc_header(io::IO) -> Union{Tuple{Int32, Int32, Int64, Float32}, Nothing}

读取 XTC 帧头。返回 (magic, natoms, step, time) 元组，或在文件末尾返回 nothing。
"""
function read_xtc_header(io::IO)
    if eof(io)
        return nothing
    end

    try
        magic = read_xdr_int(io)

        # 验证魔数
        if magic != XTC_MAGIC && magic != XTC_NEW_MAGIC
            error("Invalid XTC magic number: $magic (expected $XTC_MAGIC or $XTC_NEW_MAGIC)")
        end

        natoms = read_xdr_int(io)
        step = Int64(read_xdr_int(io))  # XTC 使用 32 位 step，但我们存储为 64 位
        time = read_xdr_float(io)

        return (magic, natoms, step, time)
    catch e
        if isa(e, EOFError)
            return nothing
        end
        rethrow(e)
    end
end

"""
    read_xtc_box(io::IO) -> Matrix{Float32}

读取 3×3 盒子矩阵。
"""
function read_xtc_box(io::IO)::Matrix{Float32}
    box = Matrix{Float32}(undef, 3, 3)
    for i in 1:3
        for j in 1:3
            box[i, j] = read_xdr_float(io)
        end
    end
    return box
end

"""
    read_xtc_frame(io::IO) -> Union{XTCFrame, Nothing}

从 IO 流读取单个 XTC 帧。到达文件末尾返回 nothing。
"""
function read_xtc_frame(io::IO)::Union{XTCFrame,Nothing}
    header = read_xtc_header(io)
    if header === nothing
        return nothing
    end

    magic, natoms, step, time = header

    # 读取盒子矩阵
    box = read_xtc_box(io)

    # 读取并解压缩坐标
    precision, coords = decompress_coords!(io, natoms, magic)

    return XTCFrame(step, time, box, natoms, precision, coords)
end

"""
    read_xtc(filename::String) -> XTCTrajectory

读取完整的 XTC 轨迹文件。

# 示例
```julia
traj = read_xtc("trajectory.xtc")
println("Number of frames: ", traj.nframes)
println("Number of atoms: ", traj.natoms)
```
"""
function read_xtc(filename::String)::XTCTrajectory
    frames = XTCFrame[]
    natoms = Int32(0)

    open(filename, "r") do io
        while true
            frame = read_xtc_frame(io)
            if frame === nothing
                break
            end

            if isempty(frames)
                natoms = frame.natoms
            end

            push!(frames, frame)
        end
    end

    return XTCTrajectory(filename, natoms, length(frames), frames)
end

"""
    eachframe(filename::String)

返回一个迭代器，逐帧读取 XTC 文件。适合处理大型轨迹文件。

# 示例
```julia
for frame in eachframe("trajectory.xtc")
    println("Step: ", frame.step, " Time: ", frame.time)
end
```
"""
function eachframe(filename::String)
    return XTCFrameIterator(filename)
end

struct XTCFrameIterator
    filename::String
end

function Base.iterate(iter::XTCFrameIterator)
    io = open(iter.filename, "r")
    return iterate(iter, io)
end

function Base.iterate(iter::XTCFrameIterator, io::IO)
    frame = read_xtc_frame(io)
    if frame === nothing
        close(io)
        return nothing
    end
    return (frame, io)
end

Base.IteratorSize(::Type{XTCFrameIterator}) = Base.SizeUnknown()
Base.eltype(::Type{XTCFrameIterator}) = XTCFrame
