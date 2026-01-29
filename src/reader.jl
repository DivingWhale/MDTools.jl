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
    read_xtc_box!(io::IO, box::Matrix{Float32})

读取 3×3 盒子矩阵到预分配的矩阵中。
"""
@inline function read_xtc_box!(io::IO, box::Matrix{Float32})
    @inbounds for i in 1:3
        for j in 1:3
            box[i, j] = read_xdr_float(io)
        end
    end
end

"""
    read_xtc_box(io::IO) -> Matrix{Float32}

读取 3×3 盒子矩阵。
"""
function read_xtc_box(io::IO)::Matrix{Float32}
    box = Matrix{Float32}(undef, 3, 3)
    read_xtc_box!(io, box)
    return box
end

"""
    read_xtc_frame!(io::IO, frame::MutableXTCFrame, buf::XTCBuffer) -> Bool

从 IO 流读取单个 XTC 帧到预分配的 MutableXTCFrame 中。
返回 true 表示成功，false 表示到达文件末尾。
"""
function read_xtc_frame!(io::IO, frame::MutableXTCFrame, buf::XTCBuffer)::Bool
    header = read_xtc_header(io)
    if header === nothing
        return false
    end

    magic, natoms, step, time = header

    frame.step = step
    frame.time = time
    frame.natoms = natoms

    # 读取盒子矩阵
    read_xtc_box!(io, frame.box)

    # 读取并解压缩坐标到预分配矩阵
    frame.precision = decompress_coords!(io, natoms, magic, frame.coords, buf)

    return true
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
此迭代器复用同一个帧对象以避免内存分配。

# 示例
```julia
for frame in eachframe("trajectory.xtc")
    println("Step: ", frame.step, " Time: ", frame.time)
end
```

# 注意
迭代器返回的帧对象在每次迭代时会被修改，如需保留数据请复制。
"""
function eachframe(filename::String)
    return XTCFrameIterator(filename)
end

struct XTCFrameIterator
    filename::String
end

# 迭代器状态：包含 IO、预分配的 Frame 和 Buffer
mutable struct XTCIteratorState
    io::IOStream
    frame::MutableXTCFrame
    buf::XTCBuffer
    initialized::Bool
end

function Base.iterate(iter::XTCFrameIterator)
    io = open(iter.filename, "r")

    # 首先读取第一帧的头部来获取原子数
    header = read_xtc_header(io)
    if header === nothing
        close(io)
        return nothing
    end

    magic, natoms, step, time = header

    # 创建预分配的帧和缓冲区
    frame = MutableXTCFrame(natoms)
    buf = XTCBuffer()

    # 填充第一帧数据
    frame.step = step
    frame.time = time
    frame.natoms = natoms
    read_xtc_box!(io, frame.box)
    frame.precision = decompress_coords!(io, natoms, magic, frame.coords, buf)

    state = XTCIteratorState(io, frame, buf, true)
    return (frame, state)
end

function Base.iterate(iter::XTCFrameIterator, state::XTCIteratorState)
    if !read_xtc_frame!(state.io, state.frame, state.buf)
        close(state.io)
        return nothing
    end
    return (state.frame, state)
end

Base.IteratorSize(::Type{XTCFrameIterator}) = Base.SizeUnknown()
Base.eltype(::Type{XTCFrameIterator}) = MutableXTCFrame
