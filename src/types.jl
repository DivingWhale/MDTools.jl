# XTC 文件格式相关类型定义

"""
    XTCFrame

表示 XTC 轨迹文件中的单个帧。

# 字段
- `step::Int64`: 帧编号（模拟步数）
- `time::Float32`: 模拟时间（单位：ps）
- `box::Matrix{Float32}`: 3×3 盒子矩阵（单位：nm）
- `natoms::Int32`: 原子数量
- `precision::Float32`: 坐标精度
- `coords::Matrix{Float32}`: 3×natoms 坐标矩阵（单位：nm）
"""
struct XTCFrame
    step::Int64
    time::Float32
    box::Matrix{Float32}
    natoms::Int32
    precision::Float32
    coords::Matrix{Float32}
end

"""
    MutableXTCFrame

可变版本的 XTCFrame，用于迭代器复用以避免内存分配。
"""
mutable struct MutableXTCFrame
    step::Int64
    time::Float32
    box::Matrix{Float32}
    natoms::Int32
    precision::Float32
    coords::Matrix{Float32}
end

function MutableXTCFrame(natoms::Integer)
    return MutableXTCFrame(
        Int64(0),
        Float32(0),
        Matrix{Float32}(undef, 3, 3),
        Int32(natoms),
        Float32(0),
        Matrix{Float32}(undef, 3, Int(natoms))
    )
end

"""
    XTCTrajectory

表示完整的 XTC 轨迹文件。

# 字段
- `filename::String`: 文件路径
- `natoms::Int32`: 原子数量
- `nframes::Int`: 帧数
- `frames::Vector{XTCFrame}`: 所有帧
"""
struct XTCTrajectory
    filename::String
    natoms::Int32
    nframes::Int
    frames::Vector{XTCFrame}
end

"""
    XTCBuffer

预分配的缓冲区，用于解压缩过程中避免重复内存分配。
"""
mutable struct XTCBuffer
    # 压缩数据缓冲区
    compressed::Vector{UInt8}
    compressed_capacity::Int

    # 临时整数数组
    minint::Vector{Int32}
    maxint::Vector{Int32}
    sizeint::Vector{Int}
    bitsizeint::Vector{Int}
    sizesmall::Vector{Int}

    # 坐标临时缓冲区
    thiscoord::Vector{Int}
    prevcoord::Vector{Int}

    # receiveints 专用缓冲区
    ri_bytes::Vector{Int}
end

function XTCBuffer(max_compressed_size::Int=1024 * 1024)
    return XTCBuffer(
        Vector{UInt8}(undef, max_compressed_size),
        max_compressed_size,
        zeros(Int32, 3),
        zeros(Int32, 3),
        zeros(Int, 3),
        zeros(Int, 3),
        zeros(Int, 3),
        zeros(Int, 3),
        zeros(Int, 3),
        zeros(Int, 32)
    )
end

# XTC 魔数常量
const XTC_MAGIC = Int32(1995)
const XTC_NEW_MAGIC = Int32(2023)

# 相当于 GROMACS 代码中约 ~300M 原子的限制
const XTC_1995_MAX_NATOMS = 298261617

# Magic integers 表，用于自适应位宽编码
const MAGICINTS = Int32[
    0, 0, 0, 0, 0, 0, 0, 0, 0, 8,
    10, 12, 16, 20, 25, 32, 40, 50, 64, 80,
    101, 128, 161, 203, 256, 322, 406, 512, 645, 812,
    1024, 1290, 1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192,
    10321, 13003, 16384, 20642, 26007, 32768, 41285, 52015, 65536, 82570,
    104031, 131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561, 832255,
    1048576, 1321122, 1664510, 2097152, 2642245, 3329021, 4194304, 5284491, 6658042, 8388607,
    10568983, 13316085, 16777216
]

# Julia 1-indexed，对应 C 代码的 FIRSTIDX = 9
const FIRSTIDX = 10
const LASTIDX = length(MAGICINTS)
