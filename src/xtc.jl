# XTC trajectory file reading module
# Merged from: types.jl, xdr.jl, bitbuffer.jl, compression.jl, reader.jl

# =============================================================================
# 常量定义
# =============================================================================

# XTC magic number constants
const XTC_MAGIC = Int32(1995)
const XTC_NEW_MAGIC = Int32(2023)

# Equivalent to ~300M atom limit in GROMACS code
const XTC_1995_MAX_NATOMS = 298261617

# Magic integers table for adaptive bit-width encoding
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

# Julia 1-indexed, corresponding to FIRSTIDX = 9 in C code
const FIRSTIDX = 10
const LASTIDX = length(MAGICINTS)

# =============================================================================
# 数据类型
# =============================================================================

"""
    XTCFrame

Represents a single frame in an XTC trajectory file.

# Fields
- `step::Int64`: Frame number (simulation step)
- `time::Float32`: Simulation time (ps)
- `box::Matrix{Float32}`: 3×3 box matrix (nm)
- `natoms::Int32`: Number of atoms
- `precision::Float32`: Coordinate precision
- `coords::Matrix{Float32}`: 3×natoms coordinate matrix (nm)
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

Mutable version of XTCFrame for iterator reuse to avoid memory allocation.
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
        Int64(0), Float32(0), Matrix{Float32}(undef, 3, 3),
        Int32(natoms), Float32(0), Matrix{Float32}(undef, 3, Int(natoms))
    )
end

"""
    XTCTrajectory

Represents a complete XTC trajectory file.

# Fields
- `filename::String`: File path
- `natoms::Int32`: Number of atoms
- `nframes::Int`: Number of frames
- `frames::Vector{XTCFrame}`: All frames
"""
struct XTCTrajectory
    filename::String
    natoms::Int32
    nframes::Int
    frames::Vector{XTCFrame}
end

"""
    XTCBuffer

Pre-allocated buffer to avoid repeated memory allocation during decompression.
"""
mutable struct XTCBuffer
    compressed::Vector{UInt8}
    compressed_capacity::Int
    minint::Vector{Int32}
    maxint::Vector{Int32}
    sizeint::Vector{Int}
    bitsizeint::Vector{Int}
    sizesmall::Vector{Int}
    thiscoord::Vector{Int}
    prevcoord::Vector{Int}
    ri_bytes::Vector{Int}
end

function XTCBuffer(max_compressed_size::Int=1024 * 1024)
    return XTCBuffer(
        Vector{UInt8}(undef, max_compressed_size), max_compressed_size,
        zeros(Int32, 3), zeros(Int32, 3), zeros(Int, 3), zeros(Int, 3),
        zeros(Int, 3), zeros(Int, 3), zeros(Int, 3), zeros(Int, 32)
    )
end

# =============================================================================
# XDR 读取函数
# =============================================================================

@inline function read_xdr_int(io::IO)::Int32
    return ntoh(read(io, Int32))
end

@inline function read_xdr_int64(io::IO)::Int64
    return ntoh(read(io, Int64))
end

@inline function read_xdr_float(io::IO)::Float32
    bytes = read(io, 4)
    return reinterpret(Float32, ntoh(reinterpret(UInt32, bytes)[1]))
end

# =============================================================================
# 位缓冲区
# =============================================================================

mutable struct BitBuffer{T<:AbstractVector{UInt8}}
    data::T
    index::Int
    lastbits::Int
    lastbyte::UInt32
end

@inline function BitBuffer(data::T) where {T<:AbstractVector{UInt8}}
    return BitBuffer{T}(data, 1, 0, UInt32(0))
end

@inline function receivebits(buffer::BitBuffer, num_of_bits::Int)::Int
    mask = (1 << num_of_bits) - 1
    num = 0

    @inbounds while num_of_bits >= 8
        buffer.lastbyte = (buffer.lastbyte << 8) | buffer.data[buffer.index]
        buffer.index += 1
        num |= (buffer.lastbyte >> buffer.lastbits) << (num_of_bits - 8)
        num_of_bits -= 8
    end

    @inbounds if num_of_bits > 0
        if buffer.lastbits < num_of_bits
            buffer.lastbits += 8
            buffer.lastbyte = (buffer.lastbyte << 8) | buffer.data[buffer.index]
            buffer.index += 1
        end
        buffer.lastbits -= num_of_bits
        num |= (buffer.lastbyte >> buffer.lastbits) & ((1 << num_of_bits) - 1)
    end

    return num & mask
end

@inline function receiveints!(buffer::BitBuffer, num_of_ints::Int, num_of_bits::Int,
    sizes::Vector{<:Integer}, nums::Vector{Int}, bytes::Vector{Int})
    num_of_bytes = 0
    bits_remaining = num_of_bits

    @inbounds while bits_remaining > 8
        num_of_bytes += 1
        bytes[num_of_bytes] = receivebits(buffer, 8)
        bits_remaining -= 8
    end
    @inbounds if bits_remaining > 0
        num_of_bytes += 1
        bytes[num_of_bytes] = receivebits(buffer, bits_remaining)
    end

    @inbounds for i in num_of_ints:-1:2
        num = 0
        for j in num_of_bytes:-1:1
            num = (num << 8) | bytes[j]
            p = num ÷ sizes[i]
            bytes[j] = p
            num = num - p * sizes[i]
        end
        nums[i] = num
    end

    @inbounds nums[1] = bytes[1] | (bytes[2] << 8) | (bytes[3] << 16) | (bytes[4] << 24)
end

# =============================================================================
# 压缩算法
# =============================================================================

@inline function sizeofint(size::Integer)::Int
    num = 1
    num_of_bits = 0
    while size >= num && num_of_bits < 32
        num_of_bits += 1
        num <<= 1
    end
    return num_of_bits
end

function sizeofints(num_of_ints::Int, sizes::Vector{<:Integer})::Int
    bytes = zeros(UInt8, 32)
    bytes[1] = 1
    num_of_bytes = 1

    @inbounds for i in 1:num_of_ints
        tmp::UInt64 = 0
        for j in 1:num_of_bytes
            tmp = UInt64(bytes[j]) * sizes[i] + tmp
            bytes[j] = tmp & 0xff
            tmp >>= 8
        end
        while tmp != 0
            num_of_bytes += 1
            bytes[num_of_bytes] = tmp & 0xff
            tmp >>= 8
        end
    end

    num = 1
    num_of_bits = 0
    @inbounds while bytes[num_of_bytes] >= num
        num_of_bits += 1
        num *= 2
    end

    return num_of_bits + (num_of_bytes - 1) * 8
end

function decompress_coords!(io::IO, natoms::Int32, magic::Int32,
    coords::Matrix{Float32}, buf::XTCBuffer)::Float32
    lsize = read_xdr_int(io)
    if lsize != natoms
        @warn "Atom count mismatch: expected $natoms, got $lsize"
    end

    # Small systems don't use compression
    if lsize <= 9
        @inbounds for i in 1:lsize
            for j in 1:3
                coords[j, i] = read_xdr_float(io)
            end
        end
        return Float32(-1)
    end

    precision = read_xdr_float(io)

    @inbounds for i in 1:3
        buf.minint[i] = read_xdr_int(io)
    end
    @inbounds for i in 1:3
        buf.maxint[i] = read_xdr_int(io)
    end
    @inbounds for i in 1:3
        buf.sizeint[i] = buf.maxint[i] - buf.minint[i] + 1
    end

    local bitsize::Int
    if any(s -> s > 0xffffff, buf.sizeint)
        @inbounds for i in 1:3
            buf.bitsizeint[i] = sizeofint(buf.sizeint[i])
        end
        bitsize = 0
    else
        @inbounds for i in 1:3
            buf.bitsizeint[i] = 0
        end
        bitsize = sizeofints(3, buf.sizeint)
    end

    smallidx = Int(read_xdr_int(io))

    local bufsize::Int64
    if magic == XTC_NEW_MAGIC
        bufsize = read_xdr_int64(io)
    else
        bufsize = Int64(read_xdr_int(io))
    end

    if bufsize > buf.compressed_capacity
        resize!(buf.compressed, Int(bufsize))
        buf.compressed_capacity = Int(bufsize)
    end

    read!(io, @view buf.compressed[1:Int(bufsize)])
    padding = (4 - bufsize % 4) % 4
    if padding > 0
        skip(io, padding)
    end

    buffer = BitBuffer(@view buf.compressed[1:Int(bufsize)])

    smaller = smallidx > (FIRSTIDX - 1) ? MAGICINTS[smallidx] ÷ 2 : 0
    smallnum = MAGICINTS[smallidx+1] ÷ 2
    @inbounds for i in 1:3
        buf.sizesmall[i] = Int(MAGICINTS[smallidx+1])
    end

    inv_precision = 1.0f0 / precision
    thiscoord = buf.thiscoord
    prevcoord = buf.prevcoord
    ri_bytes = buf.ri_bytes

    run = 0
    i = 1
    coord_idx = 1

    @inbounds while i <= lsize
        if bitsize == 0
            for k in 1:3
                thiscoord[k] = receivebits(buffer, buf.bitsizeint[k])
            end
        else
            receiveints!(buffer, 3, bitsize, buf.sizeint, thiscoord, ri_bytes)
        end

        for k in 1:3
            thiscoord[k] += buf.minint[k]
        end

        prevcoord[1] = thiscoord[1]
        prevcoord[2] = thiscoord[2]
        prevcoord[3] = thiscoord[3]

        flag = receivebits(buffer, 1)
        is_smaller = 0
        if flag == 1
            run = receivebits(buffer, 5)
            is_smaller = run % 3
            run -= is_smaller
            is_smaller -= 1
        end

        if run > 0
            for k in 0:3:(run-1)
                receiveints!(buffer, 3, smallidx, buf.sizesmall, thiscoord, ri_bytes)
                i += 1

                for j in 1:3
                    thiscoord[j] += prevcoord[j] - smallnum
                end

                if k == 0
                    for j in 1:3
                        thiscoord[j], prevcoord[j] = prevcoord[j], thiscoord[j]
                    end
                    for j in 1:3
                        coords[j, coord_idx] = prevcoord[j] * inv_precision
                    end
                    coord_idx += 1
                else
                    prevcoord[1] = thiscoord[1]
                    prevcoord[2] = thiscoord[2]
                    prevcoord[3] = thiscoord[3]
                end

                for j in 1:3
                    coords[j, coord_idx] = thiscoord[j] * inv_precision
                end
                coord_idx += 1
            end
        else
            for j in 1:3
                coords[j, coord_idx] = thiscoord[j] * inv_precision
            end
            coord_idx += 1
        end

        smallidx += is_smaller
        if is_smaller < 0
            smallnum = smaller
            if smallidx > (FIRSTIDX - 1)
                smaller = MAGICINTS[smallidx] ÷ 2
            else
                smaller = 0
            end
        elseif is_smaller > 0
            smaller = smallnum
            smallnum = MAGICINTS[smallidx+1] ÷ 2
        end
        for j in 1:3
            buf.sizesmall[j] = MAGICINTS[smallidx+1]
        end

        i += 1
    end

    return precision
end

function decompress_coords!(io::IO, natoms::Int32, magic::Int32)::Tuple{Float32,Matrix{Float32}}
    coords = Matrix{Float32}(undef, 3, Int(natoms))
    buf = XTCBuffer()
    precision = decompress_coords!(io, natoms, magic, coords, buf)
    return precision, coords
end

# =============================================================================
# 读取接口
# =============================================================================

function read_xtc_header(io::IO)
    if eof(io)
        return nothing
    end
    try
        magic = read_xdr_int(io)
        if magic != XTC_MAGIC && magic != XTC_NEW_MAGIC
            error("Invalid XTC magic number: $magic")
        end
        natoms = read_xdr_int(io)
        step = Int64(read_xdr_int(io))
        time = read_xdr_float(io)
        return (magic, natoms, step, time)
    catch e
        if isa(e, EOFError)
            return nothing
        end
        rethrow(e)
    end
end

@inline function read_xtc_box!(io::IO, box::Matrix{Float32})
    @inbounds for i in 1:3
        for j in 1:3
            box[i, j] = read_xdr_float(io)
        end
    end
end

function read_xtc_box(io::IO)::Matrix{Float32}
    box = Matrix{Float32}(undef, 3, 3)
    read_xtc_box!(io, box)
    return box
end

function read_xtc_frame!(io::IO, frame::MutableXTCFrame, buf::XTCBuffer)::Bool
    header = read_xtc_header(io)
    if header === nothing
        return false
    end
    magic, natoms, step, time = header
    frame.step = step
    frame.time = time
    frame.natoms = natoms
    read_xtc_box!(io, frame.box)
    frame.precision = decompress_coords!(io, natoms, magic, frame.coords, buf)
    return true
end

function read_xtc_frame(io::IO)::Union{XTCFrame,Nothing}
    header = read_xtc_header(io)
    if header === nothing
        return nothing
    end
    magic, natoms, step, time = header
    box = read_xtc_box(io)
    precision, coords = decompress_coords!(io, natoms, magic)
    return XTCFrame(step, time, box, natoms, precision, coords)
end

"""
    read_xtc(filename::String) -> XTCTrajectory

Reads a complete XTC trajectory file.

# Example
```julia
traj = read_xtc("trajectory.xtc")
println("Frames: ", traj.nframes)
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

Returns an iterator to read XTC file frame by frame. Reuses frame objects to avoid memory allocation.

# Example
```julia
for frame in eachframe("trajectory.xtc")
    println("Step: ", frame.step)
end
```
"""
function eachframe(filename::String)
    return XTCFrameIterator(filename)
end

struct XTCFrameIterator
    filename::String
end

mutable struct XTCIteratorState
    io::IOStream
    frame::MutableXTCFrame
    buf::XTCBuffer
end

function Base.iterate(iter::XTCFrameIterator)
    io = open(iter.filename, "r")
    header = read_xtc_header(io)
    if header === nothing
        close(io)
        return nothing
    end
    magic, natoms, step, time = header
    frame = MutableXTCFrame(natoms)
    buf = XTCBuffer()
    frame.step = step
    frame.time = time
    frame.natoms = natoms
    read_xtc_box!(io, frame.box)
    frame.precision = decompress_coords!(io, natoms, magic, frame.coords, buf)
    state = XTCIteratorState(io, frame, buf)
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
