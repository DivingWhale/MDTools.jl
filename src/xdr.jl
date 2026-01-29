# XDR (External Data Representation) 序列化读取函数
# XDR 使用大端序 (Big Endian)

"""
    read_xdr_int(io::IO) -> Int32

从 IO 流读取一个 XDR 格式的 32 位有符号整数。
"""
function read_xdr_int(io::IO)::Int32
    return ntoh(read(io, Int32))
end

"""
    read_xdr_int64(io::IO) -> Int64

从 IO 流读取一个 XDR 格式的 64 位有符号整数。
"""
function read_xdr_int64(io::IO)::Int64
    return ntoh(read(io, Int64))
end

"""
    read_xdr_float(io::IO) -> Float32

从 IO 流读取一个 XDR 格式的 32 位浮点数。
"""
function read_xdr_float(io::IO)::Float32
    # IEEE 754 浮点数的字节序也需要转换
    bytes = read(io, 4)
    return reinterpret(Float32, ntoh(reinterpret(UInt32, bytes)[1]))
end

"""
    read_xdr_opaque(io::IO, len::Integer) -> Vector{UInt8}

从 IO 流读取指定长度的不透明数据块。
XDR 要求数据块按 4 字节对齐，但这里我们只读取实际数据。
"""
function read_xdr_opaque(io::IO, len::Integer)::Vector{UInt8}
    data = read(io, len)
    # XDR 要求 4 字节对齐，跳过填充字节
    padding = (4 - len % 4) % 4
    if padding > 0
        skip(io, padding)
    end
    return data
end

"""
    read_xdr_floats(io::IO, n::Integer) -> Vector{Float32}

从 IO 流读取 n 个 XDR 格式的浮点数。
"""
function read_xdr_floats(io::IO, n::Integer)::Vector{Float32}
    result = Vector{Float32}(undef, n)
    for i in 1:n
        result[i] = read_xdr_float(io)
    end
    return result
end
