# 位缓冲区实现，用于读取压缩坐标数据

"""
    BitBuffer

位级别的数据缓冲区，用于从压缩数据中提取任意位数的整数。

# 字段
- `data::Vector{UInt8}`: 原始字节数据
- `index::Int`: 当前读取位置（1-indexed）
- `lastbits::Int`: 上一字节中剩余的有效位数
- `lastbyte::UInt32`: 缓存的字节数据
"""
mutable struct BitBuffer
    data::Vector{UInt8}
    index::Int
    lastbits::Int
    lastbyte::UInt32
end

"""
    BitBuffer(data::Vector{UInt8}) -> BitBuffer

从字节数组创建位缓冲区。
"""
function BitBuffer(data::Vector{UInt8})
    return BitBuffer(data, 1, 0, UInt32(0))
end

"""
    receivebits(buffer::BitBuffer, num_of_bits::Int) -> Int

从位缓冲区读取指定数量的位，返回对应的整数值。
"""
function receivebits(buffer::BitBuffer, num_of_bits::Int)::Int
    mask = (1 << num_of_bits) - 1
    num = 0
    
    while num_of_bits >= 8
        buffer.lastbyte = (buffer.lastbyte << 8) | buffer.data[buffer.index]
        buffer.index += 1
        num |= (buffer.lastbyte >> buffer.lastbits) << (num_of_bits - 8)
        num_of_bits -= 8
    end
    
    if num_of_bits > 0
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

"""
    receiveints!(buffer::BitBuffer, num_of_ints::Int, num_of_bits::Int, 
                 sizes::Vector{<:Integer}, nums::Vector{Int})

从位缓冲区解码多个整数到 nums 数组中。
这是 sendints 的逆操作，用于解压使用混合基数编码的整数组。
"""
function receiveints!(buffer::BitBuffer, num_of_ints::Int, num_of_bits::Int,
                      sizes::Vector{<:Integer}, nums::Vector{Int})
    bytes = zeros(Int, 32)
    num_of_bytes = 0
    
    # 从位缓冲区读取所有字节
    bits_remaining = num_of_bits
    while bits_remaining > 8
        num_of_bytes += 1
        bytes[num_of_bytes] = receivebits(buffer, 8)
        bits_remaining -= 8
    end
    if bits_remaining > 0
        num_of_bytes += 1
        bytes[num_of_bytes] = receivebits(buffer, bits_remaining)
    end
    
    # 从高位整数开始，依次除以对应的 size 来还原各个整数
    for i in num_of_ints:-1:2
        if sizes[i] == 0
            error("Cannot read trajectory, file possibly corrupted (size is zero).")
        end
        num = 0
        for j in num_of_bytes:-1:1
            num = (num << 8) | bytes[j]
            p = num ÷ sizes[i]
            bytes[j] = p
            num = num - p * sizes[i]
        end
        nums[i] = num
    end
    
    # 最后一个整数直接从剩余字节组合
    nums[1] = bytes[1] | (bytes[2] << 8) | (bytes[3] << 16) | (bytes[4] << 24)
end
