# XTC 坐标压缩/解压缩算法实现

"""
    sizeofint(size::Integer) -> Int

计算存储给定大小的整数所需的位数。
"""
@inline function sizeofint(size::Integer)::Int
    num = 1
    num_of_bits = 0
    while size >= num && num_of_bits < 32
        num_of_bits += 1
        num <<= 1
    end
    return num_of_bits
end

"""
    sizeofints(num_of_ints::Int, sizes::Vector{<:Integer}) -> Int

计算使用混合基数编码存储多个整数所需的总位数。
"""
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

"""
    decompress_coords!(io::IO, natoms::Int32, magic::Int32, 
                       coords::Matrix{Float32}, buf::XTCBuffer) -> Float32

解压缩 XTC 坐标数据到预分配的 coords 矩阵中。
返回精度值。
"""
function decompress_coords!(io::IO, natoms::Int32, magic::Int32,
    coords::Matrix{Float32}, buf::XTCBuffer)::Float32
    lsize = read_xdr_int(io)

    if lsize != natoms
        @warn "Atom count mismatch: expected $natoms, got $lsize"
    end

    size3 = Int(lsize) * 3

    # 小系统不使用压缩
    if lsize <= 9
        @inbounds for i in 1:lsize
            for j in 1:3
                coords[j, i] = read_xdr_float(io)
            end
        end
        return Float32(-1)
    end

    # 读取精度
    precision = read_xdr_float(io)

    # 读取最小/最大整数值（复用缓冲区）
    @inbounds for i in 1:3
        buf.minint[i] = read_xdr_int(io)
    end
    @inbounds for i in 1:3
        buf.maxint[i] = read_xdr_int(io)
    end

    # 计算尺寸
    @inbounds for i in 1:3
        buf.sizeint[i] = buf.maxint[i] - buf.minint[i] + 1
    end

    # 判断是否使用独立位编码
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

    # 读取 smallidx
    smallidx = Int(read_xdr_int(io))

    # 读取压缩数据大小
    local bufsize::Int64
    if magic == XTC_NEW_MAGIC
        bufsize = read_xdr_int64(io)
    else
        bufsize = Int64(read_xdr_int(io))
    end

    # 确保缓冲区足够大
    if bufsize > buf.compressed_capacity
        resize!(buf.compressed, Int(bufsize))
        buf.compressed_capacity = Int(bufsize)
    end

    # 读取压缩数据到预分配缓冲区
    read!(io, @view buf.compressed[1:Int(bufsize)])
    # XDR 4字节对齐
    padding = (4 - bufsize % 4) % 4
    if padding > 0
        skip(io, padding)
    end

    # 创建 BitBuffer（复用数据）
    buffer = BitBuffer(@view buf.compressed[1:Int(bufsize)])

    # 初始化解压缩参数
    smaller = smallidx > (FIRSTIDX - 1) ? MAGICINTS[smallidx] ÷ 2 : 0
    smallnum = MAGICINTS[smallidx+1] ÷ 2
    @inbounds for i in 1:3
        buf.sizesmall[i] = Int(MAGICINTS[smallidx+1])
    end

    inv_precision = 1.0f0 / precision

    # 使用预分配的临时缓冲区
    thiscoord = buf.thiscoord
    prevcoord = buf.prevcoord
    ri_bytes = buf.ri_bytes

    run = 0
    i = 1
    coord_idx = 1

    @inbounds while i <= lsize
        # 读取基础坐标
        if bitsize == 0
            for k in 1:3
                thiscoord[k] = receivebits(buffer, buf.bitsizeint[k])
            end
        else
            receiveints!(buffer, 3, bitsize, buf.sizeint, thiscoord, ri_bytes)
        end

        # 加上最小值偏移
        for k in 1:3
            thiscoord[k] += buf.minint[k]
        end

        prevcoord[1] = thiscoord[1]
        prevcoord[2] = thiscoord[2]
        prevcoord[3] = thiscoord[3]

        # 读取 run-length 控制位
        flag = receivebits(buffer, 1)
        is_smaller = 0
        if flag == 1
            run = receivebits(buffer, 5)
            is_smaller = run % 3
            run -= is_smaller
            is_smaller -= 1
        end

        if run > 0
            # 处理 run-length 编码的后续原子
            for k in 0:3:(run-1)
                receiveints!(buffer, 3, smallidx, buf.sizesmall, thiscoord, ri_bytes)
                i += 1

                for j in 1:3
                    thiscoord[j] += prevcoord[j] - smallnum
                end

                if k == 0
                    # 水分子优化：交换
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

        # 更新 smallidx
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

# 保留旧版本接口以保持兼容性（会分配内存）
function decompress_coords!(io::IO, natoms::Int32, magic::Int32)::Tuple{Float32,Matrix{Float32}}
    coords = Matrix{Float32}(undef, 3, Int(natoms))
    buf = XTCBuffer()
    precision = decompress_coords!(io, natoms, magic, coords, buf)
    return precision, coords
end
