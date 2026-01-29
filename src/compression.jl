# XTC 坐标压缩/解压缩算法实现

"""
    sizeofint(size::Integer) -> Int

计算存储给定大小的整数所需的位数。
"""
function sizeofint(size::Integer)::Int
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

    for i in 1:num_of_ints
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
    while bytes[num_of_bytes] >= num
        num_of_bits += 1
        num *= 2
    end

    return num_of_bits + (num_of_bytes - 1) * 8
end

"""
    decompress_coords!(io::IO, natoms::Int32, magic::Int32) -> Tuple{Float32, Matrix{Float32}}

解压缩 XTC 坐标数据。

返回 (precision, coords) 元组，其中 coords 是 3×natoms 的矩阵。
"""
function decompress_coords!(io::IO, natoms::Int32, magic::Int32)::Tuple{Float32,Matrix{Float32}}
    lsize = read_xdr_int(io)

    if lsize != natoms
        @warn "Atom count mismatch: expected $natoms, got $lsize"
    end

    size3 = Int(lsize) * 3

    # 小系统不使用压缩
    if lsize <= 9
        coords_flat = read_xdr_floats(io, size3)
        coords = reshape(coords_flat, 3, Int(lsize))
        return Float32(-1), coords
    end

    # 读取精度
    precision = read_xdr_float(io)

    # 读取最小/最大整数值
    minint = [read_xdr_int(io) for _ in 1:3]
    maxint = [read_xdr_int(io) for _ in 1:3]

    # 计算尺寸
    sizeint = [maxint[i] - minint[i] + 1 for i in 1:3]

    # 判断是否使用独立位编码
    if any(s > 0xffffff for s in sizeint)
        bitsizeint = [sizeofint(s) for s in sizeint]
        bitsize = 0
    else
        bitsizeint = zeros(Int, 3)
        bitsize = sizeofints(3, sizeint)
    end

    # 读取 smallidx - 保持原始值，因为它直接用作位数
    # 在索引 MAGICINTS 时需要 +1 转换为 Julia 1-indexed
    smallidx = Int(read_xdr_int(io))

    # 读取压缩数据大小
    if magic == XTC_NEW_MAGIC
        bufsize = read_xdr_int64(io)
    else
        bufsize = Int64(read_xdr_int(io))
    end

    # 读取压缩数据
    compressed_data = read_xdr_opaque(io, Int(bufsize))
    buffer = BitBuffer(compressed_data)

    # 初始化解压缩参数
    # 注意：smallidx 用作位数，smallidx+1 用于索引 MAGICINTS（Julia 1-indexed）
    smaller = smallidx > (FIRSTIDX - 1) ? MAGICINTS[smallidx] ÷ 2 : 0
    smallnum = MAGICINTS[smallidx+1] ÷ 2
    sizesmall = fill(Int(MAGICINTS[smallidx+1]), 3)

    inv_precision = 1.0f0 / precision

    # 输出坐标
    coords = Matrix{Float32}(undef, 3, Int(lsize))

    # 临时缓冲区
    thiscoord = zeros(Int, 3)
    prevcoord = zeros(Int, 3)

    run = 0
    i = 1
    coord_idx = 1

    while i <= lsize
        # 读取基础坐标
        if bitsize == 0
            for k in 1:3
                thiscoord[k] = receivebits(buffer, bitsizeint[k])
            end
        else
            receiveints!(buffer, 3, bitsize, sizeint, thiscoord)
        end

        # 加上最小值偏移
        for k in 1:3
            thiscoord[k] += minint[k]
        end

        prevcoord .= thiscoord

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
                receiveints!(buffer, 3, smallidx, sizesmall, thiscoord)
                i += 1

                for j in 1:3
                    thiscoord[j] += prevcoord[j] - smallnum
                end

                if k == 0
                    # 水分子优化：交换第一个后续原子和前一个原子的值
                    # prevcoord 包含主循环读取的坐标，thiscoord 包含刚读取的小坐标
                    # 交换它们的值
                    for j in 1:3
                        thiscoord[j], prevcoord[j] = prevcoord[j], thiscoord[j]
                    end
                    # 输出交换后的 prevcoord（原来的 thiscoord，即小坐标）
                    for j in 1:3
                        coords[j, coord_idx] = prevcoord[j] * inv_precision
                    end
                    coord_idx += 1
                else
                    prevcoord .= thiscoord
                end

                # 输出 thiscoord
                for j in 1:3
                    coords[j, coord_idx] = thiscoord[j] * inv_precision
                end
                coord_idx += 1
            end
        else
            # 没有 run-length 编码，直接输出
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
        sizesmall .= MAGICINTS[smallidx+1]

        i += 1
    end

    return precision, coords
end
