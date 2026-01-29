using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using XTCReader
using Printf

function run_benchmark()
    println("=== XTCReader.jl 性能基准测试 ===\n")

    # 预热
    println("预热中...")
    for frame in eachframe("test.xtc")
        break
    end

    # 测试 eachframe 迭代器（零分配版本）
    println("\n测试 eachframe 迭代器...")
    start_time = time()
    count = 0
    for frame in eachframe("test.xtc")
        count += 1
    end
    end_time = time()

    duration = end_time - start_time
    @printf("XTCReader.jl (eachframe) 耗时: %.4f 秒\n", duration)
    @printf("帧率 (FPS): %.2f\n", count / duration)
    println("总帧数: $count")

    # 获取原子数
    traj = read_xtc("test.xtc")
    println("原子数: $(traj.natoms)")

    # 再运行一次取平均
    println("\n第二次运行...")
    start_time = time()
    count = 0
    for frame in eachframe("test.xtc")
        count += 1
    end
    end_time = time()

    duration2 = end_time - start_time
    @printf("XTCReader.jl (eachframe) 耗时: %.4f 秒\n", duration2)
    @printf("帧率 (FPS): %.2f\n", count / duration2)

    avg_duration = (duration + duration2) / 2
    @printf("\n平均耗时: %.4f 秒\n", avg_duration)
    @printf("平均帧率: %.2f FPS\n", count / avg_duration)
end

run_benchmark()
