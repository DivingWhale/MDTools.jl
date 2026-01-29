using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MDTools
using Printf

function run_benchmark()
    println("=== XTCReader.jl Performance Benchmark ===\n")

    # Warm-up
    println("Warming up...")
    for frame in eachframe("test.xtc")
        break
    end

    # Test eachframe iterator (zero-allocation version)
    println("\nTesting eachframe iterator...")
    start_time = time()
    count = 0
    for frame in eachframe("test.xtc")
        count += 1
    end
    end_time = time()

    duration = end_time - start_time
    @printf("XTCReader.jl (eachframe) time: %.4f seconds\n", duration)
    @printf("Frame rate (FPS): %.2f\n", count / duration)
    println("Total frames: $count")

    # Get number of atoms
    traj = read_xtc("test.xtc")
    println("Number of atoms: $(traj.natoms)")

    # Run again to get average
    println("\nSecond run...")
    start_time = time()
    count = 0
    for frame in eachframe("test.xtc")
        count += 1
    end
    end_time = time()

    duration2 = end_time - start_time
    @printf("XTCReader.jl (eachframe) time: %.4f seconds\n", duration2)
    @printf("Frame rate (FPS): %.2f\n", count / duration2)

    avg_duration = (duration + duration2) / 2
    @printf("\nAverage time: %.4f seconds\n", avg_duration)
    @printf("Average frame rate: %.2f FPS\n", count / avg_duration)
end

run_benchmark()
