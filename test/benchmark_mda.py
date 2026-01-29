import MDAnalysis as mda
import time
import os

def run_benchmark():
    u = mda.Universe("test.gro", "test.xtc")
    print(f"Frames: {len(u.trajectory)}")
    print(f"Atoms: {len(u.atoms)}")
    
    # Warmup
    for ts in u.trajectory:
        break
    
    # Run 1
    start_time = time.time()
    count = 0
    for ts in u.trajectory:
        count += 1
    end_time = time.time()
    duration1 = end_time - start_time
    print(f"MDAnalysis Run 1: {duration1:.4f} seconds")
    
    # Run 2
    start_time = time.time()
    count = 0
    for ts in u.trajectory:
        count += 1
    end_time = time.time()
    duration2 = end_time - start_time
    print(f"MDAnalysis Run 2: {duration2:.4f} seconds")
    
    avg = (duration1 + duration2) / 2
    print(f"\nAverage time: {avg:.4f} seconds")
    print(f"Average FPS: {count / avg:.2f}")

if __name__ == "__main__":
    run_benchmark()
