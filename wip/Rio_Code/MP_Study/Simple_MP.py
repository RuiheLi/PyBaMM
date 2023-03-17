
import multiprocessing
import time
 
def cube(x):
    return x**3
 
if __name__ == "__main__":
    pool = multiprocessing.Pool(6)
    start_time = time.perf_counter()
    processes = [pool.apply_async(cube, args=(x,)) for x in range(1,1000)]
    result = [p.get() for p in processes]
    finish_time = time.perf_counter()
    print(f"Program finished in {finish_time-start_time} seconds")
    print(result)