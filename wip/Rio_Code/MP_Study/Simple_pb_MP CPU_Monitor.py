import pybamm; import os;import numpy as np
import multiprocessing
import time
from scipy.io import savemat,loadmat;

def Run_model(Temper_i ):
    model = pybamm.lithium_ion.SPMe()
    ChemistryChen=pybamm.parameter_sets.Chen2020_coupled   
    Para_0=pybamm.ParameterValues(chemistry=ChemistryChen)
    Para_0.update(   {'Ambient temperature [K]':273.15+Temper_i });
    sim = pybamm.Simulation(model,parameter_values = Para_0,)
    sim.solve([0, 360])
    solution = sim.solution
    t = solution["Time [s]"].entries
    V = solution["Terminal voltage [V]"].entries
    DataPack = [t,V];

    return DataPack

if __name__ == "__main__":
    pool = multiprocessing.Pool(6)
    start_time = time.perf_counter()
    processes = [pool.apply_async(Run_model, args=(x,)) for x in np.arange(0,1,5e-3)]
    result = [p.get() for p in processes]
    finish_time = time.perf_counter()
    print(f"Program finished in {finish_time-start_time} seconds")
    # print(result)