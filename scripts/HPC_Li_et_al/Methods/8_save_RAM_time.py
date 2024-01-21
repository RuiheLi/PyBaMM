from memory_profiler import profile

import matplotlib.pyplot as plt
import numpy as np
import pybamm
import pandas as pd

def Run_model_1st(exp_text,param,save_at_cycles,var_pts,):
    
    model = pybamm.lithium_ion.DFN(options={
        "SEI": "interstitial-diffusion limited", 
        "SEI film resistance": "distributed", 
        "SEI porosity change": "true",})
    exp = pybamm.Experiment( exp_text*save_at_cycles )
    sim = pybamm.Simulation(
        model, 
        experiment=exp,
        parameter_values=param,
        solver=pybamm.CasadiSolver(return_solution_if_failed_early=True),
        var_pts=var_pts,
    )
    sol = sim.solve(save_at_cycles=save_at_cycles,calc_esoh=False)
    
    return sol,model

def Run_model_later(model_old, sol_old, exp_text,param, save_at_cycles,var_pts,):
    model = model_old.set_initial_conditions_from(sol_old, inplace=False)
    exp = pybamm.Experiment( exp_text*save_at_cycles )
    sim = pybamm.Simulation(
        model, 
        experiment=exp,
        parameter_values=param,
        solver=pybamm.CasadiSolver(return_solution_if_failed_early=True),
        var_pts=var_pts,
    )
    sol = sim.solve(save_at_cycles=save_at_cycles,calc_esoh=False)
    return sol,model

def Post_processing(sol,cyc,step):
    time = (
        sol.cycles[cyc].steps[step]["Time [h]"].entries 
        -
        sol.cycles[cyc].steps[step]["Time [h]"].entries[0] )
    v = sol.cycles[cyc].steps[step]["Terminal voltage [V]"].entries
    return time, v


def Run_Big_Model(exp_text,param,save_at_cycles,var_pts,total_cycles,return_sol):
    start_1 = pybamm.Timer()
    Time = []; Voltage = []; Cyc_count = []; Sol_all = []
    i=1; n = total_cycles / save_at_cycles
    # Run the first N cycles:
    sol,model = Run_model_1st(exp_text,param,save_at_cycles,var_pts,)
    if return_sol == True:
        Sol_all.append(sol)
    time, v = Post_processing(sol,0,1)
    Cyc_count.append(i+1)
    Time.append(time)
    Voltage.append(v)
    time, v = Post_processing(sol,-1,1)
    Cyc_count.append(save_at_cycles)
    Time.append(time)
    Voltage.append(v)
    while i<n:
        sol,model = Run_model_later(model, sol, exp_text,param, save_at_cycles,var_pts,)
        if return_sol == True:
            Sol_all.append(sol)
        time, v = Post_processing(sol,-1,1)
        Time.append(time)
        Voltage.append(v)
        Cyc_count.append(save_at_cycles*i)
        i += 1
    print(f'running time is {start_1.time()}')
    return Time, Voltage, Cyc_count,Sol_all

exp_text = [(
    "Hold at 4.2V until C/100", 
    "Discharge at 1C until 2.5V", 
    "Rest for 1 hours", 
    "Charge at 1C until 4.2V", 
)]
param = pybamm.ParameterValues("OKane2022")
param.update({"Inner SEI lithium interstitial diffusivity [m2.s-1]": 3.5e-19})  

var_pts={"x_n": 5, "x_s": 5,"x_p": 5, "r_n": 30, "r_p": 20, }
save_at_cycles = 5; total_cycles = 1000

# Return solution 
# instantiating the decorator
@profile
def Main():
    Time, Voltage, Cyc_count,Sol_all = Run_Big_Model(
        exp_text,param,save_at_cycles,var_pts,
        total_cycles,return_sol=True)
if __name__ == "__main__":
    Main()
