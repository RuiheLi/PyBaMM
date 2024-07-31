""" Functions to actually run PyBaMM model and catch any errors """

import csv, random, os, gc
import pybamm as pb;import pandas as pd   ;import numpy as np;import os;
import matplotlib.pyplot as plt;import os;#import imageio;import timeit
from scipy.io import savemat,loadmat;from pybamm import constants,exp,sqrt;
import matplotlib as mpl; 
from multiprocessing import Queue, Process, set_start_method
from queue import Empty
import openpyxl
import traceback
import random;import time, signal
from .TaskResult import TaskResult
from .Plot import *
from .Post_process import *
from .Fun_Upgrade import *

# DEFINE my callback:
class RioCallback(pb.callbacks.Callback):
    def __init__(self, logfile=None):
        self.logfile = logfile
        self.success  = True
        if logfile is None:
            # Use pybamm's logger, which prints to command line
            self.logger = pb.logger
        else:
            # Use a custom logger, this will have its own level so set it to the same
            # level as the pybamm logger (users can override this)
            self.logger = pb.get_new_logger(__name__, logfile)
            self.logger.setLevel(pb.logger.level)
    
    def on_experiment_error(self, logs):
        self.success  = False
    def on_experiment_infeasible(self, logs):
        self.success  = False

class Experiment_error_infeasible(ValueError):
    pass



###################################################################
#############    From Patrick from Github        ##################
###################################################################
class TimeoutFunc(object):
    def __init__(self, func, timeout=None, timeout_val=None):
        assert callable(func), 'Positional argument 1 must be a callable method'
        self.func = func
        self.timeout = timeout
        self.timeout_val = timeout_val
    def _queued_f(self, *args, queue, **kwargs):
        Result_List = self.func(*args, **kwargs)
        try:
            queue.put(Result_List)
        except BrokenPipeError as exc:
            pass

    def __call__(self, *args, **kwargs):
        q = Queue(1)
        p = Process(target=self._queued_f, args=args,
            kwargs={**kwargs, 'queue': q})
        p.daemon = True
        p.start()
        try:
            Result_List = q.get(timeout=self.timeout)
        except Empty as exc: 
            # todo: need to incooperate other cases such as pb.expression_tree.exceptions.ModelError,
            #      pb.expression_tree.exceptions.SolverError
            p.terminate()
            Call_ref = RioCallback() 
            # careful! the length of Result_List should be same 
            #   as what you get from main function!
            Result_List = [         
                self.timeout_val,
                self.timeout_val,
                Call_ref,
                self.timeout_val]
        return Result_List  

# define to kill too long runs - Jianbo Huang kindly writes this
class TimeoutError(Exception):
	pass


# define the model and run break-in cycle - 
# input parameter: model_options, Experiment_Breakin, Para_0, mesh_list, submesh_strech
# output: Sol_0 , Model_0, Call_Breakin
def Run_Breakin(config, Experiment_Breakin, Para_0, Model_0):

    mesh_list = config.model_config.mesh_list
    submesh_strech = config.model_config.submesh_strech
    # update 220926 - add diffusivity and conductivity as variables:
    c_e = Model_0.variables["Electrolyte concentration [mol.m-3]"]
    T = Model_0.variables["Cell temperature [K]"]
    D_e = Para_0["Electrolyte diffusivity [m2.s-1]"]
    sigma_e = Para_0["Electrolyte conductivity [S.m-1]"]
    Model_0.variables["Electrolyte diffusivity [m2.s-1]"] = D_e(c_e, T)
    Model_0.variables["Electrolyte conductivity [S.m-1]"] = sigma_e(c_e, T)
    var = pb.standard_spatial_vars  
    var_pts = {
        var.x_n: int(mesh_list[0]),  
        var.x_s: int(mesh_list[1]),  
        var.x_p: int(mesh_list[2]),  
        var.r_n: int(mesh_list[3]),  
        var.r_p: int(mesh_list[4]),  }       
    submesh_types = Model_0.default_submesh_types
    if submesh_strech == "nan":
        pass
    else:
        particle_mesh = pb.MeshGenerator(
            pb.Exponential1DSubMesh, 
            submesh_params={"side": "right", "stretch": submesh_strech})
        submesh_types["negative particle"] = particle_mesh
        submesh_types["positive particle"] = particle_mesh 

    # generate a list of capacity increase perturbation:
    import random
    Cap_in_perturbation = []
    try_no = 4
    for k in range(try_no):
        random_number = 0
        while abs(random_number) <= 1e-6:
            random_number = random.uniform(-1e-5, 1e-5)
        Cap_in_perturbation.append(random_number)
    Cap_in_perturbation[0] = 0.0 # first time must be zero
    Cap_in_perturbation[1] = 3.9116344182112036E-4
    c_s_neg_baseline = Para_0["Initial concentration in negative electrode [mol.m-3]"]
    # update 240603 - try shift neg soc 4 times until give up 
    i_run_try = 0
    while i_run_try<try_no:
        try: 
            Sim_0    = pb.Simulation(
                Model_0,        experiment = Experiment_Breakin,
                parameter_values = Para_0,
                solver = pb.CasadiSolver(),
                var_pts=var_pts,
                submesh_types=submesh_types) 
            Call_Breakin = RioCallback()    
            Sol_0    = Sim_0.solve(calc_esoh=False,callbacks=Call_Breakin)
        except (
            pb.expression_tree.exceptions.ModelError,
            pb.expression_tree.exceptions.SolverError
            ) as e:
            Sol_0 = "Model error or solver error"
            str_err = (
                f"Fail to run break in due to {Sol_0} for the {i_run_try}th time"
                f" with perturbation of {Cap_in_perturbation[i_run_try]:.2e}Ah")
            print();print(str_err);print()
            i_run_try += 1
        else:
            str_err = (
                f"Succeed to run break in for the {i_run_try}th time "
                f"with perturbation of {Cap_in_perturbation[i_run_try]:.2e}Ah")
            print();print(str_err);print()
            break
    DeBug_List = ["Place holder"] # TODO
    Result_list_breakin = [Model_0, Sol_0, Call_Breakin, DeBug_List]

    return Result_list_breakin


def Run_Breakin_lump_with_try_catch(
        task_result,config,
        Experiment_Breakin, Para_0, Model_0, SmallTimer):
    """
    Function to wrap Run_Breakin.

    Main purpose is to wrap "try ... except"

    Parameters
    ----------
    task_result : TaskResult class
        collect all outputs for the whole task
    config : TaskConfig class
        collect all inputs for the whole task
    Experiment_Breakin : pybamm.experiment
        define break-in cycle
    Para_0 : pybamm.parameter
        define parameter, destination of most input
    Model_0 : pybamm.model
        destination of many input
    SmallTimer : pybamm.timer()
        check short time
    """

    ##########    2-1: Define model and run break-in cycle
    Scan_No = config.Scan_No
    Re_No = config.global_config.Re_No
    
    try:  
        # the following turns on for HPC only!
        if config.global_config.Timeout == True:
            timeout_RPT = TimeoutFunc(
                Run_Breakin, 
                timeout=config.global_config.Timelimit, 
                timeout_val=config.global_config.Timeout_text)
            
            Result_list_breakin  = timeout_RPT(
                config, Experiment_Breakin, Para_0, Model_0)
        else:
            Result_list_breakin  = Run_Breakin(
                config, Experiment_Breakin, Para_0, Model_0)
            
        [Model_0,Sol_0,Call_Breakin,DeBug_List_Breakin] = Result_list_breakin
        if config.global_config.Return_Sol == True:
            task_result.Sol_RPT.append(Sol_0)
        if Call_Breakin.success == False:
            print("Fail due to Experiment error or infeasible")
            1/0
        if Sol_0 == config.global_config.Timeout_text: # to do: distinguish different failure cases
            print("Fail due to Timeout")
            1/0
        if Sol_0 == "Model error or solver error":
            print("Fail due to Model error or solver error")
            1/0
    except ZeroDivisionError as e:
        str_error_Breakin = str(e)
        if config.global_config.Check_Short_Time == True:
            str_error_Breakin = f"Scan {Scan_No} Re {Re_No}: Fail break-in "
            f"cycle within {SmallTimer.time()}, need to exit the whole "
            f"scan now due to {str_error_Breakin} but do not know how!"

            print(str_error_Breakin)
            SmallTimer.reset()
        else:
            str_error_Breakin = f"Scan {Scan_No} Re {Re_No}: Fail break-in "
            "cycle, need to exit the whole scan now due to "
            f"{str_error_Breakin} but do not know how!"

            print(str_error_Breakin)
        Flag_Breakin = False 
    else:
        if config.global_config.Check_Short_Time == True:    
            print(f"Scan {Scan_No} Re {Re_No}: Finish break-in cycle within {SmallTimer.time()}")
            SmallTimer.reset()
        else:
            print(f"Scan {Scan_No} Re {Re_No}: Finish break-in cycle")
        # post-process for break-in cycle - 0.1C only
        keys_all_RPT = task_result.keys_all_RPT
        task_result.my_dict_RPT = GetSol_dict(
            task_result.my_dict_RPT,keys_all_RPT, Sol_0, 0, 
            config.exp_config.step_0p1C_CD, 
            config.exp_config.step_0p1C_CC,
            config.exp_config.step_0p1C_RE, 
            config.exp_config.step_AGE_CV   )
        # update 230517 - Get R from C/2 discharge only, discard GITT
        cap_full = Para_0["Nominal cell capacity [A.h]"]
        if config.global_config.R_from_GITT: 
            Res_midSOC,Res_full,SOC_Res = Get_0p1s_R0(Sol_0,cap_full)
        else: 
            step_0P5C_CD = Sol_0.cycles[0].steps[config.exp_config.step_0p5C_CD]
            Res_midSOC,Res_full,SOC_Res = Get_R_from_0P5C_CD(step_0P5C_CD,cap_full)
        task_result.my_dict_RPT["SOC_Res"].append(SOC_Res)
        task_result.my_dict_RPT["Res_full"].append(Res_full)
        task_result.my_dict_RPT["Res_midSOC"].append(Res_midSOC)    

        task_result.my_dict_RPT["avg_Age_T"].append(
            config.exp_config.Age_T_in_K-273.15)  # Update add 230617              
        del SOC_Res,Res_full,Res_midSOC; gc.collect() # ensure ram by forcing the garbage collector to run
        # add zero for break in
        task_result.my_dict_RPT["Cycle_RPT"].append(0)
        task_result.my_dict_RPT["Cyc_Update_Index"].append(0)
        Flag_Breakin = True

        str_error_Breakin = "nan"

        print(f"Scan {Scan_No} Re {Re_No}: Finish post-process for break-in cycle within {SmallTimer.time()}")
        SmallTimer.reset()
       
    task_result.Flag_Breakin = Flag_Breakin
    task_result.Call_Breakin = Call_Breakin
    task_result.str_error_Breakin = str_error_Breakin
    task_result.DeBug_List_Breakin = DeBug_List_Breakin
    del Flag_Breakin,Call_Breakin,str_error_Breakin,DeBug_List_Breakin; gc.collect()

    return task_result, config, Para_0, Model_0, SmallTimer, Sol_0

def Run_RPT_lump_with_try_catch(task_result,config,
    Experiment_RPT, Paraupdate, 
    Model_Dry_old, Sol_Dry_old, SmallTimer, k, i):

    Scan_No = config.Scan_No
    Re_No = config.global_config.Re_No
    Timeout = config.global_config.Timeout    
    Timelimit = config.global_config.Timelimit
    Timeout_text = config.global_config.Timeout_text

    try:
        # Timelimit = int(60*60*2)
        if Timeout == True:
            timeout_RPT = TimeoutFunc(
                Run_Model_Base_On_Last_Solution_RPT, 
                timeout=Timelimit, 
                timeout_val=Timeout_text)
            Result_list_RPT = timeout_RPT(
                config, Experiment_RPT, Paraupdate, 
                Model_Dry_old, Sol_Dry_old, 
            )
        else:
            Result_list_RPT = Run_Model_Base_On_Last_Solution_RPT(
                config, Experiment_RPT, Paraupdate, 
                Model_Dry_old, Sol_Dry_old, 
            )
        [Model_Dry_i, Sol_Dry_i,Call_RPT,DeBug_List_RPT]  = Result_list_RPT
        if config.global_config.Return_Sol  == True:
            task_result.Sol_RPT.append(Sol_Dry_i)
        #print(f"Temperature for RPT is now: {RPT_T_in_K}")  
        if Call_RPT.success == False:
            print("Fail due to Experiment error or infeasible")
            str_error_RPT = "Experiment error or infeasible"
            1/0 
        if Sol_Dry_i == Timeout_text:
            #print("Fail due to Timeout")
            str_error_RPT = "Timeout"
            1/0
        if Sol_Dry_i == "Model error or solver error":
            print("Fail due to Model error or solver error")
            str_error_RPT = "Model error or solver error"
            1/0
    except ZeroDivisionError as e:
        str_error_RPT = (
            f"""Scan {Scan_No} Re {Re_No}: Fail during No.
            {task_result.my_dict_RPT['Cyc_Update_Index'][-1]} RPT
            cycles within {SmallTimer.time()}, due to {str_error_RPT}"""
        )
        print(str_error_RPT)
        SmallTimer.reset()
        #  break   # TODO consider how to break in  this case
    else:
        # post-process for RPT
        task_result.my_dict_RPT["Cyc_Update_Index"].append(
            task_result.my_dict_AGE["Agecycle_count"])

        print(
            f"""Scan {Scan_No} Re {Re_No}: Finish No."
            {task_result.my_dict_RPT['Cyc_Update_Index'][-1]}"
            RPT cycles within {SmallTimer.time()}"""
        )
        SmallTimer.reset()

        task_result.my_dict_RPT = GetSol_dict (
            task_result.my_dict_RPT,task_result.keys_all_RPT, 
            Sol_Dry_i, 0,
            config.exp_config.step_0p1C_CD, 
            config.exp_config.step_0p1C_CC,
            config.exp_config.step_0p1C_RE, 
            config.exp_config.step_AGE_CV) 
        task_result.my_dict_RPT["Cycle_RPT"].append(k)
        # this is the average T during previous ageing cycles 
        task_result.my_dict_RPT["avg_Age_T"].append(
            np.mean(task_result.my_dict_AGE["avg_Age_T"])
        ) 
        
        # update 230517 - Get R from C/2 discharge only, discard GITT
        cap_full = Paraupdate["Nominal cell capacity [A.h]"] # 5
        if config.global_config.R_from_GITT: 
            Res_midSOC,Res_full,SOC_Res = Get_0p1s_R0(Sol_Dry_i,cap_full)
        else: 
            step_0P5C_CD = Sol_Dry_i.cycles[0].steps[config.exp_config.step_0p5C_CD]
            Res_midSOC,Res_full,SOC_Res = Get_R_from_0P5C_CD(step_0P5C_CD,cap_full)
        task_result.my_dict_RPT["SOC_Res"].append(SOC_Res)
        task_result.my_dict_RPT["Res_full"].append(Res_full)
        task_result.my_dict_RPT["Res_midSOC"].append(Res_midSOC)             
        del SOC_Res,Res_full,Res_midSOC; gc.collect() 
    
        if config.model_config.DryOut == "On":
            task_result.Update_mdic_dry()
        Para_0_Dry_old = Paraupdate    
        Model_Dry_old = Model_Dry_i     
        Sol_Dry_old = Sol_Dry_i      
        del Paraupdate,Model_Dry_i,Sol_Dry_i; 
        gc.collect() # ensure ram by forcing the garbage collector to run

        print(
            f"""Scan {Scan_No} Re {Re_No}: Finish post-process for No.
            {task_result.my_dict_RPT["Cyc_Update_Index"][-1]} 
            RPT cycles within {SmallTimer.time()}""")
        str_error_RPT = "nan"
        SmallTimer.reset()

    task_result.Call_Age = Call_RPT
    task_result.str_error_AGE = str_error_RPT
    task_result.DeBug_List_AGE = DeBug_List_RPT

    return (
        task_result, config, Para_0_Dry_old, Model_Dry_old, 
        SmallTimer, Sol_Dry_old)

def Run_Age_lump_with_try_catch(
    task_result, config, Experiment_Long, Paraupdate, 
    Model_Dry_old, Sol_Dry_old, SmallTimer, k, i):

    Scan_No = config.Scan_No
    Re_No = config.global_config.Re_No
    Timeout = config.global_config.Timeout
    try:
        if Timeout == True:
            timeout_AGE = TimeoutFunc(
                Run_Model_Base_On_Last_Solution, 
                timeout=config.global_config.Timelimit, 
                timeout_val=config.global_config.Timeout_text)
            Result_list_AGE = timeout_AGE( 
                config, Experiment_Long, Paraupdate, 
                Model_Dry_old, Sol_Dry_old)
        else:
            Result_list_AGE = Run_Model_Base_On_Last_Solution( 
                config, Experiment_Long, Paraupdate, 
                Model_Dry_old, Sol_Dry_old)
        [Model_Dry_i, Sol_Dry_i, Call_Age, DeBug_List_AGE] = Result_list_AGE
        
        if config.global_config.Return_Sol == True:
            task_result.Sol_AGE.append(Sol_Dry_i)
        if "Partially" in DeBug_List_AGE[-1]:
            Flag_partial_AGE = True
            succeed_cycs = len(Sol_Dry_i.cycles) 
            if succeed_cycs < config.exp_config.update:
                print(f"Instead of {config.exp_config.update}, "
                f"succeed only {succeed_cycs} cycles")
        else:
            Flag_partial_AGE = False
            if Call_Age.success == False:
                print("Fail due to Experiment error or infeasible")
                str_error_AGE = "Experiment error or infeasible"
                1/0
            elif Sol_Dry_i == "Model error or solver error":
                print("Fail due to Model error or solver error")
                str_error_AGE = "Model error or solver error"
                1/0
            else:
                pass
        if Sol_Dry_i == config.global_config.Timeout_text: # fail due to timeout
            print("Fail due to Timeout")
            str_error_AGE = "Timeout"
            1/0
    except ZeroDivisionError as e: # ageing cycle fails
        if config.global_config.Check_Short_Time == True:    
            str_error_AGE = f"""Scan {Scan_No} Re {Re_No}: Fail 
                during No. {task_result.my_dict_RPT["Cyc_Update_Index"][-1]} 
                ageing cycles within {SmallTimer.time()} 
                due to {str_error_AGE}"""
            print(str_error_AGE)
            SmallTimer.reset()
        else:
            str_error_AGE = f"""Scan {Scan_No} Re {Re_No}: Fail during No.
                {task_result.my_dict_RPT["Cyc_Update_Index"][-1]} ageing 
                cycles due to {str_error_AGE}"""
            print(str_error_AGE)
        Flag_AGE = False
        # break - Mark: TODO Ask Dan, how to break in this case
    else:                           # ageing cycle SUCCEED
        succeed_cycs = len(Sol_Dry_i.cycles) 
        Para_0_Dry_old = Paraupdate; Model_Dry_old = Model_Dry_i; Sol_Dry_old = Sol_Dry_i;   
        del Paraupdate,Model_Dry_i,Sol_Dry_i; gc.collect()
    
        print(f"""Scan {Scan_No} Re {Re_No}: Finish No.
              {task_result.my_dict_RPT["Cyc_Update_Index"][-1]} 
              ageing cycles within {SmallTimer.time()}""")
        SmallTimer.reset()
        
        # post-process for first ageing cycle and every -1 ageing cycle
        step_AGE_CD = config.exp_config.step_AGE_CD 
        step_AGE_CC = config.exp_config.step_AGE_CC
        step_0p1C_RE= config.exp_config.step_0p1C_RE
        step_AGE_CV = config.exp_config.step_AGE_CV 
        keys_all_AGE = task_result.keys_all_AGE
        if k==0 and i==0:    
            task_result.my_dict_AGE = GetSol_dict (
                task_result.my_dict_AGE,keys_all_AGE, Sol_Dry_old, 0, 
                step_AGE_CD, step_AGE_CC, step_0p1C_RE, step_AGE_CV)     
            task_result.my_dict_AGE["Cycle_AGE"].append(1)
        # update 240111
        if Flag_partial_AGE == True:
            try:
                task_result.my_dict_AGE = GetSol_dict (
                    task_result.my_dict_AGE,keys_all_AGE, Sol_Dry_old, -1, 
                    step_AGE_CD, step_AGE_CC, step_0p1C_RE, step_AGE_CV)    
            except IndexError:
                print("The last cycle is incomplete, try [-2] cycle")
                try:
                    task_result.my_dict_AGE = GetSol_dict (
                        task_result.my_dict_AGE,keys_all_AGE, Sol_Dry_old, -2, 
                        step_AGE_CD, step_AGE_CC, step_0p1C_RE, step_AGE_CV)   
                except IndexError:
                    print("[-2] cycle also does not work, try first one")
                    try:
                        task_result.my_dict_AGE = GetSol_dict (
                            task_result.my_dict_AGE,keys_all_AGE, Sol_Dry_old, 0, 
                            step_AGE_CD, step_AGE_CC, step_0p1C_RE, step_AGE_CV)  
                    except:
                        print("Still does not work, less than one cycle, we are in trouble")
        else:
            try:
                task_result.my_dict_AGE = GetSol_dict (
                    task_result.my_dict_AGE,keys_all_AGE, Sol_Dry_old, -1, 
                    step_AGE_CD, step_AGE_CC, step_0p1C_RE, step_AGE_CV)    
            except:
                print("GetSol_dict fail for a complete ageing set for unknown reasons!!!")
        task_result.my_dict_AGE["Agecycle_count"] +=  succeed_cycs 

        task_result.my_dict_AGE["avg_Age_T"].append(np.mean(
            Sol_Dry_old["Volume-averaged cell temperature [C]"].entries))
        
        task_result.my_dict_AGE["Cycle_AGE"].append(
            task_result.my_dict_AGE["Agecycle_count"])     
              
        task_result.my_dict_RPT["Cyc_Update_Index"].append(
            task_result.my_dict_AGE["Agecycle_count"])
        
        if config.model_config.DryOut == "On":
            task_result.Update_mdic_dry()
        Flag_AGE = True
        str_error_AGE = "nan"
   
        print(f"""Scan {Scan_No} Re {Re_No}: Finish post-process for 
              No. {task_result.my_dict_RPT["Cyc_Update_Index"][-1]} 
              ageing cycles within {SmallTimer.time()}""")
        SmallTimer.reset()

    # Para_0_Dry_old = Paraupdate; Model_Dry_old = Model_Dry_i; Sol_Dry_old = Sol_Dry_i;  
    task_result.Flag_partial_AGE = Flag_partial_AGE
    task_result.Flag_AGE = Flag_AGE
    task_result.Call_Age = Call_Age
    task_result.str_error_AGE = str_error_AGE
    task_result.DeBug_List_AGE = DeBug_List_AGE

    return (
        task_result, config, Para_0_Dry_old, Model_Dry_old, 
        SmallTimer, Sol_Dry_old)



def Run_Model_Base_On_Last_Solution_RPT( 
    config, Experiment_RPT, Para_update, Model, Sol):


    Update_Cycles = config.exp_config.update
    Age_T_in_K = config.exp_config.Age_T_in_K
    mesh_list = config.model_config.mesh_list
    submesh_strech = config.model_config.submesh_strech
    Ratio_CeLi = Para_update[
        "Ratio of Li-ion concentration change in electrolyte"
        " consider solvent consumption"]
    # print("Model is now using average EC Concentration of:",Para_update['Bulk solvent concentration [mol.m-3]'])
    # print("Ratio of electrolyte dry out in jelly roll is:",Para_update['Ratio of electrolyte dry out in jelly roll'])
    # print("Model is now using an electrode width of:",Para_update['Electrode width [m]'])

    if isinstance(Sol, pb.solvers.solution.Solution):
        list_short,dict_short = Get_Last_state(Model, Sol)
    elif isinstance(Sol, list):
        [dict_short, getSth] = Sol
        list_short = []
        for var, equation in Model.initial_conditions.items():
            list_short.append(var._name)
    else:
        print("!! Big problem, Sol here is neither solution or list")

    dict_short["Negative electrode porosity times concentration [mol.m-3]"] = (
        dict_short["Negative electrode porosity times concentration [mol.m-3]"] 
        * Ratio_CeLi )# important: update sol here!
    dict_short["Separator porosity times concentration [mol.m-3]"] = (
        dict_short["Separator porosity times concentration [mol.m-3]"] 
        * Ratio_CeLi )# important: update sol here!
    dict_short["Positive electrode porosity times concentration [mol.m-3]"] = (
        dict_short["Positive electrode porosity times concentration [mol.m-3]"] 
        * Ratio_CeLi )# important: update sol here!
    Model_new = Model.set_initial_conditions_from(dict_short, inplace=False)
    Para_update.update(   {'Ambient temperature [K]':Age_T_in_K });
    
    var = pb.standard_spatial_vars  
    var_pts = {
        var.x_n: int(mesh_list[0]),  
        var.x_s: int(mesh_list[1]),  
        var.x_p: int(mesh_list[2]),  
        var.r_n: int(mesh_list[3]),  
        var.r_p: int(mesh_list[4]),  }
    submesh_types = Model.default_submesh_types
    if submesh_strech == "nan":
        pass
    else:
        particle_mesh = pb.MeshGenerator(
            pb.Exponential1DSubMesh, 
            submesh_params={"side": "right", "stretch": int(submesh_strech)})
        submesh_types["negative particle"] = particle_mesh
        submesh_types["positive particle"] = particle_mesh 
    Simnew = pb.Simulation(
        Model_new,
        experiment = Experiment_RPT, 
        parameter_values=Para_update, 
        solver = pb.CasadiSolver(),
        var_pts = var_pts,
        submesh_types=submesh_types
    )
    Call_RPT = RioCallback()  # define callback
    # update 231208 - try 3 times until give up 
    i_run_try = 0
    while i_run_try<3:
        try:
            Sol_new = Simnew.solve(
                calc_esoh=False,
                callbacks=Call_RPT) 
            if Call_RPT.success == False:
                raise Experiment_error_infeasible("Self detect")        
        except (
            pb.expression_tree.exceptions.ModelError,
            pb.expression_tree.exceptions.SolverError
            ) as e:
            i_run_try += 1
            Sol_new = "Model error or solver error"
            str_err = f"Fail to run RPT due to {Sol_new} for the {i_run_try}th time"
            print(str_err)
            DeBug_List = [ Para_update, Update_Cycles, dict_short, str_err ]
        except Experiment_error_infeasible as custom_error:
            i_run_try += 1
            Sol_new = "Experiment error or infeasible"
            DeBug_List = [
                Model, Model_new, Call_RPT, Simnew, Sol, Sol_new, 
                Para_update, Experiment_RPT, Update_Cycles,
                Age_T_in_K ,mesh_list,submesh_strech, var_pts,
                submesh_types,  list_short, dict_short, 
            ]
            str_err = f"{Sol_new}: {custom_error} for RPT for the {i_run_try}th time"
            print(str_err)
            DeBug_List = [ Para_update, Update_Cycles, dict_short, str_err ]
        else:
            # add 230221 - update 230317 try to access Throughput capacity more than once
            i_run_try += 1
            i_try = 0
            while i_try<3:
                try:
                    getSth2 = Sol_new['Throughput capacity [A.h]'].entries[-1]
                except:
                    i_try += 1
                    print(f"Fail to read Throughput capacity for the {i_try}th time")
                else:
                    break
            i_try = 0
            while i_try<3:
                try:
                    getSth = Sol['Throughput capacity [A.h]'].entries[-1]
                except:
                    i_try += 1
                    print(f"Fail to read Throughput capacity for the {i_try}th time")
                else:
                    break
            Sol_new['Throughput capacity [A.h]'].entries += getSth
            # update 23-11-17 change method to get throughput capacity to avioid problems of empty solution:
            thr_tot = Get_ThrCap(Sol_new)
            Sol_new['Throughput capacity [A.h]'].entries[-1] = getSth + thr_tot # only the last one is true
            DeBug_List = "Empty"
            print(f"Succeed to run RPT for the {i_run_try}th time")
            break # terminate the loop of trying to solve if you can get here. 
    # print("Solved this model in {}".format(ModelTimer.time()))
    Result_List_RPT = [Model_new, Sol_new,Call_RPT,DeBug_List]
    return Result_List_RPT


# Define a function to calculate based on previous solution
def Run_Model_Base_On_Last_Solution( 
    config, Experiment_Long, Para_update, Model, Sol):
    # 
    Update_Cycles = config.exp_config.update
    Age_T_in_K = config.exp_config.Age_T_in_K
    mesh_list = config.model_config.mesh_list
    submesh_strech = config.model_config.submesh_strech

    Ratio_CeLi = Para_update[
        "Ratio of Li-ion concentration change in electrolyte"
        " consider solvent consumption"]
    if isinstance(Sol, pb.solvers.solution.Solution):
        list_short,dict_short = Get_Last_state(Model, Sol)
    elif isinstance(Sol, list):
        [dict_short, getSth] = Sol
        list_short = []
        for var, equation in Model.initial_conditions.items():
            list_short.append(var._name)
    else:
        print("!! Big problem, Sol here is neither solution or list")

    dict_short["Negative electrode porosity times concentration [mol.m-3]"] = (
        dict_short["Negative electrode porosity times concentration [mol.m-3]"] 
        * Ratio_CeLi )# important: update sol here!
    dict_short["Separator porosity times concentration [mol.m-3]"] = (
        dict_short["Separator porosity times concentration [mol.m-3]"] 
        * Ratio_CeLi )# important: update sol here!
    dict_short["Positive electrode porosity times concentration [mol.m-3]"] = (
        dict_short["Positive electrode porosity times concentration [mol.m-3]"] 
        * Ratio_CeLi )# important: update sol here!
    Model_new = Model.set_initial_conditions_from(dict_short, inplace=False)
    Para_update.update({'Ambient temperature [K]':Age_T_in_K });   # run model at 45 degree C
    
    var = pb.standard_spatial_vars  
    var_pts = {
        var.x_n: int(mesh_list[0]),  
        var.x_s: int(mesh_list[1]),  
        var.x_p: int(mesh_list[2]),  
        var.r_n: int(mesh_list[3]),  
        var.r_p: int(mesh_list[4]),  
        }
    submesh_types = Model.default_submesh_types
    if submesh_strech == "nan":
        pass
    else:
        particle_mesh = pb.MeshGenerator(
            pb.Exponential1DSubMesh, 
            submesh_params={"side": "right", "stretch": int(submesh_strech)})
        submesh_types["negative particle"] = particle_mesh
        submesh_types["positive particle"] = particle_mesh 
    Call_Age = RioCallback()  # define callback
    
    # update 231208 - try 3 times until give up 
    i_run_try = 0
    str_err = "Initialize only"
    while i_run_try<3:
        try:
            if i_run_try < 2:
                Simnew = pb.Simulation(
                    Model_new,
                    experiment = Experiment_Long, 
                    parameter_values=Para_update, 
                    solver = pb.CasadiSolver(),
                    var_pts = var_pts,
                    submesh_types=submesh_types )
                Sol_new = Simnew.solve(
                    calc_esoh=False,
                    save_at_cycles = Update_Cycles,
                    callbacks=Call_Age)
                if Call_Age.success == False:
                    raise Experiment_error_infeasible("Self detect")
            else:   # i_run_try = 2, 
                Simnew = pb.Simulation(
                    Model_new,
                    experiment = Experiment_Long, 
                    parameter_values=Para_update, 
                    solver = pb.CasadiSolver(return_solution_if_failed_early=True),
                    var_pts = var_pts,
                    submesh_types=submesh_types )
                Sol_new = Simnew.solve(
                    calc_esoh=False,
                    # save_at_cycles = Update_Cycles,
                    callbacks=Call_Age)
                Succeed_AGE_cycs = len(Sol_new.cycles)
                if Succeed_AGE_cycs < Update_Cycles:
                    str_err = f"Partially succeed to run the ageing set for {Succeed_AGE_cycs} cycles the {i_run_try+1}th time"
                    print(str_err)
                else: 
                    str_err = f"Fully succeed to run the ageing set for {Succeed_AGE_cycs} cycles the {i_run_try+1}th time"
                    print(str_err)
        except (
            pb.expression_tree.exceptions.ModelError,
            pb.expression_tree.exceptions.SolverError
            ) as e:
            i_run_try += 1
            Sol_new = "Model error or solver error"
            str_err = f"Fail to run the ageing set due to {Sol_new} for the {i_run_try}th time"
            print(str_err)
            DeBug_List = [ Para_update, Update_Cycles, dict_short, str_err ]
        except Experiment_error_infeasible as custom_error:
            i_run_try += 1
            Sol_new = "Experiment error or infeasible"
            DeBug_List = [
                Model, Model_new, Call_Age, Simnew, Sol, Sol_new, 
                Para_update, Experiment_Long, Update_Cycles,
                Age_T_in_K ,mesh_list,submesh_strech, var_pts,
                submesh_types,  list_short, dict_short, 
            ]
            str_err= f"{Sol_new}: {custom_error} for ageing set for the {i_run_try}th time"
            print(str_err)
            DeBug_List = [ Para_update, Update_Cycles, dict_short, str_err ]
        else:
            i_run_try += 1
            # add 230221 - update 230317 try to access Throughput capacity more than once
            i_try = 0
            while i_try<3:
                try:
                    getSth2 = Sol_new['Throughput capacity [A.h]'].entries[-1]
                except:
                    i_try += 1
                    print(f"Fail to read Throughput capacity for the {i_try}th time")
                else:
                    break
            # update 240110
            if isinstance(Sol, pb.solvers.solution.Solution):
                i_try = 0
                while i_try<3:
                    try:
                        getSth = Sol['Throughput capacity [A.h]'].entries[-1]
                    except:
                        i_try += 1
                        print(f"Fail to read Throughput capacity for the {i_try}th time")
                    else:
                        break
            elif isinstance(Sol, list):
                print("Must be the first run after restart, already have getSth")
            else:
                print("!! Big problem, Sol here is neither solution or list")
            # update 23-11-16 change method to get throughput capacity:
            # # (1) add old ones, as before; (2): change values of the last one
            Sol_new['Throughput capacity [A.h]'].entries += getSth
            cyc_number = len(Sol_new.cycles)
            if not Update_Cycles == 1: # the solution is imcomplete in this case
                thr_1st = np.trapz(
                    abs(Sol_new.cycles[0]["Current [A]"].entries), 
                    Sol_new.cycles[0]["Time [h]"].entries) # in A.h
                thr_end = np.trapz(
                    abs(Sol_new.cycles[-1]["Current [A]"].entries), 
                    Sol_new.cycles[-1]["Time [h]"].entries) # in A.h
                thr_tot = (thr_1st+thr_end) / 2 * cyc_number
            else:
                thr_tot = abs(
                    Sol_new['Throughput capacity [A.h]'].entries[-1] - 
                    Sol_new['Throughput capacity [A.h]'].entries[0]  )
            Sol_new['Throughput capacity [A.h]'].entries[-1] = getSth + thr_tot # only the last one is true
            DeBug_List = [ Para_update, Update_Cycles, dict_short, str_err ]
            print(f"Succeed to run the ageing set for {cyc_number} cycles the {i_run_try}th time")
            break # terminate the loop of trying to solve if you can get here. 
    
    # print("Solved this model in {}".format(ModelTimer.time()))
    Result_list = [Model_new, Sol_new,Call_Age,DeBug_List]
    return Result_list





























