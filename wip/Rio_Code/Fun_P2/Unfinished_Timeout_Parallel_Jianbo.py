import pybamm as pb;import pandas as pd   ;import numpy as np;import os;import matplotlib.pyplot as plt;import os;#import imageio
from scipy.io import savemat,loadmat;from pybamm import constants,exp;import matplotlib as mpl; fs=17; # or we can set import matplotlib.pyplot as plt then say 'mpl.rc...'
import openpyxl
import traceback
import multiprocessing
import time
from multiprocessing import Pool   
import scipy.optimize
import random;import time, signal
fs=17;
font = {'family' : 'DejaVu Sans','size'   : fs}
mpl.rc('font', **font)
from Fun_P2_Patrick import (
    recursive_scan,
    write_excel_xlsx,
    Run_P2_till_Fail,
    Run_P2_Opt_Timeout,
    Para_init,
    Initialize_mdic_dry,
    RioCallback,
    GetSol_dict,
)
""" 
try写在子任务里，子任务是有返回值的，我在例子里写的是字符串，改成你定义的结果就好了。返回结果有三种情况：
1. 正常运行，返回一个正常的结果
2. 中途出错，在except里返回一个包含错误信息的结果
3. 运行超时，没有结果
然后主程序里，根据这三种结果分别处理：
1. 子任务正常返回，就生成下一轮的参数，准备下一次仿真
2. 子任务运行异常，调整参数，准备下一次仿真（或者别的处理）
3. 子任务超时，标记为下次不再仿真 
"""
########################  Input  ########################
# all values here must be a list, even it is a single object
Para_dict_All = {
   "Total ageing cycles":[9,],
   "Ageing cycles between RPT":[3,],
   "Update cycles for ageing":[3,],
   "Cycles within RPT":[1,],
   "Ageing temperature":[25,],
   "RPT temperature":[25,],
   "Particle mesh points":[30,],   # Simon uses 30
   #"Exponential mesh stretch":[1.0],
   "Para_Set":[ "Li2023_Coupled",],
   "Model option":[
         {
            "calculate discharge energy":"true",
            "particle": "Fickian diffusion",          
            "SEI":"interstitial-diffusion limited",   
            "SEI on cracks":"true",  
            "SEI film resistance":"distributed",          
            "SEI porosity change":"true",      
            "particle mechanics":("swelling and cracking", "swelling only"), 
            "loss of active material":"stress-driven", 
            "lithium plating":"partially reversible"      },
         ],
   "Inner SEI reaction proportion":[0.5,],
   "Ratio of lithium moles to SEI moles":[2,], # I have always been using 1 for solvent consumption model
   "Initial inner SEI thickness [m]":[2.5E-9,],
   "Initial outer SEI thickness [m]":[2.5E-9,],
   "SEI growth activation energy [J.mol-1]":[38000,],
   # Solvent consumption sub-model
   "Initial electrolyte excessive amount ratio":[ 2], # set to <1 for DryOut=Off 
   "Current solvent concentration in the reservoir [mol.m-3]":[4541.0,],
   "Current electrolyte concentration in the reservoir [mol.m-3]":[1000,],
   "Ratio of Li-ion concentration change in electrolyte consider solvent consumption":[1.0,],
   # DFN parameter
   "Upper voltage cut-off [V]":[4.21,],
   "Lower voltage cut-off [V]":[2.49,],

   # interstitial-diffusion limited
   'Inner SEI lithium interstitial diffusivity [m2.s-1]':[1e-22],    
   'Lithium interstitial reference concentration [mol.m-3]':[15,],
   # ec-reaction limited
   'EC diffusivity [m2.s-1]':[1e-22,],
   'SEI kinetic rate constant [m.s-1]':[1e-12,], 
   'EC initial concentration in electrolyte [mol.m-3]':[4541.0,],
   'Typical EC concentration in electrolyte [mol.m-3]':[4541.0,], # Mark Ruihe change, act as an initial value here
   # LiP and coupling with SEI:
   "Dead lithium decay constant [s-1]":[ 1e-6,],            # default: 1e-6
   'Lithium plating kinetic rate constant [m.s-1]':[1E-12], # default: 1e-9
   # Crack model
   "Negative electrode LAM constant proportional term [s-1]":[ 2.7778e-10,], # default: 2.7778e-07
   "Positive electrode LAM constant proportional term [s-1]":[ 2.7778e-10,], # default: 2.7778e-07
   # make it simple for now,], but may want to have T dependency in the future
   "Negative electrode cracking rate":[ 3.9e-22,],   # default: function, ~3.9e-20
   "Positive electrode cracking rate":[ 3.9e-22,],   # default: function, ~3.9e-20
   "Negative electrode volume change":[ 0.0,],
   "Positive electrode volume change":[ 0.0,],
   "Initial Neg SOC":[0.850],    #list(np.linspace(0.84,0.90,6)),
   "Initial Pos SOC":[0.2705], # list(np.linspace(0.22,0.27,6)),
}
Para_dict_list = []
recursive_scan(Para_dict_list,Para_dict_All, list(Para_dict_All.keys()), {})
print(f"Total scan case is {len(Para_dict_list)}")

#BasicPath = 'D:/OneDrive - Imperial College London/SimDataSave/P2R7'; 
BasicPath=os.getcwd()
Target  = '/Test_timeout_parallel/'
if not os.path.exists(BasicPath + Target):
   os.mkdir(BasicPath + Target);
book_name_xlsx = 'Test_timeout_parallel.xlsx';
sheet_name_xlsx = 'Results';
Path_pack = [BasicPath,Target,book_name_xlsx,sheet_name_xlsx,];

# define experiments and output keys
V_max = 4.2;        
V_min = 2.5; 
exp_AGE_text = [(f"Discharge at 1 C until {V_min} V", 
        f"Charge at 0.3 C until {V_max} V", 
        f"Hold at {V_max} V until C/100"),  ]
# step index for ageing
step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2;

exp_RPT_text = [ (f"Discharge at 0.1C until {V_min} V",  
        "Rest for 1 hours",  
        f"Charge at 0.1C until {V_max} V" ) ]
# step index for RPT
step_RPT_CD = 0;  step_RPT_RE =1;   step_RPT_CC = 2;  

exp_text_list = [exp_AGE_text, exp_RPT_text,];
cycle_no = -1; 
exp_index_pack = [cycle_no,step_AGE_CD,step_AGE_CC,step_AGE_CV,
   step_RPT_CD,step_RPT_RE , step_RPT_CC ];
########################  Output  ########################
keys_loc_RPT = [ # MAY WANT TO SELECT AGEING CYCLE later
    # Default output:
    "x [m]",
    "x_n [m]",
    "x_s [m]",
    "x_p [m]",
    # default: end; 
    "CCend Porosity",
    "CCend Negative electrode interfacial current density [A.m-2]",
    "CCend Electrolyte potential [V]",
    "CCend Electrolyte concentration [mol.m-3]",
    "CCend Negative electrode reaction overpotential [V]",
    "CCend Negative particle surface concentration [mol.m-3]",
    "CCend Negative electrode roughness ratio",
    "CCend Total SEI on cracks thickness [m]",

    "CDend Porosity",
    "CDend Negative electrode interfacial current density [A.m-2]",
    "CDend Electrolyte potential [V]",
    "CDend Electrolyte concentration [mol.m-3]",
    "CDend Negative electrode reaction overpotential [V]",
    "CDend Negative particle surface concentration [mol.m-3]",
    "CDend Negative electrode roughness ratio",
    "CDend Total SEI on cracks thickness [m]",
    "REend Total SEI on cracks thickness [m]",
]
keys_tim_RPT = [
    # default: CD
    "CD Time [h]",
    "CD Terminal voltage [V]",
    "RE Terminal voltage [V]",
]
keys_cyc_RPT = [   # default: CDend
    "Discharge capacity [A.h]",
    "CDend Loss of capacity to lithium plating [A.h]",
    "CDend Loss of capacity to SEI [A.h]",
    "CDend Loss of capacity to SEI on cracks [A.h]",
    "CDend X-averaged total SEI on cracks thickness [m]",
    "CDend X-averaged negative electrode roughness ratio",
    "CDend Local ECM resistance [Ohm]",
    "CDsta Negative electrode SOC", 
    "CDend Negative electrode SOC",
    "CDsta Positive electrode SOC", 
    "CDend Positive electrode SOC",
    "CDend Negative electrode capacity [A.h]",
    "CDend Positive electrode capacity [A.h]",
    "CDend Throughput capacity [A.h]",
]

keys_loc_AGE = [ # MAY WANT TO SELECT AGEING CYCLE later
    # Default output:
    "x [m]",
    "x_n [m]",
    "x_s [m]",
    "x_p [m]",
    # default: end; 
    "CCend Porosity",
    "CCend Negative electrode interfacial current density [A.m-2]",
    "CCend Electrolyte potential [V]",
    "CCend Electrolyte concentration [mol.m-3]",
    "CCend Negative electrode reaction overpotential [V]",
    "CCend Negative particle surface concentration [mol.m-3]",
    "CCend Negative electrode roughness ratio",
    "CCend Total SEI on cracks thickness [m]",

    "CDend Porosity",
    "CDend Negative electrode interfacial current density [A.m-2]",
    "CDend Electrolyte potential [V]",
    "CDend Electrolyte concentration [mol.m-3]",
    "CDend Negative electrode reaction overpotential [V]",
    "CDend Negative particle surface concentration [mol.m-3]",
    "CDend Negative electrode roughness ratio",
    "CDend Total SEI on cracks thickness [m]",
    "CDend Electrolyte diffusivity [m2.s-1]",
    "CDend Electrolyte conductivity [S.m-1]",
]
keys_tim_AGE = [];
keys_cyc_AGE = [];
keys_all_RPT = [keys_loc_RPT,keys_tim_RPT,keys_cyc_RPT];
keys_all_AGE = [keys_loc_AGE,keys_tim_AGE,keys_cyc_AGE];
keys_all = [keys_all_RPT,keys_all_AGE];

# Write the head for excel file:
head_list = list(Para_dict_list[0].keys())
head_list.insert(0,"Index");
head_list.insert(1,"Dry out");
head_list.extend([ "exp_AGE_text", "exp_RPT_text",
   "Cap Loss","LLI to LiP",
   "LLI to SEI","LLI to sei-on-cracks",
   "LAM to Neg","LAM to Pos",
   "Vol_Elely_Tot Final", "Vol_Elely_JR Final","Width Final","Error"])
write_excel_xlsx(BasicPath + Target+book_name_xlsx, sheet_name_xlsx, [head_list])

Niall_data = loadmat( 'Extracted_all_cell.mat')

##########################################################
####  Part-1: Prepare for break-in cycle    ##############
##########################################################
# state clearly input and output for results? 
index_list = np.arange(1,len(Para_dict_list)+1,1)
CyclePack_All = []; Para_0_All = []; 
ce_EC_0_All   = []; DryOut_All = [];
Dict_RPT_All =[]; Dict_AGE_All=[]; Dict_dry_All =[];
Cyc_Update_Index_All =[]; Flag_timeout_All =[];
for Scan_i,Para_dict_i in zip(index_list,Para_dict_list):
    CyclePack,Para_0 = Para_init(Para_dict_i)
    temp_Int_ElelyExces_Ratio =  Para_0["Initial electrolyte excessive amount ratio"] 
    ce_EC_0 = Para_0['EC initial concentration in electrolyte [mol.m-3]'] 
    Flag_timeout = False
    if temp_Int_ElelyExces_Ratio < 1:
        Int_ElelyExces_Ratio = -1;
        DryOut = "Off";
    else:
        Int_ElelyExces_Ratio = temp_Int_ElelyExces_Ratio;
        DryOut = "On";
    print(f"Scan {Scan_i}: DryOut = {DryOut}")
    # initialize my_dict for outputs
    my_dict_RPT = {}
    for keys in keys_all_RPT:
        for key in keys:
            my_dict_RPT[key]=[];
    my_dict_AGE = {}; 
    for keys in keys_all_AGE:
        for key in keys:
            my_dict_AGE[key]=[];
    my_dict_RPT["Cycle_RPT"] = []; my_dict_AGE["Cycle_AGE"] = []; 
    Cyc_Update_Index     =[]; 
    if DryOut == "On":  
        mdic_dry,Para_0 = Initialize_mdic_dry(Para_0,Int_ElelyExces_Ratio)
    else:
        mdic_dry ={}
    # Sum up things:
    CyclePack_All.append(CyclePack)
    Para_0_All.append(Para_0)
    ce_EC_0_All.append(ce_EC_0)
    DryOut_All.append(DryOut)
    Dict_RPT_All.append(my_dict_RPT)
    Dict_AGE_All.append(my_dict_AGE)
    Dict_dry_All.append(mdic_dry)
    Cyc_Update_Index_All.append(Cyc_Update_Index)
    Flag_timeout_All.append(Flag_timeout)

def Run_Breakin_Timeout(exp_RPT_text, Para_0, CyclePack,Scan_i):

    [Total_Cycles,Cycle_bt_RPT,Update_Cycles,RPT_Cycles,
        Temper_i,Temper_RPT,mesh_par,submesh_strech,model_options] = CyclePack
    Experiment_Breakin = pb.Experiment( exp_RPT_text * RPT_Cycles     )
    Model_0 = pb.lithium_ion.DFN(options=model_options ) #
    # update 220926 - add diffusivity and conductivity as variables:
    c_e = Model_0.variables["Electrolyte concentration [mol.m-3]"]
    T = Model_0.variables["Cell temperature [K]"]
    D_e = Para_0["Electrolyte diffusivity [m2.s-1]"]
    sigma_e = Para_0["Electrolyte conductivity [S.m-1]"]
    Model_0.variables["Electrolyte diffusivity [m2.s-1]"] = D_e(c_e, T)
    Model_0.variables["Electrolyte conductivity [S.m-1]"] = sigma_e(c_e, T)
    var = pb.standard_spatial_vars  
    var_pts = {
        var.x_n: 20,  
        var.x_s: 10,  
        var.x_p: 20,  
        var.r_n: int(mesh_par),  
        var.r_p: int(mesh_par),  }       
    submesh_types = Model_0.default_submesh_types
    if submesh_strech == "nan":
        pass
    else:
        particle_mesh = pb.MeshGenerator(
            pb.Exponential1DSubMesh, 
            submesh_params={"side": "right", "stretch": submesh_strech})
        submesh_types["negative particle"] = particle_mesh
        submesh_types["positive particle"] = particle_mesh 
    Sim_0    = pb.Simulation(
        Model_0,        experiment = Experiment_Breakin,
        parameter_values = Para_0,
        solver = pb.CasadiSolver(),
        var_pts=var_pts,
        submesh_types=submesh_types) #mode="safe"
    Call_Breakin = RioCallback()
    str_error_Breakin = "Empty";  Flag_Breakin = True
    try:
        Sol_0    = Sim_0.solve(calc_esoh=False,callbacks=Call_Breakin)
        if Call_Breakin.success == False:
            print("Break in cycle fail due to Experiment error or infeasible")
            1/0
    except (
        pb.expression_tree.exceptions.ModelError,
        pb.expression_tree.exceptions.SolverError
        ) as e:
        Sol_0 = "Model error or solver error"
        str_error_Breakin = "Model error or solver error"
        Flag_Breakin = False
    except ZeroDivisionError as e:
        Sol_0 = "Experiment error or infeasible"
        str_error_Breakin = "Experiment error or infeasible"
        Flag_Breakin = False
    else:
        print(f"Scan {Scan_i}: Finish break-in cycle")
    Result_list_breakin = [
        Model_0,Sol_0,Call_Breakin,
        Flag_Breakin,str_error_Breakin]

    return Result_list_breakin

##########################################################
############  Part-2: Do break-in cycle    ###############
##########################################################
TASK_NUM = len(Para_dict_list)
worker_pool = Pool(processes=TASK_NUM)   
# Submit all break-in cycles 
async_Breakins = [ [] for x in range(TASK_NUM)]
for Scan_i in index_list:
    i = int(Scan_i-1)
    if Flag_timeout_All[i] == False: # Only submit this job if never timeout before
        async_Breakin = worker_pool.apply_async(
            Run_Breakin_Timeout, 
            (
                exp_RPT_text,Para_0_All[i],
                CyclePack_All[i],Scan_i)
            )
        async_Breakins[i] = async_Breakin
# Stop job submission
worker_pool.close()

# Set timeout:  Note: if set too long, timeout won't apply
#     for break-in cycle: 0.5 hour, check every 1 minutes
Timelimit = int(3600*0.5); check_i = 60; 
Timecheck = np.arange(check_i,Timelimit+check_i,check_i).tolist()
for time_check in Timecheck:
    time.sleep(time_check)
    all_ready = all([async_Breakin.ready() for async_Breakin in async_Breakins])
    if all_ready:
        break
# Now a time lengh of Timelimit has passed, try to get results, even if some are not ready yet.
for Scan_i, async_Breakin in zip(index_list,async_Breakins):
    # check whehter timeout has happen:
    i = Scan_i - 1
    if Flag_timeout_All[i] == False:
        if async_Breakin.ready() and Flag_Breakin == True:
            Result_list_breakin=async_Breakin.get()
            [
                Model_0,Sol_0,Call_Breakin,
                Flag_Breakin,str_error_Breakin,
                ] = Result_list_breakin 

        else:
            print('Task %d timeout or fail during break-in' % Scan_i)
            Flag_timeout_All[i] = True
            # TODO: # Fun_process_break_in_fail



# post-process for break-in cycle
my_dict_RPT = GetSol_dict (my_dict_RPT,keys_all_RPT, Sol_0, 
    cycle_no, step_RPT_CD , step_RPT_CC , step_RPT_RE, step_AGE_CV   )
cycle_count =0; 
my_dict_RPT["Cycle_RPT"].append(cycle_count)
Cyc_Update_Index.append(cycle_count);