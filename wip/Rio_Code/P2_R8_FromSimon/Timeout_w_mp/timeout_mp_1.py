import pybamm as pb;import pandas as pd   ;import numpy as np;import os;import matplotlib.pyplot as plt;import os;#import imageio
from scipy.io import savemat,loadmat;from pybamm import constants,exp;import matplotlib as mpl; fs=17; # or we can set import matplotlib.pyplot as plt then say 'mpl.rc...'
import openpyxl
import traceback
from multiprocessing import Pool   
import scipy.optimize
import random;import time, signal
fs=17;
font = {'family' : 'DejaVu Sans','size'   : fs}
mpl.rc('font', **font)

import sys  
import sys  
str_path_0 = os.path.abspath(os.path.join(pb.__path__[0],'..'))
str_path_1 = os.path.abspath(os.path.join(str_path_0,"wip\Rio_Code\Fun_P2"))
sys.path.append(str_path_1) 
from Fun_P2 import *
pb.set_logging_level("INFO")
########################  Input  ########################
# all values here must be a list, even it is a single object
Para_dict_All = {
   "Total ageing cycles":[4,],
   "Ageing cycles between RPT":[2,],
   "Update cycles for ageing": [2,],
   "Cycles within RPT":[1,],
   "Ageing temperature":[25,40],
   "RPT temperature":[25,],
   "Mesh list":[[5,5,5,30,20],],   # Simon uses 30
   "Para_Set":[ "OKane2023",], # Li2023_Coupled
   "Model option":[
         {
            #"calculate discharge energy":"true",
            "thermal": "lumped",
            "particle": "Fickian diffusion",          
            "SEI":"interstitial-diffusion limited",   
            "SEI on cracks":"true",  
            "SEI film resistance":"distributed",          
            "SEI porosity change":"true",      
            "particle mechanics":("swelling and cracking", "swelling only"), 
            "loss of active material":"stress-driven", 
            "lithium plating":"partially reversible"
         },
         ],
   "Inner SEI reaction proportion":[0.5,],
   "Ratio of lithium moles to SEI moles":[2,], # I have always been using 1 for solvent consumption model
   "Initial inner SEI thickness [m]":[2.5E-9,],
   "Initial outer SEI thickness [m]":[2.5E-9,],
   # Solvent consumption sub-model
   "Initial electrolyte excessive amount ratio":[ 0.5,1.2], # set to <1 for DryOut=Off 
   "Current solvent concentration in the reservoir [mol.m-3]":[4541.0,],
   "Current electrolyte concentration in the reservoir [mol.m-3]":[1000,],
   "Ratio of Li-ion concentration change in electrolyte consider solvent consumption":[1.0,],
   # DFN parameter
   "Upper voltage cut-off [V]":[4.21,],
   "Lower voltage cut-off [V]":[2.49,],

   # interstitial-diffusion limited
   'Inner SEI lithium interstitial diffusivity [m2.s-1]':[1e-19,],    
   'Lithium interstitial reference concentration [mol.m-3]':[15,],
   # ec-reaction limited
   'EC diffusivity [m2.s-1]':[1e-22,],
   'SEI kinetic rate constant [m.s-1]':[1e-12,], 
   'EC initial concentration in electrolyte [mol.m-3]':[4541.0,],
   'Typical EC concentration in electrolyte [mol.m-3]':[4541.0,], # Mark Ruihe change, act as an initial value here
   # LiP and coupling with SEI:
   "Dead lithium decay constant [s-1]":[ 1e-6,],            # default: 1e-6
   'Lithium plating kinetic rate constant [m.s-1]':[1E-10], # default: 1e-9
   # Crack model
   "Negative electrode LAM constant proportional term [s-1]":[ 1e-9], # default: 2.7778e-07
   "Positive electrode LAM constant proportional term [s-1]":[ 1e-9,], # default: 2.7778e-07
   # make it simple for now,], but may want to have T dependency in the future
   "Negative electrode cracking rate":[ 1e-22,],   # default: function, ~3.9e-20
   #"Positive electrode cracking rate":[ 1e-22,],   # default: function, ~3.9e-20
   #"Negative electrode volume change":[ 0.0,],
   #"Positive electrode volume change":[ 0.0,],
   #"Initial Neg SOC":[0.850],    #list(np.linspace(0.84,0.90,6)),
   #"Initial Pos SOC":[0.2705], # list(np.linspace(0.22,0.27,6)),
}
Para_dict_list = []
recursive_scan(Para_dict_list,Para_dict_All, list(Para_dict_All.keys()), {})
print(f"Total scan case is {len(Para_dict_list)}")

BasicPath = 'D:/OneDrive - Imperial College London/SimDataSave/P2_R8_FromSimon'; 
#BasicPath=os.getcwd()
Target  = '/timeout_mp_1/'
if not os.path.exists(BasicPath + Target):
   os.mkdir(BasicPath + Target);
book_name_xlsx = 'timeout_mp_1.xlsx';
sheet_name_xlsx = 'Results';
Path_pack = [BasicPath,Target,book_name_xlsx,sheet_name_xlsx,];

# define experiments-1 and output keys
V_max = 4.2;        
V_min = 2.5; 
charge_time_mins = 60 * 4.86491/5
exp_AGE_text = [(
    "Charge at 1C for 5 seconds",
    "Discharge at 1C for 5 seconds", 
    ),  ]
# step index for ageing
step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2;

exp_RPT_text = [ (
    "Discharge at 1C for 5 seconds", 
    "Rest for 2 seconds", 
    "Charge at 1C for 5 seconds",
    ) ]
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
    "CDsta Negative electrode stoichiometry", 
    "CDend Negative electrode stoichiometry",
    "CDsta Positive electrode stoichiometry", 
    "CDend Positive electrode stoichiometry",
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

# define global index and dict for all experiment data - prepare for read!
Exp_Path = [
    "Expt 1 - Si-based Degradation/",
    "Expt 2,2 - C-based Degradation 2/",
    "Expt 3 - Cathode Degradation and Li-Plating/",
    "Expt 4 - Drive Cycle Aging (Control)/",
    "Expt 5 - Standard Cycle Aging (Control)/",]
Exp_head = [
    "Expt 1",
    "Expt 2,2",
    "Expt 3",
    "Expt 4",
    "Expt 5",]
Exp_1_Cell = ["A","B","J","D","E","F","K","L","M"];
Exp_1_Temp = {
    "A":"10","B":"10","J":"10",
    "D":"25","E":"25","F":"25",
    "K":"40","L":"40","M":"40",}
Temp_Cell_Exp_1 = {
    "10":["A","B","J"],
    "25":["D","E","F"],
    "40":["K","L","M"],}
Exp_2_Cell = ["A","B","C","D","E","F"];
Exp_2_Temp = {
    "A":"10","B":"10",
    "C":"25","D":"25",
    "E":"40","F":"40",}
Temp_Cell_Exp_2 = {
    "10":["A","B"],
    "25":["C","D"],
    "40":["E","F"],}
Exp_3_Cell = ["A","B","C","D","E","F","G","H","I"];
Exp_3_Temp = {
    "A":"10","B":"10","C":"10",
    "D":"25","E":"25","F":"25",
    "G":"40","H":"40","I":"40"}
Temp_Cell_Exp_3 = {
    "10":["A","B","C"],
    "25":["D","E","F"],
    "40":["G","H","I"],}
Exp_4_Cell = ["A","B","C","D","E","F","G","H"];
Exp_4_Temp = {
    "A":"10","B":"10","C":"10",
    "D":"25","E":"25",
    "F":"40","G":"40","H":"40",}
Temp_Cell_Exp_4 = {
    "10":["A","B","C"],
    "25":["D","E",],
    "40":["F","G","H"],}
Exp_5_Cell = ["A","B","C","D","E","F","G","H"];
Exp_5_Temp = {
    "A":"10","B":"10","C":"10",
    "D":"25","E":"25",
    "F":"40","G":"40","H":"40",}
Temp_Cell_Exp_5 = {
    "10":["A","B","C"],
    "25":["D","E",],
    "40":["F","G","H"],}
Exp_All_Cell  = [Exp_1_Cell,Exp_2_Cell,Exp_3_Cell,Exp_4_Cell,Exp_5_Cell]
Exp_Temp_Cell = [Exp_1_Temp,Exp_2_Temp,Exp_3_Temp,Exp_4_Temp,Exp_5_Temp]
Temp_Cell_Exp_All = [Temp_Cell_Exp_1,Temp_Cell_Exp_2,Temp_Cell_Exp_3,Temp_Cell_Exp_4,Temp_Cell_Exp_5]
Mark_Cell_All = [
    {
    "A":"o","B":">","J":"v",
    "D":"o","E":">","F":"v",
    "K":"o","L":">","M":"v",},
    {
    "A":"10","B":"10",
    "C":"25","D":"25",
    "E":"40","F":"40",},
    {
    "A":"o","B":">","C":"v",
    "D":"o","E":">","F":"v",
    "G":"o","H":">","I":"v",},
    {
    "A":"o","B":">","C":"v",
    "D":"o","E":">",
    "F":"o","G":">","H":"v",},
    {
    "A":"o","B":">","C":"v",
    "D":"o","E":">",
    "F":"o","G":">","H":"v",}]
Color_Cell_All = [
    {
    "A":[2/255, 3/255, 226/255,0.7],"B":[2/255, 3/255, 226/255,0.7],"J":[2/255, 3/255, 226/255,0.7],
    "D":[0, 0, 0,0.7],"E":[0, 0, 0,0.7],"F":[0, 0, 0,0.7],
    "K":[1,0,0,0.4],"L":[1,0,0,0.4],"M":[1,0,0,0.4],},
    {
    "A":[2/255, 3/255, 226/255,0.7],"B":[2/255, 3/255, 226/255,0.7],
    "D":[0, 0, 0,0.7],"C":[0, 0, 0,0.7],
    "E":[1,0,0,0.4],"F":[1,0,0,0.4],},
    {
    "A":[2/255, 3/255, 226/255,0.7],"B":[2/255, 3/255, 226/255,0.7],"C":[2/255, 3/255, 226/255,0.7],
    "D":[0, 0, 0,0.7],"E":[0, 0, 0,0.7],"F":[0, 0, 0,0.7],
    "G":[1,0,0,0.4],"H":[1,0,0,0.4],"I":[1,0,0,0.4],},
    {
    "A":[2/255, 3/255, 226/255,0.7],"B":[2/255, 3/255, 226/255,0.7],"C":[2/255, 3/255, 226/255,0.7],
    "D":[0, 0, 0,0.7],"E":[0, 0, 0,0.7],
    "F":[1,0,0,0.4],"G":[1,0,0,0.4],"H":[1,0,0,0.4],},
    {
    "A":[2/255, 3/255, 226/255,0.7],"B":[2/255, 3/255, 226/255,0.7],"C":[2/255, 3/255, 226/255,0.7],
    "D":[0, 0, 0,0.7],"E":[0, 0, 0,0.7],
    "F":[1,0,0,0.4],"G":[1,0,0,0.4],"H":[1,0,0,0.4],}]

# Global path 
Path_NiallDMA = "D:/OneDrive - Imperial College London/SimDataSave/InputData/"
index_exp = 1  # index for experiment set, range from 1~5
Temp_Cell_Exp = Temp_Cell_Exp_All[index_exp-1] # select a cell list through temperature
Exp_Any_AllData = Read_Exp(
    Path_NiallDMA,Exp_All_Cell[index_exp-1],Exp_Path,Exp_head,Exp_Temp_Cell[index_exp-1],index_exp-1)

# Write para - 1st round:
head_list = list(Para_dict_list[0].keys())
head_list.insert(0,"Index");
head_list.insert(1,"Dry out");
head_list.extend([ "exp_AGE_text", "exp_RPT_text",
   "Cap Loss","LLI to LiP",
   "LLI to SEI","LLI to sei-on-cracks",
   "LAM to Neg","LAM to Pos",
   "Vol_Elely_Tot Final", "Vol_Elely_JR Final","Width Final","Error"])
Values_1 = [head_list,];
index_list = np.arange(1,len(Para_dict_list)+1,1)
for Scan_i,Para_dict_i in zip(index_list,Para_dict_list):
    value_list_temp = list(Para_dict_i.values())
    values = []
    for value_list_temp_i in value_list_temp:
        values.append(str(value_list_temp_i))
    values.insert(0,str(Scan_i))
    Values_1.append(values)
write_excel_xlsx(
    BasicPath + Target+book_name_xlsx, 
    sheet_name_xlsx, Values_1)   

fs = 10

#################### need to change:  ####################
####################################################
# task: save the last solution if fail in the middle, so that we can have a closer look
# build a big class of my own setting, including results
class TaskConfig:
    def __init__(self,index_xlsx, Para_dict_i,   Path_pack ,  fs,
        keys_all,   exp_text_list, exp_index_pack , 
        Exp_Any_AllData,Temp_Cell_Exp, 
        Plot_Exp,Timeout,Return_Sol):  # = true or false

        self.Scan_i = int(index_xlsx)
        CyclePack,Para_0 = Para_init(Para_dict_i) # 
        [Total_Cycles,Cycle_bt_RPT,Update_Cycles,RPT_Cycles,
            Temper_i,Temper_RPT,mesh_list,submesh_strech,model_options] = CyclePack;
        self.Total_Cycles = Total_Cycles
        self.Cycle_bt_RPT=Cycle_bt_RPT
        self.Update_Cycles=Update_Cycles
        self.RPT_Cycles=RPT_Cycles
        self.Temper_i=Temper_i
        self.Temper_RPT=Temper_RPT
        self.mesh_list=mesh_list
        self.submesh_strech=submesh_strech
        self.model_options=model_options
        [BasicPath,Target,book_name_xlsx,sheet_name_xlsx,] = Path_pack
        self.BasicPath = BasicPath
        self.Target = Target
        self.book_name_xlsx=book_name_xlsx
        self.sheet_name_xlsx=sheet_name_xlsx
        self.fs = fs
        [keys_all_RPT,keys_all_AGE] = keys_all
        self.keys_all_RPT=keys_all_RPT
        self.keys_all_AGE=keys_all_AGE
        [exp_AGE_text, exp_RPT_text,] = exp_text_list;
        str_exp_RPT_text  = str(exp_RPT_text);
        str_exp_AGE_text  = str(exp_AGE_text);
        self.str_exp_AGE_text = str_exp_AGE_text
        self.str_exp_RPT_text=str_exp_RPT_text
        [cycle_no,step_AGE_CD,step_AGE_CC,step_AGE_CV,
            step_RPT_CD,step_RPT_RE , step_RPT_CC ] = exp_index_pack;
        self.cycle_no = cycle_no
        self.step_AGE_CD=step_AGE_CD
        self.step_AGE_CC=step_AGE_CC
        self.step_AGE_CV=step_AGE_CV
        self.step_RPT_CD=step_RPT_CD
        self.step_RPT_RE=step_RPT_RE
        self.step_RPT_CC=step_RPT_CC
        self.Exp_Any_AllData = Exp_Any_AllData
        self.Temp_Cell_Exp = Temp_Cell_Exp
        self.Plot_Exp = Plot_Exp
        self.Timeout = Timeout
        self.Return_Sol = Return_Sol
        # define experiment
        Experiment_Long   = pb.Experiment( exp_AGE_text * Update_Cycles  )  
        Experiment_RPT    = pb.Experiment( exp_RPT_text * RPT_Cycles     ) 
        Experiment_Breakin= pb.Experiment( exp_RPT_text * RPT_Cycles     )

        #####  index definition ######################
        Small_Loop =  int(Cycle_bt_RPT/Update_Cycles);   
        SaveTimes = int(Total_Cycles/Cycle_bt_RPT);
        self.Experiment_Long=Experiment_Long
        self.Experiment_RPT=Experiment_RPT
        self.Experiment_Breakin=Experiment_Breakin
        self.Small_Loop=Small_Loop
        self.SaveTimes=SaveTimes
        # update 220924: merge DryOut and Int_ElelyExces_Ratio
        temp_Int_ElelyExces_Ratio =  Para_0["Initial electrolyte excessive amount ratio"] 
        ce_EC_0 = Para_0['EC initial concentration in electrolyte [mol.m-3]'] # used to calculate ce_EC_All
        if temp_Int_ElelyExces_Ratio < 1:
            Int_ElelyExces_Ratio = -1;
            DryOut = "Off";
        else:
            Int_ElelyExces_Ratio = temp_Int_ElelyExces_Ratio;
            DryOut = "On";  
        if DryOut == "On":  
            mdic_dry,Para_0 = Initialize_mdic_dry(Para_0,Int_ElelyExces_Ratio)
        else:
            mdic_dry ={}
        self.DryOut=DryOut
        print(f"Scan {self.Scan_i}: DryOut = {self.DryOut}")
        self.Para_0 = Para_0            # most important, the parameter set!!!
        self.ce_EC_0=ce_EC_0
        self.Int_ElelyExces_Ratio=Int_ElelyExces_Ratio
        self.mdic_dry=mdic_dry
        self.Timeout_text = 'I timed out'
        # define results:
        self.Sol_RPT = [];  self.Sol_AGE = [];
        my_dict_RPT = {};my_dict_AGE = {}; 
        for keys in keys_all_RPT:
            for key in keys:
                my_dict_RPT[key]=[];
        for keys in keys_all_AGE:
            for key in keys:
                my_dict_AGE[key]=[];
        my_dict_RPT["Cycle_RPT"] = []; my_dict_AGE["Cycle_AGE"] = []; 
        Cyc_Update_Index     =[]
        self.midc_merge_all=[]
        self.my_dict_RPT=my_dict_RPT
        self.my_dict_AGE=my_dict_AGE
        self.Cyc_Update_Index=Cyc_Update_Index
        # output for break-in cycle:
        self.Model_0=[];
        self.Sol_0=[];
        self.Call_Breakin=[]

        self.Flag_Breakin=True
        self.str_error_Breakin=[]
        self.Flag_AGE = True; 
        self.str_error_AGE_final = "Empty"
        self.str_error_RPT = "Empty"
        self.Call_Age=[]
        self.Call_RPT=[]

########### Most important new feature: task_configs ###########
################################################################
task_configs = [] # get all configurations
for index_i, Para_dict_i in zip(index_list,Para_dict_list):
    task_configs.append(TaskConfig(
        index_i    ,    Para_dict_i,   Path_pack, fs,
        keys_all,   exp_text_list, exp_index_pack,
        Exp_Any_AllData,Temp_Cell_Exp,
        False,True,True,  ))
# define the model and run break-in cycle - 
# input parameter: model_options, Experiment_Breakin, Para_0, mesh_list, submesh_strech
# output: Sol_0 , Model_0, Call_Breakin
def run_breakin_single(task_config):
    task_config.Model_0 = pb.lithium_ion.DFN(options=task_config.model_options ) #
    # update 220926 - add diffusivity and conductivity as variables:
    c_e = task_config.Model_0.variables["Electrolyte concentration [mol.m-3]"]
    T = task_config.Model_0.variables["Cell temperature [K]"]
    D_e = task_config.Para_0["Electrolyte diffusivity [m2.s-1]"]
    sigma_e = task_config.Para_0["Electrolyte conductivity [S.m-1]"]
    task_config.Model_0.variables["Electrolyte diffusivity [m2.s-1]"] = D_e(c_e, T)
    task_config.Model_0.variables["Electrolyte conductivity [S.m-1]"] = sigma_e(c_e, T)
    var = pb.standard_spatial_vars  
    var_pts = {
        var.x_n: int(task_config.mesh_list[0]),  
        var.x_s: int(task_config.mesh_list[1]),  
        var.x_p: int(task_config.mesh_list[2]),  
        var.r_n: int(task_config.mesh_list[3]),  
        var.r_p: int(task_config.mesh_list[4]),  }       
    submesh_types = task_config.Model_0.default_submesh_types
    if task_config.submesh_strech == "nan":
        pass
    else:
        particle_mesh = pb.MeshGenerator(
            pb.Exponential1DSubMesh, 
            submesh_params={"side": "right", "stretch": task_config.submesh_strech})
        submesh_types["negative particle"] = particle_mesh
        submesh_types["positive particle"] = particle_mesh 
    Sim_0    = pb.Simulation(
        task_config.Model_0,        experiment = task_config.Experiment_Breakin,
        parameter_values = task_config.Para_0,
        solver = pb.CasadiSolver(),
        var_pts=var_pts,
        submesh_types=submesh_types) #mode="safe"
    task_config.Call_Breakin = RioCallback()
    try:
        task_config.Sol_0    = Sim_0.solve(calc_esoh=False,callbacks=task_config.Call_Breakin)
    except (
        pb.expression_tree.exceptions.ModelError,
        pb.expression_tree.exceptions.SolverError
        ) as e:
        task_config.Sol_0 = "Model error or solver error"
    else:
        pass
    return task_config

# ,Timelimit,check_interval, # task_configs
def run_Breakin_round(configs,Timelimit,check_interval):
    worker_pool = Pool(processes=len(configs))  
    async_results = []
    # submit jobs and start multiprosessing 
    for config in configs:
        async_result = worker_pool.apply_async(run_breakin_single, (config,))
        async_results.append(async_result)
    # stop submitting jobs
    worker_pool.close()
    # set time limit = Timelimit; check every check_interval period 
    Timelimit_series = np.arange(check_interval,Timelimit+check_interval,check_interval)
    for check in Timelimit_series:
        time.sleep(check)
        print(f"Checking at {check} second")
        all_ready = all([async_result.ready() for async_result in async_results])
        if all_ready:
            break
            
    # either all finish, or timeout
    configs = [] # clear configs
    for async_result in async_results:
        if async_result.ready():
            task_result = async_result.get()
            configs.append(task_result)
        else:
            print('Task timeout')
    
    # close the pool to release resources
    worker_pool.terminate()
    
    return configs

# really run the cases with both timeout and multiprocess enabled!
Timelimit = int(60*20); # set timeout: unit: second
check_interval = int(4*60) # check interval unit: second
task_configs = run_Breakin_round(task_configs,Timelimit,check_interval)

print("Reach end of the file")