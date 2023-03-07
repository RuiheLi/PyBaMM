import pybamm as pb;import pandas as pd   ;import numpy as np;import os;import matplotlib.pyplot as plt;import os;#import imageio
from scipy.io import savemat,loadmat;from pybamm import constants,exp;import matplotlib as mpl; fs=17; # or we can set import matplotlib.pyplot as plt then say 'mpl.rc...'
import openpyxl
import traceback
import multiprocessing
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
   "Total ageing cycles":[6,],
   "Ageing cycles between RPT":[3,],
   "Update cycles for ageing": [3,],
   "Cycles within RPT":[1,],
   "Ageing temperature":[10,25,40],
   "RPT temperature":[25,],
   "Mesh list":[[5,5,5,60,20],],   # Simon uses 30
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

BasicPath = 'c:/Users/rl1120/OneDrive - Imperial College London/SimDataSave/P2_R8_FromSimon'; 
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
    f"Charge at 0.3 C for {charge_time_mins} minutes",
    f"Discharge at 1 C until {V_min} V", 
    ),  ]
# step index for ageing
step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2;

exp_RPT_text = [ (
    f"Charge at 1C until {V_max} V",  
    f"Hold at {V_max}V until C/100",
    f"Discharge at 0.1C until {V_min} V (5 minute period)",  
    "Rest for 1 hours (5 minute period)",  
    f"Charge at 0.1C until {V_max} V (5 minute period)",
    "Rest for 1 hours (5 minute period)",
    f"Discharge at 1C until {V_min} V",  
    f"Hold at {V_min} V until C/100",
    ) ]
# step index for RPT
step_RPT_CD = 2;  step_RPT_RE =3;   step_RPT_CC = 4;  

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
Path_NiallDMA = "c:/Users/rl1120/OneDrive - Imperial College London/SimDataSave/InputData/"
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


# scan:
index_list = np.arange(1,len(Para_dict_list)+1,1)

#################### original:  ####################
####################################################
# scan:
fs = 10
index_list = np.arange(1,len(Para_dict_list)+1,1)
midc_merge_all = [];Sol_RPT_all = [];Sol_AGE_all = [];
for index_i, Para_dict_i in zip(index_list,Para_dict_list):
    midc_merge,Sol_RPT,Sol_AGE = Run_P2_Opt_Timeout(
        index_i    ,    Para_dict_i,   Path_pack, fs,
        keys_all,   exp_text_list, exp_index_pack,
        Exp_Any_AllData,Temp_Cell_Exp,False,
        True,True,)  
    midc_merge_all.append(midc_merge)
    Sol_RPT_all.append(Sol_RPT)
    Sol_AGE_all.append(Sol_AGE)

#################### need to change:  ####################
####################################################
# task: save the last solution if fail in the middle, so that we can have a closer look

def Run_P2_Opt_Timeout(
    index_xlsx, Para_dict_i,   Path_pack ,  fs,
    keys_all,   exp_text_list, exp_index_pack , 
    Exp_Any_AllData,Temp_Cell_Exp, Plot_Exp,   # = true or false
    Timeout,Return_Sol ):

    ##########################################################
    ##############    Part-0: Log of the scripts    ##########
    ##########################################################
    # add 221205: if Timeout=='True', use Patrick's version, disable pool
    #             else, use pool to accelerate 
    # add Return_Sol, on HPC, always set to False, as it is useless, 
    # add 230221: do sol_new['Throughput capacity [A.h]'].entries += sol_old['Throughput capacity [A.h]'].entries 
    #             and for "Throughput energy [W.h]", when use Model.set_initial_conditions_from
    #             this is done inside the two functions Run_Model_Base_On_Last_Solution(_RPT)

    ##########################################################
    ##############    Part-1: Initialization    ##############
    ##########################################################
    font = {'family' : 'DejaVu Sans','size'   : fs}
    mpl.rc('font', **font)
    ModelTimer = pb.Timer()
    Scan_i = int(index_xlsx)
    print('Start Now! Scan %d.' % Scan_i)  
    Sol_RPT = [];  Sol_AGE = [];
    # pb.set_logging_level('INFO') # show more information!
    # set_start_method('fork') # from Patrick

    # Un-pack data:
    [cycle_no,step_AGE_CD,step_AGE_CC,step_AGE_CV,
        step_RPT_CD,step_RPT_RE , step_RPT_CC ] = exp_index_pack;
    [exp_AGE_text, exp_RPT_text,] = exp_text_list;
    [BasicPath,Target,book_name_xlsx,sheet_name_xlsx,] = Path_pack
    CyclePack,Para_0 = Para_init(Para_dict_i)
    [Total_Cycles,Cycle_bt_RPT,Update_Cycles,RPT_Cycles,
        Temper_i,Temper_RPT,mesh_list,submesh_strech,model_options] = CyclePack;
    [keys_all_RPT,keys_all_AGE] = keys_all
    str_exp_AGE_text  = str(exp_AGE_text);
    str_exp_RPT_text  = str(exp_RPT_text);

    # define experiment
    Experiment_Long   = pb.Experiment( exp_AGE_text * Update_Cycles  )  
    Experiment_RPT    = pb.Experiment( exp_RPT_text * RPT_Cycles     ) 
    Experiment_Breakin= pb.Experiment( exp_RPT_text * RPT_Cycles     )

    #####  index definition ######################
    Small_Loop =  int(Cycle_bt_RPT/Update_Cycles);   
    SaveTimes = int(Total_Cycles/Cycle_bt_RPT);   

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
            
    # update 220924: merge DryOut and Int_ElelyExces_Ratio
    temp_Int_ElelyExces_Ratio =  Para_0["Initial electrolyte excessive amount ratio"] 
    ce_EC_0 = Para_0['EC initial concentration in electrolyte [mol.m-3]'] # used to calculate ce_EC_All
    if temp_Int_ElelyExces_Ratio < 1:
        Int_ElelyExces_Ratio = -1;
        DryOut = "Off";
    else:
        Int_ElelyExces_Ratio = temp_Int_ElelyExces_Ratio;
        DryOut = "On";
    print(f"Scan {Scan_i}: DryOut = {DryOut}")
    if DryOut == "On":  
        mdic_dry,Para_0 = Initialize_mdic_dry(Para_0,Int_ElelyExces_Ratio)
    else:
        mdic_dry ={}

    ##########################################################
    ##############    Part-2: Run model         ##############
    ##########################################################
    ##########################################################
    Timeout_text = 'I timed out'
    ##########    2-1: Define model and run break-in cycle
    try:  
        Timelimit = int(3600*2)
        # the following turns on for HPC only!
        if Timeout == True:
            timeout_RPT = TimeoutFunc(
                Run_Breakin, 
                timeout=Timelimit, 
                timeout_val=Timeout_text)
            Result_list_breakin  = timeout_RPT(
                model_options, Experiment_Breakin, 
                Para_0, mesh_list, submesh_strech)
        else:
            Result_list_breakin  = Run_Breakin(
                model_options, Experiment_Breakin, 
                Para_0, mesh_list, submesh_strech)
        [Model_0,Sol_0,Call_Breakin] = Result_list_breakin
        if Return_Sol == True:
            Sol_RPT.append(Sol_0)
        if Call_Breakin.success == False:
            print("Fail due to Experiment error or infeasible")
            1/0
        if Sol_0 == Timeout_text: # to do: distinguish different failure cases
            print("Fail due to Timeout")
            1/0
        if Sol_0 == "Model error or solver error":
            print("Fail due to Model error or solver error")
            1/0
    except ZeroDivisionError as e:
        str_error_Breakin = str(e)
        print(f"Scan {Scan_i}: Fail break-in cycle, need to exit the whole scan now due to {str_error_Breakin} but do not know how!")
        
        Flag_Breakin = False
    else:
        print(f"Scan {Scan_i}: Finish break-in cycle")
        # post-process for break-in cycle
        my_dict_RPT = GetSol_dict (my_dict_RPT,keys_all_RPT, Sol_0, 
            cycle_no, step_RPT_CD , step_RPT_CC , step_RPT_RE, step_AGE_CV   )
        cycle_count =0; 
        my_dict_RPT["Cycle_RPT"].append(cycle_count)
        Cyc_Update_Index.append(cycle_count);
        Flag_Breakin = True
        
    Flag_AGE = True; str_error_AGE_final = "Empty";   str_error_RPT = "Empty";
    #############################################################
    #######   2-2: Write a big loop to finish the long experiment    
    if Flag_Breakin == True: 
        k=0
        # Para_All.append(Para_0);Model_All.append(Model_0);Sol_All_i.append(Sol_0); 
        Para_0_Dry_old = Para_0;     Model_Dry_old = Model_0  ; Sol_Dry_old = Sol_0;   del Model_0,Sol_0
        while k < SaveTimes:    
            i=0    
            while i < Small_Loop:
                if DryOut == "On":
                    Data_Pack,Paraupdate   = Cal_new_con_Update (  Sol_Dry_old,   Para_0_Dry_old )
                if DryOut == "Off":
                    Paraupdate = Para_0
                # Run aging cycle:
                try:
                    Timelimit = int(3600*2)
                    if Timeout == True:
                        timeout_AGE = TimeoutFunc(
                            Run_Model_Base_On_Last_Solution, 
                            timeout=Timelimit, 
                            timeout_val=Timeout_text)
                        Result_list_AGE = timeout_AGE( 
                            Model_Dry_old  , Sol_Dry_old , Paraupdate ,Experiment_Long, 
                            Update_Cycles,Temper_i,mesh_list,submesh_strech )
                    else:
                        Result_list_AGE = Run_Model_Base_On_Last_Solution( 
                            Model_Dry_old  , Sol_Dry_old , Paraupdate ,Experiment_Long, 
                            Update_Cycles,Temper_i,mesh_list,submesh_strech )
                    [Model_Dry_i, Sol_Dry_i , Call_Age ] = Result_list_AGE
                    if Return_Sol == True:
                        Sol_AGE.append(Sol_Dry_i)
                    #print(f"Temperature for ageing is now: {Temper_i}")  
                    if Call_Age.success == False:
                        print("Fail due to Experiment error or infeasible")
                        str_error_AGE = "Experiment error or infeasible"
                        1/0
                    if Sol_Dry_i == Timeout_text: # fail due to timeout
                        print("Fail due to Timeout")
                        str_error_AGE = "Timeout"
                        1/0
                    if Sol_Dry_i == "Model error or solver error":
                        print("Fail due to Model error or solver error")
                        str_error_AGE = "Model error or solver error"
                        1/0
                except ZeroDivisionError as e:
                    print(f"Scan {Scan_i}: Fail during No.{Cyc_Update_Index[-1]} ageing cycles due to {str_error_AGE}")
                    Flag_AGE = False
                    str_error_AGE_final = str_error_AGE
                    break
                else:
                    Para_0_Dry_old = Paraupdate;       Model_Dry_old = Model_Dry_i;      Sol_Dry_old = Sol_Dry_i;   
                    del Paraupdate,Model_Dry_i,Sol_Dry_i
                    # post-process for first ageing cycle and every -1 ageing cycle
                    if k==0 and i==0:    
                        my_dict_AGE = GetSol_dict (my_dict_AGE,keys_all_AGE, Sol_Dry_old, 
                            0, step_AGE_CD , step_AGE_CC , step_RPT_RE, step_AGE_CV   )     
                        my_dict_AGE["Cycle_AGE"].append(1)
                    my_dict_AGE = GetSol_dict (my_dict_AGE,keys_all_AGE, Sol_Dry_old, 
                        cycle_no, step_AGE_CD , step_AGE_CC , step_RPT_RE, step_AGE_CV   )    
                    cycle_count +=  Update_Cycles; 
                    my_dict_AGE["Cycle_AGE"].append(cycle_count)           
                    Cyc_Update_Index.append(cycle_count)
                    print(f"Scan {Scan_i}: Finish for No.{Cyc_Update_Index[-1]} ageing cycles")
                    if DryOut == "On":
                        mdic_dry = Update_mdic_dry(Data_Pack,mdic_dry)
                    i += 1;   
            # run RPT, and also update parameters (otherwise will have problems)
            if DryOut == "On":
                Data_Pack , Paraupdate  = Cal_new_con_Update (  Sol_Dry_old,   Para_0_Dry_old   )
            if DryOut == "Off":
                Paraupdate = Para_0     
            try:
                Timelimit = int(3600*2)
                if Timeout == True:
                    timeout_RPT = TimeoutFunc(
                        Run_Model_Base_On_Last_Solution_RPT, 
                        timeout=Timelimit, 
                        timeout_val=Timeout_text)
                    Result_list_RPT = timeout_RPT(
                        Model_Dry_old  , Sol_Dry_old ,   
                        Paraupdate,      Experiment_RPT, RPT_Cycles, 
                        Temper_RPT ,mesh_list ,submesh_strech
                    )
                else:
                    Result_list_RPT = Run_Model_Base_On_Last_Solution_RPT(
                        Model_Dry_old  , Sol_Dry_old ,   
                        Paraupdate,      Experiment_RPT, RPT_Cycles, 
                        Temper_RPT ,mesh_list ,submesh_strech
                    )
                [Model_Dry_i, Sol_Dry_i,Call_RPT]  = Result_list_RPT
                if Return_Sol == True:
                    Sol_RPT.append(Sol_Dry_i)
                #print(f"Temperature for RPT is now: {Temper_RPT}")  
                if Call_RPT.success == False:
                    print("Fail due to Experiment error or infeasible")
                    str_error_RPT = "Experiment error or infeasible"
                    1/0 
                if Sol_Dry_i == Timeout_text:
                    print("Fail due to Timeout")
                    str_error_RPT = "Timeout"
                    1/0
                if Sol_Dry_i == "Model error or solver error":
                    print("Fail due to Model error or solver error")
                    str_error_RPT = "Model error or solver error"
                    1/0
            except ZeroDivisionError as e:
                print(f"Scan {Scan_i}: Fail during No.{Cyc_Update_Index[-1]} RPT cycles, due to {str_error_RPT}")
                break
            else:
                my_dict_RPT = GetSol_dict (my_dict_RPT,keys_all_RPT, Sol_Dry_i, 
                    cycle_no, step_RPT_CD , step_RPT_CC , step_RPT_RE, step_AGE_CV   )
                my_dict_RPT["Cycle_RPT"].append(cycle_count)
                Cyc_Update_Index.append(cycle_count)
                print(f"Scan {Scan_i}: Finish for No.{Cyc_Update_Index[-1]} RPT cycles")
                if DryOut == "On":
                    mdic_dry = Update_mdic_dry(Data_Pack,mdic_dry)
                Para_0_Dry_old = Paraupdate;    Model_Dry_old = Model_Dry_i  ;     Sol_Dry_old = Sol_Dry_i    ;   
                del Paraupdate,Model_Dry_i,Sol_Dry_i
                if Flag_AGE == False:
                    break
            k += 1 
    ############################################################# 
    #########   An extremely bad case: cannot even finish breakin
    if Flag_Breakin == False: 
        value_list_temp = list(Para_dict_i.values())
        values = []
        for value_list_temp_i in value_list_temp:
            values.append(str(value_list_temp_i))
        values.insert(0,str(Scan_i));
        values.insert(1,DryOut);
        values.extend([
            str_exp_AGE_text,
            str_exp_RPT_text,
            "nan","nan",
            "nan","nan", 
            "nan","nan",
            "nan","nan",
            "nan",str_error_Breakin])
        values = [values,]
        print(str_error_Breakin)
        print("Fail in {}".format(ModelTimer.time())) 
        book_name_xlsx_seperate =   str(Scan_i)+ '_' + book_name_xlsx;
        sheet_name_xlsx =  str(Scan_i);
        write_excel_xlsx(
            BasicPath + Target+book_name_xlsx_seperate, 
            sheet_name_xlsx, values)
        
        midc_merge = {**my_dict_RPT, **my_dict_AGE,**mdic_dry}

        return midc_merge,Sol_RPT,Sol_AGE
    ##########################################################
    ##############   Part-3: Post-prosessing    ##############
    ##########################################################
    # Newly add (220517): save plots, not just a single line in excel file:     
    # Newly add (221114): make plotting as functions
    # Niall_data = loadmat( 'Extracted_all_cell.mat'
    else:
        if not os.path.exists(BasicPath + Target + str(Scan_i)):
            os.mkdir(BasicPath + Target + str(Scan_i) );
        dpi= 100;
        # Update 230221 - Add model LLI, LAM manually 
        my_dict_RPT['Throughput capacity [kA.h]'] = (
            np.array(my_dict_RPT['Throughput capacity [A.h]'])/1e3).tolist()
        my_dict_RPT['CDend SOH [%]'] = ((
            np.array(my_dict_RPT["Discharge capacity [A.h]"])
            /my_dict_RPT["Discharge capacity [A.h]"][0])*100).tolist()
        my_dict_RPT["CDend LAM_ne [%]"] = ((1-
            np.array(my_dict_RPT['CDend Negative electrode capacity [A.h]'])
            /my_dict_RPT['CDend Negative electrode capacity [A.h]'][0])*100).tolist()
        my_dict_RPT["CDend LAM_pe [%]"] = ((1-
            np.array(my_dict_RPT['CDend Positive electrode capacity [A.h]'])
            /my_dict_RPT['CDend Positive electrode capacity [A.h]'][0])*100).tolist()
        my_dict_RPT["CDend LLI [%]"] = ((1-
            np.array(my_dict_RPT["CDend Total lithium capacity in particles [A.h]"])
            /my_dict_RPT["CDend Total lithium capacity in particles [A.h]"][0])*100).tolist()
        if model_options.__contains__("SEI"):
            my_dict_RPT["CDend LLI SEI [%]"] = ((
                np.array(
                    my_dict_RPT["CDend Loss of capacity to SEI [A.h]"]-
                    my_dict_RPT["CDend Loss of capacity to SEI [A.h]"][0]
                    )
                /my_dict_RPT["CDend Total lithium capacity in particles [A.h]"][0])*100).tolist()
        if model_options.__contains__("SEI on cracks"):
            my_dict_RPT["CDend LLI SEI on cracks [%]"] = ((
                np.array(
                    my_dict_RPT["CDend Loss of capacity to SEI on cracks [A.h]"]-
                    my_dict_RPT["CDend Loss of capacity to SEI on cracks [A.h]"][0]
                    )
                /my_dict_RPT["CDend Total lithium capacity in particles [A.h]"][0])*100).tolist()
        if model_options.__contains__("lithium plating"):
            my_dict_RPT["CDend LLI lithium plating [%]"] = ((
                np.array(
                    my_dict_RPT["CDend Loss of capacity to lithium plating [A.h]"]-
                    my_dict_RPT["CDend Loss of capacity to lithium plating [A.h]"][0]
                    )
                /my_dict_RPT["CDend Total lithium capacity in particles [A.h]"][0])*100).tolist()
  
        ##########################################################
        #########      3-1: Plot cycle,location, Dryout related 
        Plot_Cyc_RPT_4(
            my_dict_RPT,
            Exp_Any_AllData,Temp_Cell_Exp, Plot_Exp ,  # =True or False,
            Scan_i,Temper_i,model_options,BasicPath, Target,fs,dpi)
        if len(my_dict_AGE["CDend Porosity"])>1:
            Plot_Loc_AGE_4(my_dict_AGE,Scan_i,model_options,BasicPath, Target,fs,dpi)
        if DryOut == "On":
            Plot_Dryout(Cyc_Update_Index,mdic_dry,ce_EC_0,Scan_i,BasicPath, Target,fs,dpi)
        ##########################################################
        #########      3-2: Save data as .mat 
        my_dict_RPT["Cyc_Update_Index"] = Cyc_Update_Index
        my_dict_RPT["SaveTimes"]    = SaveTimes
        midc_merge = {**my_dict_RPT, **my_dict_AGE,**mdic_dry}
        savemat(BasicPath + Target+ str(Scan_i) + '/' + str(Scan_i)+ '-StructDara_for_Mat.mat',midc_merge)  
        ##########################################################
        #########      3-3: Save summary to excel 
        values=Get_Values_Excel(
            model_options,my_dict_RPT,mdic_dry,
            DryOut,Scan_i,Para_dict_i,str_exp_AGE_text,
            str_exp_RPT_text,
            str_error_AGE_final,
            str_error_RPT)
        values = [values,]
        book_name_xlsx_seperate =   str(Scan_i)+ '_' + book_name_xlsx;
        sheet_name_xlsx =  str(Scan_i);
        write_excel_xlsx(
            BasicPath + Target+book_name_xlsx_seperate, 
            sheet_name_xlsx, values)
        print("Succeed doing something in {}".format(ModelTimer.time()))
        print('This is the end of No.', Scan_i, ' scan')
        return midc_merge,Sol_RPT,Sol_AGE



print("Reach end of the file")