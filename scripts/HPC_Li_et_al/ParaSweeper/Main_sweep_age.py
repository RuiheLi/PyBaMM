### Example notebook for a general script to run aging cycle

import csv, random, os,sys
import pybamm as pb;import pandas as pd;import numpy as np;
import os, json,openpyxl,traceback,multiprocessing,scipy.optimize
import matplotlib.pyplot as plt;
import imageio,timeit,random,time, signal
from scipy.io import savemat,loadmat;
from pybamm import constants,exp;import matplotlib as mpl; 
fs=17; 
font = {'family' : 'DejaVu Sans','size'   : fs}
mpl.rc('font', **font)

########################     Global settings!!!
rows_per_file = 4;  Scan_end_end = 1000;
purpose_i = "Example_Sweep_Age"

On_HPC =   True
Runshort = False                    # a long run or a quick test

if On_HPC:
    i_bundle = int(os.environ["PBS_ARRAY_INDEX"])
else:
    i_bundle = 1; # manually specify

Scan_start = (i_bundle-1)*rows_per_file+1;    
Scan_end   = min(Scan_start + rows_per_file-1, Scan_end_end)    
purpose = f"{purpose_i}_Case_{Scan_start}_{Scan_end}"
# interpetation: Simnon suggested, with cracking activation, heat transfer
para_csv = f"Bundle_{i_bundle}.csv"  # name of the random file to get parameters

# specify path:
if On_HPC:                          # Run on HPC
    Path_Data_pre = "InputData/" 
    BasicPath=os.getcwd() 
    Para_file = f"InputData/{purpose_i}/"  +  para_csv
else:
    Path_Data_pre = os.path.expanduser(
        "~/EnvPB_Linux/PyBaMM/scripts/HPC_Li_et_al/HPC_Paper_SimSave/Example_Expdata") # for Linux
    BasicPath =  os.path.expanduser(
        "~/EnvPB_Linux/PyBaMM/scripts/HPC_Li_et_al/HPC_Paper_SimSave/Results")
    Para_file = os.path.expanduser(
        "~/EnvPB_Linux/PyBaMM/scripts/HPC_Li_et_al/HPC_Paper_SimSave")+f'/Get_Input/{purpose_i}/'+para_csv

if not os.path.exists(BasicPath +"/"+ purpose):
    os.mkdir(BasicPath +"/"+ purpose);

from scripts.HPC_Li_et_al.ParaSweeper.Fun_HPC import *
parameter_names, combinations = load_combinations_from_csv(Para_file)


pool_no = len(combinations)
Indexs  = np.arange(Scan_start-1,Scan_end)
index_list = Indexs+1

# Get all para
Para_dict_list = []
# get all dictionaries
for combination in combinations:
    input_dict = {}
    for parameter_name,para_value in zip(parameter_names,combination ):
        input_dict[parameter_name] = para_value
    Para_dict_list.append(input_dict)
print(f"Total scan case is {len(Para_dict_list)}")

# try to plot experiment data
Path_Data  = Path_Data_pre+"/example_Age_SOH.csv"
my_data = pd.read_csv( Path_Data)

cap_thr   = my_data["Charge Throughput (A.h)"] / 1e3
cap_C_10  = my_data["C/10 Capacity (mA.h)"]/1e3
soh_C_10  = cap_C_10 / cap_C_10[0] * 100
# plt.plot(cap_thr,soh_C_10,"-o")
cap_thr   = np.array(cap_thr).tolist()
cap_C_10  = np.array(cap_C_10).tolist()
soh_C_10  = np.array(soh_C_10).tolist()
XY_Exp =[cap_thr,soh_C_10]

#####################################################################
########################  Define experiment  ########################
#####################################################################
V_max = 4.2;        V_min = 2.5; 
exp_AGE_text = [(
    f"Discharge at 1C until {V_min}V", 
    f"Charge at 0.3C until {V_max}V",
    f"Hold at {V_max} V until C/100",
    ),  ]  # *  78
step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2;
exp_RPT_text = [ (
    # refill
    "Rest for 1 minutes", # change this step to rest to avoid errors. 
    "Rest for 1 hours", 
    # 0.1C cycle 
    f"Discharge at 0.1C until {V_min} V",  
    "Rest for 3 hours (20 minute period)",  
    f"Charge at 0.1C until {V_max} V",
    f"Hold at {V_max}V until C/100",
    "Rest for 1 hours",
    ) ] * 1
exp_breakin_text = [ (
    # refill
    f"Hold at {V_max}V until C/100",
    "Rest for 1 hours", 
    # 0.1C cycle 
    f"Discharge at 0.1C until {V_min} V",  
    "Rest for 3 hours (20 minute period)",  
    f"Charge at 0.1C until {V_max} V",
    f"Hold at {V_max}V until C/100",
    "Rest for 1 hours",
    ) ] * 1
# step index for RPT
step_0p1C_CD = 2; step_0p1C_CC = 4;   step_0p1C_RE =3;    
cycle_no = -1; 

#####################################################################
########################  Output  ###################################
#####################################################################
#####################################################################
########################  Output  ###################################
#####################################################################
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
    #"CCend Negative electrode roughness ratio",
    #"CCend Total SEI on cracks thickness [m]",

    "CDend Porosity",
    "CDend Negative electrode interfacial current density [A.m-2]",
    "CDend Electrolyte potential [V]",
    "CDend Electrolyte concentration [mol.m-3]",
    "CDend Negative electrode reaction overpotential [V]",
    "CDend Negative particle surface concentration [mol.m-3]",
    #"CDend Negative electrode roughness ratio",
    #"CDend Total SEI on cracks thickness [m]",
    #"REend Total SEI on cracks thickness [m]",
]
keys_tim_RPT = [
    # default: CD
    "CD Time [h]",
    "CD Terminal voltage [V]",
    #"RE Terminal voltage [V]",
]
keys_cyc_RPT = [   # default: CDend
    "Discharge capacity [A.h]",
    "Throughput capacity [A.h]",
    "CDend Total lithium capacity in particles [A.h]",
    "CDend Loss of capacity to lithium plating [A.h]",
    "CDend Loss of capacity to SEI [A.h]",
    "CDend Loss of capacity to SEI on cracks [A.h]",
    #"CDend X-averaged total SEI on cracks thickness [m]",
    #"CDend X-averaged negative electrode roughness ratio",
    "CDend Local ECM resistance [Ohm]",
    "CDsta Negative electrode stoichiometry", 
    "CDend Negative electrode stoichiometry",
    "CDsta Positive electrode stoichiometry", 
    "CDend Positive electrode stoichiometry",
    "CDend Negative electrode capacity [A.h]",
    "CDend Positive electrode capacity [A.h]",
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
    #"CCend Negative electrode roughness ratio",
    #"CCend Total SEI on cracks thickness [m]",

    "CDend Porosity",
    "CDend Negative electrode interfacial current density [A.m-2]",
    "CDend Electrolyte potential [V]",
    "CDend Electrolyte concentration [mol.m-3]",
    "CDend Negative electrode reaction overpotential [V]",
    "CDend Negative particle surface concentration [mol.m-3]",
    #"CDend Negative electrode roughness ratio",
    #"CDend Total SEI on cracks thickness [m]",
    "CDend Electrolyte diffusivity [m2.s-1]",
    "CDend Electrolyte conductivity [S.m-1]",
]
keys_tim_AGE = [];
keys_cyc_AGE = [];
keys_all_RPT = [keys_loc_RPT,keys_tim_RPT,keys_cyc_RPT];
keys_all_AGE = [keys_loc_AGE,keys_tim_AGE,keys_cyc_AGE];
keys_all = [keys_all_RPT,keys_all_AGE];

# Write para - 1st round:
Values_1 = []
head_keys = list(Para_dict_list[0].keys())
head_pre = [
    "Scan No","Y or N",
    "Error %","Punish",
    "Dry out",]

head_pos = [ 
   "exp_AGE_text", "exp_RPT_text",
   "SOH [%]","LLI [%]","LLI to LiP [%]",
   "LLI to SEI [%]","LLI to sei-on-cracks [%]",
   "LLI due to LAM [%]","LAM_to_Crack_NE [%]","LAM_to_Crack_PE [%]",
   "LAM_to_Dry [%]","Error"]
Values_1 .append([*head_pre,*head_keys,*head_pos])
book_name_xlsx = f'Summary_{purpose}.xlsx';
sheet_name_xlsx = 'Output'
Target  = f'/{purpose}/'
if not os.path.exists(BasicPath +Target):
    os.mkdir(BasicPath +Target);
write_excel_xlsx(
    BasicPath + Target + book_name_xlsx, 
    sheet_name_xlsx, Values_1)   
Exp_pack =[
    exp_AGE_text,exp_RPT_text,exp_breakin_text,
    step_AGE_CD,step_AGE_CC,step_AGE_CV,
    step_0p1C_CD ,step_0p1C_CC,step_0p1C_RE,
    cycle_no,book_name_xlsx,
    ] 
if not os.path.exists(BasicPath +Target+"Mats"):
    os.mkdir(BasicPath +Target +"Mats");
if not os.path.exists(BasicPath +Target+"Plots"):
    os.mkdir(BasicPath +Target+"Plots");
if not os.path.exists(BasicPath +Target+"Excel"):
    os.mkdir(BasicPath +Target+"Excel");


midc_merge_all = [];Sol_RPT_all = [];Sol_AGE_all = [];
Plot_Exp=True;     Timeout=False;  flag_RunOneCyc = False
Return_Sol=True;   Check_Small_Time=True;
fs = 13; dpi = 100;

""" my_dict_1Cyc,Sol_RPT,Sol_0,Call_1Cyc = Run_model (
    Para_dict_list[2],BasicPath, XY_Exp, 
    purpose,    Exp_pack, keys_all,dpi,fs,
    flag_RunOneCyc,   Plot_Exp,Timeout,Return_Sol,
    Check_Small_Time
) """

if __name__ == "__main__":
    pool = multiprocessing.Pool(int(pool_no))
    processes = [
    pool.apply_async(
        Run_model, 
        args=(
            Para_dict_i,BasicPath, XY_Exp, 
            purpose,    Exp_pack, keys_all,dpi,fs,
            Runshort,   flag_RunOneCyc,   Plot_Exp,
            Timeout,Return_Sol, Check_Small_Time
        ) )
        for Para_dict_i in Para_dict_list]
    Result = [p.get() for p in processes]  

for result in Result:
    midc_merge_all.append(result[0])
    Sol_RPT_all.append(result[1])
    Sol_AGE_all.append(result[2]) 

# Write summary into excel
Index_List_succeed = index_list
for k,index_i in enumerate(Index_List_succeed):
    #print(index_i)
    try:
        old_book = str(index_i) + '_' + book_name_xlsx
        #print(old_book)
        #open excel:
        data_old = openpyxl.load_workbook(
            BasicPath + Target+ "Excel/" + old_book)   
        data_tar = openpyxl.load_workbook(
            BasicPath + Target + book_name_xlsx) 

        table_old = data_old[str(index_i)]
        nrows_old = table_old.max_row  # 获得行数
        ncolumns_old = table_old.max_column  # 获得列数

        table_tar = data_tar[sheet_name_xlsx]
        nrows_tar = table_tar.max_row # ncolumns_old + k +1 # Mark!!! Most important changes!
        ncolumns_old = table_old.max_column  # 获得列数
        list_old = [];
        #print(nrows_old,nrows_tar)
        for i in range(1,nrows_old+1):
            for j in range(1,ncolumns_old+1):
                list_old.append(table_old.cell(row=i,column=j).value)
        
        list_old = [list_old,]
        for i in range(1, len(list_old)+1):
                for j in range(1, len(list_old[i-1])+1):
                    #print(i,j,list_old[i-1][j-1]    )
                    table_tar.cell(nrows_tar+i, j).value = list_old[i-1][j-1]     
        data_tar.save(BasicPath + Target + book_name_xlsx) 
        data_tar.close()
    except:
        print(f"Something goes wrong for Scan {index_i}!")
    else:
        print(f"Successfuly write results for Scan {index_i}!") 

















