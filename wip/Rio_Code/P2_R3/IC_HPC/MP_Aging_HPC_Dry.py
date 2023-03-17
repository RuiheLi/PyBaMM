import pybamm as pb;import pandas as pd   ;import numpy as np;import os;import matplotlib.pyplot as plt;import os;#import imageio
from scipy.io import savemat,loadmat;from pybamm import constants,exp;import matplotlib as mpl; fs=17; # or we can set import matplotlib.pyplot as plt then say 'mpl.rc...'
for k in range(0,1):
    mpl.rcParams["axes.labelsize"] = fs
    mpl.rcParams["axes.titlesize"] = fs
    mpl.rcParams["xtick.labelsize"] =  fs
    mpl.rcParams["ytick.labelsize"] =  fs
    mpl.rcParams["legend.fontsize"] =  fs
    mpl.rcParams['font.sans-serif'] = ['Times New Roman']
    mpl.rcParams['font.serif'] = ['Times New Roman']
    mpl.rcParams['axes.unicode_minus'] = False

import openpyxl
import traceback

# Add funtion to sys.path
os.chdir(pb.__path__[0]+'/..')   
# pb.__path__[0] allow us to location path without concerning different machine
import sys  
sys.path.append(os.path.join(os.getcwd(),'wip\Rio_Code\Fun_P2'))  
from Fun_P2_Union import (GetScan,Run_model_w_dry_out ,Run_model_wo_Dry_out,Cal_new_con_Update,Run_Model_Base_On_Last_Solution,Run_Model_Base_On_Last_Solution_RPT,write_excel_xlsx,)
import multiprocessing

# Run model with dry-out:
for i in range(0,1):  
    Total_Cycles = 104; Cycle_bt_RPT = 13; Update_Cycles = 26;  
    CyclePack = [Total_Cycles,Cycle_bt_RPT,Update_Cycles];
    # Key scan parameters:
    Ratio_excess = [1.02,1.04,1.10,];
    cs_Neg_Init = [29866,]; Diff_SEI = [1.7e-20,];    # default: 29866
    R_SEI = [2E5,];   Bulk_Sol_Con =[ 4541.0,];
    D_Li_inSEI = [3e-19,3e-18,3e-17,  3e-16,3e-15,3e-14,3e-13,];    # default: 1e-20 
    c_Li_inte_ref = [15,];    # default: 15
    Couple_SEI_LiP = [1e-6,]; # default: 1e-6
    k_LiP = [1E-10];         # default: 1e-10
    Temper = [10,25,45];
    MESH_PAR = [30,];

    (TotalScan, DatePack_scan) = GetScan(Ratio_excess,cs_Neg_Init,Diff_SEI,
        R_SEI,Bulk_Sol_Con,D_Li_inSEI,c_Li_inte_ref,Couple_SEI_LiP,k_LiP,Temper,MESH_PAR);
    BasicPath = 'D:/OneDrive - Imperial College London/SimDataSave/P2R3/'; 
    # BasicPath=os.getcwd() # for HPC
    Target  = '/HPC_Dry_220725/'
    if not os.path.exists(BasicPath + Target):
        os.mkdir(BasicPath + Target);
    book_name_xlsx = 'Inte_SEI_104cycles.xlsx';

    sheet_name_xlsx = 'Results'
    value3 = [
        ["Index", "Ratio_ex","cs_Neg_Init", "Diff_SEI", "R_SEI", 
        "Bulk_Sol_Con","D_Li_inSEI", "c_Li_inte_ref",
        "Couple_SEI_LiP","k_LiP","Temper_i","mesh_par",
        "Cap Loss","LLI to LiP",
        "LLI to SEI","Vol_Elely_Tot Final", "Vol_Elely_JR Final","Width Final","Error"],
        ]
    write_excel_xlsx(BasicPath + Target+book_name_xlsx, sheet_name_xlsx, value3)



if __name__ == "__main__":
    pool = multiprocessing.Pool(63)
    processes = [pool.apply_async(
        Run_model_w_dry_out, 
        args=(
            CyclePack , DatePack_scan_jj, BasicPath , Target, book_name_xlsx,
        )
            ) 
            for DatePack_scan_jj in DatePack_scan]
    result = [p.get() for p in processes]

