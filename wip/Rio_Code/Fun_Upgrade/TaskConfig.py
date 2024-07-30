# define class for all inputs and customized settings:
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
#from .Post_process import *
#from .Fun_Upgrade import *


# Function to read exp
def Read_Exp(BasicPath_Save,Exp_Any_Cell,Exp_Path,Exp_head,Exp_Any_Temp,i):
    Exp_Any_AllData  = {}
    for cell in Exp_Any_Cell:
        Exp_Any_AllData[cell] = {} # one cell
        # For extracted directly measured capacity, resistance, etc.
        Exp_Any_AllData[cell]["Extract Data"] = pd.read_csv(
            BasicPath_Save+Exp_Path[i]+
            f"{Exp_head[i]} - cell {cell} ({Exp_Any_Temp[cell]}degC) - Extracted Data.csv", 
            index_col=0)
        Exp_Any_AllData[cell]["Extract Data"].loc[     # update 230617- fill this with avg T
            0, 'Age set average temperature (degC)'] = float(Exp_Any_Temp[cell])
        # Read for DMA results, further a dictionary
        Exp_Any_AllData[cell]["DMA"] = {}
        Exp_Any_AllData[cell]["DMA"]["Cap_Offset"]=pd.read_csv(
            BasicPath_Save+Exp_Path[i]+ "DMA Output/" + f"cell {cell}/" + 
            f"{Exp_head[i]} - Capacity and offset data from OCV-fitting for cell {cell}.csv", 
            index_col=0)
        Exp_Any_AllData[cell]["DMA"]["LLI_LAM"]=pd.read_csv(
            BasicPath_Save+Exp_Path[i]+ "DMA Output/" + f"cell {cell}/" + 
            f"{Exp_head[i]} - DM data from OCV-fitting for cell {cell}.csv", 
            index_col=0)
        Exp_Any_AllData[cell]["DMA"]["Fit_SOC"]=pd.read_csv(
            BasicPath_Save+Exp_Path[i]+ "DMA Output/" + f"cell {cell}/" + 
            f"{Exp_head[i]} - fitting parameters from OCV-fitting for cell {cell}.csv", 
            index_col=0)
        Exp_Any_AllData[cell]["DMA"]["RMSE"]=pd.read_csv(
            BasicPath_Save+Exp_Path[i]+ "DMA Output/" + f"cell {cell}/" + 
            f"{Exp_head[i]} - RMSE data from OCV-fitting for cell {cell}.csv", 
            index_col=0)
        # update 230726: read 0.1C voltage curve, discharge only
        Exp_Any_AllData[cell]["0.1C voltage"] = {}
        for m in range(16):
            try:
                C_10_curve_temp = pd.read_csv(
                    BasicPath_Save+Exp_Path[i]+ "0.1C Voltage Curves/"+ f"cell {cell}/" +
                    f"{Exp_head[i]} - cell {cell} - RPT{m} - 0.1C discharge data.csv", )    
            except:
                pass # print(f"Exp-{i+1} - Cell {cell} doesn't have RPT {m}")
            else:
                C_10_curve_temp["Time (h)"] = (
                    C_10_curve_temp["Time (s)"] - 
                    C_10_curve_temp["Time (s)"].iloc[0]) / 3600
                Exp_Any_AllData[cell]["0.1C voltage"][f"RPT{m}"] = C_10_curve_temp
                # print(f"Read Exp-{i+1} - Cell {cell} RPT {m}")
    print("Finish reading Experiment!")
    return Exp_Any_AllData


class Exp_config:
    def __init__(self, V_max = 4.2 , V_min = 2.5, 
            exp_AGE_text=None, step_AGE_CD=None, step_AGE_CC=None, step_AGE_CV=None,
            exp_breakin_text=None, exp_RPT_text=None, exp_GITT_text=None, 
            exp_refill=None, exp_adjust_before_age=None,
            step_0p1C_CD=None, step_0p1C_CC=None, step_0p1C_RE=None, step_0p5C_CD=None,
            Exp_No=None,Age_T=None,tot_cyc=None,Cycle_bt_RPT=None,update=None,
            Runs_bt_RPT=None,RPT_num=None,
            RPT_Cycles=None,Age_T_in_K=None,RPT_T_in_K=None):

        self.V_max = V_max 
        self.V_min = V_min

        self.exp_AGE_text = exp_AGE_text
        self.step_AGE_CD = step_AGE_CD
        self.step_AGE_CC = step_AGE_CC
        self.step_AGE_CV = step_AGE_CV
        self.exp_breakin_text = exp_breakin_text
        self.exp_RPT_text = exp_RPT_text
        self.exp_GITT_text = exp_GITT_text
        self.exp_refill = exp_refill
        self.exp_adjust_before_age = exp_adjust_before_age
        self.step_0p1C_CD = step_0p1C_CD
        self.step_0p1C_CC = step_0p1C_CC
        self.step_0p1C_RE = step_0p1C_RE
        self.step_0p5C_CD = step_0p5C_CD

        self.Exp_No = Exp_No
        self.Age_T = Age_T
        self.tot_cyc = tot_cyc
        self.Cycle_bt_RPT = Cycle_bt_RPT
        self.update = update
        self.Runs_bt_RPT = Runs_bt_RPT
        self.RPT_num = RPT_num
        self.RPT_Cycles = RPT_Cycles
        self.Age_T_in_K=Age_T_in_K
        self.RPT_T_in_K=RPT_T_in_K

class Para_config:
    def __init__(self, Para_dict_i=None):
        self.Para_dict_i = Para_dict_i

class Model_config:
    def __init__(self,mesh_list=None,
            submesh_strech=None,model_options=None,DryOut=None ):
        self.mesh_list=mesh_list
        self.submesh_strech=submesh_strech
        self.model_options=model_options
        self.DryOut=DryOut

class ExpData_config:
    def __init__(self,cap_0 = None, 
                Temp_Cell_Exp = None, 
                Exp_Any_AllData = None,
                XY_pack = None):
        self.cap_0 = cap_0 # set initial capacity to standardize SOH and get initial SEI thickness 
        self.Temp_Cell_Exp = Temp_Cell_Exp
        self.Exp_Any_AllData = Exp_Any_AllData
        self.XY_pack = XY_pack
        

class Global_config:
    def __init__(self, On_HPC=False, Runshort="GEM-2", 
            Plot_Exp=True, Timeout=True,Return_Sol=True, 
            Check_Short_Time=True, R_from_GITT=True,
            fs=13, dpi=100, Re_No=0, Timelimit=int(3600*48) ,
            Timeout_text = 'I timed out',
            colormap = "cool"):
        self.On_HPC = On_HPC
        self.Runshort = Runshort
        self.Plot_Exp = Plot_Exp
        self.Timeout = Timeout
        self.Return_Sol = Return_Sol
        self.Check_Short_Time = Check_Short_Time
        if Runshort == "GEM-2":
            self.R_from_GITT = R_from_GITT
        else:   # overwrite 
            self.R_from_GITT = False 
        
        self.fs = fs
        self.dpi = dpi
        self.Re_No = Re_No
        self.Timelimit = Timelimit
        self.Timeout_text = Timeout_text
        self.colormap = colormap
       

class Path_config:
    def __init__(
            self, On_HPC = None, Path_Input = None, BasicPath_Save = None,
            purpose_i = None, case_no = None, 
            rows_per_file = None, Scan_end_end = None ):
        
        Scan_start = (case_no-1)*rows_per_file+1  
        Scan_end   = min(Scan_start + rows_per_file-1, Scan_end_end)    
        purpose = f"{purpose_i}_Case_{Scan_start}_{Scan_end}"
        Target  = f'/{purpose}/'
        para_csv = f"Bundle_{case_no}.csv"  
        # Path setting:
        if On_HPC:                          # Run on HPC
            Path_csv = f"InputData/{purpose_i}/" 
            Para_file = Path_csv +  para_csv
        else:
            Para_file = Path_Input+f'{purpose_i}/'+para_csv
        self.BasicPath_Save = BasicPath_Save
        self.Path_Input = Path_Input
        self.Target = Target
        self.purpose = purpose
        self.purpose_i = purpose_i
        self.Para_file = Para_file
        self.rows_per_file = rows_per_file
        self.Scan_end_end = Scan_end_end
        self.case_no = case_no

class TaskConfig:
    def __init__(
            self, path_config, global_config, exp_config,
            para_config,model_config,expData_config,Scan_No=None,): 

        self.path_config = path_config  
        self.global_config = global_config
        self.exp_config = exp_config
        self.para_config = para_config
        self.model_config = model_config 
        self.expData_config = expData_config
        self.Scan_No = Scan_No
        

    # do further initialization:
    def Get_Para_dict_i(self,Para_dict_i):
        self.para_config.Para_dict_i = Para_dict_i

    def Get_Exp_config(self):
        self.para_config.Para_dict_i["Scan No"] = int(self.para_config.Para_dict_i["Scan No"])
        self.para_config.Para_dict_i["Exp No."] = int(self.para_config.Para_dict_i["Exp No."])
        Scan_No = self.para_config.Para_dict_i["Scan No"]

        Exp_No = self.para_config.Para_dict_i["Exp No."]
        Age_T = self.para_config.Para_dict_i["Ageing temperature"]  

        Runshort = self.global_config.Runshort
        self.exp_config.Exp_No = Exp_No 
        self.exp_config.Age_T = Age_T
        self.Scan_No = Scan_No
        print(f"Aging T is {Age_T}degC")
        if Runshort == "GEM-2":
            if Exp_No == 1: # 0~30% SOC cycling
                tot_cyc = 3084; Cycle_bt_RPT = 257; update = 257; 
            if Exp_No == 2: # 70~85% SOC cycling                 
                tot_cyc = 6168; Cycle_bt_RPT = 514; update = 514;  # actually 516 on cycler
            if Exp_No == 3: # 70~85% SOC cycling
                tot_cyc = 6168; Cycle_bt_RPT = 514; update = 514;  # actually 516 on cycler
            if Exp_No == 5: # 0~100% SOC cycling
                tot_cyc = 1155; Cycle_bt_RPT = 77; update = 77;    # actually 77 on cycler
            if Exp_No == 6:
                tot_cyc = 1170*10; Cycle_bt_RPT = 77; update = 77;
            if Exp_No == 7:
                tot_cyc = 1170*10; Cycle_bt_RPT = 77; update = 77;
            if Exp_No == 8:
                tot_cyc = 1170*10; Cycle_bt_RPT = 77; update = 77;
            if Exp_No == 9:
                tot_cyc = 1170*10; Cycle_bt_RPT = 77; update = 77;
        elif Runshort == "Reservoir":
            pass
        else:  # really running short    
            if Exp_No in list(np.arange(1,10)):
                tot_cyc = 2; Cycle_bt_RPT = 1; update =1           # 2 1 1 
            else:
                print(f"Exp {Exp_No} Not yet implemented!")
        self.exp_config.tot_cyc = int(tot_cyc)
        self.exp_config.Cycle_bt_RPT = int(Cycle_bt_RPT)
        self.exp_config.update = int(update)
        self.exp_config.Runs_bt_RPT =  int(Cycle_bt_RPT/update) 
        self.exp_config.RPT_num = int(tot_cyc/Cycle_bt_RPT)



    # # Get information for Experiment - LG M50 Ageing dataset
    def Get_Exp_Pack(self):
        Exp_No = self.exp_config.Exp_No
        Age_T = self.exp_config.Age_T
        Path_Input = self.path_config.Path_Input
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
        Temp_Cell_Exp_All = [
            Temp_Cell_Exp_1,Temp_Cell_Exp_2,Temp_Cell_Exp_3,
            Temp_Cell_Exp_4,Temp_Cell_Exp_5]
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
            "A":[2/255, 3/255, 226/255,0.7],"B":[2/255, 3/255, 226/255,0.7],
            "J":[2/255, 3/255, 226/255,0.7],
            "D":[0, 0, 0,0.7],"E":[0, 0, 0,0.7],"F":[0, 0, 0,0.7],
            "K":[1,0,0,0.4],"L":[1,0,0,0.4],"M":[1,0,0,0.4],},
            {
            "A":[2/255, 3/255, 226/255,0.7],"B":[2/255, 3/255, 226/255,0.7],
            "D":[0, 0, 0,0.7],"C":[0, 0, 0,0.7],
            "E":[1,0,0,0.4],"F":[1,0,0,0.4],},
            {
            "A":[2/255, 3/255, 226/255,0.7],"B":[2/255, 3/255, 226/255,0.7],
            "C":[2/255, 3/255, 226/255,0.7],
            "D":[0, 0, 0,0.7],"E":[0, 0, 0,0.7],"F":[0, 0, 0,0.7],
            "G":[1,0,0,0.4],"H":[1,0,0,0.4],"I":[1,0,0,0.4],},
            {
            "A":[2/255, 3/255, 226/255,0.7],"B":[2/255, 3/255, 226/255,0.7],
            "C":[2/255, 3/255, 226/255,0.7],
            "D":[0, 0, 0,0.7],"E":[0, 0, 0,0.7],
            "F":[1,0,0,0.4],"G":[1,0,0,0.4],"H":[1,0,0,0.4],},
            {
            "A":[2/255, 3/255, 226/255,0.7],"B":[2/255, 3/255, 226/255,0.7],
            "C":[2/255, 3/255, 226/255,0.7],
            "D":[0, 0, 0,0.7],"E":[0, 0, 0,0.7],
            "F":[1,0,0,0.4],"G":[1,0,0,0.4],"H":[1,0,0,0.4],}]
        # Exp_No should be 1~5 to really have experimental data
        if Exp_No in list(np.arange(1,6)) and \
            int(Age_T) in [10,25,40]:

            Temp_Cell_Exp = Temp_Cell_Exp_All[Exp_No-1] 
            Exp_Any_AllData = Read_Exp(
                Path_Input,Exp_All_Cell[Exp_No-1],
                Exp_Path,Exp_head,Exp_Temp_Cell[Exp_No-1],
                Exp_No-1)
        else:
            Temp_Cell_Exp = "nan"
            Exp_Any_AllData = "nan"
        self.expData_config.Temp_Cell_Exp = Temp_Cell_Exp 
        self.expData_config.Exp_Any_AllData = Exp_Any_AllData
    # initialize expeirment text:
    def Initialize_exp_text(self):
        Exp_No = self.exp_config.Exp_No
        V_max = self.exp_config.V_max
        V_min = self.exp_config.V_min
        if Exp_No ==1:
            charge_time_mins = 1 * 60 
            exp_AGE_text = [(
                f"Discharge at 1C until {V_min}V", 
                f"Charge at 0.3C for {charge_time_mins} minutes or until {V_max}V",
                ),  ]  # *  setting on cycler is 516, rather than 514 in wiki
            step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2;
        
        elif Exp_No ==2:
            discharge_time_mins = 0.15* 60 * 4.86491/5
            charge_time_mins = 0.5* 60 * 4.86491/5
            exp_AGE_text = [(
                f"Discharge at 1C for {discharge_time_mins} minutes or until {V_min}V", 
                f"Charge at 0.3C for {charge_time_mins} minutes or until {V_max}V",
                ),  ]  # *  setting on cycler is 516, rather than 514 in wiki
            step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2;
        elif Exp_No ==3:
            discharge_time_mins = 0.15* 60 * 4.86491/5
            charge_time_mins = 0.5* 60 * 4.86491/5
            exp_AGE_text = [(
                f"Discharge at 1C for {discharge_time_mins} minutes or until {V_min}V", 
                f"Charge at 0.3C until {V_max}V",
                f"Hold at {V_max} V until C/100",
                ),  ]   # *  setting on cycler is 515, rather than 514 in wiki
            step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2;
        elif Exp_No ==5:
            exp_AGE_text = [(
                f"Discharge at 1C until {V_min}V", 
                f"Charge at 0.3C until {V_max}V",
                f"Hold at {V_max} V until C/100",
                ),  ]  # *  78
            step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2
        elif Exp_No ==6:
            exp_AGE_text = [(
                f"Discharge at 1C until {V_min}V", 
                f"Charge at 1.2C until {V_max}V",
                f"Hold at {V_max} V until C/100",
                ),  ] 
            step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2
        elif Exp_No ==7:
            exp_AGE_text = [(
                f"Discharge at 0.5C until {V_min}V",  
                f"Charge at 0.3C until {V_max}V",
                f"Hold at {V_max} V until C/100",
                ),  ] 
            step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2
        elif Exp_No ==8:
            exp_AGE_text = [(
                f"Discharge at 2C until {V_min}V", 
                f"Charge at 0.3C until {V_max}V",
                f"Hold at {V_max} V until C/100",
                ),  ]
            step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2
        elif Exp_No ==9:
            exp_AGE_text = [(
                f"Discharge at 1C until {V_min}V", 
                f"Charge at 0.3C until {V_max}V",
                f"Hold at {V_max} V until C/100",
                ),  ] 
            step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2
        else:
            print("Not yet implemented!")

        # now for RPT: 
        exp_RPT_Need_TopUp = [ (
            f"Charge at 0.3C until {V_max}V",
            f"Hold at {V_max}V until C/100",
            "Rest for 1 hours",
            "Rest for 10 s",   # add here to correct values of step_0p1C_CD 
            # 0.1C cycle 
            f"Discharge at 0.1C until {V_min} V",  
            "Rest for 3 hours",  
            f"Charge at 0.1C until {V_max} V",
            f"Hold at {V_max}V until C/100",
            "Rest for 1 hours",
            # 0.5C cycle 
            f"Discharge at 0.5C until {V_min} V",  
            "Rest for 3 hours",
            f"Charge at 0.5C until {V_max} V",
            f"Hold at {V_max}V until C/100",
            # Update 23-11-17: add one more 0.5C cycle to increase throughput capacity
            f"Discharge at 0.5C until {V_min} V",  
            "Rest for 3 hours",
            f"Charge at 0.5C until {V_max} V",
            f"Hold at {V_max}V until C/100",   
            "Rest for 3 hours",  
            ) ] 
        exp_RPT_No_TopUp = [ (
            "Rest for 10 s",   # add here to correct values of step_0p1C_CD
            "Rest for 10 s",   # add here to correct values of step_0p1C_CD
            "Rest for 1 hours", 
            "Rest for 10 s",   # add here to correct values of step_0p1C_CD
            # 0.1C cycle 
            f"Discharge at 0.1C until {V_min} V",  
            "Rest for 3 hours",  
            f"Charge at 0.1C until {V_max} V",
            f"Hold at {V_max}V until C/100",
            "Rest for 1 hours",
            # 0.5C cycle 
            f"Discharge at 0.5C until {V_min} V",  
            "Rest for 3 hours",
            f"Charge at 0.5C until {V_max} V",
            f"Hold at {V_max}V until C/100",
            # Update 23-11-17: add one more 0.5C cycle to increase throughput capacity
            f"Discharge at 0.5C until {V_min} V",  
            "Rest for 3 hours",
            f"Charge at 0.5C until {V_max} V",
            f"Hold at {V_max}V until C/100",   
            "Rest for 3 hours",  
            ) ] 
        exp_breakin_text = [ (
            # refill
            #"Rest for 10 s",   # add here to correct values of step_0p1C_CD
            #f"Hold at {V_max}V until C/100",
            #"Rest for 10 s", # Mark Ruihe change ad hoc setting for LFP 
            f"Discharge at 0.5C until {V_max-0.2}V", # start from discharge as it is easier for unbalanced cells
            f"Charge at 0.3C until {V_max}V",
            f"Hold at {V_max}V until C/100",
            "Rest for 1 hours", 
            # 0.1C cycle 
            f"Discharge at 0.1C until {V_min} V",  
            "Rest for 3 hours",  
            f"Charge at 0.1C until {V_max} V",
            f"Hold at {V_max}V until C/100",
            "Rest for 1 hours",
            # 0.5C cycle 
            f"Discharge at 0.5C until {V_min} V",  
            "Rest for 3 hours",
            f"Charge at 0.5C until {V_max} V",
            f"Hold at {V_max}V until C/100",
            # Update 23-11-17: add one more 0.5C cycle to increase throughput capacity
            f"Discharge at 0.5C until {V_min} V",  
            "Rest for 3 hours",
            f"Charge at 0.5C until {V_max} V",
            f"Hold at {V_max}V until C/100",   
            "Rest for 3 hours",  
            ) ] 
        
        # now for RPT - short version 
        exp_RPT_Need_TopUp_short = [ (
            f"Charge at 0.3C until {V_max}V",
            f"Hold at {V_max}V until C/100",
            "Rest for 10 s",
            "Rest for 10 s",   # add here to correct values of step_0p1C_CD 
            # 0.1C cycle 
            f"Discharge at 0.1C for 10 s",  
            "Rest for 10 s",  
            f"Charge at 0.1C for 10 s",
            f"Rest for 10 s",
            "Rest for 10 s",
            # 0.5C cycle 
            f"Discharge at 0.5C for 10 s",  
            "Rest for 10 s",
            f"Charge at 0.5C for 10 s",
            f"Rest for 10 s",
            ) ] 
        exp_RPT_No_TopUp_short = [ (
            "Rest for 10 s",   # add here to correct values of step_0p1C_CD
            "Rest for 10 s",   # add here to correct values of step_0p1C_CD
            "Rest for 10 s",
            "Rest for 10 s",   # add here to correct values of step_0p1C_CD 
            # 0.1C cycle 
            f"Discharge at 0.1C for 10 s",  
            "Rest for 10 s",  
            f"Charge at 0.1C for 10 s",
            f"Rest for 10 s",
            "Rest for 10 s",
            # 0.5C cycle 
            f"Discharge at 0.5C for 10 s",  
            "Rest for 10 s",
            f"Charge at 0.5C for 10 s",
            f"Rest for 10 s",
            ) ] 
        exp_breakin_text_short = [ (
            f"Discharge at 0.5C until {V_max-0.2}V", # start from discharge as it is easier for unbalanced cells
            f"Charge at 0.3C until {V_max}V",
            f"Hold at {V_max}V until C/100",
            "Rest for 1 hours", 
            # 0.1C cycle 
            f"Discharge at 0.1C for 10 s",  
            "Rest for 10 s",  
            f"Charge at 0.1C for 10 s",
            f"Rest for 10 s",
            "Rest for 10 s",
            # 0.5C cycle 
            f"Discharge at 0.5C for 10 s",  
            "Rest for 10 s",
            f"Charge at 0.5C for 10 s",
            f"Rest for 10 s",
            ) ] 
        
        exp_GITT_text = [ (
            "Rest for 5 minutes (1 minute period)",  
            "Rest for 1.2 seconds (0.1 second period)",  
            f"Discharge at C/2 for 4.8 minutes or until {V_min}V (0.1 second period)",
            "Rest for 1 hour", # (5 minute period)  
            ) ]
        exp_refill = [ (
            f"Charge at 0.3C until {V_max}V",
            f"Hold at {V_max}V until C/100",
            "Rest for 1 hours", 
            ) ] 
        if Exp_No == 1: 
            # adjust to 30% SOC, SOC before this step must be 100%
            charge_time_mins = 1 * 60 
            exp_adjust_before_age = [ (
                f"Discharge at 1C until {V_min} V",
                f"Hold at {V_min}V until C/100",
                "Rest for 4 hours", 
                f"Charge at 0.3C for {charge_time_mins} minutes or until {V_max}V",
                ) ] 
            # Update: 231219 for Exp-2, need to top up the cell after ageing and before C/10
            if self.global_config.Runshort == "GEM-2":
                exp_RPT_text = exp_RPT_Need_TopUp 
            else:
                exp_RPT_text = exp_RPT_Need_TopUp_short
        elif Exp_No == 2:
            # adjust to 85% SOC, SOC before this step must be 100%
            discharge_time_mins = 0.15* 60 * 4.86491/5
            exp_adjust_before_age = [ (
                f"Discharge at 1C for {discharge_time_mins} minutes or until {V_min}V", # discharge for 15%SOC
                "Rest for 3 hours", 
                ) ] 
            # Update: 231219 for Exp-2, need to top up the cell after ageing and before C/10
            if self.global_config.Runshort == "GEM-2":
                exp_RPT_text = exp_RPT_Need_TopUp 
            else:
                exp_RPT_text = exp_RPT_Need_TopUp_short
        elif Exp_No ==3:
            exp_adjust_before_age = [ (
                # just a place holder for now TODO 
                "Rest for 1 hours", 
                ) ] 
            if self.global_config.Runshort == "GEM-2":
                exp_RPT_text = exp_RPT_No_TopUp 
            else:
                exp_RPT_text = exp_RPT_No_TopUp_short
        else:
            exp_adjust_before_age = [ (
                # just a place holder for now TODO 
                "Rest for 1 hours", 
                ) ] 
            if self.global_config.Runshort == "GEM-2":
                exp_RPT_text = exp_RPT_No_TopUp 
            else:
                exp_RPT_text = exp_RPT_No_TopUp_short

        if self.global_config.Runshort == "GEM-2":
            pass
        else:
            exp_breakin_text = exp_breakin_text_short
        
        # step index for RPT
        step_0p1C_CD = 4; step_0p1C_CC = 6;   
        step_0p1C_RE = 5; step_0p5C_CD = 9;

        self.exp_config.exp_AGE_text = exp_AGE_text
        self.exp_config.step_AGE_CD = step_AGE_CD
        self.exp_config.step_AGE_CC = step_AGE_CC
        self.exp_config.step_AGE_CV = step_AGE_CV
        self.exp_config.exp_breakin_text = exp_breakin_text
        self.exp_config.exp_RPT_text = exp_RPT_text
        self.exp_config.exp_GITT_text = exp_GITT_text
        self.exp_config.exp_refill = exp_refill
        self.exp_config.exp_adjust_before_age = exp_adjust_before_age
        self.exp_config.step_0p1C_CD = step_0p1C_CD
        self.exp_config.step_0p1C_CC = step_0p1C_CC
        self.exp_config.step_0p1C_RE = step_0p1C_RE
        self.exp_config.step_0p5C_CD = step_0p5C_CD    
    
    

    # 
    def Create_folders(self):
        BasicPath_Save = self.path_config.BasicPath_Save
        Target = self.path_config.Target
        if not os.path.exists(BasicPath_Save + Target):
            os.mkdir(BasicPath_Save + Target)
        if not os.path.exists(BasicPath_Save + Target+"Mats"):
            os.mkdir(BasicPath_Save + Target +"Mats");
        if not os.path.exists(BasicPath_Save + Target+"Plots"):
            os.mkdir(BasicPath_Save + Target+"Plots");
        if not os.path.exists(BasicPath_Save + Target+"Excel"):
            os.mkdir(BasicPath_Save + Target+"Excel")

    #
    ## Update 23-05-18
    # Function to get cell average index from several cells in one T and one Exp
    def Get_Cell_Mean_1T_1Exp(self):
        Exp_No = self.exp_config.Exp_No
        Age_T = self.exp_config.Age_T

        Exp_Any_AllData = self.expData_config.Exp_Any_AllData
        Temp_Cell_Exp = self.expData_config.Temp_Cell_Exp
        Exp_temp_i_cell = Temp_Cell_Exp[str(int(Age_T))]

        # Get the cell with smallest X "Charge Throughput (A.h)" - to interpolate
        X_1_last = [];    X_5_last = [];
        for cell in Exp_temp_i_cell:
            df = Exp_Any_AllData[cell]["Extract Data"]
            xtemp = np.array(df["Charge Throughput (A.h)"])/1e3
            X_1_last.append(xtemp[-1])
            index_Res = df[df['0.1s Resistance (Ohms)'].le(10)].index
            xtemp2 = np.array(df["Charge Throughput (A.h)"][index_Res])/1e3
            X_5_last.append(xtemp2[-1])
        index_X1 = np.argmin(X_1_last) # find index of the smallest value
        index_X5 = np.argmin(X_5_last)
        # Get the interpolate ones
        Y_1_st = []; Y_2_st = []; Y_3_st = []; Y_4_st = []; Y_5_st = [];Y_6_st = [];
        df_t1 = Exp_Any_AllData[Exp_temp_i_cell[index_X1]]["Extract Data"]
        X_1_st = np.array(df_t1["Charge Throughput (A.h)"])/1e3
        df_t5 = Exp_Any_AllData[Exp_temp_i_cell[index_X5]]["Extract Data"]
        index_Res_5 = df_t5[df_t5['0.1s Resistance (Ohms)'].le(10)].index
        X_5_st = np.array(df_t5["Charge Throughput (A.h)"][index_Res_5])/1e3
        for cell in Exp_temp_i_cell:
            df = Exp_Any_AllData[cell]["Extract Data"]
            chThr_temp = np.array(df["Charge Throughput (A.h)"])/1e3
            df_DMA = Exp_Any_AllData[cell]["DMA"]["LLI_LAM"]
            Y_1_st.append( list(np.interp(X_1_st,chThr_temp,
                            np.array(df_DMA["SoH"])*100)) )
            Y_2_st.append( list(np.interp(X_1_st,chThr_temp,
                            np.array(df_DMA["LLI"])*100)) )
            Y_3_st.append( list(np.interp(X_1_st,chThr_temp,
                            np.array(df_DMA["LAM NE_tot"])*100)) )
            Y_4_st.append( list(np.interp(X_1_st,chThr_temp,
                            np.array(df_DMA["LAM PE"])*100)) )
            Y_6_st.append(list(
                np.interp(X_1_st, chThr_temp, 
                df["Age set average temperature (degC)"].astype(float))))

            index_Res = df[df['0.1s Resistance (Ohms)'].le(10)].index
            Y_5_st.append( list(np.interp(
                X_5_st,np.array(df["Charge Throughput (A.h)"][index_Res])/1e3,
                np.array(df["0.1s Resistance (Ohms)"][index_Res])*1e3, )) )
        # Get the mean values of the interpolate ones:
        def Get_Mean(X_1_st,Y_1_st):
            Y_1_st_avg = []
            for i in range(len(X_1_st)):
                sum=0
                for j in range(len(Y_1_st)):
                    sum += Y_1_st[j][i]
                Y_1_st_avg.append(sum/len(Y_1_st))
            return Y_1_st_avg   
        
        Y_1_st_avg = Get_Mean(X_1_st,Y_1_st)
        Y_2_st_avg = Get_Mean(X_1_st,Y_2_st)
        Y_3_st_avg = Get_Mean(X_1_st,Y_3_st)
        Y_4_st_avg = Get_Mean(X_1_st,Y_4_st)
        Y_5_st_avg = Get_Mean(X_5_st,Y_5_st)
        Y_6_st_avg = Get_Mean(X_1_st,Y_6_st)
        XY_pack = [
            X_1_st, X_5_st, Y_1_st_avg, Y_2_st_avg,
            Y_3_st_avg, Y_4_st_avg, Y_5_st_avg, Y_6_st_avg]
        self.expData_config.XY_pack = XY_pack

