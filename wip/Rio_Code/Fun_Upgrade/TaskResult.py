""" Define TaskResult class to get output """

import os
import numpy as np

class TaskResult:
    def __init__(
            self, midc_merge=None,
            Sol_RPT=[],Sol_AGE=[],

            Call_Breakin=None, Flag_Breakin = None,
            str_error_Breakin = None,
            DeBug_List_Breakin=None,
            
            Call_Age=None,      Flag_AGE = True, 
            str_error_AGE = "Empty", 
            DeBug_List_AGE = "Empty",

            Call_RPT = None,    Flag_partial_AGE = False,
            str_error_RPT = "Empty", 
            DeBug_List_RPT = "Empty", 

            DryOut_List = None,
            mpe_all = None,
            Pass_Fail = None,
            Keys_error = None,
            
            ):
        # get output keys - customized    
        self.midc_merge = midc_merge 
        self.Sol_RPT = Sol_RPT 
        self.Sol_AGE = Sol_AGE  

        self.Call_Breakin = Call_Breakin 
        self.Flag_Breakin = Flag_Breakin 
        self.str_error_Breakin = str_error_Breakin 
        self.DeBug_List_Breakin = DeBug_List_Breakin

        self.Call_Age = Call_Age 
        self.Flag_AGE = Flag_AGE 
        self.str_error_AGE = str_error_AGE
        self.DeBug_List_AGE = DeBug_List_AGE

        self.Call_RPT = Call_RPT 
        self.Flag_partial_AGE = Flag_partial_AGE 
        self.str_error_RPT = str_error_RPT
        self.DeBug_List_RPT = DeBug_List_RPT

        self.DryOut_List = DryOut_List
        self.mpe_all = mpe_all
        self.Pass_Fail = Pass_Fail
        self.Keys_error = Keys_error




    def Get_Output_Keys(self,model_options):
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
            "CD Anode potential [V]",    # self defined
            "CD Cathode potential [V]",  # self defined
            "CC Time [h]",
            "CC Terminal voltage [V]",
            "CC Anode potential [V]",    # self defined
            "CC Cathode potential [V]",  # self defined
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
            "CCend Negative electrode surface potential difference [V]",
            "CCend SEI film overpotential [V]",
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
            "CDend Negative electrode surface potential difference [V]",
            "CDend SEI film overpotential [V]",
            "CDend Electrolyte diffusivity [m2.s-1]",
            "CDend Electrolyte conductivity [S.m-1]",
        ]
        keys_tim_AGE = [
            # default: CD
            "CD Time [h]",
            "CD Terminal voltage [V]",
            "CD Anode potential [V]",    # self defined
            "CD Cathode potential [V]",  # self defined
            
            "CC Time [h]",
            "CC Terminal voltage [V]",
            "CC Anode potential [V]",    # self defined
            "CC Cathode potential [V]",  # self defined
        ]
        keys_Cycle_bt_RPT = []
        keys_all_RPT = [keys_loc_RPT,keys_tim_RPT,keys_cyc_RPT]
        keys_all_AGE = [keys_loc_AGE,keys_tim_AGE,keys_Cycle_bt_RPT]
        self.keys_all_RPT = keys_all_RPT
        self.keys_all_AGE = keys_all_AGE

    def Initialize_my_dict(self):
        keys_all_RPT = self.keys_all_RPT
        keys_all_AGE = self.keys_all_AGE
        my_dict_RPT = {}
        for keys in keys_all_RPT:
            for key in keys:
                my_dict_RPT[key]=[];
        my_dict_AGE = {}; 
        for keys in keys_all_AGE:
            for key in keys:
                my_dict_AGE[key]=[]
        my_dict_RPT["Cycle_RPT"] = []
        my_dict_RPT["Res_full"] = []
        my_dict_RPT["Res_midSOC"] = []
        my_dict_RPT["SOC_Res"] = []
        my_dict_AGE["Cycle_AGE"] = []
        my_dict_RPT["avg_Age_T"] = [] 
        my_dict_AGE["avg_Age_T"] = [] 
        my_dict_RPT["Cyc_Update_Index"] = []
        # one number that adds on during calculation
        my_dict_AGE["Agecycle_count"] = 0 
        self.my_dict_RPT = my_dict_RPT
        self.my_dict_AGE = my_dict_AGE
    
    def Add_Dryout_Dict(self,mdic_dry):
        self.mdic_dry = mdic_dry

    def Update_mdic_dry(self):
        DryOut_List = self.DryOut_List
        mdic_dry = self.mdic_dry
        [
            Vol_EC_consumed, 
            Vol_Elely_need, 
            Test_V, 
            Test_V2, 
            Vol_Elely_add, 
            Vol_Elely_Tot_new, 
            Vol_Elely_JR_new, 
            Vol_Pore_tot_new, 
            Vol_Pore_decrease, 
            c_e_r_new, c_EC_r_new,
            Ratio_Dryout, Ratio_CeEC_JR, 
            Ratio_CeLi_JR,
            Width_new, ]= DryOut_List
        mdic_dry["Vol_Elely_Tot_All"].append(Vol_Elely_Tot_new*1e6)          
        mdic_dry["Vol_Elely_JR_All"].append(Vol_Elely_JR_new*1e6)     
        mdic_dry["Vol_Pore_tot_All"].append(Vol_Pore_tot_new*1e6)           
        mdic_dry["Ratio_CeEC_All"].append(Ratio_CeEC_JR);                      
        mdic_dry["Ratio_CeLi_All"].append(Ratio_CeLi_JR);             
        mdic_dry["Ratio_Dryout_All"].append(Ratio_Dryout);
        mdic_dry["Vol_EC_consumed_All"].append(Vol_EC_consumed*1e6)        
        mdic_dry["Vol_Elely_need_All"].append(Vol_Elely_need*1e6)    
        mdic_dry["Width_all"].append(Width_new)
        mdic_dry["Vol_Elely_add_All"].append(Vol_Elely_add*1e6)           
        mdic_dry["Vol_Pore_decrease_All"].append(Vol_Pore_decrease*1e6)
        mdic_dry["Test_V_All"].append(Test_V*1e6)
        mdic_dry["Test_V2_All"].append(Test_V2*1e6)
        mdic_dry["c_e_r_new_All"].append(c_e_r_new)
        mdic_dry["c_EC_r_new_All"].append(c_EC_r_new)

        self.mdic_dry = mdic_dry

   

