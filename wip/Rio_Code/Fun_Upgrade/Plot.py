"""Functions for plotting"""


import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

def adjust_category_names(categories):
    return [
        '\n'.join(category.split()) if len(category) > 12 else category
        for category in categories
    ]

# plot inside the function:
def Plot_Cyc_RPT_4(task_result, config):
    # Unpack:
    my_dict_RPT = task_result.my_dict_RPT

    Exp_Any_AllData = config.expData_config.Exp_Any_AllData
    Temp_Cell_Exp = config.expData_config.Temp_Cell_Exp
    XY_pack = config.expData_config.XY_pack
    Exp_No = config.exp_config.Exp_No
    Plot_Exp = config.global_config.Plot_Exp
    R_from_GITT = config.global_config.R_from_GITT
    Scan_No = config.Scan_No
    Re_No = config.global_config.Re_No
    Age_T_in_K = config.global_config.Re_No
    model_options = config.model_config.model_options
    BasicPath_Save = config.path_config.BasicPath_Save
    Target = config.path_config.Target
    fs = config.global_config.fs
    dpi = config.global_config.dpi


    Num_subplot = 5;
    fig, axs = plt.subplots(2,3, figsize=(15,7.8),tight_layout=True)
    axs[0,0].plot(
        my_dict_RPT['Throughput capacity [kA.h]'], 
        my_dict_RPT['CDend SOH [%]'],     
        '-o', label="Scan=" + str(Scan_No) )
    axs[0,1].plot(
        my_dict_RPT['Throughput capacity [kA.h]'], 
        my_dict_RPT["CDend LLI [%]"],'-o', label="total LLI")
    if model_options.__contains__("lithium plating"):
        axs[0,1].plot(
            my_dict_RPT['Throughput capacity [kA.h]'], 
            my_dict_RPT["CDend LLI lithium plating [%]"],'--o', label="LiP")
    if model_options.__contains__("SEI"):
        axs[0,1].plot(
            my_dict_RPT['Throughput capacity [kA.h]'], 
            my_dict_RPT["CDend LLI SEI [%]"] ,'--o', label="SEI")
    if model_options.__contains__("SEI on cracks"):
        axs[0,1].plot(
            my_dict_RPT['Throughput capacity [kA.h]'], 
            my_dict_RPT["CDend LLI SEI on cracks [%]"] ,
            '--o', label="SEI-on-cracks")
    axs[0,2].plot(
        my_dict_RPT["Throughput capacity [kA.h]"], 
        my_dict_RPT["CDend LAM_ne [%]"],     '-o', ) 
    axs[1,0].plot(
        my_dict_RPT["Throughput capacity [kA.h]"], 
        my_dict_RPT["CDend LAM_pe [%]"],     '-o',  ) 
    axs[1,1].plot(
        my_dict_RPT["Throughput capacity [kA.h]"], 
        np.array(my_dict_RPT["Res_midSOC"]),     '-o', ) 
    axs[1,2].plot(
        my_dict_RPT["Throughput capacity [kA.h]"][1:], 
        np.array(my_dict_RPT["avg_Age_T"][1:]),     '-o', ) 
    # Plot Charge Throughput (A.h) vs SOH
    color_exp     = [0, 0, 0, 0.3]; marker_exp     = "v";
    color_exp_Avg = [0, 0, 0, 0.7]; marker_exp_Avg = "s";
    if Exp_No in list(np.arange(1,6)) and int(Age_T_in_K- 273.15) in [10,25,40]:
        Exp_temp_i_cell = Temp_Cell_Exp[str(int(Age_T_in_K- 273.15))]
    else:
        Exp_temp_i_cell = "nan"
        Plot_Exp = False

    if Plot_Exp == True:
        for cell in Exp_temp_i_cell:
            df = Exp_Any_AllData[cell]["Extract Data"]
            chThr_temp = np.array(df["Charge Throughput (A.h)"])/1e3
            df_DMA = Exp_Any_AllData[cell]["DMA"]["LLI_LAM"]
            axs[0,0].plot(
                chThr_temp,np.array(df_DMA["SoH"])*100,
                color=color_exp,marker=marker_exp,label=f"Cell {cell}") 
            axs[0,1].plot(
                chThr_temp,np.array(df_DMA["LLI"])*100,
                color=color_exp,marker=marker_exp,label=f"Cell {cell}")  
            axs[0,2].plot(
                chThr_temp,np.array(df_DMA["LAM NE_tot"])*100,
                color=color_exp,marker=marker_exp, )
            axs[1,0].plot(
                chThr_temp,np.array(df_DMA["LAM PE"])*100,
                color=color_exp,marker=marker_exp,)
            # update 230312- plot resistance here
            # Exp_1_AllData["A"]["Extract Data"]["0.1s Resistance (Ohms)"]
            index_Res = df[df['0.1s Resistance (Ohms)'].le(10)].index
            axs[1,1].plot(
                #df["Days of degradation"][index_Res],
                np.array(df["Charge Throughput (A.h)"][index_Res])/1e3,
                np.array(df["0.1s Resistance (Ohms)"][index_Res])*1e3,
                color=color_exp,marker=marker_exp)
            axs[1,2].plot(
                chThr_temp[1:],
                np.array(df["Age set average temperature (degC)"][1:]).astype(float),
                color=color_exp,marker=marker_exp,)
        # Update 230518: Plot Experiment Average - at 1 expeirment and 1 temperature
        [X_1_st,X_5_st,Y_1_st_avg,Y_2_st_avg,
            Y_3_st_avg,Y_4_st_avg,Y_5_st_avg,Y_6_st_avg]  = XY_pack
        axs[0,0].plot(
            X_1_st,Y_1_st_avg,color=color_exp_Avg,
            marker=marker_exp_Avg,label=f"Exp-Avg") 
        axs[0,1].plot(
            X_1_st,Y_2_st_avg,color=color_exp_Avg,
            marker=marker_exp_Avg,label=f"Exp-Avg")  
        axs[0,2].plot(
            X_1_st,Y_3_st_avg,color=color_exp_Avg,
            marker=marker_exp_Avg, )
        axs[1,0].plot(
            X_1_st,Y_4_st_avg,
            color=color_exp_Avg,marker=marker_exp_Avg,)
        axs[1,1].plot(
            X_5_st,Y_5_st_avg,
            color=color_exp_Avg,marker=marker_exp_Avg)
        axs[1,2].plot(
            X_1_st[1:],Y_6_st_avg[1:],
            color=color_exp_Avg,marker=marker_exp_Avg,)
    axs[0,0].set_ylabel("SOH %")
    axs[0,1].set_ylabel("LLI %")
    axs[0,2].set_ylabel("LAM NE %")
    axs[1,0].set_ylabel("LAM PE %")
    axs[1,1].set_ylabel(r"Lump resistance [m$\Omega$]")
    axs[1,2].set_ylabel(r"Avg age T [$^\circ$C]")
    axs[0,2].set_xlabel("Charge Throughput (kA.h)")
    axs[1,2].set_xlabel("Charge Throughput (kA.h)")
    axf = axs.flatten()
    for i in range(0,6):
        labels = axf[i].get_xticklabels() + axf[i].get_yticklabels(); 
        [label.set_fontname('DejaVu Sans') for label in labels]
        axf[i].tick_params(labelcolor='k', labelsize=fs, width=1);del labels
    axs[1,1].ticklabel_format(style='sci', axis='x', scilimits=(-1e-2,1e-2))
    axs[0,0].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)
    axs[0,1].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)
    fig.suptitle(
        f"Scan_{Scan_No}-Exp-{Exp_No}-{str(int(Age_T_in_K- 273.15))}"
        +r"$^\circ$C - Summary", fontsize=fs+2)
    plt.savefig(
        BasicPath_Save + Target+    "Plots/" +  
        f"0_Scan_{Scan_No}_Re_{Re_No}-Exp-{Exp_No}-{str(int(Age_T_in_K- 273.15))}degC Summary.png", dpi=dpi)
    plt.close()  # close the figure to save RAM

    if model_options.__contains__("SEI on cracks"):
        """ Num_subplot = 2;
        fig, axs = plt.subplots(1,Num_subplot, figsize=(12,4.8),tight_layout=True)
        axs[0].plot(my_dict_RPT['Throughput capacity [kA.h]'], my_dict_RPT["CDend X-averaged total SEI on cracks thickness [m]"],     '-o', label="Scan=" + str(Scan_No) )
        axs[1].plot(my_dict_RPT['Throughput capacity [kA.h]'], my_dict_RPT["CDend X-averaged negative electrode roughness ratio"],'-o', label="Scan=" + str(Scan_No) )
        axs[0].set_ylabel("SEI on cracks thickness [m]")
        axs[1].set_ylabel("Roughness ratio")
        for i in range(0,Num_subplot):
            axs[i].set_xlabel("Charge Throughput (kA.h)")
            labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
            axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
            axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)
        axs[0].set_title("X-avg tot Neg SEI on cracks thickness",   fontdict={'family':'DejaVu Sans','size':fs+1})
        axs[1].set_title("X-avg Neg roughness ratio",   fontdict={'family':'DejaVu Sans','size':fs+1})
        plt.savefig(BasicPath_Save + Target+"Plots/" +
            f"Scan_{Scan_No}_Re_{Re_No}-Exp-{Exp_No}-{str(int(Age_T_in_K- 273.15))}degC - Cracks related_Scan.png", dpi=dpi)
        plt.close()  # close the figure to save RAM """

        Num_subplot = 2;
        fig, axs = plt.subplots(1,Num_subplot, figsize=(12,4.8),tight_layout=True)
        axs[0].plot(my_dict_RPT['Throughput capacity [kA.h]'], 
            my_dict_RPT["CDend Negative electrode capacity [A.h]"][0]
            -
            my_dict_RPT["CDend Negative electrode capacity [A.h]"],'-o',label="Neg Scan=" + str(Scan_No))
        axs[1].plot(my_dict_RPT['Throughput capacity [kA.h]'], 
            my_dict_RPT["CDend Positive electrode capacity [A.h]"][0]
            -
            my_dict_RPT["CDend Positive electrode capacity [A.h]"],'-^',label="Pos Scan=" + str(Scan_No))
        """ axs[0].plot(
            my_dict_RPT['Throughput capacity [kA.h]'], 
            my_dict_RPT["CDend X-averaged total SEI on cracks thickness [m]"],                  
            '-o',label="Scan="+ str(Scan_No)) """
        for i in range(0,2):
            axs[i].set_xlabel("Charge Throughput (kA.h)")
            labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
            axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
            axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)
        #axs[0].set_ylabel("SEI on cracks thickness [m]")
        #axs[0].set_title("CDend X-avg tot SEI on cracks thickness",   fontdict={'family':'DejaVu Sans','size':fs+1})
        for i in range(0,2):
            axs[i].set_xlabel("Charge Throughput (kA.h)")
            axs[i].set_ylabel("Capacity [A.h]")
            labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
            axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
            axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)
            axs[i].set_title("LAM of Neg and Pos",   fontdict={'family':'DejaVu Sans','size':fs+1})
        plt.savefig(BasicPath_Save + Target+"Plots/" +
            f"Scan_{Scan_No}_Re_{Re_No}-Exp-{Exp_No}-{str(int(Age_T_in_K- 273.15))}degC LAM-IR.png", dpi=dpi)
        plt.close()  # close the figure to save RAM
    Num_subplot = 2;
    fig, axs = plt.subplots(1,Num_subplot, figsize=(8,3.2),tight_layout=True)
    axs[0].plot(my_dict_RPT['Throughput capacity [kA.h]'], my_dict_RPT["CDsta Positive electrode stoichiometry"] ,'-o',label="Start" )
    axs[0].plot(my_dict_RPT['Throughput capacity [kA.h]'], my_dict_RPT["CDend Positive electrode stoichiometry"] ,'-^',label="End" )
    axs[1].plot(my_dict_RPT['Throughput capacity [kA.h]'], my_dict_RPT["CDsta Negative electrode stoichiometry"],'-o',label="Start" )
    axs[1].plot(my_dict_RPT['Throughput capacity [kA.h]'], my_dict_RPT["CDend Negative electrode stoichiometry"],'-^',label="End" )
    for i in range(0,2):
        axs[i].set_xlabel("Charge Throughput (kA.h)")
        axs[i].set_ylabel("Stoichiometry")
        labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
        axs[i].ticklabel_format(style='sci', axis='x', scilimits=(-1e-2,1e-2))
        axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
        axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)     
    axs[0].set_title("Neg Sto. range (Dis)",   fontdict={'family':'DejaVu Sans','size':fs+1})
    axs[1].set_title("Pos Sto. range (Dis)",   fontdict={'family':'DejaVu Sans','size':fs+1})
    plt.savefig(BasicPath_Save + Target+"Plots/"+
        f"Scan_{Scan_No}_Re_{Re_No}-Exp-{Exp_No}-{str(int(Age_T_in_K- 273.15))}degC SOC_RPT_dis.png", dpi=dpi) 
    plt.close()  # close the figure to save RAM
    # update 230518: plot resistance in C/2 discharge:
    N_RPT = len(my_dict_RPT["Res_full"])
    colormap_i = mpl.cm.get_cmap("gray", 14) 
    fig, axs = plt.subplots(figsize=(4,3.2),tight_layout=True)
    for i in range(N_RPT):
        axs.plot(
            my_dict_RPT["SOC_Res"][i], my_dict_RPT["Res_full"][i] ,
            color=colormap_i(i),marker="o",  label=f"RPT {i}" )
    if R_from_GITT: 
        axs.set_xlabel("SOC-GITT %")
        axs.set_ylabel(r'Res GITT (m$\Omega$)')
        # axs.legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)     
        axs.set_title("Res during GITT Dis",   fontdict={'family':'DejaVu Sans','size':fs+1})
    else:
        axs.set_xlabel("SOC-C/2 %")
        axs.set_ylabel(r'Res C/2 (m$\Omega$)')
        # axs.legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)     
        axs.set_title("Res during C/2 Dis",   fontdict={'family':'DejaVu Sans','size':fs+1})
    plt.savefig(BasicPath_Save + Target+ "Plots/"+
        f"Scan_{Scan_No}_Re_{Re_No}-Exp-{Exp_No}-{str(int(Age_T_in_K- 273.15))}degC Res_full.png", dpi=dpi) 
    plt.close()  # close the figure to save RAM

    return



def Plot_DMA_Dec(task_result, config):
    # Unpack:
    my_dict_RPT = task_result.my_dict_RPT
    Scan_No = config.Scan_No
    Re_No = config.global_config.Re_No
    Age_T_in_K = config.global_config.Re_No
    BasicPath_Save = config.path_config.BasicPath_Save
    Target = config.path_config.Target
    fs = config.global_config.fs
    dpi = config.global_config.dpi


    fig, ax = plt.subplots(figsize=(6,5),tight_layout=True) 
    categories = [
        'SEI', 'SEI on cracks', 'Li plating', 
        'LAM (cracking and dry-out)']
    adjusted_categories = adjust_category_names(categories)
    values = [
        my_dict_RPT["CDend LLI SEI [%]"][-1], 
        my_dict_RPT["CDend LLI SEI on cracks [%]"][-1], 
        my_dict_RPT["CDend LLI lithium plating [%]"][-1] , 
        my_dict_RPT["CDend LLI due to LAM [%]"][-1]   ]
    plt.bar(adjusted_categories, values ,width=0.45 )
    plt.ylabel('LLI %')
    fig.suptitle(
        f"Scan_{Scan_No}_Re_{Re_No}-{str(int(Age_T_in_K- 273.15))}"
        +r"$^\circ$C - LLI break down", fontsize=fs+2)
    plt.savefig(
        BasicPath_Save + Target+    "Plots/" +  
        f"0_Scan_{Scan_No}_Re_{Re_No}-"
        f"{str(int(Age_T_in_K- 273.15))}degC "
        "LLI break down.png", dpi=dpi)
    plt.close()  # close the figure to save RAM

    fig, axs = plt.subplots(1,2, figsize=(12,4.5),tight_layout=True) 
    width = 0.45
    categories = ['Dry out', 'Stress',]
    values_ne = [
        my_dict_RPT["LAM_to_Dry [%] end"],
        my_dict_RPT["LAM_to_Crack_NE [%] end"] ]
    values_pe = [
        my_dict_RPT["LAM_to_Dry [%] end"],
        my_dict_RPT["LAM_to_Crack_PE [%] end"] ]
    axs[0].bar(categories, values_ne, width )
    axs[1].bar(categories, values_pe, width )
    axs[0].set_ylabel("LAM NE %")
    axs[1].set_ylabel("LAM PE %")
    fig.suptitle(
        f"Scan_{Scan_No}_Re_{Re_No}-{str(int(Age_T_in_K- 273.15))}"
        +r"$^\circ$C - LAM break down", fontsize=fs+2)
    plt.savefig(
        BasicPath_Save + Target+    "Plots/" +  
        f"0_Scan_{Scan_No}_Re_{Re_No}-"
        f"{str(int(Age_T_in_K- 273.15))}degC"
        " LAM break down.png", dpi=dpi)
    plt.close()  # close the figure to save RAM
    return



def Plot_HalfCell_V(task_result, config):
    # Unpack:
    my_dict_RPT = task_result.my_dict_RPT
    my_dict_AGE = task_result.my_dict_AGE
    colormap = config.global_config.colormap
    Exp_No = config.exp_config.Exp_No
    Scan_No = config.Scan_No
    Re_No = config.global_config.Re_No
    Age_T_in_K = config.global_config.Re_No
    BasicPath_Save = config.path_config.BasicPath_Save
    Target = config.path_config.Target
    fs = config.global_config.fs
    dpi = config.global_config.dpi

    #~~~~~~~~~~~~~~~~~ plot RPT 
    def inFun_Plot(my_dict,str_jj):
        fig, axs = plt.subplots(3,2, figsize=(12,10),tight_layout=True)
        Str_Front= ["CD ","CC ",]
        Str_Back = ["Cathode potential [V]","Terminal voltage [V]","Anode potential [V]",]
        for j in range(2): 
            for i in range(3):
                Time = my_dict[Str_Front[j]+"Time [h]"]
                Num_Lines = len(Time)
                cmap = mpl.cm.get_cmap(colormap, Num_Lines) # cmap(i)
                for k in range(Num_Lines):
                    Y = my_dict[Str_Front[j]+Str_Back[i]] [k]
                    axs[i,j].plot( Time[k] , Y, color = cmap(k)   )
                if j == 0 :
                    axs[i,j].set_ylabel(
                        Str_Back[i])
            axs[2,j].set_xlabel("Time [h]")
        axs[0,0].set_title("During Discharge",   fontdict={'family':'DejaVu Sans','size':fs+1})
        axs[0,1].set_title("During Charge",   fontdict={'family':'DejaVu Sans','size':fs+1})
        fig.suptitle(
            f"Scan {str(Scan_No)}-Exp-{Exp_No}-{str(int(Age_T_in_K-273.15))}"
            +r"$^\circ$C"+f" - Half cell Potential ({str_jj})", fontsize=fs+2)
        plt.savefig(BasicPath_Save + Target+"Plots/" +
            f"Scan_{Scan_No}_Re_{Re_No}-Exp-{Exp_No}-{str(int(Age_T_in_K- 273.15))}degC"
            f" Half cell Potential ({str_jj}).png", dpi=dpi) 
        plt.close() 
    inFun_Plot(my_dict_RPT,"RPT")
    inFun_Plot(my_dict_AGE,"AGE")
    return 



def Plot_Loc_AGE_4(task_result,config):

    my_dict_AGE = task_result.my_dict_AGE
    Exp_No = config.exp_config.Exp_No
    Scan_No = config.Scan_No
    Re_No = config.global_config.Re_No
    Age_T_in_K = config.global_config.Re_No
    BasicPath_Save = config.path_config.BasicPath_Save
    Target = config.path_config.Target
    fs = config.global_config.fs
    dpi = config.global_config.dpi
    
    Num_subplot = 2;
    fig, axs = plt.subplots(1,Num_subplot, figsize=(8,3.2),tight_layout=True)
    axs[0].plot(my_dict_AGE["x [m]"], my_dict_AGE["CDend Porosity"][0],'-o',label="First")
    axs[0].plot(my_dict_AGE["x [m]"], my_dict_AGE["CDend Porosity"][-1],'-^',label="Last"  )
    axs[1].plot(
        my_dict_AGE["x_n [m]"],
        my_dict_AGE["CDend Negative electrode reaction overpotential [V]"][0],'-o',label="First" )
    axs[1].plot(
        my_dict_AGE["x_n [m]"],
        my_dict_AGE["CDend Negative electrode reaction overpotential [V]"][-1],'-^',label="Last" )

    axs[0].set_xlabel("Dimensional Cell thickness")
    axs[1].set_xlabel("Dimensional Neg thickness")
    axs[0].set_title("Porosity",   fontdict={'family':'DejaVu Sans','size':fs+1})
    axs[1].set_title("Neg electrode reaction overpotential",   fontdict={'family':'DejaVu Sans','size':fs+1})
    axs[0].set_ylabel("Porosity")
    axs[1].set_ylabel("Overpotential [V]")
    for i in range(0,2):
        labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
        axs[i].ticklabel_format(style='sci', axis='x', scilimits=(-1e-2,1e-2))
        axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
        axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)    
    plt.savefig(BasicPath_Save + Target+"Plots/" +
        f"Scan_{Scan_No}_Re_{Re_No}-Exp-{Exp_No}-{str(int(Age_T_in_K- 273.15))}degC Por Neg_S_eta.png", dpi=dpi) 
    plt.close()  # close the figure to save RAM

    """ update: disable cracking related things just because when sei-on-crack is false it will easily give errors; 
    if model_options.__contains__("SEI on cracks"):
        Num_subplot = 2;
        fig, axs = plt.subplots(1,Num_subplot, figsize=(8,3.2),tight_layout=True)
        axs[0].plot(my_dict_AGE["x_n [m]"], my_dict_AGE["CDend Negative electrode roughness ratio"][0],'-o',label="First")
        axs[0].plot(my_dict_AGE["x_n [m]"], my_dict_AGE["CDend Negative electrode roughness ratio"][-1],'-^',label="Last"  )
        axs[1].plot(my_dict_AGE["x_n [m]"], my_dict_AGE["CDend Total SEI on cracks thickness [m]"][0],'-o',label="First" )
        axs[1].plot(my_dict_AGE["x_n [m]"], my_dict_AGE["CDend Total SEI on cracks thickness [m]"][-1],'-^',label="Last" )
        axs[0].set_title("Tot Neg SEI on cracks thickness",   fontdict={'family':'DejaVu Sans','size':fs+1})
        axs[1].set_title("Neg roughness ratio",   fontdict={'family':'DejaVu Sans','size':fs+1})
        axs[0].set_ylabel("SEI on cracks thickness [m]")
        axs[1].set_ylabel("Roughness ratio")
        for i in range(0,2):
            axs[i].set_xlabel("Dimensional Neg thickness")
            labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
            axs[i].ticklabel_format(style='sci', axis='x', scilimits=(-1e-2,1e-2))
            axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
            axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)    
        plt.savefig(BasicPath_Save + Target+"Plots/" +
            f"Scan_{Scan_No}_Re_{Re_No}-Exp-{Exp_No}-{str(int(Age_T_in_K- 273.15))}degC Cracks related spatial.png", dpi=dpi) 
        plt.close()  # close the figure to save RAM """
    Num_subplot = 2;
    fig, axs = plt.subplots(1,Num_subplot, figsize=(8,3.2),tight_layout=True)
    axs[0].plot(my_dict_AGE["x [m]"], my_dict_AGE["CDend Electrolyte concentration [mol.m-3]"][0],'-o',label="First")
    axs[0].plot(my_dict_AGE["x [m]"], my_dict_AGE["CDend Electrolyte concentration [mol.m-3]"][-1],'-^',label="Last"  )
    axs[1].plot(my_dict_AGE["x [m]"], my_dict_AGE["CDend Electrolyte potential [V]"][0],'-o',label="First" )
    axs[1].plot(my_dict_AGE["x [m]"], my_dict_AGE["CDend Electrolyte potential [V]"][-1],'-^',label="Last" )
    axs[0].set_title("Electrolyte concentration",   fontdict={'family':'DejaVu Sans','size':fs+1})
    axs[1].set_title("Electrolyte potential",   fontdict={'family':'DejaVu Sans','size':fs+1})
    axs[0].set_ylabel("Concentration [mol.m-3]")
    axs[1].set_ylabel("Potential [V]")
    for i in range(0,2):
        axs[i].set_xlabel("Dimensional Cell thickness")
        labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
        axs[i].ticklabel_format(style='sci', axis='x', scilimits=(-1e-2,1e-2))
        axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
        axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)    
    plt.savefig(BasicPath_Save + Target+"Plots/" +
        f"Scan_{Scan_No}_Re_{Re_No}-Exp-{Exp_No}-{str(int(Age_T_in_K- 273.15))}degC Electrolyte concentration and potential.png", dpi=dpi)
    plt.close()  # close the figure to save RAM
    Num_subplot = 2;
    fig, axs = plt.subplots(1,Num_subplot, figsize=(8,3.2),tight_layout=True)
    axs[0].plot(my_dict_AGE["x [m]"], my_dict_AGE["CDend Electrolyte diffusivity [m2.s-1]"][0],'-o',label="First")
    axs[0].plot(my_dict_AGE["x [m]"], my_dict_AGE["CDend Electrolyte diffusivity [m2.s-1]"][-1],'-^',label="Last"  )
    axs[1].plot(my_dict_AGE["x [m]"], my_dict_AGE["CDend Electrolyte conductivity [S.m-1]"][0],'-o',label="First" )
    axs[1].plot(my_dict_AGE["x [m]"], my_dict_AGE["CDend Electrolyte conductivity [S.m-1]"][-1],'-^',label="Last" )
    axs[0].set_title("Electrolyte diffusivity",   fontdict={'family':'DejaVu Sans','size':fs+1})
    axs[1].set_title("Electrolyte conductivity",   fontdict={'family':'DejaVu Sans','size':fs+1})
    axs[0].set_ylabel("Diffusivity [m2.s-1]")
    axs[1].set_ylabel("Conductivity [S.m-1]")
    for i in range(0,2):
        axs[i].set_xlabel("Dimensional Cell thickness")
        labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
        axs[i].ticklabel_format(style='sci', axis='x', scilimits=(-1e-2,1e-2))
        axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
        axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)    
    plt.savefig(BasicPath_Save + Target+"Plots/" +
        f"Scan_{Scan_No}_Re_{Re_No}-Exp-{Exp_No}-{str(int(Age_T_in_K- 273.15))}degC Electrolyte diffusivity and conductivity.png", dpi=dpi)
    plt.close()  # close the figure to save RAM
    return

def Plot_Dryout(task_result, config):
    # Unpack:
    mdic_dry = task_result.mdic_dry
    ce_EC_0 = mdic_dry["c_EC_r_new_All"][0] # be careful
    Cyc_Update_Index = task_result.my_dict_RPT["Cyc_Update_Index"] 
    Exp_No = config.exp_config.Exp_No
    Scan_No = config.Scan_No
    Re_No = config.global_config.Re_No
    Age_T_in_K = config.global_config.Re_No
    BasicPath_Save = config.path_config.BasicPath_Save
    Target = config.path_config.Target
    fs = config.global_config.fs
    dpi = config.global_config.dpi
    
    CeEC_All =np.full(np.size(mdic_dry["Ratio_CeEC_All"]),ce_EC_0); 
    for i in range(1,np.size(mdic_dry["Ratio_CeEC_All"])):
        for k in range(0,i):
            CeEC_All[i] *= mdic_dry["Ratio_CeEC_All"][k]
    mdic_dry["CeEC_All"]= CeEC_All.tolist()

    Num_subplot = 3;
    fig, axs = plt.subplots(1,Num_subplot, figsize=(12,3.2),tight_layout=True)
    axs[0].plot(Cyc_Update_Index, mdic_dry["Vol_EC_consumed_All"],'-o',label="EC consumed")
    axs[0].plot(Cyc_Update_Index, mdic_dry["Vol_Elely_need_All"],'-.',label="Elely needed")
    axs[0].plot(Cyc_Update_Index, mdic_dry["Vol_Elely_add_All"],'-s',label="Elely added")
    axs[0].plot(Cyc_Update_Index, mdic_dry["Vol_Pore_decrease_All"],'--',label="Pore decreased")
    axs[1].plot(Cyc_Update_Index, mdic_dry["Vol_Elely_Tot_All"],     '-o', label="Total electrolyte in cell" )
    axs[1].plot(Cyc_Update_Index, mdic_dry["Vol_Elely_JR_All"],     '--', label="Total electrolyte in JR" )
    axs[1].plot(Cyc_Update_Index, mdic_dry["Vol_Pore_tot_All"],     '-s', label="Total pore in JR" )
    axs[2].plot(Cyc_Update_Index, mdic_dry["Ratio_Dryout_All"],     '-s', label="Dry out ratio" )
    axs[2].set_ylabel("Ratio")
    axs[2].set_xlabel("Cycle number")
    axs[2].set_title("Dry out ratio")
    for i in range(0,2):
        axs[i].set_xlabel("Cycle number")
        axs[i].set_ylabel("Volume [mL]")
        axs[i].set_title("Volume",   fontdict={'family':'DejaVu Sans','size':fs+1})
        labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
        axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
        axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)    
    plt.savefig(BasicPath_Save + Target+"Plots/" +
        f"Scan_{Scan_No}_Re_{Re_No}-Exp-{Exp_No}-{str(int(Age_T_in_K- 273.15))}degC Volume_total.png", 
        dpi=dpi)
    plt.close()  # close the figure to save RAM
    return

