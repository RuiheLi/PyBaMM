import pybamm;import pandas as pd   ;import numpy as np;import os;
import matplotlib.pyplot as plt;import os;#import imageio
from scipy.io import savemat,loadmat;from pybamm import constants,exp;
import matplotlib as mpl; fs=17; # or we can set import matplotlib.pyplot as plt then say 'mpl.rc...'

import openpyxl
import traceback

from pybamm import tanh,exp,sqrt
def electrolyte_conductivity_Valoen2005Constant(c_e,c_EC, T):# Mark Ruihe change
    # mol/m3 to molar
    c_e = c_e / 1000
    sigma = (c_e <= 4.5) * (
        (1e-3 / 1e-2) * (
        c_e
        * (
            (-10.5 + 0.0740 * T - 6.96e-5 * T ** 2)
            + c_e * (0.668 - 0.0178 * T + 2.80e-5 * T ** 2)
            + c_e ** 2 * (0.494 - 8.86e-4 * T)
        )
        ** 2
    )) + (c_e > 4.5) *  (
        (1e-3 / 1e-2) * (
        4.5
        * (
            (-10.5 + 0.0740 * T - 6.96e-5 * T ** 2)
            + 4.5 * (0.668 - 0.0178 * T + 2.80e-5 * T ** 2)
            + 4.5 ** 2 * (0.494 - 8.86e-4 * T)
        )
        ** 2
    ))
    # mS/cm to S/m
    return sigma
def electrolyte_conductivity_Valoen2005Constant_EC_Haya(c_e,c_EC, T):
    # mol/m3 to molar
    c_e = c_e / 1000
    sigma = (c_e <= 4.5) * (
        (1e-3 / 1e-2) * (
        c_e
        * (
            (-10.5 + 0.0740 * T - 6.96e-5 * T ** 2)
            + c_e * (0.668 - 0.0178 * T + 2.80e-5 * T ** 2)
            + c_e ** 2 * (0.494 - 8.86e-4 * T)
        )
        ** 2
    )) + (c_e > 4.5) *  (
        (1e-3 / 1e-2) * (
        4.5
        * (
            (-10.5 + 0.0740 * T - 6.96e-5 * T ** 2)
            + 4.5 * (0.668 - 0.0178 * T + 2.80e-5 * T ** 2)
            + 4.5 ** 2 * (0.494 - 8.86e-4 * T)
        )
        ** 2
    ))
    a=1.092; b=-6.497e-6; c=-0.7877; d=-0.0004808
    ratio= (
        a*exp(b*c_EC)+c*exp(d*c_EC) )
    return sigma*ratio
def electrolyte_conductivity_Valoen2005Constant_ECtanh500_1(c_e,c_EC, T):# Mark Ruihe change
    # mol/m3 to molar
    c_e = c_e / 1000
    sigma = (c_e <= 4.5) * (
        (1e-3 / 1e-2) * (
        c_e
        * (
            (-10.5 + 0.0740 * T - 6.96e-5 * T ** 2)
            + c_e * (0.668 - 0.0178 * T + 2.80e-5 * T ** 2)
            + c_e ** 2 * (0.494 - 8.86e-4 * T)
        )
        ** 2
    )) + (c_e > 4.5) *  (
        (1e-3 / 1e-2) * (
        4.5
        * (
            (-10.5 + 0.0740 * T - 6.96e-5 * T ** 2)
            + 4.5 * (0.668 - 0.0178 * T + 2.80e-5 * T ** 2)
            + 4.5 ** 2 * (0.494 - 8.86e-4 * T)
        )
        ** 2
    ))
    coff = 1
    ratio = ( (1-coff)+ coff/2 + coff/2 *  tanh((c_EC-4500*0.5)/500))
    return sigma*ratio


# DEFINE my callback:
class RioCallback(pybamm.callbacks.Callback):
    def __init__(self, logfile=None):
        self.logfile = logfile
        self.success  = True
        if logfile is None:
            # Use pybamm's logger, which prints to command line
            self.logger = pybamm.logger
        else:
            # Use a custom logger, this will have its own level so set it to the same
            # level as the pybamm logger (users can override this)
            self.logger = pybamm.get_new_logger(__name__, logfile)
            self.logger.setLevel(pybamm.logger.level)
    
    def on_experiment_error(self, logs):
        self.success  = False
    def on_experiment_infeasible(self, logs):
        self.success  = False

def write_excel_xlsx(path, sheet_name, value):
    import numpy as np
    index = len(value)
    workbook = openpyxl.Workbook()  # 新建工作簿（默认有一个sheet？）
    sheet = workbook.active  # 获得当前活跃的工作页，默认为第一个工作页
    sheet.title = sheet_name  # 给sheet页的title赋值
    for i in range(0, index):
        for j in range(0, len(value[i])):
            sheet.cell(row=i + 1, column=j + 1, value=str(value[i][j]))  # 行，列，值 这里是从1开始计数的
    workbook.save(path)  # 一定要保存
    print("Successfully create a excel file")

# valid only for single parameter scan:
def PlotDynamics(Sol,str,Para_scan,BasicPath , Target,Save,fs):
    font = {'family' : 'DejaVu Sans','size': fs}
    mpl.rc('font', **font)
    # plot overall figures
    label = [str+f"={Para_i}" for Para_i in Para_scan]
    output_variables1 = [
        "Terminal voltage [V]",   
        "Discharge capacity [A.h]",
        "EC concentration [mol.m-3]",
        "Electrolyte concentration [mol.m-3]",
        "Li+ flux [mol.m-2.s-1]",
        "EC flux [mol.m-2.s-1]",
    ]
    quick_plot = pybamm.QuickPlot(
        [sol for sol in Sol], 
        output_variables1,label,
        variable_limits='fixed',time_unit='hours',n_rows=2,
        figsize = (12,8)) #     spatial_unit='mm',
    quick_plot.dynamic_plot();
    if Save == 'True':
        quick_plot.create_gif(
            number_of_images=10, duration=2,
            output_filename=BasicPath + Target+"All cases - overview.gif")

    output_variables2 = [
        "Electrolyte concentration",
        "Minus div Li+ flux",
        "Li+ source term",
        "Minus div Li+ flux by diffusion",
        "Minus div Li+ flux by migration",
        "Minus div Li+ flux by solvent",
        #"Li+ source term refill",
    ]
    quick_plot = pybamm.QuickPlot(
        [sol for sol in Sol], 
        output_variables2,label,
        variable_limits='tight',time_unit='hours',n_rows=2,
        figsize = (12,9)) #     spatial_unit='mm',
    quick_plot.dynamic_plot();
    if Save == 'True':
        quick_plot.create_gif(
            number_of_images=10, duration=2,
            output_filename=BasicPath + Target+"All cases - c(Li) and breakdown.gif")
    else:
        pass

    output_variables3 = [
        "EC concentration",
        "EC source term (SEI)",
        "Minus div EC flux",
        "Minus div EC flux by diffusion",
        "Minus div EC flux by migration",
        "Minus div EC flux by Li+",
        #"Li+ source term refill",
    ]
    quick_plot = pybamm.QuickPlot(
        [sol for sol in Sol], 
        output_variables3,label,
        variable_limits='fixed',time_unit='hours',n_rows=2,
        figsize = (12,9)) #     spatial_unit='mm',
    quick_plot.dynamic_plot();
    if Save == 'True':
        quick_plot.create_gif(
            number_of_images=10, duration=2,
            output_filename=BasicPath + Target+"All cases - c(EC) and breakdown.gif")
    else:
        pass
    return 


# Plot a pair of loc dependent varibles - different cycles
def Plot_Loc_Var( key_all, my_dict,colormap ,fs): # for my_dict only
    font = {'family' : 'DejaVu Sans','size': fs}
    mpl.rc('font', **font)
    Num_subplot = len(key_all); # must have 2+ keys
    fig, axs = plt.subplots(1,Num_subplot, figsize=(7*Num_subplot,5),tight_layout=True)
    for i in range(0,Num_subplot):
        cmap_i = mpl.cm.get_cmap(colormap, len(my_dict[ key_all[i]] ) ) 
        if 'Negative' in key_all[i] or 'negative' in key_all[i]:
            x_loc = "x_n [m]";
        elif 'Positive' in key_all[i] or 'positive' in key_all[i]:
            x_loc = "x_p [m]";
        elif 'Seperator' in key_all[i] or 'seperator' in key_all[i]:
            x_loc = "x_s [m]";
        else:
            x_loc = "x [m]";
        for j in range(0,len(my_dict[ key_all[i] ])):
            axs[i].plot(
                my_dict[x_loc], 
                my_dict[ key_all[i] ][j],'-',
                color=cmap_i(j),)
            axs[i].set_title(key_all[i] ,   fontdict={'family':'DejaVu Sans','size':fs-1})
            #axs[1].set_ylabel("Potential [V]",   fontdict={'family':'DejaVu Sans','size':fs})
            axs[i].set_xlabel(x_loc,   fontdict={'family':'DejaVu Sans','size':fs})
            
            labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
            
            axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
            axs[i].ticklabel_format( axis='x', style='sci',scilimits=[-0.01,0.01], useOffset=None, useLocale=None, useMathText=None)
            #axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)  
    return fig, axs 

# Plot a pair of loc dependent varibles - within one step
def Plot_Loc_Var_sol( sol,x_loc_all, key_all, cycle, step,colormap,fs): # for initial solution object
    font = {'family' : 'DejaVu Sans','size': fs}
    mpl.rc('font', **font)
    Num_subplot = len(key_all); # must have 2+ keys
    fig, axs = plt.subplots(1,Num_subplot, figsize=(4*Num_subplot,5),tight_layout=True)
    for i in range(0,Num_subplot):
        x_loc=x_loc_all[i]; key=key_all[i];
        LinesNmum = len(sol.cycles[cycle].steps[step][key].entries[0,:] )
        cmap_i = mpl.cm.get_cmap(colormap, LinesNmum) 
        for j in range(0,LinesNmum):
            axs[i].plot(
                sol.cycles[cycle].steps[step][x_loc].entries[:,0], 
                sol.cycles[cycle].steps[step][key].entries[:,j], '-',
                color=cmap_i(j),)
            axs[i].set_title(key ,   fontdict={'family':'DejaVu Sans','size':fs-1})
            #axs[1].set_ylabel("Potential [V]",   fontdict={'family':'DejaVu Sans','size':fs})
            axs[i].set_xlabel(x_loc,   fontdict={'family':'DejaVu Sans','size':fs})
            
            labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
            
            axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
            axs[i].ticklabel_format( axis='x', style='sci',scilimits=[-0.01,0.01], useOffset=None, useLocale=None, useMathText=None)
            #axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)  
    return  fig, axs 

def Plot_Loc_Var_sol_6( sol,x_loc_all, key_all, cycle, step,colormap,fs): # for initial solution object
    font = {'family' : 'DejaVu Sans','size'   : fs}
    mpl.rc('font', **font)
    Num_subplot = len(key_all); # must have 2+ keys
    if Num_subplot <= 3:
        fig, axs = plt.subplots(1,Num_subplot, figsize=(2.5*Num_subplot,5),tight_layout=True)
        for i in range(0,Num_subplot):
            x_loc=x_loc_all[i]; key=key_all[i];
            LinesNmum = len(sol.cycles[cycle].steps[step][key].entries[0,:] )
            cmap_i = mpl.cm.get_cmap(colormap, LinesNmum) 
            for j in range(0,LinesNmum):
                axs[i].plot(
                    sol.cycles[cycle].steps[step][x_loc].entries[:,0], 
                    sol.cycles[cycle].steps[step][key].entries[:,j], '-',
                    color=cmap_i(j),)
                axs[i].set_title(key ,   fontdict={'family':'DejaVu Sans','size':fs-1})
                #axs[1].set_ylabel("Potential [V]",   fontdict={'family':'DejaVu Sans','size':fs})
                axs[i].set_xlabel(x_loc,   fontdict={'family':'DejaVu Sans','size':fs})
                
                labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
                
                axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
                axs[i].ticklabel_format( axis='x', style='sci',
                scilimits=[-0.01,0.01], useOffset=None, 
                useLocale=None, useMathText=None)
                #axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)  
    elif Num_subplot >3:
        fig, axs = plt.subplots(2,3, figsize=(2.5*Num_subplot,8.5),tight_layout=True)
        Plot_Count = 0
        for m in range(0,2):
            for n in range(0,3):
                x_loc=x_loc_all[Plot_Count]; 
                key=key_all[Plot_Count]; Plot_Count +=1;
                temp_L_Num= len(sol.cycles[cycle].steps[step][key].entries[0,:] )
                xx = np.arange(0,temp_L_Num,int(np.rint(temp_L_Num/100)));
                xx = xx.tolist()
                if not xx[-1]==temp_L_Num-1:
                    xx.append(temp_L_Num-1)
                cmap_i = mpl.cm.get_cmap(colormap, len(xx)) 
                for index_j,j in zip(  range(0,len(xx)), xx):
                    len_1 = len(sol.cycles[cycle].steps[step][x_loc].entries[:,0])
                    len_2 = len(sol.cycles[cycle].steps[step][key].entries[:,j])
                    if len_1>len_2:
                        axs[m,n].plot(
                            sol.cycles[cycle].steps[step][x_loc].entries[0:len_2,0], 
                            sol.cycles[cycle].steps[step][key].entries[:,j], '-',
                            color=cmap_i(index_j),)
                    elif len_1<len_2:
                        axs[m,n].plot(
                            sol.cycles[cycle].steps[step][x_loc].entries[:,0], 
                            sol.cycles[cycle].steps[step][key].entries[0:len_1,j], '-',
                            color=cmap_i(index_j),)
                    else:    
                        axs[m,n].plot(
                            sol.cycles[cycle].steps[step][x_loc].entries[:,0], 
                            sol.cycles[cycle].steps[step][key].entries[:,j], '-',
                            color=cmap_i(index_j),)
                    axs[m,n].set_title(key ,   fontdict={'family':'DejaVu Sans','size':fs-1})
                    #axs[1].set_ylabel("Potential [V]",   fontdict={'family':'DejaVu Sans','size':fs})
                    axs[m,n].set_xlabel(x_loc,   fontdict={'family':'DejaVu Sans','size':fs})
                    
                    labels = axs[m,n].get_xticklabels() + axs[m,n].get_yticklabels(); 
                    [label.set_fontname('DejaVu Sans') for label in labels]
                    
                    axs[m,n].tick_params(
                        labelcolor='k', labelsize=fs, width=1) ;  del labels;
                    axs[m,n].ticklabel_format( 
                        axis='x', style='sci',
                        scilimits=[-0.01,0.01], useOffset=None, 
                        useLocale=None, useMathText=None)
    return  fig, axs 

# Fig. 1
def Plot_Fig_1(Full_cycle,my_dict_AGE,
    BasicPath, Target,   Scan_i,  fs,  dpi):
    font = {'family' : 'DejaVu Sans','size': fs}
    mpl.rc('font', **font)
    fig, axs = plt.subplots(2,3, figsize=(15,8.5),tight_layout=True)
    axs[0,0].plot(Full_cycle, my_dict_AGE["Discharge capacity [A.h]"] ,'-o',)
    axs[0,0].set_ylabel("Capacity [A.h]",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[0,0].set_title("Discharge capacity",   fontdict={'family':'DejaVu Sans','size':fs+1})

    axs[0,1].plot(Full_cycle, my_dict_AGE["CDend Loss of capacity to SEI [A.h]"] ,'-o',)
    axs[0,1].set_ylabel("Capacity [A.h]",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[0,1].set_title("Loss of capacity to SEI",   fontdict={'family':'DejaVu Sans','size':fs+1})

    axs[0,2].plot(Full_cycle, my_dict_AGE["CDend Local ECM resistance [Ohm]"] ,'-o',)
    axs[0,2].set_ylabel("Resistance [Ohm]",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[0,2].set_title("Local ECM resistance",   fontdict={'family':'DejaVu Sans','size':fs+1})

    axs[1,0].plot(Full_cycle, my_dict_AGE["CDsta Positive electrode SOC"] ,'-o',label="Start" )
    axs[1,0].plot(Full_cycle, my_dict_AGE["CDend Positive electrode SOC"] ,'-^',label="End" )
    axs[1,1].plot(Full_cycle, my_dict_AGE["CDsta Negative electrode SOC"],'-o',label="Start" )
    axs[1,1].plot(Full_cycle, my_dict_AGE["CDend Negative electrode SOC"],'-^',label="End" )
    axs[1,0].set_ylabel("SOC",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[1,0].set_title("Pos SOC range (Dis)",   fontdict={'family':'DejaVu Sans','size':fs+1})
    axs[1,1].set_ylabel("SOC",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[1,1].set_title("Neg SOC range (Dis)",   fontdict={'family':'DejaVu Sans','size':fs+1})


    axs[1,2].plot(Full_cycle, my_dict_AGE["CDend Negative electrode capacity [A.h]"],'-o',label="Neg" )
    axs[1,2].plot(Full_cycle, my_dict_AGE["CDend Positive electrode capacity [A.h]"],'-^',label="Pos" )
    axs[1,2].set_ylabel("Capacity [A.h]",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[1,2].set_title("Electrode capacity",   fontdict={'family':'DejaVu Sans','size':fs+1})

    for i in range(0,2):
        for j in range(0,3):
            axs[i,j].set_xlabel("Cycle numbers",   fontdict={'family':'DejaVu Sans','size':fs})
            
            labels = axs[i,j].get_xticklabels() + axs[i,j].get_yticklabels(); 
            [label.set_fontname('DejaVu Sans') for label in labels]
            axs[i,j].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
            axs[i,j].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)     
    plt.savefig(
        BasicPath + Target+f"{Scan_i}th Scan/" + 
        "Fig. 1 - Time based overall change.png", dpi=dpi)

        
def Plot_Loc_Var_2( Full_cycle, key_all, my_dict,fs): # for my_dict only
    font = {'family' : 'DejaVu Sans','size': fs}
    mpl.rc('font', **font)
    Num_subplot = len(key_all); # must have 2+ keys
    fig, axs = plt.subplots(3,3, figsize=(15,12),tight_layout=True)
    count_i = 0; 
    for i in range(0,3):
        for j in range(0,3):
            key_ij = key_all[count_i]; count_i += 1
            if 'Negative' in key_ij or 'negative' in key_ij:
                x_loc = "x_n [m]";
            elif 'Positive' in key_ij or 'positive' in key_ij:
                x_loc = "x_p [m]";
            elif 'Seperator' in key_ij or 'seperator' in key_ij:
                x_loc = "x_s [m]";
            else:
                x_loc = "x [m]";
            X_Len = min(len(my_dict[x_loc]),len(my_dict[ key_ij ][0]))
            #print(x_loc,X_Len)
            axs[i,j].plot(my_dict[x_loc][0:X_Len], my_dict[ key_ij ][0][0:X_Len],'-o',label="1st cycle")
            axs[i,j].plot(my_dict[x_loc][0:X_Len], my_dict[ key_ij ][-1][0:X_Len],'-^',label=f"{Full_cycle[-1]}th cycle")
            axs[i,j].set_title(key_ij ,   fontdict={'family':'DejaVu Sans','size':fs-3})
            #axs[1].set_ylabel("Potential [V]",   fontdict={'family':'DejaVu Sans','size':fs})
            axs[i,j].set_xlabel(x_loc,   fontdict={'family':'DejaVu Sans','size':fs})
            
            labels = axs[i,j].get_xticklabels() + axs[i,j].get_yticklabels(); 
            [label.set_fontname('DejaVu Sans') for label in labels]
            
            axs[i,j].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
            axs[i,j].ticklabel_format( 
                axis='x', style='sci',scilimits=[-0.01,0.01], 
                useOffset=None, useLocale=None, useMathText=None)
            axs[i,j].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)  
    return fig, axs 


def Plot_Last_Single_Step(
    sol,cycle,step,BasicPath, Target,Scan_i,
    index_cyc,Save,colormap,fs,dpi):

    font = {'family' : 'DejaVu Sans','size': fs}
    mpl.rc('font', **font)
    # 1st plot Li+ and flux
    i =Scan_i;
    fig, axs= Plot_Loc_Var_sol_6(
        sol,
        ["x [m]","x [m]","x [m]","x [m]","x [m]","x [m]",], 
        ["Electrolyte concentration",
        "Minus div Li+ flux",
        "Li+ source term",
        "Minus div Li+ flux by diffusion",
        "Minus div Li+ flux by migration",
        "Minus div Li+ flux by solvent",], 
        cycle,step,colormap,fs)
    if Save == 'True':
        plt.savefig(
            BasicPath + 
            Target +     f"{i}th Scan/" +
            f"after {index_cyc}th cycle-cycle={cycle}, step={step}, c(Li+) and flux.png", dpi=dpi)
    else:
        pass
    
    # 2nd plot: for EC and flux
    fig, axs = Plot_Loc_Var_sol_6(
        sol,
        ["x [m]","x [m]","x [m]","x [m]","x [m]","x [m]",], 
        ["EC concentration",
        "c(EC) over c(Li+)",
        "Minus div EC flux",
        "Minus div EC flux by diffusion",
        "Minus div EC flux by migration",
        "Minus div EC flux by Li+",], 
        cycle,step,colormap,fs)
    if Save == 'True':
        plt.savefig(
            BasicPath + 
            Target+ f"{i}th Scan/" +
            f"after {index_cyc}th cycle-cycle={cycle}, step={step}, c(EC) and flux.png", dpi=dpi)
    else:
        pass
    # 3rd plot: porosity, potential
    fig, axs = Plot_Loc_Var_sol_6(
        sol,
        ["x_n [m]","x_p [m]","x [m]","x [m]","x [m]","x [m]",], 
        [
            "Negative electrode porosity",
            "Positive electrode potential [V]",
            "Electrolyte current density [A.m-2]",
            "Electrolyte potential [V]",
            "Electrolyte diffusivity [m2.s-1]",
            "Electrolyte conductivity [S.m-1]",
        ], 
        cycle,step,colormap,fs)
    if Save == 'True':
        plt.savefig(
            BasicPath + 
            Target+f"{i}th Scan/" +
            f"after {index_cyc}th cycle-cycle={cycle}, step={step}, Porosity, elely and potential.png", dpi=dpi)
    else:
        pass

    return


def Plot_Single_Static(Sol,str,cycle, step, 
    Para_scan,BasicPath , Target,Save,colormap,fs,dpi):

    font = {'family' : 'DejaVu Sans','size': fs}
    mpl.rc('font', **font)
    for sol, Para_i in zip(Sol,Para_scan):
        # 1st plot
        fig, axs= Plot_Loc_Var_sol(
            sol,
            ["x [m]","x [m]","x [m]",], 
            ["Electrolyte concentration",
            "Minus div Li+ flux",
            "Li+ source term",], 
            cycle,step,colormap,fs)
        if Save == 'True':
            plt.savefig(
                BasicPath + 
                Target +    str +
                f"={Para_i} - cycle={cycle}, step={step}, c(Li+) and flux.png", dpi=dpi)
        else:
            pass
        
        # 2nd plot
        fig, axs= Plot_Loc_Var_sol(
            sol,
            ["x [m]","x [m]","x [m]",], 
            ["Minus div Li+ flux by diffusion",
            "Minus div Li+ flux by migration",
            "Minus div Li+ flux by solvent",], 
            cycle,step,colormap,fs)
        if Save == 'True':
            plt.savefig(
                BasicPath + 
                Target+str+
                f"={Para_i} - cycle={cycle}, step={step}, Li+ flux breakdown.png", dpi=dpi)
        else:
            pass
        

        fig, axs = Plot_Loc_Var_sol(
            sol,
            ["x [m]","x [m]","x [m]",], 
            ["EC concentration",
            "c(EC) over c(Li+)",
            "Minus div EC flux",], 
            cycle,step,colormap,fs)
        if Save == 'True':
            plt.savefig(
                BasicPath + 
                Target+str+
                f"={Para_i} - cycle={cycle}, step={step}, c(EC) and flux.png", dpi=dpi)
        else:
            pass

        fig, axs = Plot_Loc_Var_sol(
            sol,
            ["x [m]","x [m]","x [m]",], 
            ["Minus div EC flux by diffusion",
            "Minus div EC flux by migration",
            "Minus div EC flux by Li+",], 
            cycle,step,colormap,fs)
        if Save == 'True':
            plt.savefig(
                BasicPath + 
                Target+str+
                f"={Para_i} - cycle={cycle}, step={step}, EC flux breakdown.png", dpi=dpi)
        else:
            pass
    return


import random
'''
有一定概率出错的函数
'''
def may_cause_error():
	if (random.random() > 0.8):
	    return 1 / 0
	else:
	    return 0

def recursive_scan(mylist,kvs, key_list, acc):
    # 递归终止条件
    if len(key_list) == 0:
        mylist.append(acc.copy())   # copy.deepcopy(acc) 如果value是个list，就要deep copy了
        return mylist
    # 继续递归
    k = key_list[0]
    for v in kvs[k]:
        acc[k] = v
        # print(v)
        recursive_scan(mylist,kvs, key_list[1:], acc)

# this function is to initialize the para with a known dict
def Para_init(Para_dict):
    Para_dict_used = Para_dict.copy();
    Para_0=pybamm.ParameterValues(Para_dict_used["Para_Set"]  )
    Para_dict_used.pop("Para_Set")

    if Para_dict_used.__contains__("Total ageing cycles"):
        Total_Cycles = Para_dict_used["Total ageing cycles"]  
        Para_dict_used.pop("Total ageing cycles")
    if Para_dict_used.__contains__("SaveAsList"):
        SaveAsList = Para_dict_used["SaveAsList"]  
        Para_dict_used.pop("SaveAsList")
    if Para_dict_used.__contains__("Ageing cycles between RPT"):
        Cycle_bt_RPT = Para_dict_used["Ageing cycles between RPT"]  
        Para_dict_used.pop("Ageing cycles between RPT")
    if Para_dict_used.__contains__("Update cycles for ageing"):
        Update_Cycles = Para_dict_used["Update cycles for ageing"]  
        Para_dict_used.pop("Update cycles for ageing")
    if Para_dict_used.__contains__("Cycles within RPT"):
        RPT_Cycles = Para_dict_used["Cycles within RPT"]  
        Para_dict_used.pop("Cycles within RPT")
    if Para_dict_used.__contains__("Ageing temperature"):
        Temper_i = Para_dict_used["Ageing temperature"]  
        Para_dict_used.pop("Ageing temperature")
    if Para_dict_used.__contains__("RPT temperature"):
        Temper_RPT = Para_dict_used["RPT temperature"]  
        Para_dict_used.pop("RPT temperature")
    if Para_dict_used.__contains__("Particle mesh points"):
        mesh_par = Para_dict_used["Particle mesh points"]  
        Para_dict_used.pop("Particle mesh points")
    if Para_dict_used.__contains__("Exponential mesh stretch"):
        submesh_strech = Para_dict_used["Exponential mesh stretch"]  
        Para_dict_used.pop("Exponential mesh stretch")
    else:
        submesh_strech = "nan";
    if Para_dict_used.__contains__("Model option"):
        model_options = Para_dict_used["Model option"]  
        Para_dict_used.pop("Model option")
    if Para_dict_used.__contains__("Func Electrolyte conductivity [S.m-1]"):
        Para_0.update({
            "Electrolyte conductivity [S.m-1]": 
            eval(Para_dict_used["Func Electrolyte conductivity [S.m-1]"])})  
        Para_dict_used.pop("Func Electrolyte conductivity [S.m-1]")
    if Para_dict_used.__contains__("Func Electrolyte diffusivity [m2.s-1]"):
        Para_0.update({
            "Electrolyte diffusivity [m2.s-1]": 
            eval(Para_dict_used["Func Electrolyte diffusivity [m2.s-1]"])})
        Para_dict_used.pop("Func Electrolyte diffusivity [m2.s-1]")

    if Para_dict_used.__contains__("Initial Neg SOC"):
        c_Neg1SOC_in = (
            Para_0["Maximum concentration in negative electrode [mol.m-3]"]
            *Para_dict_used["Initial Neg SOC"]  )
        Para_0.update(
            {"Initial concentration in negative electrode [mol.m-3]":
            c_Neg1SOC_in})
        Para_dict_used.pop("Initial Neg SOC")
    if Para_dict_used.__contains__("Initial Pos SOC"):    
        c_Pos1SOC_in = (
            Para_0["Maximum concentration in positive electrode [mol.m-3]"]
            *Para_dict_used["Initial Pos SOC"]  )
        Para_0.update(
            {"Initial concentration in positive electrode [mol.m-3]":
            c_Pos1SOC_in})
        Para_dict_used.pop("Initial Pos SOC")

    CyclePack = [ 
        Total_Cycles,SaveAsList,Temper_i,model_options];
    
    for key, value in Para_dict_used.items():
        # risk: will update parameter that doesn't exist, 
        # so need to make sure the name is right 
        Para_0.update({key: value},check_already_exists=False)
    return CyclePack,Para_0

def GetSol_dict (my_dict, keys_all, Sol, 
    cycle_no,step_CD, step_CC , step_RE, step_CV ):

    [keys_loc,keys_tim,keys_cyc]  = keys_all;
    # get time_based variables:
    if len(keys_tim): 
        for key in keys_tim:
            step_no = eval("step_{}".format(key[0:2]))
            if key[3:] == "Time [h]":
                my_dict[key].append  (
                    Sol.cycles[cycle_no].steps[step_no][key[3:]].entries
                    -
                    Sol.cycles[cycle_no].steps[step_no][key[3:]].entries[0])
            else:
                #print("Solve up to Step",step_no)
                my_dict[key].append  (
                    Sol.cycles[cycle_no].steps[step_no][key[3:]].entries)
    # get cycle_step_based variables:
    if len(keys_cyc): 
        for key in keys_cyc:
            if key in ["Discharge capacity [A.h]"]:
                step_no = step_CD;
                my_dict[key].append  (
                    Sol.cycles[cycle_no].steps[step_no][key].entries[-1]
                    - 
                    Sol.cycles[cycle_no].steps[step_no][key].entries[0])
            elif key[0:5] in ["CDend","CCend","CVend","REend",
                            "CDsta","CCsta","CVsta","REsta",]:
                step_no = eval("step_{}".format(key[0:2]))
                if key[2:5] == "sta":
                    my_dict[key].append  (
                        Sol.cycles[cycle_no].steps[step_no][key[6:]].entries[0])
                elif key[2:5] == "end":
                    my_dict[key].append  (
                        Sol.cycles[cycle_no].steps[step_no][key[6:]].entries[-1])
                
    # get location_based variables:
    if len(keys_loc): 
        for key in keys_loc:
            if key in ["x_n [m]","x [m]","x_s [m]","x_p [m]"]:
                #print("These variables only add once")
                if not len(my_dict[key]):   # special: add only once
                    my_dict[key] = Sol[key].entries[:,-1]
            elif key[0:5] in ["CDend","CCend","CVend","REend",
                            "CDsta","CCsta","CVsta","REsta",]:      
                #print("These variables add multiple times")
                step_no = eval("step_{}".format(key[0:2]))
                if key[2:5] == "sta":
                    my_dict[key].append  (
                        Sol.cycles[cycle_no].steps[step_no][key[6:]].entries[:,0])
                elif key[2:5] == "end":
                    my_dict[key].append  (
                        Sol.cycles[cycle_no].steps[step_no][key[6:]].entries[:,-1])
    return my_dict 


# Run model   # mark:
def Run_P3_model(
    index_xlsx, Para_dict_i,   Path_pack , 
    keys_all_AGE,   Exp_AGE_List, exp_index_pack ):

    ModelTimer = pybamm.Timer()
    Para_dict_old = Para_dict_i.copy();
    count_i = int(index_xlsx);

    # Un-pack data:
    [cycle_no,step_AGE_CD,
        step_AGE_CC,step_AGE_CV, ] = exp_index_pack
    [exp_AGE,exp_AGE,exp_AGE_2, # may be deleted later
        exp_AGE_CD,exp_AGE_CC,exp_AGE_CV] = Exp_AGE_List
    [BasicPath,Target,
        book_name_xlsx,sheet_name_xlsx,] = Path_pack

    ##### Initialise Para_0 and model 
    CyclePack,Para_0 = Para_init(Para_dict_i)
    [Total_Cycles,SaveAsList,
        Temper_i,model_options] = CyclePack
    model_0 = pybamm.lithium_ion.DFN(options=model_options)
    str_model_options = str(model_options)
    # add electrolyte properties as variables, 
    # only when they are not constant
    c_e = model_0.variables["Electrolyte concentration [mol.m-3]"]
    T = model_0.variables["Cell temperature [K]"]
    c_EC = model_0.variables["EC concentration [mol.m-3]"]
    model_0.variables["c(EC) over c(Li+)"] = c_EC / c_e
    if not Para_dict_i.__contains__(
        "Electrolyte conductivity [S.m-1]"):
        model_0.variables["Electrolyte conductivity [S.m-1]"] =(
            Para_0['Electrolyte conductivity [S.m-1]'](c_e,c_EC, T))
    if not Para_dict_i.__contains__(
        "Electrolyte diffusivity [m2.s-1]"):
        model_0.variables["Electrolyte diffusivity [m2.s-1]"] =(
            Para_0['Electrolyte diffusivity [m2.s-1]'](c_e,c_EC, T))
    

    # Define experiment- for ageing only, NOT RPT
    str_exp_All = [];Experiment_All = [];
    for SaveAs_i, exp_i in zip(SaveAsList,Exp_AGE_List):
        str_exp_All.append(str(exp_i))
        Experiment_All.append(
            pybamm.Experiment( exp_i * int(SaveAs_i)  ) )


    #####Important: index canot be pre-determined anymore! ######


    # initialize my_dict for outputs
    my_dict_AGE = {}; 
    for keys in keys_all_AGE:
        for key in keys:
            my_dict_AGE[key]=[];	
    ###### Biggest loop:  #####
    Succ_Cyc = []; Sol_All = []; 
    i = 0;step_switch = 0 # initialize the loop
    while (i < Total_Cycles):
        print('try to run %d cycles' % SaveAsList[step_switch])
        try:
            if i==0: # 1st time or never succeed, run from scratch:
                rioCall = RioCallback()  # define callback
                sim_0    = pybamm.Simulation(
                    model_0,
                    experiment = Experiment_All[step_switch],
                    parameter_values = Para_0,
                    solver = pybamm.CasadiSolver(),
                    #var_pts=var_pts,
                    #submesh_types=submesh_types
                    ) 
                sol_0    = sim_0.solve(
                    calc_esoh=False,
                    save_at_cycles = SaveAsList[step_switch],
                    callbacks=rioCall)
            else: # succeeded before, 
                model_new = model_old.set_initial_conditions_from(
                    sol_old, inplace=False)
                rioCall = RioCallback()  # define callback
                sim_new    = pybamm.Simulation(
                    model_new,        
                    experiment =  Experiment_All[step_switch],
                    parameter_values = Para_0,
                    solver = pybamm.CasadiSolver(),
                    #var_pts=var_pts,
                    #submesh_types=submesh_types
                    ) 
                sol_new    = sim_new.solve(
                    calc_esoh=False,
                    save_at_cycles = SaveAsList[step_switch],
                    callbacks=rioCall)   
            if rioCall.success == False:
                1/0
        except:
            print('Failed and shorten cycles')
            step_switch += 1
            if (step_switch >= len(SaveAsList)):
                print('Exit as no options left')
                str_error = traceback.format_exc() # final errors 
                print('Finally finish %d cycles' % i)  
                break
        else:        
            if i == 0: 
                model_old = model_0; sol_old = sol_0    
            else: 
                model_old = model_new; sol_old = sol_new
                del model_new,sol_new
            ### Should do early post-prosessing and delete to 
            ### save RAM in the near future, not now 
            # because how to post-processing depends on step_switch
            Sol_All.append(sol_old)
            Succ_Cyc.append(SaveAsList[step_switch])
            i += SaveAsList[step_switch]
            print('Succeed! Now it is the %dth cycles' % i)  
            if step_switch <= 2:
                step_switch += 0
            elif step_switch > 2 and step_switch < len(SaveAsList)-1:
                step_switch += 1
                print('Succeed a single step and switch to next step normally')
            else:
                print('Finish last single step and Exit as no options left')
                print('Finally finish %d cycles' % i)  
                break
    return Sol_All,Succ_Cyc


