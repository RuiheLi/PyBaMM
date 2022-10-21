import pybamm;import pandas as pd   ;import numpy as np;import os;
import matplotlib.pyplot as plt;import os;#import imageio
from scipy.io import savemat,loadmat;from pybamm import constants,exp;
import matplotlib as mpl; fs=17; # or we can set import matplotlib.pyplot as plt then say 'mpl.rc...'

import openpyxl
import traceback

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

def Plot_Single_Static(Sol,str,cycle, step, Para_scan,BasicPath , Target,Save,colormap,fs,dpi):
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
