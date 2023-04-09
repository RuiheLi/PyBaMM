# Concentration, LJP and overpotential 
Crate_index = -1
fig, axs = plt.subplots(1,2, figsize=(9.3,3.2),tight_layout=True)
step_sd = SD_Crate['Sol_All'][Crate_index].cycles[0].steps[1]
axs[0].plot(
    step_sd['Time [s]'].entries-step_sd['Time [s]'].entries[0], 
    step_sd["X-averaged battery concentration overpotential [V]"].entries,
    color=Colors[0],linestyle=LS[0],label=r"Single-High $D_\times$") 
axs[1].plot(
    step_sd['Time [s]'].entries-step_sd['Time [s]'].entries[0], 
    step_sd["X-averaged EC concentration overpotential [V]"].entries,
    color=Colors[0],linestyle=LS[0],label=r"Single-High $D_\times$") 

step_DD_LDx = DD_LDx_Crate['Sol_All'][Crate_index].cycles[0].steps[1]
axs[0].plot(
    step_DD_LDx['Time [s]'].entries-step_DD_LDx['Time [s]'].entries[0], 
    step_DD_LDx["X-averaged battery concentration overpotential [V]"].entries,
    color=Colors[3],linestyle=LS[3],label=r"Double-Low $D_\times$") 
axs[1].plot(
    step_DD_LDx['Time [s]'].entries-step_DD_LDx['Time [s]'].entries[0], 
    step_DD_LDx["X-averaged EC concentration overpotential [V]"].entries,
    color=Colors[3],linestyle=LS[3],label=r"Double-Low $D_\times$") 

step_DD_HDx = DD_HDx_Crate['Sol_All'][Crate_index].cycles[0].steps[1]
axs[0].plot(
    step_DD_HDx['Time [s]'].entries-step_DD_HDx['Time [s]'].entries[0], 
    step_DD_HDx["X-averaged battery concentration overpotential [V]"].entries,
    color=Colors[2],linestyle=LS[2],label=r"Double-High $D_\times$") 
axs[1].plot(
    step_DD_HDx['Time [s]'].entries-step_DD_HDx['Time [s]'].entries[0], 
    step_DD_HDx["X-averaged EC concentration overpotential [V]"].entries,
    color=Colors[2],linestyle=LS[2],label=r"Double-High $D_\times$") 

axs[0].set_ylabel("Potential [V]",fontsize=fs)
axs[1].set_xlabel("Time [s]",fontsize=fs)
axs[0].set_xlabel("Time [s]",fontsize=fs)
fig.suptitle(f"EC and Li+ overpotential - {Rate_Dis_All[Crate_index]}C Discharge", fontsize=fs+1)