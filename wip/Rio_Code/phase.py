sol_dd = SD_Dis_All[0]
sol_dd = DD_Dis_All[0]
t_seconds = sol_dd["Time [s]"].entries
t_hours = (t_seconds - 60) / 3600
I = sol_dd["Current [A]"].entries
Q = sol_dd["Discharge capacity [A.h]"].entries
V = sol_dd["Terminal voltage [V]"].entries

# No zoom in:
fig, ax = plt.subplots()
ax.plot(t_exp-t_exp[0],V_exp,color='k',linewidth=2,  linestyle='-',label='Experiment')
ax.plot(t_hours,V,color='r',linewidth=2,  linestyle='--',label='Simulation')

ax.set_xlabel('Time [h]')
ax.set_ylabel('Terminal voltage [V]')
ax.set_title('25x1C pulses at 298K, comparison')
ax.legend(loc='best',frameon=False)
fig.suptitle(f'Scan = {index_i}')
plt.show()