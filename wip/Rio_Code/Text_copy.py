V_RMSE_All = {}
for i,case in enumerate(Cases):
    V_RMSE_All[Str_cases[i]] = {}
    for T_deg in T_deg_All:
        V_RMSE_All[Str_cases[i]] [T_deg] = []

for i,case in enumerate(Cases):
    for T_deg in T_deg_All:
        for i_RPT in range(13):
            V_RMSE_All[Str_cases[i]][T_deg].append(Calculate_0P1C_V_RMSE(
                index_exp,Exp_2_AllData,
                Cases[i],T_deg,i_RPT)  )
i=0; T_deg=10
#print(f"{Str_cases[i]}, T={T_deg}degC")
#plt.plot( V_RMSE_All[Str_cases[i]][T_deg])
fig, axs = plt.subplots(1,3, figsize=(16,4),gridspec_kw={'top': 0.9, 'bottom': 0.15})
for i_ax,T_deg in enumerate(T_deg_All):
    for j in range(len(Str_cases)):
        axs[i_ax].plot(V_RMSE_All[Str_cases[j]][T_deg],  
            '-', color=Default_Colors_Alpha[j],  
            #linestyle=LS[0],
            label=Str_cases[j]) 
        axs[i_ax].set_xlabel(r"RPT No.") 
        axs[i_ax].set_ylim([0.2,1.3])
        axs[i_ax].set_xticks([0,3,6,9,11])
        axs[i_ax].set_title(f"T={T_deg}"+r"$^\circ$C")
    if i_ax == 0:
        axs[i_ax].legend(prop={'family':'DejaVu Sans','size':fs},loc='best',frameon=False)
        axs[i_ax].set_ylabel("RMSE (%)")
plt.savefig(target_folder + f"/RMSE All 3Ts.png", dpi=800)
plt.savefig(target_folder + f"/RMSE All 3Ts.svg")
for i_ax,T_deg in enumerate(T_deg_All):
    for j in range(len(Str_cases)):
        print(f"T={T_deg},{Str_cases[j]}:{np.mean(V_RMSE_All[Str_cases[j]][T_deg]):.2f}")