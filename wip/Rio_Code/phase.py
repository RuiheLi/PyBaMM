Targets = [
    "Andrew_int=6e-19_2000_chi","Valoen_int=6e-19_2000_chi","Constant_int=6e-19_2000_chi",
    "Andrew_sol=3e-20_2000_chi","Valoen_sol=3e-20_2000_chi","Constant_sol=3e-20_2000_chi",
    ]


my_dict = HPC_Age_230105[target][str(26)]
keys_time = "CD Time [h]"
time_1 = my_dict[keys_time][0][10][0]
vol_1  = my_dict["CD Terminal voltage [V]"][0][10][0]
plt.plot(time_1,vol_1)