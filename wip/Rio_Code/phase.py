N=100
sto = np.linspace(0, 1, N); T_0 = 298.15; lw=2
aa=np.array(electrolyte_diffusivity_EC_DMC_1_1_Landesfeind2019_Con(sto,4000,T_0))
fig, axs = plt.subplots(1,2, figsize=(18/2.54,6/2.54),tight_layout=True) # 

axs[0].plot(
    sto/1e3, 
    graphite_LGM50_diffusivity_ORegan2022(sto,4000,T_0), 
    lw=lw, ) 
axs[1].plot(
    sto/1e3, 
    nmc_LGM50_diffusivity_ORegan2022(sto,4000,T_0), 
    lw=lw, ) 
for i in range(3):
    axs[i].set_xlabel(r"c(Li$^+$) [M]")

axs[0].set_ylabel("$\kappa_\mathrm{e}$ [S/m]")
axs[1].set_ylabel("D$_\mathrm{e}$ [m$^\mathrm{2}$/s]") # 
axs[2].set_ylabel("$\mathit{t}_\mathrm{0}^\mathrm{+}$")
plt.savefig(BasicPath +  Target+ 
    f"Landesfeind2019 EC_DMC_1_1.png", dpi=1000)
plt.savefig(BasicPath +  Target+ 
    f"Landesfeind2019 EC_DMC_1_1.svg")  