def Plot_All_8_Para(D_0_e_EC,Save_fig):
    """ Function to look: 
    # electrolyte_conductivity_Wang2021(c_e,c_EC)
    # Fun_t_0plus_Wang2021(c_e,c_EC)
    # dLJP_2_Species_dc_e_np(c_e,c_EC)
    # dLJP_2_Species_dc_EC_np(c_e,c_EC)
    # dLJP_1_Specie_dc_e_np(c_e)
    # Fun_Xi_tidle(c_e,c_EC)
    # diff_constant_3E_10(c_e, c_EC )
    # EC_diffusivity_5E_10(c_e, c_EC)
    # Fun_D_0_EC_e """
    fig, axss = plt.subplots(
        3,3, figsize=(18/2.54,14/2.54),tight_layout=True)
    axs   = axss.flatten()
    # get c range:  # When c_e=1000 mol/m3, c_EC=5622.86, c_T=13245.72 mol/m3
    c_EC_1_specie = 6209.49; C_EC = np.linspace(1E3,1E4,300).tolist();   
    C_e = np.linspace(100,3.5E3,100) # .tolist()
    # first: dLJP_2_Species_dc_e_np and dLJP_1_Specie_dc_e_np
    axs[0].plot(C_e/1e3,dLJP_2_Species_dc_e_np(C_e,c_EC_1_specie),"-",label="Double")
    axs[0].plot(C_e/1e3,dLJP_1_Specie_dc_e_np(C_e),"--",label="Single")
    axs[0].set_ylabel(r'$\partial \Delta U / \partial c_\mathrm{e}$ / $\mathrm{V/(mol/m^3)}$')
    # second: dLJP_2_Species_dc_EC_np
    axs[1].plot(C_e/1e3,dLJP_2_Species_dc_EC_np(C_e,c_EC_1_specie),"-",label="Double")
    axs[1].set_ylabel(r'$\partial \Delta U / \partial c_\mathrm{EC}$ / $\mathrm{V/(mol/m^3)}$')
    for i in range(2):
        axs[i].ticklabel_format( 
            axis='y', style='sci',scilimits=[-0.01,0.01], 
            useOffset=None, useLocale=None, useMathText=None)
    # 3rd: electrolyte_conductivity_Wang2021
    axs[2].plot(C_e/1e3, electrolyte_conductivity_Wang2021(C_e,c_EC_1_specie), ) 
    axs[2].set_ylabel("$\kappa_\mathrm{e}$ / S/m")
    # 4th: Fun_t_0plus_Wang2021
    axs[3].plot(C_e/1e3, Fun_t_0plus_Wang2021(C_e,c_EC_1_specie), ) 
    axs[3].set_ylabel(r'$t_{+}^{0}$')
    # 5rd: Fun_Xi_tidle
    axs[4].plot(C_e/1e3, Fun_Xi_tidle(C_e,c_EC_1_specie), ) 
    axs[4].set_ylabel(r"$\widetilde{\Xi}$") # $$
    # 6-9: D_ee, D_EC,EC, D_e,EC, D_EC,e 
    axs[5].plot(C_e/1e3, np.ones(np.size(C_e))*3e-10,"-")
    axs[5].set_ylabel(r"$D_\mathrm{ee}^\mathrm{0}$ / $\mathrm{m}^2/\mathrm{s}$")
    axs[6].plot(C_e/1e3, np.ones(np.size(C_e))*5e-10,"-")
    axs[6].set_ylabel(r"$D_\mathrm{EC,EC}^\mathrm{0}$ / $\mathrm{m}^2/\mathrm{s}$")
    D_0_e_EC = 1.5e-11
    axs[7].plot(C_e/1e3, np.ones(np.size(C_e))*D_0_e_EC,"-")
    axs[7].set_ylabel(r"$D_\mathrm{e,EC}^\mathrm{0}$ / $\mathrm{m}^2/\mathrm{s}$")
    axs[8].plot(C_e/1e3, Fun_D_0_EC_e(C_e,c_EC_1_specie,D_0_e_EC),"-")
    axs[8].set_ylabel(r"$D_\mathrm{EC,e}^\mathrm{0}$ / $\mathrm{m}^2/\mathrm{s}$")


    return 
Plot_All_8_Para(0,Save_fig)