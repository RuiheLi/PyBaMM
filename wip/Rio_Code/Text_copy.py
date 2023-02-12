def Wrap_OCV_Fit(Cell_i_Dis_1RPT, str_cell, RPT_num, Num_try, Input_pack):  
    
    # Unpack input parameters:
    [guess_values_0, Min_guess_0, Max_guess_0,
        NE_Gr_OCV, NE_Si_OCV, PE_OCV,Path, Save] = Input_pack
    fs=17; lw=2.5; ms=6

    Guess_Values = Guess_0_Random(Num_try,guess_values_0,Min_guess_0,Max_guess_0)

    Cell_i_Dis_1RPT = Cell_i_Dis_1RPT[
        Cell_i_Dis_1RPT[
        'Current (mA)']<0].loc[
        :, ['Charge (mA.h)', 'Voltage (V)']]
    Cell_i_Dis_1RPT.reset_index(inplace=True, drop=True)
    Cell_i_Dis_1RPT['SOC (%)'] =(
        1 - 
        Cell_i_Dis_1RPT['Charge (mA.h)']
        /Cell_i_Dis_1RPT['Charge (mA.h)'].max()   )
    z_1Cell1RPT_Alltries = []; Err_1Cell1RPT_Alltries=[]; Fit_1Cell1RPT_Alltries = [];
    for guess_values in Guess_Values:
        z_out, _, _ = dmc.stoich_OCV_fit_multi_comp(
            NE_Gr_OCV, NE_Si_OCV, PE_OCV, Cell_i_Dis_1RPT,
            z_guess=guess_values)
        Fit_1Cell1RPT, _, _, _ = dmc.calc_full_cell_OCV_multi_standalone(
            NE_Gr_OCV, NE_Si_OCV, PE_OCV, Cell_i_Dis_1RPT['SOC (%)'], 
            *z_out)
        fitted_Gr_fract = z_out[4].round(4)
        fitted_Si_fract = 1-fitted_Gr_fract
        err_BoL = dmc.DM_error_check(
            NE_Gr_OCV, NE_Si_OCV, PE_OCV, Cell_i_Dis_1RPT, 
            z_out)
        # print(f"RMSE for {guess_values} is {(err_BoL*1e3).round(2)} mV")
        # Append
        z_1Cell1RPT_Alltries.append(z_out); Err_1Cell1RPT_Alltries.append(err_BoL); 
        Fit_1Cell1RPT_Alltries.append(Fit_1Cell1RPT)
    Err_1Cell1RPT_Alltries = np.array(Err_1Cell1RPT_Alltries)    
    Min_index = Find_min(Err_1Cell1RPT_Alltries,8)
    # History: for 20 rounds of optimization, 64s is needed. 
    import matplotlib as mpl; import copy; 
    font = {'family' : 'DejaVu Sans','size': fs}
    mpl.rc('font', **font)
    Index_Plot = {
        "0":[0,1],"1":[0,2],"2":[1,0],"3":[1,1],
        "4":[1,2],"5":[2,0],"6":[2,1],"7":[2,2] }

    # Start a plot
    Fig_DMA_1cell1RPT, ax = plt.subplots(3,3,  figsize=(16,12))
    Err_1Cell1RPT_Alltries_temp = copy.deepcopy(Err_1Cell1RPT_Alltries)
    Err_1Cell1RPT_Alltries_temp = abs( np.sort(-Err_1Cell1RPT_Alltries_temp)); 
    Temp_List = np.arange(1,len(Err_1Cell1RPT_Alltries_temp)+1)
    ax[0,0].plot(Temp_List,  Err_1Cell1RPT_Alltries_temp*1e3,  
        color='grey',linewidth=lw,  linestyle='-',marker = 'o',markersize = ms)
    ax[0,0].plot(Temp_List[-8:] ,  Err_1Cell1RPT_Alltries_temp[-8:]*1e3,  
        color='b',linewidth=lw,  linestyle='none',marker = 'o',markersize = ms+2)
    ax[0,0].set_ylabel('RMSE (mV)')
    ax[0,0].set_xlabel('Try number (sort)')
    # plot things:
    for min_index,i in zip(Min_index,range(len(Min_index))):
        SOC_OCV_Temp = Fit_1Cell1RPT_Alltries[min_index].to_numpy()
        SOC_OCV_Exp  = Cell_i_Dis_1RPT.to_numpy()
        m,n = Index_Plot[str(i)]
        ax[m,n].plot(SOC_OCV_Exp[:,2]*100,SOC_OCV_Exp[:,1],linewidth=lw, ls='-',color='k',label='Exp', ) 
        ax[m,n].plot(SOC_OCV_Temp[:,0]*100,SOC_OCV_Temp[:,1],linewidth=lw, ls='--',color='b',label='Fit', ) 
        ax[m,n].set_xlabel('SOC %')
        ax[m,n].set_ylabel('Voltage (V)')
        ax[m,n].text(35,2.7,
            f'Try Num: {min_index}\nRMSE: {format(Err_1Cell1RPT_Alltries[min_index]*1e3,".2f")} mV',
            horizontalalignment='left',
            verticalalignment='bottom', fontsize=fs)
    ax[0,1].legend(frameon=False,)

    Fig_DMA_1cell1RPT.suptitle(
        f'Fit result of Cell {str_cell}, RPT-{RPT_num}', 
        fontsize=fs+4)
    Fig_DMA_1cell1RPT.tight_layout()
    if Save == True:
        Fig_DMA_1cell1RPT.savegif(
            BasicPath + Target+ str(Scan_i)+"/Cracks related.png", dpi=dpi)
            plt.savefig(BasicPath + Target+ str(Scan_i)+"/Cracks related.png", dpi=dpi)
        

    return z_1Cell1RPT_Alltries,Err_1Cell1RPT_Alltries,Fit_1Cell1RPT_Alltries,Min_index