    axs[0,0].plot(
        All_Scans[str(scan)]['Throughput capacity [kA.h]'], 
        All_Scans[str(scan)]['CDend SOH [%]'],     
        '-o', label="Scan=" + str(scan) )
    axs[0,1].plot(
        All_Scans[str(scan)]['Throughput capacity [kA.h]'], 
        All_Scans[str(scan)]["CDend LLI [%]"],'-o', label="total LLI")
    if model_options.__contains__("lithium plating"):
        axs[0,1].plot(
            All_Scans[str(scan)]['Throughput capacity [kA.h]'], 
            All_Scans[str(scan)]["CDend LLI lithium plating [%]"],'--o', label="LiP")
    if model_options.__contains__("SEI"):
        axs[0,1].plot(
            All_Scans[str(scan)]['Throughput capacity [kA.h]'], 
            All_Scans[str(scan)]["CDend LLI SEI [%]"] ,'--o', label="SEI")
    if model_options.__contains__("SEI on cracks"):
        axs[0,1].plot(
            All_Scans[str(scan)]['Throughput capacity [kA.h]'], 
            All_Scans[str(scan)]["CDend LLI SEI on cracks [%]"] ,'--o', label="SEI-on-cracks")
    axs[0,2].plot(
        All_Scans[str(scan)]["Throughput capacity [kA.h]"], 
        All_Scans[str(scan)]["CDend LAM_ne [%]"],     '-o', ) 
    axs[1,0].plot(
        All_Scans[str(scan)]["Throughput capacity [kA.h]"], 
        All_Scans[str(scan)]["CDend LAM_pe [%]"],     '-o',  ) 
    axs[1,1].plot(
        All_Scans[str(scan)]["Throughput capacity [kA.h]"], 
        np.array(All_Scans[str(scan)]["Res_0p5C_50SOC"]),     '-o', ) 
    axs[1,2].plot(
        All_Scans[str(scan)]["Throughput capacity [kA.h]"][1:], 
        np.array(All_Scans[str(scan)]["avg_Age_T"][1:]),     '-o', ) 