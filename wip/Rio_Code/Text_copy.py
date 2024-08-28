LLI_to_LAM_Exp3_5couple.append(
        np.interp(
        x_0, 
        Full_Exp3_GoodFit[str(T_deg)][0][x_key],
        Full_Exp3_GoodFit[str(T_deg)][0]['CDend LLI [%]'])
        - 
        np.interp(
        x_0, 
        Full_Exp3_GoodFit[str(T_deg)][0][x_key],
        Full_Exp3_GoodFit[str(T_deg)][0]['CDend LLI SEI [%]'])
        - 
        np.interp(
        x_0, 
        Full_Exp3_GoodFit[str(T_deg)][0][x_key],
        Full_Exp3_GoodFit[str(T_deg)][0]['CDend LLI SEI on cracks [%]'])
        - 
        np.interp(
        x_0, 
        Full_Exp3_GoodFit[str(T_deg)][0][x_key],
        Full_Exp3_GoodFit[str(T_deg)][0]['CDend LLI lithium plating [%]'])
        )