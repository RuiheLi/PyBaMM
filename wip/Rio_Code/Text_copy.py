df_10degC = pd.read_excel(BasicPath+"Summary of 200 cases-10oC.xlsx")
# Display the DataFrame
# print(df_10degC.head())
column_mapping = {
    "Inner SEI lithium interstitial diffusivity [m2.s-1]":"Dint",
    "Dead lithium decay constant [s-1]":"Decay",
    "Lithium plating kinetic rate constant [m.s-1]"	:"k_LiP",
    "Negative electrode LAM constant proportional term [s-1]"	:"pLAM_Ne",
    "Negative electrode cracking rate":"k_Ne_cr",
    "Outer SEI partial molar volume [m3.mol-1]":"V_SEI",	
}
# Rename the columns using the mapping   - Error_1~6;   Para_1~8
df_10degC = df_10degC.rename(columns=column_mapping)
df_10degC_s = df_10degC[df_10degC["Error Tot%"] != "Unknown"].copy()
df_10degC_f = df_10degC[df_10degC["Error Tot%"] == "Unknown"].copy()
df_10degC_f = df_10degC_f[~df_10degC_f["Scan No"].isin(df_10degC_s["Scan No"])].copy()