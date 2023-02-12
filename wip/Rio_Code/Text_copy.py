SoH = Cap_1CellAllRPTs_df['Cell Capacity']/Cap_1CellAllRPTs_df['Cell Capacity'][0]
LAM_pe = 1 - (Cap_1CellAllRPTs_df['PE Capacity']/Cap_1CellAllRPTs_df['PE Capacity'][0])
LAM_ne_tot = 1 - (Cap_1CellAllRPTs_df['NE(tot) Capacity']/Cap_1CellAllRPTs_df['NE(tot) Capacity'][0])
LAM_ne_Gr = 1 - (Cap_1CellAllRPTs_df['NE(Gr) Capacity']/Cap_1CellAllRPTs_df['NE(Gr) Capacity'][0])
LAM_ne_Si = 1 - (Cap_1CellAllRPTs_df['NE(Si) Capacity']/Cap_1CellAllRPTs_df['NE(Si) Capacity'][0])
LLI = (Cap_1CellAllRPTs_df['PE Capacity'][0] - Cap_1CellAllRPTs_df['PE Capacity'] - (Cap_1CellAllRPTs_df['Offset'][0]-Cap_1CellAllRPTs_df['Offset']))/Cap_1CellAllRPTs_df['Cell Capacity'][0]

# Compile the DM parameters into a dataframe
DM_df = pd.DataFrame(
    data={
    'SoH':SoH, 'LAM PE':LAM_pe, 'LAM NE_tot':LAM_ne_tot, 
    'LAM NE_Gr':LAM_ne_Gr, 'LAM NE_Si':LAM_ne_Si, 'LLI':LLI})