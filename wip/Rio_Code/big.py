# int - Dec=3e-10; single vs double 
Model_DFN_int_lowDec  = pybamm.lithium_ion.DFN(options={
    "SEI":"interstitial-diffusion limited",             
    "SEI film resistance":"distributed",          
    "SEI porosity change":"true",   
    "solvent diffusion": "single no consume wo refill",
    "electrolyte conductivity": "full"  ,}) 
Model_DD_int_lowDec  = pybamm.lithium_ion.DFN(options={
    "SEI":"interstitial-diffusion limited",             
    "SEI film resistance":"distributed",          
    "SEI porosity change":"true",   
    "solvent diffusion": "double no consume wo refill",
    "electrolyte conductivity": "sol full"  ,}) 
Model_All_int_lowDec =[ Model_DFN_int_lowDec ,   Model_DD_int_lowDec ]
Str_model_int_lowDec =[ 'Model_DFN_int_lowDec', 'Model_DD_int_lowDec' ]
# Para_All[0]: for single; Para_All[1]: for double
para=pybamm.ParameterValues("Li2023_ECdrag")
Para_All = []
for i in range(0,2):
    para=pybamm.ParameterValues("Li2023_ECdrag")
    para.update({"EC initial concentration in electrolyte [mol.m-3]":3500})
    para.update({"Typical EC concentration [mol.m-3]":3500})
    para.update({'EC diffusivity in SEI [m2.s-1]':5e-20})
    para.update({"EC diffusivity in electrolyte [m2.s-1]":EC_diffusivity_3E_10})
    para.update({'Inner SEI lithium interstitial diffusivity [m2.s-1]':3e-19})
    para.update({"EC Lithium ion cross diffusivity [m2.s-1]":Cross_diffusivity_1p5E_12})
    para.update({"Cation transference number":electrolyte_transference_number_EC_EMC_3_7_Landesfeind2019_Con})
    para.update({"Electrolyte conductivity [S.m-1]":electrolyte_conductivity_EC_EMC_3_7_Landesfeind2019_Con})
    para.update({"Electrolyte diffusivity [m2.s-1]":electrolyte_diffusivity_EC_EMC_3_7_Landesfeind2019_Con})
    Para_All.append(para)
Para_All[0].update({"Measured dLJP_dce":dLJP_One_Specie_dce_Jung2023})
# check 
print(Para_All[0]["Measured dLJP_dce"])
print(Para_All[1]["Measured dLJP_dce"])
Sol_all_int_lowDec =[];
for model,para in zip(Model_All_int_lowDec,Para_All):
    c_e = model.variables["Electrolyte concentration [mol.m-3]"]
    c_EC= model.variables["EC concentration [mol.m-3]"]
    T = model.variables["Cell temperature [K]"]
    D_e = para["Electrolyte diffusivity [m2.s-1]"]
    D_EC= para["EC diffusivity in electrolyte [m2.s-1]"]
    sigma_e = para["Electrolyte conductivity [S.m-1]"]
    Xi = para["EC transference number"]
    model.variables["Electrolyte diffusivity [m2.s-1]"] = D_e(c_e,c_EC, T)
    model.variables["EC diffusivity in electrolyte [m2.s-1]"] = D_EC(c_e,c_EC, T)
    model.variables["Electrolyte conductivity [S.m-1]"] = sigma_e(c_e,c_EC, T)
    model.variables["EC transference number"] = Xi(c_e,c_EC, T)
    model.variables["c(EC) over c(Li+)"] = c_EC / c_e
    t_0plus = para["Cation transference number"]
    model.variables["Cation transference number"] = t_0plus(c_e,c_EC, T)
    var_pts = {
        "x_n": 10,  # negative electrode
        "x_s": 5,  # separator 
        "x_p": 10,  # positive electrode
        "r_n": 30,  # negative particle
        "r_p": 20,  # positive particle
    }
    sim    = pybamm.Simulation(
        model, experiment = exp,
        parameter_values = para,
        solver = pybamm.CasadiSolver(return_solution_if_failed_early=True),
        var_pts=var_pts,)       
    Sol_all_int_lowDec.append(sim.solve())