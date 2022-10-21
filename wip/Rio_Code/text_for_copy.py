model_2 = pybamm.lithium_ion.DFN()
c_e = model_2.variables["Electrolyte concentration [mol.m-3]"]
T = model_2.variables["Cell temperature [K]"]
model_2.variables["Electrolyte conductivity [S.m-1]"] =(
    Para_0['Electrolyte conductivity [S.m-1]'](c_e, T))
model_2.variables["Electrolyte diffusivity [m2.s-1]"] =(
    Para_0['Electrolyte diffusivity [m2.s-1]'](c_e, T))

sim_2 = pybamm.sim_2ulation(
    model_2, experiment = Experiment_Long,
    parameter_values = Para_0,
    solver = pybamm.CasadiSolver(),
    var_pts=var_pts,)  
try:
    sol_2 = sim_2.solve(save_at_cycles=save_at_cycles,);
    print(sol_2.cycles[-1].steps[-1]);  # a way to check whether the solution is finalized 
except:
    print('Fail for electrolyte: Exp')
else:
    Sol.append(sol_2)   
    print('Succeed for electrolyte: Exp')