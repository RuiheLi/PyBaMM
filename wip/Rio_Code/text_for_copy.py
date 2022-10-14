
c_e = Model.variables["Electrolyte concentration [mol.m-3]"]
T = Model.variables["Cell temperature [K]"]
c_EC = Model.variables["EC concentration [mol.m-3]"]
Model.variables["c(EC) over c(Li+)"] = c_EC / c_e
Model.variables["Electrolyte conductivity [S.m-1]"] =(
    Para_0['Electrolyte conductivity [S.m-1]'](c_e,c_EC, T))
Model.variables["Electrolyte diffusivity [m2.s-1]"] =(
    Para_0['Electrolyte diffusivity [m2.s-1]'](c_e,c_EC, T))