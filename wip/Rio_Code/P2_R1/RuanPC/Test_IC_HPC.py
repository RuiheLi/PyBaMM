import pybamm; import os;
from scipy.io import savemat,loadmat;
model = pybamm.lithium_ion.SPMe()
ChemistryChen=pybamm.parameter_sets.Chen2020_coupled   
Para_0=pybamm.ParameterValues(chemistry=ChemistryChen)
sim = pybamm.Simulation(model,parameter_values = Para_0,)
sim.solve([0, 360])

BasicPath=os.getcwd()
Target  = '/SaveData/'

if not os.path.exists(BasicPath + Target):
    os.mkdir(BasicPath + Target);


solution = sim.solution

t = solution["Time [s]"].entries
V = solution["Terminal voltage [V]"].entries

mdic = {
    "Time": t,
    "Terminal_voltage": V,
}

savemat(BasicPath + Target+'StructDara_for_Mat.mat',mdic)

print("I did it!")