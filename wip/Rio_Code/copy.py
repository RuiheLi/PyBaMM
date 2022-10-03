label = ["0.3 d","0.3 s","3e-4 d","3e-4 s",] 
output_variables3 = [
    "Terminal voltage [V]",   
    "Current [A]",
    "EC concentration [mol.m-3]",
    "Electrolyte concentration [mol.m-3]",
    #"Loss of capacity to SEI [A.h]",
    "Electrolyte flux [mol.m-2.s-1]",
    "EC flux [mol.m-2.s-1]",
    #"Porosity times EC concentration",
    #"Porosity times concentration",  0p3 3e_4_
    "Cation transference number 3e_4",
    "Cation transference number 0p3",
]
quick_plot = pybamm.QuickPlot([
    Sol_0p3_ddiff,Sol_0p3_sdiff,Sol_3e_4_ddiff,Sol_3e_4_sdiff
    ], output_variables3,label,
    variable_limits='fixed',time_unit='hours',n_rows=2,
    ) #     spatial_unit='mm',figsize = (330,140)
quick_plot.dynamic_plot();
BasicPath = 'D:/OneDrive - Imperial College London/SimDataSave/P3R4/'; 
Target  = 'a0_case3_doublediffusion_noSEI/'
if not os.path.exists(BasicPath + Target):
    os.mkdir(BasicPath + Target);
output_filename = BasicPath + Target + '/t_0+_func.gif'

quick_plot.create_gif(
    number_of_images=10, duration=1,output_filename=output_filename)