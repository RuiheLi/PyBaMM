Para_wEC=pybamm.ParameterValues(chemistry=ChemistryChen);
Para_wEC['EC transference number'] =    Xi
Para_wEC['Cation transference number'] =     t_0plus
Para_wEC['EC Lithium ion cross diffusivity [m2.s-1]'] = Dcross
Para_wEC['Typical EC Lithium ion cross diffusivity [m2.s-1]'] =  Dcross
Para_wEC['Electrolyte diffusivity [m2.s-1]'] =  De
Para_wEC['EC diffusivity in electrolyte [m2.s-1]'] =  Dec
Para_wEC['Ratio of lithium moles to SEI moles'] =  1
Para_wEC.update({
    'Electrolyte conductivity [S.m-1]':
    electrolyte_conductivity_Nyman2008Exp_wEC}, )