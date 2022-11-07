axs[2,0].plot(
    Full_cycle, 
    my_dict_AGE["CDend Total EC in electrolyte [mol]"] ,
    '-o',label="In elely" )
axs[2,0].plot(
    Full_cycle, 
    my_dict_AGE["CDend Total EC in electrolyte and SEI [mol]"] ,
    '-^',label="In elely and SEI"  )
axs[2,1].plot(
    Full_cycle, 
    my_dict_AGE["CDend Total lithium in electrolyte [mol]"],
    '-o', )

axs[2,2].plot(
    Full_cycle, 
    my_dict_AGE["CDend Total lithium in particles [mol]"],
    '-o',)

axs[2,0].set_ylabel("Quantity [mol]",   fontdict={'family':'DejaVu Sans','size':fs})
axs[2,0].set_title("Total EC",   fontdict={'family':'DejaVu Sans','size':fs+1})
axs[2,1].set_ylabel("Quantity [mol]",   fontdict={'family':'DejaVu Sans','size':fs})
axs[2,1].set_title("Total lithium in electrolyte",   fontdict={'family':'DejaVu Sans','size':fs+1})
axs[2,2].set_ylabel("Quantity [mol]",   fontdict={'family':'DejaVu Sans','size':fs})
axs[2,2].set_title("Total lithium in particles",   fontdict={'family':'DejaVu Sans','size':fs+1})
