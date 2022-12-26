    axs[1,1].plot(
        Full_cycle, 
        my_dict_AGE["CC LLI to SEI per Ah"],
        my_dict_AGE["CC LLI to SEI per Ah"],
        '--^',label="CC" )
    axs[1,1].plot(
        Full_cycle, 
        my_dict_AGE["CV LLI to SEI per Ah"],
        my_dict_AGE["CV LLI to SEI per Ah"],
        '--^',label="CV" )
    axs[1,1].plot(
        Full_cycle, 
        my_dict_AGE["Cha LLI to SEI per Ah"],
        my_dict_AGE["Cha LLI to SEI per Ah"],
        '-s',label="Charge" )
    axs[1,1].plot(
        Full_cycle, 
        my_dict_AGE["Dis LLI to SEI per Ah"],
        my_dict_AGE["Dis LLI to SEI per Ah"],
        '-o',label="Discharge" )