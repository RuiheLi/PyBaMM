# This is ongoing work, still tidying up the script to make it more readable. 



# add 221205: if Timeout=='True', use Patrick's version, disable pool
#             else, use pool to accelerate 
# add Return_Sol, on HPC, always set to False, as it is useless, 
# add 230221: do sol_new['Throughput capacity [A.h]'].entries += sol_old['Throughput capacity [A.h]'].entries 
#             and for "Throughput energy [W.h]", when use Model.set_initial_conditions_from
#             this is done inside the two functions run_age_model_base_on_last_solution(_RPT)
# change 230621: add ability to scan one parameter set at different temperature 
#                and at Exp-2,3,5
# 240429: tidy up the script
# idea to tidy up: create a class: 
# 240723: tidy up, define class to take in everything


# 240828
idea: seperate files that are specific for GEM-2 and general for ParaSweeper