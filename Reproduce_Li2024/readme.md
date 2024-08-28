## This folder contains scripts that can reproduce results in the paper: Li, Ruihe, Kirkaldy, Niall D., Oehler, Fabian, Marinescu, Monica, Offer, Gregory J., & O'Kane, Simon E. J. (2023). Lithium-ion battery degradation: using degradation mode analysis to validate lifetime prediction modelling. arXiv, [arXiv:2311.05482v2]. doi:10.48550/arXiv.2311.05482, which is currently under review and will be published soon. 

## Note that currently the base parameter set is called OKane2023, which doesn't associate to one specific parameter set, latter when tidied up, it should called Li2024. 

## This folder contains the following files:

1. Fun_NC.py contains necessary functions needed for all other files.

2. Run_long_Full_Exp1235.py contains scripts to run long simulation (ageing tests). To make it run, you will need to ensure you have generated the input files under folders "InputData/Full_Exp1235_NC". If you want to run other cases, just change this to other input files such as "Full_Exp23_Paper_11_fine", "SEI_Dry_Exp23_Paper_11_fine".

3. Get_Input_files.ipynb creates input files like "Full_Exp23_Paper_11_fine", "SEI_Dry_Exp23_Paper_11_fine" to be used in Run_long_Full_Exp1235.py to run long simulation.

4. Reload_Exp2.ipynb reloads simulation results for Experiment 2 (Fig. 2~5, Table 2 and 3 in main text, and Parametrization results section in SI)

5. Reload_Exp3.ipynb reloads simulation results for Experiment 3 (Fig. 6 in main text and Validation results section in SI)

6. 1C_GITT.ipynb and 2C_GITT.ipynb fit to 1C and 2C GITT results (Fig. S12 in the paper)

7. Compare_hysteresis.ipynb compares 0.1C charge/discharge voltage curves without and with current sigmoid. (Fig. S19 in SI)

8. Reload_Compare_Exp_2_3.ipynb reload simulation results to compare Experiment 2 and Experiment 3 (Fig. S17 and S18 in SI)

9. Reload_Exp1_5.ipynb reload simulation results on Experiment 1 and 5 (Fig. S16 in SI)

10. Full_Exp23_Paper_11_fine, SEI_Dry_Exp23_Paper_11_fine and SEI_Exp23_Paper_11_fine sweep over SEI growth activation energy (E_act^SEI) and the negative electrode diffusivity activation energy (E_act^(D_(s,n) )) and generate Fig. S10 in SI




