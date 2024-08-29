# ParaSweeper is a script developed by the ESE group at Imperial College London that aims at sweeping over variables (parameters) of a complicated physics-based models. The obtained input-output pairs can then be used to do sensitivity analysis, train suggorate models, etc. The main features/advantages of this tool lays in that it handles inputs/outputs and errors properly. It facilitates high-throughput computing for those "heavy", computational expensive models. 

## In its examples, we showcase how to use this model to run solvent-consumption model, and reproduce the Li2024 paper (GEM-2 paper).

## Define the terminology:
- sweep_task: a series of cases to be run, the different cases are defined slightly in a way that it sweep over one, two or more "variables". 

- case: running a case can be regarded as running some experiments (can be pybamm.experiments) on a pre-defined virtual cells (defined by pybamm.parameter_values) with pre-defined physics (pybamm.models), and get outputs (pybamm.solution)

- config_exp: configuration to define pybamm.experiments

- config_para: configuration to define pybamm.parameter_values

- config_model: configuration to define pybamm.model

- config_sol: container to store customized results of pybamm solution

- config_data: configuration to define experimentally measured data

- solvent_consumption_setting: 

- gem_setting:



