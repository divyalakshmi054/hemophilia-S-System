# hemophilia-S-System
# BST Model for Hemophilia and Treatment

To run simulations: Keep visit number = 4 in `sample_ensemble.jl`, then `include("sample_ensemble.jl")`

To estimate parameters by reducing differences: Specify FVIIa and FII levels in the learn_routine() function in `estimate_model_parameters.jl`, then change p_best dump and SIM data dump filenames to include info about FVIIa added and FII percentage. Finally, run `include("estimate_model_parameters.jl")`

To use estimated parameters to generate other cases: Pick a PSET, and change FVIIa and FII levels in `construction.jl`. Then,  `include("construction.jl")`. This file has a plotting routine to generate FIIa v. Time *for each patient* and store to `figs` in pwd.

To generate average plots for N patients: Open the appropriate experimental data file and in a loop running from 1 to N, open simulation file(s) in `plot_thrombin.jl` to compute mean and standard error. Then, `include("plot_thrombin.jl")`

To run Morris sensitivity analysis: Edit number of samples in gsa method in `sensitivity.jl`, then `include("sensitivity.jl")`

To plot Sobol indices: `include("plot_sensitivity.jl")`; you may need to update the CSV filename or the path to file containing the Morris results
