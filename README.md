# Binary Activity Thermodynamics Model
Water plays an essential role in environmental chemistry, but when thousands of organic species are involved highly-detailed thermodynamic modeling can become problematic if not altogether impractical, due to limited species information.
The BAT model was built to correct for this limitation as it is a water-sensitive organic species thermodynamic model that can adapt to different levels of chemical information. 
The BAT modeling framework describes a mixture of organics with water as binary mixtures of one organic and associated water. 

For gas-particle partitioning, we use a non-ideal Volatility Basis Set (VBS) formulation that includes activity coefficients and accounts for liquid-liquid equilibrium. 
Our model is capable of using both detailed molecular structures as well as bulk properties like O:C, molar mass, and component volatility classes.

If you use the model please cite:
Gorkowski, K., Preston, T. C. and Zuend, A.: RH-dependent organic aerosol thermodynamics via an efficient reduced-complexity model, Atmos. Chem. Phys. Discuss., (June), 1–37, https://doi.org/10.5194/acp-2019-495, 2019.

## Code and File Overview
The source code was written in MATLAB 2018b, though it will likely run with past and future versions. 
The source code is in [Matlab_source_code](https://github.com/Gorkowski/Binary_Activity_Thermodynamics_Model/tree/master/Matlab_source_code), which contains [update_Matlab_paths.m](https://github.com/Gorkowski/Binary_Activity_Thermodynamics_Model/blob/master/Matlab_source_code/update_Matlab_paths.m) function that should be run to add the source folder to your MATLAB path.

In addition to the source code, a standalone executable file has been compiled and can run the model if you have no MATLAB license.

## Standalone Executable Steps
The Windows installation file is [BAT_installer_v1_web.exe](https://github.com/Gorkowski/Binary_Activity_Thermodynamics_Model/tree/master/Matlab_runtime/BAT_Model/for_redistribution), which needs the freely available MATLAB Runtime R2018b engine. 
The installation file will download the respective runtime engine. 
You can also install the runtime engine independently from [MathWorks]( https://www.mathworks.com/products/compiler/matlab-runtime.html).

Then run the installed program, and a file dialog box will appear. Then select one of the example input files, e.g., [simple_input_single_run.txt](https://github.com/Gorkowski/Binary_Activity_Thermodynamics_Model/blob/master/Input_examples/simple_input_single_run.txt).

### Input File for Simulations

The executable is expecting is a comma delimited text file, with one or more simulation input blocks. 
Each input block starts with the `run name`, which is used when writing the output file.
>run name=aPsoa_t1 

The next line contains the refinement method for the BAT model, which is different from the VBS model. The `BAT_refinement_mode` can only be interpolate, different methods can be added. Note, this is only used when the water activity is above `VBSBAT_options.BAT_refinement_aw=0.9` and the error in water activity is above `VBSBAT_options.BAT_refinement_tolerance=1e-06`
>BAT_refinement_mode=interpolate

Next, are the VBS+BAT options, which have predefined run modes. You have to start with a predefined run mode, and then declare any changes in subsequent lines. The run mode options are `default` (used in [Gorkowski et al. (2019)](https://doi.org/10.5194/acp-2019-495)), `robust` which takes longer and provides a check on how well the neural networks are doing, `NN only` uses only the neural networks without refinement, and `beta only` which forces a single organic-rich phase.
>VBSBAT_options.run_mode_used=default ,{default,robust,NN only,beta only} 

Below this, you can add your modifications, like turning off the graphs.
>VBSBAT_options.plot_PM=no ,{yes, no : option to make simple output graph}

The complete option list and settings used will be generated after each simulation run in a separate text file, like `VBSBAT_input_used_aPsoat1.txt`. 

Then set the 'dry' water activity used to calculate the Csat_j, if only C_liquid and Cstar are given.
>calculate_Csat_j_with_aw=0

Now `water activity` is used as a trigger word to indicate the next line is a list of water activities used to run this mixture and that the `VBSBAT_options` have ended. The `water activity` values can have non-uniform spacing, but should be monatomic. 
>water activity
>0.9999,0.99,0.985,0.98,0.97,0.96,0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.15,0.050.01,0.001,

The next trigger word is `system properties`, which is used to indicate the reading in of the organic composition. Note, currently the `BAT functional group` **must have a space** after and before the separating comma. The `optional Name` **IS NOT OPTIONAL** and does nothing currently, but it will be used in the future for more detailed output files. We are working on improving the input file reader.
>system properties,
>M (g/mol), O:C, H:C  , eff. Csat_j (ug/m3), Ctotal_j (ug/m3) , BAT functional group, optional Name, optional Cstar (ug/m3), optional C^liquid (ug/m3),

>2.00E+02,4.00E-01,1.60E+00,5.74E+03,8.79E+00, hydroperoxideSOA ,C107OOH,8620.171693,0.007057093,

>1.88E+02,4.44E-01,1.78E+00,0.00E+00,3.98E+00, hydroperoxideSOA ,C97OOH,522.7659518,0.052085067,

>2.16E+02,5.00E-01,1.60E+00,0.00E+00,1.13E+00, hydroperoxideSOA ,C108OOH,231.757194,0.032911505,

An empty line is used to separate simulation input blocks, and next simulation input block is triggered by a line containing `run name`


## Data and Figures From Paper
The figures presented in [Gorkowski et al. (2019)](https://doi.org/10.5194/acp-2019-495) can be reproduced using the source code in [paper_figures](https://github.com/Gorkowski/Binary_Activity_Thermodynamics_Model/tree/master/Matlab_source_code/paper_figures). The data shown in the paper and the fit parameters are compiled in the [Excel workbook file](https://github.com/Gorkowski/Binary_Activity_Thermodynamics_Model/tree/master/Matlab_source_code/paper_figures/Figure_data).
