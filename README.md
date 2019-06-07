# Binary Activity Thermodynamics Model
Water plays an essential role in environmental chemistry, but when thousands of organic species are involved highly-detailed thermodynamic modeling can become problematic if not altogether impractical, due to limited species information.
The BAT model was built to correct for this limitation as it is a water-sensitive organic species thermodynamic model that can adapt to different levels of chemical information. 
The BAT modeling framework describes a mixture of organics with water as binary mixtures of one organic and associated water. 

For gas-particle partitioning, we use a non-ideal Volatility Basis Set (VBS) formulation that includes activity coefficients and accounts for liquid-liquid equilibrium. 
Our model is capable of using both detailed molecular structures as well as bulk properties like O:C, molar mass, and component volatility classes.

If you use the model please cite:
Gorkowski, K., Preston, T. C. and Zuend, A.: RH-dependent organic aerosol thermodynamics via an efficient reduced-complexity model, Atmos. Chem. Phys. Discuss., (June), 1–37, https://doi.org/10.5194/acp-2019-495, 2019.

## Code and File Overview
The source code is written in MATLAB 2018b, though it will likely run with past and future versions. 
The source code is in [Matlab_source_code](https://github.com/Gorkowski/Binary_Activity_Thermodynamics_Model/tree/master/Matlab_source_code), which contains [update_Matlab_paths.m](https://github.com/Gorkowski/Binary_Activity_Thermodynamics_Model/blob/master/Matlab_source_code/update_Matlab_paths.m) function that should be run to add the source folder to your MATLAB path.

In addition to the source code, a standalone executable file has been compiled and can run the model if you have no MATLAB license.

## Standalone Executable Steps
The Windows installation file is [BAT_installer_v1_web.exe](https://github.com/Gorkowski/Binary_Activity_Thermodynamics_Model/tree/master/Matlab_runtime/BAT_Model/for_redistribution), which needs the freely available MATLAB Runtime R2018b engine. 
The installation file will download the respective runtime engine. 
You can also install the runtime engine independently from [MathWorks]( https://www.mathworks.com/products/compiler/matlab-runtime.html).

Then run the installed program, and a file dialog box will appear. Then select one of the example input files, e.g., [simple_input_single_run.txt](https://github.com/Gorkowski/Binary_Activity_Thermodynamics_Model/blob/master/Input_examples/simple_input_single_run.txt).

## Data and Figures From Paper
The figures presented in [Gorkowski et al. (2019)](https://doi.org/10.5194/acp-2019-495) can be reproduced using the source code in [paper_figures](https://github.com/Gorkowski/Binary_Activity_Thermodynamics_Model/tree/master/Matlab_source_code/paper_figures). The data shown in the paper and the fit parameters are compiled in the [Excel workbook file](https://github.com/Gorkowski/Binary_Activity_Thermodynamics_Model/tree/master/Matlab_source_code/paper_figures/Figure_data).
