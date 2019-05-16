# Binary Activity Thermodynamics Model
Water plays an essential role in aerosol chemistry, gas-to-particle partitioning, and particle viscosity, but it is not typically included in organic aerosol thermodynamic models.  
We have developed a water-sensitive organic aerosol thermodynamic model that can adapt to different levels of chemical information. 
Our model is capable of using both detailed molecular structures as well as bulk properties like O:C and component volatility classes. 
The first development stage is accounting for non-ideal interactions of organic aerosol components with water. 
With the non-ideal water-organic interactions accounted for we can realistically simulate organic aerosol using ambient measurements. 
This capability uses our new Binary Activity Thermodynamics (BAT) modeling framework, which describes a mixture of organics with water as binary mixtures of one organic and associated water. 
For gas-particle partitioning, we use a non-ideal Volatility Basis Set (VBS) formulation that includes activity coefficients and accounts for liquid-liquid equilibrium. 

<!-- If you use the model please cite/reference our published paper XXX -->

## Code and File Overview
The source code in written in MATLAB 2018b, though it will likely run with past and future versions. 
The source code is in "\Matlab_source_code", which contains "Matlab_source_code\update_Matlab_paths.m" function that should be run to add the source folder to your MATLAB path.

In addition to the source code, a standalone executable file has be compiled and can run the model if you have no MATLAB license.

## Standalone Executable Steps
The installation file is "Matlab_runtime\BAT_Model\for_redistribution\BAT_installer_v1_web.exe", which needs the freely available MATLAB Runtime R2018b engine. 
The installation file will download the 

## Generate Figures From Paper

<!-- ## Data Repository TBD-->
