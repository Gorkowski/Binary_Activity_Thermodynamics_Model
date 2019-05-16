function [Csat_ugPm3, Csat_ugPm3_log10] = Csat_calculate_v1(Molar_mass_gPmol,saturation_vapor_pressure_Pa, Temperature_K, optional_convert_Psat_from_mmHg)
% [Csat_ugPm3, Csat_ugPm3_log10] = Csat_calculate_v1(Molar_mass_gPmol,saturation_vapor_pressure_Pa, Temperature_K, optional_convert_Psat_from_mmHg)
%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2019-Feb-05 10:53 AM
% Copyright 2019 Kyle Gorkowski 
%% 
if not(exist('optional_convert_Psat_from_mmHg'))
    % set default options
    optional_convert_Psat_from_mmHg='no';
end

if strcmpi(optional_convert_Psat_from_mmHg, 'yes')
    saturation_vapor_pressure_Pa= saturation_vapor_pressure_Pa.*133.322365; % convert  1 mmHg = 133.322365 pa
end
    

gas_constant_m3Pa_per_Kmol=8.3144598; %8.3144598(48)	m3?Pa / K?1?mol?1
Csat_gPm3 = Molar_mass_gPmol.*saturation_vapor_pressure_Pa./(gas_constant_m3Pa_per_Kmol.*Temperature_K);

Csat_ugPm3 = Csat_gPm3.*10^6;

Csat_ugPm3_log10=log10_zeropass(Csat_ugPm3);

end

