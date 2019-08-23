function [mass_frac1, mass_frac2] = mole_to_mass_fraction(mole_frac1, mole_frac2, molar_mass1, molar_mass2)
%%
% Created by Kyle Gorkowski [GORKOWFALCON] on 2019-Aug-18  6:23 AM
% Copyright 2019 Kyle Gorkowski 
%% 

mass_frac1=mole_frac1.*molar_mass1./(mole_frac1.*molar_mass1+mole_frac2.*molar_mass2);
mass_frac2=mole_frac2.*molar_mass2./(mole_frac1.*molar_mass1+mole_frac2.*molar_mass2);


end

