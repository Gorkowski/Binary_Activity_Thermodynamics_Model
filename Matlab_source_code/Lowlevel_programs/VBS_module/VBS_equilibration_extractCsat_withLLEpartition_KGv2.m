function [Csat_approx]=VBS_equilibration_extractCsat_withLLEpartition_KGv2(Cp_j_VBSold, Cstar_j_VBSold, ...
     aw_measurment, molecular_weight, O2C_values, H2C_values,  McGlashan_mode, BAT_refinement_mode,N2C_values_denistyOnly )
%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2018-Jul-13 11:21 AM
% Copyright 2018 Kyle Gorkowski 
%% 
if not(exist('N2C_values_denistyOnly'))
    N2C_values_denistyOnly=molecular_weight.*0;
end

S=size(Cp_j_VBSold);

Molar_mass_ratios=18.016./molecular_weight;

aw_vec=ones(S(1,1),1).*aw_measurment;

McGlashan_fit_tolerance=10^-8;

[~, mole_frac_org_beta] = inverted_NNMcGlashan_v8(O2C_values,H2C_values, Molar_mass_ratios, aw_vec, McGlashan_mode);


mass_fraction_water_beta=molecular_weight.*0;
activity_coefficient_beta=mass_fraction_water_beta;
for i=1:S(1,1)
    
           mole_frac_bounds_beta=[0,1];

        if strcmpi(BAT_refinement_mode,'interpolate')

              [~, func2, ycalc1, activity_coefficient_beta(i,1), activity_calc2_water, activity_calc2_beta, mass_fraction_water_beta(i,1), mass_fraction2,~,~, ~ ]=...
                BAT_activity_calc_with_refinement_v1(mole_frac_org_beta(i,1), O2C_values(i,1), H2C_values(i,1),...
                Molar_mass_ratios(i,1),McGlashan_mode,[], [BAT_refinement_mode,'beta'], aw_vec(i,1), N2C_values_denistyOnly);
  
        else
            [~, func2, ycalc1, activity_coefficient_beta(i,1), activity_calc2_water, activity_calc2_beta, mass_fraction_water_beta(i,1), mass_fraction2,~,~, ~ ]=...
                BAT_activity_calc_with_refinement_v1(mole_frac_org_beta(i,1), O2C_values(i,1), H2C_values(i,1),...
                Molar_mass_ratios(i,1),McGlashan_mode,[], BAT_refinement_mode,  aw_vec(i,1), N2C_values_denistyOnly);
        end
end



q_beta=1;

%% New Coa j
Coa_j=Cp_j_VBSold;


% beta 
Coa_j_beta=Coa_j.*q_beta;
Caq_beta=sum(Coa_j_beta.*mass_fraction_water_beta./(1-mass_fraction_water_beta));
Coaaq_beta=sum(Coa_j_beta)+Caq_beta;

Coaaq_via_Ej_guess=Coaaq_beta;


% C* via beta phase
massweighted_molar_weight_beta=sum(Coa_j_beta./(molecular_weight))+Caq_beta./18.015;

    %Cstar_j_VBSold=Cstar_dry.*activity_coefficient_beta.*q_beta.*Coaaq_via_Ej_guess./(molecular_weight.*massweighted_molar_weight_beta); % calc new Cstar

    Csat_approx=Cstar_j_VBSold./(activity_coefficient_beta.*q_beta.*Coaaq_via_Ej_guess./(molecular_weight.*massweighted_molar_weight_beta));




