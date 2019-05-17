function [partition_coefficents_AB, Coa_j_AB, Caq_j_AB, Cstar_j, Coa_AB,...
    Caq_AB, Coa, q_alpha_water, fit_exit_flag, error_out]=VBS_equilibration_withLLEpartition_KGv2(...
    guess_C_OAalpha_ugPm3,guess_C_OAbeta_ugPm3,guess_partition_coefficents, ...
    C_OM_ugPm3, Cstar_dry, activity_coefficient_AB, q_alpha_molefrac_phase_split_org, ...
    mass_fraction_water_AB, molecular_weight, a_water, O2C_values, BAT_functional_group, ...
    VBSBAT_options)
%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2018-Apr-10 11:48 AM
% Copyright 2018 Kyle Gorkowski
%%
% Solves for VSB partition coefficients Ej

%% get Ej guess
Molar_mass_raitos=18.015./molecular_weight;
if guess_C_OAbeta_ugPm3>0
    guess_C_oaaqbeta_ugPm3=guess_C_OAbeta_ugPm3./(1-mean(mass_fraction_water_AB(:,2)));
elseif guess_C_OAalpha_ugPm3>0
    guess_C_oaaqbeta_ugPm3=guess_C_OAalpha_ugPm3./(1-mean(mass_fraction_water_AB(:,2)));
else
    guess_C_oaaqbeta_ugPm3=1./(1-mean(mass_fraction_water_AB(:,2)));
end

if sum(guess_partition_coefficents)==0 || strcmpi(VBSBAT_options.optimization.independent_aw,'yes')
    
    if strcmpi(VBSBAT_options.VBSBAT_NN_options.use_NN_for_VBS_initial_guess,'no') % guess for Ej using approximation
        
        Cstar_j_guess=Cstar_dry.*activity_coefficient_AB(:,2)./(1+mass_fraction_water_AB(:,2).*(1./Molar_mass_raitos-1)); %not right but fine for guess
        Ej_guess=(1+Cstar_j_guess./guess_C_oaaqbeta_ugPm3).^-1;
        
    elseif strcmpi(VBSBAT_options.VBSBAT_NN_options.use_NN_for_VBS_initial_guess,'yes') % guess for Ej using neural network
        
        % network trained on hydroxyl so need to convert to OH first
        [O2C_temp,Mratio_temp] = convert_chemical_structure_to_OH_eqv_v3(O2C_values,Molar_mass_raitos, BAT_functional_group );
        
        [Ej_guess] =  VBSBAT_nerual_network_estimate_v1(...
            Cstar_dry, C_OM_ugPm3, O2C_temp, 18.015./Mratio_temp, mass_fraction_water_AB(:,2)...
            ,a_water, VBSBAT_options.VBSBAT_NN_options);
    end
    
    
elseif strcmpi(VBSBAT_options.optimization.independent_aw,'no') % use previous Ej_guess results
    Ej_guess=guess_partition_coefficents;
else
    error(' Select VBSBAT_options.optimization.independent_aw option either yes or no')
end

% error in guess runs the cost funciton to see if there is an error in
% guess.
guess_error=nested_opt_cost_sub1(Ej_guess, C_OM_ugPm3, Cstar_dry, activity_coefficient_AB,...
    q_alpha_molefrac_phase_split_org, mass_fraction_water_AB, molecular_weight);
% guess_CoaEj=[Ej_guess];

if VBSBAT_options.optimization.guess_refinement_threshold<=guess_error % check if refinement is needed
    
    S_om=size(C_OM_ugPm3);
    bounds_lower=[ones(S_om(1,1),1).*10^-10]; % Coa qAlpha and Ej bounds
    bounds_upper=[ones(S_om(1,1),1)];
    
    %options for fmincon optimization
    options=optimoptions('fmincon');
    options.Display='off';
    options.Algorithm=VBSBAT_options.optimization.fit_Algorithm; %'interior-point'; % 'interior-point' (default) 'trust-region-reflective' 'sqp' 'active-set'
    options.TolX=VBSBAT_options.optimization.fit_tolerance;
    options.MaxIter=VBSBAT_options.optimization.MaxIter;
    
    problem = createOptimProblem('fmincon','objective',...
        @(guess_CoaEj)nested_opt_cost_sub1(guess_CoaEj, C_OM_ugPm3, Cstar_dry, activity_coefficient_AB,...
        q_alpha_molefrac_phase_split_org, mass_fraction_water_AB, molecular_weight),...
        'x0',Ej_guess,'lb',bounds_lower,'ub',bounds_upper,'options',options);
    
    % possible different optimization models
    if strcmpi(VBSBAT_options.optimization.opt_method,'fmincon') % a single initialization point which is default
        [fit_CoaEj,~,fit_exit_flag] = fmincon(problem);
        
    elseif strcmpi(VBSBAT_options.optimization.opt_method,'global') % multiple start points taken at random
        gs = GlobalSearch('Display', 'off','TolX',VBSBAT_options.optimization.fit_tolerance,...
            'NumTrialPoints', VBSBAT_options.optimization.global_NumTrialPoints,'NumStageOnePoints',...
            VBSBAT_options.optimization.global_NumStageOnePoints);
        [fit_CoaEj,~,fit_exit_flag]  = run(gs,problem);
        
    elseif strcmpi(VBSBAT_options.optimization.opt_method,'none')
        fit_CoaEj=Ej_guess;
        fit_exit_flag=NaN;
    else
        error('Invalid VBSBAT_options.optimization.opt_method... possible options fmincon, global, none')
    end
    
else % guess has a lower error than guess_refinement_threshold
    fit_CoaEj=Ej_guess;
    fit_exit_flag=NaN;
end

% runs a final calcuation for output
[partition_coefficents_AB, Coa_j_AB, Caq_j_AB, Cstar_j,  Coa_AB, Caq_AB, Coa, ...
    q_alpha_water, error_out]= VBS_equilibration_withLLEpartition_objFun_KGv3(fit_CoaEj, ...
    C_OM_ugPm3, Cstar_dry, activity_coefficient_AB, q_alpha_molefrac_phase_split_org, mass_fraction_water_AB, molecular_weight);

end


% cost function for optimization 
function [error_out]= nested_opt_cost_sub1(guess_CoaEj, ...
    C_OM_ugPm3, Cstar_dry, activity_coefficient_AB, q_alpha_molefrac_phase_split_org, mass_fraction_water_AB, molecular_weight)


[partition_coefficents_AB, Coa_j_AB, Caq_j_AB, Cstar_j, Coa_AB, Caq_AB, Coa, q_alpha_water, error_out]= VBS_equilibration_withLLEpartition_objFun_KGv3(guess_CoaEj, ...
    C_OM_ugPm3, Cstar_dry, activity_coefficient_AB, q_alpha_molefrac_phase_split_org, mass_fraction_water_AB, molecular_weight);

end



