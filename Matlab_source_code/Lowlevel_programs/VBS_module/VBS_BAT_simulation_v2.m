function [C_OA_PM, Caq_PM, kappaHGF, details]=VBS_BAT_simulation_v2(Csat_j_value, C_OM_ugPm3, ...
    O2C_values, H2C_values,  Molecular_weight, aw_series, BAT_functional_group,...
    BAT_refinement_mode, VBSBAT_options, sim_name, N2C_values_denistyOnly)
% % % %
% % % %%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2018-Apr-10  1:28 PM
% Copyright 2018 Kyle Gorkowski
% Updated by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2018-Oct-26  3:38 PM
% Updated by Kyle Gorkowski [GORKOWFALCON] on 2018-Dec-16  1:38 PM
% 
% Purpose: This is the main program used in the VBS+BAT simulations. There are 4 main program steps
% Step 1: calculates the a_w separation point
% Step 2: calculates the activity coefficients at a given a_w/RH using the BAT model
% Step 3: calculates the equilibrium aerosol mass using VBS+BAT
% Step 4: calculates kappa HGF and collects/writes outputs 
%
% program starts*************

if not(exist('sim_name'))
    sim_name=['Simulation_' datestr(now,30)];
end
if not(exist('N2C_values_denistyOnly'))
    N2C_values_denistyOnly=O2C_values.*0;
end

tic % cpu timer starts

S=size(Csat_j_value);
special_options=[];

[~,min_aw_i]=min(aw_series); % min aw index
S_full=size(aw_series);

Molar_mass_ratios=18.016./Molecular_weight; % molecular weight ratios
[BAT_functional_group] = check_BAT_functional_group_inputs_v1(O2C_values, BAT_functional_group);

%%***Step 1***** check for phase separation conditions and calc q alpha values
[O2C_temp,Mratio_temp] = convert_chemical_structure_to_OH_eqv_v3(O2C_values,Molar_mass_ratios, BAT_functional_group );
O2C_single_phase_cross_point=O2C_temp<single_phase_O2C_point_KGv3(Mratio_temp);
mean_prop_mask=O2C_single_phase_cross_point;

if strcmpi(VBSBAT_options.force_phase.onePhase, 'no') % includes two phase options
    if sum(O2C_single_phase_cross_point)>0 % has a separation activity
        
        % get a_w_sep_point
        [a_w_sep_point] = biphasic_to_single_phase_RH_master_v4(O2C_values, H2C_values, Molar_mass_ratios, BAT_functional_group);
        
        % isolates miscible
        a_w_sep_point=a_w_sep_point.*O2C_single_phase_cross_point;
        
        % create a_w_sep_point matrix to match aw_series length and
        % chemical species width
        a_w_sep_matrix=repmat(a_w_sep_point',S_full(1,1),1);
        aw_series_matrix=repmat(aw_series,1,S(1,1));

        %calculate q_alpha as a function of a_w
        [q_alpha_vsRH_values] = q_alpha_transfer_vs_aw_calc_v1(a_w_sep_matrix, aw_series_matrix, VBSBAT_options);
        
        % threshold start and end q_alpha values
        max_q_alpha=q_alpha_vsRH_values>VBSBAT_options.q_alpha.q_alpha_bounds(1,1);
        min_q_alpha=q_alpha_vsRH_values<VBSBAT_options.q_alpha.q_alpha_bounds(1,2);
        
        % applies q_alpha limits 
        q_alpha_vsRH_values=max_q_alpha.*q_alpha_vsRH_values;
        q_alpha_vsRH_values=min_q_alpha.*q_alpha_vsRH_values+not(min_q_alpha);
        
    else
        q_alpha_vsRH_values=ones(S_full(1,1),S(1,1)); % fully miscible
    end
elseif strcmpi(VBSBAT_options.force_phase.onePhase, 'beta')
    q_alpha_vsRH_values=zeros(S_full(1,1),S(1,1)); % organic phase only
elseif strcmpi(VBSBAT_options.force_phase.onePhase, 'alpha')
    q_alpha_vsRH_values=ones(S_full(1,1),S(1,1)); % water rich phase only
else
    error('select VBSBAT_options.force_phase.onePhase')
end

% VBS calc start
%% matrix and vectors to save output data
partition_coefficients_out=zeros(S_full(1,1),S(1,1));
Cstar_j_out=partition_coefficients_out;
Coa_j_alpha=partition_coefficients_out;
Coa_j_beta=partition_coefficients_out;
Caq_j_alpha=partition_coefficients_out;
Caq_j_beta=partition_coefficients_out;
gamma_alpha=partition_coefficients_out;
gamma_beta=partition_coefficients_out;
mass_fraction_water_alpha_out=partition_coefficients_out;
mass_fraction_water_beta_out=partition_coefficients_out;

C_OA_out=zeros(S_full(1,1),1);
Coa_alpha=C_OA_out;
Coa_beta=C_OA_out;
q_alpha_water=C_OA_out;
q_alpha_water_out=C_OA_out;
weight_q_alpha=C_OA_out;
a_w_sep_point_of_meanPM=C_OA_out;

ycalc_org_beta_temp=C_OA_out;
mass_fraction_water_beta_temp=C_OA_out;
ycalc_org_alpha_temp=C_OA_out;
mass_fraction_water_alpha_temp=C_OA_out;
fit_exit_flag_save=C_OA_out;
error_save=C_OA_out;

VBSBAT_fit_counter=0; % some timing counters
toc_BAT_cumlitive=0;
toc_VBSBAT_cumlitive=0;
BAT_fit_counter=0;

if sum(C_OM_ugPm3)==0
    error('zero concentration of all species')
end

for s_i=1:S_full(1,1) % iterates through aw values

    aw=aw_series(s_i,1); % selects aw
    
    aw_vec=ones(S(1,1),1).*aw; % duplicates aw to match organic numb
    
    %%***Step 2*****  
    %nn version gets mole fraction at aw
    [mole_frac_org_alpha, mole_frac_org_beta] = inverted_NNMcGlashan_v8(O2C_values, H2C_values, Molar_mass_ratios, aw_vec, BAT_functional_group);
    
    mass_fraction_water_alpha=Molecular_weight.*0;
    ycalc_org_alpha=mass_fraction_water_alpha;
    mass_fraction_water_beta=Molecular_weight.*0;
    ycalc_org_beta=mass_fraction_water_alpha;
    for i=1:S(1,1)
        mole_frac_bounds_alpha=[0,1];
        mole_frac_bounds_beta=[0,1];
        
        if VBSBAT_options.BAT_refinement_aw >= aw % sets when the refinement should start
            BAT_refinement_mode_temp='none';
        else
            BAT_refinement_mode_temp=BAT_refinement_mode;
        end
        
        toc_BAT_start=toc;
        if strcmpi(BAT_refinement_mode_temp,'interpolate')
            
            [~, ~, ~, ycalc_org_alpha(i,1), ~, activity_calc2_alpha, mass_fraction_water_alpha(i,1), ~,~,~, ~ ]=...
                BAT_activity_calc_with_refinement_v1(mole_frac_org_alpha(i,1), O2C_values(i,1), H2C_values(i,1), Molar_mass_ratios(i,1),BAT_functional_group(i,1),special_options,...
                [BAT_refinement_mode_temp,'alpha'],  aw, N2C_values_denistyOnly(i,1));
            
            [~, ~, ~, ycalc_org_beta(i,1), ~, activity_calc2_beta, mass_fraction_water_beta(i,1), ~,~,~, ~ ]=...
                BAT_activity_calc_with_refinement_v1(mole_frac_org_beta(i,1), O2C_values(i,1), H2C_values(i,1), Molar_mass_ratios(i,1),BAT_functional_group(i,1),special_options,...
                [BAT_refinement_mode_temp,'beta'],  aw, N2C_values_denistyOnly(i,1));
            
        else
            
            [~, ~, ~, ycalc_org_alpha(i,1), ~, activity_calc2_alpha, mass_fraction_water_alpha(i,1), ~,~,~, ~ ]=...
                BAT_activity_calc_with_refinement_v1(mole_frac_org_alpha(i,1), O2C_values(i,1), H2C_values(i,1), Molar_mass_ratios(i,1),BAT_functional_group(i,1),special_options,...
                BAT_refinement_mode_temp,  aw, N2C_values_denistyOnly(i,1));
            
            [~, ~, ~, ycalc_org_beta(i,1), ~, activity_calc2_beta, mass_fraction_water_beta(i,1), ~,~,~, ~ ]=...
                BAT_activity_calc_with_refinement_v1(mole_frac_org_beta(i,1), O2C_values(i,1), H2C_values(i,1), Molar_mass_ratios(i,1),BAT_functional_group(i,1),special_options,...
                BAT_refinement_mode_temp,  aw, N2C_values_denistyOnly(i,1));
        end
        BAT_fit_counter=BAT_fit_counter+2;
        toc_BAT_cumlitive=toc_BAT_cumlitive+toc-toc_BAT_start;
        
    end
    
    activity_coefficient_AB=[ycalc_org_alpha,ycalc_org_beta];
    mass_fraction_water_AB=[mass_fraction_water_alpha,mass_fraction_water_beta];
    
    %%***Step 3***** 

    % guess options for VBS start point 
    % not used if the NN is used for the initial guess (in VBSBAT_options)
    % but we still need to pass a value to it.  
    if s_i>1
        guess_C_OAalpha_ugPm3=Coa_alpha(s_i-1,:);
        guess_C_OAbeta_ugPm3=Coa_beta(s_i-1,:);
        guess_partition_coefficients=partition_coefficients_out(s_i-1,:)';
%         guess_q_alpha_water=q_alpha_water_out(s_i-1,:);
    else
        guess_C_OAalpha_ugPm3=.1;
        guess_C_OAbeta_ugPm3=.1;
        guess_partition_coefficients=partition_coefficients_out(1,:)';
%         guess_q_alpha_water=q_alpha_water_out(1,:);
    end
    
    toc_VBSBAT_start=toc;
    if sum(q_alpha_vsRH_values(s_i,:))==0      % single phase calculation BETA
        
        %         %beta phase calc only
        [partition_coefficients_temp, Coa_j_AB, Caq_j_AB, Cstar_j, Coa_AB, Caq_AB,...
            C_OA_PM, q_alpha_water, fit_exit_flag, error_out]=VBS_equilibration_withLLEpartition_KGv2(...
            guess_C_OAalpha_ugPm3,guess_C_OAbeta_ugPm3,guess_partition_coefficients,...
            C_OM_ugPm3, Csat_j_value, activity_coefficient_AB, q_alpha_vsRH_values(s_i,:)', mass_fraction_water_AB, ...
            Molecular_weight, aw, O2C_values, BAT_functional_group,...
            VBSBAT_options);
        
        partition_coefficients=partition_coefficients_temp; % outputs
        weight_q_alpha(s_i,1)=q_alpha_vsRH_values(s_i,1);
        VBSBAT_fit_counter=VBSBAT_fit_counter+1;
        
    elseif sum(q_alpha_vsRH_values(s_i,:))==S(1,1)     % single phase calculation ALPHA
        
        % alpha phase calc only
        [partition_coefficients_temp, Coa_j_AB, Caq_j_AB,Cstar_j, Coa_AB, Caq_AB,...
            C_OA_PM, q_alpha_water, fit_exit_flag, error_out]=VBS_equilibration_withLLEpartition_KGv2(...
            guess_C_OAalpha_ugPm3,guess_C_OAbeta_ugPm3,guess_partition_coefficients,...
            C_OM_ugPm3, Csat_j_value, activity_coefficient_AB, q_alpha_vsRH_values(s_i,:)', mass_fraction_water_AB, ...
            Molecular_weight,aw,O2C_values, BAT_functional_group,...
            VBSBAT_options);
        
        partition_coefficients=partition_coefficients_temp; % outputs
        weight_q_alpha(s_i,1)=q_alpha_vsRH_values(s_i,1);
        VBSBAT_fit_counter=VBSBAT_fit_counter+1;
        
    else % q alpha has some value
        
        % beta phase calc first which is used for average q_alpha calc and
        % used as a check to see if the two-phase simulation is realistic 
        
        guess_C_OA_ugPm3=max(guess_C_OAbeta_ugPm3,guess_C_OAalpha_ugPm3); % not used if the NN is used for the initial guess (in VBSBAT_options)
        
        [partition_coefficients_temp, Coa_j_AB_temp, Caq_j_AB_temp, Cstar_j_temp, Coa_AB_temp,...
            Caq_AB_temp, C_OA_temp, q_alpha_water_temp, fit_exit_flag, error_out]=VBS_equilibration_withLLEpartition_KGv2(...
            guess_C_OA_ugPm3,guess_C_OA_ugPm3,guess_partition_coefficients,...
            C_OM_ugPm3, Csat_j_value, activity_coefficient_AB, q_alpha_vsRH_values(s_i,:)'.*0, mass_fraction_water_AB, ...
            Molecular_weight,aw,O2C_values, BAT_functional_group,...
            VBSBAT_options);
        
        VBSBAT_fit_counter=VBSBAT_fit_counter+1;
        
        %% old mass weighted averages 
%         mass_inPM_temp=partition_coefficients_temp.*C_OM_ugPm3.*mean_prop_mask; % mask for phase sep org add .*mean_prop_mask
%         mass_fraction_inPM=(mass_inPM_temp)./C_OA_temp;
%         
%         %% mean prop used in mean q_alpha calc
%          Mratio_temp=18.015./sum(mass_fraction_inPM.* Molecular_weight);
%          O2C_temp=sum(mass_fraction_inPM.* O2C_values);
%         H2C_temp=sum(mass_fraction_inPM.* H2C_values);
        
        %molar based averages
        [avg_M_gPmol, O2C_temp, H2C_temp, ~, ~] = ...
    molar_based_means(Coa_j_AB_temp(:,2)', Molecular_weight, O2C_values, H2C_values, N2C_values_denistyOnly);
        Mratio_temp=18.015./avg_M_gPmol;

        [a_w_sep_point_of_meanPM(s_i,1)] = biphasic_to_single_phase_RH_master_v4(O2C_temp, H2C_temp, Mratio_temp, VBSBAT_options.mean_BAT_functional_group);
        
        [weight_q_alpha(s_i,1)] = q_alpha_transfer_vs_aw_calc_v1(a_w_sep_point_of_meanPM(s_i,1), aw, VBSBAT_options);
        
        % threshold q_alpha value for mean
        if weight_q_alpha(s_i,1)>VBSBAT_options.q_alpha.q_alpha_bounds_mean(1,2)
            weight_q_alpha(s_i,1)=1;
        end
        
        % make sure q_alpha mean is greater than mean min, to proceed with 2 phase
        % calc. as this must pass before a 2 phase calc is used with
        % individual q alphas
        if weight_q_alpha(s_i,1)>VBSBAT_options.q_alpha.q_alpha_bounds_mean(1,1) 
            
            % select whigh q alpha to use in two phase transition calc
            if strcmpi(VBSBAT_options.q_alphaVBS.method_to_use,'mean_prop')
                q_alpha_molefrac_phase_split_org=ones(S(1,1),1).*weight_q_alpha(s_i,1);
                
            elseif strcmpi(VBSBAT_options.q_alphaVBS.method_to_use,'individual')
                q_alpha_molefrac_phase_split_org=q_alpha_vsRH_values(s_i,:)';
            else
                error('select q alpha method for transition')
            end
            
            % two phases calculation
            activity_coefficient_AB=[ycalc_org_alpha,ycalc_org_beta];
            mass_fraction_water_AB=[mass_fraction_water_alpha,mass_fraction_water_beta];
            
            [partition_coefficients, Coa_j_AB, Caq_j_AB, Cstar_j, Coa_AB, Caq_AB, C_OA_PM, q_alpha_water, fit_exit_flag, error_out]=VBS_equilibration_withLLEpartition_KGv2(...
                guess_C_OAalpha_ugPm3,guess_C_OAbeta_ugPm3,guess_partition_coefficients,...
                C_OM_ugPm3, Csat_j_value, activity_coefficient_AB, q_alpha_molefrac_phase_split_org, mass_fraction_water_AB,...
                Molecular_weight,aw,O2C_values, BAT_functional_group,...
                VBSBAT_options);
            
            VBSBAT_fit_counter=VBSBAT_fit_counter+1;
            
            
            % checks if the two phase mass calculation is realistic
            if strcmpi(VBSBAT_options.q_alphaVBS.Option_to_overwrite_two_phase_calculation , 'yes') ...
                    && (C_OA_PM-C_OA_temp)/C_OA_temp <VBSBAT_options.q_alphaVBS.overwrite_threshold_two_phase_fraction_difference
                
                
                % 2 phase calc difference is too low calc
                % the mean of two phase and beta only
                partition_coefficients=(partition_coefficients+partition_coefficients_temp)./2; % outputs
                Coa_j_AB=(Coa_j_AB_temp+Coa_j_AB)./2;
                Caq_j_AB=(Caq_j_AB_temp+Caq_j_AB)./2;
                Cstar_j=(Cstar_j_temp+Cstar_j)./2;
                Coa_AB=(Coa_AB+Coa_AB_temp)./2;
                C_OA_PM=(C_OA_PM+C_OA_temp)./2;
                q_alpha_water=(q_alpha_water+q_alpha_water_temp)./2;
                %disp('force_change')
            end
            
        else 
            % 2 phase calculation not needed
            partition_coefficients=partition_coefficients_temp; % outputs
            Coa_j_AB=Coa_j_AB_temp;
            Caq_j_AB=Caq_j_AB_temp;
            Cstar_j=Cstar_j_temp;
            Coa_AB=Coa_AB_temp;
            C_OA_PM=C_OA_temp;
            q_alpha_water=q_alpha_water_temp;
            weight_q_alpha(s_i,1)=0;
        end
        
    end
       % VBSBAT time cumulative sum
        toc_VBSBAT_cumlitive=toc_VBSBAT_cumlitive+toc-toc_VBSBAT_start;

    
    % outputs for this aw iteration
    partition_coefficients_out(s_i,:)=partition_coefficients(:,1);
    Cstar_j_out(s_i,:)=Cstar_j;
    Coa_j_alpha(s_i,:)=Coa_j_AB(:,1);
    Coa_j_beta(s_i,:)=Coa_j_AB(:,2);
    Caq_j_alpha(s_i,:)=Caq_j_AB(:,1);
    Caq_j_beta(s_i,:)=Caq_j_AB(:,2);
    gamma_alpha(s_i,:)=ycalc_org_alpha';
    gamma_beta(s_i,:)=ycalc_org_beta';
    mass_fraction_water_alpha_out(s_i,:)=mass_fraction_water_alpha';
    mass_fraction_water_beta_out(s_i,:)=mass_fraction_water_beta';

    q_alpha_water_out(s_i,1)=q_alpha_water;
    C_OA_out(s_i,:)=C_OA_PM;
    Coa_alpha(s_i,:)=Coa_AB(:,1);
    Coa_beta(s_i,:)=Coa_AB(:,2);
    
    error_save(s_i,:)=error_out;
    fit_exit_flag_save(s_i,1)=fit_exit_flag;
    
    
end

%%***Step 4***** 
% bulk properties
Coa_j_PM=Coa_j_alpha+Coa_j_beta;
Caq_j_PM=Caq_j_alpha+Caq_j_beta;

C_OA_PM=C_OA_out;
Caq_PM=sum(Caq_j_PM,2);
C_OA_ratio=C_OA_PM./C_OA_PM(min_aw_i,1);


% organic mass weighted O2C, MW, H:C
mass_fraction_inPM=Coa_j_PM./C_OA_PM;
[avg_M_gPmol, molar_avg_O2C, molar_avg_H2C, molar_avg_N2C, mass_weighted_avg_O2C] = ...
    molar_based_means(Coa_j_PM, Molecular_weight, O2C_values, H2C_values, N2C_values_denistyOnly);

%old method has a high O2C bias
% mass_weighted_avg_MW=sum(mass_fraction_inPM.* repmat(Molecular_weight',S_full(1,1),1),2);
% mass_weighted_avg_O2C=sum(mass_fraction_inPM.* repmat(O2C_values',S_full(1,1),1),2);
% mass_weighted_avg_H2C=sum(mass_fraction_inPM.* repmat(H2C_values',S_full(1,1),1),2);


% calc kappa
[kappaHGF,kappa, growth] = bulk_kappa_v1(O2C_values, H2C_values,...
    Molecular_weight, aw_series, Coa_j_PM, Caq_PM, VBSBAT_options.kappa_CCN_settings);

timer_end=toc;

%% collect outputs

% outputs that could be of use species specific
details.species_specific.partition_coefficients=partition_coefficients_out;
% details.species_specific.partition_coefficients_beta=partition_coefficients_beta;
details.species_specific.Cstar_j=Cstar_j_out;
details.species_specific.Coa_j_alpha=Coa_j_alpha;
details.species_specific.Coa_j_beta=Coa_j_beta;
details.species_specific.Caq_j_alpha=Caq_j_alpha;
details.species_specific.Caq_j_beta=Caq_j_beta;
details.species_specific.Coa_j_PM=Coa_j_PM;
details.species_specific.Caq_j_PM=Caq_j_PM;
details.species_specific.gamma_alpha=gamma_alpha;
details.species_specific.gamma_beta=gamma_beta;
details.species_specific.mass_fraction_water_alpha=mass_fraction_water_alpha_out;
details.species_specific.mass_fraction_water_beta=mass_fraction_water_beta_out;
details.species_specific.q_alpha_water=q_alpha_water_out;
details.species_specific.q_alpha_org_values=q_alpha_vsRH_values;
details.species_specific.q_alpha_org_used=weight_q_alpha;
details.species_specific.mass_fraction_inPM=mass_fraction_inPM;
% usefull calcs
Molecular_weight_withwater=[Molecular_weight;18.015];
details.species_specific.postEquilb.mole_fraction_water_free_alpha=(Coa_j_alpha./Molecular_weight')./sum(Coa_j_alpha./Molecular_weight',2);
details.species_specific.postEquilb.mole_fraction_water_free_beta=(Coa_j_beta./Molecular_weight')./sum(Coa_j_beta./Molecular_weight',2);
details.species_specific.postEquilb.mole_fraction_alpha=([Coa_j_alpha,sum(Caq_j_alpha,2)]./Molecular_weight_withwater')./sum([Coa_j_alpha,sum(Caq_j_alpha,2)]./Molecular_weight_withwater',2);
details.species_specific.postEquilb.mole_fraction_beta=([Coa_j_beta,sum(Caq_j_beta,2)]./Molecular_weight_withwater')./sum([Coa_j_beta,sum(Caq_j_beta,2)]./Molecular_weight_withwater',2);



% the totals
details.totals.C_OA=C_OA_out;
details.totals.Coa_alpha=Coa_alpha;
details.totals.Coa_beta=Coa_beta;
details.totals.C_OA_PM=C_OA_PM;
details.totals.Caq_PM=Caq_PM;
details.totals.Caq_PM_alpha=sum(details.species_specific.Caq_j_alpha,2);
details.totals.Caq_PM_beta=sum(details.species_specific.Caq_j_beta,2);
details.totals.C_OA_ratio=C_OA_ratio;
details.totals.avg_M_gPmol=avg_M_gPmol;
details.totals.mass_weighted_avg_O2C=mass_weighted_avg_O2C;
details.totals.molar_avg_O2C=molar_avg_O2C;
details.totals.molar_avg_H2C=molar_avg_H2C;
details.totals.molar_avg_N2C=molar_avg_N2C;
details.kappas.kappaHGF=kappaHGF;
details.kappas.kappa=kappa;
details.growth=growth;

% some fit outputs
details.fit.error_VBSBAT=error_save;
details.fit.VBSBAT_fit_exit_flag=fit_exit_flag_save;
details.fit.mean_time_per_fit=timer_end./VBSBAT_fit_counter;
details.fit.total_fits_sec=VBSBAT_fit_counter;
details.fit.total_time_sec=timer_end;
details.fit.VBSBAT_total=toc_VBSBAT_cumlitive;
details.fit.VBSBAT_mean_time_per_fit=toc_VBSBAT_cumlitive./VBSBAT_fit_counter;
details.fit.BAT_total=toc_BAT_cumlitive;
details.fit.BAT_mean_time_per_fit=toc_BAT_cumlitive./BAT_fit_counter;



%save inputs in output
details.inputs.Csat_j_value=Csat_j_value;
details.inputs.C_OM_ugPm3=C_OM_ugPm3;
details.inputs.O2C_values=O2C_values;
details.inputs.H2C_values=H2C_values;
details.inputs.Molecular_weight=Molecular_weight;
details.inputs.aw_series=aw_series;
details.inputs.BAT_functional_group=BAT_functional_group;
details.inputs.BAT_refinement_mode=BAT_refinement_mode;
details.inputs.VBSBAT_options=VBSBAT_options;

%% graph
if strcmpi(VBSBAT_options.plot_PM,'yes')
   plot_VBSBAT_simple_graph(details, VBSBAT_options.plot_PM, sim_name)
end

% save output files
if strcmpi(VBSBAT_options.write_simple_output,'yes')
    write_VBSBAT_output_files(sim_name,details)
end


