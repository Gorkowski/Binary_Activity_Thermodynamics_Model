function [func1, func2, ycal_water, ycalc_org, activity_water, activity_org, mass_fraction1, mass_fraction2,Gibbs_RT, dGibbs_RTdx2, mole_frac_fit, error_out ]=...
    BAT_activity_calc_with_refinement_v1(mole_frac, O2C_values, H2C_values, molarmass_ratio, BAT_functional_group, McGlashan_special_options, ...
    refinement_mode, aw_desired, N2C_values_denistyOnly);

%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2018-Jul-27  1:48 PM
% Copyright 2018 Kyle Gorkowski
%% BAT optional solver
if not(exist('N2C_values_denistyOnly'))
    N2C_values_denistyOnly=molarmass_ratio.*0;
end

%interpolate number of steps steps the resolution interpolate uses, this could be improved
interpolate_step_numb=500; 

if strcmpi(refinement_mode,'none')
    [func1, func2, ycal_water, ycalc_org, activity_water, activity_org, mass_fraction1, mass_fraction2,Gibbs_RT, dGibbs_RTdx2]...
        =BAT_properties_calculation_v1(mole_frac, O2C_values, H2C_values, molarmass_ratio,...
        BAT_functional_group,McGlashan_special_options, N2C_values_denistyOnly);
    
    error_out=mole_frac.*NaN;
    mole_frac_fit=mole_frac;
    

elseif strcmpi(refinement_mode,'interpolatebeta')
    
    mole_frac=[0:1/interpolate_step_numb:1]';
    
    [func1, func2, ycal_water, ycalc_org, activity_water, activity_org, mass_fraction1, mass_fraction2,Gibbs_RT, dGibbs_RTdx2]...
        =BAT_properties_calculation_v1(mole_frac,O2C_values, H2C_values, molarmass_ratio, BAT_functional_group,McGlashan_special_options);

    [phaseSep_via_activity,phaseSep_via_activity_curvature,index_phase_sep_starts,index_phase_sep_end]=finds_PhaseSep_and_activity_curve_dips_v2(activity_water);
%     [phase_sep_check,lower_a_w_sep_index,upper_a_w_sep_index, matching_Upper_a_w_sep_index]=finds_PhaseSep_w_and_org(activity_water,activity_org);
% check how the order of this is working 
    if phaseSep_via_activity_curvature==1
        if index_phase_sep_end<length(activity_water) % fixed error for always sep
            mole_frac_fit=interp1(activity_water(index_phase_sep_end:end,1), mole_frac(index_phase_sep_end:end,1),aw_desired,'linear');
        else
            mole_frac_fit=mole_frac(end,1);
        end
    else
        
        mole_frac_fit=interp1(activity_water(1:end,1), mole_frac(1:end,1),aw_desired,'linear');
    end
    
    if isnan(mole_frac_fit)
        mole_frac_fit=mole_frac(index_phase_sep_end,1);
    end
    [~,aw_index]=min(abs(mole_frac-mole_frac_fit));
    
    %final output
    [func1, func2, ycal_water, ycalc_org, activity_water, activity_org, mass_fraction1, mass_fraction2,Gibbs_RT, dGibbs_RTdx2]...
        =BAT_properties_calculation_v1(mole_frac_fit, O2C_values, H2C_values, molarmass_ratio, BAT_functional_group,McGlashan_special_options);
    error_out=abs(aw_desired-activity_water)./aw_desired;
    
    
elseif strcmpi(refinement_mode,'interpolatealpha')
    mole_frac=[0:1/interpolate_step_numb:1]';
    
    [func1, func2, ycal_water, ycalc_org, activity_water, activity_org, mass_fraction1, mass_fraction2,Gibbs_RT, dGibbs_RTdx2]...
        =BAT_properties_calculation_v1(mole_frac, O2C_values, H2C_values, molarmass_ratio,...
        BAT_functional_group,McGlashan_special_options, N2C_values_denistyOnly);
    
    [phaseSep_via_activity,phaseSep_via_activity_curvature,index_phase_sep_starts,index_phase_sep_end]=finds_PhaseSep_and_activity_curve_dips_v2(activity_water);
    
    if phaseSep_via_activity_curvature==1
        if index_phase_sep_end<length(activity_water) % fixed error for always sep
            mole_frac_fit=interp1(activity_water(1:index_phase_sep_starts,1), mole_frac(1:index_phase_sep_starts,1),aw_desired,'linear');
        else
            mole_frac_fit=mole_frac(end,1);
        end
    else
        mole_frac_fit=interp1(activity_water(1:end,1), mole_frac(1:end,1),aw_desired,'linear');
    end
    
    if isnan(mole_frac_fit)
        mole_frac_fit=mole_frac(index_phase_sep_starts,1);
    end
    [~,aw_index]=min(abs(mole_frac-mole_frac_fit));
    
    %final output
    [func1, func2, ycal_water, ycalc_org, activity_water, activity_org, mass_fraction1, mass_fraction2,Gibbs_RT, dGibbs_RTdx2]...
        =BAT_properties_calculation_v1(mole_frac_fit, O2C_values, H2C_values, molarmass_ratio,...
        BAT_functional_group,McGlashan_special_options, N2C_values_denistyOnly);
    
    error_out=abs(aw_desired-activity_water)./aw_desired;
    
else
    error('pick McGlashan refinement method')
end



end

