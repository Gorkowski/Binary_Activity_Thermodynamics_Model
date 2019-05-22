 function [phaseSep_via_activity,phaseSep_via_activity_curvature,index_phase_sep_starts,index_phase_sep_end]=finds_PhaseSep_and_activity_curve_dips_v2(activity_data)

%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2018-Oct-11 11:01 AM
% Copyright 2018 Kyle Gorkowski
%%

% % EXAMPLE
% clear
% 
% fix_molarmass_ratio_temp=    0.01;
% 
% H2C=0;
% BAT_functional_group='initial fitting for biphasic transition'; % 'SOA chemicals' 'standard chemicals' 'initial fitting for biphasic transition'
% 
% 
% O2C_values_temp =0.2;
% 
% mole_frac_grid=[1:-0.0001:0]';
% 
% [func1, func2, ycal_water, ycalc_org, activity_water, activity_org, mass_fraction1, mass_fraction2,Gibbs_RT, dGibbs_RTdx2]...
%     =mcglashan_activity_calc_KGv8(mole_frac_grid, O2C_values_temp, H2C, fix_molarmass_ratio_temp, BAT_functional_group,[]);
% activity_data=activity_org;
% %activity_data=activity_water;


%% Program


activity_diff=diff(activity_data);
L_m=length(activity_data);

% figure
% plot(activity_data);
% figure
% plot(activity_diff);

if L_m>3

[min_value,min_index]=min(activity_diff);
[max_value,max_index]=max(activity_diff);
mean_sign=sign(mean(activity_diff));

% if the max and min signs are both the same then there is no dip in the
% activity curve
if sign(min_value)==sign(max_value)
    phaseSep_via_activity_curvature=0;
    
            index_phase_sep_starts=NaN;
        index_phase_sep_end=NaN;
        
elseif sign(min_value)~=sign(max_value)
    phaseSep_via_activity_curvature=1;

    % get mole fraction value at start or phasesep and end of phasesep

    % curve change
    activity_calc1_diff_sign_change=sign([activity_diff(1,1);activity_diff])~=sign(activity_diff(1,1));
    
    index_start=find(activity_calc1_diff_sign_change,1);
    back_index=index_start-1+find(not(activity_calc1_diff_sign_change(index_start:end,1)),1);
    
    if back_index<L_m
        [~,activity_data_gap]=min(abs(activity_data(back_index:end,1)-activity_data(index_start,1)));
        restart_match_index=activity_data_gap+back_index-1;
    else
        restart_match_index=L_m;
    end
    
    % activity change, greater than 1
    if sum(activity_data>1)
        
        % get min RH in high mole fraction region.
        
        [min_value_Idilute,min_index_Idilute]=min(activity_data(index_start:end));
        min_index_Idilute=min_index_Idilute+index_start-1;
        
        % git mole fraction in low mole fraction region
        [~,activity_data_gap_start]=min(abs(activity_data(1:index_start,1)-activity_data(min_index_Idilute,1)));
        
        
        %output
        % check which region starts earlier and ends earlier
        if activity_data_gap_start<index_start
            index_phase_sep_starts=activity_data_gap_start;
        else
            index_phase_sep_starts=index_start;
        end
        
        if min_index_Idilute<restart_match_index
            index_phase_sep_end=min_index_Idilute;
        else
            index_phase_sep_end=restart_match_index;
        end
        
    else
        %output
        index_phase_sep_starts=index_start;
        
        index_phase_sep_end=restart_match_index;
        
    end
    
    
end

else
    phaseSep_via_activity=activity_data;
    phaseSep_via_activity_curvature=0;
    
    index_phase_sep_starts=NaN;
    index_phase_sep_end=NaN;
        

end

if sum(activity_data>1)
    phaseSep_via_activity=1;
else
    phaseSep_via_activity=0;
end



