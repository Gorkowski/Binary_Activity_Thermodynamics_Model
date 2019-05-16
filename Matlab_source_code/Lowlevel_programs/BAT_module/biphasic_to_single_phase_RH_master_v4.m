function [RH_cross_point] = biphasic_to_single_phase_RH_master_v4(O2C, H2C, Mratio, BAT_functional_group, McGlashan_refinement_mode)

% Updated by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2018-Oct-26  3:38 PM
% Updated by Kyle Gorkowski [GORKOWFALCON] on 2018-Dec-16 10:38 AM
% always does interpolation method

RH_cross_point=O2C.*0;

interpolate_step_numb=500; % interpolation points
mole_frac=[0:1/interpolate_step_numb:1]';

for i=1:length(O2C)

% calculate activities 
    [func1, func2, ycal_water, ycalc_org, activity_water, activity_org, mass_fraction1, mass_fraction2,Gibbs_RT, dGibbs_RTdx2]...
        =BAT_properties_calculation_v1(mole_frac, O2C(i,1),  H2C(i,1), Mratio(i,1),BAT_functional_group(i,:),[]);
    
    % check for phase seperation
%     [~,phaseSep_via_activity_curvature_w,~,index_phase_sep_end_w]=...
%         finds_PhaseSep_and_activity_curve_dips_v2(activity_water);
%     [~,phaseSep_via_activity_curvature_org,~,index_phase_sep_end_org,]=...
%         finds_PhaseSep_and_activity_curve_dips_v2(activity_org);
%     a_w_sep_i=min(index_phase_sep_end_w,index_phase_sep_end_org);
    
    [phase_sep_check,~,upper_a_w_sep_index, ~]=finds_PhaseSep_w_and_org(activity_water,activity_org);

    
    if phase_sep_check==1
        RH_cross_point(i,1)=activity_water(upper_a_w_sep_index); % save phase sep RH
    else
        RH_cross_point(i,1)=0;
    end
end
% Checks outputs with in physical limits 
%round to zero
lower_round=RH_cross_point>0;
RH_cross_point=RH_cross_point.*lower_round+not(lower_round);%.*10^-10;
% round max to 1
upper_round=RH_cross_point>1;
RH_cross_point=not(upper_round).*RH_cross_point+upper_round;


end
