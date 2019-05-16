function [phase_sep_check,lower_a_w_sep_index,upper_a_w_sep_index, matching_Upper_a_w_sep_index]=finds_PhaseSep_w_and_org(activity_water,activity_org);

%%
% Created by Kyle Gorkowski [GORKOWFALCON] on 2019-May-04  2:42 PM
% Copyright 2019 Kyle Gorkowski 
%% 


% check for phase seperation in each activity curve
[~,phaseSep_via_activity_curvature_w,index_phase_sep_starts_w,index_phase_sep_end_w]=...
    finds_PhaseSep_and_activity_curve_dips_v2(activity_water);
[~,phaseSep_via_activity_curvature_org,index_phase_sep_starts_org,index_phase_sep_end_org]=...
    finds_PhaseSep_and_activity_curve_dips_v2(activity_org);

indexes=[index_phase_sep_starts_w,index_phase_sep_end_w,index_phase_sep_starts_org,index_phase_sep_end_org]; % get start and end indexes
if phaseSep_via_activity_curvature_w==1
    phase_sep_check=1;
    if activity_water(1)<activity_water(end) % increasing a_w with index
        lower_a_w_sep_index=min(indexes);
        upper_a_w_sep_index=max(indexes);
        % find matching a_w point for upper_a_w_sep_index, as this is not
        % necessarily the same as lower_a_w_sep_index, and is likely on a metastable
        % phase-separation branch in Gibbs mix.
        mid_sep_index=floor((lower_a_w_sep_index+upper_a_w_sep_index)/2);
        activity_water_beta=activity_water(1:mid_sep_index);
        match_a_w=activity_water(upper_a_w_sep_index);
        [match_index_prime]=find((activity_water_beta-match_a_w)>0);
        if isempty(match_index_prime)
           [~, match_index_prime]=max(activity_water_beta-match_a_w);
        end
        matching_Upper_a_w_sep_index=match_index_prime(1,1)-1;
    else
        lower_a_w_sep_index=max(indexes);
        upper_a_w_sep_index=min(indexes);
        % find matching a_w point for upper_a_w_sep_index, as this is not
        % necessarily the same as lower_a_w_sep_index, and is likely on a metastable
        % phase-separation branch in Gibbs mix.
        mid_sep_index=floor((lower_a_w_sep_index+upper_a_w_sep_index)/2);
        activity_water_beta=activity_water(mid_sep_index:end);
        match_a_w=activity_water(upper_a_w_sep_index);
        match_index_prime=find(activity_water_beta<=match_a_w);
        matching_Upper_a_w_sep_index=mid_sep_index+match_index_prime(1,1)-1;
    end
    


    
else
    lower_a_w_sep_index=1; % no phase sep
    upper_a_w_sep_index=2;
    matching_Upper_a_w_sep_index=2;
    phase_sep_check=0;% no phase sep
end


end