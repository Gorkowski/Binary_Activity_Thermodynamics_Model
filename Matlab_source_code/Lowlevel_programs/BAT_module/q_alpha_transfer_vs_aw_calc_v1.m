function [q_alpha_vsRH_values] = q_alpha_transfer_vs_aw_calc_v1(RH_cross_point,aw_series,VBSBAT_options)
%%
% Created by Kyle Gorkowski [GORKOWFALCON] on 2018-Dec-16  8:43 AM
% Copyright 2018 Kyle Gorkowski 
%% makes a squeezed logistic funciton to transfer for q_alpha ~0 to q_alpha ~1, desicribed in VBSBAT_options.q_alpha

    S=size(RH_cross_point);
    S_full=size(aw_series);
    
    % spread in transfer from 50/50 point
    RH_delta_to_one=1-RH_cross_point;
    org_1phase_spread_in_aw=min(Removes_zero_rows(RH_delta_to_one))+VBSBAT_options.q_alpha.org_1phase_shift_in_aw;
    if isempty(org_1phase_spread_in_aw)
       org_1phase_spread_in_aw= VBSBAT_options.q_alpha.min_spread_in_aw;
    end
        
    min_spread_points=RH_delta_to_one<=VBSBAT_options.q_alpha.min_spread_in_aw; % very high aw transistion to 1 phase
    RH_delta_to_one=RH_delta_to_one+min_spread_points.*VBSBAT_options.q_alpha.min_spread_in_aw;
    max_spread_points=RH_delta_to_one>=1; % no transition
%     RH_delta_to_one=RH_delta_to_one+max_spread_points.*(org_1phase_spread_in_aw-1);
    RH_delta_to_one=RH_delta_to_one+max_spread_points.*(0);

    RH_delta_to_one=VBSBAT_options.q_alpha.scale_transfer_range.*RH_delta_to_one;
    
    q_alpha_transfer_constant=log(1/(1-VBSBAT_options.q_alpha.q_alpha_at_1phase_aw)-1)./RH_delta_to_one; % squeeze the transfer function
    
    % aw at q_alpha 50/50 point
    min_1phase_RH=min(Removes_zero_rows(RH_cross_point))-org_1phase_spread_in_aw;
    
    
    q_alpha_RH_transfer_center=RH_cross_point+max_spread_points.*min_1phase_RH;
    q_alpha_RH_transfer_center=q_alpha_RH_transfer_center-RH_delta_to_one;
    
    % calc q_alpha vs aw
    q_alpha_vsRH_values=1-1./(1+exp((-repmat(q_alpha_RH_transfer_center',S_full(1,1),1)+repmat(aw_series,1,S(1,1) )).*repmat(q_alpha_transfer_constant',S_full(1,1),1)));

    
% t=1;

end

