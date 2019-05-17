function [q_alpha_value] = q_alpha_transfer_vs_aw_calc_v1(a_w_sep,aw_series,VBSBAT_options)
%%
% Created by Kyle Gorkowski [GORKOWFALCON] on 2018-Dec-16  8:43 AM
% Copyright 2018 Kyle Gorkowski
%% makes a squeezed logistic funciton to transfer for q_alpha ~0 to q_alpha ~1, desicribed in VBSBAT_options.q_alpha

S=size(a_w_sep);
S_full=size(aw_series);

mask_of_miscible_points=a_w_sep==0; % values held for correction at the end

% spread in transfer from 50/50 point
delta_a_w_sep=1-a_w_sep;

% check min value allowed
above_min_delta_a_w_sep_value=delta_a_w_sep>VBSBAT_options.q_alpha.min_spread_in_aw;
delta_a_w_sep=delta_a_w_sep.*above_min_delta_a_w_sep_value+...
    not(above_min_delta_a_w_sep_value).*VBSBAT_options.q_alpha.min_spread_in_aw;

% calculate curve parameter of simgmoid
sigmoid_curve_parameter=ln_zeropass(1./(1-VBSBAT_options.q_alpha.q_alpha_at_1phase_aw)-1)./delta_a_w_sep;

% calculate q_alpha value
q_alpha_value=1-1./(1+exp_wlimiter(sigmoid_curve_parameter.*(aw_series-a_w_sep+delta_a_w_sep)));


% apply mask for complete misciblity, turns misicible organics to
% q_alpha=1 for all a_w
q_alpha_value=q_alpha_value.*not(mask_of_miscible_points)+mask_of_miscible_points;


end