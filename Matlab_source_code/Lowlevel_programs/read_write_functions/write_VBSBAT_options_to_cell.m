function [options_2write]=write_VBSBAT_options_to_cell(VBSBAT_options);
%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2019-Mar-22  3:32 PM
% Copyright 2019 Kyle Gorkowski 
%% write cell of VBSBAT_options=default_VBSBAT_options('robust');%'robust' 'default'

% clear
% 
% VBSBAT_options=default_VBSBAT_options('default');%'robust' 'default'


options_2write=cell(24,1);
k=1;
options_2write(k,1)={['VBSBAT_options.run_mode_used='  VBSBAT_options.run_mode_used, ...
    ',{default,robust,NN only,beta only} ']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.BAT_refinement_aw='  num2str(VBSBAT_options.BAT_refinement_aw), ...
    ',{any water activity value above which refinement starts 0 to 1}']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.BAT_refinement_tolerance='  num2str(VBSBAT_options.BAT_refinement_tolerance), ...
    ',{not required in interpolation mode}']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.plot_PM=' VBSBAT_options.plot_PM ',{yes, no : option to make simple output graph}']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.mean_BAT_functional_group=' VBSBAT_options.mean_BAT_functional_group, ...
    ',{used in mean beta phase qalpha calc.: hydroxyl;carboxyl;ether;ketone;ester;hydroperoxide;hydroperoxideSOA;PEG}']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.q_alpha.org_1phase_shift_in_aw=' num2str(VBSBAT_options.q_alpha.org_1phase_shift_in_aw), ...
    ',{optional offset in a_w_sep}']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.q_alpha.min_spread_in_aw=' num2str(VBSBAT_options.q_alpha.min_spread_in_aw), ...
    ',{min delta in 1-a_w_sep}']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.q_alpha.scale_transfer_range=' num2str(VBSBAT_options.q_alpha.scale_transfer_range), ...
    ',{scale delta(a_w_sep) by this multiple}']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.q_alpha.q_alpha_bounds=' num2str(VBSBAT_options.q_alpha.q_alpha_bounds(1,1)), ',' num2str(VBSBAT_options.q_alpha.q_alpha_bounds_mean(1,2)), ...
    ',{bounds of q_alpha, above round to 1 and below round to 0}']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.q_alpha.q_alpha_bounds_mean=' num2str(VBSBAT_options.q_alpha.q_alpha_bounds_mean(1,1)), ',' num2str(VBSBAT_options.q_alpha.q_alpha_bounds_mean(1,2)), ...
    ',{bounds of mean beta phase q_alpha, above round to 1 and below round to 0}']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.q_alphaVBS.method_to_use=' VBSBAT_options.q_alphaVBS.method_to_use, ...
    ',{q_alpha method use in q_alpha calc.: mean_prop,individual}']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.q_alphaVBS.Option_to_overwrite_two_phase_calculation=' VBSBAT_options.q_alphaVBS.Option_to_overwrite_two_phase_calculation, ...
    ',{check for non-physical changes at q_alpha transitions: yes, no}']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.q_alphaVBS.overwrite_threshold_two_phase_fraction_difference=' num2str(VBSBAT_options.q_alphaVBS.overwrite_threshold_two_phase_fraction_difference), ...
    ',{fractional difference needed to trigger an overwrite}']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.force_phase.onePhase=' VBSBAT_options.force_phase.onePhase, ...
    ',{force single phase calc: no, beta, alpha }']};

% VBS fitting options
k=k+1;
options_2write(k,1)={['VBSBAT_options.optimization.opt_method=' VBSBAT_options.optimization.opt_method, ...
    ',{VBS equil. solver: none, fmincon, global }']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.optimization.fit_Algorithm=' VBSBAT_options.optimization.fit_Algorithm, ...
    ',{VBS equil. algorithm: sqp, interior-point, trust-region-reflective, active-set}']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.optimization.fit_tolerance=' num2str(VBSBAT_options.optimization.fit_tolerance), ...
    ',{solver exit tolerance}']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.optimization.global_NumTrialPoints=' num2str(VBSBAT_options.optimization.global_NumTrialPoints), ...
    ',{total trial points in global solver}']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.optimization.global_NumStageOnePoints=' num2str(VBSBAT_options.optimization.global_NumStageOnePoints), ...
    ',{stage one trial points in global solver}']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.optimization.MaxIter=' num2str(VBSBAT_options.optimization.MaxIter), ...
    ',{max iterations in solver}']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.optimization.independent_aw=' VBSBAT_options.optimization.independent_aw, ...
    ',{treat each a_w independently}']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.optimization.guess_refinement_threshold=' num2str(VBSBAT_options.optimization.guess_refinement_threshold), ...
    ',{if the error in the guess is below, then no refinement is done}']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.VBSBAT_NN_options.use_NN_for_VBS_initial_guess=' VBSBAT_options.VBSBAT_NN_options.use_NN_for_VBS_initial_guess, ...
    ',{use NN to get initial guess}']};
k=k+1;
options_2write(k,1)={['VBSBAT_options.VBSBAT_NN_options.NN_type=' VBSBAT_options.VBSBAT_NN_options.NN_type, ...
    ',{NN options: mean_properties, individual_properties }']};




