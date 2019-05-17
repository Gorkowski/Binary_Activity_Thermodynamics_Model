function [VBSBAT_options]=default_VBSBAT_options(run_mode)
%%
% Created by Kyle Gorkowski [GORKOWFALCON] on 2018-Dec-16 10:38 AM
% Copyright 2018 Kyle Gorkowski
%%
% generates the options structure/type for VBS_BAT_simulation_vX

if not(exist('run_mode'))
    % set default options
    run_mode='default';
end

% q_alpha options
VBSBAT_options.q_alpha.min_spread_in_aw=10^-6;
VBSBAT_options.q_alpha.q_alpha_at_1phase_aw=0.99;%1-10^-6;

VBSBAT_options.q_alpha.q_alpha_bounds=[10^-6,1];
VBSBAT_options.q_alpha.q_alpha_bounds_mean=[10^-6,1];

VBSBAT_options.q_alphaVBS.Option_to_overwrite_two_phase_calculation='yes';
VBSBAT_options.q_alphaVBS.overwrite_threshold_two_phase_fraction_difference=0;
VBSBAT_options.q_alphaVBS.method_to_use='individual'; %
VBSBAT_options.q_alphaVBS.method_to_use_options={'mean_prop', 'individual'};
% VBSBAT_options.q_alphaVBS.two_phase_Cstar_average='mean';

% phase calc isolation
VBSBAT_options.force_phase.onePhase='no';
VBSBAT_options.force_phase.onePhase_options={'alpha','beta','no'};

% optimization options
VBSBAT_options.optimization.fit_tolerance=10^-5;
VBSBAT_options.optimization.guess_refinement_threshold=10^-5;
VBSBAT_options.optimization.global_NumTrialPoints=100;
VBSBAT_options.optimization.global_NumStageOnePoints=10;
VBSBAT_options.optimization.MaxIter=250;
VBSBAT_options.optimization.opt_method='fmincon'; % 'global' or 'fmincon' or 'none'
VBSBAT_options.optimization.independent_aw='yes';
% if the C OA mass vs aw doesn't look consitent at low RHs change
% fit_Algorithm to 'interior-point' as it is more robust
VBSBAT_options.optimization.fit_Algorithm='sqp'; % 'interior-point' (default) 'sqp' 'active-set'

% using a neural network for the inital guess is not computationally faster
% when using 'interior-point' fit_Algorithm above.
VBSBAT_options.VBSBAT_NN_options.use_NN_for_VBS_initial_guess='yes'; % yes or no
VBSBAT_options.VBSBAT_NN_options.NN_type='individual_properties'; % mean_properties individual_properties



%% BAT
VBSBAT_options.BAT_refinement_aw=0.90; % due to errors in NN fit
VBSBAT_options.BAT_refinement_tolerance=10^-6;

%% plots
VBSBAT_options.plot_PM='no';

VBSBAT_options.mean_BAT_functional_group='hydroxyl';
VBSBAT_options.BAT_functional_group_options={'hydroxyl';'carboxyl';'ether';'ketone';'ester';'hydroperoxide';'hydroperoxideSOA';'PEG'};

VBSBAT_options.run_mode_used=run_mode;
VBSBAT_options.run_mode_possible_options={{'default', 'robust', 'NN only', 'beta only'}; {' Add your own in default_VBSBAT_options.m lines 61-80'}};
if strcmpi(run_mode, 'default')
    % changes nothing
elseif strcmpi(run_mode, 'robust')
    VBSBAT_options.optimization.opt_method='fmincon'; % 'global' or 'fmincon' or 'none'
    VBSBAT_options.optimization.fit_Algorithm='interior-point'; % 'interior-point' (default) 'sqp' 'active-set'
    VBSBAT_options.BAT_refinement_aw=0;
    VBSBAT_options.optimization.fit_tolerance=10^-6;
    VBSBAT_options.optimization.guess_refinement_threshold=0;
    
elseif strcmpi(run_mode, 'NN only')
    
    VBSBAT_options.optimization.opt_method='none'; % 'global' or 'fmincon' or 'none'
    VBSBAT_options.VBSBAT_NN_options.use_NN_for_VBS_initial_guess='yes';
    
elseif strcmpi(run_mode, 'beta only')
    VBSBAT_options.force_phase.onePhase='beta';
        VBSBAT_options.BAT_refinement_aw=1;

else
    
    warning('no run_mode selected, so using default run_mode')

end



end


