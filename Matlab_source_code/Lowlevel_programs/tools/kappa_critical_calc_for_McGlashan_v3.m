function [kappa_SatCritical, kappa_a_w, SatCritical, Dp_critical, V_w_critical, alternate]= kappa_critical_calc_for_McGlashan_v3(...
    density_org_g_cm3, mass_fraction_water, mass_fraction_org, activity_water, activity_org, kappa_CCN_settings)

%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2019-Jan-17 10:25 AM
% Copyright 2019 Kyle Gorkowski
% added beta and alpha phase distinction

if not(exist('kappa_CCN_settings'))
    % set default options
    kappa_CCN_settings=get_default_kappa_CCN_settings;
end
%% program
Diameter_vol_eqv_org_nm=kappa_CCN_settings.Diameter_vol_eqv_org_nm;
density_water_g_cm3=kappa_CCN_settings.density_water_g_cm3;
Mw=kappa_CCN_settings.Mw;
R=kappa_CCN_settings.R;
Temp=kappa_CCN_settings.Temp;
a_w=activity_water;

% check to flip data so aw is increaseing
if min(a_w)==length(a_w)
    mass_fraction_water=flip(mass_fraction_water);
    mass_fraction_org=flip(mass_fraction_org);
    a_w=flip(a_w);
end



V_dp_m=4./3.*pi.*(Diameter_vol_eqv_org_nm.*10^-9./2).^3; % volume dp particle
w_org=density_org_g_cm3.*1000.*V_dp_m;  % estimate mass of org

V_w=1./(density_water_g_cm3.*1000).*w_org.*(mass_fraction_water./mass_fraction_org);

V_total=V_dp_m+V_w;

if strcmpi(kappa_CCN_settings.surface_tension_method, 'volume mix')
    sigma_droplet_N_m=(kappa_CCN_settings.sigmaW_droplet_N_m.*V_w+kappa_CCN_settings.sigmaOrg_droplet_N_m.*V_dp_m)./V_total;
elseif strcmpi(kappa_CCN_settings.surface_tension_method, 'water')
    sigma_droplet_N_m=kappa_CCN_settings.sigmaW_droplet_N_m;
elseif strcmpi(kappa_CCN_settings.surface_tension_method, 'organic')
    sigma_droplet_N_m=kappa_CCN_settings.sigmaOrg_droplet_N_m;
end


Diameter_total=real(2*(V_total.*3/(4*pi)).^(1/3));

SatRatio=a_w.*exp_wlimiter(4.*sigma_droplet_N_m.*Mw./(R.*Temp.*density_water_g_cm3.*1000.*Diameter_total));
SatRatio_percent=(SatRatio-1).*100;

kappa=(1./a_w-1).*(V_w./V_dp_m);

% check for phase seperation
% [phaseSep_via_activity,phaseSep_via_activity_curvature,index_phase_sep_starts,index_phase_sep_end]=finds_PhaseSep_and_activity_curve_dips_v2(a_w);
[phase_sep_check,~,index_phase_sep_end, index_phase_sep_starts]=finds_PhaseSep_w_and_org(activity_water,activity_org);




fixed_aw_val=kappa_CCN_settings.fixed_aw_kappaCCN;
alternate.beta.phase=0;
alternate.fixed_aw.fixed_aw_val=fixed_aw_val; %save
alternate.fixed_aw.beta.phase=0;

if phase_sep_check==0 % no phase sep.
    [SatCritical,Sc_index]=max(SatRatio);
    
    kappa_SatCritical=kappa(Sc_index,1);
    kappa_a_w=a_w(Sc_index,1);
    
    Dp_critical=Diameter_total(Sc_index,1);
    V_w_critical=V_w(Sc_index,1);
    
    % kappa at fixed aw alpha phase
    [~,aw_index]=min(abs(a_w-fixed_aw_val));
    alternate.fixed_aw.alpha.aw=a_w(aw_index,1);
    
    alternate.fixed_aw.alpha.SatCritical=SatRatio(aw_index,1);
    
    alternate.fixed_aw.alpha.kappa_SatCritical=kappa(aw_index,1);
    alternate.fixed_aw.alpha.kappa_a_w=a_w(aw_index,1);
    
    alternate.fixed_aw.alpha.Dp_critical=Diameter_total(aw_index,1);
    alternate.fixed_aw.alpha.V_w_critical=V_w(aw_index,1);
    
    %save alpha
    alternate.alpha.kappa_SatCritical=kappa_SatCritical;
    alternate.alpha.kappa_a_w=kappa_a_w;
    alternate.alpha.SatCritical=SatCritical;
    alternate.alpha.Dp_critical=Dp_critical;
    alternate.alpha.V_w_critical=V_w_critical;
    
else % has phase sep so pulls kappa at onset of phase sep.
    %% fixed aw
    % kappa at fixed aw alpha phase
    [~,aw_index]=min(abs(a_w(1:index_phase_sep_starts,:)-fixed_aw_val));
    alternate.fixed_aw.beta.aw=a_w(aw_index,1);
    
    alternate.fixed_aw.beta.SatCritical=SatRatio(aw_index,1);
    
    alternate.fixed_aw.beta.kappa_SatCritical=kappa(aw_index,1);
    alternate.fixed_aw.beta.kappa_a_w=a_w(aw_index,1);
    
    alternate.fixed_aw.beta.Dp_critical=Diameter_total(aw_index,1);
    alternate.fixed_aw.beta.V_w_critical=V_w(aw_index,1);
    
    % kappa at fixed aw beta phase
    alternate.fixed_aw.beta.phase=1;
    
    [~,aw_index]=min(abs(a_w(index_phase_sep_end:end,:)-fixed_aw_val));
    aw_index=aw_index+index_phase_sep_end-1;
    alternate.fixed_aw.alpha.aw=a_w(aw_index,1);
    
    alternate.fixed_aw.alpha.SatCritical=SatRatio(aw_index,1);
    
    alternate.fixed_aw.alpha.kappa_SatCritical=kappa(aw_index,1);
    alternate.fixed_aw.alpha.kappa_a_w=a_w(aw_index,1);
    
    alternate.fixed_aw.alpha.Dp_critical=Diameter_total(aw_index,1);
    alternate.fixed_aw.alpha.V_w_critical=V_w(aw_index,1);
    
    %% beta
    alternate.beta.phase=1;
    
    
    % kappa beta
    Sc_index=index_phase_sep_starts;
    
    alternate.beta.SatCritical=SatRatio(Sc_index,1);
    
    alternate.beta.kappa_SatCritical=kappa(Sc_index,1);
    alternate.beta.kappa_a_w=a_w(Sc_index,1);
    
    alternate.beta.Dp_critical=Diameter_total(Sc_index,1);
    alternate.beta.V_w_critical=V_w(Sc_index,1);
    
    % alpha
    [~,Sc_index]=max(SatRatio(index_phase_sep_end:end,:));
    Sc_index=Sc_index+index_phase_sep_end-1;
    alternate.alpha.SatCritical=SatRatio(Sc_index,1);
    
    alternate.alpha.kappa_SatCritical=kappa(Sc_index,1);
    alternate.alpha.kappa_a_w=a_w(Sc_index,1);
    
    alternate.alpha.Dp_critical=Diameter_total(Sc_index,1);
    alternate.alpha.V_w_critical=V_w(Sc_index,1);
    
    
    
    %output uses beta
    kappa_SatCritical=alternate.beta.kappa_SatCritical;
    kappa_a_w=alternate.beta.kappa_a_w;
    SatCritical=alternate.beta.SatCritical;
    Dp_critical=alternate.beta.Dp_critical;
    V_w_critical=alternate.beta.V_w_critical;
    
end




alternate.settings=kappa_CCN_settings;


