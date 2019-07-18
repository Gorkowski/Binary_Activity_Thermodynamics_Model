function [kappa_CCN_settings] = get_default_kappa_CCN_settings
%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2019-Jan-18 11:14 AM
% Copyright 2019 Kyle Gorkowski 
%% 

kappa_CCN_settings.Diameter_vol_eqv_org_nm=100;
% kappa_CCN_settings.sigma_droplet_N_m=0.072;
kappa_CCN_settings.density_water_g_cm3=1;
kappa_CCN_settings.Mw=18.016.*10^-3;
kappa_CCN_settings.R=8.31446;
kappa_CCN_settings.Temp=293;

kappa_CCN_settings.sigmaW_droplet_N_m=0.072;
kappa_CCN_settings.sigmaOrg_droplet_N_m=0.03;
kappa_CCN_settings.surface_tension_method='volume mix'; % water, organic, volume mix
kappa_CCN_settings.surface_tension_method_options={'water', 'organic', 'volume mix'};
% fixed aw kappa
kappa_CCN_settings.fixed_aw_kappaCCN=0.99;

end

