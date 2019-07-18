function [kappaHGF,kappa, growth] = bulk_kappa_v1(O2C_values, H2C_values,...
    Molecular_weight, aw_series, Coa_j_PM, Caq_PM, kappa_CCN_settings)
%%
% Created by Kyle Gorkowski [GORKOWFALCON] on 2018-Dec-16  1:38 PM
% Copyright 2018 Kyle Gorkowski 
%% 
if not(exist('kappa_CCN_settings'))
    % set default options
    kappa_CCN_settings=get_default_kappa_CCN_settings;
end

%% calc kappa
%Petters, M. D. and Kreidenweis, S. M.: A single parameter representation
%of hygroscopic growth and cloud condensation nucleus activity – Part 2:
%Including solubility, Atmos. Chem. Phys., 7(8), 1961–1971, doi:10.5194/acp-7-1961-2007, 2007.

[min_aw_val,min_aw_i]=min(aw_series);
% zero aw limit
trip_fix_0=false;
if min_aw_val==0 % replaces 0aw with small number to advoide singularity
    aw_series(min_aw_i)=10^-10;
    trip_fix_0=true;
end
S_full=size(aw_series);

C_OA_PM=sum(Coa_j_PM,2);

% get density
[densityEst]=Org_density_Estimate_KGv1(Molecular_weight, O2C_values, H2C_values);
densityEst_ugPm3=densityEst.*10^12; % covert from g/cm3
densityEst_water_ugPm3=0.997.*10^12;

% calc total mass
Coa_j_dry=repmat(Coa_j_PM(min_aw_i,:),S_full(1,1),1);

% mean density of OA mass weighted
Mass_fractions=Coa_j_PM./C_OA_PM;
Mass_fractions_dry=Coa_j_dry./sum(Coa_j_dry,2);

mean_densityEst_ugPm3=sum(repmat(densityEst_ugPm3',S_full(1,1),1).*Mass_fractions,2);
mean_densityEst_ugPm3_dry=sum(repmat(densityEst_ugPm3',S_full(1,1),1).*Mass_fractions_dry,2);

% convert to
Vs_dry=sum(Coa_j_dry,2)./mean_densityEst_ugPm3_dry;
Vs=sum(Coa_j_PM,2)./mean_densityEst_ugPm3;
Vw=sum(Caq_PM,2)./densityEst_water_ugPm3;

% bulk kappa methods
kappaHGF=(1./aw_series-1).*(Vw+Vs-Vs_dry)./Vs_dry; % accounts for co-condensation
kappa=(1./aw_series-1).*(Vw)./Vs; % with out co-condensation

growth.kappa.kappaHGF=kappaHGF;
growth.kappa.kappa=kappa;
growth.kappa.note=' kappa HGF accounts for co-condensation [kappaHGF=(1./aw_series-1).*(Vw+Vs-Vs_dry)./Vs_dry] and standard kappa does not [kappa=(1./aw_series-1).*(Vw)./Vs]'; 

%Organic volume growth factor
growth.Organic.volume=Vs./Vs_dry;
growth.Organic.radius=growth.Organic.volume.^(1/3);
growth.Organic.mass=C_OA_PM./C_OA_PM(min_aw_i);
growth.Organic.mean_density_ugPm3=mean_densityEst_ugPm3;

growth.PM.volume_HGF=(Vw+Vs-Vs_dry)./Vs_dry;
growth.PM.radius_HGF=growth.PM.volume_HGF.^(1/3);
growth.PM.mass_HGF=(C_OA_PM+Caq_PM)./(C_OA_PM(min_aw_i)+Caq_PM(min_aw_i));

growth.PM.volume=(Vw+Vs)./Vs_dry;
growth.PM.radius=growth.PM.volume.^(1/3);
growth.PM.mass=(C_OA_PM+Caq_PM)./(C_OA_PM(min_aw_i)+Caq_PM(min_aw_i));


%% Satruation calculation, direct from VBS-BAT with default size and surface tension mixing
mass_fraction_water=Caq_PM./(Caq_PM+C_OA_PM);
mass_fraction_org=C_OA_PM./(Caq_PM+C_OA_PM);

% calc saturation ratios
%% program
Diameter_vol_eqv_org_nm=kappa_CCN_settings.Diameter_vol_eqv_org_nm;
density_water_g_cm3=kappa_CCN_settings.density_water_g_cm3;
Mw=kappa_CCN_settings.Mw;
R=kappa_CCN_settings.R;
Temp=kappa_CCN_settings.Temp;
a_w=aw_series;


V_dp_m=4./3.*pi.*(Diameter_vol_eqv_org_nm.*10^-9./2).^3; % volume dp particle
w_org=mean_densityEst_ugPm3.*10^-12.*1000.*V_dp_m;  % estimate mass of org
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

[sat_max, sat_max_i]=max(SatRatio);



% [kappa_SatCritical, ~, ~, ~, ~, alternate]= kappa_critical_calc_for_McGlashan_v3(...
%     mean_densityEst_ugPm3.*10^-12, mass_fraction_water, mass_fraction_org, aw_series, aw_series, kappa_CCN_settings);

growth.CCN_direct.diamter_total=Diameter_total;
growth.CCN_direct.SatRatio=SatRatio;
growth.CCN_direct.kappa_with_cocondensation=kappaHGF(sat_max_i);
growth.CCN_direct.kappa_without_cocondensation=kappa(sat_max_i);
growth.CCN_direct.SatRatio_max=sat_max;


end

