function [kappaHGF,kappa, growth] = bulk_kappa_v1(O2C_values, H2C_values,...
    Molecular_weight, aw_series, Coa_j_PM, Caq_PM)
%%
% Created by Kyle Gorkowski [GORKOWFALCON] on 2018-Dec-16  1:38 PM
% Copyright 2018 Kyle Gorkowski 
%% 


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
kappa=(1./aw_series-1).*(Vw)./Vs; % with co-condensation

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


%% kappa CCN direct with defaults
mass_fraction_water=Caq_PM./(Caq_PM+C_OA_PM);
mass_fraction_org=C_OA_PM./(Caq_PM+C_OA_PM);
kappa_CCN_settings=get_default_kappa_CCN_settings;

[~, ~, ~, ~, ~, alternate]= kappa_critical_calc_for_McGlashan_v3(...
    mean_densityEst_ugPm3.*10^-12, mass_fraction_water, mass_fraction_org, aw_series, aw_series, kappa_CCN_settings);

growth.kappaCCN_direct=alternate;

end

