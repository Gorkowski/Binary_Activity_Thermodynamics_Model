function [data_header, data_matrix, input_header, input_matrix, input_functional_groups, input_aw]=gen_simple_output_data_from_VBSBAT_sim(sim_details);
%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2019-Mar-25 10:06 AM
% Copyright 2019 Kyle Gorkowski 
%% simple output data


% clear
% load('VBSBAT_calc.mat','details_aP')
% 
% name='ap_simulation';
% sim_details=details_aP;

Saw=length(sim_details.inputs.aw_series);

data_matrix=zeros(Saw,13);
data_header=cell(13,1);

k=1;
data_header(k,1)={'water activity'};
data_matrix(:,k)=sim_details.inputs.aw_series;
k=k+1;
data_header(k,1)={'C_water_PM (ug/m3)'};
data_matrix(:,k)=sim_details.totals.Caq_PM;
k=k+1;
data_header(k,1)={'Calpha_water_PM (ug/m3)'};
data_matrix(:,k)=sim_details.totals.Caq_PM_alpha;
k=k+1;
data_header(k,1)={'Cbeta_water_PM (ug/m3)'};
data_matrix(:,k)=sim_details.totals.Caq_PM_beta;
k=k+1;
data_header(k,1)={'C_OA_PM (ug/m3)'};
data_matrix(:,k)=sim_details.totals.C_OA_PM;
k=k+1;
data_header(k,1)={'Calpha_OA_PM (ug/m3)'};
data_matrix(:,k)=sim_details.totals.Coa_alpha;
k=k+1;
data_header(k,1)={'Cbeta_OA_PM (ug/m3)'};
data_matrix(:,k)=sim_details.totals.Coa_beta;
k=k+1;
data_header(k,1)={'PM_VGF'};
data_matrix(:,k)=sim_details.growth.PM.volume;
k=k+1;
data_header(k,1)={'PM_DGF'};
data_matrix(:,k)=sim_details.growth.PM.radius;
k=k+1;
data_header(k,1)={'PM_MGF'};
data_matrix(:,k)=sim_details.growth.PM.mass;
k=k+1;
data_header(k,1)={'PM_OA_VGF'};
data_matrix(:,k)=sim_details.growth.Organic.volume;
k=k+1;
data_header(k,1)={'PM_OA_DGF'};
data_matrix(:,k)=sim_details.growth.Organic.radius;
k=k+1;
data_header(k,1)={'PM_OA_MGF'};
data_matrix(:,k)=sim_details.growth.Organic.mass;
k=k+1;
data_header(k,1)={'kappa_bulk'};
data_matrix(:,k)=sim_details.growth.kappa.kappa;
k=k+1;
data_header(k,1)={'kappa_HGF_bulk'};
data_matrix(:,k)=sim_details.growth.kappa.kappaHGF;
k=k+1;
data_header(k,1)={'OA_mean_density(ug/m3)'};
data_matrix(:,k)=sim_details.growth.Organic.mean_density_ugPm3;


%% system input prop
k=0;
k=k+1;
input_header(k,1)={'M (g/mol)'};
input_matrix(:,k)=sim_details.inputs.Molecular_weight;
k=k+1;
input_header(k,1)={'O:C'};
input_matrix(:,k)=sim_details.inputs.O2C_values;
k=k+1;
input_header(k,1)={'H:C'};
input_matrix(:,k)=sim_details.inputs.H2C_values;
k=k+1;
input_header(k,1)={'eff. Csat_j (ug/m3)'};
input_matrix(:,k)=sim_details.inputs.Csat_j_value;
k=k+1;
input_header(k,1)={'Ctotal_j (ug/m3)'};
input_matrix(:,k)=sim_details.inputs.C_OM_ugPm3;

input_functional_groups=sim_details.inputs.BAT_functional_group;
input_aw=sim_details.inputs.aw_series;











