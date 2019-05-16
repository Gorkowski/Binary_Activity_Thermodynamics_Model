function [partition_coefficent_estimate] =  VBSBAT_nerual_network_estimate_v1(...
    Cj_sat_ugPm3, Cj_ugPm3, O2C, molecular_weight_gPmol, mass_fraction_water_beta ,a_water, VBSBAT_NN_options)
%%
% Created by Kyle Gorkowski [GORKOWFALCON] on 2019-Jan-19  2:01 PM
% Copyright 2019 Kyle Gorkowski 
%% one data set at a time

% example
% clear
% Cj_sat_ugPm3=10.^(rand(1,11).*diff([-10,2])-6)';
% Cj_ugPm3=[0.00728379183006649,0.00126729755041757,3.58418117500636e-05,0.00198483412363864,0.000698811733797980,0.00217387638180601,0.00515682007180860,0.00236906334566364,0.00139239290447069,0.00517271660215104,0.00187109754882274]';
% O2C=[0.749846747723861,0.513221697853408,1.03257931967642,0.504915917234067,1.19123018306739,0.403169421013742,0.290684071535326,0.636224306583813,0.407313104603885,0.879798318676640,0.412229873520999]';
% molecular_weight_gPmol=[576.578927511129,200.891273298969,558.718754454970,326.220430769759,511.207564680670,193.603693811755,163.951798522873,440.790633917725,414.712310110913,403.197257011539,209.965783352877]';
% mass_fraction_water_beta =[0.178805971613121,0.121823615843382,0.190549070223948,0.122145069616407,0.189092949464273,0.0911283724125909,0.0528819270150390,0.154222994422637,0.0999865996220898,0.181684410401761,0.0944107467283075]';
% a_water=[0.878329135065060]';


if not(exist('VBSBAT_NN_options'))
    % set default options
    VBSBAT_NN_options.NN_type='individual_properties'; %'mean_properties' 'individual_properties'
end

% blank output
Sc=size(Cj_ugPm3);
partition_coefficent_estimate=zeros(Sc);

% bining of inputs 
Csat_j_bins=[0;10.^(0.5+[-6:1:4-.5]');Inf]; % shift to bin edges and extend to full numberline 0 to inf
index_to_bin = discretize(Cj_sat_ugPm3,Csat_j_bins); % bin indexes

bins_to_use=unique(index_to_bin);
indexes=[1:length(Cj_sat_ugPm3)]';

input_Cj=zeros(length(Csat_j_bins)-1,1);
input_O2C=input_Cj;
input_M=input_Cj;
input_mf_water=input_Cj;

for i=1:length(bins_to_use) % sort data into bins
    
    input_i=bins_to_use(i);
    k = find(index_to_bin==input_i);
    
    input_Cj(input_i)=sum(Cj_ugPm3(k));
    input_O2C(input_i)=mean(O2C(k));
    input_M(input_i)=mean(molecular_weight_gPmol(k));
    input_mf_water(input_i)=mean(mass_fraction_water_beta(k));

end

%% nerual network calc         
if strcmpi(VBSBAT_NN_options.NN_type, 'mean_properties')
    input=[input_Cj; mean(input_O2C); mean(input_M); input_mf_water ;a_water(1,1)]; % could of only used 1 input for water mass fraction, may update NN later
    
    pc_temp=NN_VBSBAT_layers24_meanProp_v1(input);
    
elseif strcmpi(VBSBAT_NN_options.NN_type, 'individual_properties')
    

    
    input=[input_Cj; input_O2C; input_M; input_mf_water ;a_water(1,1)];
    [pc_temp] = NN_VBSBAT_layers28_individProp_v1(input);
    
%         input=[log10_zeropass(input_Cj); input_O2C; input_M; input_mf_water ;a_water(1,1)];
% 
%     [pc_temp] = NN_VBSBAT_layers26_individProp_logC_v2(input);
%     
elseif strcmpi(VBSBAT_NN_options.NN_type, 'individual_properties_v2')

    O2C_miscibility_line=single_phase_O2C_point_KGv3(18.016./input_M);
    O2C_miscibility_delta=input_O2C-O2C_miscibility_line;

    input_M=replace_data_A_to_B_KGv1(input_M,0,18.016./2);
    
        input=[log10_zeropass(input_Cj);O2C_miscibility_delta;18.016./input_M;input_mf_water;a_water(1,1)];

        [pc_temp] = NN_VBSBAT_layers20_individProp_v2(input);
else
    error('select VBSBAT_NN_options.NN_type')
end
[pc_temp] = round_extremes_v1(pc_temp,0, 1);

% resort to match inputs

partition_coefficent_estimate=Cj_sat_ugPm3.*0;
for i=1:length(bins_to_use) % sort data into bins
   
        input_i=bins_to_use(i);
    k = find(index_to_bin==input_i);

    
    partition_coefficent_estimate(k,:)=pc_temp(input_i,1);
end



% 
% 
% 
% end

