  function [partition_coefficents, Coa_j_AB, Caq_j_AB, Cstar_j, Coa_AB, Caq_AB, Coa, q_alpha_water, error_out]= VBS_equilibration_withLLEpartition_objFun_KGv3(...
      guess_CoaEj, C_OM_ugPm3, Cstar_dry, activity_coefficient_AB, q_alpha_molefrac_phase_split_org, mass_fraction_water_AB, molecular_weight)
% % %%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2018-May-22  2:41 PM
% Copyright 2018 Kyle Gorkowski 
% Updated by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2018-May-30  2:20 PM
% rework of 2phase vbs with correction to C*
%% 

% %% test extremes
% clear
% guess_CoaEj=  [0.0269482661652176;0.356015485468693;0.161961397410055;0.999984266018315;0.344628565184500;0.992535073686157;0.998023977276474;0.999984160920101;0.124264815131966;0.999492046478144];;
% C_OM_ugPm3=[9.40227657770280;4.26067371914800;1.21313760354800;4.35217753885218;0.671761394400000;0.982360452692740;0.819566367988052;1.08804960692435;0.427227910730400;0.334933036705977];
% Cstar_dry=[268.759058670634;12.6541277225292;41.5753100838800;2.29192903909092e-07;13.1614064362187;0.0570212237834359;0.0142702364468685;2.73383721711946e-07;41.4409242242485;0.00377919930805804];
% activity_coefficient_AB=[54.9574380195832,1;26.7816625770900,1;97.4448787188580,1;8109.21338153931,1;23.6448809091187,1;1,1;1,1;8109.76843586726,1;9.31882845259539,1;1,1];
% q_alpha_molefrac_phase_split_org=[0,0,0,0,0,0,0,0,0,0]';
% mass_fraction_water_AB=[0.881777636853347,0;0.843664973477244,0;0.999999880034427,0;0.949288859923453,0;0.828120435885848,0;0,0;0,0;0.949278861472932,0;0.668541900636056,0;0,0];
% molecular_weight=[200.166000000000;188.174000000000;216.130000000000;368.298000000000;186.166000000000;204.182000000000;195.172000000000;368.306000000000;158.174000000000;206.138000000000];
% 

% guess_CoaAlpha_ugPm3=guess_CoaEj(1,1);
% guess_CoaBeta_ugPm3=guess_CoaEj(2,1);
% guess_q_alpha_water=guess_CoaEj(3,1);
guess_Ej=guess_CoaEj(1:end,1);


% get terms 
mass_fraction_water_alpha=mass_fraction_water_AB(:,1);
mass_fraction_water_beta=mass_fraction_water_AB(:,2);
activity_coefficient_alpha=activity_coefficient_AB(:,1);
activity_coefficient_beta=activity_coefficient_AB(:,2);

q_alpha=q_alpha_molefrac_phase_split_org;
q_beta=1-q_alpha_molefrac_phase_split_org;

% tiny_0=10^-32;
% tiny_almost1=1-10^-32;
% q_alpha=replace_data_A_to_B_KGv1(q_alpha, 0, tiny_0);
% q_alpha=replace_data_A_to_B_KGv1(q_alpha, 1, tiny_almost1);
% q_beta=replace_data_A_to_B_KGv1(q_beta, 0, tiny_0);
% q_beta=replace_data_A_to_B_KGv1(q_beta, 1, tiny_almost1);

%% checks that denomintor is not Inf or NaN
mass_fraction_water_alpha_denominator=(1-mass_fraction_water_alpha);
mass_fraction_water_beta_denominator=(1-mass_fraction_water_beta);
alpha_good_denominators=mass_fraction_water_alpha_denominator>0; % check points that have physical realistic values of mass fraciton water
beta_good_denominators=mass_fraction_water_beta_denominator>0;
% relpace bad points, x/(1-x) is ~x for x~=0... bellow set (1-x)==1
mass_fraction_water_alpha_denominator=abs(mass_fraction_water_alpha_denominator).*alpha_good_denominators+not(alpha_good_denominators);
mass_fraction_water_beta_denominator=abs(mass_fraction_water_beta_denominator).*beta_good_denominators+not(beta_good_denominators);

        
%% New Coa j
Coa_j=guess_Ej.*C_OM_ugPm3;
Coa_guess_viaEj=sum(Coa_j);

% alpha
Coa_j_alpha=Coa_j.*q_alpha;
Caq_alpha=sum(Coa_j_alpha.*mass_fraction_water_alpha./mass_fraction_water_alpha_denominator);
Coaaq_alpha=sum(Coa_j_alpha)+Caq_alpha;

% beta 
Coa_j_beta=Coa_j.*q_beta;
Caq_beta=sum(Coa_j_beta.*mass_fraction_water_beta./mass_fraction_water_beta_denominator);
Coaaq_beta=sum(Coa_j_beta)+Caq_beta;

Coaaq_via_Ej_guess=Coaaq_beta+Coaaq_alpha;

% C* via alpha phase
massweighted_molar_weight_alpha=sum(Coa_j_alpha./(molecular_weight))+Caq_alpha./18.015;

if massweighted_molar_weight_alpha>0
    Cstar_j_via_alpha=Cstar_dry.*activity_coefficient_alpha.*q_alpha.*Coaaq_via_Ej_guess./(molecular_weight.*massweighted_molar_weight_alpha); % calc new Cstar
else
    Cstar_j_via_alpha=q_alpha.*0;
end

% C* via beta phase
massweighted_molar_weight_beta=sum(Coa_j_beta./(molecular_weight))+Caq_beta./18.015;

if massweighted_molar_weight_beta>0
    Cstar_j_via_beta=Cstar_dry.*activity_coefficient_beta.*q_beta.*Coaaq_via_Ej_guess./(molecular_weight.*massweighted_molar_weight_beta); % calc new Cstar
else
    Cstar_j_via_beta=q_beta.*0;
end

%% select C* to use 
% could have different methods used here
% weight average by q_alpha
 Cstar_j_used=Cstar_j_via_alpha.*q_alpha+Cstar_j_via_beta.*q_beta;
% Cstar_j_used=Cstar_j_via_alpha.^q_alpha.*Cstar_j_via_beta.^q_beta;
% weighted by mass fraction in each phase
%Cstar_j_used=Cstar_j_via_alpha.*Coa_j_alpha./Coa_j+Cstar_j_via_beta.*Coa_j_beta./Coa_j;
% % use mean ov Cstar
% if massweighted_molar_weight_beta==0
%     Cstar_j_used=Cstar_j_via_alpha;
% elseif massweighted_molar_weight_alpha==0
%     Cstar_j_used=Cstar_j_via_beta;
% else
%     Cstar_j_used=(Cstar_j_via_alpha+Cstar_j_via_beta)./2;
% end

Ej_new=(1+Cstar_j_used./Coaaq_via_Ej_guess).^-1;

% calculate new mass values
Coa_j_new=Ej_new.*C_OM_ugPm3;
Coa_new_viaEj=sum(Coa_j);

% alpha
Coa_j_alpha_new=Coa_j.*q_alpha;
Caq_j_alpha_new=(Coa_j_alpha_new.*mass_fraction_water_alpha./mass_fraction_water_alpha_denominator);
Caq_alpha_new=sum(Caq_j_alpha_new);
Coa_alpha_new=sum(Coa_j_alpha_new);

% beta
Coa_j_beta_new=Coa_j_new.*q_beta;
Caq_j_beta_new=(Coa_j_beta_new.*mass_fraction_water_beta./mass_fraction_water_beta_denominator);
Caq_beta_new=sum(Caq_j_beta_new);
Coa_beta_new=sum(Coa_j_beta_new);

q_alpha_water_new=Caq_alpha_new./(Caq_alpha_new+Caq_beta_new);

% error for convergence of VBS
% error_out=(guess_CwOA_ugPm3-CwOA_out).^2+sum((Ei-Ei_out).^2);
error_out=sum((guess_Ej-Ej_new).^2+(Coa_guess_viaEj-Coa_new_viaEj).^2);


% collect the outputs
partition_coefficents=[Ej_new];
q_alpha_water=q_alpha_water_new;
Coa_j_AB=[Coa_j_alpha_new,Coa_j_beta_new];
Caq_j_AB=[Caq_j_alpha_new,Caq_j_beta_new];
Coa_AB=[Coa_alpha_new,Coa_beta_new]; 
Caq_AB=[Caq_alpha_new,Caq_beta_new];
Coa=sum(Coa_AB);
Cstar_j=Cstar_j_used;

end



