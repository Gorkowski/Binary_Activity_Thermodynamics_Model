 function [avg_M_gPmol, molar_avg_O2C, molar_avg_H2C, molar_avg_N2C, mass_weighted_avg_O2C] = molar_based_means(C_ug_m3, M_gPmol, O2C, H2C, N2C)
%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2019-Aug-12  3:11 PM
% Copyright 2019 Kyle Gorkowski 
%% 
% calculated average properties base on molar abundance. 

% example
% clear
% 
% O2C=[.5;.8;1];
% M_gPmol=[180,200,350]';
% C_ug_m3=[[5,3,8];[5,3,3]];

% program
if not(exist('N2C'))
    N2C=M_gPmol.*0;
end
if not(exist('H2C'))
    H2C=M_gPmol.*0;
end

MC = 12.01D0; MO = 16.0D0; MH = 1.008D0; MN = 14.0067D0; %the molar masses of the carbon, oxygen and hydrogen atoms in [g/mol]

%1) estimate the H2C value if not provided from input
if (H2C < 0.1D0) 
    %estimate H2C assuming an aliphatic compound with H2C = 2 in absence of oxygen functional groups, then correct for oxygen content assuming a -1 slope (Van Krevelen Diagram of typical SOA).
    H2Cest = 2.0D0 -O2C;
else
    H2Cest = H2C;
end

%2) compute the approx. number of atoms per organic molecule
NC = M_gPmol./(MC + H2Cest.*MH + O2C.*MO + N2C.*MN); %carbon
NO = O2C.*NC;
NH = H2Cest.*NC;
NN = N2C.*NC;


% 3) molar concentration
shape_C_ug_m3=size(C_ug_m3);

molar_C_ug_m3=C_ug_m3./repmat(M_gPmol',shape_C_ug_m3(1,1),1);


% new ratios based on molar concentrations
molar_NC=sum(molar_C_ug_m3.* repmat(NC',shape_C_ug_m3(1,1),1),2);
molar_avg_O2C=sum(molar_C_ug_m3.* repmat(NO',shape_C_ug_m3(1,1),1),2)./molar_NC;
molar_avg_H2C=sum(molar_C_ug_m3.* repmat(NH',shape_C_ug_m3(1,1),1),2)./molar_NC;
molar_avg_N2C=sum(molar_C_ug_m3.* repmat(NN',shape_C_ug_m3(1,1),1),2)./molar_NC;


% avg molar weight
mass_fraction_inPM=C_ug_m3./sum(C_ug_m3,2);
avg_M_gPmol=sum(mass_fraction_inPM.* repmat(M_gPmol',shape_C_ug_m3(1,1),1),2);
mass_weighted_avg_O2C=sum(mass_fraction_inPM.* repmat(O2C',shape_C_ug_m3(1,1),1),2);



end

