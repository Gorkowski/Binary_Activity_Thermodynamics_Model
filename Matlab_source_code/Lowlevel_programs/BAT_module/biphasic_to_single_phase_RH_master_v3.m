function [RH_cross_point] = biphasic_to_single_phase_RH_master_v3(O2C, H2C, Mratio, mode)
    error('change to v4')

% Updated by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2018-Oct-26  3:38 PM
% Updated by Kyle Gorkowski [GORKOWFALCON] on 2018-Dec-16 10:38 AM

[O2C,Mratio] = convert_chemical_structure_to_OH_eqv_v3(O2C, Mratio, mode);

x1=[O2C, Mratio];
RH_cross_point=O2C.*0;

[org_mole_fraction] = biphasic_to_single_phase_molfrac_org_NN_v5(x1')';

for i=1:length(O2C)
        % calcs activty

    [~, ~, ~, ~, RH_cross_point(i,1), ~, ~, ~, ~, ~]=BAT_properties_calculation_v1(...
        org_mole_fraction(i,1), O2C(i,1), H2C(i,1), Mratio(i,1), mode,[]); 
end
% Checks outputs with in physical limits 
%round to zero
lower_round=RH_cross_point>0;
RH_cross_point=RH_cross_point.*lower_round+not(lower_round);%.*10^-10;
% round max to 1
upper_round=RH_cross_point>1;
RH_cross_point=not(upper_round).*RH_cross_point+upper_round;


end
