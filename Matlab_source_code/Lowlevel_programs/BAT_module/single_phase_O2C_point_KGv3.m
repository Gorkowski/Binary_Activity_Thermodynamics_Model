function O2C_single_phase_cross_point=single_phase_O2C_point_KGv3(molarmass_ratio);

% sigmoid eqiv.
O2C_single_phase_cross_point=.205./(1+exp(26.6.*(molarmass_ratio-.12))).^.843+.225; % changed from 0.23 to 0.22
