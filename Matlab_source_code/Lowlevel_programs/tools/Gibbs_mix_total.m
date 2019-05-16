function [ideal_Gibbs_RT,ideal_dGibbs_dxRT, Gibbs_total_RT,dGibbs_total_dxRT] = Gibbs_mix_total(org_mole_frac_scan, GibbsEx_RT, dGibbsEx_RTdx)
%%
% Created by Kyle Gorkowski [GORKOWFALCON] on 2019-May-05 10:16 AM
% Copyright 2019 Kyle Gorkowski 
% calculates the total Gibbs energy of mixing, non-ideal from BAT plus
% ideal part from mole fractions

ideal_Gibbs_RT=org_mole_frac_scan.*ln_zeropass(org_mole_frac_scan)+(1-org_mole_frac_scan).*ln_zeropass(1-org_mole_frac_scan);
ideal_dGibbs_dxRT=ln_zeropass(org_mole_frac_scan)-ln_zeropass(1-org_mole_frac_scan);

Gibbs_total_RT=GibbsEx_RT+ideal_Gibbs_RT;
dGibbs_total_dxRT=dGibbsEx_RTdx+ideal_dGibbs_dxRT;


end

