function O2C_single_phase_cross_point=single_phase_O2C_point_KGv3(molarmass_ratio);

x=molarmass_ratio;

%        p1 =      0.9483 ;% (0.5704, 1.326)
%        p2 =     -0.5388 ;% (-0.7592, -0.3185)
%        p3 =      0.2442 ;% (0.1098, 0.3787)
%        p4 =    0.004689 ;% (0.002557, 0.006822)
%        q1 =      0.4936 ;% (-0.131, 1.118)
%        q2 =      0.4185  ;%(0.2092, 0.6278)
%        q3 =    0.002564 ;% (0.0014, 0.003727)
% O2C_single_phase_cross_point = (p1.*x.^3 + p2.*x.^2 + p3.*x + p4)./(x.^3 + q1.*x.^2 + q2.*x + q3);


% new fit with OH sim Chain

%  Linear model Poly5:
%      f(x) = p1*x^5 + p2*x^4 + p3*x^3 + p4*x^2 + p5*x + p6
% Coefficients (with 95% confidence bounds):
%        p1 =      -444.7  (-606.9, -282.4)
%        p2 =       480.7  (326, 635.4)
%        p3 =      -181.1  (-235.8, -126.5)
%        p4 =       28.65  (19.9, 37.41)
%        p5 =      -2.462  (-3.081, -1.843)
%        p6 =        0.52  (0.5051, 0.5348)
% 
% Goodness of fit:
%   SSE: 0.0001161
%   R-square: 0.9991
%   Adjusted R-square: 0.9989
%   RMSE: 0.002073



% sigmoid eqiv.
O2C_single_phase_cross_point=.205./(1+exp(26.6.*(x-.12))).^.843+.225; % changed from 0.23 to 0.22

%O2C_single_phase_cross_point=0.1756./(1+exp(26.6.*(x-.12))).^.843+.225; %3-11-2019
