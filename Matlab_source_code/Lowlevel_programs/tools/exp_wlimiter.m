function [y] = exp_wlimiter(x)
%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2018-Oct-31  4:58 PM
% Copyright 2018 Kyle Gorkowski 
%% 
ulimiter=690.7755;%log(10^300);
llimiter=-690.7755;%log(10^-300);

% force upper limit
upperlimit=x>ulimiter;
x=x.*not(upperlimit)+upperlimit.*ulimiter;

% force lower limit
lowerlimit=x<llimiter;
x=x.*not(lowerlimit)+lowerlimit.*llimiter;


y=exp(x);

end

