function [data] = round_extremes_v1(data,lower_value, upper_value)
%%
% Created by Kyle Gorkowski [GORKOWFALCON] on 2019-Jan-19  3:38 PM
% Copyright 2019 Kyle Gorkowski 
%% 
lower_round=data>lower_value;
data=data.*lower_round+not(lower_round).*lower_value;
% round max to 1
upper_round=data>upper_value;
data=not(upper_round).*data+upper_round.*upper_value;

end

