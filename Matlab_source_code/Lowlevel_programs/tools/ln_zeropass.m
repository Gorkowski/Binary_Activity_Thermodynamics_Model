function [output] = ln_zeropass(input)
%%
% Created by Kyle Gorkowski [GORKOWFALCON] on 2019-Jan-20  4:46 PM
% Copyright 2019 Kyle Gorkowski 
%% 


input=replace_data_A_to_B_KGv1(input, 0,NaN);
output=log(input);
output=replace_data_A_to_B_KGv1(output, NaN,0);

end

