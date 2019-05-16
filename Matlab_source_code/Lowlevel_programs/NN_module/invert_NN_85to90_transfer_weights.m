function [transfer,weight_lower,weight_higher] = invert_NN_85to90_transfer_weights(a_w)

%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2018-Oct-30 10:08 AM
% Copyright 2018 Kyle Gorkowski 
%% 

upper_stop=.90;
lower_stop=.85;



if a_w<lower_stop;
    transfer=0;
    weight_lower=1;
    weight_higher=0;
elseif a_w>upper_stop;
        transfer=0;
    weight_lower=0;
    weight_higher=1;
else
    
normalize_x_value=(a_w-lower_stop)./(upper_stop-lower_stop);
    
zero_shift=0.5;
scale_factor=10;
x_span=0.5;
min_val=[0.00669285092428486];%1/(1+exp(-scale_factor.*((zero_shift-x_span)-zero_shift)));
max_val=[0.986614298151430];%1/(1+exp(-scale_factor.*((zero_shift+x_span)-zero_shift)))-min_val;

 scaled_sigmoid=(1.0./(1+exp(-scale_factor.*(normalize_x_value-zero_shift)))-min_val)./max_val;

         transfer=1;
    weight_lower=1-scaled_sigmoid;
    weight_higher=scaled_sigmoid;
 
end



end

