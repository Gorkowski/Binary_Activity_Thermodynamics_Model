function [shift_method] = check_BAT_functional_group_inputs_v1(O2C, shift_method)
%%
% Created by Kyle Gorkowski [GORKOWFALCON] on 2018-Dec-17  7:51 AM
% Copyright 2018 Kyle Gorkowski 
%% make the sizes of the O2C and shift method equal

[max_dim, di]=max(size(O2C));
Smethod=size(shift_method,di);
O2C_eqv=zeros(max_dim,1);
molarmass_ratio_eqv=zeros(max_dim,1);

if Smethod==1
    if iscell(shift_method)
        shift_method=repmat(shift_method,max_dim,1);
    else
        shift_method=repmat({shift_method},max_dim,1);
    end
elseif Smethod<max_dim
    error(['shift_method has less points than O2C ' num2str(Smethod) ' vs '  num2str(max_dim)])
end

end

