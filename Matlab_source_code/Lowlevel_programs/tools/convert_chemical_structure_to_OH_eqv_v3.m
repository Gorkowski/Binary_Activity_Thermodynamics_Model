function [O2C_eqv,molarmass_ratio_eqv] = convert_chemical_structure_to_OH_eqv_v3(O2C, molarmass_ratio, shift_method )
%%
% Created by Kyle Gorkowski [GORKOWFALCON] on 2018-Oct-27 10:33 AM
% Updated by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2018-Nov-21  3:35 PM
% Updated by Kyle Gorkowski [GORKOWFALCON] on 2018-Dec-17  6:46 AM
%  % now includes independent species, if shift method is only 1 group that
%  is applied to all.
% Copyright 2018 Kyle Gorkowski
%%
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

for sp_i=1:max_dim
    calc_shift=true;
    if  strcmpi('hydroxyl', shift_method(sp_i)) || strcmpi('carboxyl', shift_method(sp_i))
        % no change
        calc_shift=false;
        
    elseif strcmpi('hydroperoxideSOA', shift_method(sp_i)) || strcmpi('SOA chemicals', shift_method(sp_i));
        % SOA shift
        %fit_shift=[0.000590960272608028;2.93659435788714e-05;0.793966315974577;0.314333201045729];  % old
        fit_shift=[0.000149021455554774;0.00473627706154738;0.869057801919811;0.564783492434737];% updated to only include hydroperoxide containing SOA products
    elseif strcmpi('PEG', shift_method(sp_i));
        %PEG shift
        fit_shift=[0.00544768078879267;3.86433605822883;-0.267168022244528;0.255486696379870];
        
    elseif strcmpi('ketone', shift_method(sp_i))
        
        % ketones of little change for CH2CO
        fit_shift=[0.00453425820008037;0.000648450309348991;0.138144408029760;0.352454367906330];
        
    elseif strcmpi('hydroperoxide', shift_method(sp_i))
        % hydroperoxide
        fit_shift=[8.17160348517250e-06;4.53182593743281e-07;0.966089559236154;0.459433193460024];
        
    elseif strcmpi('ether', shift_method(sp_i))
        % ether
        fit_shift=[2.44336043644347e-05;0.000158316167479487;0.284974095922167;0.229338647030993];
        
        
    elseif strcmpi('ester', shift_method(sp_i))
        % ester,  very hydrophobic
        fit_shift=[-1.29324633442629;0.00108128314380665;1.24051435479678;0.405354156019572];
        
    else
        warning('no O2C and MW system selected')
        calc_shift=false;
        
    end
    
    if calc_shift
        % new function form
        %replace O2C value
        O2C_eqv(sp_i,1)=O2C(sp_i)./(1+fit_shift(3,1).*exp(-(O2C(sp_i))*fit_shift(1,1)));
        % change mass ratio
        MW_old=(18.016./molarmass_ratio(sp_i));
        MW_new=MW_old./(1+fit_shift(4,1).*exp(-(MW_old)*fit_shift(2,1)));
        
        molarmass_ratio_eqv(sp_i,1)=18.016./MW_new;
    else
        O2C_eqv(sp_i,1)=O2C(sp_i);
        molarmass_ratio_eqv(sp_i,1)=molarmass_ratio(sp_i);
    end
    
end


end

