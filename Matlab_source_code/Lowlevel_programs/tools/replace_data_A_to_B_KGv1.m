function output=replace_data_A_to_B_KGv1(data_in, A ,B)
% Created by Kyle Gorkowski on 2012-Jan-24  4.20 PM
% Updated by Kyle Gorkowski [GORKOWFALCON] on 2016-Sep-10  1:54 PM
% Updated by Kyle Gorkowski [GORKOWFALCON] on 2016-Oct-08 11:38 AM
% works best with intigers as you can have floating point problems for
% 16bit data

%% 

% s=size(data_in);

if isnan(A)
    mask=not(isnan(data_in)); % makes mask
    [row,col]=find(+mask);
    
    data_size=size(data_in);
    data_in2=zeros(data_size);
    if data_size(1,2)==1
        col=1;
    end
    if data_size(1,1)==1
        row=1;
    end
    data_in2(row,col)=data_in(row,col); % pulls out all non nan values
    data_in=data_in2;
else
    mask=not(data_in==A);
end
 % data that satisfies conditions 
full_mask=+mask;

good_data=full_mask.*data_in; % saves all the good data and replaces everything else with zeros

if isnan(B)
    replacement_data=full_mask./full_mask-1;
else
    replacement_data=not(full_mask).*B;
end


output=good_data+replacement_data;
