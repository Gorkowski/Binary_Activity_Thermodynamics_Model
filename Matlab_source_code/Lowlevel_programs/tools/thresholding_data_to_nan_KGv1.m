function output=thresholding_data_to_nan_KGv1(data_in, lower_bound, upper_bound)
%output=thresholding_data_to_nan_KGv1(data_in, lower_bound, upper_bound)
% Created by Kyle Gorkowski on 2012-Jan-24  4.20 PM
%% 

s=size(data_in);

for i=1:1:s(1,1)
    for k=1:1:s(1,2)
        if data_in(i,k)>upper_bound || data_in(i,k)<lower_bound
            data_in(i,k)=NaN;
        end
    end
end
output=data_in;
