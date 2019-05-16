function [output_time, output_data]=sort_and_remove_duplicates_KGv1(time, data)
% Created by Kyle Gorkowski on 2012-Jun-19 12:48 PM
%[output_time, output_data]=sort_and_remove_duplicates_KGv1(time, data)
D=size(time);

both_sets=sortrows([time, data], 1);
B=size(both_sets);
time=both_sets(:,1);
data=both_sets(:,2:B(1,2));

for i=1:1:D(1,1)
    time_min=min(abs([time(1:i-1,1);time(i+1:D(1,1),1)]-ones(D(1,1)-1,1).*time(i,1)));
    if time_min==0
        time(i,1)=NaN;
        
    end
end


[output_time, output_data]=Removes_nan_rows_dual(time, data);