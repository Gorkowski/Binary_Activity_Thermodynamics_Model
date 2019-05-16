function [data1, data2]=Removes_nan_rows_dual(data1_nan, data2_other)

% Created by Kyle Gorkowski on 2012-Jan-24  4.20 PM
%% [data1, data2]=Removes_nan_rows_dual(data1_nan, data2_other)
nan_rows=sum(isnan(data1_nan),2);

s=size(data1_nan);
q=0;
for i=1:1:s(1,1)

        if nan_rows(i,1)==0
            q=q+1;
            data1_nan(q,:)=data1_nan(i,:);
            data2_other(q,:)=data2_other(i,:);

        end

end
data1=data1_nan(1:q,:);
data2=data2_other(1:q,:);
