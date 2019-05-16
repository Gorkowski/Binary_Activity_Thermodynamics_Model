function output=Removes_zero_rows(data_in)

% Created by Kyle Gorkowski on 2012-Jan-24  4.20 PM
%% output=Removes_nan_rows(data_in)
zero_rows=any(data_in,2);

s=size(data_in);
q=0;
for i=1:1:s(1,1)

        if zero_rows(i,1)==1
            q=q+1;
            data_in(q,:)=data_in(i,:);
        end

end
output=data_in(1:q,:);
