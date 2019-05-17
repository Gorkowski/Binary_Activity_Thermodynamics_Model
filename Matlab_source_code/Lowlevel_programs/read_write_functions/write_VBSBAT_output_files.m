function write_VBSBAT_output_files(name,sim_details)
%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2019-Mar-25 11:26 AM
% Copyright 2019 Kyle Gorkowski 
%% 
% write text files
% 
% clear
% load('VBSBAT_calc.mat','details_aP')
% 
% name='ap_simulation';
% sim_details=details_aP;

[options_2write]=write_VBSBAT_options_to_cell(sim_details.inputs.VBSBAT_options);
% sim_details.inputs.BAT_refinement_mode

[data_header, data_matrix, input_header, input_matrix, input_functional_groups, input_aw]=gen_simple_output_data_from_VBSBAT_sim(sim_details);


%% write input file used
input_file_name=['VBSBAT_input_used_' name '.txt'];
input_file = fopen(input_file_name,'w');
% model name 
% 
fprintf(input_file,'%s\n', ['run name=' name]);
fprintf(input_file,'%s\n', ['BAT_refinement_mode=' sim_details.inputs.BAT_refinement_mode]);

s_options=size(options_2write);
for i=1:s_options(1,1)
    fprintf(input_file,'%s\n', string(options_2write(i,1)));
end

fprintf(input_file,'%s\n', ['water activity']);
fprintf(input_file,'%f,', transpose(input_aw));
fprintf(input_file,'\n%s\n',['system properties']);

input_header=pad(string(input_header),18,'both');
input_header=[input_header;'BAT functional group']; 
fprintf(input_file,'%s,', input_header);
fprintf(input_file,'\n');

% save species input
for i=1:size(input_matrix,1)
    fprintf(input_file,'%18.6E,', input_matrix(i,:));
    fprintf(input_file,'%s\n', pad(string(input_functional_groups(i,1)), 18,'both'));
end
fclose(input_file);


%% write simple output file

output_file_name=['VBSBAT_simple_output_' name '.txt'];
output_file = fopen(output_file_name,'w');

fprintf(output_file,'%s,', pad(string(data_header), 18,'both'));
fprintf(output_file,'\n');

%save output
for i=1:size(data_matrix,1)
    fprintf(output_file,'%18.6E,', data_matrix(i,:));
    fprintf(output_file,'\n');

end
fclose(output_file);





