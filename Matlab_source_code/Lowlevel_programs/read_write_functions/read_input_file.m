function simulation_input=read_input_file(filename)
%clear
%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2019-Mar-26 11:38 AM
% Copyright 2019 Kyle Gorkowski 
%% 
% import file for VBSBAT calc

% testing exmaple
% clear
% %% Initialize variables.
%   filename = 'C:\Users\kkgor\OneDrive\Professional Current DnD\McGill_Docs\AIOMFAC\VBSBAT_standalone\Copy_of_VBSBAT_test_aPsoa.txt';
%       filename = 'C:\Users\kkgor\OneDrive\Professional Current DnD\McGill_Docs\AIOMFAC\VBSBAT_standalone\VBSBAT_test_aPsoa.csv';

delimiter = {''};

%% Format for each line of text:
formatSpec = '%s%[^\n\r]';
fileID = fopen(filename,'r','n','UTF-8');
% fseek(fileID, 3, 'bof');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);
input_file = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;


%% generate input structures




run_start=false;
settings_start=false;
system_start=false;


sim_i=1;
for i=1:size(input_file,1)
    texLine=char(input_file(i,:));
    texLength1=length(texLine)+1;
    if contains(texLine, 'run name=') % sets run name
        if system_start
            sim_i=sim_i+1;
        end
        
        equal_i=strfind(texLine,'=');
        first_comma=strfind(texLine,',');
        if isempty(first_comma)
            first_comma(1,1)=texLength1;
        end
        simulation_input(sim_i).run_name=texLine(equal_i+1:first_comma(1,1)-1);
        run_start=true;
        settings_start=false;
        system_start=false;
    end
    
    if run_start && contains(texLine, 'McGlashan_refinement_mode=') % sets refinement mode
        equal_i=strfind(texLine,'=');
        first_comma=strfind(texLine,',');
        if isempty(first_comma)
            first_comma(1,1)=texLength1;
        end
        simulation_input(sim_i).McGlashan_refinement_mode=texLine(equal_i+1:first_comma(1,1)-1);
    end
    
    if run_start && contains(texLine, 'VBSBAT_options') % put in options

        equal_i=strfind(texLine,'=');
        first_comma=strfind(texLine,',');
        if isempty(first_comma)
            first_comma(1,1)=texLength1;
        end
        
        if not(settings_start) % initiate VBSBAT.options structure
            if contains(texLine, 'VBSBAT_options.run_mode_used=')
                simulation_input(sim_i).VBSBAT_options= default_VBSBAT_options(texLine(equal_i+1:first_comma(1,1)-1));
            else
                simulation_input(sim_i).VBSBAT_options= default_VBSBAT_options('default');
            end
            settings_start=true;
        
        else % change options
            start_i=strfind(texLine,'VBSBAT_options');
            option_tochange=texLine(start_i:equal_i-1);
            
            if contains(texLine, 'VBSBAT_options.q_alpha.q_alpha_bounds') || contains(texLine, 'VBSBAT_options.q_alpha.q_alpha_bounds_mean')
                %two inputs      
                bounds(1,1)=str2double(texLine( equal_i+1:first_comma(1,1)-1));
                bounds(1,2)=str2double(texLine( first_comma(1,1)+1:first_comma(1,2)-1));
                eval(['simulation_input(sim_i).' option_tochange '=bounds ;'])
            else % single input
                numeric_entry=regexp(texLine(equal_i+1:first_comma(1,1)-1),'\d');
                if isempty(numeric_entry) % string
                    eval(['simulation_input(sim_i).' option_tochange '=texLine( equal_i+1:first_comma(1,1)-1);'])
                else % number
                    eval(['simulation_input(sim_i).' option_tochange '=str2double(texLine( equal_i+1:first_comma(1,1)-1)) ;'])
                end
            end
        end 
    end
    
    if contains(texLine, 'water activity') %import water activity
        if not(settings_start) % check if settings has not be initiated. 
           simulation_input(sim_i).VBSBAT_options=default_VBSBAT_options('default');
        end
        
        aw_line=char(input_file(i+1,:));
        aw_input_cell = textscan(aw_line,'%f,');
        simulation_input(sim_i).water_activity = cell2mat(aw_input_cell(1,1));
    end
    
    if contains(texLine, 'M (g/mol)') %system properties reader switch
        system_start=true;
        run_start=false;
        settings_start=false;
        system_i=1;
    end
                    
    % imports system properties
    if system_start && not(contains(texLine, 'M (g/mol)')) && sum(regexp(texLine,'\d'))>0 % past the header line and contains some numbers
        first_comma=strfind(texLine,',');
        % something funny with textscan, as I think I only need to call texscan once,
        % but it kept on messing up on the second half of the numbers.
        system_input_cell = [textscan(texLine(1:first_comma(1,7)-1),'%f,%f,%f,%f,%f,%s,%s,'),textscan(texLine(first_comma(1,7)+1:first_comma(1,9)-1),'%f,%f')];
        
        simulation_input(sim_i).system.Molecular_weight(system_i,1)=cell2mat(system_input_cell(1,1));
        simulation_input(sim_i).system.O2C_values(system_i,1)=cell2mat(system_input_cell(1,2));
        simulation_input(sim_i).system.H2C_values(system_i,1)=cell2mat(system_input_cell(1,3));
        simulation_input(sim_i).system.Csat_j_value(system_i,1)=cell2mat(system_input_cell(1,4));
        simulation_input(sim_i).system.C_OM_ugPm3(system_i,1)=cell2mat(system_input_cell(1,5));
        simulation_input(sim_i).system.BAT_functional_group(system_i,1)=system_input_cell{1,6};
        simulation_input(sim_i).system.speciesName(system_i,1)=system_input_cell{1,7};
        simulation_input(sim_i).system.optional_Cstar_ugPm3(system_i,1)=cell2mat(system_input_cell(1,8));
        simulation_input(sim_i).system.optional_Cliquid_ugPm3(system_i,1)=cell2mat(system_input_cell(1,9));

        system_i=system_i+1;
    end
    
    if run_start && contains(texLine, 'calculate_Csat_j_with_aw=')
%         calculate_Csat_j_with_aw=0
        equal_i=strfind(texLine,'=');
        first_comma=strfind(texLine,',');
        if isempty(first_comma)
            first_comma(1,1)=texLength1;
        end
        simulation_input(sim_i).system.calculate_Csat_j_with_aw=str2double(texLine( equal_i+1:first_comma(1,1)-1));
    end
end






