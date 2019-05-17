function Run_VBSBAT_from_input_file(file_full_path)
%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2019-Mar-26 12:55 PM
% Copyright 2019 Kyle Gorkowski 
%% 


if not(exist('file_full_path'))
    % set default options
    [file,path] = uigetfile('*.*','Select VBSBAT input file');
    file_full_path=[path file];
end
% move to file location
back_slashes=strfind(file_full_path,'\');
cd(file_full_path(1:back_slashes(end)))

% read input file
simulation_input=read_input_file(file_full_path);

% make output folder
mkdir('Outputs')
cd('Outputs')


for i=1:size(simulation_input,2)
    
    
    
    % check if effective Csat needs to be calculated
    if sum(simulation_input(i).system.Csat_j_value)==0
        if sum(simulation_input(i).system.optional_Cliquid_ugPm3)>0 % checks to make sure those values are non-zero
            
            if isfield(simulation_input(i).system,'calculate_Csat_j_with_aw')
                aw_to_convert_at=simulation_input(i).system.calculate_Csat_j_with_aw;
            else
                aw_to_convert_at=0;
            end
            % estimates effective Csat value at dry conditions
            [Csat_approx]=VBS_equilibration_extractCsat_withLLEpartition_KGv2(simulation_input(i).system.optional_Cliquid_ugPm3, simulation_input(i).system.optional_Cstar_ugPm3, ...
                aw_to_convert_at, simulation_input(i).system.Molecular_weight, simulation_input(i).system.O2C_values, simulation_input(i).system.H2C_values,...
                simulation_input(i).system.BAT_functional_group, simulation_input(i).BAT_refinement_mode);
        else
            error('needs input of eff. Csat_j (ug/m3) or Cstar (ug/m3) and C^liquid (ug/m3), [if you want only the condensed phase set Csat_j to small value i.e. 10^-10]')
        end

    else
        Csat_approx=simulation_input(i).system.Csat_j_value;
    end
    
    
    [C_OA_PM, Caq_PM, kappaHGF, details]=VBS_BAT_simulation_v2(...
    Csat_approx, simulation_input(i).system.C_OM_ugPm3, ...
    simulation_input(i).system.O2C_values, simulation_input(i).system.H2C_values, simulation_input(i).system.Molecular_weight, ...
    simulation_input(i).water_activity, simulation_input(i).system.BAT_functional_group, simulation_input(i).BAT_refinement_mode, simulation_input(i).VBSBAT_options,simulation_input(i).run_name );
    
    % save out simulation as matlab file
    simulation_input_settings=simulation_input(i);
    save(['VBSBAT_sim_' simulation_input(i).run_name '.mat'],...
        'C_OA_PM', 'Caq_PM', 'kappaHGF', 'details','Csat_approx',...
        'simulation_input_settings')

end

cd ../