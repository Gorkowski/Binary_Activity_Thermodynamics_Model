%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2018-Jun-01  4:16 PM
% Copyright 2018 Kyle Gorkowski 
%% 

% plot of data used in McGlashan fit vs Mr and O2C phase diagram
clear
load('binary_fit_data_imported_SOA')

% error O2C for hydroxyl average O2C error 0.040859 max error 0.10004
error_O2C_value=0.040859;

data1_i=24; % data mixture run

mole_frac_scan=[1:-0.0001:0]';
mode='hydroperoxideSOA'; %    {'hydroxyl';'carboxyl';'ether';'ketone';'ester';'hydroperoxide';'hydroperoxideSOA';'PEG'}
% run BAT sim for hydroperoxideSOA
[~, ~, ycal_water1, ycalc_org1, activity_water1, activity_org1, ~, ~,Gibbs_RT, dGibbs_RTdx2]...
    =BAT_properties_calculation_v1(mole_frac_scan, O2C(data1_i),  H2C(data1_i), Mr(data1_i),mode,[]);
% run BAT sim for hydroxyl
[~, ~, ~, ~, activity_water1_OH, activity_org1_OH, ~, ~,Gibbs_RT, dGibbs_RTdx2]...
    =BAT_properties_calculation_v1(mole_frac_scan, O2C(data1_i),  H2C(data1_i), Mr(data1_i),'hydroxyl',[]);

% get AIOMFAC data for this run
mole_org1=molef_comp2(:,data1_i);
mole_water1=molef_comp1(:,data1_i);
AIOMFAC_act_org1=yin2(:,data1_i).*mole_org1;
AIOMFAC_act_water1=yin1(:,data1_i).*mole_water1;
data1_name=char(compName(data1_i));

% load second data set to compare translation method
load('binary_fit_data_imported_PEG')

data2_i=2;
mode='PEG'; %    {'hydroxyl';'carboxyl';'ether';'ketone';'ester';'hydroperoxide';'hydroperoxideSOA';'PEG'}
% simulation BAT with PEG transform
[~, ~, ycal_water1, ycalc_org2, activity_water2, activity_org2, ~, ~,Gibbs_RT, dGibbs_RTdx2]...
    =BAT_properties_calculation_v1(mole_frac_scan, O2C(data2_i),  H2C(data2_i), Mr(data2_i),mode,[]);
% simulation BAT with hydroxyl
[~, ~, ~, ~, activity_water2_OH, activity_org2_OH, ~, ~,Gibbs_RT, dGibbs_RTdx2]...
    =BAT_properties_calculation_v1(mole_frac_scan, O2C(data2_i),  H2C(data2_i), Mr(data2_i),'hydroxyl',[]);

% get AIOMFAC mixture data
mole_org2=molef_comp2(:,data2_i);
mole_water2=molef_comp1(:,data2_i);
AIOMFAC_act_org2=yin2(:,data2_i).*mole_org2;
AIOMFAC_act_water2=yin1(:,data2_i).*mole_water2;
data2_name=char(compName(data2_i));


% caculate phase line
shift_list={'hydroxyl';'carboxyl';'ether';'ketone';'ester';'hydroperoxide';'hydroperoxideSOA';'PEG'};
O2C_simple=[0:.01:2]';
MW_simple=[50:10:500]';
fgroups=size(shift_list);
for i=1:fgroups(1,1)
    
[O2C_eqv(:,i),~] = convert_chemical_structure_to_OH_eqv_v3(O2C_simple, O2C_simple.*0, char(shift_list(i,1)) );
[~,molarmass_ratio_eqv(:,i)] = convert_chemical_structure_to_OH_eqv_v3(MW_simple.*0, 18.016./MW_simple, char(shift_list(i,1)) );

% get OH phase sep O2C
O2C_single_phase_OH=single_phase_O2C_point_KGv3(molarmass_ratio_eqv(:,i));
O2C_single_phase_strucutre(:,i)=interp1(O2C_eqv(:,i),O2C_simple,O2C_single_phase_OH,'linear');
O2C_single_phase_strucutre_error_upper(:,i)=interp1(O2C_eqv(:,i)+error_O2C_value,O2C_simple,O2C_single_phase_OH,'linear');
O2C_single_phase_strucutre_error_lower(:,i)=interp1(O2C_eqv(:,i)-error_O2C_value,O2C_simple,O2C_single_phase_OH,'linear');

end
MW_eqv=18.016./molarmass_ratio_eqv;
limts_O2C=max(max(O2C_eqv));


%% figure
plot_name=['paper figure functional groups phaseline'];
 paper_postion=[0, 0, 6.5, 3].*1;

figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', paper_postion+.5,'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]


axes1 = axes('Parent',figure1,...
    'Position',[0.09 0.19 0.439761904761905 0.780238095238095],...
    'LineWidth',1.75,...
    'FontSize',12);  %'CLim',[80 300],
set(axes1,'FontSize',12,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% view(axes1,[175 35]);
grid(axes1,'on');
hold(axes1,'all');
% phase line
hold on


fill([MW_simple;flip(MW_simple)],[O2C_single_phase_strucutre_error_upper(:,1);...
    flip(O2C_single_phase_strucutre_error_lower(:,1))],'r','FaceAlpha',0.25,...
    'FaceColor',[0.3 0.3 0.3],...
    'EdgeColor','none');

fill([MW_simple;flip(MW_simple)],[O2C_single_phase_strucutre_error_upper(:,6);...
    flip(O2C_single_phase_strucutre_error_lower(:,6))],'r','FaceAlpha',0.25,...
    'FaceColor',[0.3 0.3 0.3],...
    'EdgeColor','none');

plot( MW_simple,   O2C_single_phase_strucutre(:,1),'-k','LineWidth', 4,'DisplayName', char(shift_list(1,1)))

for i=2:fgroups(1,1)
    
    if strcmpi(char(shift_list(i,1)),'ester')
        % ester always phase seperate
        plot( MW_simple,   MW_simple./MW_simple,':', 'DisplayName', char(shift_list(i,1)),'LineWidth', 2)
    else
plot( MW_simple,   O2C_single_phase_strucutre(:,i), 'DisplayName', char(shift_list(i,1)),'LineWidth', 2)
    end
end
 
 ylim([0,1.1])
  xlim([50,500])

xlabel('structure MW')
ylabel('structure O:C')
%legend1 = legend(axes1,'show');

% legend('show')
axes2 = axes('Parent',figure1,...
    'Position',[0.7 0.601190476190476 0.28 0.363095238095236],...
    'LineWidth',1.75,...
    'FontSize',12,'Color',[0.729411780834198 0.831372559070587 0.95686274766922]);  %'CLim',[80 300],
set(axes2,'FontSize',12,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75,'Xcolor', 'none');
hold on
grid on
disp( data1_name)
plot(1-mole_frac_scan,activity_water1,'-','DisplayName', data1_name,...
    'LineWidth',1,'Color',[0 0.450980392156863 0.741176470588235])
plot(1-mole_frac_scan,activity_water1_OH,'-','DisplayName', 'OH',...
    'LineWidth',1,'Color',[.6 0.6 0.6])
plot(mole_water1,AIOMFAC_act_water1,'--','LineWidth',1,'DisplayName', 'AIOMFAC',...
    'Color',[0 0.450980392156863 0.741176470588235])
plot(1-mole_frac_scan,activity_org1,'-','DisplayName', data1_name,'LineWidth',1,...
    'Color',[0.149019607843137 0.149019607843137 0.149019607843137])
plot(1-mole_frac_scan,activity_org1_OH,'-','DisplayName', 'OH',...
    'LineWidth',1,'Color',[.6 0.6 0.6])
plot(mole_water1,AIOMFAC_act_org1,'--','LineWidth',1,'DisplayName', 'AIOMFAC',...
    'Color',[0.149019607843137 0.149019607843137 0.149019607843137])
ylabel('Activity')

axes3 = axes('Parent',figure1,...
    'Position',[0.7 0.19 0.28 0.37547619047619],...
    'LineWidth',1.75,...
    'FontSize',12,'Color',[0.756862759590149 0.866666674613953 0.776470601558685]);  %'CLim',[80 300],
set(axes3,'FontSize',12,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
hold on
grid on
disp( data2_name)
plot(1-mole_frac_scan,activity_water2,'-','DisplayName', data2_name,...
    'LineWidth',1,'Color',[0 0.450980392156863 0.741176470588235])
plot(1-mole_frac_scan,activity_water2_OH,'-','DisplayName', 'OH',...
    'LineWidth',1,'Color',[.6 0.6 0.6])
plot(mole_water2,AIOMFAC_act_water2,'--','LineWidth',1,'DisplayName', 'AIOMFAC',...
    'Color',[0 0.450980392156863 0.741176470588235])
plot(1-mole_frac_scan,activity_org2,'-','DisplayName', data2_name,'LineWidth',1,...
    'Color',[0.149019607843137 0.149019607843137 0.149019607843137])
plot(1-mole_frac_scan,activity_org2_OH,'-','DisplayName', 'OH',...
    'LineWidth',1,'Color',[.6 0.6 0.6])
plot(mole_water2,AIOMFAC_act_org2,'--','LineWidth',1,'DisplayName', 'AIOMFAC',...
    'Color',[0.149019607843137 0.149019607843137 0.149019607843137])
ylabel('Activity')
xlabel('Water Mole Fraction')
print(figure1,plot_name,'-dpng','-r600')
saveas(figure1, [plot_name '.fig'], 'fig')






