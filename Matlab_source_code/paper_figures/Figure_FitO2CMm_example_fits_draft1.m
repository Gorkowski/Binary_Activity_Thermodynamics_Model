%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2018-Jun-01  4:16 PM
% Copyright 2018 Kyle Gorkowski 
%% 

% plot of data used in McGlashan fit vs Mr and O2C phase diagram

clear
load('binary_fit_data_imported.mat') % load AIOMFAC binary mixture data

Molar_mass_range=[45:1:600]'; % mass range to plot
water_mass=18.015;

data_start=1; % data range to plot
data_end=160;

S=size(O2C);
indexes_data=[1:S(1,2)]';
compound_M=water_mass./Mr;

%get data points
O2C_notHP=[O2C(data_start:data_end)];
Mmass_notHP=[compound_M(data_start:data_end)];

% phase sep. line
O2C_phase_sep_comp=single_phase_O2C_point_KGv3(water_mass./Molar_mass_range);

% mixtures to use in graph
data1_i=5;
data2_i=25;

% run BAT comparison line
mole_frac_scan=[1:-0.0001:0]';
mode='hydroxyl'; %    'hydroxyl'='carboxyl', 'hydroperoxideSOA', 'PEG', 'hydroperoxide'

[~, ~, ycal_water1, ycalc_org1, activity_water1, activity_org1, ~, ~,Gibbs_RT, dGibbs_RTdx2]...
    =BAT_properties_calculation_v1(mole_frac_scan, O2C(data1_i),  H2C(data1_i), Mr(data1_i),mode,[]);

[~, ~, ycal_water1, ycalc_org2, activity_water2, activity_org2, ~, ~,Gibbs_RT, dGibbs_RTdx2]...
    =BAT_properties_calculation_v1(mole_frac_scan, O2C(data2_i),  H2C(data2_i), Mr(data2_i),mode,[]);

% get AIOMFAC comparison line
mole_org1=molef_comp2(:,data1_i);
mole_water1=molef_comp1(:,data1_i);

AIOMFAC_act_org1=yin2(:,data1_i).*mole_org1;
AIOMFAC_act_water1=yin1(:,data1_i).*mole_water1;

mole_org2=molef_comp2(:,data2_i);
mole_water2=molef_comp1(:,data2_i);

AIOMFAC_act_org2=yin2(:,data2_i).*mole_org2;
AIOMFAC_act_water2=yin1(:,data2_i).*mole_water2;


%% figure
plot_name=['paper figure input data and BAT fit'];
 paper_postion=[0, 0, 6.5, 3].*1;

figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', paper_postion+.5,'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]


axes1 = axes('Parent',figure1,...
    'Position',[0.09 0.19 0.439761904761905 0.780238095238095],...
    'LineWidth',1.75,...
    'FontSize',12,'Color',[0.729411780834198 0.831372559070587 0.95686274766922]);  %'CLim',[80 300],
set(axes1,'FontSize',12,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% view(axes1,[175 35]);
grid(axes1,'on');
hold(axes1,'all');
hold on
area(Molar_mass_range,O2C_phase_sep_comp,'DisplayName', 'Middle O2C model ', 'FaceColor',[0.756862759590149 0.866666674613953 0.776470601558685]);
area(Molar_mass_range,O2C_phase_sep_comp.*.2,'DisplayName', 'Lower O2C model','FaceColor',[0.23137255012989 0.443137258291245 0.337254911661148]);


plot(Mmass_notHP,O2C_notHP, 'DisplayName','hydroxyl functionality','MarkerFaceColor',[0 0 0],...
    'MarkerEdgeColor',[1 1 1],...
    'MarkerSize',4,...
    'Marker','square',...
    'LineStyle','none',...
    'Color',[1 1 1]);

xlim(axes1,[50 500])
ylim([0,2.1])
ylabel('O:C')
xlabel('Molar Mass (g/mole)')
% legend('show')
axes2 = axes('Parent',figure1,...
    'Position',[0.7 0.601190476190476 0.28 0.363095238095236],...
    'LineWidth',1.75,...
    'FontSize',12,'Color',[0.729411780834198 0.831372559070587 0.95686274766922]);  %'CLim',[80 300],
set(axes2,'FontSize',12,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75,'Xcolor', 'none');
hold on
grid on
disp( char(compName(data1_i)))
plot(1-mole_frac_scan,activity_water1,'-','DisplayName', ['BAT a_w'],...
    'LineWidth',1,'Color',[0 0.450980392156863 0.741176470588235])
plot(mole_water1,AIOMFAC_act_water1,'--','LineWidth',1,'DisplayName', 'AIOMFAC a_w',...
    'Color',[0 0.450980392156863 0.741176470588235])
plot(1-mole_frac_scan,activity_org1,'-','DisplayName', [char(compName(data1_i)) ' a_{org} BAT'],'LineWidth',1,...
    'Color',[0.149019607843137 0.149019607843137 0.149019607843137])
plot(mole_water1,AIOMFAC_act_org1,'--','LineWidth',1,'DisplayName', [char(compName(data1_i)) ' a_{org} AIOMFAC'],...
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
disp( char(compName(data2_i)))
plot(1-mole_frac_scan,activity_water2,'-','DisplayName', ['BAT a_w'],...
    'LineWidth',1,'Color',[0 0.450980392156863 0.741176470588235])
plot(mole_water2,AIOMFAC_act_water2,'--','LineWidth',1,'DisplayName', 'AIOMFAC a_w',...
    'Color',[0 0.450980392156863 0.741176470588235])
plot(1-mole_frac_scan,activity_org2,'-','DisplayName', [char(compName(data2_i)) ' a_{org} BAT'],'LineWidth',1,...
    'Color',[0.149019607843137 0.149019607843137 0.149019607843137])
plot(mole_water2,AIOMFAC_act_org2,'--','LineWidth',1,'DisplayName', [char(compName(data2_i)) ' a_{org} AIOMFAC'],...
    'Color',[0.149019607843137 0.149019607843137 0.149019607843137])
ylabel('Activity')
xlabel('Water Mole Fraction')
print(figure1,plot_name,'-dpng','-r600')
saveas(figure1, [plot_name '.fig'], 'fig')






