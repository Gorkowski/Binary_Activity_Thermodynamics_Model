%%
% Created by Kyle Gorkowski [GORKOWFALCON] on 2018-Dec-18  6:44 AM
% Copyright 2018 Kyle Gorkowski 
%%
% loads AIOMFAC data and generates comparison plots. Using direct calls to
% VBSBAT and not using input file format
% outputs 

clear
BAT_refinement_mode='interpolate'; % 'perfect water activity' 'none' 'interpolate'

% individual species simulation VBSBAT 
VBSBAT_options=default_VBSBAT_options('default');%'robust' 'default'
% VBSBAT_options.q_alpha.min_spread_in_aw=10^-6;
% VBSBAT_options.q_alphaVBS.method_to_use='individual';%mean_prop or individual
% mean prop species simulation VBSBAT
VBSBAT_options_mean=VBSBAT_options;
VBSBAT_options_mean.VBSBAT_NN_options.NN_type='individual_properties'; 
VBSBAT_options_mean.mean_BAT_functional_group='hydroperoxideSOA';

extension_name='_final';

load('sim_database_in_matlab.mat') % loads AIOMFAC sim

aw_series=data_systems(1).basic.aw_dataset; % get water activity points
[~,min_aw_i]=min(aw_series);
Saw=size(aw_series);

aw_series_VBSBAT=[0.9999999; 0.999999; 0.99995; aw_series];
added_aw=4;

%% alpha pinene SOA system
index=3; % index in AIOMFAC simulaiton database
%AIOMFAC information
AIOMFAC_aPsoa_run=data_systems(index).comp(1).nd;
AIOMFAC_aPsoa=data_systems(index).Org_mass_ugPm3.Org_Conc_massPM;
AIOMFAC_aPsoa_water=data_systems(index).mass_ugPm3.massPM_water;
AIOMFAC_aPsoa_dry=repmat(AIOMFAC_aPsoa(min_aw_i),Saw(1,1),1);

AIOMFAC_aPsoa_run_NHSO4=data_systems(1).comp(1).nd;
AIOMFAC_aPsoa_NHSO4=data_systems(1).Org_mass_ugPm3.Org_Conc_massPM;
AIOMFAC_aPsoa_water_NHSO4=data_systems(1).mass_ugPm3.massPM_water;
AIOMFAC_aPsoa_kappaHGF=data_systems(index).kappa.kappaHGF;

% make inputs for VBSBAT simulaiton
% names=                         ["C107OOH";      "C97OOH";      "C108OOH";   "ALDOL_dimer_C19H28O7";  "PINIC";       "C921OOH";         "C812OOH";   "ESTER_dimer";      "C811OH";         "C813OOH"]
%BAT_functional_group={'hydroperoxide'; 'hydroperoxide'; 'hydroperoxide'; 'hydroperoxide'; 'carboxyl';  'hydroperoxide'; 'hydroperoxide'; 'ester' ; 'hydroperoxide' ; 'hydroperoxide'};
BAT_functional_group={'hydroperoxideSOA'; 'hydroperoxideSOA'; 'hydroperoxideSOA'; 'hydroperoxideSOA'; 'carboxyl';  'hydroperoxideSOA'; 'hydroperoxideSOA'; 'ester' ; 'hydroperoxideSOA' ; 'hydroperoxideSOA'};
VBSBAT_options.mean_BAT_functional_group='hydroperoxideSOA';% {'hydroxyl';'carboxyl';'ether';'ketone';'ester';'hydroperoxide';'hydroperoxideSOA';'PEG'};

H2C_values=[16/10; 16/9; 16/10; 28/19; 14/9; 16/9; 14/8;  28/18;  14/8; 14/8] ;
O2C_values=data_systems(index).Org_mass_ugPm3.Org_O2C';
MolarMass=data_systems(index).Org_mass_ugPm3.Org_MW_gMole';
sim_name='aPsoa';

% calc effective Csat
[Csat_approx]=VBS_equilibration_extractCsat_withLLEpartition_KGv2(data_systems(index).Org_mass_ugPm3.Org_Conc_massPM_comp(min_aw_i,:)', data_systems(index).Org_mass_ugPm3.Org_Cstar_comp(min_aw_i,:)', ...
    aw_series(min_aw_i,1), MolarMass, O2C_values, H2C_values,  BAT_functional_group, BAT_refinement_mode);

% run VBSBAT simulation
[VBSBAT_aPsoa_C_OA_PM, VBSBAT_aPsoa_Caq_PM, VBSBAT_aPsoa_kappaHGF, details_aP]=VBS_BAT_simulation_v2(...
    Csat_approx, data_systems(index).Org_mass_ugPm3.Org_Conc_OM(min_aw_i,:)', ...
    O2C_values, H2C_values,  MolarMass, aw_series_VBSBAT, BAT_functional_group, BAT_refinement_mode, VBSBAT_options,sim_name );
disp(['fit per aw ' num2str(details_aP.fit.mean_time_per_fit) ' sec'])

%% alpha pinene SOA system, using mean properties
% calculate average dry prop for VBSBAT simulation
mass_fraction=data_systems(index).Org_mass_ugPm3.Org_Conc_massPM_comp(min_aw_i,:)'./sum(data_systems(index).Org_mass_ugPm3.Org_Conc_massPM_comp(min_aw_i,:)');
H2C_values_dry_mean=repmat(sum(H2C_values.*mass_fraction),size(H2C_values,1),1);
O2C_values_dry_mean=repmat(sum(O2C_values.*mass_fraction),size(O2C_values,1),1);
MolarMass_dry_mean=repmat(sum(MolarMass.*mass_fraction),size(MolarMass,1),1);
BAT_functional_group='hydroperoxideSOA';
sim_name='aPsoa_meanProp';

% calc. new effective Csat
[Csat_approx]=VBS_equilibration_extractCsat_withLLEpartition_KGv2(data_systems(index).Org_mass_ugPm3.Org_Conc_massPM_comp(min_aw_i,:)', data_systems(index).Org_mass_ugPm3.Org_Cstar_comp(min_aw_i,:)', ...
    aw_series(min_aw_i,1), MolarMass_dry_mean, O2C_values_dry_mean, H2C_values_dry_mean,  BAT_functional_group, BAT_refinement_mode);
% calc. VBSBAT simulation
[VBSBAT_aPsoa_Mean_C_OA_PM, VBSBAT_aPsoa_Mean_Caq_PM, VBSBAT_aPsoa_Mean_kappaHGF, details_aPsoa_mean]=VBS_BAT_simulation_v2(Csat_approx, data_systems(index).Org_mass_ugPm3.Org_Conc_OM(min_aw_i,:)', ...
    O2C_values_dry_mean, H2C_values_dry_mean,  MolarMass_dry_mean, aw_series_VBSBAT, BAT_functional_group, BAT_refinement_mode, VBSBAT_options_mean, sim_name);
disp(['fit per aw ' num2str(details_aPsoa_mean.fit.mean_time_per_fit) ' sec'])


%% Isoprene SOA system
index=4;
% get data from AIOMFAC database
AIOMFAC_IPsoa_run=data_systems(index).comp(1).nd;
AIOMFAC_IPsoa=data_systems(index).Org_mass_ugPm3.Org_Conc_massPM;
AIOMFAC_IPsoa_water=data_systems(index).mass_ugPm3.massPM_water;
AIOMFAC_IPsoa_dry=repmat(AIOMFAC_IPsoa(min_aw_i),Saw(1,1),1);
AIOMFAC_IPsoa_run_NHSO4=data_systems(2).comp(1).nd; % salt case
AIOMFAC_IPsoa_NHSO4=data_systems(2).Org_mass_ugPm3.Org_Conc_massPM;
AIOMFAC_IPsoa_water_NHSO4=data_systems(2).mass_ugPm3.massPM_water;
AIOMFAC_IPsoa_kappaHGF=data_systems(index).kappa.kappaHGF;

% enter VBSBAT data
%                     ["IEB1OOH";        "IEB2OOH";     "C59OOH";        "IEC1OOH";        "C58OOH";         "IEPOXA";   "C57OOH";       "IEPOXC";   "HIEB1OOH";       "INDOOH";          "IEACO3H";     "C525OOH";         "HIEB2OOH";     "IEC2OOH";         "INAOOH";        "C510OOH";    "INB1OOH";         "IECCO3H";         "INCOOH";        "INB2OOH";    "2-Methyltetrol_dimer"          ]
%BAT_functional_group={'hydroperoxide'; 'hydroperoxide'; 'hydroperoxide';  'hydroperoxide'; 'hydroperoxide'; 'hydroxyl'; 'hydroperoxide';  'hydroxyl'; 'hydroperoxide';  'hydroperoxide';  'hydroperoxide';  'hydroperoxide'; 'hydroperoxide'; 'hydroperoxide';  'hydroperoxide'; 'hydroperoxide'; 'hydroperoxide'; 'hydroperoxide'; 'hydroperoxide'; 'hydroperoxide'; 'hydroxyl'         } ;
BAT_functional_group={'hydroperoxideSOA'; 'hydroperoxideSOA'; 'hydroperoxideSOA';  'hydroperoxideSOA'; 'hydroperoxideSOA'; 'hydroxyl'; 'hydroperoxideSOA';  'hydroxyl'; 'hydroperoxideSOA';  'hydroperoxideSOA';  'hydroperoxideSOA';  'hydroperoxideSOA'; 'hydroperoxideSOA'; 'hydroperoxideSOA';  'hydroperoxideSOA'; 'hydroperoxideSOA'; 'hydroperoxideSOA'; 'hydroperoxideSOA'; 'hydroperoxideSOA'; 'hydroperoxideSOA'; 'hydroxyl'         } ;
H2C_values=[              10/5;              12/5;             10/5;            10/5;              10/5;      10/5;        10/5;          10/5;      10/5;               11/5;              8/5;              10/5;           10/5;           8/5;               11/5;            9/5;              11/5;            8/5;            11/5;             11/5;             23/10];
% {'hydroxyl';'carboxyl';'ether';'ketone';'ester';'hydroperoxide';'hydroperoxideSOA';'PEG'};
VBSBAT_options.mean_BAT_functional_group='hydroperoxideSOA';

O2C_values=data_systems(index).Org_mass_ugPm3.Org_O2C';
MolarMass=data_systems(index).Org_mass_ugPm3.Org_MW_gMole';
sim_name='IPsoa';

%Calculate effective Csat
[Csat_approx]=VBS_equilibration_extractCsat_withLLEpartition_KGv2(data_systems(index).Org_mass_ugPm3.Org_Conc_massPM_comp(min_aw_i,:)', data_systems(index).Org_mass_ugPm3.Org_Cstar_comp(min_aw_i,:)', ...
    aw_series(min_aw_i,1), MolarMass, O2C_values, H2C_values,  BAT_functional_group, BAT_refinement_mode);
% VBSBAT simulation
[VBSBAT_IPsoa_C_OA_PM, VBSBAT_IPsoa_Caq_PM, VBSBAT_IPsoa_kappaHGF, details_Isoa]=VBS_BAT_simulation_v2(Csat_approx,...
    data_systems(index).Org_mass_ugPm3.Org_Conc_OM(min_aw_i,:)', ...
    O2C_values, H2C_values,  MolarMass, aw_series_VBSBAT, BAT_functional_group, BAT_refinement_mode, VBSBAT_options, sim_name);
disp(['fit per aw ' num2str(details_Isoa.fit.mean_time_per_fit) ' sec'])

%% Isoprene SOA system, using average properties in VBSBAT
% average dry prop
mass_fraction=data_systems(index).Org_mass_ugPm3.Org_Conc_massPM_comp(min_aw_i,:)'./sum(data_systems(index).Org_mass_ugPm3.Org_Conc_massPM_comp(min_aw_i,:)');
H2C_values_dry_mean=repmat(sum(H2C_values.*mass_fraction),size(H2C_values,1),1);
O2C_values_dry_mean=repmat(sum(O2C_values.*mass_fraction),size(O2C_values,1),1);
MolarMass_dry_mean=repmat(sum(MolarMass.*mass_fraction),size(MolarMass,1),1);
BAT_functional_group='hydroperoxideSOA';
sim_name='IPsoa_meanProp';

[Csat_approx]=VBS_equilibration_extractCsat_withLLEpartition_KGv2(data_systems(index).Org_mass_ugPm3.Org_Conc_massPM_comp(min_aw_i,:)', data_systems(index).Org_mass_ugPm3.Org_Cstar_comp(min_aw_i,:)', ...
    aw_series(min_aw_i,1), MolarMass_dry_mean, O2C_values_dry_mean, H2C_values_dry_mean,  BAT_functional_group, BAT_refinement_mode);
[VBSBAT_IPsoa_Mean_C_OA_PM, VBSBAT_IPsoa_Mean_Caq_PM, VBSBAT_IPsoa_Mean_kappaHGF, details_IPsoa_mean]=VBS_BAT_simulation_v2(Csat_approx, data_systems(index).Org_mass_ugPm3.Org_Conc_OM(min_aw_i,:)', ...
    O2C_values_dry_mean, H2C_values_dry_mean,  MolarMass_dry_mean, aw_series_VBSBAT, BAT_functional_group, BAT_refinement_mode, VBSBAT_options_mean, sim_name);
disp(['fit per aw ' num2str(details_IPsoa_mean.fit.mean_time_per_fit) ' sec'])

% save data outputs
save('VBSBAT_calc.mat')








%% graphs   **************************************************************************************

% salt efflorescence for org mass
AIOMFAC_aPsoa_NHSO4(92:end)=AIOMFAC_aPsoa(92:end);
AIOMFAC_IPsoa_NHSO4(92:end)=AIOMFAC_IPsoa(92:end);



line_width=1;
fontsize=12;
paper_mult=1;

%% org mass
plot_name=['C_OA mass plot', extension_name];
paper_postion=[0, 0, 3.5, 3.5*1.25].*paper_mult;
figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', paper_postion.*1+.5,'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]


axes1 = axes('Parent',figure1,...
    'Position',[0.14 0.47 0.8 0.50],...
    'LineWidth',1.75,...
    'FontSize',fontsize, 'Xcolor', 'none');  %'CLim',[80 300],
set(axes1,'FontSize',fontsize,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% view(axes1,[175 35]);
grid(axes1,'on');
hold(axes1,'all');
hold on

%base line
% plot([0;100],zeros(2,1),'Color',[0 0 0],'LineWidth',1.75)

color_plot=[0 0.498039215803146 0]; %alpha pinene
p1=plot(100.*aw_series_VBSBAT,VBSBAT_aPsoa_Mean_C_OA_PM,'Color',color_plot.*.8, 'DisplayName' ,...
    ['VBSBAT aPsoa MeanProp run ' num2str(AIOMFAC_aPsoa_run)],'LineStyle','-','LineWidth',line_width,'Visible','on');
p1=plot(100.*aw_series,AIOMFAC_aPsoa_dry,'Color',color_plot, 'DisplayName' ,...
    ['AIOMFAC aPsoa Dry run ' num2str(AIOMFAC_aPsoa_run)],'LineStyle','-.','LineWidth',line_width,'Visible','on');
p1=plot(100.*aw_series,AIOMFAC_aPsoa_NHSO4,'Color',color_plot, 'DisplayName' ,...
    ['AIOMFAC aPsoa with (NH_4)_2SO_4 run ' num2str(AIOMFAC_aPsoa_run_NHSO4)],'LineStyle',':','LineWidth',line_width,'Visible','on');
p1=plot(100.*aw_series,AIOMFAC_aPsoa,'Color',color_plot, 'DisplayName' ,...
    ['AIOMFAC aPsoa run ' num2str(AIOMFAC_aPsoa_run)],'LineStyle','--','LineWidth',line_width,'Visible','on');

p2=plot(100.*aw_series_VBSBAT,VBSBAT_aPsoa_C_OA_PM,'Color',color_plot, 'DisplayName' ,['VBSBAT ' num2str(AIOMFAC_aPsoa_run)],...
    'LineStyle','-','LineWidth',line_width+1,'Visible','on');

color_plot=[1 0 1]; % Isoprene SOA
p1=plot(100.*aw_series_VBSBAT,VBSBAT_IPsoa_Mean_C_OA_PM,'Color',color_plot.*.8, 'DisplayName' ,...
    ['VBSBAT aPsoa MeanProp run ' num2str(AIOMFAC_IPsoa_run)],'LineStyle','-','LineWidth',line_width,'Visible','on');
p1=plot(100.*aw_series,AIOMFAC_IPsoa_dry,'Color',color_plot, 'DisplayName' ,...
    ['AIOMFAC IPsoa Dry run ' num2str(AIOMFAC_aPsoa_run)],'LineStyle','-.','LineWidth',line_width,'Visible','on');
p1=plot(100.*aw_series,AIOMFAC_IPsoa_NHSO4,'Color',color_plot, 'DisplayName' ,...
    ['AIOMFAC IPsoa with (NH_4)_2SO_4 run ' num2str(AIOMFAC_IPsoa_run_NHSO4)],'LineStyle',':','LineWidth',line_width,'Visible','on');
p1=plot(100.*aw_series,AIOMFAC_IPsoa,'Color',color_plot, 'DisplayName' ,...
    ['AIOMFAC IPsoa run ' num2str(AIOMFAC_IPsoa_run)],'LineStyle','--','LineWidth',line_width,'Visible','on');
p2=plot(100.*aw_series_VBSBAT,VBSBAT_IPsoa_C_OA_PM,'Color',color_plot, 'DisplayName' ,['VBSBAT ' num2str(AIOMFAC_IPsoa_run)],...
    'LineStyle','-','LineWidth',line_width+1,'Visible','on');

ylim([6,18.5])
ylabel('Organic PM Mass (\mug/m^3)')
% xlabel('Bulk Relative Humidity (%)  (a_w x100) ')

axes2 = axes('Parent',figure1,...
    'Position',[0.14 0.14 0.8 0.30],...
    'LineWidth',1.75,...
    'FontSize',fontsize);  %'CLim',[80 300],
set(axes2,'FontSize',fontsize,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% view(axes1,[175 35]);
grid(axes2,'on');
hold(axes2,'all');
hold on

%base line
plot([0;100],zeros(2,1),'Color',[0 0 0],'LineWidth',line_width,'LineStyle','--')

color_plot=[0 0.498039215803146 0]; %alpha pinene
p1=plot(100.*aw_series,100*(VBSBAT_aPsoa_Mean_C_OA_PM(added_aw:end)-AIOMFAC_aPsoa)./AIOMFAC_aPsoa,'Color',color_plot.*.8,...
    'DisplayName' ,['aPsoa Mean Org. run ' num2str(AIOMFAC_aPsoa_run)],...
    'LineStyle','-','LineWidth',line_width,'Visible','on');
p1=plot(100.*aw_series,100*(AIOMFAC_aPsoa_dry-AIOMFAC_aPsoa)./AIOMFAC_aPsoa,'Color',color_plot,...
    'DisplayName' ,['aPsoa with Dry run ' num2str(AIOMFAC_aPsoa_run)],...
    'LineStyle','-.','LineWidth',line_width,'Visible','on');
% p1=plot(100.*aw_series,100*(VBSBAT_aPsoa_C_OA_PM-AIOMFAC_aPsoa_NHSO4)./AIOMFAC_aPsoa_NHSO4,'Color',color_plot, ...
%     'DisplayName' ,['aPsoa with (NH_4)_2SO_4 run ' num2str(AIOMFAC_aPsoa_run)],...
%     'LineStyle',':','LineWidth',line_width,'Visible','off');
p2=plot(100.*aw_series,100*(VBSBAT_aPsoa_C_OA_PM(added_aw:end)-AIOMFAC_aPsoa)./AIOMFAC_aPsoa,'Color',color_plot, ...
    'DisplayName' ,['aPsoa Org. run ' num2str(AIOMFAC_aPsoa_run)],...
    'LineStyle','-','LineWidth',line_width+1,'Visible','on');

color_plot=[1 0 1]; % Isoprene SOA
p1=plot(100.*aw_series,100*(VBSBAT_IPsoa_Mean_C_OA_PM(added_aw:end)-AIOMFAC_IPsoa)./AIOMFAC_IPsoa,'Color',color_plot.*.8, ...
    'DisplayName' ,['IPsoa Mean Org. run ' num2str(AIOMFAC_aPsoa_run)],...
    'LineStyle','-','LineWidth',line_width,'Visible','on');
p1=plot(100.*aw_series,100*(AIOMFAC_IPsoa_dry-AIOMFAC_IPsoa)./AIOMFAC_IPsoa,'Color',color_plot, ...
    'DisplayName' ,['IPsoa with Dry run ' num2str(AIOMFAC_aPsoa_run)],...
    'LineStyle','-.','LineWidth',line_width,'Visible','on');
% p1=plot(100.*aw_series,100*(VBSBAT_IPsoa_C_OA_PM-AIOMFAC_IPsoa_NHSO4)./AIOMFAC_IPsoa_NHSO4,'Color',color_plot, ...
%     'DisplayName' ,['IPsoa with (NH_4)_2SO_4 run ' num2str(AIOMFAC_IPsoa_run)],...
%     'LineStyle',':','LineWidth',line_width,'Visible','off');
p2=plot(100.*aw_series,100*(VBSBAT_IPsoa_C_OA_PM(added_aw:end)-AIOMFAC_IPsoa)./AIOMFAC_IPsoa,'Color',color_plot, ...
    'DisplayName' ,['IPsoa Org. run ' num2str(AIOMFAC_IPsoa_run)],...
    'LineStyle','-','LineWidth',line_width+1,'Visible','on');

ylabel('Percent Difference (%)')
xlabel('Bulk Relative Humidity (RH% = a_w x100) ')
ylim([-10,10])

linkaxes([axes2, axes1],'x')

saveas(figure1, [plot_name '.png'], 'png')
saveas(figure1, [plot_name '.fig'], 'fig')




%% kappa

plot_name=['water and kappa plot', extension_name];
paper_postion=[0, 0, 3.5, 3.5*1.25].*paper_mult;
figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', paper_postion.*1+.5,'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]


axes1 = axes('Parent',figure1,...
    'Position',[0.14 0.58 0.8 0.4],...
    'LineWidth',1.75,...
    'FontSize',fontsize, 'Xcolor', 'none');  %'CLim',[80 300],
set(axes1,'FontSize',fontsize,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% view(axes1,[175 35]);
grid(axes1,'on');
hold(axes1,'all');
hold on

%base line
% plot([0;100],zeros(2,1),'Color',[0 0 0],'LineWidth',1.75)

color_plot=[0 0.498039215803146 0]; %alpha pinene
p1=plot(100.*aw_series_VBSBAT,VBSBAT_aPsoa_Mean_Caq_PM,'Color',color_plot.*.8, 'DisplayName' ,...
    ['VBSBAT aPsoa Mean run ' num2str(AIOMFAC_aPsoa_run)],...
    'LineStyle','-','LineWidth',line_width,'Visible','on');
p1=plot(100.*aw_series,AIOMFAC_aPsoa_water_NHSO4,'Color',color_plot, 'DisplayName' ,...
    ['AIOMFAC aPsoa with (NH_4)_2SO_4 run ' num2str(AIOMFAC_aPsoa_run_NHSO4)],...
    'LineStyle',':','LineWidth',line_width,'Visible','on');
p1=plot(100.*aw_series,AIOMFAC_aPsoa_water,'Color',color_plot, 'DisplayName' ,...
    ['AIOMFAC aPsoa run ' num2str(AIOMFAC_aPsoa_run)],'LineStyle','--','LineWidth',line_width,'Visible','on');

p2=plot(100.*aw_series_VBSBAT,VBSBAT_aPsoa_Caq_PM,'Color',color_plot, 'DisplayName' ,['VBSBAT ' num2str(AIOMFAC_aPsoa_run)],...
    'LineStyle','-','LineWidth',line_width+1,'Visible','on');

color_plot=[1 0 1]; % Isoprene SOA
p1=plot(100.*aw_series_VBSBAT,VBSBAT_IPsoa_Mean_Caq_PM,'Color',color_plot.*.8, 'DisplayName' ,...
    ['VBSBAT IPsoa Mean run ' num2str(AIOMFAC_IPsoa_run)],...
    'LineStyle','-','LineWidth',line_width,'Visible','on');
p1=plot(100.*aw_series,AIOMFAC_IPsoa_water_NHSO4,'Color',color_plot, 'DisplayName' ,...
    ['AIOMFAC IPsoa with (NH_4)_2SO_4 run ' num2str(AIOMFAC_IPsoa_run_NHSO4)],...
    'LineStyle',':','LineWidth',line_width,'Visible','on');
p1=plot(100.*aw_series,AIOMFAC_IPsoa_water,'Color',color_plot, 'DisplayName' ,...
    ['AIOMFAC IPsoa run ' num2str(AIOMFAC_IPsoa_run)],'LineStyle','--',...
    'LineWidth',line_width,'Visible','on');
p2=plot(100.*aw_series_VBSBAT,VBSBAT_IPsoa_Caq_PM,'Color',color_plot, 'DisplayName' ,['VBSBAT ' num2str(AIOMFAC_IPsoa_run)],...
    'LineStyle','-','LineWidth',line_width+1,'Visible','on');
ylim([0,30])

ylabel('Water PM Mass (\mug/m^3)')

axes2 = axes('Parent',figure1,...
    'Position',[0.14 0.14 0.8 0.4],...
    'LineWidth',1.75,...
    'FontSize',fontsize);  %'CLim',[80 300],
set(axes2,'FontSize',fontsize,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% view(axes1,[175 35]);
grid(axes2,'on');
hold(axes2,'all');
hold on

color_plot=[0 0.498039215803146 0]; %alpha pinene
p1=plot(100.*aw_series_VBSBAT,VBSBAT_aPsoa_Mean_kappaHGF,'Color',color_plot.*.8, 'DisplayName' ,...
    ['VBSBAT aPsoa Mean run ' num2str(AIOMFAC_aPsoa_run)],'LineStyle','-','LineWidth',line_width,'Visible','on');
p1=plot(100.*aw_series,AIOMFAC_aPsoa_kappaHGF,'Color',color_plot, 'DisplayName' ,...
    ['AIOMFAC aPsoa run ' num2str(AIOMFAC_aPsoa_run)],'LineStyle','--','LineWidth',line_width,'Visible','on');

p2=plot(100.*aw_series_VBSBAT,VBSBAT_aPsoa_kappaHGF,'Color',color_plot, 'DisplayName' ,['VBSBAT ' num2str(AIOMFAC_aPsoa_run)],...
    'LineStyle','-','LineWidth',line_width+1,'Visible','on');

color_plot=[1 0 1]; % Isoprene SOA
p1=plot(100.*aw_series_VBSBAT,VBSBAT_IPsoa_Mean_kappaHGF,'Color',color_plot.*.8, 'DisplayName' ,...
    ['VBSBAT IPsoa Mean run ' num2str(AIOMFAC_IPsoa_run)],'LineStyle','-','LineWidth',line_width,'Visible','on');
p1=plot(100.*aw_series,AIOMFAC_IPsoa_kappaHGF,'Color',color_plot, 'DisplayName' ,...
    ['AIOMFAC IPsoa run ' num2str(AIOMFAC_IPsoa_run)],'LineStyle','--','LineWidth',line_width,'Visible','on');
p2=plot(100.*aw_series_VBSBAT,VBSBAT_IPsoa_kappaHGF,'Color',color_plot, 'DisplayName' ,['VBSBAT ' num2str(AIOMFAC_IPsoa_run)],...
    'LineStyle','-','LineWidth',line_width+1,'Visible','on');


ylabel('\kappa_{HGF}')
xlabel('Bulk Relative Humidity (RH% = a_w x100) ')

linkaxes([axes2, axes1],'x')

saveas(figure1, [plot_name '.png'], 'png')
saveas(figure1, [plot_name '.fig'], 'fig')

%% zoom out on mass error
%% org mass
plot_name=['SI C_OA error', extension_name];
paper_postion=[0, 0, 3.5, 3.5].*paper_mult;
figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', paper_postion.*1+.5,'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]


axes2 = axes('Parent',figure1,...
    'Position',[0.16 0.16 0.75 0.75],...
    'LineWidth',1.75,...
    'FontSize',fontsize ); %'CLim',[80 300],
set(axes2,'FontSize',fontsize,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% view(axes1,[175 35]);
grid(axes2,'on');
hold(axes2,'all');
hold on

%base line
plot([0;100],zeros(2,1),'Color',[0 0 0],'LineWidth',line_width,'LineStyle','--')

color_plot=[0 0.498039215803146 0]; %alpha pinene
p1=plot(100.*aw_series,100*(VBSBAT_aPsoa_Mean_C_OA_PM(added_aw:end)-AIOMFAC_aPsoa)./AIOMFAC_aPsoa,'Color',color_plot.*.8,...
    'DisplayName' ,['aPsoa Mean Org. run ' num2str(AIOMFAC_aPsoa_run)],...
    'LineStyle','-','LineWidth',line_width,'Visible','on');
p1=plot(100.*aw_series,100*(AIOMFAC_aPsoa_dry-AIOMFAC_aPsoa)./AIOMFAC_aPsoa,'Color',color_plot,...
    'DisplayName' ,['aPsoa with Dry run ' num2str(AIOMFAC_aPsoa_run)],...
    'LineStyle','-.','LineWidth',line_width,'Visible','on');
% p1=plot(100.*aw_series,100*(VBSBAT_aPsoa_C_OA_PM-AIOMFAC_aPsoa_NHSO4)./AIOMFAC_aPsoa_NHSO4,'Color',color_plot, ...
%     'DisplayName' ,['aPsoa with (NH_4)_2SO_4 run ' num2str(AIOMFAC_aPsoa_run)],...
%     'LineStyle',':','LineWidth',line_width,'Visible','off');
p2=plot(100.*aw_series,100*(VBSBAT_aPsoa_C_OA_PM(added_aw:end)-AIOMFAC_aPsoa)./AIOMFAC_aPsoa,'Color',color_plot, ...
    'DisplayName' ,['aPsoa Org. run ' num2str(AIOMFAC_aPsoa_run)],...
    'LineStyle','-','LineWidth',line_width+1,'Visible','on');

color_plot=[1 0 1]; % Isoprene SOA
p1=plot(100.*aw_series,100*(VBSBAT_IPsoa_Mean_C_OA_PM(added_aw:end)-AIOMFAC_IPsoa)./AIOMFAC_IPsoa,'Color',color_plot.*.8, ...
    'DisplayName' ,['IPsoa Mean Org. run ' num2str(AIOMFAC_aPsoa_run)],...
    'LineStyle','-','LineWidth',line_width,'Visible','on');
p1=plot(100.*aw_series,100*(AIOMFAC_IPsoa_dry-AIOMFAC_IPsoa)./AIOMFAC_IPsoa,'Color',color_plot, ...
    'DisplayName' ,['IPsoa with Dry run ' num2str(AIOMFAC_aPsoa_run)],...
    'LineStyle','-.','LineWidth',line_width,'Visible','on');
% p1=plot(100.*aw_series,100*(VBSBAT_IPsoa_C_OA_PM-AIOMFAC_IPsoa_NHSO4)./AIOMFAC_IPsoa_NHSO4,'Color',color_plot, ...
%     'DisplayName' ,['IPsoa with (NH_4)_2SO_4 run ' num2str(AIOMFAC_IPsoa_run)],...
%     'LineStyle',':','LineWidth',line_width,'Visible','off');
p2=plot(100.*aw_series,100*(VBSBAT_IPsoa_C_OA_PM(added_aw:end)-AIOMFAC_IPsoa)./AIOMFAC_IPsoa,'Color',color_plot, ...
    'DisplayName' ,['IPsoa Org. run ' num2str(AIOMFAC_IPsoa_run)],...
    'LineStyle','-','LineWidth',line_width+1,'Visible','on');

ylabel('Percent Difference (%)')
xlabel('Bulk Relative Humidity (RH% = a_w x100) ')
ylim([-70,70])

% linkaxes([axes2, axes1],'x')

saveas(figure1, [plot_name '.png'], 'png')
saveas(figure1, [plot_name '.fig'], 'fig')

