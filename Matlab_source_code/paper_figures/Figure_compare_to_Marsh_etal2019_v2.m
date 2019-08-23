 %%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2019-Feb-05 10:32 AM
% Copyright 2019 Kyle Gorkowski 
%%
% compare data to
%Marsh, A., Rovelli, G., Miles, R. E. H. and Reid, J. P.: The Complexity of Measuring and Representing the Hygroscopicity of Mixed Component Aerosol, J. Phys. Chem. A, acs.jpca.8b11623, doi:10.1021/acs.jpca.8b11623, 2019.
%
clear
 BAT_refinement_mode='interpolate'; % 'perfect water activity' 'none' 'interpolate'

VBSBAT_options=default_VBSBAT_options('robust');%('robust');
% VBSBAT_options.BAT_refinement_aw=0;
 VBSBAT_options.plot_PM='no';

% 
paper_mult=1;
fontsize=10;

Mlength=50;
MFS_aw=zeros(Mlength,7);
MFS_OT=zeros(Mlength,7);
MFS_VBSBAT=zeros(Mlength,7);

%% mixture1 calc
extension_name='mixture1';

aw_series=[0.5: 0.001: 1]';

mixture1_compounds={'glycine', 'lysine', 'glutaric acid', 'malonic acid'}';
mixture1_functional_group={'carboxyl', 'carboxyl','carboxyl', 'carboxyl'}'; % note carboxyl is identical to hydroxyl in the BAT model
mixture1_total_moles= [1, 1, 1, 1]';
mixture1_Psat= mixture1_total_moles.*10^-9;

% {'hydroxyl';'carboxyl';'ether';'ketone';'ester';'hydroperoxide';'hydroperoxideSOA';'PEG'};

mixture1_H2C=[5/2, 14/6, 8/5, 4/3 ]' ;
mixture1_O2C=[2/2, 2/6, 4/5, 4/3]' ;
mixture1_N2C=[1/2, 2/6, 0/5, 0/3]' ;
mixture1_MolarMass=[75, 146.19, 132.12, 104.06 ]' ;

[Csat_ugPm3, Csat_ugPm3_log10] = Csat_calculate_v1(mixture1_MolarMass,mixture1_Psat, 298, 'no');

% mixture normalized mole fractions
Org_Conc_OM=mixture1_total_moles.*mixture1_MolarMass./(sum(mixture1_total_moles.*mixture1_MolarMass));

[VBSBAT_mix1_C_OA_PM, VBSBAT_mix1_Caq_PM, VBSBAT_mix1_kappaHGF, details_mix1]=VBS_BAT_simulation_v2(...
    Csat_ugPm3, Org_Conc_OM, ...
    mixture1_O2C, mixture1_H2C,  mixture1_MolarMass,...
    aw_series, mixture1_functional_group, BAT_refinement_mode, VBSBAT_options, 'mixture 1', mixture1_N2C);

disp(['fit per aw ' num2str(details_mix1.fit.mean_time_per_fit) ' sec'])



mixture_data=[0.5448	0.00158	0.00189	0.79691	0.0015
0.56824	0.00944	0.01128	0.76741	0.0055
0.58983	0.00145	0.00176	0.78525	0.00106
0.61273	0.00537	0.00647	0.75369	0.00535
0.63077	0.00151	0.00184	0.76666	6.29E-04
0.65171	0.00426	0.00517	0.74069	0.00482
0.67405	0.00134	0.00164	0.73727	5.18E-04
0.69153	0.00404	0.00494	0.71552	0.007
0.71652	0.00116	0.00144	0.70497	7.21E-04
0.7309	0.00335	0.00413	0.68904	0.0041
0.7506	0.00402	0.00497	0.67359	0.00276
0.77141	0.00626	0.00777	0.65536	0.00645
0.79414	7.04E-04	8.92E-04	0.6274	6.90E-04
0.80673	0.001	5.58E-04	0.60803	7.59E-04
0.83199	0.00205	0.0017	0.5723	0.00299
0.85703	6.58E-04	4.17E-04	0.52364	9.44E-04
0.87007	0.00111	9.88E-04	0.49849	0.0027
0.8916	0.0014	0.00131	0.44744	0.00372
0.91212	0.00137	0.00122	0.39599	0.00622
0.93242	0.00126	0.00104	0.33353	0.00638
0.95313	0.00102	8.77E-04	0.24709	0.00484
0.97007	0.00101	6.97E-04	0.18316	0.00403];

kappa_data=[0.85632	0.00188	0.20684	0.00444
0.86436	0.00263	0.20471	0.00626
0.87566	0.00228	0.20511	0.00614
0.8847	0.00266	0.20673	0.00513
0.8951	0.00297	0.20761	0.00431
0.90552	0.00307	0.20869	0.00402
0.9157	0.00254	0.20996	0.00399
0.92529	0.00299	0.21171	0.00385
0.93554	0.00307	0.21269	0.00452
0.94573	0.00291	0.21295	0.00457
0.95577	0.00308	0.21491	0.00558
0.96636	0.0029	0.21067	0.00631
0.97255	0.00115	0.20794	0.00586];

%save MFS matrix
MFS_aw(1:size(mixture_data,1),1)=mixture_data(:,1);
MFS_OT(1:size(mixture_data,1),1)=mixture_data(:,4);
MFS_VBSBAT(1:size(mixture_data,1),1)=interp1(aw_series,VBSBAT_mix1_C_OA_PM./(VBSBAT_mix1_Caq_PM+VBSBAT_mix1_C_OA_PM),mixture_data(:,1),'linear');
% figure
% plot(MFS_OT,MFS_VBSBAT)

%plot
plot_name=['water MFS compare ', extension_name];
paper_postion=[0, 0, 3.5, 3.5*1.25].*paper_mult;
figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', paper_postion.*1+.5,'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]


axes1 = axes('Parent',figure1,...
    'Position',[0.16 0.16 0.75 0.75],...
    'LineWidth',1.75,...
    'FontSize',fontsize );  %'CLim',[80 300],
set(axes1,'FontSize',fontsize,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% view(axes1,[175 35]);
grid(axes1,'on');
hold(axes1,'all');
hold on

plot(aw_series,VBSBAT_mix1_C_OA_PM./(VBSBAT_mix1_Caq_PM+VBSBAT_mix1_C_OA_PM),'LineWidth',2)

fill(axes1,[mixture_data(:,1)+mixture_data(:,2);flip(mixture_data(:,1)-mixture_data(:,3))],[mixture_data(:,4)+mixture_data(:,5);flip(mixture_data(:,4)-mixture_data(:,5))],'r','FaceAlpha',0.25,...
    'FaceColor',[0.63921568627451 0.07843137254902 0.180392156862745],...
    'EdgeColor','none');
    
plot(mixture_data(:,1),mixture_data(:,4),'or')

plot(aw_series, VBSBAT_mix1_kappaHGF,'k','LineWidth',2)
plot(kappa_data(:,1),kappa_data(:,3),'ok')
ylabel('MFS or \kappa')
xlabel('water activity (a_w)')
title(extension_name)
annotation(figure1,'textbox',...
    [0.214541666666667 0.385897435897436 0.347958333333333 0.235042735042737],...
    'String',{char(mixture1_compounds)},...
    'FitBoxToText','off','LineStyle','none');
xlim([0.5,1]) 
ylim([0,1])
saveas(figure1, [plot_name '.png'], 'png')
saveas(figure1, [plot_name '.fig'], 'fig')


%% %% mixture2 calc
extension_name='mixture2';

aw_series=[0.5: 0.001: 1]';

mixture2_compounds={  'methyl succinic acid', 'arginine','glutaric acid', 'citric acid'}';
mixture2_functional_group={'carboxyl', 'carboxyl','carboxyl', 'carboxyl'}'; % note carboxyl is identical to hydroxyl in the BAT model
mixture2_total_moles= [1, 1, 1, 1]';
mixture2_Psat= mixture1_total_moles.*10^-9;

% {'hydroxyl';'carboxyl';'ether';'ketone';'ester';'hydroperoxide';'hydroperoxideSOA';'PEG'};

mixture2_H2C=[8/5, 14/6, 8/5, 8/6 ]' ;
mixture2_O2C=[4/5, 2/6, 4/5, 7/6]' ;
mixture2_N2C=[0/5, 4/6, 0/5, 0/6]' ;
mixture2_MolarMass=[132.12, 174.2, 132.12, 192.12 ]' ;

[Csat_ugPm3, Csat_ugPm3_log10] = Csat_calculate_v1(mixture2_MolarMass,mixture2_Psat, 298, 'no');

% mixture normalized mole fractions
Org_Conc_OM=mixture2_total_moles.*mixture2_MolarMass./(sum(mixture2_total_moles.*mixture2_MolarMass));

[VBSBAT_mix2_C_OA_PM, VBSBAT_mix2_Caq_PM, VBSBAT_mix2_kappaHGF, details_mix2]=VBS_BAT_simulation_v2(...
    Csat_ugPm3, Org_Conc_OM, ...
    mixture2_O2C, mixture2_H2C,  mixture2_MolarMass, aw_series,...
    mixture2_functional_group, BAT_refinement_mode, VBSBAT_options,extension_name, mixture2_N2C );

% disp(['fit per aw ' num2str(details_mix1.fit.mean_time_per_fit) ' sec'])

mixture_data=[0.54895	0.00196	0.00234	0.87207	0.00319
0.58698	0.00505	0.00609	0.85121	0.01173
0.60151	0.00245	0.00296	0.85376	0.00535
0.64	0.0025	0.00305	0.84976	0.00679
0.65144	0.00257	0.00313	0.8464	0.00444
0.69133	0.0023	0.00284	0.84717	0.00477
0.70286	0.00246	0.00304	0.83422	0.00838
0.72897	0.00593	0.00729	0.78434	0.02006
0.75906	0.00696	0.00855	0.75759	0.01401
0.7817	0.01008	0.01243	0.74846	0.01335
0.81649	0.00389	0.002	0.72851	0.00178
0.83553	0.00233	0.0017	0.70511	0.00248
0.85482	0.00137	0.00109	0.68312	0.0022
0.88236	0.00333	0.0033	0.63106	0.00849
0.90755	6.28E-04	6.59E-04	0.61525	0.00414
0.93285	0.00187	0.00187	0.50221	0.01232
0.95965	0.00135	0.00134	0.34044	0.01191
0.98491	7.05E-04	7.03E-04	0.14535	0.00559
0.99812	0.01197	0.01163	0.09018	0.00515];

kappa_data=[0.90595	0.00132	0.08486	0.00278
0.9159	0.00256	0.08487	0.00279
0.92648	0.0032	0.08536	0.00212
0.93519	0.00234	0.08613	0.00171
0.94521	0.0032	0.08725	0.00205
0.95617	0.00282	0.08969	0.002
0.96581	0.00288	0.09505	0.0027
0.97615	0.00292	0.10078	0.00313
0.98623	0.00281	0.10593	0.009
0.99245	0.00147	0.09508	0.01635];

%save MFS matrix
MFS_aw(1:size(mixture_data,1),2)=mixture_data(:,1);
MFS_OT(1:size(mixture_data,1),2)=mixture_data(:,4);
MFS_VBSBAT(1:size(mixture_data,1),2)=interp1(aw_series,VBSBAT_mix2_C_OA_PM./(VBSBAT_mix2_Caq_PM+VBSBAT_mix2_C_OA_PM),mixture_data(:,1),'linear');
%
plot_name=['water MFS compare ', extension_name];
paper_postion=[0, 0, 3.5, 3.5*1.25].*paper_mult;
figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', paper_postion.*1+.5,'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]


axes1 = axes('Parent',figure1,...
    'Position',[0.16 0.16 0.75 0.75],...
    'LineWidth',1.75,...
    'FontSize',fontsize );  %'CLim',[80 300],
set(axes1,'FontSize',fontsize,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% view(axes1,[175 35]);
grid(axes1,'on');
hold(axes1,'all');
hold on

plot(aw_series,VBSBAT_mix2_C_OA_PM./(VBSBAT_mix2_Caq_PM+VBSBAT_mix2_C_OA_PM),'LineWidth',2)

fill(axes1,[mixture_data(:,1)+mixture_data(:,2);flip(mixture_data(:,1)-mixture_data(:,3))],[mixture_data(:,4)+mixture_data(:,5);flip(mixture_data(:,4)-mixture_data(:,5))],'r','FaceAlpha',0.25,...
    'FaceColor',[0.63921568627451 0.07843137254902 0.180392156862745],...
    'EdgeColor','none');
    
plot(mixture_data(:,1),mixture_data(:,4),'or')


plot(aw_series, VBSBAT_mix2_kappaHGF,'k','LineWidth',2)
plot(kappa_data(:,1),kappa_data(:,3),'ok')
ylabel('MFS or \kappa')
xlabel('water activity (a_w)')
title(extension_name)
annotation(figure1,'textbox',...
    [0.214541666666667 0.385897435897436 0.347958333333333 0.235042735042737],...
    'String',{char(mixture2_compounds)},...
    'FitBoxToText','off','LineStyle','none');
xlim([0.5,1]) 
ylim([0,1])
saveas(figure1, [plot_name '.png'], 'png')
saveas(figure1, [plot_name '.fig'], 'fig')

%% mixture3 calc
extension_name='mixture3';

aw_series=[0.5: 0.001: 1]';

mixture3_compounds={'oxalic acid', 'malonic acid'}';
mixture3_functional_group={'carboxyl', 'carboxyl'}'; % note carboxyl is identical to hydroxyl in the BAT model
mixture3_total_moles= [1, 1]';
mixture3_Psat= mixture3_total_moles.*10^-9;

% {'hydroxyl';'carboxyl';'ether';'ketone';'ester';'hydroperoxide';'hydroperoxideSOA';'PEG'};

mixture3_H2C=[2/2, 4/3 ]' ;
mixture3_O2C=[4/2, 4/3]' ;
mixture3_N2C=[0/2, 0/3]' ;
mixture3_MolarMass=[90.03, 104.06]' ;

[Csat_ugPm3, Csat_ugPm3_log10] = Csat_calculate_v1(mixture3_MolarMass,mixture3_Psat, 298, 'no');

% mixture normalized mole fractions
Org_Conc_OM=mixture3_total_moles.*mixture3_MolarMass./(sum(mixture3_total_moles.*mixture3_MolarMass));

[VBSBAT_mix3_C_OA_PM, VBSBAT_mix3_Caq_PM, VBSBAT_mix3_kappaHGF, details_mix3]=VBS_BAT_simulation_v2(...
    Csat_ugPm3, Org_Conc_OM, ...
    mixture3_O2C, mixture3_H2C,  mixture3_MolarMass, aw_series, mixture3_functional_group, ...
    BAT_refinement_mode, VBSBAT_options, extension_name, mixture3_N2C );

% disp(['fit per aw ' num2str(details_mix1.fit.mean_time_per_fit) ' sec'])

mixture_data=[0.81668	0.00123	6.04E-04	0.47931	5.06E-04
0.82871	0.00676	0.00332	0.46542	0.00189
0.83934	0.00499	0.00245	0.45156	0.00289
0.85029	0.00676	0.00332	0.43463	0.00383
0.85893	0.01012	0.00497	0.42159	0.00677
0.86655	0.00796	0.00391	0.40876	0.00312
0.87857	0.00676	0.00332	0.38807	0.0033
0.89137	0.00749	0.00368	0.36511	0.00222
0.90113	0.01012	0.00497	0.34507	0.00766
0.90654	0.00922	0.00453	0.33328	0.0069
0.91586	0.0071	0.00349	0.31151	0.00321
0.92842	0.00542	0.00266	0.28017	0.00603
0.9394	0.00513	0.00252	0.24816	0.00553
0.94997	0.00559	0.00275	0.21661	0.0048
0.95967	0.00388	0.00191	0.18498	0.00481
0.96991	0.00348	0.00171	0.14691	0.00369
0.9798	0.00268	0.00132	0.10946	0.0031
0.98601	0.00487	0.00239	0.0882	0.00143];

kappa_data=[0.81632	0.00128	0.3836	0.00326
0.8242	0.00322	0.38076	0.00551
0.83637	0.002	0.36756	0.00231
0.84393	0.00133	0.36453	0.00233
0.85393	0.00167	0.35807	0.00212
0.86587	0.00168	0.35073	0.00209
0.878	0.00129	0.34292	0.00175
0.88683	0.00417	0.33817	0.00608
0.89264	0.00143	0.334	0.0021
0.90453	0.00177	0.3258	0.0024
0.91567	0.0015	0.3194	0.00282
0.92607	0.0011	0.31273	0.00246
0.9346	8.28E-04	0.30793	0.00301
0.94537	0.00353	0.30167	0.00532
0.95673	0.00252	0.29207	0.00478
0.96698	0.00248	0.28361	0.0133
0.97641	0.00285	0.26608	0.01082
0.98397	0.00202	0.24103	0.01582];
%save MFS matrix
MFS_aw(1:size(mixture_data,1),3)=mixture_data(:,1);
MFS_OT(1:size(mixture_data,1),3)=mixture_data(:,4);
MFS_VBSBAT(1:size(mixture_data,1),3)=interp1(aw_series,VBSBAT_mix3_C_OA_PM./(VBSBAT_mix3_Caq_PM+VBSBAT_mix3_C_OA_PM),mixture_data(:,1),'linear');
%
plot_name=['water MFS compare ', extension_name];
paper_postion=[0, 0, 3.5, 3.5*1.25].*paper_mult;
figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', paper_postion.*1+.5,'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]


axes1 = axes('Parent',figure1,...
    'Position',[0.16 0.16 0.75 0.75],...
    'LineWidth',1.75,...
    'FontSize',fontsize );  %'CLim',[80 300],
set(axes1,'FontSize',fontsize,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% view(axes1,[175 35]);
grid(axes1,'on');
hold(axes1,'all');
hold on

plot(aw_series,VBSBAT_mix3_C_OA_PM./(VBSBAT_mix3_Caq_PM+VBSBAT_mix3_C_OA_PM),'LineWidth',2)

fill(axes1,[mixture_data(:,1)+mixture_data(:,2);flip(mixture_data(:,1)-mixture_data(:,3))],[mixture_data(:,4)+mixture_data(:,5);flip(mixture_data(:,4)-mixture_data(:,5))],'r','FaceAlpha',0.25,...
    'FaceColor',[0.63921568627451 0.07843137254902 0.180392156862745],...
    'EdgeColor','none');
    
plot(mixture_data(:,1),mixture_data(:,4),'or')

plot(aw_series, VBSBAT_mix3_kappaHGF,'k','LineWidth',2)
plot(kappa_data(:,1),kappa_data(:,3),'ok')
ylabel('MFS or \kappa')
xlabel('water activity (a_w)')
title(extension_name)
annotation(figure1,'textbox',...
    [0.214541666666667 0.385897435897436 0.347958333333333 0.235042735042737],...
    'String',{char(mixture3_compounds)},...
    'FitBoxToText','off','LineStyle','none');
xlim([0.5,1]) 
ylim([0,1])
saveas(figure1, [plot_name '.png'], 'png')
saveas(figure1, [plot_name '.fig'], 'fig')


%% mixture4 calc
extension_name='mixture4';

aw_series=[0.5: 0.001: 1]';

mixture4_compounds={'glycine', 'lysine'}';
mixture4_functional_group={'carboxyl', 'carboxyl'}'; % note carboxyl is identical to hydroxyl in the BAT model
mixture4_total_moles= [1, 1]';
mixture4_Psat= mixture4_total_moles.*10^-9;

% {'hydroxyl';'carboxyl';'ether';'ketone';'ester';'hydroperoxide';'hydroperoxideSOA';'PEG'};

mixture4_H2C=[5/2, 14/6 ]' ;
mixture4_O2C=[2/2, 2/6]' ;
mixture4_N2C=[1/2, 2/6]' ;

% mixture4_H2C=[5/3, 14/8 ]' ; % convert n to C
% mixture4_O2C=[2/3, 2/8]' ;
% mixture4_H2C=[5/2, 14/6 ]' ;% convert nx1=O
% mixture4_O2C=[3/2, 4/6]' ;
mixture4_MolarMass=[75.07, 146.19]' ;

[Csat_ugPm3, Csat_ugPm3_log10] = Csat_calculate_v1(mixture4_MolarMass,mixture4_Psat, 298, 'no');

% mixture normalized mole fractions
Org_Conc_OM=mixture4_total_moles.*mixture4_MolarMass./(sum(mixture4_total_moles.*mixture4_MolarMass));

[VBSBAT_mix4_C_OA_PM, VBSBAT_mix4_Caq_PM, VBSBAT_mix4_kappaHGF, details_mix3]=VBS_BAT_simulation_v2(...
    Csat_ugPm3, Org_Conc_OM, ...
    mixture4_O2C, mixture4_H2C,  mixture4_MolarMass, aw_series, mixture4_functional_group, ...
    BAT_refinement_mode, VBSBAT_options, extension_name, mixture4_N2C);

% disp(['fit per aw ' num2str(details_mix1.fit.mean_time_per_fit) ' sec'])

mixture_data=[0.83885	0.00122	6.75E-04	0.44529	6.43E-04
0.85074	0.00403	0.00223	0.42811	0.00188
0.86121	0.00523	0.00289	0.41403	0.00292
0.87116	0.00522	0.00288	0.39794	0.00411
0.88239	0.00573	0.00316	0.37969	0.00482
0.89382	0.00818	0.00451	0.36062	0.00567
0.90216	0.00573	0.00316	0.34347	0.00646
0.91324	0.00643	0.00355	0.32039	0.00603
0.92331	0.00573	0.00316	0.29846	0.00529
0.93241	0.00643	0.00355	0.27503	0.00158
0.94039	0.00605	0.00334	0.2538	0.00388
0.95133	0.00425	0.00235	0.2271	0.00496
0.96115	0.00329	0.00181	0.18808	0.0049
0.97138	0.00228	0.00126	0.1379	0.00461
0.97969	0.00209	0.00115	0.0943	0.00285];

kappa_data=[0.83807	0.00105	0.32636	0.00291
0.84371	0.003	0.3246	0.00331
0.8551	0.00304	0.3152	0.00297
0.86512	0.00262	0.30718	0.00178
0.87533	0.0013	0.29992	0.00151
0.88458	0.00162	0.29325	0.00142
0.89558	0.00116	0.28317	0.0014
0.90525	0.00136	0.27817	0.00147
0.91525	0.00129	0.27025	0.00171
0.92475	0.00114	0.26225	0.00249
0.93577	0.00374	0.25491	0.0038
0.94671	0.00258	0.242	0.00464
0.95704	0.00222	0.23164	0.00519
0.96609	0.00298	0.23614	0.00684
0.9755	0.00304	0.25834	0.01415
0.9815	5.99E-04	0.27258	0.01613];
    %save MFS matrix
MFS_aw(1:size(mixture_data,1),4)=mixture_data(:,1);
MFS_OT(1:size(mixture_data,1),4)=mixture_data(:,4);
MFS_VBSBAT(1:size(mixture_data,1),4)=interp1(aw_series,VBSBAT_mix4_C_OA_PM./(VBSBAT_mix4_Caq_PM+VBSBAT_mix4_C_OA_PM),mixture_data(:,1),'linear');
%
plot_name=['water MFS compare ', extension_name];
paper_postion=[0, 0, 3.5, 3.5*1.25].*paper_mult;
figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', paper_postion.*1+.5,'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]


axes1 = axes('Parent',figure1,...
    'Position',[0.16 0.16 0.75 0.75],...
    'LineWidth',1.75,...
    'FontSize',fontsize );  %'CLim',[80 300],
set(axes1,'FontSize',fontsize,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% view(axes1,[175 35]);
grid(axes1,'on');
hold(axes1,'all');
hold on

plot(aw_series,VBSBAT_mix4_C_OA_PM./(VBSBAT_mix4_Caq_PM+VBSBAT_mix4_C_OA_PM),'LineWidth',2)

fill(axes1,[mixture_data(:,1)+mixture_data(:,2);flip(mixture_data(:,1)-mixture_data(:,3))],[mixture_data(:,4)+mixture_data(:,5);flip(mixture_data(:,4)-mixture_data(:,5))],'r','FaceAlpha',0.25,...
    'FaceColor',[0.63921568627451 0.07843137254902 0.180392156862745],...
    'EdgeColor','none');
    
plot(mixture_data(:,1),mixture_data(:,4),'or')


plot(aw_series, VBSBAT_mix4_kappaHGF,'k','LineWidth',2)
plot(kappa_data(:,1),kappa_data(:,3),'ok')
ylabel('MFS or \kappa')
xlabel('water activity (a_w)')
title(extension_name)
annotation(figure1,'textbox',...
    [0.214541666666667 0.385897435897436 0.347958333333333 0.235042735042737],...
    'String',{char(mixture4_compounds)},...
    'FitBoxToText','off','LineStyle','none');
xlim([0.5,1]) 
ylim([0,1])
saveas(figure1, [plot_name '.png'], 'png')
saveas(figure1, [plot_name '.fig'], 'fig')

%% mixture5 calc
extension_name='mixture5';

aw_series=[0.5: 0.001: 1]';

mixture5_compounds={'glutaric acid', 'malonic acid'}';
mixture5_functional_group={'carboxyl', 'carboxyl'}'; % note carboxyl is identical to hydroxyl in the BAT model
mixture5_total_moles= [1, 1]';
mixture5_Psat= mixture5_total_moles.*10^-9;

% {'hydroxyl';'carboxyl';'ether';'ketone';'ester';'hydroperoxide';'hydroperoxideSOA';'PEG'};

mixture5_H2C=[8/5, 4/3 ]' ;
mixture5_O2C=[4/2, 4/3]' ;
mixture5_N2C=[0/2, 0/3]' ;
mixture5_MolarMass=[132.12, 104.06]' ;

[Csat_ugPm3, Csat_ugPm3_log10] = Csat_calculate_v1(mixture5_MolarMass,mixture5_Psat, 298, 'no');

% mixture normalized mole fractions
Org_Conc_OM=mixture5_total_moles.*mixture5_MolarMass./(sum(mixture5_total_moles.*mixture5_MolarMass));

[VBSBAT_mix5_C_OA_PM, VBSBAT_mix5_Caq_PM, VBSBAT_mix5_kappaHGF, details_mix3]=VBS_BAT_simulation_v2(...
    Csat_ugPm3, Org_Conc_OM, ...
    mixture5_O2C, mixture5_H2C,  mixture5_MolarMass, aw_series, mixture5_functional_group, ...
    BAT_refinement_mode, VBSBAT_options,  extension_name, mixture5_N2C);

% disp(['fit per aw ' num2str(details_mix1.fit.mean_time_per_fit) ' sec'])

mixture_data=[0.84163	0.00134	7.56E-04	0.55996	5.99E-04
0.85092	0.00354	0.002	0.54442	0.00248
0.86222	0.0062	0.00349	0.5284	0.0038
0.87216	0.00619	0.00349	0.50977	0.00553
0.88081	0.00786	0.00443	0.49303	0.00316
0.89212	0.00663	0.00373	0.46975	0.00326
0.90413	0.00663	0.00373	0.44226	0.00363
0.9155	0.00663	0.00373	0.41198	0.00352
0.92585	0.00663	0.00373	0.38148	0.00357
0.93449	0.00663	0.00373	0.35126	0.00307
0.94418	0.00466	0.00262	0.31257	0.00835
0.95366	0.00466	0.00262	0.27119	0.00875
0.96464	0.00312	0.00176	0.2197	0.008
0.97474	0.00207	0.00117	0.1483	0.00617
0.98127	0.00233	0.00132	0.10012	0.00307];

kappa_data=[0.83987	3.39E-04	0.19577	0.00245
0.84299	0.00266	0.1943	0.00235
0.85453	0.00324	0.19259	0.00272
0.86536	0.00327	0.18743	0.0015
0.87517	0.00338	0.18567	0.00115
0.884	0.00339	0.18322	8.33E-04
0.8935	0.00201	0.1801	0.00166
0.90467	0.0021	0.17742	0.00138
0.916	0.00191	0.17483	0.00175
0.92625	0.0016	0.17292	0.00178
0.93523	0.00196	0.17262	0.00185
0.94548	0.0032	0.17287	0.00234
0.95609	0.00292	0.17352	0.00304
0.96639	0.00228	0.17576	0.00822
0.97618	0.00279	0.2062	0.0189
0.98226	9.53E-04	0.23313	0.01108];

%save MFS matrix
MFS_aw(1:size(mixture_data,1),5)=mixture_data(:,1);
MFS_OT(1:size(mixture_data,1),5)=mixture_data(:,4);
MFS_VBSBAT(1:size(mixture_data,1),5)=interp1(aw_series,VBSBAT_mix5_C_OA_PM./(VBSBAT_mix5_Caq_PM+VBSBAT_mix5_C_OA_PM),mixture_data(:,1),'linear');
%

plot_name=['water MFS compare ', extension_name];
paper_postion=[0, 0, 3.5, 3.5*1.25].*paper_mult;
figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', paper_postion.*1+.5,'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]


axes1 = axes('Parent',figure1,...
    'Position',[0.16 0.16 0.75 0.75],...
    'LineWidth',1.75,...
    'FontSize',fontsize );  %'CLim',[80 300],
set(axes1,'FontSize',fontsize,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% view(axes1,[175 35]);
grid(axes1,'on');
hold(axes1,'all');
hold on

plot(aw_series,VBSBAT_mix5_C_OA_PM./(VBSBAT_mix5_Caq_PM+VBSBAT_mix5_C_OA_PM),'LineWidth',2)

fill(axes1,[mixture_data(:,1)+mixture_data(:,2);flip(mixture_data(:,1)-mixture_data(:,3))],[mixture_data(:,4)+mixture_data(:,5);flip(mixture_data(:,4)-mixture_data(:,5))],'r','FaceAlpha',0.25,...
    'FaceColor',[0.63921568627451 0.07843137254902 0.180392156862745],...
    'EdgeColor','none');
    
plot(mixture_data(:,1),mixture_data(:,4),'or')

plot(aw_series, VBSBAT_mix5_kappaHGF,'k','LineWidth',2)
plot(kappa_data(:,1),kappa_data(:,3),'ok')
ylabel('MFS or \kappa')
xlabel('water activity (a_w)')
title(extension_name)
annotation(figure1,'textbox',...
    [0.214541666666667 0.385897435897436 0.347958333333333 0.235042735042737],...
    'String',{char(mixture5_compounds)},...
    'FitBoxToText','off','LineStyle','none');
xlim([0.5,1]) 
ylim([0,1])
saveas(figure1, [plot_name '.png'], 'png')
saveas(figure1, [plot_name '.fig'], 'fig')

%% mixture6 calc
extension_name='mixture6';

aw_series=[0.5: 0.001: 1]';

mixture6_compounds={'glycine', 'lysine', 'glutaric acid', 'malonic acid'}';
mixture6_functional_group={'carboxyl', 'carboxyl','carboxyl', 'carboxyl'}'; % note carboxyl is identical to hydroxyl in the BAT model
mixture6_total_moles= [2, 2, 1, 1]';
mixture6_Psat= mixture6_total_moles.*10^-9;

% {'hydroxyl';'carboxyl';'ether';'ketone';'ester';'hydroperoxide';'hydroperoxideSOA';'PEG'};

mixture6_H2C=[5/2, 14/6, 8/5, 4/3 ]' ;
mixture6_O2C=[2/2, 2/6, 4/5, 4/3]' ;
mixture6_N2C=[1/2, 2/6, 0/5, 0/3]' ;

%  mixture6_O2C=[3/2, 4/6, 4/5, 4/3]' ; %% mixture4_O2C=[4/2, 6/6]' ;


mixture6_MolarMass=[75.07, 146.19, 132.12, 104.06]' ;

[Csat_ugPm3, Csat_ugPm3_log10] = Csat_calculate_v1(mixture6_MolarMass,mixture6_Psat, 298, 'no');

% mixture normalized mole fractions
Org_Conc_OM=mixture6_total_moles.*mixture6_MolarMass./(sum(mixture6_total_moles.*mixture6_MolarMass));

[VBSBAT_mix6_C_OA_PM, VBSBAT_mix6_Caq_PM, VBSBAT_mix6_kappaHGF, details_mix3]=VBS_BAT_simulation_v2(...
    Csat_ugPm3, Org_Conc_OM, ...
    mixture6_O2C, mixture6_H2C,  mixture6_MolarMass, aw_series, mixture6_functional_group,...
    BAT_refinement_mode, VBSBAT_options, extension_name, mixture6_N2C);

% disp(['fit per aw ' num2str(details_mix1.fit.mean_time_per_fit) ' sec'])

mixture_data=[0.85049	0.00111	6.55E-04	0.57391	7.41E-04
0.86289	0.00306	0.00181	0.55213	0.00217
0.87145	0.0034	0.00201	0.5333	0.00302
0.88392	0.00366	0.00216	0.50919	0.00386
0.8946	0.00568	0.00336	0.48854	0.00206
0.90365	0.00387	0.00229	0.46361	0.00476
0.9124	0.00507	0.003	0.43796	0.00547
0.92221	0.00366	0.00216	0.40784	0.00512
0.93415	0.00366	0.00216	0.37012	0.005
0.94475	0.00366	0.00216	0.33408	0.00483
0.95386	0.00255	0.00151	0.28661	0.00616
0.96296	0.00245	0.00145	0.23491	0.00506];

kappa_data=[0.84894	4.23E-04	0.16996	6.24E-04
0.85551	0.00245	0.16796	0.00214
0.86674	0.00319	0.16694	0.00159
0.87386	0.00199	0.16814	0.00156
0.88404	0.00258	0.16496	9.28E-04
0.89639	0.00255	0.16172	0.00136
0.90556	0.00333	0.16006	0.00151
0.91587	0.00283	0.16225	9.89E-04
0.92567	6.51E-04	0.15883	0.00103
0.93425	0.00303	0.15792	0.00164
0.94597	0.0029	0.15423	0.00161
0.95584	0.00308	0.16088	0.0038
0.96403	0.00187	0.16716	0.00247];
%save MFS matrix
MFS_aw(1:size(mixture_data,1),6)=mixture_data(:,1);
MFS_OT(1:size(mixture_data,1),6)=mixture_data(:,4);
MFS_VBSBAT(1:size(mixture_data,1),6)=interp1(aw_series,VBSBAT_mix6_C_OA_PM./(VBSBAT_mix6_Caq_PM+VBSBAT_mix6_C_OA_PM),mixture_data(:,1),'linear');
%
plot_name=['water MFS compare ', extension_name];
paper_postion=[0, 0, 3.5, 3.5*1.25].*paper_mult;
figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', paper_postion.*1+.5,'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]


axes1 = axes('Parent',figure1,...
    'Position',[0.16 0.16 0.75 0.75],...
    'LineWidth',1.75,...
    'FontSize',fontsize );  %'CLim',[80 300],
set(axes1,'FontSize',fontsize,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% view(axes1,[175 35]);
grid(axes1,'on');
hold(axes1,'all');
hold on

plot(aw_series,VBSBAT_mix6_C_OA_PM./(VBSBAT_mix6_Caq_PM+VBSBAT_mix6_C_OA_PM),'LineWidth',2)

fill(axes1,[mixture_data(:,1)+mixture_data(:,2);flip(mixture_data(:,1)-mixture_data(:,3))],[mixture_data(:,4)+mixture_data(:,5);flip(mixture_data(:,4)-mixture_data(:,5))],'r','FaceAlpha',0.25,...
    'FaceColor',[0.63921568627451 0.07843137254902 0.180392156862745],...
    'EdgeColor','none');
    
plot(mixture_data(:,1),mixture_data(:,4),'or')

plot(aw_series, VBSBAT_mix6_kappaHGF,'k','LineWidth',2)
plot(kappa_data(:,1),kappa_data(:,3),'ok')
ylabel('MFS or \kappa')
xlabel('water activity (a_w)')
title(extension_name)
annotation(figure1,'textbox',...
    [0.214541666666667 0.385897435897436 0.347958333333333 0.235042735042737],...
    'String',{char(mixture6_compounds)},...
    'FitBoxToText','off','LineStyle','none');
xlim([0.5,1]) 
ylim([0,1])
saveas(figure1, [plot_name '.png'], 'png')
saveas(figure1, [plot_name '.fig'], 'fig')


%% mixture7 calc
extension_name='mixture7';

aw_series=[0.5: 0.001: 1]';

mixture7_compounds={'3,3-Dimethyl Glutaric acid', '2,2-Dimethyl Glutaric acid', 'pimelic acid'}';
mixture7_functional_group={'carboxyl', 'carboxyl', 'carboxyl'}'; % note carboxyl is identical to hydroxyl in the BAT model
mixture7_total_moles= [1, 1, 1]';
mixture7_Psat= mixture7_total_moles.*10^-9;

% {'hydroxyl';'carboxyl';'ether';'ketone';'ester';'hydroperoxide';'hydroperoxideSOA';'PEG'};

mixture7_H2C=[12/7, 12/7, 12/7 ]' ;
mixture7_O2C=[4/7, 4/7, 4/7]' ;
mixture7_N2C=[0/7, 0/7, 0/7]' ;
mixture7_MolarMass=[160.17, 160.17, 160.17]' ;

[Csat_ugPm3, Csat_ugPm3_log10] = Csat_calculate_v1(mixture7_MolarMass,mixture7_Psat, 298, 'no');

% mixture normalized mole fractions
Org_Conc_OM=mixture7_total_moles.*mixture7_MolarMass./(sum(mixture7_total_moles.*mixture7_MolarMass));

[VBSBAT_mix7_C_OA_PM, VBSBAT_mix7_Caq_PM, VBSBAT_mix7_kappaHGF, details_mix3]=VBS_BAT_simulation_v2(...
    Csat_ugPm3, Org_Conc_OM, ...
    mixture7_O2C, mixture7_H2C,  mixture7_MolarMass, aw_series, mixture7_functional_group, ...
    BAT_refinement_mode, VBSBAT_options,  extension_name, mixture7_N2C);

% disp(['fit per aw ' num2str(details_mix1.fit.mean_time_per_fit) ' sec'])

mixture_data=[0.82077	0.0014	7.01E-04	0.74813	0.00244
0.85862	0.00619	0.00508	0.72532	0.01535
0.89681	0.01096	0.00687	0.6714	0.0214
0.91536	0.00761	0.00604	0.61181	0.0252
0.9541	0.0052	0.00417	0.47687	0.02239
0.98874	0.00141	0.00107	0.09738	0.0124
0.99752	0.00493	0.00247	0.03811	0.00286];

kappa_data=[0.81932	7.72E-04	0.09162	0.00637
0.82214	0.00205	0.08518	0.00839
0.8345	0.00404	0.07825	0.00988
0.8455	0.00311	0.076	0.00739
0.85612	0.00344	0.07475	0.00717
0.86467	0.00208	0.07233	0.00751
0.89525	0.0034	0.07075	0.00486
0.9066	0.00351	0.0708	0.00449
0.9145	0.00432	0.0695	0.00404
0.9462	0.00239	0.0654	0.0023
0.95489	0.00293	0.06344	0.00246
0.96575	0.00381	0.06187	0.00264
0.97607	0.00266	0.0637	0.00579
0.98659	0.00285	0.1018	0.03286
0.99469	0.00217	0.11322	0.03319];

%save MFS matrix
MFS_aw(1:size(mixture_data,1),7)=mixture_data(:,1);
MFS_OT(1:size(mixture_data,1),7)=mixture_data(:,4);
MFS_VBSBAT(1:size(mixture_data,1),7)=interp1(aw_series,VBSBAT_mix7_C_OA_PM./(VBSBAT_mix7_Caq_PM+VBSBAT_mix7_C_OA_PM),mixture_data(:,1),'linear');
%

plot_name=['water MFS compare ', extension_name];
paper_postion=[0, 0, 3.5, 3.5*1.25].*paper_mult;
figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', paper_postion.*1+.5,'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]


axes1 = axes('Parent',figure1,...
    'Position',[0.16 0.16 0.75 0.75],...
    'LineWidth',1.75,...
    'FontSize',fontsize );  %'CLim',[80 300],
set(axes1,'FontSize',fontsize,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% view(axes1,[175 35]);
grid(axes1,'on');
hold(axes1,'all');
hold on

plot(aw_series,VBSBAT_mix7_C_OA_PM./(VBSBAT_mix7_Caq_PM+VBSBAT_mix7_C_OA_PM),'LineWidth',2)

fill(axes1,[mixture_data(:,1)+mixture_data(:,2);flip(mixture_data(:,1)-mixture_data(:,3))],[mixture_data(:,4)+mixture_data(:,5);flip(mixture_data(:,4)-mixture_data(:,5))],'r','FaceAlpha',0.25,...
    'FaceColor',[0.63921568627451 0.07843137254902 0.180392156862745],...
    'EdgeColor','none');
    
plot(mixture_data(:,1),mixture_data(:,4),'or')

plot(aw_series, VBSBAT_mix7_kappaHGF,'k','LineWidth',2)
plot(kappa_data(:,1),kappa_data(:,3),'ok')
ylabel('MFS or \kappa')
xlabel('water activity (a_w)')
title(extension_name)
annotation(figure1,'textbox',...
    [0.214541666666667 0.385897435897436 0.347958333333333 0.235042735042737],...
    'String',{char(mixture7_compounds)},...
    'FitBoxToText','off','LineStyle','none');

xlim([0.5,1]) 
ylim([0,1])
saveas(figure1, [plot_name '.png'], 'png')
saveas(figure1, [plot_name '.fig'], 'fig')



%%

plot_name=['MFS mode vs measurement '];
paper_postion=[0, 0, 3.5, 3.5*1.25].*paper_mult;
figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', paper_postion.*1+.5,'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]


axes1 = axes('Parent',figure1,...
    'Position',[0.16 0.16 0.75 0.75],...
    'LineWidth',1.75,...
    'FontSize',fontsize );  %'CLim',[80 300],
set(axes1,'FontSize',fontsize,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% view(axes1,[175 35]);
grid(axes1,'on');
hold(axes1,'all');
hold on


MFS_aw=replace_data_A_to_B_KGv1(MFS_aw,0,NaN);
MFS_OT=replace_data_A_to_B_KGv1(MFS_OT,0,NaN);
MFS_VBSBAT=replace_data_A_to_B_KGv1(MFS_VBSBAT,0,NaN);

plot([0,1],[0,1],'k','LineWidth',1.5)
% plot([0,1],[0,1.1],':k','LineWidth',1.5)
% plot([0,1],[0,.9],':k','LineWidth',1.5)
% plot([0,1],[0.1,1],':b','LineWidth',1.5)
% plot([0,1],[-0.1,1],':b','LineWidth',1.5)

% error bars
fill([0,1,1,0],[0.1,1,1,-0.1],'r','FaceAlpha',0.25,...
        'FaceColor',[0.00,0.45,0.74],...
        'EdgeColor','none')
fill([0,1,1,0],[0,1.1,.9,0],'r','FaceAlpha',0.25,...
        'FaceColor',[0.65,0.65,0.65],...
        'EdgeColor','none')

    
xlabel('MFS measured')
ylabel('MFS BAT')

m_size=4;
% Create plot
plot(MFS_OT(:,1),MFS_VBSBAT(:,1),...
    'MarkerFaceColor',[0.698039215686274 0.694117647058824 0.850980392156863],...
    'MarkerEdgeColor',[0.411764705882353 0.325490196078431 0.611764705882353],...
    'MarkerSize',m_size,...
    'Marker','diamond',...
    'LineStyle','none', 'DisplayName', 'Mixture 1');

% Create plot
plot(MFS_OT(:,2),MFS_VBSBAT(:,2),...
    'MarkerFaceColor',[0.925490196078431 0.92156862745098 0.964705882352941],...
    'MarkerEdgeColor',[0.411764705882353 0.325490196078431 0.611764705882353],...
    'MarkerSize',m_size,...
    'Marker','diamond',...
    'LineStyle','none', 'DisplayName', 'Mixture 2');

% Create plot
plot(MFS_OT(:,3),MFS_VBSBAT(:,3),...
    'MarkerFaceColor',[0.627450980392157 0.803921568627451 0.309803921568627],...
    'MarkerEdgeColor',[0.250980392156863 0.388235294117647 0.145098039215686],...
    'MarkerSize',m_size,...
    'Marker','o',...
    'LineStyle','none', 'DisplayName', 'Mixture 3');

% Create plot
plot(MFS_OT(:,4),MFS_VBSBAT(:,4),...
    'MarkerFaceColor',[0.4 0.176470588235294 0.568627450980392],...
    'MarkerEdgeColor',[0.533333333333333 0.490196078431373 0.729411764705882],...
    'MarkerSize',m_size,...
    'Marker','v',...
    'LineStyle','none', 'DisplayName', 'Mixture 4');

% Create plot
plot(MFS_OT(:,5),MFS_VBSBAT(:,5),...
    'MarkerFaceColor',[0.250980392156863 0.388235294117647 0.145098039215686],...
    'MarkerEdgeColor',[0.27843137254902 0.36078431372549 0.196078431372549],...
    'MarkerSize',m_size,...
    'Marker','o',...
    'LineStyle','none', 'DisplayName', 'Mixture 5');

% Create plot
plot(MFS_OT(:,6),MFS_VBSBAT(:,6),...
    'MarkerFaceColor',[0.533333333333333 0.490196078431373 0.729411764705882],...
    'MarkerEdgeColor',[0.411764705882353 0.325490196078431 0.611764705882353],...
    'MarkerSize',m_size,...
    'Marker','diamond',...
    'LineStyle','none', 'DisplayName', 'Mixture 6');

% Create plot
plot(MFS_OT(:,7),MFS_VBSBAT(:,7),...
    'MarkerFaceColor',[0.909803921568627 0.690196078431373 0.129411764705882],...
    'MarkerEdgeColor',[0.603921568627451 0.450980392156863 0.168627450980392],...
    'MarkerSize',m_size,...
    'Marker','square',...
    'LineStyle','none', 'DisplayName', 'Mixture 7');


xlim([0,1]) 
ylim([0,1])
saveas(figure1, [plot_name '.png'], 'png')
saveas(figure1, [plot_name '.fig'], 'fig')