%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2018-Jul-23 11:11 AM
% Copyright 2018 Kyle Gorkowski
%%
clear

%% settings
mode='hydroxyl'; %    'hydroxyl'='carboxyl', 'hydroperoxideSOA', 'PEG', 'hydroperoxide'
molar_mass=750;

O2C_values=[0:0.005:2]';

a_w_specific=[0.999, 0.70 ,0.40]';

fix_molarmass1=18.016/molar_mass;
Dp_sim_nm=200;
sigma_droplet_N_m=0.072;


%% program
mole_frac_scan=[1:-0.0005:0]';

fit_tolerance=10^-6;

ycal_org_specific_beta=zeros(length(a_w_specific),length(O2C_values));
ycal_org_specific_alpha=ycal_org_specific_beta;
ycal_water_specific_beta=ycal_org_specific_beta;
ycal_water_specific_alpha=ycal_org_specific_beta;
activity_org_specific=ycal_org_specific_beta;
mole_fraction_org_beta=ycal_org_specific_beta;
mole_fraction_org_alpha=ycal_org_specific_beta;

H2C=0;
progressbartext('activity space calc')
for i=1:length(O2C_values)
    
    [func1, func2, ycal_water(:,i), ycalc_org(:,i), activity_water(:,i), activity_org(:,i), mass_fraction1(:,i), mass_fraction2(:,i),Gibbs_RT, dGibbs_RTdx2]...
        =BAT_properties_calculation_v1(mole_frac_scan, O2C_values(i,1),  O2C_values(i,1).*H2C, fix_molarmass1,mode,[]);
%     
%     density_org_g_cm3=Org_density_Estimate_KGv1(18.016/fix_molarmass1, O2C_values(i), O2C_values(i).*0);
%     
%     [kappa_SatCritical(i,1), kappa_a_w(i,1), SatCritical(i,1), Dp_critical(i,1), V_w_critical(i,1)]= kappa_critical_calc_for_McGlashan_v1(Dp_sim_nm, sigma_droplet_N_m,...
%     density_org_g_cm3, mass_fraction1(:,i),  mass_fraction2(:,i), activity_water(:,i));
%     
%     
    [ycal_waterLw(:,i), ycalc_orgL(:,i), activity_waterLw(:,i), activity_orgL(:,i), ycal_waterU(:,i), ...
        ycalc_orgU(:,i), activity_waterU(:,i), activity_orgU(:,i), mole_frac_fit_lowRH_water(:,i),...
        mole_frac_fit_upperRH_water(:,i) ]=mcglashan_activity_calc_SRH_edges_KGv8(0.5,...
        O2C_values(i,1),  O2C_values(i,1).*H2C, fix_molarmass1,mode, 'a_water SRH', fit_tolerance);
    
    [ycal_waterL(:,i), ycalc_orgL(:,i), activity_waterL(:,i), activity_orgL(:,i), ycal_waterU(:,i), ...
        ycalc_orgU(:,i), activity_waterU(:,i), activity_orgU(:,i), mole_frac_fit_lowRH_org(:,i),...
        mole_frac_fit_upperRH_org(:,i) ]=mcglashan_activity_calc_SRH_edges_KGv8(0.5,...
        O2C_values(i,1),  O2C_values(i,1).*H2C, fix_molarmass1,mode, 'a_org SRH', fit_tolerance);
    
    
    
    
    % interperate to specific a_w
    [~,p_i]=max(activity_water(:,i));
    
    [phaseSep_via_activity,phaseSep_via_activity_curvature,index_phase_sep_starts,index_phase_sep_end]=finds_PhaseSep_and_activity_curve_dips_v2(activity_water(:,i));
    index_phase_sep_starts=p_i;
    % index_phase_sep_end=p_i;
    
    if phaseSep_via_activity_curvature==1
        % beta organic rich side
        ycal_org_specific_beta(:,i)=interp1(activity_water(1:index_phase_sep_starts,i), ycalc_org(1:index_phase_sep_starts,i),a_w_specific,'linear');
        ycal_water_specific_beta(:,i)=interp1(activity_water(1:index_phase_sep_starts,i), ycal_water(1:index_phase_sep_starts,i),a_w_specific,'linear');
        mole_fraction_org_beta(:,i)=interp1(activity_water(1:index_phase_sep_starts,i), mole_frac_scan(1:index_phase_sep_starts,1),a_w_specific,'linear');
        
        
        % alpha water rich side
        ycal_org_specific_alpha(:,i)=interp1(activity_water(index_phase_sep_end-10:end,i), ycalc_org(index_phase_sep_end-10:end,i),a_w_specific,'linear');
        ycal_water_specific_alpha(:,i)=interp1(activity_water(index_phase_sep_end-10:end,i), ycal_water(index_phase_sep_end-10:end,i),a_w_specific,'linear');
        mole_fraction_org_alpha(:,i)=interp1(activity_water(index_phase_sep_end-10:end,i), mole_frac_scan(index_phase_sep_end-10:end,1),a_w_specific,'linear');
        
        %ycalc_org_specific(:,i)=interp1(activity_water(1:p_i,i), ycalc_org(1:p_i,i),a_w_specific,'linear');
        %ycal_water_specific(:,i)=interp1(activity_water(1:p_i,i), ycal_water(1:p_i,i),a_w_specific,'linear');
        %activity_org_specific(:,i)=interp1(activity_water(1:p_i,i), activity_org(1:p_i,i),a_w_specific,'linear');
        %mole_fraction_org(:,i)=interp1(activity_water(1:p_i,i), mole_frac_scan(1:p_i,1),a_w_specific,'linear');
        
    else
        index_phase_sep_starts=length(mole_frac_scan);
        % beta organic rich side
        ycal_org_specific_beta(:,i)=interp1(activity_water(1:index_phase_sep_starts,i), ycalc_org(1:index_phase_sep_starts,i),a_w_specific,'linear');
        ycal_water_specific_beta(:,i)=interp1(activity_water(1:index_phase_sep_starts,i), ycal_water(1:index_phase_sep_starts,i),a_w_specific,'linear');
        mole_fraction_org_beta(:,i)=interp1(activity_water(1:index_phase_sep_starts,i), mole_frac_scan(1:index_phase_sep_starts,1),a_w_specific,'linear');

    end
    
    
    
    progressbartext( i/length(O2C_values))
end

% d_ideal=(activity_water-repmat(1-mole_frac_scan,1,length(O2C_values)));
% [d_ideal_line,d_ideal_i]=min(d_ideal,[],2);
% 
% ideal_mole_frac=mole_frac_scan(d_ideal_i);
% figure
% hold on
% plot(mole_frac_scan,activity_org')
% xlabel('mole fraction')
% ylabel('a_organic')

% figure
% hold on
% surf(1-mole_frac_scan,O2C_values,log10(ycalc_org'),'LineStyle','none')
% view([0 90]);
% xlabel('x_{water}')
% ylabel('O:C')
% zlabel('\gamma organic')
% colorbar('Direction','reverse');
mole_frac_fit_upperRH_water=replace_data_A_to_B_KGv1(mole_frac_fit_upperRH_water,0,NaN);
mole_frac_fit_upperRH_org=replace_data_A_to_B_KGv1(mole_frac_fit_upperRH_org,0,NaN);
mole_frac_fit_lowRH_water=replace_data_A_to_B_KGv1(mole_frac_fit_lowRH_water,0,NaN);
mole_frac_fit_lowRH_org=replace_data_A_to_B_KGv1(mole_frac_fit_lowRH_org,0,NaN);

combined_SRH_water=[mole_frac_fit_lowRH_water,mole_frac_fit_upperRH_water]';
combined_SRH_org=[mole_frac_fit_lowRH_org,mole_frac_fit_upperRH_org]';
[combined_SRH_water, O2C_SRH_water]=Removes_nan_rows_dual(combined_SRH_water,[O2C_values;O2C_values]);
[combined_SRH_org, O2C_SRH_org]=Removes_nan_rows_dual(combined_SRH_org,[O2C_values;O2C_values]);

[combined_SRH_water, O2C_SRH_water]=sort_and_remove_duplicates_KGv1(combined_SRH_water, O2C_SRH_water);
[combined_SRH_org, O2C_SRH_org]=sort_and_remove_duplicates_KGv1(combined_SRH_org, O2C_SRH_org);


%% mole fraction
plot_name=['Mole fraction- water activity vs O2C SI MW' num2str(molar_mass)];
 paper_postion=[0, 0, 3.5, 3.5].*1;

figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', paper_postion+.5,'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]
axes1 = axes('Parent',figure1,...
    'Position',[0.1846875 0.1484375 0.78015842495637 0.815],...
    'LineWidth',1.75,...
    'FontSize',12);  %'CLim',[80 300],
set(axes1,'FontSize',12,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% reverse color map
%fliped_color=[0.976899981498718,0.983900010585785,0.0804999992251396;0.969114303588867,0.960628569126129,0.106361903250217;0.966412723064423,0.936758756637573,0.122133336961269;0.963711082935333,0.912888884544373,0.137904763221741;0.961009502410889,0.889019072055817,0.153676196932793;0.966782510280609,0.864744484424591,0.167588099837303;0.972555518150330,0.840469837188721,0.181500002741814;0.978328585624695,0.816195249557495,0.195411905646324;0.984101593494415,0.791920661926270,0.209323808550835;0.989874601364136,0.767646014690399,0.223235711455345;0.995647609233856,0.743371427059174,0.237147614359856;0.972757160663605,0.729771435260773,0.239442855119705;0.939257144927979,0.728785693645477,0.209957137703896;0.903542876243591,0.733028590679169,0.177414283156395;0.863433361053467,0.740599989891052,0.159633338451385;0.820314288139343,0.749814271926880,0.153528571128845;0.773833334445953,0.759804785251617,0.164609521627426;0.724200010299683,0.769842863082886,0.191028565168381;0.671985685825348,0.779271423816681,0.222699999809265;0.608423769474030,0.782274603843689,0.272054761648178;0.544861912727356,0.785277783870697,0.321409523487091;0.481299996376038,0.788280963897705,0.370764285326004;0.417738080024719,0.791284143924713,0.420119047164917;0.354176163673401,0.794287323951721,0.469473809003830;0.290614277124405,0.797290503978729,0.518828570842743;0.246957138180733,0.791795253753662,0.556742846965790;0.216085717082024,0.784300029277802,0.592299997806549;0.201290473341942,0.778520643711090,0.613347589969635;0.186495244503021,0.772741317749023,0.634395241737366;0.171700000762939,0.766961932182312,0.655442833900452;0.140771433711052,0.758400022983551,0.684157133102417;0.0963952392339706,0.750000000000000,0.712038099765778;0.0433285720646381,0.741095244884491,0.739409506320953;0.00665714265778661,0.731214284896851,0.766014277935028;0.00357142859138548,0.720266640186310,0.791700005531311;0.0296666659414768,0.708166658878326,0.816333353519440;0.0688714310526848,0.694771409034729,0.839357137680054;0.0952095240354538,0.679828584194183,0.859780967235565;0.109480954706669,0.662428557872772,0.872284114360809;0.123752377927303,0.645028591156006,0.884787321090698;0.138023808598518,0.627628564834595,0.897290468215942;0.146011903882027,0.608914256095886,0.909545242786407;0.153999999165535,0.590200006961823,0.921800017356873;0.168742850422859,0.570261895656586,0.935871422290802;0.176438093185425,0.549904763698578,0.952019035816193;0.179921418428421,0.528638124465942,0.965907096862793;0.183404758572578,0.507371425628662,0.979795217514038;0.196333333849907,0.484719038009644,0.989152371883392;0.220642849802971,0.460257142782211,0.997285723686218;0.244033336639404,0.435833334922791,0.998833358287811;0.260242849588394,0.412328571081161,0.995157122612000;0.262822449207306,0.389977544546127,0.981430590152740;0.265402019023895,0.367626518011093,0.967704057693481;0.267981618642807,0.345275491476059,0.953977525234222;0.270561218261719,0.322924494743347,0.940251052379608;0.273140817880630,0.300573468208313,0.926524519920349;0.275720387697220,0.278222441673279,0.912797987461090;0.278299987316132,0.255871415138245,0.899071455001831;0.275114297866821,0.234238088130951,0.870985686779022;0.270647615194321,0.214676186442375,0.836371421813965;0.263535708189011,0.198607146739960,0.792353570461273;0.256423801183701,0.182538092136383,0.748335719108582;0.249311909079552,0.166469037532806,0.704317867755890;0.242200002074242,0.150399997830391,0.660300016403198];
fliped_color=flip(artistic_color);
set(gca,'Colormap',fliped_color)

hold on

Cmatrix=contourf(1-mole_frac_scan,O2C_values,thresholding_data_to_nan_KGv1(activity_water,0,1)','LineWidth',1,'LineStyle','none',...
    'LineColor',[0.9501960784313726 0.9501960784313726 0.9501960784313726],...
    'levelStep', 0.025);
contour(1-mole_frac_scan,O2C_values,thresholding_data_to_nan_KGv1(activity_water,0,1.01)','LineWidth',1,...
    'LineColor',[0.9501960784313726 0.9501960784313726 0.9501960784313726],...
    'levels', 0.1)
contour(1-mole_frac_scan,O2C_values,thresholding_data_to_nan_KGv1(activity_water,0,1)','LineWidth',1,...
    'LineColor',[0.9501960784313726 0.9501960784313726 0.9501960784313726],...
    'levels', 0.99)
% plot(1-combined_SRH_water, O2C_SRH_water,'LineWidth',2,'Color',[0 0 0])
% plot(1-combined_SRH_org, O2C_SRH_org,'--','LineWidth',2,'Color',[0 0 0])
area(1-combined_SRH_water, O2C_SRH_water,'LineWidth',2,'FaceColor',[0.24705882370472 0.24705882370472 0.24705882370472],  'EdgeColor',[0.24705882370472 0.24705882370472 0.24705882370472])
area(1-combined_SRH_org, O2C_SRH_org,'LineWidth',2,'FaceAlpha',.7,'FaceColor',[0.800000011920929 0.800000011920929 0.800000011920929],...
    'EdgeColor',[0.800000011920929 0.800000011920929 0.800000011920929])
%plot(1-ideal_mole_frac,O2C_values,'LineWidth',10)
set(axes1 ,'Layer', 'Top')

view([0 90]);
xlabel('Water mole fraction')
ylabel('O:C (hydroxyl)')
zlabel('a_w')
% colorbar(axes1,'eastoutside','Position',...
%     [0.930440820918603 0.125 0.0193236714975846 0.815]);
% linkaxes([axes1,axes2],'xy')
%colorbar('Direction','reverse');
saveas(figure1, [plot_name '.png'], 'png')
saveas(figure1, [plot_name '.fig'], 'fig')

%% mass fraction
plot_name=['Mass fraction- water activity vs O2C SI MW' num2str(molar_mass)];
 paper_postion=[0, 0, 3.5, 3.5].*1;

figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', paper_postion+.5,'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]
axes1 = axes('Parent',figure1,...
    'Position',[0.1846875 0.1484375 0.78015842495637 0.815],...
    'LineWidth',1.75,...
    'FontSize',12);  %'CLim',[80 300],
set(axes1,'FontSize',12,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% reverse color map
%fliped_color=[0.976899981498718,0.983900010585785,0.0804999992251396;0.969114303588867,0.960628569126129,0.106361903250217;0.966412723064423,0.936758756637573,0.122133336961269;0.963711082935333,0.912888884544373,0.137904763221741;0.961009502410889,0.889019072055817,0.153676196932793;0.966782510280609,0.864744484424591,0.167588099837303;0.972555518150330,0.840469837188721,0.181500002741814;0.978328585624695,0.816195249557495,0.195411905646324;0.984101593494415,0.791920661926270,0.209323808550835;0.989874601364136,0.767646014690399,0.223235711455345;0.995647609233856,0.743371427059174,0.237147614359856;0.972757160663605,0.729771435260773,0.239442855119705;0.939257144927979,0.728785693645477,0.209957137703896;0.903542876243591,0.733028590679169,0.177414283156395;0.863433361053467,0.740599989891052,0.159633338451385;0.820314288139343,0.749814271926880,0.153528571128845;0.773833334445953,0.759804785251617,0.164609521627426;0.724200010299683,0.769842863082886,0.191028565168381;0.671985685825348,0.779271423816681,0.222699999809265;0.608423769474030,0.782274603843689,0.272054761648178;0.544861912727356,0.785277783870697,0.321409523487091;0.481299996376038,0.788280963897705,0.370764285326004;0.417738080024719,0.791284143924713,0.420119047164917;0.354176163673401,0.794287323951721,0.469473809003830;0.290614277124405,0.797290503978729,0.518828570842743;0.246957138180733,0.791795253753662,0.556742846965790;0.216085717082024,0.784300029277802,0.592299997806549;0.201290473341942,0.778520643711090,0.613347589969635;0.186495244503021,0.772741317749023,0.634395241737366;0.171700000762939,0.766961932182312,0.655442833900452;0.140771433711052,0.758400022983551,0.684157133102417;0.0963952392339706,0.750000000000000,0.712038099765778;0.0433285720646381,0.741095244884491,0.739409506320953;0.00665714265778661,0.731214284896851,0.766014277935028;0.00357142859138548,0.720266640186310,0.791700005531311;0.0296666659414768,0.708166658878326,0.816333353519440;0.0688714310526848,0.694771409034729,0.839357137680054;0.0952095240354538,0.679828584194183,0.859780967235565;0.109480954706669,0.662428557872772,0.872284114360809;0.123752377927303,0.645028591156006,0.884787321090698;0.138023808598518,0.627628564834595,0.897290468215942;0.146011903882027,0.608914256095886,0.909545242786407;0.153999999165535,0.590200006961823,0.921800017356873;0.168742850422859,0.570261895656586,0.935871422290802;0.176438093185425,0.549904763698578,0.952019035816193;0.179921418428421,0.528638124465942,0.965907096862793;0.183404758572578,0.507371425628662,0.979795217514038;0.196333333849907,0.484719038009644,0.989152371883392;0.220642849802971,0.460257142782211,0.997285723686218;0.244033336639404,0.435833334922791,0.998833358287811;0.260242849588394,0.412328571081161,0.995157122612000;0.262822449207306,0.389977544546127,0.981430590152740;0.265402019023895,0.367626518011093,0.967704057693481;0.267981618642807,0.345275491476059,0.953977525234222;0.270561218261719,0.322924494743347,0.940251052379608;0.273140817880630,0.300573468208313,0.926524519920349;0.275720387697220,0.278222441673279,0.912797987461090;0.278299987316132,0.255871415138245,0.899071455001831;0.275114297866821,0.234238088130951,0.870985686779022;0.270647615194321,0.214676186442375,0.836371421813965;0.263535708189011,0.198607146739960,0.792353570461273;0.256423801183701,0.182538092136383,0.748335719108582;0.249311909079552,0.166469037532806,0.704317867755890;0.242200002074242,0.150399997830391,0.660300016403198];
fliped_color=flip(artistic_color);
set(gca,'Colormap',fliped_color)

hold on
[mass_frac_water, mass_frac_org] = mole_to_mass_fraction(1-mole_frac_scan, mole_frac_scan, 18.016, molar_mass);

Cmatrix=contourf(mass_frac_water,O2C_values,thresholding_data_to_nan_KGv1(activity_water,0,1)','LineWidth',1,'LineStyle','none',...
    'LineColor',[0.9501960784313726 0.9501960784313726 0.9501960784313726],...
    'levelStep', 0.025);
contour(mass_frac_water,O2C_values,thresholding_data_to_nan_KGv1(activity_water,0,1.01)','LineWidth',1,...
    'LineColor',[0.9501960784313726 0.9501960784313726 0.9501960784313726],...
    'levels', 0.1)
contour(mass_frac_water,O2C_values,thresholding_data_to_nan_KGv1(activity_water,0,1)','LineWidth',1,...
    'LineColor',[0.9501960784313726 0.9501960784313726 0.9501960784313726],...
    'levels', 0.99)
% plot(1-combined_SRH_water, O2C_SRH_water,'LineWidth',2,'Color',[0 0 0])
% plot(1-combined_SRH_org, O2C_SRH_org,'--','LineWidth',2,'Color',[0 0 0])
[combined_SRH_water_mass, ~] = mole_to_mass_fraction(1-combined_SRH_water, combined_SRH_water, 18.016, molar_mass);
area(combined_SRH_water_mass, O2C_SRH_water,'LineWidth',2,'FaceColor',[0.24705882370472 0.24705882370472 0.24705882370472],  'EdgeColor',[0.24705882370472 0.24705882370472 0.24705882370472])
[combined_SRH_org_mass, ~] = mole_to_mass_fraction(1-combined_SRH_org, combined_SRH_org, 18.016, molar_mass);
area(combined_SRH_org_mass, O2C_SRH_org,'LineWidth',2,'FaceAlpha',.7,'FaceColor',[0.800000011920929 0.800000011920929 0.800000011920929],...
    'EdgeColor',[0.800000011920929 0.800000011920929 0.800000011920929])
%plot(1-ideal_mole_frac,O2C_values,'LineWidth',10)
set(axes1 ,'Layer', 'Top')

view([0 90]);
xlabel('Water mass fraction')
ylabel('O:C (hydroxyl)')
zlabel('a_w')
% colorbar(axes1,'eastoutside','Position',...
%     [0.930440820918603 0.125 0.0193236714975846 0.815]);
% linkaxes([axes1,axes2],'xy')
%colorbar('Direction','reverse');
saveas(figure1, [plot_name '.png'], 'png')
saveas(figure1, [plot_name '.fig'], 'fig')


%
% plot_name=['organic activity vs O2C v8'];
% paper_postion=[.5, .5, 11, 6.5];
%
% figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', [1, 1, 11.5, 7],'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]
% axes1 = axes('Parent',figure1,...
%     'Position',[0.13 0.135 0.81584249563699825 0.815],...
%     'LineWidth',1.75,...
%     'FontSize',16);  %'CLim',[80 300],
% set(axes1,'FontSize',13,'XMinorTick','on','YMinorTick','on',...
%     'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% hold on
% %surf(1-mole_frac_scan,O2C_values,thresholding_data_to_nan_KGv1(activity_org,0,1)','LineStyle','none')
% contourf(1-mole_frac_scan,O2C_values,thresholding_data_to_nan_KGv1(activity_org,0,1)','LineWidth',2,...
%     'LineColor',[0.901960784313726 0.901960784313726 0.901960784313726],...
%     'LevelStep',0.1)
%  area(1-combined_SRH_water, O2C_SRH_water,'LineWidth',2,'FaceColor',[0.24705882370472 0.24705882370472 0.24705882370472])
%  area(1-combined_SRH_org, O2C_SRH_org,'LineWidth',2,'FaceAlpha',0.5,'FaceColor',[1 1 1],...
%     'EdgeColor',[0 0.498039215803146 0])
%
% view([0 90]);
% xlabel('Water mole fraction')
% ylabel('O:C')
% zlabel('a_w')
% % colorbar('Direction','reverse');
% saveas(figure1, [plot_name '.png'], 'png')
% saveas(figure1, [plot_name '.fig'], 'fig')

% ycal_org_specific_alpha=replace_data_A_to_B_KGv1(ycal_org_specific_alpha,0,NaN);
% ycal_water_specific_alpha=replace_data_A_to_B_KGv1(ycal_water_specific_alpha,0,NaN);
% mole_fraction_org_alpha=replace_data_A_to_B_KGv1(mole_fraction_org_alpha,0,NaN);

% 
% 
% plot_name=['activity coefficients vs O2C'];
% paper_postion=[.5, .5, 11, 6.5];
% 
% figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', [1, 1, 11.5, 7],'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]
% axes1 = axes('Parent',figure1,...
%     'Position',[0.13 0.135 0.81584249563699825 0.815],...
%     'LineWidth',1.75,...
%     'FontSize',16);  %'CLim',[80 300],
% set(axes1,'FontSize',13,'XMinorTick','on','YMinorTick','on',...
%     'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% hold on
% grid on
% colorOrder = get(gca, 'ColorOrder');
% lines_totoal=length(a_w_specific);
% colorOrder=repmat(colorOrder, ceil(lines_totoal./length(colorOrder)),1);
% 
% for p_i=1:length(a_w_specific)
%     
%     if not(isnan(ycal_org_specific_alpha(p_i,1)))
%         plot(O2C_values,ycal_org_specific_alpha(p_i,:)', 'DisplayName', ['\gamma_{org}^{\alpha}  a_w=' num2str(a_w_specific(p_i))],...
%             'Color',colorOrder(p_i,:), 'LineWidth', 4)
%     end
%     plot(O2C_values,ycal_org_specific_beta(p_i,:)','DisplayName', ['\gamma_{org}^{\beta}  a_w=' num2str(a_w_specific(p_i))],...
%         'Color',colorOrder(p_i,:), 'LineWidth', 2)
%     
%     if not(isnan(ycal_water_specific_alpha(p_i,1)))
%         plot(O2C_values,ycal_water_specific_alpha(p_i,:)','--',  'DisplayName', ['\gamma_{w}^{\alpha}  a_w=' num2str(a_w_specific(p_i))],...
%             'Color',colorOrder(p_i,:), 'LineWidth', 4)
%     end
%     plot(O2C_values,ycal_water_specific_beta(p_i,:)','--', 'DisplayName', ['\gamma_{w}^{\beta}  a_w=' num2str(a_w_specific(p_i))],...
%         'Color',colorOrder(p_i,:), 'LineWidth', 2)
%     
% end
% ylim([0,3])
% xlabel('O:C_{hydroxyl}')
% ylabel('\gamma')
% legend(axes1,'show');
% saveas(figure1, [plot_name '.png'], 'png')
% saveas(figure1, [plot_name '.fig'], 'fig')
% 
% %
% plot_name=['mole fraction vs O2C'];
% paper_postion=[.5, .5, 11, 6.5];
%
% figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', [1, 1, 11.5, 7],'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]
% axes1 = axes('Parent',figure1,...
%     'Position',[0.13 0.135 0.81584249563699825 0.815],...
%     'LineWidth',1.75,...
%     'FontSize',16);  %'CLim',[80 300],
% set(axes1,'FontSize',13,'XMinorTick','on','YMinorTick','on',...
%     'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% hold on
% grid on
% colorOrder = get(gca, 'ColorOrder');
% lines_totoal=length(a_w_specific);
% colorOrder=repmat(colorOrder, ceil(lines_totoal./length(colorOrder)),1);
%
% for p_i=1:length(a_w_specific)
%
%     if not(isnan(ycal_org_specific_alpha(p_i,1)))
%         plot(O2C_values,(1-mole_fraction_org_alpha(p_i,:))', 'DisplayName', ['x_{org}^{\alpha}  a_w=' num2str(a_w_specific(p_i))],...
%             'Color',colorOrder(p_i,:), 'LineWidth', 4)
%     end
%     plot(O2C_values,(1-mole_fraction_org_beta(p_i,:))','DisplayName', ['x_{org}^{\beta}  a_w=' num2str(a_w_specific(p_i))],...
%         'Color',colorOrder(p_i,:), 'LineWidth', 2)
%
%
% end
% ylim([0,1.1])
% xlabel('O:C_{hydroxyl}')
% ylabel('water mole fraction')
% legend(axes1,'show');
% saveas(figure1, [plot_name '.png'], 'png')
% saveas(figure1, [plot_name '.fig'], 'fig')


