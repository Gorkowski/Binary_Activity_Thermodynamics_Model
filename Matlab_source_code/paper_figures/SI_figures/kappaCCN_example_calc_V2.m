clear
% kappa test calc

plot_name=['figure kappaCCN test O2C'];
% 
% load('binary_fit_data_imported.mat')
% index_id=8;

% test system
O2C_value=.380;
Molar_mass=200;
Dp_sim_nm=100;
density_water_g_cm3=1;
% density_org_g_cm3=Org_density_Estimate_KGv1(Molar_mass,O2C_value, O2C_value.*0);
% sigma_droplet_N_m=0.072;
sigmaW_droplet_N_m=0.072;
sigmaOrg_droplet_N_m=0.03;
surface_tension_method='water'; % water, organic, volume mix

Mw=18.016.*10^-3;
R=8.31446;
Temp=293;
mole_frac_scan=[1:-0.0001:0]';
fix_molarmass1=18.016./Molar_mass;
H2C=0;
    mode='hydroxyl';
    [func1, func2, ycal_water, ycalc_org, activity_water, activity_org, mass_fraction1, mass_fraction2,Gibbs_RT, dGibbs_RTdx2]...
        =BAT_properties_calculation_v1(mole_frac_scan, O2C_value,  O2C_value.*H2C, fix_molarmass1,mode,[]);
    
     density_org_g_cm3=Org_density_Estimate_KGv1(18.016/fix_molarmass1, O2C_value, O2C_value.*0);

mass_fraction_water=mass_fraction1;
mass_fraction_org=mass_fraction2;
a_w=activity_water;


V_dp_m=4./3.*pi.*(Dp_sim_nm.*10^-9./2).^3; % volume organics

w_org=density_org_g_cm3.*1000.*V_dp_m;  % estimate mass of org
V_w=1./(density_water_g_cm3.*1000).*w_org.*(mass_fraction_water./mass_fraction_org);% volume water

V_total=V_dp_m+V_w;

if strcmpi(surface_tension_method, 'volume mix')
    sigma_droplet_N_m=(sigmaW_droplet_N_m.*V_w+sigmaOrg_droplet_N_m.*V_dp_m)./V_total;
elseif strcmpi(surface_tension_method, 'water')
    sigma_droplet_N_m=sigmaW_droplet_N_m;
elseif strcmpi(surface_tension_method, 'organic')
    sigma_droplet_N_m=sigmaOrg_droplet_N_m;
end


Diameter_total=real(2*(V_total.*3/(4*pi)).^(1/3));

SatRatio=a_w.*exp_wlimiter(4.*sigma_droplet_N_m.*Mw./(R.*Temp.*density_water_g_cm3.*1000.*Diameter_total));
SatRatio_percent=(SatRatio-1).*100;

kappa=(1./a_w-1).*(V_w./V_dp_m);

% check for phase seperation
% [phaseSep_via_activity,phaseSep_via_activity_curvature,index_phase_sep_starts,index_phase_sep_end]=finds_PhaseSep_and_activity_curve_dips_v2(a_w);
[phaseSep_via_activity_curvature,~,index_phase_sep_end, index_phase_sep_starts]=finds_PhaseSep_w_and_org(activity_water,activity_org);

disp(activity_water(index_phase_sep_end)-activity_water(index_phase_sep_starts))
disp(activity_water(index_phase_sep_end))
disp(activity_water([index_phase_sep_starts-1,index_phase_sep_starts,index_phase_sep_starts+1]))

fixed_aw_val=0.99;
alternate.beta.phase=0;
alternate.fixed_aw.fixed_aw_val=fixed_aw_val; %save
alternate.fixed_aw.beta.phase=0;

if phaseSep_via_activity_curvature==0 % no phase sep.
[SatCritical,Sc_index]=max(SatRatio);

kappa_SatCritical=kappa(Sc_index,1);
kappa_a_w=a_w(Sc_index,1);

Dp_critical=Diameter_total(Sc_index,1);
V_w_critical=V_w(Sc_index,1);

% kappa at fixed aw alpha phase
[~,aw_index]=min(abs(a_w-fixed_aw_val));
alternate.fixed_aw.alpha.aw=a_w(aw_index,1);

alternate.fixed_aw.alpha.SatCritical=SatRatio(aw_index,1);

alternate.fixed_aw.alpha.kappa_SatCritical=kappa(aw_index,1);
alternate.fixed_aw.alpha.kappa_a_w=a_w(aw_index,1);

alternate.fixed_aw.alpha.Dp_critical=Diameter_total(aw_index,1);
alternate.fixed_aw.alpha.V_w_critical=V_w(aw_index,1);

index_phase_sep_starts=length(a_w);
index_phase_sep_end=1;

else % has phase sep so pulls kappa at onset of phase sep.
    
    %% fixed aw
   % kappa at fixed aw alpha phase
[~,aw_index]=min(abs(a_w(1:index_phase_sep_starts,:)-fixed_aw_val));
alternate.fixed_aw.beta.aw=a_w(aw_index,1);

alternate.fixed_aw.beta.SatCritical=SatRatio(aw_index,1);

alternate.fixed_aw.beta.kappa_SatCritical=kappa(aw_index,1);
alternate.fixed_aw.beta.kappa_a_w=a_w(aw_index,1);

alternate.fixed_aw.beta.Dp_critical=Diameter_total(aw_index,1);
alternate.fixed_aw.beta.V_w_critical=V_w(aw_index,1); 

   % kappa at fixed aw beta phase
   alternate.fixed_aw.beta.phase=1;

[~,aw_index]=min(abs(a_w(index_phase_sep_end:end,:)-fixed_aw_val));
    aw_index=aw_index+index_phase_sep_end-1;
alternate.fixed_aw.alpha.aw=a_w(aw_index,1);

alternate.fixed_aw.alpha.SatCritical=SatRatio(aw_index,1);

alternate.fixed_aw.alpha.kappa_SatCritical=kappa(aw_index,1);
alternate.fixed_aw.alpha.kappa_a_w=a_w(aw_index,1);

alternate.fixed_aw.alpha.Dp_critical=Diameter_total(aw_index,1);
alternate.fixed_aw.alpha.V_w_critical=V_w(aw_index,1); 

    
%% max sat 
alternate.beta.phase=1;
    % kappa beta
    Sc_index=index_phase_sep_starts;

    alternate.beta.SatCritical=SatRatio(Sc_index,1);
    
    alternate.beta.kappa_SatCritical=kappa(Sc_index,1);
    alternate.beta.kappa_a_w=a_w(Sc_index,1);
    
    alternate.beta.Dp_critical=Diameter_total(Sc_index,1);
    alternate.beta.V_w_critical=V_w(Sc_index,1);
    
    % alpha
    [SatCritical2,Sc_index]=max(SatRatio(index_phase_sep_end:end,:));
    Sc_index=Sc_index+index_phase_sep_end-1;
    SatCritical=SatRatio(Sc_index,1);
    
    kappa_SatCritical=kappa(Sc_index,1);
    kappa_a_w=a_w(Sc_index,1);
    
    Dp_critical=Diameter_total(Sc_index,1);
    V_w_critical=V_w(Sc_index,1);

    

    
end

disp(['kappa (at Sc) = ' num2str(kappa_SatCritical) ' Dp cloud ' num2str(Dp_critical.*10^6) ' um'])

text_settings=['V eqv. diameter ' num2str(Dp_sim_nm) ' nm, O:C ' num2str(O2C_value) ', M ' num2str(Molar_mass) ' g/mol'];

disp(text_settings)
% figure
% plot(a_w,SatRatio,'.')
% xlim([.98,1.02])
% %ylim([.98,1.06])
% grid on
% ylabel('Sat ratio')
% xlabel('a_w')

save(['DATA for ' plot_name '.mat'])

paper_postion=[0, 0, 6.5, 3].*1;

figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', paper_postion+.5,'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]


axes1 = axes('Parent',figure1,...
    'Position',[0.114583333333333 0.19 0.355178571428571 0.78],...
    'LineWidth',1.75,...
    'FontSize',12);  %'CLim',[80 300],
set(axes1,'FontSize',12,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
hold on
grid on

% plot(Diameter_total.*10^6,SatRatio,'DisplayName','Unphysical','Parent',axes1,'LineWidth',2,...
%     'Color',[0.800000011920929 0 0.200000002980232])
plot(Diameter_total(1:index_phase_sep_starts,:).*10^6,SatRatio(1:index_phase_sep_starts,:),...
    'DisplayName','\beta-phase','Parent',axes1,'LineWidth',2,...
    'Color',[0.87058824300766 0.490196079015732 0])
plot(Diameter_total(index_phase_sep_end:end,:).*10^6,SatRatio(index_phase_sep_end:end,:),...
    'DisplayName','\alpha-phase','Parent',axes1,'LineWidth',2,...
    'Color',[0 0.447058826684952 0.74117648601532])

plot(Dp_critical.*10^6,SatCritical,'DisplayName','\alpha-phase critical','Parent',axes1,'MarkerSize',30,'Marker','.','LineStyle','none',...
    'Color',[0 0.447058826684952 0.74117648601532])
if alternate.beta.phase==1


    plot(alternate.beta.Dp_critical.*10^6,alternate.beta.SatCritical,'DisplayName','\beta-phase critical','Parent',axes1,'MarkerSize',30,'Marker','.','LineStyle','none',...
    'Color',[0.87058824300766 0.490196079015732 0])

    plot(Diameter_total([index_phase_sep_starts,index_phase_sep_end],:).*10^6,SatRatio([index_phase_sep_starts,index_phase_sep_end],:),'LineStyle',':',...
        'DisplayName','Miscibility gap','Parent',axes1,'LineWidth',2,...
     'Color',[0.800000011920929 0 0.200000002980232])
end



xlim([0,1])
%  ylim([SatCritical.*.99,SatCritical.*1.01])
 ylim([.995,1.01])

grid on
ylabel('Saturation ratio')
xlabel('d (\mum)')


axes2 = axes('Parent',figure1,...
    'Position',[0.623511904761905 0.19 0.352678571428569 0.78],...
    'LineWidth',1.75,...
    'FontSize',12);  %'CLim',[80 300],
set(axes2,'FontSize',12,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
hold on
grid on


% plot(a_w,kappa,'DisplayName','Unphysical','Parent',axes2,'LineWidth',2,...
%     'Color',[0.800000011920929 0 0.200000002980232])
plot(a_w(1:index_phase_sep_starts,:),kappa(1:index_phase_sep_starts,:),...
    'DisplayName','\beta-phase','Parent',axes2,'LineWidth',2,...
    'Color',[0.87058824300766 0.490196079015732 0])
plot(a_w(index_phase_sep_end:end,:),kappa(index_phase_sep_end:end,:),...
    'DisplayName','\alpha-phase','Parent',axes2,'LineWidth',2,...
    'Color',[0 0.447058826684952 0.74117648601532])

plot(kappa_a_w,kappa_SatCritical,'DisplayName','\alpha-phase critical','Parent',axes2,'MarkerSize',30,'Marker','.','LineStyle','none',...
    'Color',[0 0.447058826684952 0.74117648601532])
if alternate.beta.phase==1
    plot(alternate.beta.kappa_a_w,alternate.beta.kappa_SatCritical,'DisplayName','\beta-phase critical','Parent',axes2,'MarkerSize',30,'Marker','.','LineStyle','none',...
    'Color',[0.87058824300766 0.490196079015732 0])

    plot(a_w([index_phase_sep_starts,index_phase_sep_end],:),kappa([index_phase_sep_starts,index_phase_sep_end],:),'LineStyle',':',...
        'DisplayName','Miscibility gap','Parent',axes2,'LineWidth',2,...
     'Color',[0.800000011920929 0 0.200000002980232])
end





xlim([.98,1])
ylim([0,.12])
grid on
ylabel('\kappa_{CCN}')
xlabel('a_w')


annotation(figure1,'textbox',...
    [0.377763605442176 0.0267857142857143 0.32015306122449 0.08965934065934],...
    'String',{text_settings},...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off');


saveas(figure1, [plot_name '.png'], 'png')
saveas(figure1, [plot_name '.fig'], 'fig')