%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2018-Nov-26  9:23 AM
% Copyright 2018 Kyle Gorkowski
%%

clear

%% settings
O2C_values=[0:0.0025:1.2]';
O2C_values_plot=[0.35, 0.9];

a_w_specific=[0.999, 0.70 ,0.40]';
molar_masses_plot=[125,150,200,300]';
molar_masses=[80:5:600]';

kappa_CCN_settings=get_default_kappa_CCN_settings;

%% program
mole_frac_scan=[1:-0.0001:0]';

fit_tolerance=10^-6;


H2C=0;
progressbartext('activity space calc')
for i=1:length(O2C_values)
    for m_i=1:length(molar_masses)
        fix_molarmass1=18.016/molar_masses(m_i,1);
        % hydroxyl
        mode1='hydroxyl';
        [func1, func2, ycal_water, ycalc_org, activity_water, activity_org, mass_fraction1, mass_fraction2,Gibbs_RT, dGibbs_RTdx2]...
            =BAT_properties_calculation_v1(mole_frac_scan, O2C_values(i,1),  O2C_values(i,1).*H2C, fix_molarmass1,mode1,[]);
        
        density_org_g_cm3=Org_density_Estimate_KGv1(18.016/fix_molarmass1, O2C_values(i), O2C_values(i).*0);
        
        [kappa_SatCritical_OH(i,m_i), kappa_a_w_OH(i,m_i), SatCritical_OH(i,m_i), Dp_critical(i,m_i), V_w_critical(i,m_i), alternate]= ...
            kappa_critical_calc_for_McGlashan_v3(...
    density_org_g_cm3,  mass_fraction1,  mass_fraction2, activity_water, kappa_CCN_settings);

        if alternate.beta.phase==0
            kappa_SatCritical_OH_beta(i,m_i)=NaN;
                        kappa_SatCritical_OH_alpha(i,m_i)=NaN;

        else
            kappa_SatCritical_OH_beta(i,m_i)=alternate.beta.kappa_SatCritical;
            kappa_SatCritical_OH_alpha(i,m_i)=alternate.alpha.kappa_SatCritical;
        end
        
        % hydroperoxide sim
        mode2='hydroperoxide';
        [func1, func2, ycal_water, ycalc_org, activity_water, activity_org, mass_fraction1, mass_fraction2,Gibbs_RT, dGibbs_RTdx2]...
            =BAT_properties_calculation_v1(mole_frac_scan, O2C_values(i,1),  O2C_values(i,1).*H2C, fix_molarmass1,mode2,[]);
        
        density_org_g_cm3=Org_density_Estimate_KGv1(18.016/fix_molarmass1, O2C_values(i), O2C_values(i).*0);
        
%         [kappa_SatCritical_OOH(i,m_i), kappa_a_w_OOH(i,m_i), SatCritical_OOH(i,m_i), Dp_critical(i,m_i), V_w_critical(i,m_i), alternate]= ...
%             kappa_critical_calc_for_McGlashan_v2(Dp_sim_nm, sigma_droplet_N_m,...
%             density_org_g_cm3, mass_fraction1,  mass_fraction2, activity_water);
                [kappa_SatCritical_OOH(i,m_i), kappa_a_w_OOH(i,m_i), SatCritical_OOH(i,m_i), Dp_critical(i,m_i), V_w_critical(i,m_i), alternate]= ...
            kappa_critical_calc_for_McGlashan_v3(...
    density_org_g_cm3,  mass_fraction1,  mass_fraction2, activity_water, kappa_CCN_settings);

        if alternate.beta.phase==0
            kappa_SatCritical_OOH_beta(i,m_i)=NaN;
            kappa_SatCritical_OOH_alpha(i,m_i)=NaN;
        else
            kappa_SatCritical_OOH_beta(i,m_i)=alternate.beta.kappa_SatCritical;
            kappa_SatCritical_OOH_alpha(i,m_i)=alternate.alpha.kappa_SatCritical;

        end
        
        
        
        kappa_intrinsic(i,m_i)=density_org_g_cm3./1.*(fix_molarmass1);
        
    end
    progressbartext( i/length(O2C_values))
    
end

%% data
alpha_pinene_data=[[0.0102100000000000,0.429930000000000;0.104530000000000,0.385170000000000;0.101560000000000,0.407490000000000;0.105860000000000,0.411330000000000;0.0999100000000000,0.438530000000000;0.103710000000000,0.501070000000000;0.101140000000000,0.389680000000000;0.0991900000000000,0.408210000000000;0.104100000000000,0.423430000000000;0.0983900000000000,0.442300000000000;0.106740000000000,0.499970000000000;0.0954400000000000,0.394540000000000;0.0939800000000000,0.409670000000000;0.0961600000000000,0.427120000000000;0.0968200000000000,0.447960000000000;0.104920000000000,0.495780000000000;0.0792600000000000,0.399270000000000;0.0839200000000000,0.413720000000000;0.0865200000000000,0.432690000000000;0.0953600000000000,0.449840000000000;0.104500000000000,0.502600000000000;0.103960000000000,0.338020000000000;0.106930000000000,0.340810000000000;0.108120000000000,0.349170000000000;0.111680000000000,0.352890000000000;0.106340000000000,0.354750000000000;0.104550000000000,0.358160000000000;0.104550000000000,0.361260000000000;0.105150000000000,0.364360000000000;0.111680000000000,0.360330000000000;0.114650000000000,0.360020000000000;0.109310000000000,0.366840000000000;0.111680000000000,0.373660000000000;0.101580000000000,0.372730000000000;0.106340000000000,0.373350000000000;0.108120000000000,0.377070000000000;0.108120000000000,0.383570000000000;0.108710000000000,0.386360000000000;0.102770000000000,0.387290000000000;0.0974300000000000,0.389150000000000;0.105740000000000,0.393490000000000;0.105150000000000,0.404340000000000;0.100990000000000,0.402480000000000;0.0974300000000000,0.401860000000000;0.0944600000000000,0.404960000000000;0.0962400000000000,0.408680000000000;0.0938600000000000,0.410540000000000;0.100400000000000,0.412090000000000;0.105150000000000,0.411780000000000;0.105150000000000,0.408370000000000;0.112280000000000,0.408680000000000;0.0968300000000000,0.415190000000000;0.0998000000000000,0.417050000000000;0.110500000000000,0.425100000000000;0.106930000000000,0.422310000000000;0.102180000000000,0.421380000000000;0.103960000000000,0.423550000000000;0.0974300000000000,0.425720000000000;0.102770000000000,0.435950000000000;0.106340000000000,0.444940000000000;0.112280000000000,0.452070000000000;0.110500000000000,0.454860000000000;0.106340000000000,0.460120000000000;0.107530000000000,0.466940000000000;0.110500000000000,0.469110000000000;0.112280000000000,0.467560000000000;0.123560000000000,0.476550000000000;0.119410000000000,0.471900000000000;0.118810000000000,0.473760000000000;0.115840000000000,0.472520000000000;0.104550000000000,0.476240000000000;0.115250000000000,0.482130000000000;0.122970000000000,0.483680000000000;0.125350000000000,0.486160000000000;0.118810000000000,0.486470000000000;0.118220000000000,0.484300000000000;0.113470000000000,0.485850000000000;0.119410000000000,0.489570000000000;0.125940000000000,0.493280000000000;0.124160000000000,0.496070000000000;0.108120000000000,0.488330000000000;0.116440000000000,0.496380000000000;0.117030000000000,0.501960000000000;0.120000000000000,0.480000000000000;0.100000000000000,0.550000000000000;0.100000000000000,0.520000000000000;0.100000000000000,0.500000000000000;0.100000000000000,0.500000000000000;0.110000000000000,0.530000000000000;0.110000000000000,0.570000000000000;0.100000000000000,0.440000000000000;0.100000000000000,0.590000000000000;0.110000000000000,0.610000000000000;0.100000000000000,0.480000000000000;0.100000000000000,0.440000000000000;0.110000000000000,0.410000000000000;0.100000000000000,0.400000000000000;0.110000000000000,0.380000000000000;0.150080000000000,0.405830000000000;0.153290000000000,0.467060000000000;0.184920000000000,0.692710000000000;0.223510000000000,0.839650000000000;0.134530000000000,0.398830000000000;0.149540000000000,0.495040000000000;0.210640000000000,0.729450000000000;0.250840000000000,0.878130000000000]];

%plot
plot_name=['paper figure kappaCCN'];
paper_postion=[0, 0, 6.5, 3].*1;

figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', paper_postion+.5,'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]


axes1 = axes('Parent',figure1,...
    'Position',[0.114583333333333 0.19 0.465178571428571 0.78],...
    'LineWidth',1.75,...
    'FontSize',12);  %'CLim',[80 300],
set(axes1,'FontSize',12,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
hold on
grid on
colorOrder = get(gca, 'ColorOrder');
lines_totoal=length(molar_masses_plot);
colorOrder=repmat(colorOrder, ceil(lines_totoal./length(colorOrder)),1);

for m_i=1:length(molar_masses_plot)
    [~,plot_index]=min(abs(molar_masses_plot(m_i)-molar_masses));
    plot(O2C_values,kappa_SatCritical_OH(:,plot_index)','DisplayName', [mode1 ' ' num2str(molar_masses_plot(m_i)) ' g/mol' ],...
        'Color',colorOrder(m_i,:), 'LineWidth', 2)
        plot(O2C_values,kappa_SatCritical_OH_alpha(:,plot_index)','DisplayName', [mode1 ' beta ' num2str(molar_masses_plot(m_i)) ' g/mol' ],...
        'Color',colorOrder(m_i,:), 'LineWidth', 1)
    
    % miscibility gap fill.
    [plot_top,xvalues]=Removes_nan_rows_dual(kappa_SatCritical_OH_alpha(:,plot_index),O2C_values);
        [plot_bottom,~]=Removes_nan_rows_dual(kappa_SatCritical_OH_beta(:,plot_index),O2C_values);

    fill(axes1,[xvalues;flip(xvalues)],[plot_top;flip(plot_bottom)],'r','FaceAlpha',0.25,...
        'FaceColor',[0.63921568627451 0.07843137254902 0.180392156862745],...
        'EdgeColor','none');
    
    %
    if m_i==4 %|| m_i==2
        plot(O2C_values,kappa_SatCritical_OOH(:,plot_index)','--', 'DisplayName', [mode2 ' ' num2str(molar_masses_plot(m_i)) ' g/mol' ],...
            'Color',colorOrder(m_i,:), 'LineWidth', 2)
        plot(O2C_values,kappa_SatCritical_OOH_alpha(:,plot_index)','--', 'DisplayName', [mode2 ' beta ' num2str(molar_masses_plot(m_i)) ' g/mol' ],...
            'Color',colorOrder(m_i,:), 'LineWidth', 1)
        
        plot(O2C_values,kappa_intrinsic(:,plot_index)',':', 'DisplayName', ['\kappa intr ' num2str(molar_masses_plot(m_i)) ' g/mol' ],...
            'Color',colorOrder(m_i,:), 'LineWidth', 2)
        
            % miscibility gap fill.
    [plot_top,xvalues]=Removes_nan_rows_dual(kappa_SatCritical_OOH_alpha(:,plot_index),O2C_values);
        [plot_bottom,~]=Removes_nan_rows_dual(kappa_SatCritical_OOH_beta(:,plot_index),O2C_values);

    fill(axes1,[xvalues;flip(xvalues)],[plot_top;flip(plot_bottom)],'r','FaceAlpha',0.25,...
        'FaceColor',[0.63921568627451 0.07843137254902 0.180392156862745],...
        'EdgeColor','none');
    
    
    end
end
point_size=20;

scatter(alpha_pinene_data(:,2),alpha_pinene_data(:,1),point_size,'DisplayName','\alpha Pinene SOA',...
    'MarkerFaceColor',[1 0.074509803921569 0.650980392156863],...
    'MarkerEdgeColor',[0 0 0],...
    'LineWidth',2);


xlim([.0,1.2])
ylim([0,.25])

xlabel('O:C atmoic ratio')
ylabel('\kappa_{CCN}')
% legend(axes1,'show');

axes2 = axes('Parent',figure1,...
    'Position',[0.623511904761905 0.19 0.352678571428569 0.78],...
    'LineWidth',1.75,...
    'FontSize',12);  %'CLim',[80 300],
set(axes2,'FontSize',12,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
hold on
grid on
for m_i=1:length(O2C_values_plot)
    [~,plot_index]=min(abs(O2C_values_plot(m_i)-O2C_values));
    
    plot(molar_masses,kappa_SatCritical_OH(plot_index,:)','DisplayName', [mode1 ' ' num2str(O2C_values_plot(m_i)) ' O:C' ],...
        'Color',colorOrder(m_i,:), 'LineWidth', 2)
    plot(molar_masses,kappa_SatCritical_OH_alpha(plot_index,:)','DisplayName', [mode1 ' beta ' num2str(O2C_values_plot(m_i)) ' O:C' ],...
        'Color',colorOrder(m_i,:), 'LineWidth', 1)
    
    plot(molar_masses,kappa_SatCritical_OOH(plot_index,:)','--','DisplayName', [mode2 ' ' num2str(O2C_values_plot(m_i)) ' O:C' ],...
        'Color',colorOrder(m_i,:), 'LineWidth', 2)
    plot(molar_masses,kappa_SatCritical_OOH_alpha(plot_index,:)','--','DisplayName', [mode2 ' beta ' num2str(O2C_values_plot(m_i)) ' O:C' ],...
        'Color',colorOrder(m_i,:), 'LineWidth', 1)
    
    plot(molar_masses,kappa_intrinsic(plot_index,:)',':', 'DisplayName', ['\kappa intr ' num2str(O2C_values_plot(m_i)) ' O:C' ],...
        'Color',colorOrder(m_i,:), 'LineWidth', 2)
    
    
        % miscibility gap fill.
    [plot_top,xvalues]=Removes_nan_rows_dual(kappa_SatCritical_OH_alpha(plot_index,:)',molar_masses);
        [plot_bottom,~]=Removes_nan_rows_dual(kappa_SatCritical_OH_beta(plot_index,:)',molar_masses);

    fill(axes2,[xvalues;flip(xvalues)],[plot_top;flip(plot_bottom)],'r','FaceAlpha',0.25,...
        'FaceColor',[0.63921568627451 0.07843137254902 0.180392156862745],...
        'EdgeColor','none');
    
        
        % miscibility gap fill.
    [plot_top,xvalues]=Removes_nan_rows_dual(kappa_SatCritical_OOH_alpha(plot_index,:)',molar_masses);
        [plot_bottom,~]=Removes_nan_rows_dual(kappa_SatCritical_OOH_beta(plot_index,:)',molar_masses);

    fill(axes2,[xvalues;flip(xvalues)],[plot_top;flip(plot_bottom)],'r','FaceAlpha',0.25,...
        'FaceColor',[0.63921568627451 0.07843137254902 0.180392156862745],...
        'EdgeColor','none');
    
end
xlabel('Molar Mass')

ylim([0,.25])
xlim([80,600])

saveas(figure1, [plot_name '.png'], 'png')
saveas(figure1, [plot_name '.fig'], 'fig')





