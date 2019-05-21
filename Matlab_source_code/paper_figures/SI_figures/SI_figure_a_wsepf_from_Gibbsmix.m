clear 

O2C=.2250;
H2C=0;
Morg=[100]; 
Mr=18.016./Morg;
mole_frac_scan=[1:-0.00001:0]';


mode1='hydroxyl';
[func1, func2, ycal_water, ycalc_org, activity_water, activity_org, mass_fraction1, mass_fraction2,Gibbs_RT, dGibbs_RTdx2]...
    =BAT_properties_calculation_v1(mole_frac_scan, O2C,  H2C, Mr,mode1,[]);

[ideal_Gibbs_RT,ideal_dGibbs_dxRT, Gibbs_total_RT,dGibbs_total_dxRT] = Gibbs_mix_total(mole_frac_scan, Gibbs_RT, dGibbs_RTdx2);

% check for phase seperation
[phaseSep_via_activity_w,phaseSep_via_activity_curvature_w,index_phase_sep_starts_w,index_phase_sep_end_w]=...
    finds_PhaseSep_and_activity_curve_dips_v2(activity_water);
[phaseSep_via_activity_org,phaseSep_via_activity_curvature_org,index_phase_sep_starts_org,index_phase_sep_end_org]=...
    finds_PhaseSep_and_activity_curve_dips_v2(activity_org);

[phase_sep_check,lower_a_w_sep_index,upper_a_w_sep_index, matching_Upper_a_w_sep_index]=finds_PhaseSep_w_and_org(activity_water,activity_org);

disp(activity_water(upper_a_w_sep_index)-activity_water(matching_Upper_a_w_sep_index))
disp(activity_water(upper_a_w_sep_index))
disp(activity_water([matching_Upper_a_w_sep_index-1,matching_Upper_a_w_sep_index,matching_Upper_a_w_sep_index+1]))
disp(activity_water(lower_a_w_sep_index))

if phaseSep_via_activity_curvature_org>0
    % index_phase_sep_end_org=index_phase_sep_end_org-2;
    a_org_sep_i=[index_phase_sep_starts_org;index_phase_sep_end_org];
    a_w_sep_i=[index_phase_sep_starts_w;index_phase_sep_end_w];
else
    a_org_sep_i=[1;2];
    a_w_sep_i=[1;2];
end

G_sep(1,1)=lower_a_w_sep_index;
G_sep(2,1)=upper_a_w_sep_index;

x_variable=1-mole_frac_scan;



%% plot 
line_width=1;
fontsize=12;
paper_mult=1;

% org mass
plot_name=['SI-GibbsTotal_sep'];
paper_postion=[0, 0, 6.5, 4].*paper_mult;
figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', paper_postion.*1+.5,'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]


axes1 = axes('Parent',figure1,...
    'Position',[0.14 0.58 0.8 0.40],...
    'LineWidth',1.75,...
    'FontSize',fontsize, 'Xcolor', 'none');  %'CLim',[80 300],
set(axes1,'FontSize',fontsize,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% view(axes1,[175 35]);
grid(axes1,'on');
hold(axes1,'all');
hold on

plot(x_variable,Gibbs_total_RT,'DisplayName', 'Gibbs','LineWidth', line_width,...
    'Color',[0 0 0])

data_match=Gibbs_total_RT;
% plot(x_variable(a_w_sep_i,:),data_match(a_w_sep_i,:),':k', 'DisplayName', 'a_{w} sep')
% plot(x_variable(a_org_sep_i,:),data_match(a_org_sep_i,:),'k', 'DisplayName', 'a_{org} sep')
plot(x_variable(G_sep,:),data_match(G_sep,:),'r',  'DisplayName', 'G sep','LineWidth', line_width,...
    'LineStyle','--',...
    'Color',[1 0 0])

plot(x_variable(matching_Upper_a_w_sep_index,:),data_match(matching_Upper_a_w_sep_index,:),'square', ...
    'Color',[0 0.498039215803146 0], 'DisplayName', 'metastable','LineWidth', line_width)
plot(x_variable(index_phase_sep_end_w,:),data_match(index_phase_sep_end_w,:),'square',  ...
    'Color',[0 0.447058826684952 0.74117648601532], 'DisplayName', 'metastable','LineWidth', line_width)

plot(x_variable(lower_a_w_sep_index,:),data_match(lower_a_w_sep_index,:),'o',  ...
    'Color',[1 0 0], 'DisplayName', 'beta end','LineWidth', line_width)
plot(x_variable(upper_a_w_sep_index,:),data_match(upper_a_w_sep_index,:),'o',...
    'Color',[1 0 0],  'DisplayName', 'alpha end','LineWidth', line_width)


% xlabel('x_w')
ylabel('Gibbs/RT mol')

axes2 = axes('Parent',figure1,...
    'Position',[0.14 0.14 0.8 0.40],...
    'LineWidth',1.75,...
    'FontSize',fontsize);  %'CLim',[80 300],
set(axes2,'FontSize',fontsize,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% view(axes1,[175 35]);
grid(axes2,'on');
hold(axes2,'all');
hold on

plot(x_variable,activity_water,'LineWidth', line_width, 'Color',[0 0.447058826684952 0.74117648601532])
plot(x_variable,activity_org,'LineWidth', line_width,'Color',[0 0.498039215803146 0])
plot(x_variable(a_w_sep_i,:),activity_water(a_w_sep_i,:),'LineStyle','--',...
    'Color',[0 0.447058826684952 0.74117648601532])
plot(x_variable(a_org_sep_i,:),activity_org(a_org_sep_i,:),'LineStyle','--',...
    'Color',[0 0.498039215803146 0])

xlabel('x')
ylabel('a_w or a_{org}')

linkaxes([axes2, axes1],'x')

saveas(figure1, [plot_name '.png'], 'png')
saveas(figure1, [plot_name '.fig'], 'fig')
