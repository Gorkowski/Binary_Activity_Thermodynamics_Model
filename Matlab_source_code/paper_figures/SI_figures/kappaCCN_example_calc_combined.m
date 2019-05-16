clear
% kappa test calc

plot_name=['figure kappaCCN combinded'];
% 
dots_size=20;

sets=3;
% set 1
load('DATA for figure kappaCCN mid O2C.mat', 'Dp_sim_nm', ...
    'Molar_mass', 'O2C_value', 'Diameter_total', 'SatRatio', 'kappa', 'a_w', 'alternate',...
    'index_phase_sep_starts', 'index_phase_sep_end','SatCritical','Dp_critical', 'kappa_a_w', 'kappa_SatCritical')
labelO2C(1)=O2C_value;
Diameter_total_combo(:,1)=Diameter_total;
SatRatio_combo(:,1)=SatRatio;
kappa_combo(:,1)=kappa;
a_w_combo(:,1)=a_w;
alternate_combo(1)=alternate;
index_phase_sep_starts_combo(1)=index_phase_sep_starts;
index_phase_sep_end_combo(1)=index_phase_sep_end;
SatCritical_combo(1)=SatCritical;
Dp_critical_combo(1)=Dp_critical;
kappa_a_w_combo(1)=kappa_a_w;
kappa_SatCritical_combo(1)=kappa_SatCritical;
text_settings=['V eqv. diameter ' num2str(Dp_sim_nm) ' nm, M ' num2str(Molar_mass) ' g/mol'];

% set 2
load('DATA for figure kappaCCN lower O2C.mat', 'Dp_sim_nm', ...
    'Molar_mass', 'O2C_value', 'Diameter_total', 'SatRatio', 'kappa', 'a_w', 'alternate',...
    'index_phase_sep_starts', 'index_phase_sep_end', 'SatCritical','Dp_critical', 'kappa_a_w', 'kappa_SatCritical')
labelO2C(:,2)=O2C_value;
Diameter_total_combo(:,2)=Diameter_total;
SatRatio_combo(:,2)=SatRatio;
kappa_combo(:,2)=kappa;
a_w_combo(:,2)=a_w;
alternate_combo(2)=alternate;
index_phase_sep_starts_combo(2)=index_phase_sep_starts;
index_phase_sep_end_combo(2)=index_phase_sep_end;
SatCritical_combo(2)=SatCritical;
Dp_critical_combo(2)=Dp_critical;
kappa_a_w_combo(2)=kappa_a_w;
kappa_SatCritical_combo(2)=kappa_SatCritical;
% set 3
load('DATA for figure kappaCCN zero O2C.mat', 'Dp_sim_nm', ...
    'Molar_mass', 'O2C_value', 'Diameter_total', 'SatRatio', 'kappa', 'a_w', 'alternate',...
    'index_phase_sep_starts', 'index_phase_sep_end','SatCritical','Dp_critical', 'kappa_a_w', 'kappa_SatCritical')
labelO2C(:,3)=O2C_value;
Diameter_total_combo(:,3)=Diameter_total;
SatRatio_combo(:,3)=SatRatio;
kappa_combo(:,3)=kappa;
a_w_combo(:,3)=a_w;
alternate_combo(3)=alternate;
index_phase_sep_starts_combo(3)=index_phase_sep_starts;
index_phase_sep_end_combo(3)=index_phase_sep_end;
SatCritical_combo(3)=SatCritical;
Dp_critical_combo(3)=Dp_critical;
kappa_a_w_combo(3)=kappa_a_w;
kappa_SatCritical_combo(3)=kappa_SatCritical;
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


for i=1:sets

plot(Diameter_total_combo(1:index_phase_sep_starts_combo(i),i).*10^6,SatRatio_combo(1:index_phase_sep_starts_combo(i),i),...
    'DisplayName','\beta-phase','Parent',axes1,'LineWidth',2,...
    'Color',[0.87058824300766 0.490196079015732 0])
plot(Diameter_total_combo(index_phase_sep_end_combo(i):end,i).*10^6,SatRatio_combo(index_phase_sep_end_combo(i):end,i),...
    'DisplayName',['\alpha-phase O2C' num2str(labelO2C)],'Parent',axes1,'LineWidth',2,...
    'Color',[0 0.447058826684952 0.74117648601532])

plot(Dp_critical_combo(i).*10^6,SatCritical_combo(i),'DisplayName','\alpha-phase critical','Parent',axes1,'MarkerSize',dots_size,'Marker','.','LineStyle','none',...
    'Color',[0 0.447058826684952 0.74117648601532])

if alternate_combo(i).beta.phase==1

    plot(alternate_combo(i).beta.Dp_critical.*10^6,alternate_combo(i).beta.SatCritical,'DisplayName','\beta-phase critical',...
        'Parent',axes1,'MarkerSize',dots_size,'Marker','.','LineStyle','none',...
    'Color',[0.87058824300766 0.490196079015732 0])

    plot(Diameter_total_combo([index_phase_sep_starts_combo(i),index_phase_sep_end_combo(i)],i).*10^6,SatRatio_combo([index_phase_sep_starts_combo(i),index_phase_sep_end_combo(i)],i),'LineStyle',':',...
        'DisplayName','Miscibility gap','Parent',axes1,'LineWidth',2,...
     'Color',[0.800000011920929 0 0.200000002980232])
end

end

xlim([0,1])
%  ylim([SatCritical.*.99,SatCritical.*1.01])
 ylim([.998,1.006])

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

for i=1:sets
% plot(a_w_combo(:,i),kappa_combo(:,i),'DisplayName','Unphysical','Parent',axes2,'LineWidth',2,...
%      'Color',[0.800000011920929 0 0.200000002980232])
plot(a_w_combo(1:index_phase_sep_starts_combo(i),i),kappa_combo(1:index_phase_sep_starts_combo(i),i),...
    'DisplayName','\beta-phase','Parent',axes2,'LineWidth',2,...
    'Color',[0.87058824300766 0.490196079015732 0])
plot(a_w_combo(index_phase_sep_end_combo(i):end,i),kappa_combo(index_phase_sep_end_combo(i):end,i),...
    'DisplayName','\alpha-phase','Parent',axes2,'LineWidth',2,...
    'Color',[0 0.447058826684952 0.74117648601532])

plot(kappa_a_w_combo(:,i),kappa_SatCritical_combo(i),'DisplayName','\alpha-phase critical','Parent',axes2,'MarkerSize',dots_size,'Marker','.','LineStyle','none',...
    'Color',[0 0.447058826684952 0.74117648601532])
if alternate_combo(i).beta.phase==1
    plot(alternate_combo(i).beta.kappa_a_w,alternate_combo(i).beta.kappa_SatCritical,'DisplayName','\beta-phase critical',...
        'Parent',axes2,'MarkerSize',dots_size,'Marker','.','LineStyle','none',...
    'Color',[0.87058824300766 0.490196079015732 0])

    plot(a_w_combo([index_phase_sep_starts_combo(i),index_phase_sep_end_combo(i)],i),kappa_combo([index_phase_sep_starts_combo(i),index_phase_sep_end_combo(i)],i),'LineStyle',':',...
        'DisplayName','Miscibility gap','Parent',axes2,'LineWidth',2,...
     'Color',[0.800000011920929 0 0.200000002980232])
end
end

xlim([.98,1])
ylim([-.005,.12])
grid on
ylabel('\kappa')
xlabel('a_w')


annotation(figure1,'textbox',...
    [0.377763605442176 0.0267857142857143 0.32015306122449 0.08965934065934],...
    'String',{text_settings},...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off');


saveas(figure1, [plot_name '.png'], 'png')
saveas(figure1, [plot_name '.fig'], 'fig')