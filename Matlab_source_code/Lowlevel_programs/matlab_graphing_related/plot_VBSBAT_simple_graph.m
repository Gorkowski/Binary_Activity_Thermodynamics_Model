function plot_VBSBAT_simple_graph(details, plot_mode, sim_name)
%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2019-Mar-27 10:56 AM
% Copyright 2019 Kyle Gorkowski 
%% 
% clear
% load('VBSBAT_sim_aPsoa_t4.mat')
% sim_name='test';
% plot_mode='yes';

%% plot
paper_postion=[0, 0, 6.5, 4].*1;
aw_series=details.inputs.aw_series;
line_width=2;

if contains(plot_mode, 'only') 
    figure1 = figure('Name',['VBSBAT: ' sim_name], 'Units', 'inches', ...
       'Position', paper_postion+.5,'Color',[1 1 1],'visible','off'); %'Color',[0.9999 0.9999 0.9999]
else
    figure1 = figure('Name',['VBSBAT: ' sim_name], 'Units', 'inches',...
        'Position', paper_postion+.5,'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]
end

axes1 = axes('Parent',figure1,...
    'Position',[0.114583333333333 0.19 0.355178571428571 0.60],...
    'LineWidth',1.75,...
    'FontSize',12);  %'CLim',[80 300],
set(axes1,'FontSize',12,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% view(axes1,[175 35]);
grid(axes1,'on');
hold on
plot(axes1,aw_series,details.totals.C_OA_PM,'DisplayName', 'OA total liquid','LineWidth',line_width+2,...
    'Color',[0.074509803921569 0.666666666666667 0.372549019607843]) % ,C_wOA_total,C_OA_out+C_wOA_total
plot(axes1,aw_series,details.totals.Coa_alpha,'DisplayName', 'OA \alpha phase','LineWidth',line_width,'LineStyle','--',...
    'Color',[0.564705882352941 0.325490196078431 0.631372549019608]) % ,C_wOA_total,C_OA_out+C_wOA_total
plot(axes1,aw_series,details.totals.Coa_beta,'DisplayName', 'OA \beta phase','LineWidth',line_width,'LineStyle','-.',...
    'Color',[0.913725490196078 0.615686274509804 0.231372549019608]) % ,C_wOA_total,C_OA_out+C_wOA_total
xlabel('Water activity (a_w)')
ylabel('Organic PM mass (\mug/m^3)')
legend1=legend(axes1,'show');
set(legend1, 'Position',[0.186755954221423 0.836033949892922 0.19940475822382 0.1296296261113],...
    'EdgeColor',[0.901960784313726 0.901960784313726 0.901960784313726],...
    'FontSize',10);

% water axes
axes2 = axes('Parent',figure1,...
    'Position',[0.623511904761905 0.19 0.352678571428569 0.60],...
    'LineWidth',1.75,...
    'FontSize',12);  %'CLim',[80 300],
set(axes2,'FontSize',12,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
grid(axes2,'on');
hold on
plot(axes2,aw_series,details.totals.Caq_PM,'DisplayName', 'Water total liquid','LineWidth',line_width+2,...
    'Color',[0.105882352941176 0.458823529411765 0.737254901960784]) % ,C_wOA_total,C_OA_out+C_wOA_total
plot(axes2,aw_series,sum(details.species_specific.Caq_j_alpha,2),'DisplayName', 'Water \alpha phase','LineWidth',line_width,'LineStyle','--',...
    'Color',[0.564705882352941 0.325490196078431 0.631372549019608]) % ,C_wOA_total,C_OA_out+C_wOA_total
plot(axes2,aw_series,sum(details.species_specific.Caq_j_beta,2),'DisplayName', 'Water \beta phase','LineWidth',line_width,'LineStyle','-.',...
    'Color',[0.913725490196078 0.615686274509804 0.231372549019608]) % ,C_wOA_total,C_OA_out+C_wOA_total
xlabel('Water activity (a_w)')
ylabel('Water PM mass (\mug/m^3)')

ylim([min(details.totals.Caq_PM),max(max(details.totals.Coa_beta),max(details.totals.Coa_alpha))])
legend2=legend(axes2,'show');
set(legend2,'Position',[0.687500002217434 0.8395061523876 0.224702376517512 0.1296296261113],...
    'EdgeColor',[0.901960784313726 0.901960784313726 0.901960784313726],...
    'FontSize',10);


annotation(figure1,'textbox',...
    [0.442964285714286 0.868055555555555 0.195428571428571 0.115740740740741],...
    'VerticalAlignment','middle',...
    'String',{sim_name},...
    'Interpreter','none',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'EdgeColor','none');


if not(contains(plot_mode, 'only pdf'))
    saveas(figure1,['PM_organic_Mass_' sim_name '.png'],'png')
end
if not(contains(plot_mode, 'only png'))
    saveas(figure1,['PM_organic_Mass_' sim_name '.pdf'],'pdf')
end

if not(contains(plot_mode, 'only'))
    saveas(figure1,['PM_organic_Mass_' sim_name '.fig'],'fig')
else
    close(figure1)
end








