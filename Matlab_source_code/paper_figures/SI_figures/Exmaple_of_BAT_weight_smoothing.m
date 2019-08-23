%%
% Created by Kyle Gorkowski [GORKOWFALCON] on 2019-Aug-17  7:26 PM
% Copyright 2019 Kyle Gorkowski
%% plot for BAT coefficent weights to show how the smoothing works


molar_mass=500;
O2C_values=[0:0.001:0.6]';
Mr=18.016/molar_mass;


tran_lowO2C_fractionOne_phase=[0.189974476118418];
tran_lowO2C_sigmoid_bend=[79.2606902175984];
tran_lowO2C_sigmoid_shift=[0.0604293454322489];
tran_midO2C_sigmoid_bend=[75.0159268221068];
tran_midO2C_sigmoid_shift=[0.000947111285750515];

Onephase_O2C=single_phase_O2C_point_KGv3(Mr);
weights=zeros(length(O2C_values),3);

for i=1:length(O2C_values)
    O2C=O2C_values(i,1);
    % get region transition properties
    mid_transition=Onephase_O2C.*0.75;
    if O2C < mid_transition % lower to mid O2C region
        
        % data point trasfer weight
        O2C_1phase_delta=O2C-Onephase_O2C.*tran_lowO2C_fractionOne_phase;
        weight_1phase=1/(1+exp(-tran_lowO2C_sigmoid_bend.*(O2C_1phase_delta-tran_lowO2C_sigmoid_shift))); % logistic transfer function 1/(1+e^-(75*x))
        
        % normalize to end point so at mid_transition weight 2 is 1.
        O2C_1phase_delta_norm=O2C-mid_transition.*tran_lowO2C_fractionOne_phase;
        weight_1phase_norm=1/(1+exp(-tran_lowO2C_sigmoid_bend.*(O2C_1phase_delta_norm-tran_lowO2C_sigmoid_shift))); % logistic transfer function 1/(1+e^-(75*x))
        
        weight_1phase=weight_1phase./weight_1phase_norm;
        
        weights(i,2)=weight_1phase;
        
        weights(i,1)=1-weight_1phase;
        
        fitpar_name1='midO2C';
        fitpar_name2='lowO2C';
        
    elseif O2C < Onephase_O2C*2 % mid to high O2C region
        
        O2C_1phase_delta=O2C-Onephase_O2C;
        weight_1phase=1/(1+exp(-tran_midO2C_sigmoid_bend.*(O2C_1phase_delta-tran_midO2C_sigmoid_shift))); % logistic transfer function 1/(1+e^-(75*x))
        
                weights(i,3)=weight_1phase;

        weights(i,2)=1-weight_1phase;
        
        fitpar_name1='highO2C';
        fitpar_name2='midO2C';
        
    else % high only region
        
        
        fitpar_name1='highO2C';
        fitpar_name2=NaN;
        
        weights(i,2)=0;
        weights(i,3)=1;
    end
end

% figure
% plot(O2C_values, weights)


paper_mult=1;
fontsize=10;
plot_name=['weights example'];
paper_postion=[0, 0, 3.5, 3.5].*paper_mult;
figure1 = figure('Units', 'inches', 'PaperPosition', paper_postion, 'Position', paper_postion.*1+.5,'Color',[1 1 1]); %'Color',[0.9999 0.9999 0.9999]


axes1 = axes('Parent',figure1,...
    'Position',[0.16 0.16 0.75 0.80],...
    'LineWidth',1.75,...
    'FontSize',fontsize );  %'CLim',[80 300],
set(axes1,'FontSize',fontsize,'XMinorTick','on','YMinorTick','on',...
    'TickLength',[0.02 0.04],'TickDir','out','LineWidth',1.75);
% view(axes1,[175 35]);
grid(axes1,'on');
hold(axes1,'all');
hold on
% plot([0;1],[0;1],'k','DisplayName', '1:1', 'LineStyle', '-')
% 
% fill([0,1+mean_fractional_kappa_error,1-mean_fractional_kappa_error,0],[0,1,1,0],'r','FaceAlpha',0.25,...
%         'FaceColor',[0.65,0.65,0.65],...
%         'EdgeColor','none')
%     
    plot(O2C_values, weights(:,1),'DisplayName', 'Low O:C weights, \varpi_{low}',...
            'LineWidth',2,...
    'LineStyle','-',...
    'Color',[0.850980392156863 0.329411764705882 0.101960784313725])
    plot(O2C_values, weights(:,2),'DisplayName', 'Mid O:C weights, \varpi_{mid}',...
            'LineWidth',2,...
    'LineStyle','--',...
    'Color',[0.470588235294118 0.670588235294118 0.188235294117647])
    plot(O2C_values, weights(:,3),'DisplayName', 'High O:C weights, \varpi_{high}',...
            'LineWidth',2,...
    'LineStyle','-',...
    'Color',[0.074509803921569 0.623529411764706 1])

%     plot(ratio_data_AIOMFAC(:,1),ratio_data_AIOMFAC(:,2),'DisplayName', Rname_AIOMFAC,...
%             'LineWidth',1,...
%     'LineStyle',':',...
%     'Color',[0.149019607843137 0.149019607843137 0.149019607843137])
% 
%     % data
%     plot(measured_kappa,sorted_kappa_AIOMFAC(:,4),...
%      'MarkerFaceColor',[0.25 0.25 0.25],...
%     'MarkerEdgeColor',[0.0 0.0 0.0],...
%     'Marker','o',...
%         'MarkerSize',6,...
%     'LineStyle','none',...
%     'LineWidth',1, 'DisplayName', 'BAT')
% 
%     plot(measured_kappa,sorted_kappa_AIOMFAC(:,3),...
%      'MarkerFaceColor',[1 1 1],...
%     'MarkerEdgeColor',[0.0 0.0 0.0],...
%     'Marker','o',...
%     'MarkerSize',4,...
%     'LineStyle','none',...
%     'LineWidth',1, 'DisplayName', 'AIOMFAC')
% 

ylabel('Weight')
 xlabel('O:C')
legend1 = legend(axes1,'show');
set(legend1,'Location','northoutside', 'edgeColor', 'none');
saveas(figure1, [plot_name '.png'], 'png')
saveas(figure1, [plot_name '.fig'], 'fig')

