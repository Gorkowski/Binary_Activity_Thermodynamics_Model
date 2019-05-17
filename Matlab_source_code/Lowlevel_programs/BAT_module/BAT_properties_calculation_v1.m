function [ln_func1, ln_func2, ycalc1, ycalc2, activity_calc1, activity_calc2, mass_fraction1, mass_fraction2, Gibbs_RT, dGibbs_RTdx2]=BAT_properties_calculation_v1(...
    org_mole_fraction, O2C, H2C, molarmass_ratio, BAT_functional_group, special_options, N2C_values_denistyOnly)

%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2018-Mar-20  9:41 AM
% Copyright 2018 Kyle Gorkowski
% % based off of Andi's mcglashan code
% model types {'hydroxyl';'carboxyl';'ether';'ketone';'ester';'hydroperoxide';'hydroperoxideSOA';'PEG'};
% Updated by Kyle Gorkowski [GORKOWFALCON] on 2018-Dec-16 10:38 AM
% % re-named and cleaned up for BAT v1

% % default O2C low to mid to high transition logistic function coefficients
% tran_lowO2C_fractionOne_phase=0.219;
% tran_lowO2C_sigmoid_bend=80;
% tran_lowO2C_sigmoid_shift=0.0472066;
% tran_midO2C_sigmoid_bend=75;
% tran_midO2C_sigmoid_shift=0.00094;
% default v2 O2C low to mid to high transition logistic function coefficients
tran_lowO2C_fractionOne_phase=[0.189974476118418];
tran_lowO2C_sigmoid_bend=[79.2606902175984];
tran_lowO2C_sigmoid_shift=[0.0604293454322489];
tran_midO2C_sigmoid_bend=[75.0159268221068];
tran_midO2C_sigmoid_shift=[0.000947111285750515];
% % default v2 O2C low to mid to high transition logistic function
% coefficients 3-11-2019
% tran_lowO2C_fractionOne_phase=[0.315947592234889];
% tran_lowO2C_sigmoid_bend=[79.982357891392310];
% tran_lowO2C_sigmoid_shift=[0.033217567954200];
% tran_midO2C_sigmoid_bend=[13.193243720479420];
% tran_midO2C_sigmoid_shift=[0.048246892288302];

if not(exist('N2C_values_denistyOnly'))
    N2C_values_denistyOnly=molarmass_ratio.*0;
end


if isempty(special_options)
    
elseif strcmpi(special_options.fit_option,'lowO2C_bound_fit' ) % for transfer function fitting of logistic function coefficients
    
    tran_lowO2C_fractionOne_phase=special_options.tran_lowO2C_fractionOne_phase;
    tran_lowO2C_sigmoid_bend=special_options.tran_lowO2C_simoid_bend;
    tran_lowO2C_sigmoid_shift=special_options.tran_lowO2C_simoid_shift;
    
elseif strcmpi(special_options.fit_option,'midO2C_bound_fit' )
    
    tran_midO2C_sigmoid_bend=special_options.tran_midO2C_simoid_bend;
    tran_midO2C_sigmoid_shift=special_options.tran_midO2C_simoid_shift;
    
end

% shift for equlivlent O2C and MW, for activity calc
Mr_massfrac_final=molarmass_ratio;% used in final mass fraction calc and not for activity coefficent calc.
[O2C,molarmass_ratio] = convert_chemical_structure_to_OH_eqv_v3(O2C,molarmass_ratio, BAT_functional_group);

Smole=size(org_mole_fraction);

%% force org mole fraction to be 1 or less
x2=org_mole_fraction;
mole_fraction_mask=x2>1;
x2=x2.*not(mole_fraction_mask)+mole_fraction_mask;
% greater than zero force
mole_fraction_mask=x2<0;
x2=x2.*not(mole_fraction_mask);

%% replace infinite dilution to of zero to small value
x2=replace_data_A_to_B_KGv1(x2,0,10^-20); % infinite dilution limit
Mr=molarmass_ratio; % molar1/molar2

% properties
[densityEst]=Org_density_Estimate_KGv1(18.016/Mr, O2C, H2C, N2C_values_denistyOnly);
Onephase_O2C=single_phase_O2C_point_KGv3(Mr);

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
    
    weight_2phase=1-weight_1phase;
    
    fitpar_name1='midO2C';
    fitpar_name2='lowO2C';

elseif O2C < Onephase_O2C*2 % mid to high O2C region
    
    O2C_1phase_delta=O2C-Onephase_O2C;
    weight_1phase=1/(1+exp(-tran_midO2C_sigmoid_bend.*(O2C_1phase_delta-tran_midO2C_sigmoid_shift))); % logistic transfer function 1/(1+e^-(75*x))
    
    weight_2phase=1-weight_1phase;
    
    fitpar_name1='highO2C';
    fitpar_name2='midO2C';
    
else % high only region
    
    
    fitpar_name1='highO2C';
    fitpar_name2=NaN;
    
    weight_2phase=0;
    weight_1phase=1;
end


% % fit properties data
fitpar_lowO2C=[7.089476E+00    -6.226781E-01    -7.711860E+00    -1.000000E+02    -3.885941E+01     3.081244E-09    -1.000000E+02     6.188812E+01    -5.988895E+00     6.940689E+00];
fitpar_midO2C=[ 5.872214E+00    -9.740486E-01    -4.535007E+00    -1.000000E+02    -5.129327E+00     2.109751E+00    -2.809232E+01    -2.367683E+01    -1.219164E+00     4.742729E+00];
fitpar_highO2C= [ 5.921550E+00    -1.000000E+02    -2.528295E+00    -1.000000E+02    -3.883017E+00     1.353916E+00    -7.898128E+00    -1.160145E+01    -7.868187E-02     3.650860E+00];

%% fit properties data, after I cleaned out the couple non-hydroxyl molecules (resulted in overfitting of data set) 3-11-2019
% fitpar_lowO2C=[ 6.722799E+00    -1.503450E+00    -5.956378E+00    -1.000000E+02    -2.755372E+01     4.837127E+00    -7.415255E+01    -4.627769E+01    -1.886588E+00     5.524567E+00];%Cumulative obj. func. val. of fit data subset  :    4.533449E+03
% fitpar_midO2C=[   6.293654E+00    -1.498720E+00    -5.615401E+00    -1.000000E+02    -5.672471E+00     2.349137E+00    -3.456994E+01    -3.003757E+01    -1.485205E+00     4.796483E+00];%Cumulative obj. func. val. of fit data subset  :    2.547947E+03
% fitpar_highO2C= [ 5.392588E+00     2.689195E+00    -2.658357E+00    -1.244450E+00    -3.993080E+00    -4.628953E-01    -8.389237E+00     3.151726E+00     8.374270E-02     4.773770E+00]; %Cumulative obj. func. val. of fit data subset  :    7.776732E+02



% get fit parameters in into correct phases phase 1
if strcmpi(fitpar_name1, 'highO2C')
    fitpar_1phase=fitpar_highO2C;
elseif strcmpi(fitpar_name1, 'midO2C')
    fitpar_1phase=fitpar_midO2C;
elseif strcmpi(fitpar_name1, 'lowO2C')
    fitpar_1phase=fitpar_lowO2C;
end

% get fit parameters in into correct phases phase 2
if strcmpi(fitpar_name2, 'highO2C')
    fitpar_2phase=fitpar_highO2C;
elseif strcmpi(fitpar_name2, 'midO2C')
    fitpar_2phase=fitpar_midO2C;
    
elseif strcmpi(fitpar_name2, 'lowO2C')
    fitpar_2phase=fitpar_lowO2C;
end

% for biphasic line calcs
if strcmpi(BAT_functional_group(1,:),'initial fitting for biphasic transition')
    
    %     %inital fit to bihpasic line
         fitpar_1phase=[ 5.885109E+00    -9.849013E-01    -4.731250E+00    -6.227207E+00    -5.201652E+00     2.320286E+00    -3.082297E+01    -2.584037E+01    -1.237227E+00     4.069905E+00];
    %
    %inital fit to fing bihpasic line with cleaned up input chemicals 3-11-2019
    %removed a few non-hydroxyl molecules.  obj. func. val. of fit data subset  :    2.957367E+03
%     fitpar_1phase=[6.396595E+00    -1.351389E+00    -5.117799E+00    -1.000000E+02    -3.790306E+00     1.232127E+00    -2.121294E+01    -1.176860E+01    -1.309129E+00     4.924529E+00];

    weight_2phase=0;
    weight_1phase=1;
end

gemix=zeros(Smole(1,1),2);
dgemix2dx=zeros(Smole(1,1),2);
do_calc=false;

for loop_i=1:2
    
    if loop_i==1
        
        if weight_1phase==0
            do_calc=false;
        else
            fitpar=fitpar_1phase;
            do_calc=true;
        end
    else
        
        if weight_2phase==0
            do_calc=false;
        else
            fitpar=fitpar_2phase;
            do_calc=true;
        end
        
    end
    
    
    if do_calc
        n=length(fitpar);
        
        rhor=0.997./densityEst; % assumes water 
        
        scaledMr = Mr*fitpar(n).*(1.0D0 + O2C ).^fitpar(n-1);  %scaledMr is the scaled molar mass ratio of this mixture's components.
        
        %phi2 is a scaled volume fraction, which is determined from mole fractions, scaledMr and an estimate of water/organic density ratio (rhor).
        phi2 = x2./(x2 + (1.0D0-x2).*scaledMr./rhor);  %and phi1 = 1 - phi2
        dphi2dx2 = (scaledMr./rhor)*(phi2./x2).^2; %the derivative of phi2 with respect to x2
        
        %compute O:C ratio dependent coefficients (c_i):
        % k1 = fitpar(n-2);
        % k2 = fitpar(n-3);
        n1 = (n-2)/4; %(n)/2
        coeff(1:n1) = fitpar(1:n1).*exp(fitpar(n1+1:2*n1).*O2C) + fitpar(2*n1+1:3*n1).*exp(fitpar(3*n1+1:4*n1).*Mr);
                
        sum1 = 0.0D0;
        for i = 1:n1
            sum1 = sum1 + coeff(i).*(1.0D0-2.0D0*phi2).^(i-1);
        end
        sum2 = 0.0D0;
        for i = 2:n1 % first derivative 
            sum2 = sum2 + 2.0D0.*real(i-1).*coeff(i).*(1.0D0-2.0D0*phi2).^(i-2);
        end

        gemix(:,loop_i) = phi2.*(1.0D0-phi2).*sum1;
        dgemix2dx(:,loop_i) = ( (1.0D0-2.0D0*phi2).*sum1 +phi2.*(1.0D0-phi2).*(-sum2) ).*dphi2dx2;
        
    end
    
    
end

gemix=gemix(:,1).*weight_1phase+gemix(:,2).*weight_2phase;
dgemix2dx=dgemix2dx(:,1).*weight_1phase+dgemix2dx(:,2).*weight_2phase;


%calculate the function value funcx1 (= y1(x2)) at point with w2:
ln_func1 = gemix -x2.*dgemix2dx; %the func value for component 1 = LOG(activity coeff. water)
ln_func2 = gemix +(1.0D0-x2).*dgemix2dx; %the func value of the component 2 = LOG(activity coeff. of the organic)

% ycalc1 = exp(ln_func1);
% ycalc2 = exp(ln_func2);
ycalc1 = exp_wlimiter(ln_func1); % exp with limiter for numerical overflow
ycalc2 = exp_wlimiter(ln_func2);

activity_calc1=ycalc1.*(1.0D0-x2);
activity_calc2=ycalc2.*(x2);

mass_fraction1=(1.0D0-x2).*Mr_massfrac_final./((1.0D0-x2).*(Mr_massfrac_final-1)+1);
mass_fraction2=1-mass_fraction1;
Gibbs_RT=gemix;
dGibbs_RTdx2=dgemix2dx;




