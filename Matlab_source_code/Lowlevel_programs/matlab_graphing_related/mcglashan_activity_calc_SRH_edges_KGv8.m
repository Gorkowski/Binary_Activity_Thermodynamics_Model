function [ycal_waterL, ycalc_orgL, activity_waterL, activity_orgL,ycal_waterU, ycalc_orgU, activity_waterU, activity_orgU,mole_frac_fit_lowRH, mole_frac_fit_upperRH ]=...
    mcglashan_activity_calc_SRH_edges_KGv8(mole_frac, O2C_values, H2C_values, molarmass_ratio, McGlashan_mode, refinement_mode, fit_tolerance);

%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2018-Jul-27  1:48 PM
% Copyright 2018 Kyle Gorkowski
%% McGlashan optional solver

mole_frac_grid=[1:-0.0001:0]';
[func1, func2, ycal_water, ycalc_org, activity_water, activity_org, mass_fraction1, mass_fraction2,Gibbs_RT, dGibbs_RTdx2]...
    =BAT_properties_calculation_v1(mole_frac_grid, O2C_values, H2C_values, molarmass_ratio, McGlashan_mode,[]);

if strcmpi(refinement_mode,'a_water SRH')
    activity_to_match=activity_water;
elseif strcmpi(refinement_mode,'a_org SRH')
    activity_to_match=activity_org;
end

Lg=length(mole_frac_grid);

    
    [lower_index]=find(activity_to_match>1);
    if isempty(lower_index)
        lower_index=Lg-1;
        upper_index=Lg-1;
    elseif lower_index(1,1)==Lg
        lower_index=Lg-1;
        upper_index=Lg-1;
        
    else
        [upper_index]=find((activity_to_match(lower_index(1,1)+1:end,1)<1));
%         disp(num2str([isempty(upper_index),isempty(lower_index)]))
        if isempty(upper_index) || upper_index(1,1)==Lg-lower_index(1,1)
            upper_index=Lg-1;
        else
            upper_index=upper_index(1,1)+lower_index(1,1);
        end
    end
    
  
    
    mole_frac_guess(1,1)=mole_frac_grid(lower_index(1,1),1);
    mole_frac_guess(2,1)=mole_frac_grid(upper_index(1,1),1);
    
    mole_frac_bounds_guess(1,:)=[mole_frac_grid(lower_index(1,1)+1,1),mole_frac_grid(lower_index(1,1)-1,1)];
    mole_frac_bounds_guess(2,:)=[mole_frac_grid(upper_index(1,1)+1,1),mole_frac_grid(upper_index(1,1)-1,1)];
    
%         disp(num2str([activity_water([upper_index(1,1)+1,upper_index(1,1),upper_index(1,1)-1],1)]))

            mole_frac_fit_lowRH=mole_frac.*NaN;
        mole_frac_fit_upperRH=mole_frac.*NaN;
                error=mole_frac.*NaN;
    for i=1:2
        

        

        mole_frac_bounds=mole_frac_bounds_guess(i,:);
        mole_frac_guess_temp=mole_frac_guess(i,1);
        
        % fmincon opt
        options=optimoptions('fmincon');
        options.Display='off';
        options.Algorithm='interior-point'; % 'interior-point' (default) 'trust-region-reflective' 'sqp' 'active-set'
        options.TolX=fit_tolerance;
        options.MaxIter=1000;
        
        
        
        if strcmpi(refinement_mode,'a_water SRH')
            problem = createOptimProblem('fmincon','objective',...
                @(mole_frac)nested_opt_cost_sub2(mole_frac, O2C_values, H2C_values, molarmass_ratio, McGlashan_mode),...
                'x0',mole_frac_guess_temp,'lb',mole_frac_bounds(1,1),'ub',mole_frac_bounds(1,2),'options',options);
        elseif strcmpi(refinement_mode,'a_org SRH')
            problem = createOptimProblem('fmincon','objective',...
                @(mole_frac)nested_opt_cost_sub3(mole_frac, O2C_values, H2C_values, molarmass_ratio, McGlashan_mode),...
                'x0',mole_frac_guess_temp,'lb',mole_frac_bounds(1,1),'ub',mole_frac_bounds(1,2),'options',options);
        end

        [mole_frac_fit,error(i,1),~] = fmincon(problem);
        
        if i==1
            mole_frac_fit_lowRH(1,1)=mole_frac_fit;
        else
            mole_frac_fit_upperRH(1,1)=mole_frac_fit;
        end
        
        %         Ssol=size(allmins);
        %         mole_frac_temps=zeros(Ssol(1,2),1);
        %         for i_sol=1:Ssol(1,2)
        %             mole_frac_temps(i_sol,1)=allmins(1,i_sol).X;
        %         end
        %         mole_frac_temps=sort(mole_frac_temps,'descend');
        %
        %         if Ssol(1,2)==1
        %             mole_frac_fit_lowRH(i,1)=mole_frac_temps(1,1);
        %         else
        %             mole_frac_fit_lowRH(i,1)=mole_frac_temps(1,1);
        %             mole_frac_fit_upperRH(i,1)=mole_frac_temps(2,1);
        %         end
        %         mole_frac_fit_count(i,1)=Ssol(1,2);
    end
    
    
    if mole_frac_fit_upperRH==mole_frac_fit_lowRH
        mole_frac_fit_lowRH=0;
        mole_frac_fit_upperRH=0;
    end
    
    
    [func1L, func2L, ycal_waterL, ycalc_orgL, activity_waterL, activity_orgL, mass_fraction1L, mass_fraction2L,Gibbs_RTL, dGibbs_RTdx2L]...
        =BAT_properties_calculation_v1(mole_frac_fit_lowRH, O2C_values, H2C_values, molarmass_ratio, McGlashan_mode,[]);
    
    [func1U, func2U, ycal_waterU, ycalc_orgU, activity_waterU, activity_orgU, mass_fraction1U, mass_fraction2U,Gibbs_RTU, dGibbs_RTdx2U]...
        =BAT_properties_calculation_v1(mole_frac_fit_upperRH, O2C_values, H2C_values, molarmass_ratio, McGlashan_mode,[]);
    
    %  mole_frac_fit_upperRH=mole_frac_fit_count;




end


function [error_out]= nested_opt_cost_sub1(mole_frac, O2C_values, H2C_values, molarmass_ratio, McGlashan_mode, aw_desired_temp);
% disp(num2str(mole_frac));

[func1, func2, ycal_water, ycalc_org, activity_water, activity_org, mass_fraction1, mass_fraction2,Gibbs_RT, dGibbs_RTdx2]...
    =BAT_properties_calculation_v1(mole_frac, O2C_values, H2C_values, molarmass_ratio, McGlashan_mode,[]);


error_out=(activity_water-aw_desired_temp)^2;

end


%% water phase sep
function [error_out]= nested_opt_cost_sub2(mole_frac, O2C_values, H2C_values, molarmass_ratio, McGlashan_mode);
% disp(num2str(mole_frac));

[func1, func2, ycal_water, ycalc_org, activity_water, activity_org, mass_fraction1, mass_fraction2,Gibbs_RT, dGibbs_RTdx2]...
    =BAT_properties_calculation_v1(mole_frac, O2C_values, H2C_values, molarmass_ratio, McGlashan_mode,[]);

error_out=(activity_water-1)^2;

end


%% organic phase sep
function [error_out]= nested_opt_cost_sub3(mole_frac, O2C_values, H2C_values, molarmass_ratio, McGlashan_mode);
% disp(num2str(mole_frac));

[func1, func2, ycal_water, ycalc_org, activity_water, activity_org, mass_fraction1, mass_fraction2,Gibbs_RT, dGibbs_RTdx2]...
    =BAT_properties_calculation_v1(mole_frac, O2C_values, H2C_values, molarmass_ratio, McGlashan_mode,[]);

error_out=(activity_org-1)^2;

end


