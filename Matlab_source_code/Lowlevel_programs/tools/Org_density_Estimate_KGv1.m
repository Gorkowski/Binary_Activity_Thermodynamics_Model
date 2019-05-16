function [densityEst]=Org_density_Estimate_KGv1(M, O2C, H2C, N2C)
%****************************************************************************************
%*   :: Purpose ::                                                                      *
%*   Function to estimate the density of organic compounds based on the simple model by *
%*   Girolami (1994). The input parameters include molar mass, O:C and H:C ratios.      *
%*   If the H:C ratio is unknown at input, enter a negative value. The actual H:C will  *
%*   then be estimated based on an initial assumption of H:C = 2. The model also        *
%*   estimates the number of carbon atoms per molecular structure based on molar mass,  *
%*   O:C and H:C. The density is then approximated by the formula of Girolami.          *
%*                                                                                      *
%*   :: Author & Copyright ::                                                           *
%*   Andi Zuend,                                                                        *
%*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
%*   converted to Matlab by Kyle Gorkowski                                                                                   *
%*   -> created:        2018/03/01                                                      *
%*   -> latest changes: 2018/03/01                                                      *
%*  Girolami, G. S.: A Simple “Back of the Envelope” Method for Estimating the Densities and Molecular Volumes of Liquids and Solids, J. Chem. Educ., 71(11), 962, doi:10.1021/ed071p962, 1994.                                                                                    *
%**************************************************************************************** 

if not(exist('N2C'))
    N2C=M.*0;
end
if not(exist('H2C'))
    H2C=M.*0;
end


% densityEst(M, O2C, H2C) %H:C input is optional (when known, otherwise enter a negative value)
MC = 12.01D0; MO = 16.0D0; MH = 1.008D0; MN = 14.0067D0; %the molar masses of the carbon, oxygen and hydrogen atoms in [g/mol]

% M=18.016D0/Mr;


%1) estimate the H2C value if not provided from input
if (H2C < 0.1D0) 
    %estimate H2C assuming an aliphatic compound with H2C = 2 in absence of oxygen functional groups, then correct for oxygen content assuming a -1 slope (Van Krevelen Diagram of typical SOA).
    H2Cest = 2.0D0 -O2C;
else
    H2Cest = H2C;
end

%2) compute the approx. number of carbon atoms per organic molecule
NC = M./(MC + H2Cest.*MH + O2C.*MO + N2C.*MN);

%3) compute density estimate based on method by Girolami (1994)
%   here no correction is applied for rings and aromatic compounds (due to limited info at input)
rho1 = M./(5.0D0.*NC.*(2.0D0 + H2Cest + O2C.*2.0D0 + N2C.*2.0D0));
densityEst = rho1.*(1.0D0 + min(NC.*O2C.*0.1D0+NC.*N2C.*0.1D0, 0.3D0)) ;%density in [g/cm^3]; here it is scaled assuming that most of the oxygen atoms are able to make H-bonds (donor or acceptor).





