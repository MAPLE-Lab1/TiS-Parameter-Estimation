function dy = LiTank_Dyn(t,y,p)

global F R T alphaaa alphaac alphaca alphacc ratetest

% Required Expressions
iapp = p(24)*ratetest; 
ap = 3/p(4)*(1-p(12)-p(14)); 
an = 3/p(3)*(1-p(11)-p(13)); 

% Interfacial Concentration in terms of Tank Averages
c12 = (p(14)^p(6)/p(19)*y(1)+p(15)^p(7)/p(20)*y(2))/(p(14)^p(6)/p(19)+p(15)^p(7)/p(20)); 
c23 = (p(13)^p(5)/p(18)*y(3)+p(15)^p(7)/p(20)*y(2))/(p(13)^p(5)/p(18)+p(15)^p(7)/p(20)); 

% Electrolyte Diffusivities at Interfacial Concentration
Dl_c12 = p(25);  % Electrolyte_Diffusivity(p,T,c12);
Dl_c23 = p(25);  % Electrolyte_Diffusivity(p,T,c23); 

% Ionic Conductivities in terms of Interfacial Concentration
kappa_c12 = Ionic_Conductivity(p,T,c12);
kappa_c23 = Ionic_Conductivity(p,T,c23); 

% Mass Tansport expressions in Tank Model (m/s)
K12 = 2*Dl_c12/(p(19)/(p(14)^p(6))+p(20)/(p(15)^p(7))); 
K23 = 2*Dl_c23/(p(18)/(p(13)^p(5))+p(20)/(p(15)^p(7))); 

% Electrolyte Volume expressions in Tank Model (per unit area)
V_l1 = p(14)*p(19); 
V_l2 = p(15)*p(20); 
V_l3 = p(13)*p(18); 

% Transport expressions in Electrolyte Current Balance for Tank Model (S/m^2)
S12 = 2*kappa_c12/(p(19)/(p(14)^p(6))+p(20)/(p(15)^p(7))); 
S23 = 2*kappa_c23/(p(18)/(p(13)^p(5))+p(20)/(p(15)^p(7))); 

% Open Circuit Potential for Positive Electrode
Up = OCV_pos(y(12));

% Open Circuit Potential for Negative Electrode
Un = OCV_neg(y(13));

% Over Potential for Positive and Negative Electrode
eta_p=y(10)-log(y(7))-Up; 
eta_n=y(11)-log(y(9))-Un; 

% Tank Model Equations (DAE's)
dy(1,1) = K12/V_l1*(y(2)-y(1))+ap*(1-p(23))*p(17)*(p(8)*y(1))^alphacc*p(10)*(1-y(12))^alphaca*y(12)^alphacc*(exp(alphaca*F/R/T*eta_p)-exp(-alphacc*F/R/T*eta_p))/p(14)/p(8); 
dy(2,1) = (-K12*(y(2)-y(1))+K23*(y(3)-y(2)))/V_l2; 
dy(3,1) = -K23/V_l3*(y(3)-y(2))+an*(1-p(23))*p(16)*(p(8)*y(3))^alphaaa*p(9)*(1-y(13))^alphaaa*y(13)^alphaac*(exp(alphaaa*F/R/T*eta_n)-exp(-alphaac*F/R/T*eta_n))/p(13)/p(8); 
dy(4,1) = -3*iapp/ap/p(19)/F/p(10)/p(4); 
dy(5,1) = 3*iapp/an/p(18)/F/p(9)/p(3); 
dy(6,1) = 1/3600*abs(iapp); 
dy(7,1) = -S12*(log(y(8))-log(y(7)))+2*R*T*(1-p(23))/F*S12*(y(2)-y(1))/c12-iapp; 
dy(8,1) = (p(14)^p(6)/p(19)*log(y(7))+p(15)^p(7)/p(20)*log(y(8)))/(p(14)^p(6)/p(19)+p(15)^p(7)/p(20)); 
dy(9,1) = -S23*(log(y(9))-log(y(8)))+2*R*T*(1-p(23))/F*S23*(y(3)-y(2))/c23-iapp; 
dy(10,1) = p(17)*(p(8)*y(1))^alphacc*p(10)*(1-y(12))^alphaca*y(12)^alphacc*(exp(alphaca*F/R/T*eta_p)-exp(-alphacc*F/R/T*eta_p))-iapp/ap/p(19)/F; 
dy(11,1) = p(16)*(p(8)*y(3))^alphaaa*p(9)*(1-y(13))^alphaaa*y(13)^alphaac*(exp(alphaaa*F/R/T*eta_n)-exp(-alphaac*F/R/T*eta_n))+iapp/an/p(18)/F; 
dy(12,1) = p(2)*p(10)*(y(12)-y(4))/p(4)+1/5*iapp/ap/p(19)/F; 
dy(13,1) = p(1)*p(9)*(y(13)-y(5))/p(3)-1/5*iapp/an/p(18)/F; 
dy(14,1) = y(14)-y(10)+y(11); 
end

function Dl = Electrolyte_Diffusivity(p,T,c)
Dl = 1e-4*10^(-4.43-54/(T-(229+0.005*p(8)^2*c))-2.2e-4*p(8)^2*c); 
end

function kappa = Ionic_Conductivity(p,T,c)
kappa = 0.1*0.001*p(8)*c*((-10.5+0.0740*T-6.96e-5*T^2)+0.001*p(8)*c*(0.668-0.0178*T+2.8e-5*T^2)+1e-6*(p(8)*c)^2*(0.494-8.86e-4*T))^2;
end

function Up = OCV_pos(socp)
Up = -10.72*socp^4+23.88*socp^3-16.77*socp^2+2.595*socp+4.563;
end

function Un = OCV_neg(socn)
Un = .1493+.8493*exp(-61.79*socn)+.3824*exp(-665.8*socn)-exp(39.42*socn-41.92)-.3131e-1*atan(25.59*socn-4.099)-.9434e-2*atan(32.49*socn-15.74);
end



% y(1) : Electrolyte Concentrationn in positive electrode (c1)
% y(2) : Electrolyte Concentrationn in Separator          (c2)
% y(3) : Electrolyte Concentrationn in negative electrode (c3) 
% y(4) : Average solid phase concentration in positive electrode (c1,s,avg)
% y(5) : Average solid phase concentration in negative electrode (c3,s,avg)
% y(6) : Equation for charge withdrawn while discharging (or discharge capacity)
% y(7) : Electrolyte potential in positive electrode (phi,l,1)
% y(8) : Electrolyte potential in separator (phi,l,2)
% y(9) : Electrolyte potential in negative electrode (phi,l,3)
% y(10): Electrode potential in positive electrode (phi,s,1) 
% y(11): Electrode potential in positive electrode (phi,s,3)
% y(12): Solid phase surface concentration in positive electrode (csurf,s,1) 
% y(13): Solid phase surface concentration in positive electrode (csurf,s,3) 
% y(14): Cell potential (pot(t)) not given in paper.
