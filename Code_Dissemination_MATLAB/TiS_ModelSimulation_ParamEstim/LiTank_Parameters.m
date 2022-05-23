function p = LiTank_Parameters

global F R T alphaca alphacc alphaaa alphaac socpmax socnmin
global Qnominal

% DOI :  https://doi.org/10.1115/1.4028154
% Article Reference : A Temperature Dependent, Single Particle, Lithium Ion Cell Model Including Electrolyte Diffusion
% T. R. Tanim, C. D. Rahn, and C. Y. Wang, J. Dyn. Syst. Meas. Control. Trans. ASME, 137, 011005 (2015)

% Cell Capacity 
Qnominal = 1.78;

% Constants
F = 96487;       % Faraday's Constant
R = 8.314472;    % Universal Gas Constant
T = 273.15+25; % Surface Temperature of the Cell

alphaca=0.5;
alphacc=0.5;
alphaaa=0.5;
alphaac=0.5;

socpmax = 0.955473;
socpmin = 0.359749;
socnmin = 0.01;
socnmax = 0.790813;

% Model parameters
p(1) = log(1.4E-14);      % Solid phase diffusion coefficient in negative electrode (m^2/s)        (Dsn) 
p(2) = log(2E-14);   	  % Solid phase diffusion coefficient in positive electrode (m^2/s)        (Dsp)
p(3) = 5E-6; 			  % Radius of the particle in negative electrode (m)                       (Rpn)
p(4) = 5E-6; 			  % Radius of the particle in positive electrode (m)                       (Rpp)
p(5) = 1.5; 			  % bruggemann coefficient negative electrode                              (brugn)
p(6) = 1.5; 			  % bruggemann coefficient positive electrode                              (brugp)
p(7) = 1.5; 			  % bruggemann coefficient separator                                       (brugs)
p(8) = 1200; 			  % Initial electrolyte concentration           (mol/m^3)                  (c0)
p(9) = 31080;             % Maximum particle phase concentration in negative electrode(mol/cm^3)   (ctn)
p(10) = 51830;            % Maximum particle phase concentration in positive electrode(mol/m^3)    (ctp)
p(11) = 0.038;   		  % Filler fraction in negative electrode                                  (efn)
p(12) = 0.12;   		  % Filler fraction in positive electrode                                  (efp)
p(13) = 0.30; 			  % Porosity in negative electrode                                         (en)
p(14) = 0.30; 			  % Porosity in positive electrode                                         (ep)
p(15) = 0.40; 			  % porosity in separator                                                  (es)
p(16) = 6.626E-10;	      % Reaction rate constant in negative electrode (m^2.5/mol^0.5 s)         (kn)
p(17) = 2.405E-10; 		  % Reaction rate constant in positive electrode (m^2.5/mol^0.5 s)         (kp)
p(18) = 40e-6;   		  % Thickness of negative electrode              (m)                       (ln)
p(19) = 36.55e-6;   	  % Thickness of positive electrode              (m)                       (lp)
p(20) = 25e-6;    		  % Thickness of separator                       (m)                       (ls) 
p(21) = 0.79079;          % Initial particle phase concentration (scaled) negative electrode       (socn)
p(22) = 0.35973;          % Initial particle phase concentration (scaled) positive electrode       (socp)
p(23) = .38; 			  % transference number                                                    (t+)
p(24) = 17.54;            % Current density at 1C-rate (A/m^2)                                     (capacity)
p(25) = 0.28241E-9;       % Electrolyte Diffusivity (m^2/s)                                        (Dl)

end